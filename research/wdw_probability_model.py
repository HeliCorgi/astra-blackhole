"""Bounded radius/sphere-curvature observables of the existing WDW model.

No 4D Kretschmann operator, changed Hamiltonian, tail repair or state clipping.
T is the chosen intrinsic clock, not proper time. Probabilities refer to the
already chosen L2(dx) representation. Numerical walls remain regulators.
"""
from __future__ import annotations
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import ndtr, expit, k0
from wdw_model import WDWModel, Packet


def spectral_states(model: WDWModel, packet: Packet, times) -> tuple:
    """Unchanged propagator with exact initial-vector recovery using expm1."""
    t = np.asarray(times, float)
    if t.ndim != 1 or not len(t) or not np.all(np.isfinite(t)) or np.any(t < 0):
        raise ValueError('Nonnegative finite one-dimensional times required.')
    if packet.h != model.h:
        raise ValueError('State and operator h differ.')
    v = packet.initial(model.x)
    c = model.U.T @ v
    states = v[:, None] + model.U @ (c[:, None] * np.expm1(
        -1j * model.energy[:, None] * t[None, :] / model.h))
    return states, v, c


class BoxDistribution:
    """Nonnegative piecewise-linear density, zero at the two Dirichlet walls.

    Input masses are |sqrt(dx)*psi|^2. No normalization is performed.
    All sub-grid threshold integrals are analytic for this interpolation;
    convergence of this interpolation to the continuum is tested separately.
    """
    def __init__(self, x, probability, left: float, right: float):
        x, p = np.asarray(x, float), np.asarray(probability, float)
        if x.ndim != 1 or p.shape != x.shape or len(x) < 2:
            raise ValueError('Matching one-dimensional grid and masses required.')
        if not np.all(np.isfinite(x)) or not np.all(np.isfinite(p)) or np.any(p < 0):
            raise ValueError('Finite grid and nonnegative finite masses required.')
        dx = x[1] - x[0]
        if dx <= 0 or not np.allclose(np.diff(x), dx, rtol=1e-10, atol=1e-12):
            raise ValueError('Uniform increasing grid required.')
        if not np.isclose(x[0] - left, dx) or not np.isclose(right - x[-1], dx):
            raise ValueError('Walls must be one grid step beyond the samples.')
        self.x = np.r_[left, x, right]
        self.f = np.r_[0., p / dx, 0.]
        self.cells = np.diff(self.x) * (self.f[:-1] + self.f[1:]) / 2
        self.survival = np.r_[np.cumsum(self.cells[::-1])[::-1], 0.]
        self.total = float(self.survival[0])
        if self.total <= 0:
            raise ValueError('Nonzero mass required.')

    def sf(self, threshold: float) -> float:
        """Integral over x >= threshold, computed from right to avoid 1-CDF loss."""
        if not np.isfinite(threshold):
            raise ValueError('Finite threshold required.')
        if threshold <= self.x[0]: return self.total
        if threshold >= self.x[-1]: return 0.
        i = int(np.searchsorted(self.x, threshold, side='right') - 1)
        width = self.x[i+1] - self.x[i]
        f_at = self.f[i] + (self.f[i+1]-self.f[i])*(threshold-self.x[i])/width
        return float(self.survival[i+1] + (self.x[i+1]-threshold)*(f_at+self.f[i+1])/2)

    def x_for_survival(self, probability: float) -> float:
        """Solve sf(x)=p without renormalizing the numerical state."""
        if not np.isfinite(probability) or not 0 < probability < self.total:
            raise ValueError('Probability must be strictly inside (0,total mass).')
        return float(brentq(lambda x:self.sf(x)-probability,
                           self.x[0], self.x[-1], xtol=2e-12))


def radius_threshold_x(radius_ratio: float, time: float, x0=2.) -> float:
    if not np.isfinite(radius_ratio) or radius_ratio <= 0 or not np.isfinite(time):
        raise ValueError('Positive finite radius ratio and finite time required.')
    return float(x0-time-np.log(radius_ratio))


def capped_sphere_response(x, probability, threshold_x: float) -> float:
    """E[S/(S+S_*)] where S=2/r^2 and S_*=2/r(threshold_x)^2.

    This is a bounded response, not a claim of a detector model or curvature
    saturation. The response is applied to an unchanged state.
    """
    x, p = np.asarray(x, float), np.asarray(probability, float)
    if x.ndim != 1 or not np.all(np.isfinite(x)) or x.shape != p.shape or np.any(p < 0) or not np.all(np.isfinite(p)) or not np.isfinite(threshold_x):
        raise ValueError('Finite matching data required.')
    return float(p @ expit(2*(x-threshold_x)))


def classical_radius_probability(packet: Packet, time: float, radius_ratio: float,
                                 epsabs=2e-15, epsrel=2e-10) -> tuple[float,float]:
    """Continuum Gaussian-Wigner Liouville comparison, reduced to a 1-D integral.

    x(T)=xi+log(cosh T+v sinh T), v=p/sqrt(p^2+exp(-2xi)).
    Conditional momentum is Gaussian. Integrate its analytic survival function
    over xi; no discontinuous indicator on Gauss-Hermite nodes is used.
    This is classical Hamilton flow, not a second quantum solver.
    """
    if not np.isfinite(time) or time < 0:
        raise ValueError('Finite nonnegative time required.')
    a = radius_threshold_x(radius_ratio, time, packet.x0)
    if time == 0:
        return float(ndtr((packet.x0-a)/packet.sigma)), 0.
    # Reparameterize the bounded xi interval by v=tanh(u). This removes
    # the endpoint square-root stiffness that triggered a QUADPACK warning
    # in the first direct-xi implementation at very small T.
    logsinh_t = time + np.log(-np.expm1(-2*time)) - np.log(2.)
    def logcosh(u): return np.logaddexp(u,-u)-np.log(2.)
    sp = packet.h/(2*packet.sigma)
    def integrand(u):
        lc,lt=logcosh(u),logcosh(u+time)
        xi=a-lt+lc
        logjac=logsinh_t-lc-lt
        z=(xi-packet.x0)/packet.sigma
        logdensity=-z*z/2-np.log(packet.sigma)-.5*np.log(2*np.pi)
        if logjac+logdensity < -745: return 0.  # float underflow, no tail renormalization
        if u==0: pcrit=0.
        else:
            au=abs(u); lp=-xi+au+np.log(-np.expm1(-2*au))-np.log(2.)
            if lp>700: return 0. if u>0 else float(np.exp(logdensity+logjac))
            pcrit=np.sign(u)*np.exp(lp)
        pm=packet.p0+packet.chirp*(xi-packet.x0)
        return float(np.exp(logdensity+logjac)*ndtr((pm-pcrit)/sp))
    val=0.;error=0.
    for lo,hi in [(-np.inf,-time),(-time,0.),(0.,np.inf)]:
        v,e=quad(integrand,lo,hi,epsabs=epsabs/3,epsrel=epsrel,limit=400)
        val+=v;error+=e
    return float(val+ndtr((packet.x0-a-time)/packet.sigma)),float(error)


def continuum_overlap(packet: Packet, width=12., epsabs=2e-13) -> tuple[complex,float]:
    """Independent K0 threshold-overlap integral, not fitted to propagated tails."""
    lo,hi=packet.x0-width*packet.sigma,packet.x0+width*packet.sigma
    def f(x):
        z=x-packet.x0
        psi=(2*np.pi*packet.sigma**2)**(-.25)*np.exp(-z*z/(4*packet.sigma**2)
                 +1j*(packet.p0*z+packet.chirp*z*z/2)/packet.h)
        return k0(np.exp(-x)/packet.h)*psi
    re,er=quad(lambda x:float(f(x).real),lo,hi,epsabs=epsabs,epsrel=2e-12,limit=300)
    im,ei=quad(lambda x:float(f(x).imag),lo,hi,epsabs=epsabs,epsrel=2e-12,limit=300)
    return complex(re,im),float(er+ei)


def asymptotic_survival(x: float, time: float, h: float, overlap: complex) -> float:
    """Leading A/(5X^5), conditional on the earlier fixed-T density asymptotic.

    This is NOT an error bound or a replacement for the missing numerical tail.
    """
    X=x+np.log(2*h)-np.euler_gamma
    if X <= 0 or h <= 0 or time < 0: raise ValueError('Outside asymptotic coordinates.')
    return float(16*time*time*abs(overlap)**2/(5*np.pi**2*X**5))
