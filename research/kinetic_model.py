"""Homogeneous Kantowski--Sachs Einstein--Vlasov benchmark, c=8*pi*G=1.

This is a symmetry-reduced classical model of a trapped, contracting region;
it is NOT a matched, asymptotically flat black-hole formation simulation.
Initial isotropic momentum distributions are smooth positive mixtures of
lognormal radial-number distributions. Their shape is transported exactly
along conserved covariant radial momenta and squared angular momenta.
Only the momentum integrals and Einstein ODEs are numerically approximated.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from numpy.polynomial.hermite import hermgauss
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp
from scipy.optimize import linprog

@dataclass
class MomentumBasis:
    centers: np.ndarray
    width: float = .10
    mass: float = 1.
    radial_order: int = 12
    angular_order: int = 32

    def __post_init__(self):
        self.centers = np.asarray(self.centers, dtype=float)
        if np.any(self.centers <= 0) or self.width <= 0 or self.mass < 0:
            raise ValueError('Centers/width positive; rest mass nonnegative.')
        x, rw = hermgauss(self.radial_order)
        mu, aw = leggauss(self.angular_order)
        self.p = self.centers[:, None] * np.exp(np.sqrt(2.) * self.width * x)
        self.rw = rw / np.sqrt(np.pi)
        self.mu = mu
        self.aw = aw / 2.
        p2 = self.p[:, :, None]**2
        self.radial_conserved2 = p2 * mu[None, None, :]**2
        self.angular_conserved2 = p2 * (1 - mu[None, None, :]**2)
        self.quadw = self.rw[:, None] * self.aw[None, :]

    def moment_rows(self, maximum: int = 5) -> np.ndarray:
        """Rows: N, M0=rho, M1=3p, M2=integral p^4/E^3, ... at a=b=1."""
        e = np.sqrt(self.mass**2 + self.p**2)
        rows = [np.ones(len(self.centers))]
        for k in range(maximum + 1):
            rows.append((self.p**(2*k) / e**(2*k-1)) @ self.rw)
        return np.array(rows)

    def tensors(self, y: np.ndarray, weights: np.ndarray, order: int = 1) -> dict:
        la, lb, ha, hb = y
        x = self.radial_conserved2 * np.exp(-2*la)
        z = self.angular_conserved2 * np.exp(-2*lb)
        e = np.sqrt(self.mass**2 + x + z)
        weight = weights[:, None, None] * self.quadw[None, :, :] * np.exp(-la-2*lb)
        ans = {}
        for k in range(order + 1):
            for i in range(k + 1):
                j = k-i
                ans[i,j] = np.sum(weight * x**i * z**j / e**(2*k-1))
        return ans

    def rhs(self, t: float, y: np.ndarray, weights: np.ndarray) -> np.ndarray:
        la, lb, ha, hb = y
        c = self.tensors(y, weights)
        pr, pt = c[1,0], c[0,1]/2
        dhb = -.5*(3*hb*hb + np.exp(-2*lb) + pr)
        dha = -pt - dhb - ha*ha - hb*hb - ha*hb
        return np.array([ha,hb,dha,dhb])

    def observe(self, y: np.ndarray, weights: np.ndarray) -> dict:
        la, lb, ha, hb = y
        c = self.tensors(y, weights, 2)
        rho, pr, pt = c[0,0], c[1,0], c[0,1]/2
        dy = self.rhs(0., y, weights)
        aa, bb = dy[2]+ha**2, dy[3]+hb**2
        curvature = 4*(aa**2+2*bb**2+2*(ha*hb)**2+(hb**2+np.exp(-2*lb))**2)
        scale = abs(2*ha*hb)+hb*hb+np.exp(-2*lb)+abs(rho)
        constraint = 2*ha*hb+hb**2+np.exp(-2*lb)-rho
        dpr=-(3*ha+2*hb)*pr+ha*c[2,0]+hb*c[1,1]
        dpt=-(ha+4*hb)*pt+.5*(ha*c[1,1]+hb*c[0,2])
        return dict(a=np.exp(la),b=np.exp(lb),Ha=ha,Hb=hb,rho=rho,pr=pr,pt=pt,
                    n=np.sum(weights)*np.exp(-la-2*lb),K=curvature,
                    dpr=dpr,dpt=dpt,constraint=constraint,
                    scaled_constraint=constraint/scale,
                    outgoing_expansion=2*hb,ingoing_expansion=2*hb)

    def initial_geometry(self, weights: np.ndarray, hb: float = -1.) -> np.ndarray:
        rho = self.moment_rows(0)[1] @ weights
        ha = (rho-hb*hb-1.)/(2*hb)
        return np.array([0.,0.,ha,hb])

    def evolve(self, weights: np.ndarray, start: float, end: float,
               y0: np.ndarray | None = None, rtol: float = 2e-11,
               atol: float = 2e-13, method: str = 'DOP853'):
        if np.min(weights)<0 or not np.all(np.isfinite(weights)):
            raise ValueError('Mixture weights must be nonnegative and finite.')
        if y0 is None: y0 = self.initial_geometry(weights)
        def stop(t,y): return y[1]-np.log(.10)
        stop.terminal=True; stop.direction=-1
        sol=solve_ivp(lambda t,y:self.rhs(t,y,weights),(start,end),y0,
                      rtol=rtol,atol=atol,method=method,dense_output=True,
                      events=stop,max_step=.02)
        if not sol.success: raise RuntimeError(sol.message)
        if sol.t[-1] < end-1e-10:
            raise RuntimeError('Requested time passes the finite-radius safety stop.')
        return sol


def matched_pair(basis: MomentumBasis, retained: int, density: float=.8,
                 fraction: float=.8, seed: int|None=None) -> dict:
    """Match N and M0..M_retained; change the FIRST omitted moment.

    LP maximizes the initial next-moment gap in a specified positive box.
    It never sees any future solution. This constructs a counterexample,
    not a representative population and not an optimum over all distributions.
    """
    rows=basis.moment_rows(max(5,retained+1))
    base=1/np.sqrt(1+basis.centers**2)
    if seed is not None:
        base *= np.exp(np.random.default_rng(seed).normal(0,.35,len(base)))
    base *= density/(rows[1]@base)
    constraints=rows[:retained+2]
    scaled=constraints/np.linalg.norm(constraints,axis=1)[:,None]
    target=rows[retained+2]
    result=linprog(-target/np.linalg.norm(target), A_eq=scaled,
                   b_eq=np.zeros(len(scaled)),
                   bounds=[(-fraction*v,fraction*v) for v in base],method='highs',
                   options={'dual_feasibility_tolerance':1e-9,
                            'primal_feasibility_tolerance':1e-9})
    if not result.success: raise RuntimeError(result.message)
    delta=result.x
    # Orthogonal correction removes equality roundoff, not a physical adjustment.
    delta-=np.linalg.lstsq(scaled,scaled@delta,rcond=1e-13)[0]
    wp,wm=base+delta,base-delta
    if min(wp.min(),wm.min())<=0: raise RuntimeError('Positivity failed.')
    mismatch=np.max(np.abs(constraints@(wp-wm)) / np.maximum(np.abs(constraints@base),1e-30))
    if mismatch>1e-10: raise RuntimeError('Moment matching failed.')
    return dict(plus=wp,minus=wm,base=base,retained=retained,
                relative_match_residual=mismatch,
                initial_moments_plus=rows@wp,initial_moments_minus=rows@wm,
                first_omitted_gap=float(target@(wp-wm)))


def massless_pair(basis: MomentumBasis, density: float=.8) -> dict:
    if basis.mass!=0: raise ValueError('This control requires mass=0.')
    rows=basis.moment_rows(1)
    base=1/(1+basis.centers)
    base*=density/(rows[1]@base)
    C=rows[:2]; C=C/np.linalg.norm(C,axis=1)[:,None]
    target=(basis.p**2)@basis.rw
    result=linprog(-target/np.linalg.norm(target),A_eq=C,b_eq=np.zeros(2),
                  bounds=[(-.8*v,.8*v) for v in base],method='highs')
    if not result.success:raise RuntimeError(result.message)
    d=result.x;d-=np.linalg.lstsq(C,C@d,rcond=None)[0]
    return dict(plus=base+d,minus=base-d)
