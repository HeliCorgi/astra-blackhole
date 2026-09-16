"""Massless Einstein--Boltzmann relaxation benchmark, c=8*pi*G=1.

Homogeneous Kantowski--Sachs, even axisymmetric one-particle distribution.
Anderson--Witting-type, momentum-INDEPENDENT proper-time rate; equilibrium
matches number and energy. This is an RTA ansatz, not a binary collision
integral or a measured black-hole material model. Only energy-integrated
angular brightness and number are evolved. Quantum effects are absent.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
import numpy as np
from numpy.polynomial import legendre as leg
from scipy.integrate import solve_ivp


def transport_matrix(lmax: int) -> np.ndarray:
    """Exact polynomial projection up to floating arithmetic, even Legendre modes.

    A P_l = mu*(1-mu**2)*P_l' - 4*mu**2*P_l.
    Drop only output modes above lmax; convergence must be tested.
    """
    if not isinstance(lmax, int) or lmax < 8 or lmax % 2:
        raise ValueError('lmax must be even and >= 8.')
    out = np.zeros((lmax // 2 + 1,) * 2)
    for j, ell in enumerate(range(0, lmax + 1, 2)):
        diagonal = F(ell**2, 2*ell-1) - (ell+4)*(F((ell+1)**2,(2*ell+1)*(2*ell+3)) + F(ell**2,(2*ell+1)*(2*ell-1)))
        out[j,j] = float(diagonal)
        if j > 0:
            out[j-1,j] = float(F(ell*(ell-1)*(ell-3),(2*ell-1)*(2*ell+1)))
        if j+1 < len(out):
            out[j+1,j] = float(F(-(ell+4)*(ell+1)*(ell+2),(2*ell+1)*(2*ell+3)))
    return out


@dataclass
class RelaxationModel:
    lmax: int = 128
    rate: float = 0.

    def __post_init__(self):
        if not np.isfinite(self.rate) or self.rate < 0:
            raise ValueError('rate must be finite and nonnegative.')
        self.A = transport_matrix(self.lmax)
        self.ells = np.arange(0, self.lmax + 1, 2)

    def initial(self, ell=4, amplitude=.8, density=.8, number=.8/3):
        if ell not in self.ells or ell < 4 or not np.isfinite(amplitude) or abs(amplitude) >= 1:
            raise ValueError('Even retained ell>=4 and |amplitude|<1 required.')
        if not all(np.isfinite(v) and v > 0 for v in (density, number)):
            raise ValueError('Positive density and number required.')
        c = np.zeros(len(self.ells)); c[0] = density; c[ell//2] = density*amplitude
        hb = -1.; ha = (density - hb*hb - 1.)/(2*hb)
        return np.r_[0., 0., ha, hb, number, c]

    def rhs(self, t, y):
        la, lb, ha, hb, n = y[:5]; c = y[5:]
        rho = c[0]; pr = rho/3 + 2*c[1]/15; pt = (rho-pr)/2
        dhb = -.5*(3*hb*hb + np.exp(-2*lb) + pr)
        dha = -pt-dhb-ha*ha-hb*hb-ha*hb
        dc = -4*hb*c + (ha-hb)*(self.A @ c)
        dc[1:] -= self.rate*c[1:]  # c0 = energy, exactly collision invariant.
        return np.r_[ha, hb, dha, dhb, -(ha+2*hb)*n, dc]

    def observe(self, y):
        la, lb, ha, hb, n = y[:5]; c = y[5:]
        rho = c[0]; pr = rho/3+2*c[1]/15; pt = (rho-pr)/2
        dy = self.rhs(0., y)
        aa = dy[2]+ha*ha; bb = dy[3]+hb*hb; q = np.exp(-2*lb)
        K = 4*(aa*aa+2*bb*bb+2*(ha*hb)**2+(hb*hb+q)**2)
        C = 2*ha*hb+hb*hb+q-rho
        return dict(a=float(np.exp(la)), b=float(np.exp(lb)), Ha=float(ha), Hb=float(hb),
                    rho=float(rho), pr=float(pr), pt=float(pt), n=float(n), K=float(K),
                    scaled_constraint=float(C/(abs(2*ha*hb)+hb*hb+q+abs(rho))),
                    pressure_anisotropy=float((pr-pt)/rho),
                    comoving_number=float(n*np.exp(la+2*lb)),
                    deformation_over_rate=None if self.rate==0 else float(max(abs(ha),abs(hb))/self.rate))

    def brightness(self, y, mu):
        coeff = np.zeros(self.lmax+1); coeff[::2] = y[5:]
        return leg.legval(mu, coeff)

    def evolve(self, ell=4, amplitude=.8, end=.45, rtol=2e-10, atol=2e-12,
               method='Radau', y0=None):
        if not np.isfinite(end) or end <= 0:
            raise ValueError('end must be finite and positive.')
        if y0 is None: y0 = self.initial(ell, amplitude)
        if np.shape(y0) != (5+len(self.ells),) or not np.all(np.isfinite(y0)):
            raise ValueError('Invalid initial state.')
        def stop(t,y): return y[1]-np.log(.1)
        stop.terminal=True; stop.direction=-1
        sol = solve_ivp(self.rhs,(0,end),y0,method=method,rtol=rtol,atol=atol,
                        max_step=.02,dense_output=True,events=stop)
        if not sol.success or sol.t[-1] < end-1e-10:
            raise RuntimeError(sol.message+' or finite-radius stop reached.')
        return sol


def fluid_evolve(end=.45, rtol=2e-11, atol=2e-13, density=.8, number=.8/3):
    """Separate perfect-radiation reference, not an imposed RTA solution."""
    def rhs(t,y):
        la,lb,ha,hb,n,rho=y
        pr=pt=rho/3
        dhb=-.5*(3*hb*hb+np.exp(-2*lb)+pr)
        dha=-pt-dhb-ha*ha-hb*hb-ha*hb
        return np.array([ha,hb,dha,dhb,-(ha+2*hb)*n,-4*(ha+2*hb)*rho/3])
    sol=solve_ivp(rhs,(0,end),[0,0,(2-density)/2,-1,number,density],rtol=rtol,atol=atol,
                  method='DOP853',dense_output=True,max_step=.01)
    if not sol.success: raise RuntimeError(sol.message)
    return sol


def symmetric_gap(a,b):
    return 2*abs(a-b)/(abs(a)+abs(b))
