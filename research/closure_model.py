"""Zero-tail moment closures of the existing classical RTA model.
Copyright 2026 HeliCorgi. Apache-2.0.

The independent comparison target is the high-order collision_model solver.
Low-order positive-moment feasibility does not imply that its zero-tail
polynomial is positive, nor that its time evolution approximates nature.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from numpy.polynomial import legendre as leg
from scipy.integrate import solve_ivp
from collision_model import transport_matrix


def initial_state(lmax: int, density: float, hb: float, coefficients: dict) -> np.ndarray:
    """Shared, constraint-satisfying initial geometry; only resolved modes copied."""
    if lmax < 0 or lmax % 2 or not np.isfinite(density) or density <= 0:
        raise ValueError('Nonnegative even order and positive finite density required.')
    if not np.isfinite(hb) or hb >= 0:
        raise ValueError('This audit requires finite Hb<0.')
    modes = {int(k): float(v) for k, v in coefficients.items()}
    if any(k < 4 or k % 2 or not np.isfinite(v) for k, v in modes.items()):
        raise ValueError('Only finite even modes >=4 are admitted initially.')
    if sum(abs(v) for v in modes.values()) >= 1:
        raise ValueError('The chosen sufficient initial positivity bound must be strict.')
    c = np.zeros(lmax//2+1); c[0] = density
    for ell, value in modes.items():
        if ell <= lmax: c[ell//2] = density*value
    ha = (density-hb*hb-1)/(2*hb)
    if ha+2*hb >= 0:
        raise ValueError('Initial volume must contract.')
    return np.r_[0., 0., ha, hb, density/3, c]


def moment_diagnostics(c: np.ndarray) -> dict:
    """Finite Hausdorff moment checks for x=mu^2 in [0,1].

    Up to c4: E[x]^2 <= E[x^2] <= E[x], 0<=E[x]<=1.
    These characterize a positive measure for these retained moments;
    they do NOT prove a smooth density or a positive zero-tail polynomial.
    """
    c = np.asarray(c, dtype=float)
    if c.ndim != 1 or len(c) not in (1,2,3) or not np.all(np.isfinite(c)) or c[0] <= 0:
        raise ValueError('Need positive density and 1,2,3 finite even coefficients.')
    z = c/c[0]
    m1 = 1/3+(2*z[1]/15 if len(c)>1 else 0.)
    margins = [m1, 1-m1]
    m2 = None
    if len(c) == 3:
        m2 = 1/5+4*z[1]/35+8*z[2]/315
        margins += [m2-m1*m1, m1-m2]
    poly = np.zeros(2*len(c)-1); poly[::2] = z
    candidates = [-1.,1.]
    if len(c)>1:
        for r in leg.legroots(leg.legder(poly)):
            if abs(complex(r).imag)<1e-9 and -1 <= float(np.real(r)) <= 1:
                candidates.append(float(np.real(r)))
    minimum = float(np.min(leg.legval(candidates,poly)))
    return dict(m1=float(m1), m2=None if m2 is None else float(m2),
                moment_margin=float(min(margins)), polynomial_minimum=minimum)


@dataclass
class ClosureModel:
    lmax: int
    rate: float

    def __post_init__(self):
        if self.lmax not in (0,2,4) or not np.isfinite(self.rate) or self.rate<0:
            raise ValueError('Use lmax 0,2,4 and a finite nonnegative rate.')
        size = self.lmax//2+1
        self.A = transport_matrix(8)[:size,:size].copy()

    def rhs(self, t: float, y: np.ndarray) -> np.ndarray:
        la,lb,ha,hb,n = y[:5]; c = y[5:]
        rho = c[0]; c2 = c[1] if len(c)>1 else 0.
        pr = rho/3+2*c2/15; pt=(rho-pr)/2
        dhb = -.5*(3*hb*hb+np.exp(-2*lb)+pr)
        dha = -pt-dhb-ha*ha-hb*hb-ha*hb
        dc = -4*hb*c+(ha-hb)*(self.A@c)
        dc[1:] -= self.rate*c[1:]
        return np.r_[ha,hb,dha,dhb,-(ha+2*hb)*n,dc]

    def evolve(self,y0: np.ndarray,end: float,method='Radau',rtol=2e-10,atol=2e-12):
        y0=np.asarray(y0,dtype=float)
        if y0.shape!=(6+self.lmax//2,) or not np.all(np.isfinite(y0)) or y0[5]<=0:
            raise ValueError('Invalid state shape or density.')
        if not np.isfinite(end) or end<=0: raise ValueError('Positive horizon required.')
        def stop(t,y): return y[1]-np.log(.1)
        stop.terminal=True; stop.direction=-1
        sol=solve_ivp(self.rhs,(0,end),y0,method=method,rtol=rtol,atol=atol,
                      dense_output=True,events=stop,max_step=.02)
        if not sol.success or sol.t[-1]<end-1e-10:
            raise RuntimeError(sol.message+' or finite-radius stop reached.')
        return sol


def observe_series(model,solution,times: np.ndarray) -> dict:
    """Common observables, usable for both the reference and low-order models."""
    y=solution.sol(times); la,lb,ha,hb,n=y[:5]; c=y[5:]
    rho=c[0]; c2=c[1] if len(c)>1 else np.zeros_like(rho)
    pr=rho/3+2*c2/15; pt=(rho-pr)/2; q=np.exp(-2*lb)
    dhb=-.5*(3*hb*hb+q+pr); dha=-pt-dhb-ha*ha-hb*hb-ha*hb
    K=4*((dha+ha*ha)**2+2*(dhb+hb*hb)**2+2*(ha*hb)**2+(hb*hb+q)**2)
    con=(2*ha*hb+hb*hb+q-rho)/(abs(2*ha*hb)+hb*hb+q+abs(rho))
    num=n*np.exp(la+2*lb)
    return dict(y=y,c=c,K=K,rho=rho,Pi=pr-pt,Ha=ha,Hb=hb,
                constraint=con,number_drift=num/num[0]-1,b=np.exp(lb),
                epsilon=None if model.rate==0 else np.maximum.reduce([abs(ha),abs(hb),abs(ha-hb)])/model.rate)
