"""Unbounded-observable and tail audit of the existing reduced WDW model.
Copyright 2026 HeliCorgi. Apache-2.0.

J_q = <(r_ref/r)^q>, r_ref=exp(-x0)/4. q=2 measures the intrinsic
2-sphere scalar curvature in reference units. q=6 is a normalized vacuum
curvature-per-mass-squared diagnostic, NOT a full quantum 4D Kretschmann
operator. No fixed-mass superselection, quantum mass ordering, or
spacetime singularity resolution is asserted.
"""
from __future__ import annotations
import numpy as np
from scipy.integrate import quad
from scipy.special import k0, logsumexp
from numpy.polynomial.hermite import hermgauss
from wdw_model import Packet, WDWModel, classical_path


def states_at(model: WDWModel, packet: Packet, times):
    """Same spectral propagator, using expm1 to preserve the exact input at T=0.

    U exp(-iET/h) U^T psi and this expression agree up to the computed
    orthogonality error. No time-dependent renormalization or tail clipping.
    """
    ts=np.asarray(times,dtype=float)
    if ts.ndim!=1 or not len(ts) or not np.all(np.isfinite(ts)):
        raise ValueError('Finite one-dimensional times required.')
    if model.h!=packet.h: raise ValueError('Mismatched h.')
    v=packet.initial(model.x);coef=model.U.T@v
    states=v[:,None]+model.U@(coef[:,None]*np.expm1(-1j*model.energy[:,None]*ts[None,:]/model.h))
    direct=model.U@(coef[:,None]*np.exp(-1j*model.energy[:,None]*ts[None,:]/model.h))
    return states,dict(max_norm_error=float(np.max(abs(np.sum(abs(states)**2,axis=0)-1))),
                       max_direct_state_l2_difference=float(np.max(np.linalg.norm(states-direct,axis=0))),
                       initial_roundtrip=model.U@coef)


def log_moment(x, probability, q: float, t: float, x0: float=2., cutoff=None, cap=False):
    """Natural log of a non-renormalized positive expectation / partial integral.

    probability contains quadrature masses. cutoff restricts x<=cutoff;
    cap instead clips the observable, never the state or its probability.
    A grid point exactly on a sharp cutoff gets half weight (trapezoidal
    endpoint convention). Return -inf for a zero integral. Box sums are regulators.
    """
    x=np.asarray(x,float);p=np.asarray(probability,float)
    if x.ndim!=1 or p.shape!=x.shape or not np.all(np.isfinite(x)) or not np.all(np.isfinite(p)) or np.any(p<0):
        raise ValueError('Finite matching arrays and nonnegative masses required.')
    if not np.isfinite(q) or q<=0 or not np.isfinite(t) or not np.isfinite(x0):
        raise ValueError('Positive q and finite clock/scale required.')
    if cutoff is not None and not np.isfinite(cutoff):raise ValueError('Finite cutoff required.')
    if cap and cutoff is None:raise ValueError('Cap needs a cutoff.')
    mask=np.ones(x.shape,bool) if cutoff is None or cap else x<=cutoff+1e-10
    xx=np.minimum(x,cutoff) if cap else x
    weights=p[mask].copy()
    if cutoff is not None and not cap:
        weights[np.isclose(x[mask],cutoff,atol=1e-10,rtol=0)]*=0.5
    with np.errstate(divide='ignore'):lp=np.log(weights)
    return float(logsumexp(lp+q*(xx[mask]+t-x0)))


def continuum_packet(packet: Packet,x):
    z=np.asarray(x)-packet.x0
    return (2*np.pi*packet.sigma**2)**(-.25)*np.exp(-z*z/(4*packet.sigma**2)+1j*(packet.p0*z+packet.chirp*z*z/2)/packet.h)


def threshold_overlap(packet: Packet, width=12., epsabs=2e-13):
    """B=<K_0(exp(-x)/h),psi_0>, the continuum low-energy overlap.

    Real/imaginary adaptive quadratures are independent of the finite box.
    Tail beyond 12 Gaussian sigmas is negligible for these specified packets;
    repeat with width=14 and tighter tolerance, not a certified interval bound.
    """
    lo=packet.x0-width*packet.sigma;hi=packet.x0+width*packet.sigma
    def f(x):return k0(np.exp(-x)/packet.h)*continuum_packet(packet,x)
    re,er=quad(lambda x:float(f(x).real),lo,hi,epsabs=epsabs,epsrel=2e-12,limit=300)
    im,ei=quad(lambda x:float(f(x).imag),lo,hi,epsabs=epsabs,epsrel=2e-12,limit=300)
    return complex(re,im),float(er+ei)


def tail_density(x,t,h,B):
    """Leading continuum large-x asymptotic, NOT a fitted tail.

    KL modes phi_k=sqrt(2*k*sinh(pi*k))/pi*K_{ik}(exp(-x)/h).
    a(k)=sqrt(2/pi)*B*k+O(k^3), phi_k~sqrt(2/pi)*sin(k*X+O(k^3)).
    X=x+log(2h)-EulerGamma. Then psi~4i*B*T/(pi*X^3).
    Assumes smooth rapidly decaying spectral overlap and fixed finite T.
    It is NOT used to replace any numerical density or moment sum.
    """
    X=np.asarray(x,float)+np.log(2*h)-np.euler_gamma
    if np.any(X<=0) or h<=0:raise ValueError('Asymptotic coordinate must be positive.')
    return 16*t*t*abs(B)**2/(np.pi**2*X**6)


def classical_log_moments(packet: Packet,times,q=6.,order=96):
    """Same positive initial Gaussian Wigner ensemble; NOT quantum evolution.

    |dx/dT|<=1 implies Jq(0)<=Jq(T)<=exp(2qT)Jq(0), T>=0.
    Thus all these classical moments stay finite at finite T.
    """
    z,w=hermgauss(order);z=np.sqrt(2)*z;w=w/np.sqrt(np.pi)
    x=packet.x0+packet.sigma*z[:,None]+np.zeros((1,order))
    p=packet.p0+packet.chirp*packet.sigma*z[:,None]+packet.h/(2*packet.sigma)*z[None,:]
    ts=np.asarray(times,float)
    evolved,_=classical_path(x.ravel()[:,None],p.ravel()[:,None],ts[None,:])
    lw=np.log((w[:,None]*w[None,:]).ravel())
    return logsumexp(lw[:,None]+q*(evolved+ts[None,:]-packet.x0),axis=0)
