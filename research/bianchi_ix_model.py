"""Vacuum Bianchi IX / Mixmaster classical bridge audit.

Misner variables with s=-alpha increasing toward the singular regime.
Potential normalization follows the convention
C=1/2 e^-3alpha(-p_alpha^2+p_+^2+p_-^2)+e^alpha V=0,
H=-p_alpha=sqrt(p_+^2+p_-^2+2 e^(4 alpha)V).

This module is a homogeneous classical minisuperspace model, not a full
black-hole interior, stochastic law, or quantum Bianchi IX theory.
"""
from __future__ import annotations
import numpy as np
from scipy.integrate import solve_ivp

SQ3=np.sqrt(3.0)

def _exp(x):
    return np.exp(np.clip(x,-745.0,700.0))

def wall_term_and_gradient(s: float,bp: float,bm: float):
    """Return W=2 exp(-4s)V and derivatives dW/dbeta_+, dW/dbeta_-.

    Algebraically expands cosh terms before exponentiation to postpone
    overflow near steep Mixmaster walls. No clipping of the physical W is
    applied; exponent clipping is only a float-range guard and accepted runs
    stay far below it.
    """
    A=-4*s-8*bp
    Bp=-4*s+4*bp+4*SQ3*bm; Bm=-4*s+4*bp-4*SQ3*bm
    B0=-4*s+4*bp
    C1=-4*s-2*bp+2*SQ3*bm; C2=-4*s-2*bp-2*SQ3*bm
    eA,eBp,eBm,eB0,eC1,eC2=map(_exp,(A,Bp,Bm,B0,C1,C2))
    bracket=eA+eBp+eBm-2*eB0-2*eC1-2*eC2
    W=bracket/3.0
    d_bp=(-8*eA+4*eBp+4*eBm-8*eB0+4*eC1+4*eC2)/3.0
    d_bm=(4*SQ3*(eBp-eBm)-4*SQ3*(eC1-eC2))/3.0
    return float(W),float(d_bp),float(d_bm)

def potential(bp: float,bm: float):
    return float((_exp(-8*bp)+2*_exp(4*bp)*(np.cosh(4*SQ3*bm)-1)
                  -4*_exp(-2*bp)*np.cosh(2*SQ3*bm))/6.0)

def rhs(s,y):
    bp,bm,pp,pm=np.asarray(y,float)
    W,gp,gm=wall_term_and_gradient(float(s),bp,bm)
    rad=pp*pp+pm*pm+W
    if not np.isfinite(rad) or rad<=1e-14:
        raise RuntimeError('Reduced Hamiltonian left positive branch.')
    H=np.sqrt(rad)
    return np.array([-pp/H,-pm/H,0.5*gp/H,0.5*gm/H])

def evolve(y0,s0=2.0,s1=25.0,t_eval=None,rtol=2e-10,atol=2e-12,max_step=.01,method='DOP853'):
    y0=np.asarray(y0,float)
    if y0.shape!=(4,) or not np.all(np.isfinite(y0)) or not (s1>s0):
        raise ValueError('Finite four-state and increasing s interval required.')
    sol=solve_ivp(rhs,(s0,s1),y0,t_eval=t_eval,rtol=rtol,atol=atol,
                  max_step=max_step,method=method,dense_output=t_eval is None)
    if not sol.success: raise RuntimeError(sol.message)
    return sol

def kasner_from_state(s,y):
    bp,bm,pp,pm=np.asarray(y,float)
    W,_,_=wall_term_and_gradient(float(s),bp,bm)
    H=np.sqrt(pp*pp+pm*pm+W)
    vp,vm=pp/H,pm/H # d beta/d alpha; s=-alpha changes coordinate velocity sign only
    p=np.array([(1+vp+SQ3*vm)/3,(1+vp-SQ3*vm)/3,(1-2*vp)/3])
    ordered=np.sort(p)
    u=float(ordered[2]/ordered[1]) if ordered[1]>0 else np.nan
    wall_ratio=float(abs(W)/(pp*pp+pm*pm))
    return {'p':p,'u':u,'wall_ratio':wall_ratio,
            'sum':float(p.sum()),'sum_squares':float(p@p)}

def bkl_step(u: float):
    u=float(u)
    if not np.isfinite(u) or u<=1: raise ValueError('BKL u must be >1.')
    return u-1 if u>=2 else 1/(u-1)

def bkl_orbit(u,n):
    out=[float(u)]
    for _ in range(int(n)): out.append(bkl_step(out[-1]))
    return np.asarray(out)

def first_bkl_divergence(u,delta=1e-12,threshold=.1,max_steps=200):
    a,b=float(u),float(u+delta)
    rows=[]
    for n in range(1,max_steps+1):
        a,b=bkl_step(a),bkl_step(b); gap=abs(a-b); rows.append((n,a,b,gap))
        if gap>=threshold:return rows
    return rows
