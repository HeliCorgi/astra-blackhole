"""Vacuum KS minisuperspace, one deparametrized Wheeler--DeWitt branch.
Copyright 2026 HeliCorgi. Apache-2.0.

C=p_T^2-p_x^2-exp(-2x)=0; i*h*d_T psi=sqrt(-h^2*d_x^2+exp(-2x))*psi.
T=log(a) is an intrinsic clock, NOT proper time. x=sqrt(3)*Omega-log(4)
and r=(1/4)*exp(-x-T) in the fixed Misner normalization. This is not a
controlled approximation to full quantum gravity, a curvature operator,
or an extension of the classical matter/RTA model.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from numpy.polynomial.hermite import hermgauss
from scipy.linalg import eig_banded


@dataclass(frozen=True)
class Packet:
    h: float
    sigma: float
    chirp: float = 0.
    x0: float = 2.
    p0: float = -1.

    def __post_init__(self):
        if not all(np.isfinite(v) for v in (self.h,self.sigma,self.chirp,self.x0,self.p0)):
            raise ValueError('Finite packet parameters required.')
        if self.h <= 0 or self.sigma <= 0:
            raise ValueError('Positive h and sigma required.')

    def initial(self,x):
        z=np.asarray(x)-self.x0
        v=np.exp(-z*z/(4*self.sigma**2)+1j*(self.p0*z+self.chirp*z*z/2)/self.h)
        norm=np.linalg.norm(v)
        if norm < 1e-100: raise ValueError('Packet lies outside the numerical box.')
        return v/norm # Euclidean mass amplitudes, equivalent to sqrt(dx)*psi(x).

    def continuum_moments(self):
        return dict(mean_x=self.x0,mean_p=self.p0,var_x=self.sigma**2,
                    var_p=self.h**2/(4*self.sigma**2)+self.chirp**2*self.sigma**2,
                    covariance_xp=self.chirp*self.sigma**2)


def classical_path(x0,p0,t):
    """Exact Hamilton trajectories of H=sqrt(p^2+exp(-2x)); vectorized."""
    x0,p0,t=np.asarray(x0),np.asarray(p0),np.asarray(t)
    energy=np.sqrt(p0*p0+np.exp(-2*x0));u0=np.arcsinh(p0*np.exp(x0));u=u0+t
    return np.logaddexp(u,-u)-np.log(2.)-np.log(energy),energy*np.tanh(u)


def classical_ensemble(packet,times,order=96):
    """Liouville control with the EXACT positive Gaussian initial Wigner density.

    Classical Hamilton flow of its principal symbol, not a quantum solver.
    This control separates a point-trajectory approximation from finite spread.
    """
    if order<4:raise ValueError('Quadrature order >=4 required.')
    z,w=hermgauss(order);z=np.sqrt(2.)*z;w=w/np.sqrt(np.pi)
    xx=packet.x0+packet.sigma*z[:,None]+np.zeros((1,order))
    pp=packet.p0+packet.chirp*packet.sigma*z[:,None]+packet.h/(2*packet.sigma)*z[None,:]
    weight=(w[:,None]*w[None,:]).ravel()
    x,p=classical_path(xx.ravel()[:,None],pp.ravel()[:,None],np.asarray(times)[None,:])
    mean=weight@x;variance=weight@(x*x)-mean*mean
    return dict(mean_x=mean,var_x=variance,mean_log_radius_ratio=packet.x0-mean-times)


class WDWModel:
    """Symmetric positive finite-box discretization, spectral time evolution.

    L2=-d_x^2 is the second-order Dirichlet matrix. L4=L2+dx^2 L2^2/12
    is positive and fourth-order in the interior. H=sqrt(h^2 L4+V).
    Dirichlet walls are NUMERICAL regulators; compare expanded domains.
    """
    def __init__(self,h=.2,left=-4.,right=16.,dx=.01):
        if not all(np.isfinite(v) for v in (h,left,right,dx)) or h<=0 or dx<=0 or right<=left:
            raise ValueError('Invalid grid or h.')
        intervals=int(round((right-left)/dx))
        if intervals<20 or intervals>6500 or not np.isclose(intervals*dx,right-left,atol=1e-10):
            raise ValueError('Require 20..6500 aligned intervals.')
        self.h,self.left,self.right,self.dx=h,left,right,dx
        self.x=left+dx*np.arange(1,intervals);n=len(self.x)
        self.band=np.zeros((3,n));q=h*h/(12*dx*dx)
        self.band[0]=30*q+np.exp(-2*self.x);self.band[0,[0,-1]]-=q
        self.band[1,:-1]=-16*q;self.band[2,:-2]=q
        ev,U=eig_banded(self.band,lower=True,check_finite=True)
        if ev[0]<=0:raise RuntimeError('Nonpositive operator: no clipping of eigenvalues.')
        self.eigenvalues,self.energy,self.U=ev,np.sqrt(ev),U

    def operator(self,v):
        one=v.ndim==1;w=v[:,None] if one else v
        ans=self.band[0,:,None]*w
        for j in (1,2):
            b=self.band[j,:-j,None]
            ans[j:]+=b*w[:-j];ans[:-j]+=b*w[j:]
        return ans[:,0] if one else ans

    def momentum(self,v):
        z=np.pad(v,(2,2));der=(-z[4:]+8*z[3:-1]-8*z[1:-3]+z[:-4])/(12*self.dx)
        return -1j*self.h*der

    def evolve(self,packet,times,keep_density=False):
        if not np.isclose(packet.h,self.h):raise ValueError('Packet and operator h must match.')
        times=np.asarray(times,dtype=float)
        if times.ndim!=1 or not len(times) or not np.all(np.isfinite(times)):
            raise ValueError('Finite 1-D time array required.')
        psi=packet.initial(self.x);coef=self.U.T@psi
        amplitudes=coef[:,None]*np.exp(-1j*self.energy[:,None]*times[None,:]/self.h)
        states=self.U@amplitudes;prob=abs(states)**2
        norm=prob.sum(0);mean=self.x@prob;variance=(self.x**2)@prob-mean**2
        energy_state=self.U@(self.energy[:,None]*amplitudes)
        energy_mean=np.sum(states.conj()*energy_state,axis=0).real
        square=self.U@(self.eigenvalues[:,None]*amplitudes)
        residual=np.linalg.norm(self.operator(states)-square,axis=0)/np.linalg.norm(square,axis=0)
        edge=(self.x<self.left+1)|(self.x>self.right-1)
        ppsi=self.momentum(psi);px=np.vdot(psi,ppsi).real
        initial=dict(mean_x=float(self.x@abs(psi)**2),mean_p=float(px),
                     var_x=float((self.x-packet.x0)**2@abs(psi)**2),
                     var_p=float(np.vdot(ppsi,ppsi).real-px*px),
                     covariance_xp=float(np.vdot((self.x-packet.x0)*psi,ppsi).real),
                     mean_H=float(np.sum(abs(coef)**2*self.energy)))
        out=dict(time=times,mean_x=mean,var_x=variance,
                 mean_log_radius_ratio=packet.x0-mean-times,
                 geometric_mean_radius_ratio=np.exp(packet.x0-mean-times),
                 norm_error=float(np.max(abs(norm-1))),
                 energy_relative_drift=float(np.max(abs(energy_mean/initial['mean_H']-1))),
                 discrete_constraint_relative_residual=float(np.max(residual)),
                 max_edge_mass=float(np.max(prob[edge].sum(0))),initial=initial)
        if keep_density:out['density']=prob/self.dx;out['x']=self.x.copy()
        return out
