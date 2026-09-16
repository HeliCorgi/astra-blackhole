"""Reproduce the flat-vacuum sector of Chiba et al., arXiv:2510.13175v2.

Eq. (14)/(2.14), characteristic data (27)/(3.1), appendix A diamond stencil,
and signed Klein--Gordon current (30)/(4.3). No probability interpretation
of a signed current, new quantization, singularity claim or fitted parameters.
The paper's corner data are incompatible by exp(-8); it is recorded rather
than silently changing the Gaussian. Default corner is set to zero on U=UB.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.signal import lfilter
from scipy.special import j0, j1
from scipy.interpolate import RegularGridInterpolator
from scipy.integrate import simpson

@dataclass(frozen=True)
class Data:
    kappa: float = .1
    sigma: float = .05
    rh: float = .3
    ub: float = -.5

    def __post_init__(self):
        if not all(np.isfinite(x) for x in (self.kappa,self.sigma,self.rh,self.ub)):
            raise ValueError('Finite parameters required.')
        if min(self.kappa,self.sigma,self.rh)<=0 or self.ub>=-self.rh:
            raise ValueError('Positive scales and UB<-rh required.')

    def boundary(self,u):
        z=np.asarray(u)+self.rh
        return np.exp(-z*z/(2*self.sigma**2)-1j*z/self.kappa)

    def derivative(self,u):
        return self.boundary(u)*(-(np.asarray(u)+self.rh)/self.sigma**2-1j/self.kappa)


def green(data: Data,u,v,order=192):
    """Independent continuous Volterra/Bessel solution and analytic derivatives.

    psi=g(U)-V/kappa^2 int_UB^U [2 J1(z)/z]g(s)ds,
    z=2 sqrt((U-s)V)/kappa. The finite-UB right trace is g(UB),
    unlike the paper's zero side boundary. Compare only away from this
    tiny corner mismatch and repeat with a more distant UB.
    Quadrature refines a fixed boundary, not a new state.
    """
    u,v=np.broadcast_arrays(np.asarray(u,float),np.asarray(v,float))
    if np.any(u<data.ub) or np.any(v<0) or not np.all(np.isfinite(u+v)) or order<8:
        raise ValueError('Use U>=UB,V>=0, finite points and order>=8.')
    shape=u.shape; uf=u.ravel();vf=v.ravel()
    nodes,weights=leggauss(order)
    s=data.ub+(uf[:,None]-data.ub)*(nodes+1)/2
    w=weights*(uf[:,None]-data.ub)/2
    z=2*np.sqrt((uf[:,None]-s)*vf[:,None])/data.kappa
    small=abs(z)<1e-5
    k1=np.empty_like(z)
    k1[small]=1-z[small]**2/8+z[small]**4/192
    k1[~small]=2*j1(z[~small])/z[~small]
    f=data.boundary(s);a=1/data.kappa**2;g=data.boundary(uf)
    psi=g-a*vf*np.sum(w*f*k1,axis=1)
    zb=2*np.sqrt((uf-data.ub)*vf)/data.kappa
    kb=np.ones_like(zb);active=abs(zb)>1e-5
    kb[active]=2*j1(zb[active])/zb[active]
    kb[~active]=1-zb[~active]**2/8+zb[~active]**4/192
    # Integrate by parts to avoid subtracting O((aV)^2) terms.
    du=data.derivative(uf)-a*vf*(data.boundary(data.ub)*kb+
                                   np.sum(w*data.derivative(s)*k1,axis=1))
    dv=-a*np.sum(w*f*j0(z),axis=1)
    return tuple(x.reshape(shape) for x in (psi,du,dv))


def currents(psi,du,dv):
    """qT=i psi* <-> dT psi, jX=-i psi* <-> dX psi; generally signed."""
    qu=-2*np.imag(np.conj(psi)*du);qv=-2*np.imag(np.conj(psi)*dv)
    return qu+qv,qu-qv


class Diamond:
    """Paper's phase-factored, second-order null stencil, Lambda=0,k=1."""
    def __init__(self,data: Data,step=.001,umax=1.3,vmax=1.8):
        if not all(np.isfinite(x) for x in (step,umax,vmax)) or step<=0 or umax<=data.ub or vmax<=0:
            raise ValueError('Invalid domain.')
        nu=int(round((umax-data.ub)/step));nv=int(round(vmax/step))
        if not (np.isclose(nu*step,umax-data.ub) and np.isclose(nv*step,vmax)) or min(nu,nv)<4:
            raise ValueError('Aligned grid with at least four intervals required.')
        if (nu+1)*(nv+1)>65_000_000:raise ValueError('Audit memory bound exceeded.')
        self.data,self.step=data,step
        self.u=data.ub+step*np.arange(nu+1);self.v=step*np.arange(nv+1)
        f=np.zeros((nv+1,nu+1),complex)
        f[0]=np.exp(-(self.u+data.rh)**2/(2*data.sigma**2));f[0,0]=0
        a=1/(1-1j*step/data.kappa);b=(1+1j*step/data.kappa)*a
        for j in range(1,nv+1):
            f[j,1:]=lfilter([1.],[1.,-a],a*f[j-1,1:]-b*f[j-1,:-1])
        self.envelope=f
        # Interpolate phase-factored derivatives, not an under-resolved rapid phase.
        self.interp=RegularGridInterpolator((self.v,self.u),f,bounds_error=True)

    def evaluate(self,u,v):
        u,v=np.broadcast_arrays(np.asarray(u,float),np.asarray(v,float));d=self.step
        if np.any(u<self.u[0]+2*d) or np.any(u>self.u[-1]-2*d) or np.any(v<2*d) or np.any(v>self.v[-1]-2*d):
            raise ValueError('Derivative evaluation needs a two-cell interior margin.')
        def f(uu,vv):return self.interp(np.c_[np.ravel(vv),np.ravel(uu)]).reshape(u.shape)
        f0=f(u,v)
        fu=(-f(u+2*d,v)+8*f(u+d,v)-8*f(u-d,v)+f(u-2*d,v))/(12*d)
        fv=(-f(u,v+2*d)+8*f(u,v+d)-8*f(u,v-d)+f(u,v-2*d))/(12*d)
        phase=np.exp(-1j*(u+v+self.data.rh)/self.data.kappa)
        return phase*f0,phase*(fu-1j*f0/self.data.kappa),phase*(fv-1j*f0/self.data.kappa)

    def slice(self,t,step=None):
        d=self.step;spacing=d if step is None else step
        lo=max(self.u[0]+3*d,2*t-self.v[-1]+3*d)
        hi=min(self.u[-1]-3*d,2*t-3*d)
        if hi<=lo:raise ValueError('Clock slice outside computed domain.')
        n=int(np.ceil((hi-lo)/spacing));u=np.linspace(lo,hi,n+1);v=2*t-u
        return slice_observables(t,u,*self.evaluate(u,v))


def slice_observables(t,u,psi,du,dv):
    x=t-u;q,j=currents(psi,du,dv);norm=simpson(q,x=u)
    if abs(norm)<1e-12:raise ValueError('Signed norm too close to zero.')
    negative=simpson(np.maximum(-q,0),x=u);absolute=simpson(abs(q),x=u)
    mean=simpson(x*q,x=u)/norm
    mask=u<=0
    interior=simpson(q[mask],x=u[mask]) if np.count_nonzero(mask)>2 else 0.
    return dict(t=float(t),norm=float(norm),mean_X=float(mean),mean_U=float(t-mean),
                mean_V=float(t+mean),negative_charge=float(negative),
                negative_over_absolute=float(negative/absolute),interior_charge=float(interior),
                signed_exterior_fraction=float(1-interior/norm))


def exact_slice(data,t,order=192,samples=2001):
    u=np.linspace(data.ub,2*t,samples)
    return slice_observables(t,u,*green(data,u,2*t-u,order))


def pullback(u,v):
    u,v=np.broadcast_arrays(np.asarray(u,float),np.asarray(v,float))
    if np.any(u>=0) or np.any(v<=0):raise ValueError('Internal chart requires U<0,V>0.')
    return -.5*np.log(-16*u*v),.5*np.log(-v/u)


def same_slice_flux(data,t,upper=-1e-4,lower=-.49,order=384,green_order=256):
    """Same finite interior segment, two independent quadrature coordinates.

    Constant Cartesian time t is NOT a constant boost-time tau surface.
    Its current is q_tau-j_x*d_tau/dx, not |psi|^2.
    """
    nodes,weights=leggauss(order)
    u=lower+(upper-lower)*(nodes+1)/2;v=2*t-u
    psi,du,dv=green(data,u,v,green_order);q,_=currents(psi,du,dv)
    uv=weights@q*(upper-lower)/2
    uv_numerator=weights@((t-u)*q)*(upper-lower)/2
    xlo,_=pullback(upper,2*t-upper);xhi,_=pullback(lower,2*t-lower)
    # x decreases with -U on a t>=0 slice.
    low,high=min(xlo,xhi),max(xlo,xhi)
    x=low+(high-low)*(nodes+1)/2
    X=np.sqrt(t*t+np.exp(-2*x)/16);u=t-X;v=t+X
    psi,du,dv=green(data,u,v,green_order)
    dtau=-u*du+v*dv;dx=-u*du-v*dv
    qt=-2*np.imag(np.conj(psi)*dtau);jx=2*np.imag(np.conj(psi)*dx)
    transformed=weights@(qt-jx*t/X)*(high-low)/2
    chart_numerator=weights@(X*(qt-jx*t/X))*(high-low)/2
    return dict(t=t,UV_flux=float(uv),chart_flux=float(transformed),
                relative_difference=float(abs(uv-transformed)/abs(uv)),
                UV_mean_X=float(uv_numerator/uv),
                chart_mean_X=float(chart_numerator/transformed),
                mean_X_absolute_difference=float(abs(uv_numerator/uv-chart_numerator/transformed)))
