"""Past-only local curvature forecasts. No model parameters or future samples."""
from __future__ import annotations
import numpy as np

def forecast(time, curvature, horizon):
    t=np.asarray(time,dtype=float); K=np.asarray(curvature,dtype=float)
    if t.ndim!=1 or t.shape!=K.shape or len(t)<8 or np.any(np.diff(t)<=0):
        raise ValueError("Need >=8 strictly time-ordered samples.")
    if np.any(K<=0) or not np.all(np.isfinite(K)) or horizon<=0:
        raise ValueError("Curvature must be positive and horizon positive.")
    width=t[-1]-t[0]; x=(t-t[-1])/width
    q=(K/K[-1])**(-.25)
    out={}
    for degree in (1,2,3):
        A=np.column_stack([x**j for j in range(1,degree+1)])
        c=np.linalg.lstsq(A,q-1,rcond=None)[0]
        xp=horizon/width
        qp=1+sum(c[j-1]*xp**j for j in range(1,degree+1))
        if qp<=0:
            out[degree]=None
        else:
            out[degree]=float(K[-1]*qp**(-4))
    return out
