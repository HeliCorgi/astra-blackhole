"""Convenience interface to the PREVIOUS curvature forecasters.

These are model-benchmark predictors, not validated black-hole observations.
The fixed-horizon source is under legacy/. Optional neural checkpoints and
historical studies can be restored from the original conversation archive.
"""
from __future__ import annotations
from pathlib import Path
import importlib.util
import json
import sys
import numpy as np
ROOT=Path(__file__).resolve().parent


def load_module(name: str, path: Path):
    spec=importlib.util.spec_from_file_location(name,path)
    if spec is None or spec.loader is None: raise ImportError(str(path))
    mod=importlib.util.module_from_spec(spec)
    sys.modules[name]=mod
    spec.loader.exec_module(mod)
    return mod


def checked_history(time, curvature):
    t=np.asarray(time,dtype=float);k=np.asarray(curvature,dtype=float)
    if t.ndim!=1 or k.shape!=t.shape or len(t)<8:
        raise ValueError('Need at least 8 matching, one-dimensional samples.')
    if not np.all(np.isfinite(t)) or not np.all(np.isfinite(k)):
        raise ValueError('Samples must be finite.')
    if np.any(np.diff(t)<=0) or np.any(k<=0):
        raise ValueError('Strictly increasing times and positive K required.')
    return t,k


def polynomial(time, curvature, horizon: float | None=None, tolerance: float=.01):
    t,k=checked_history(time,curvature)
    if tolerance<=0 or not np.isfinite(tolerance):raise ValueError('Positive tolerance required.')
    if horizon is None:
        m=load_module('_old_adaptive',ROOT/'adaptive_forecast.py')
        p=m.predict(t,k,tolerance)
        return dict(method='legacy_adaptive_quadratic',horizon=p['adaptive_h'],
                    predicted_K=float(k[-1]*p['adaptive_ratio']),
                    nominal_horizon=p['nominal_h'],halvings=p['halvings'],
                    disagreement=p['disagreement'],
                    warning='Disagreement tolerance is NOT a true error bound.')
    if not np.isfinite(horizon) or horizon<=0:raise ValueError('Positive finite horizon required.')
    m=load_module('_old_fixed',ROOT/'legacy/collapse_extensions/forecast.py')
    p=m.forecast(t,k,horizon)
    return dict(method='legacy_fixed_horizon_polynomials',horizon=horizon,
                predicted_K_by_degree={str(a):b for a,b in p.items()},
                warning='Null means the transformed polynomial is nonpositive. No error guarantee.')


def power_law(k0: float, dk_dt: float, horizon: float):
    if not all(np.isfinite(v) and v>0 for v in (k0,dk_dt,horizon)):
        raise ValueError('K0, dK/dt and horizon must be positive finite values.')
    remaining=4*k0/dk_dt
    if horizon>=remaining:raise ValueError('Forecast reaches the assumed pole.')
    return dict(method='assumed_inverse_fourth_power',horizon=horizon,
                predicted_K=k0/(1-horizon/remaining)**4,remaining_estimate=remaining)


def neural(time, curvature, dlogk_dt: float, robust: bool=True):
    """Frozen network; supports ONLY the original offset grid and model units.

    Does not receive future data. A larger network is NOT a different physics
    model. This interface intentionally refuses other sampling grids.
    """
    t,k=checked_history(time,curvature)
    if len(t)!=41 or not np.allclose(t-t[-1],np.linspace(-.30,0.,41),atol=1e-12,rtol=0):
        raise ValueError('Network requires 41 offsets evenly spaced from -0.30 to 0.00 in original LTB model units.')
    if not np.isfinite(dlogk_dt):raise ValueError('Finite d(log K)/dt required.')
    folder=ROOT/'legacy/latent_state_audit'
    prefix='robust_' if robust else ''
    needed=[folder/'run_audit.py',folder/(prefix+'encoder_predictor.pt'),folder/(prefix+'normalization.json')]
    if not all(path.is_file() for path in needed):
        raise FileNotFoundError('Optional neural artifacts are absent. Run tools/restore_legacy.py with the original ZIP, then install requirements-neural.txt.')
    import torch
    m=load_module('_old_latent',folder/'run_audit.py')
    saved=torch.load(folder/(prefix+'encoder_predictor.pt'),map_location='cpu',weights_only=True)
    network=m.PredictiveEncoder(1)
    network.load_state_dict(saved['state_dict']);network.eval()
    norm=json.loads((folder/(prefix+'normalization.json')).read_text())
    norm={a:{b:np.asarray(v,dtype=float) for b,v in d.items()} for a,d in norm.items()}
    current=np.array([[np.log(k[-1]),dlogk_dt]])
    history=(np.log(k/k[-1])-dlogk_dt*(t-t[-1]))[None,:-1]
    x=torch.tensor((current-norm['x']['mean'])/norm['x']['std'],dtype=torch.float32)
    h=torch.tensor((history-norm['h']['mean'])/norm['h']['std'],dtype=torch.float32)
    with torch.no_grad():
        z=network.encode(h); yy=network.decode(x,z).numpy()
    logratio=yy*norm['y']['std']+norm['y']['mean']
    return dict(method=prefix+'legacy_neural_latent1',future_offsets=[.1,.2,.3,.4],
                predicted_K=(k[-1]*np.exp(logratio[0])).tolist(),latent_z=float(z[0,0]),
                warning='Trained only on the specified synthetic LTB family. The derivative input was exact in prior tests; no real-world validity or calibrated uncertainty.')
