"""Audit whether the current reduced model can exhibit deterministic chaos."""
from __future__ import annotations
import argparse, json, platform, hashlib
from pathlib import Path
import numpy as np
from trajectory_chaos_model import (
    ks_energy, ks_path, finite_time_lyapunov_ks,
    gauss_orbit, gauss_finite_lyapunov, gauss_ensemble_lyapunov,
    gauss_exact_lyapunov,
)

ROOT = Path(__file__).resolve().parent


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def run(out: Path):
    out.mkdir(parents=True, exist_ok=True)
    times = np.geomspace(1.0, 1e4, 81)
    perturbations = [(1e-8,0.0),(0.0,1e-8),(1e-8,1e-8)]
    ks_rows=[]
    for pert in perturbations:
        lam, sep = finite_time_lyapunov_ks(2.0,-1.0,pert,times)
        ks_rows.append({
            'perturbation': list(pert),
            'lambda_at_10': float(lam[np.argmin(abs(times-10))]),
            'lambda_at_100': float(lam[np.argmin(abs(times-100))]),
            'lambda_at_1000': float(lam[np.argmin(abs(times-1000))]),
            'lambda_at_10000': float(lam[-1]),
            'max_abs_energy_relative_drift': float(max(
                np.max(abs(ks_energy(*ks_path(2.0,-1.0,times))/ks_energy(2.0,-1.0)-1)),
                np.max(abs(ks_energy(*ks_path(2.0+pert[0],-1.0+pert[1],times))/ks_energy(2.0+pert[0],-1.0+pert[1])-1)),
            )),
            'separation_at_10000': float(sep[-1]),
        })
    x0=np.pi-3.0
    gauss_steps=[10,100,1000,10000]
    gauss_rows=[{'steps':n,'finite_lyapunov':gauss_finite_lyapunov(x0,n)} for n in gauss_steps]
    exact=gauss_exact_lyapunov(); ensemble=gauss_ensemble_lyapunov()
    a=gauss_orbit(x0,40); b=gauss_orbit(x0+1e-12,40); same=gauss_orbit(x0,40)
    first_diff=int(np.argmax(abs(a-b)>1e-6)) if np.any(abs(a-b)>1e-6) else None
    result={
      'schema':1,
      'scope':'Model-internal determinism/chaos audit. Not evidence that a physical Schwarzschild interior is random.',
      'reduced_ks':{
        'hamiltonian':'H=sqrt(p^2+exp(-2x))',
        'degrees_of_freedom':1,
        'analytic_solution_used':True,
        'finite_time_lyapunov':ks_rows,
        'interpretation':'The tested asymptotic finite-time exponents tend toward zero; this reduced autonomous one-degree-of-freedom model is integrable, so it cannot test Mixmaster chaos.'
      },
      'bkl_gauss_control':{
        'map':'x -> frac(1/x)',
        'initial':x0,
        'finite_time':gauss_rows,
        'invariant_measure_ensemble':ensemble,
        'exact_invariant_measure_value':exact,
        'ensemble_relative_error':abs(ensemble/exact-1),
        'nearby_orbit_first_difference_gt_1e-6':first_diff,
        'identical_input_bitwise_repeat':bool(np.array_equal(a,same)),
        'interpretation':'Positive Lyapunov exponent is deterministic chaos, not fundamental randomness.'
      },
      'quantum_interpretation':{
        'statement':'The WDW wavefunction evolves linearly/deterministically after clock, branch, inner product and initial state are chosen. A single classical trajectory is not part of that state description; probabilities or signed KG currents are different objects. This is distinct from deterministic chaos.'
      },
      'next_model_requirement':'To test the user hypothesis inside gravity rather than as a control, relax the KS symmetry to at least two anisotropy degrees of freedom (e.g. Bianchi IX / cosmological billiard) or an inhomogeneous black-hole interior known to support BKL-like dynamics.',
      'environment':{'python':platform.python_version(),'numpy':np.__version__},
      'sources':{name:sha(ROOT/name) for name in ['trajectory_chaos_model.py','run_trajectory_chaos_audit.py']}
    }
    (out/'summary.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    curves={'ks_times':times.tolist(),'ks':[], 'gauss_first_40':a.tolist(),'gauss_nearby_first_40':b.tolist()}
    for pert in perturbations:
        lam,sep=finite_time_lyapunov_ks(2,-1,pert,times)
        curves['ks'].append({'perturbation':list(pert),'lyapunov':lam.tolist(),'separation':sep.tolist()})
    (out/'curves.json').write_text(json.dumps(curves,separators=(',',':')))
    return result


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,default=ROOT/'trajectory_chaos_results');args=ap.parse_args()
    r=run(args.out);print(json.dumps(r,indent=2,ensure_ascii=False))
if __name__=='__main__':main()
