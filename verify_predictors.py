"""Small repeatability check; does not rerun all archived studies or retrain models."""
from __future__ import annotations
from pathlib import Path
import json,numpy as np
from predictors import ROOT,load_module,polynomial,neural

def main():
    out=ROOT/'verification';out.mkdir(exist_ok=True)
    m=load_module('_old_gen_check',ROOT/'legacy/latent_state_audit/run_audit.py')
    pair=m.pair_data(1,999,k_bounds=(1.,1.),w_bounds=(3.6,3.6))
    data={'time':m.PAST.tolist(),'curvature':np.exp(pair['hist'][0]).tolist(),
          'dlogk_dt':float(pair['x'][0,1]),'note':'Synthetic LTB trajectory A; model units only.'}
    (ROOT/'example_input.json').write_text(json.dumps(data,indent=2))
    output={'scope':'Subset smoke/inference checks; historical full suites not rerun.',
            'clean_predictions':[],'robust_predictions':[]}
    ref=json.loads((ROOT/'legacy/latent_state_audit/results.json').read_text())
    for i in (0,1):
        k=np.exp(pair['hist'][i]);d=pair['x'][i,1]
        p=neural(m.PAST,k,d,False);q=neural(m.PAST,k,d,True)
        output['clean_predictions'].append(p);output['robust_predictions'].append(q)
        np.testing.assert_allclose(p['predicted_K'],ref['example']['latent_predicted_K'][i],rtol=3e-6,atol=1e-6)
    output['fixed_horizon']=polynomial(data['time'],data['curvature'],.1)
    output['adaptive']=polynomial(data['time'],data['curvature'])
    sch=load_module('_old_sch',ROOT/'legacy/singularity_forecast_check.py')
    output['schwarzschild_recomputed']=[sch.benchmark(x) for x in (.1,.01,.001)]
    output['pass']=True
    (out/'predictor_smoke_test.json').write_text(json.dumps(output,indent=2))
    print('Predictor smoke tests passed; both checkpoints and polynomial variants execute.')
if __name__=='__main__':main()
