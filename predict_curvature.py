"""Run a saved forecaster. Example: python predict_curvature.py example_input.json --method neural"""
from __future__ import annotations
import argparse,json
from pathlib import Path
from predictors import polynomial,neural,power_law

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('input',type=Path)
    p.add_argument('--method',choices=['polynomial','adaptive','neural','neural-clean','power'],default='polynomial')
    p.add_argument('--horizon',type=float,default=.1)
    p.add_argument('--output',type=Path)
    a=p.parse_args();data=json.loads(a.input.read_text())
    t=data['time'];k=data['curvature']
    if a.method in ['polynomial','adaptive']:
        out=polynomial(t,k,None if a.method=='adaptive' else a.horizon)
    elif a.method.startswith('neural'):
        out=neural(t,k,data['dlogk_dt'],a.method=='neural')
    else:
        out=power_law(k[-1],data['dlogk_dt']*k[-1],a.horizon)
    text=json.dumps(out,ensure_ascii=False,indent=2)
    if a.output:
        a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(text)
    print(text)
if __name__=='__main__':main()
