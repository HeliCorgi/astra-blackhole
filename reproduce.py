"""Reproduce the classical state audits; neural/archive checks are opt-in."""
from pathlib import Path
import argparse, os, subprocess, sys


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--with-legacy',action='store_true',help='Verify restored neural checkpoints; no retraining.')
    parser.add_argument('--out',type=Path,help='Separate generated results from committed baselines.')
    args=parser.parse_args()
    root=Path(__file__).resolve().parent
    env=os.environ.copy()
    env.setdefault('OPENBLAS_NUM_THREADS','1');env.setdefault('OMP_NUM_THREADS','1')
    commands=[[sys.executable,'-m','unittest','discover','-s','tests','-v'],
              [sys.executable,'research/verify_geometry.py']]
    radial=[sys.executable,'research/run_research.py']
    angular=[sys.executable,'research/run_angular_audit.py']
    collision=[sys.executable,'research/run_collision_audit.py']
    if args.out:
        out=args.out.resolve()
        radial+=['--out',str(out/'radial')]
        angular+=['--out',str(out/'angular'),'--check-against',str(root/'research/angular_results/results.json')]
        collision+=['--out',str(out/'collision'),'--check-against',str(root/'research/collision_results/results.json')]
    commands += [radial,angular,collision]
    if args.with_legacy:commands.append([sys.executable,'verify_predictors.py'])
    for command in commands:subprocess.run(command,cwd=root,env=env,check=True)
    print('Requested checks completed. Historical full suites were not rerun.')

if __name__=='__main__':main()
