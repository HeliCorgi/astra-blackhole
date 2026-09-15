"""Reproduce the published classical audit; neural/archive checks are opt-in."""
from pathlib import Path
import argparse, os, subprocess, sys

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--with-legacy',action='store_true',help='Also verify restored neural checkpoints; does not retrain.')
    args=parser.parse_args()
    root=Path(__file__).resolve().parent
    env=os.environ.copy()
    env.setdefault('OPENBLAS_NUM_THREADS','1');env.setdefault('OMP_NUM_THREADS','1')
    commands=[[sys.executable,'-m','unittest','discover','-s','tests','-v']]
    commands += [[sys.executable,s] for s in ['research/verify_geometry.py','research/run_research.py']]
    if args.with_legacy: commands.append([sys.executable,'verify_predictors.py'])
    for command in commands: subprocess.run(command,cwd=root,env=env,check=True)
    print('Requested checks completed. Historical full suites were not rerun.')
if __name__=='__main__':main()
