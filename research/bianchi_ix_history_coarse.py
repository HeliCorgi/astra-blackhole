"""Derived coarse-graining audit for the validated Bianchi IX history matrix."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
LABELS=('A','B+','B-')
PARTITIONS={'fine_3':((0,),(1,),(2,)),'A_vs_Bs':((0,),(1,2)),'Bplus_vs_rest':((1,),(0,2)),'Bminus_vs_rest':((2,),(0,1)),'ignore_second':((0,1,2),)}
def blocks(data):
 out={}
 for k in LABELS:
  a=np.array([[complex(z['re'],z['im']) for z in row] for row in data['third_sector_blocks'][k]],complex)
  if a.shape!=(3,3) or np.max(abs(a-a.conj().T))>1e-10:raise ValueError('invalid Hermitian 3x3 block')
  out[k]=a
 return out
def coarse_block(B,groups):
 C=np.zeros((len(groups),len(groups)),complex)
 for a,g in enumerate(groups):
  for b,h in enumerate(groups):C[a,b]=sum(B[i,j] for i in g for j in h)
 return C
def metrics(name,groups,bs,delta,regrid):
 Cs={k:coarse_block(bs[k],groups) for k in LABELS};diag=np.concatenate([np.diag(Cs[k]).real for k in LABELS])
 offvals=[];realvals=[];coh=[];rcoh=[];defects={}
 for k in LABELS:
  C=Cs[k];d=np.diag(C).real;defects[k]=float(C.sum().real-d.sum())
  for i in range(len(groups)):
   for j in range(i):
    z=C[i,j];offvals.append(abs(z));realvals.append(abs(z.real));den=np.sqrt(max(d[i],0)*max(d[j],0))
    if den>1e-15:coh.append(abs(z)/den);rcoh.append(abs(z.real)/den)
 maxoff=float(max(offvals,default=0));maxreal=float(max(realvals,default=0));maxprod=max((len(g)*len(h) for i,g in enumerate(groups) for j,h in enumerate(groups) if i!=j),default=0);bound=float(maxprod*delta)
 return {'name':name,'groups':[[LABELS[i] for i in g] for g in groups],'family_size':3*len(groups),'trace_diagonal':float(diag.sum()),'max_offdiagonal_abs':maxoff,'max_offdiagonal_real_abs':maxreal,'max_normalized_pair_coherence':float(max(coh,default=0)),'max_normalized_real_pair':float(max(rcoh,default=0)),'final_additivity_defects':defects,'max_final_additivity_defect_abs':float(max(map(abs,defects.values()))),'propagated_main_control_offdiag_bound':bound,'strict_signal_to_propagated_control_bound':None if bound==0 else float(maxoff/bound),'max_real_to_recorded_regrid_gram_change':float(maxreal/regrid) if regrid else None}
def audit(data):
 bs=blocks(data);delta=float(data['fine_max_offdiag_main_control_difference']);regrid=float(data['recorded_second_branch_gram_relative_change']);rows=[metrics(n,g,bs,delta,regrid) for n,g in PARTITIONS.items()];by={x['name']:x for x in rows}
 return {'schema':1,'method':'All five set partitions of A/B+/B- at the second time, using exact block sums of the validated fine decoherence matrix.','partitions':rows,'headline':{'fine_max_offdiagonal_abs':by['fine_3']['max_offdiagonal_abs'],'A_vs_Bs_max_offdiagonal_abs':by['A_vs_Bs']['max_offdiagonal_abs'],'A_vs_Bs_normalized_pair_coherence':by['A_vs_Bs']['max_normalized_pair_coherence'],'A_vs_Bs_signal_to_control_bound':by['A_vs_Bs']['strict_signal_to_propagated_control_bound'],'Bplus_vs_rest_signal_to_control_bound':by['Bplus_vs_rest']['strict_signal_to_propagated_control_bound'],'Bminus_vs_rest_signal_to_control_bound':by['Bminus_vs_rest']['strict_signal_to_propagated_control_bound']},'interpretation':['A versus {B+,B-} reduces strict complex interference, but the residual is only about 1.26 times the conservatively propagated 102-vs-114 control difference; strict decoherence is unresolved for this binary family.','B+ versus rest and B- versus rest retain strict complex interference more than five times the propagated two-resolution control scale; strict decoherence is not accepted for those families.','Ignoring the second-time alternative leaves only orthogonal final-time projectors and is diagonal by construction; this is algebraic coarse graining, not evidence of dynamical decoherence.','Weak/real consistency remains limited by the recorded cubic-regrid systematic.'],'limitations':['Derived from one validated finite-time D matrix; no new time evolution.','The propagated control scale compares two block-Krylov dimensions and is not a rigorous continuum error bar.','Clock, wall sectors, first-A conditioning, finite boxes and regrids are inherited from the parent audit.','No environmental decoherence, black-hole observation, singularity resolution, or unique quantum-gravity claim.']}
def signature(d):
 v=[]
 for x in d['partitions']:
  v += [x['trace_diagonal'],x['max_offdiagonal_abs'],x['max_offdiagonal_real_abs'],x['max_normalized_pair_coherence'],x['max_final_additivity_defect_abs']]
  if x['strict_signal_to_propagated_control_bound'] is not None:v.append(x['strict_signal_to_propagated_control_bound'])
 return np.asarray(v,float)
def main(argv=None):
 ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--check-against',type=Path);a=ap.parse_args(argv);r=audit(json.loads(a.input.read_text()));a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
 if a.check_against:
  old=json.loads(a.check_against.read_text());x=signature(r);y=signature(old)
  if x.shape!=y.shape or not np.allclose(x,y,rtol=2e-12,atol=2e-14):raise RuntimeError('coarse audit regression mismatch')
 print(json.dumps(r['headline'],indent=2))
if __name__=='__main__':main()
