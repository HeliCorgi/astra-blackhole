from __future__ import annotations
import argparse,json,platform,time,sys
from pathlib import Path
import numpy as np, scipy
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from bianchi_ix_model import evolve, kasner_from_state, bkl_step
from quantum_bianchi_ix_model import QuantumBianchiIX, PacketSpec
from quantum_bianchi_ix_long_model import positive_wall_terms, extended_observables, propagate_state, regrid_state, fidelity, branch_interference, classical_ensemble, WALL_NAMES

def md(a,b): return float(np.linalg.norm(np.asarray(a,float)-np.asarray(b,float)))
def tv(a,b): return 0.5*float(np.sum(abs(np.asarray(a,float)-np.asarray(b,float))))
def compact_obs(o): return {k:o[k] for k in ('mean','rms_width','norm','edge_mass','capped_region_mass','edge_mass_w1','wall_channels') if k in o}
def qgrid(bp,bm,nx,ny,cap=30): return QuantumBianchiIX(hbar=.2,bp_bounds=bp,bm_bounds=bm,nx=nx,ny=ny,potential_cap=cap)

def regrid_record(v,old,new,s):
    w,n1=regrid_state(v,old,new,3);back,n2=regrid_state(w,new,old,3)
    return w,{'s':s,'from':[old.bp_bounds,old.bm_bounds,old.nx,old.ny],'to':[new.bp_bounds,new.bm_bounds,new.nx,new.ny],'raw_forward_norm':n1,'raw_back_norm':n2,'roundtrip_fidelity':fidelity(v,back),'new_observables':compact_obs(extended_observables(new,w,s))}

def run(out:Path):
    out.mkdir(parents=True,exist_ok=True);t0=time.time();packet=PacketSpec(hbar=.2,sigma=.3)
    g1=qgrid((-5,5),(-8,5),100,130);g2=qgrid((-12,12),(-18,10),144,168);g3=qgrid((-18,18),(-24,24),144,192);g4=qgrid((-32,32),(-36,36),214,240)
    v=g1.initial(packet)
    seg1=propagate_state(g1,v,2,6,.05,48,[2,4.6,6]);v,r1=regrid_record(seg1['state'],g1,g2,6)
    seg2=propagate_state(g2,v,6,12,.1,52,[8,9.5,10.385,11,12]);v,r2=regrid_record(seg2['state'],g2,g3,12)
    seg3=propagate_state(g3,v,12,20,.15,52,[14,16,18,20]);v,r3=regrid_record(seg3['state'],g3,g4,20)
    seg4=propagate_state(g4,v,20,36,.2,52,[24,28,30,32,32.46,34,36])
    qrows={}
    for seg in (seg1,seg2,seg3,seg4):
        for s,row in seg['rows'].items(): qrows[float(s)]=row
    times=sorted(set([2,4.6,6,10.385,11,14,20,24,30,32.46,34,36]))
    ens=classical_ensemble(packet,times,n=4096,seed=20260917,s0=2,ds=.01)
    sol=evolve(np.array([0.,0.,1.,.31]),2,36,t_eval=np.asarray(times),rtol=5e-10,atol=5e-12,max_step=.01)
    central={float(s):{'mean':sol.y[:2,i].tolist(),'state':sol.y[:,i].tolist()} for i,s in enumerate(times)}
    comparisons=[]
    for s in times:
        qo=qrows[s]['obs'];eo=ens[s];co=central[s]
        comparisons.append({'s':s,'quantum_to_single_trajectory':md(qo['mean'],co['mean']),'quantum_to_classical_ensemble_mean':md(qo['mean'],eo['mean']),'quantum_rms_width':qo['rms_width'],'classical_ensemble_rms_width':eo['rms_width'],'quantum_wall_sector_mass':qo['wall_channels'],'classical_ensemble_wall_sector_mass':eo['wall_channels'],'wall_sector_total_variation':tv(qo['wall_channels'],eo['wall_channels'])})
    bounce_info=[]
    for idx,s in [(2,10.385),(3,32.46)]:
        y=np.asarray(central[s]['state']);terms=positive_wall_terms(s,np.array([[y[0]]]),np.array([[y[1]]]))[:,0,0];lab=int(np.argmax(terms));c=next(x for x in comparisons if x['s']==s)
        bounce_info.append({'bounce_number':idx,'s':s,'central_wall':WALL_NAMES[lab],'quantum_sector_mass':c['quantum_wall_sector_mass'],'classical_ensemble_sector_mass':c['classical_ensemble_wall_sector_mass'],'sector_total_variation':c['wall_sector_total_variation']})
    inter2=branch_interference(g2,seg2['rows'][9.5]['state'],9.5,11,.05,60,None)
    inter3=branch_interference(g4,seg4['rows'][30.0]['state'],30,34,.2,60,None)
    gf2=qgrid((-12,12),(-18,10),168,196);vf2,_=regrid_state(seg2['rows'][9.5]['state'],g2,gf2,3);back2,_=regrid_state(vf2,gf2,g2,3);fine2=propagate_state(gf2,vf2,9.5,11,.075,56,[11])['rows'][11]['obs'];main11=seg2['rows'][11]['obs']
    second_refine={'start_regrid_fidelity':fidelity(seg2['rows'][9.5]['state'],back2),'mean_distance':md(main11['mean'],fine2['mean']),'wall_sector_tv':tv(main11['wall_channels'],fine2['wall_channels']),'rms_abs_difference':abs(main11['rms_width']-fine2['rms_width'])}
    v30=seg4['rows'][30.0]['state'];main34=seg4['rows'][34.0]['obs'];controls=[]
    for name,model,ds,k in [('dt01',g4,.1,52),('k60',g4,.2,60),('cap50',qgrid((-32,32),(-36,36),214,240,50),.2,52)]:
        oo=propagate_state(model,v30,30,34,ds,k,[34])['rows'][34]['obs'];controls.append({'case':name,'mean_distance':md(main34['mean'],oo['mean']),'wall_sector_tv':tv(main34['wall_channels'],oo['wall_channels']),'rms_abs_difference':abs(main34['rms_width']-oo['rms_width'])})
    gf4=qgrid((-32,32),(-36,36),240,270);vf4,_=regrid_state(v30,g4,gf4,3);back4,_=regrid_state(vf4,gf4,g4,3);f4=propagate_state(gf4,vf4,30,34,.2,52,[34])['rows'][34]['obs'];controls.append({'case':'grid_240x270','start_regrid_fidelity':fidelity(v30,back4),'mean_distance':md(main34['mean'],f4['mean']),'wall_sector_tv':tv(main34['wall_channels'],f4['wall_channels']),'rms_abs_difference':abs(main34['rms_width']-f4['rms_width'])})
    k_times=np.array([2.,6.,15.,36.]);ks=evolve(np.array([0.,0.,1.,.31]),2,36,t_eval=k_times,rtol=5e-10,atol=5e-12,max_step=.01);u=[kasner_from_state(s,ks.y[:,i])['u'] for i,s in enumerate(k_times)];bkl=[{'from_u':u[i],'predicted':bkl_step(u[i]),'observed_next':u[i+1],'relative_error':abs(bkl_step(u[i])-u[i+1])/u[i+1]} for i in range(3)]
    allobs=[row['obs'] for seg in (seg1,seg2,seg3,seg4) for row in seg['rows'].values()]
    result={'scope':'Same reduced vacuum Bianchi IX principal symbol and same direct square-root quantization as the first-wall audit. Staged finite-box enlargement extends the comparison through the third central classical wall encounter.','wall_sector_definition':'A/B+/B- is whichever positive asymptotic Bianchi IX wall exponential is largest at beta. Sector mass is a one-time position probability, not a full next-wall history probability.','regrids':[r1,r2,r3],'regrid_cumulative_roundtrip_fidelity_product':float(np.prod([r1['roundtrip_fidelity'],r2['roundtrip_fidelity'],r3['roundtrip_fidelity']])),'comparisons':comparisons,'central_bounces':bounce_info,'interference':{'second_wall_split_9p5_to_11':inter2,'third_wall_split_30_to_34':inter3},'numerical_controls':{'second_wall_grid_refinement':second_refine,'third_wall_local_controls':controls,'max_saved_edge_mass_w1':max(o['edge_mass_w1'] for o in allobs),'max_saved_capped_region_mass':max(o['capped_region_mass'] for o in allobs)},'kasner_sequence':{'sample_s':k_times.tolist(),'u':u,'bkl_checks':bkl},'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},'elapsed_s':time.time()-t0}
    compact={'scope':result['scope'],'wall_sector_definition':result['wall_sector_definition'],'regrids':result['regrids'],'regrid_cumulative_roundtrip_fidelity_product':result['regrid_cumulative_roundtrip_fidelity_product'],'central_bounces':result['central_bounces'],'selected_comparisons':[x for x in comparisons if x['s'] in (4.6,10.385,11,14,20,24,30,32.46,34,36)],'interference':result['interference'],'numerical_controls':result['numerical_controls'],'kasner_sequence':result['kasner_sequence']}
    (out/'summary.json').write_text(json.dumps(compact,ensure_ascii=False,indent=2)+'\n');(out/'details.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n');return result

def sig(d):
    vals=[]
    for b in d['central_bounces']: vals+=b['quantum_sector_mass']+b['classical_ensemble_sector_mass']+[b['sector_total_variation']]
    for k in ('second_wall_split_9p5_to_11','third_wall_split_30_to_34'):
        z=d['interference'][k];vals += [z['recombination_fidelity'],z['interference_density_L1_half'],z['channel_total_variation']]
    vals += [d['regrid_cumulative_roundtrip_fidelity_product'],d['numerical_controls']['max_saved_edge_mass_w1'],d['numerical_controls']['max_saved_capped_region_mass']];vals += d['kasner_sequence']['u'];return np.asarray(vals,float)
def check(result,path):
    old=json.loads(Path(path).read_text());a=sig(result);b=sig(old)
    if a.shape!=b.shape or not np.allclose(a,b,rtol=2e-4,atol=3e-6): raise RuntimeError(f'regression mismatch {np.max(abs(a-b)) if a.shape==b.shape else "shape"}')
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,default=ROOT/'quantum_bianchi_ix_long_results');ap.add_argument('--check-against',type=Path);a=ap.parse_args();r=run(a.out)
    if a.check_against:check(r,a.check_against)
    print(json.dumps({'bounces':r['central_bounces'],'interference':r['interference'],'elapsed_s':r['elapsed_s']},indent=2))
if __name__=='__main__':main()
