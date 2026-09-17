"""Long-horizon audit helpers for the existing reduced quantum Bianchi IX model.

Extends the first-wall calculation without changing the Hamiltonian.  The only
new numerical device is staged enlargement of the finite beta-plane box, with
cubic interpolation of the physical wavefunction and explicit round-trip
fidelity/norm diagnostics.  Channel labels A/B+/B- are asymptotic wall-sector
projectors, not a theorem about quantum histories or a detector model.
"""
from __future__ import annotations
import numpy as np
from scipy.interpolate import RectBivariateSpline
from quantum_bianchi_ix_model import QuantumBianchiIX, PacketSpec, wall_grid, _rhs_vec, initial_wigner_ensemble

SQ3=np.sqrt(3.0)
WALL_NAMES=('A','B+','B-')

def positive_wall_terms(s,bp,bm):
    return np.stack([
        np.exp(np.clip(-4*s-8*bp,-745.,700.)),
        np.exp(np.clip(-4*s+4*bp+4*SQ3*bm,-745.,700.)),
        np.exp(np.clip(-4*s+4*bp-4*SQ3*bm,-745.,700.)),
    ],axis=0)

def wall_masks(model: QuantumBianchiIX,s: float):
    labels=np.argmax(positive_wall_terms(s,model.BP,model.BM),axis=0)
    return [labels==i for i in range(3)]

def channel_probabilities(model: QuantumBianchiIX,v,s: float):
    p=abs(np.asarray(v).reshape(model.nx,model.ny))**2
    return [float(p[m].sum()) for m in wall_masks(model,s)]

def extended_observables(model: QuantumBianchiIX,v,s: float,edge_width=1.0):
    out=model.observables(v,s)
    p=abs(np.asarray(v).reshape(model.nx,model.ny))**2
    edge=(model.BP<model.bp[0]+edge_width)|(model.BP>model.bp[-1]-edge_width)|(model.BM<model.bm[0]+edge_width)|(model.BM>model.bm[-1]-edge_width)
    out['edge_mass_w1']=float(p[edge].sum())
    out['wall_channels']=channel_probabilities(model,v,s)
    return out

def propagate_state(model,v,s0,s1,ds,krylov_dim,save_times=()):
    save=np.asarray(save_times,float);n=int(round((s1-s0)/ds));ds=(s1-s0)/n;s=float(s0);v=np.asarray(v,complex).copy();rows={};idx=0;min_ritz=np.inf
    if len(save) and abs(save[0]-s0)<1e-10:
        rows[float(save[0])]={'state':v.copy(),'obs':extended_observables(model,v,float(save[0]))};idx=1
    for _ in range(n):
        v,r=model.step(v,s+ds/2,ds,krylov_dim);min_ritz=min(min_ritz,r);s+=ds
        while idx<len(save) and s>=save[idx]-ds/2:
            t=float(save[idx]);rows[t]={'state':v.copy(),'obs':extended_observables(model,v,t)};idx+=1
    return {'state':v,'rows':rows,'min_ritz':float(min_ritz),'ds':float(ds)}

def regrid_state(v,old:QuantumBianchiIX,new:QuantumBianchiIX,order=3):
    """Interpolate continuum psi, then restore sqrt(dx dy) mass amplitudes."""
    psi=np.asarray(v).reshape(old.nx,old.ny)/np.sqrt(old.dx*old.dy)
    sr=RectBivariateSpline(old.bp,old.bm,psi.real,kx=order,ky=order,s=0)
    si=RectBivariateSpline(old.bp,old.bm,psi.imag,kx=order,ky=order,s=0)
    val=(sr(new.bp,new.bm)+1j*si(new.bp,new.bm))*np.sqrt(new.dx*new.dy)
    w=val.ravel();raw=float(np.linalg.norm(w))
    if raw<=0: raise RuntimeError('regrid lost all norm')
    return w/raw,raw

def fidelity(a,b):
    a=np.asarray(a);b=np.asarray(b);return float(abs(np.vdot(a,b))**2/(np.vdot(a,a).real*np.vdot(b,b).real))

def regrid_audit(v,old,new,order=3):
    fwd,n1=regrid_state(v,old,new,order);back,n2=regrid_state(fwd,new,old,order)
    return fwd,{'raw_forward_norm':n1,'raw_back_norm':n2,'roundtrip_fidelity':fidelity(v,back),
                'new_edge_mass_w1':extended_observables(new,fwd,0.0)['edge_mass_w1']}

def branch_interference(model,v_split,s_split,s_end,ds,krylov_dim,full_end=None):
    """Compare coherent propagation with an incoherent mixture of 3 wall sectors."""
    branches=[];weights=[]
    for mask in wall_masks(model,s_split):
        comp=np.where(mask.ravel(),v_split,0)
        weights.append(float(np.vdot(comp,comp).real))
        branches.append(propagate_state(model,comp,s_split,s_end,ds,krylov_dim)['state'])
    coherent=sum(branches)
    if full_end is None:
        full_end=propagate_state(model,v_split,s_split,s_end,ds,krylov_dim)['state']
    dens_full=abs(full_end.reshape(model.nx,model.ny))**2
    dens_mix=sum(abs(b.reshape(model.nx,model.ny))**2 for b in branches)
    qch=[float(dens_full[m].sum()) for m in wall_masks(model,s_end)]
    mch=[float(dens_mix[m].sum()) for m in wall_masks(model,s_end)]
    return {'split_weights':weights,'recombination_fidelity':fidelity(full_end,coherent),
            'interference_density_L1_half':0.5*float(np.sum(abs(dens_full-dens_mix))),
            'coherent_channels':qch,'incoherent_channels':mch,
            'channel_total_variation':0.5*float(np.sum(abs(np.asarray(qch)-np.asarray(mch))))}

def classical_ensemble(packet:PacketSpec,save_times,n=4096,seed=20260917,s0=2.,ds=.01):
    Y=initial_wigner_ensemble(packet,n,seed);save=np.asarray(save_times,float);out={};s=float(s0);idx=0
    def rec(t):
        xy=Y[:2];mean=xy.mean(1);cov=np.cov(xy);terms=positive_wall_terms(t,Y[0][:,None],Y[1][:,None])[:,:,0];labels=np.argmax(terms,axis=0)
        out[float(t)]={'mean':mean.tolist(),'covariance':cov.tolist(),'rms_width':float(np.sqrt(np.trace(cov))),
                       'wall_channels':[float(np.mean(labels==i)) for i in range(3)]}
    if len(save) and abs(save[0]-s0)<1e-10:rec(save[0]);idx=1
    nstep=int(round((save[-1]-s0)/ds));ds=(save[-1]-s0)/nstep
    for _ in range(nstep):
        k1=_rhs_vec(s,Y);k2=_rhs_vec(s+ds/2,Y+ds*k1/2);k3=_rhs_vec(s+ds/2,Y+ds*k2/2);k4=_rhs_vec(s+ds,Y+ds*k3)
        Y+=ds*(k1+2*k2+2*k3+k4)/6;s+=ds
        while idx<len(save) and s>=save[idx]-ds/2: rec(save[idx]);idx+=1
    return out
