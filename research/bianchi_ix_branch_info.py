"""Derived wall-sector information audit for the staged quantum Bianchi IX run.

This module does not evolve a new quantum state. It summarizes the already
recorded A/B+/B- one-time sector masses with bounded information measures.
Sector masses are configuration-space diagnostics, not consistent-histories
probabilities for which wall was hit.
"""
from __future__ import annotations
import json, math
from pathlib import Path
import numpy as np


def _prob(p):
    p=np.asarray(p,float)
    if p.shape!=(3,) or np.any(p<0) or not np.all(np.isfinite(p)):
        raise ValueError('three finite nonnegative sector weights required')
    s=float(p.sum())
    if s<=0: raise ValueError('nonzero sector weight required')
    return p/s


def entropy_bits(p):
    p=_prob(p); q=p[p>0]
    return float(-np.sum(q*np.log2(q)))


def effective_branches(p):
    return float(2.0**entropy_bits(p))


def total_variation(p,q):
    p=_prob(p);q=_prob(q)
    return 0.5*float(np.sum(abs(p-q)))


def js_bits(p,q):
    p=_prob(p);q=_prob(q);m=(p+q)/2
    def kl(a,b):
        mask=a>0
        return float(np.sum(a[mask]*np.log2(a[mask]/b[mask])))
    return 0.5*(kl(p,m)+kl(q,m))


def sector_metrics(p):
    p=_prob(p)
    return {
        'probability':p.tolist(),
        'entropy_bits':entropy_bits(p),
        'effective_branch_count':effective_branches(p),
        'dominant_sector_probability':float(p.max()),
        'maximum_entropy_fraction':float(entropy_bits(p)/math.log2(3.0)),
    }


def _interference_block(data, base):
    table=data['interference']
    for key in (base, base+'_k60'):
        if key in table:
            return table[key]
    raise KeyError(f'missing interference block {base!r} (or _k60 variant)')


def _field(block, published, local):
    if published in block: return block[published]
    if local in block: return block[local]
    raise KeyError(f'missing {published!r}/{local!r}')


def audit(data):
    rows={float(r['s']):r for r in data['selected_comparisons']}
    selected=[]
    for s in (4.6,10.385,32.46):
        r=rows[s];q=r['quantum_wall_sector_mass'];c=r['classical_ensemble_wall_sector_mass']
        selected.append({
            's':s,
            'quantum':sector_metrics(q),
            'classical_ensemble':sector_metrics(c),
            'quantum_classical_tv':total_variation(q,c),
            'quantum_classical_js_bits':js_bits(q,c),
        })
    third_controls=data['numerical_controls']['third_wall_local_controls']
    third_num=max(float(x['wall_sector_tv']) for x in third_controls)
    second_num=float(data['numerical_controls']['second_wall_grid_refinement']['wall_sector_tv'])
    inter2=_interference_block(data,'second_wall_split_9p5_to_11')
    inter3=_interference_block(data,'third_wall_split_30_to_34')
    second=next(x for x in selected if x['s']==10.385)
    third=next(x for x in selected if x['s']==32.46)
    inter2_tv=float(_field(inter2,'channel_total_variation','channel_tv'))
    inter3_tv=float(_field(inter3,'channel_total_variation','channel_tv'))
    inter2_l1=float(_field(inter2,'interference_density_L1_half','l1'))
    inter3_l1=float(_field(inter3,'interference_density_L1_half','l1'))
    return {
        'schema':1,
        'scope':'Derived information audit of saved A/B+/B- one-time wall-sector masses; no new dynamics.',
        'selected':selected,
        'interference':{
            'second_wall':{
                'coherent':sector_metrics(inter2['coherent_channels']),
                'incoherent':sector_metrics(inter2['incoherent_channels']),
                'sector_tv':inter2_tv,
                'full_density_l1_half':inter2_l1,
            },
            'third_wall':{
                'coherent':sector_metrics(inter3['coherent_channels']),
                'incoherent':sector_metrics(inter3['incoherent_channels']),
                'sector_tv':inter3_tv,
                'full_density_l1_half':inter3_l1,
            },
        },
        'numerical_scale_comparison':{
            'second_wall_quantum_classical_tv':second['quantum_classical_tv'],
            'second_wall_grid_control_tv':second_num,
            'second_wall_ratio_to_grid_control':float(second['quantum_classical_tv']/second_num),
            'second_wall_interference_sector_tv':inter2_tv,
            'second_wall_interference_ratio_to_grid_control':float(inter2_tv/second_num),
            'third_wall_quantum_classical_tv':third['quantum_classical_tv'],
            'third_wall_max_local_control_tv':third_num,
            'third_wall_ratio_to_max_local_control':float(third['quantum_classical_tv']/third_num),
            'third_wall_interference_sector_tv':inter3_tv,
            'third_wall_interference_ratio_to_max_local_control':float(inter3_tv/third_num),
            'note':'Ratios are scale comparisons, not statistical significance or certified error bars.'
        },
        'limitations':[
            'Sector entropy is not entropy of the full quantum state.',
            'One-time sector mass is not a consistent-histories probability of which wall is hit next.',
            'The three sector projectors are a chosen asymptotic partition of beta-space.',
            'Classical-ensemble spread and quantum coherence are distinct; low-order similarity does not make the states equivalent.',
            'Numerical-control ratios compare observed scales only and are not sigma levels or rigorous continuum error bounds.'
        ]
    }


def main(argv=None):
    import argparse
    ap=argparse.ArgumentParser()
    ap.add_argument('--input',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args(argv)
    d=json.loads(a.input.read_text())
    out=audit(d)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(out['selected'],ensure_ascii=False,indent=2))

if __name__=='__main__': main()
