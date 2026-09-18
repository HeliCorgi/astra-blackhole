from __future__ import annotations
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
d=json.loads((ROOT/'research/bianchi_ix_ever_entered_results/summary.json').read_text())
rows={x['case']:x for x in d['quantum_scan']}
assert len(rows)==14
assert d['decision']['strict_decoherence_accepted'] is False
assert d['decision']['probability_interpretation_accepted'] is False
assert abs(d['classical_control']['everB_given_A']-0.9987293519695044)<1e-12
assert rows['v006_k78']['offdiag_abs']>rows['v0075_k78']['offdiag_abs']
assert rows['v0085_k78']['offdiag_abs']>rows['v0075_k78']['offdiag_abs']
assert rows['v010_k78']['offdiag_abs']>rows['v0075_k78']['offdiag_abs']
assert rows['v0075_w020_k78']['offdiag_real']<0<rows['v0075_w045_k78']['offdiag_real']
assert abs(rows['v005_k66']['offdiag_abs']-rows['v005_k78']['offdiag_abs'])<0.003
assert abs(rows['v0075_k78']['offdiag_abs']-rows['v0075_k90']['offdiag_abs'])<0.002
print(json.dumps({'cases':len(rows),'decision':d['decision'],'classical_everB':d['classical_control']['everB_given_A']},ensure_ascii=False,indent=2))
