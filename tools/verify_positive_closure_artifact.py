"""Check result structure, source fingerprints, and summary counts.
Copyright 2026 HeliCorgi. Apache-2.0.
This does not replace re-running the numerical experiment.
"""
from pathlib import Path
import argparse, hashlib, json, sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'research'))
from run_positive_closure_audit import summary_of


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('results',type=Path)
    a=p.parse_args();s=json.loads((a.results/'summary.json').read_text())
    rows=json.loads((a.results/'cases.json').read_text());q=json.loads((a.results/'reference_quality.json').read_text())
    protocol=json.loads((ROOT/'research/positive_closure_protocol.json').read_text())
    assert hashlib.sha256((ROOT/'research/positive_closure_protocol.json').read_bytes()).hexdigest()==s['protocol_sha256']
    for filename,expected in s['source_sha256'].items():
        assert hashlib.sha256((ROOT/filename).read_bytes()).hexdigest()==expected,filename
    assert len(rows)==450 and len({(r['case'],r['method']) for r in rows})==450
    assert len(q)==90 and all(x['valid'] for x in q)
    recomputed=summary_of(rows,q,protocol)
    for k,v in recomputed.items():assert s[k]==v,k
    assert all(x['integration_success'] for x in rows)
    assert all(x['verification_pass'] for x in rows)
    print('450 unique result rows, 90 reference checks, summary counts and source hashes verified.')
    print('This is structural/source verification, NOT a new numerical run.')
if __name__=='__main__':main()
