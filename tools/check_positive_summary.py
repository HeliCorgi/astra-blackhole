"""Regress stored aggregates; detailed row comparison is a separate audit.
Copyright 2026 HeliCorgi. Apache-2.0.
"""
from pathlib import Path
import argparse
import copy
import json
import math

EXACT = ('schema', 'scope', 'base_commit', 'protocol_sha256', 'source_sha256',
         'horizon', 'time_samples', 'cases', 'reference_trajectories',
         'reduced_trajectories', 'positive_verification_trajectories',
         'all_references_valid', 'all_reduced_succeeded',
         'all_positive_verifications_passed', 'paired_comparisons',
         'smallest_accepted')
NUMERIC = {'worst_K_error', 'worst_Pi_error'}


def compare(actual: dict, expected: dict) -> None:
    for key in EXACT:
        if actual[key] != expected[key]:
            raise AssertionError(f'Changed summary field: {key}')
    aa, bb = actual['aggregates'], expected['aggregates']
    if len(aa) != len(bb):
        raise AssertionError('Changed aggregate count')
    for a, b in zip(aa, bb):
        if a.keys() != b.keys():
            raise AssertionError('Changed aggregate schema')
        for key in a:
            if key in NUMERIC and a[key] is not None and b[key] is not None:
                if not math.isclose(a[key], b[key], rel_tol=5e-5, abs_tol=2e-8):
                    raise AssertionError(f'Changed aggregate: {key}')
            elif a[key] != b[key]:
                raise AssertionError(f'Changed aggregate: {key}')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('actual', type=Path)
    parser.add_argument('expected', type=Path)
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args()
    a = json.loads(args.actual.read_text())
    b = json.loads(args.expected.read_text())
    compare(a, b)
    if args.self_test:
        for change in ('count', 'numeric', 'missing'):
            bad = copy.deepcopy(b)
            if change == 'count': bad['aggregates'][0]['accepted'] += 1
            elif change == 'numeric': bad['aggregates'][0]['worst_K_error'] += 1
            else: bad['aggregates'].pop()
            try:
                compare(bad, b)
            except AssertionError:
                continue
            raise AssertionError(f'Guard failed to reject {change}')
        print('Three deliberately altered summaries rejected.')
    print('Aggregate/source regression passed. This is NOT per-row regression.')


if __name__ == '__main__':
    main()
