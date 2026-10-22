import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'python'))

TOL = 0.005            # 0.5 % relative tolerance against the reference


def load_eosfit():
    with open(os.path.join(ROOT, 'reference', 'eosfit_reference.json')) as f:
        return json.load(f)


def load_abaqus():
    import csv
    with open(os.path.join(ROOT, 'reference', 'abaqus_results.csv')) as f:
        rows = [l for l in f if not l.startswith('#')]
    return {r['case']: r for r in csv.DictReader(rows)}
