from pathlib import Path
import re
from collections import defaultdict
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / 'data/raw'
OUT = ROOT / 'data/processed'
SHA = '4eaac889809dd489e3bae468c45d5c0e02b9147a'

def composition(formula):
    """Parse decimal stoichiometry and nested parentheses; reject unparsed syntax."""
    formula = ''.join(formula.split()).replace('[','(').replace(']',')')
    if '*' in formula:
        result = defaultdict(float)
        for part in formula.split('*'):
            for el,n in composition(part).items(): result[el] += n
        return dict(result)
    tokens = re.findall(r'[A-Z][a-z]?|\d*\.\d+|\d+|[()]', formula)
    if ''.join(tokens) != formula.replace(' ', ''):
        raise ValueError(formula)
    stack = [defaultdict(float)]
    i = 0
    while i < len(tokens):
        t = tokens[i]
        if t == '(':
            stack.append(defaultdict(float)); i += 1
        elif t == ')':
            part = stack.pop(); i += 1
            mult = float(tokens[i]) if i < len(tokens) and tokens[i][0] in '.0123456789' else 1.
            i += int(mult != 1. or (i < len(tokens) and tokens[i][0] in '.0123456789'))
            for el, n in part.items(): stack[-1][el] += n * mult
        elif t[0].isupper():
            i += 1
            has_num = i < len(tokens) and tokens[i][0] in '.0123456789'
            mult = float(tokens[i]) if has_num else 1.
            i += int(has_num); stack[-1][t] += mult
        else:
            raise ValueError(formula)
    if len(stack) != 1 or not stack[0]: raise ValueError(formula)
    return dict(stack[0])

def composition_key(formula):
    c = composition(formula); total = sum(c.values())
    return '|'.join(f'{el}:{n/total:.10f}' for el,n in sorted(c.items()) if n > 0)

def chemistry(formula):
    e = set(composition(formula))
    # Exclusive, explicitly heuristic classes; mixed anions are kept separately.
    if 'O' in e and 'S' in e: return 'mixed O/S'
    if 'O' in e: return 'oxide / oxyanion'
    if 'S' in e: return 'sulfide / thioanion'
    if e & {'F','Cl','Br','I'}: return 'halide (no O/S)'
    if 'N' in e: return 'nitride (no O/S/halogen)'
    if 'H' in e: return 'hydride / other H'
    return 'other'

def connected_groups(df):
    """Transitive union of exact composition and individual publication DOIs."""
    parent = list(range(len(df)))
    def find(i):
        while parent[i] != i: parent[i] = parent[parent[i]]; i = parent[i]
        return i
    seen = {}
    for i, row in enumerate(df.itertuples()):
        keys = ['c:'+row.composition_key] + ['d:'+d for d in str(row.doi_key).split(';') if d and d != 'nan']
        for k in keys:
            if k in seen: parent[find(i)] = find(seen[k])
            else: seen[k] = i
    return np.array([find(i) for i in range(len(df))])
