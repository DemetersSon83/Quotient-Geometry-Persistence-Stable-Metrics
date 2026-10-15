#!/usr/bin/env python3
"""Exact universal rank certificate for every supplied tree on >=2 vertices.

Choose a leaf at the north pole and all remaining vertices at distinct
rational points in the first-quadrant equatorial arc. A fixed quotient
gauge gives a (2n-3)-square minor whose determinant is, up to sign, the
product of equatorial-edge angular sines. This program checks that block
identity using exact rational arithmetic, for representative trees.

The mathematical proof is all-n and is not an inference from these runs.
No third-party packages required. CLI custom input: --edges-json PATH,
where PATH contains either a list of edges or {"n": N, "edges": [...] }.
"""
from fractions import Fraction as Q
import argparse
import heapq
import json
from pathlib import Path
import random

HERE = Path(__file__).resolve().parent


def dot(x, y):
    return sum(a*b for a, b in zip(x, y))


def determinant(a):
    a = [[Q(x) for x in row] for row in a]
    det = Q(1)
    for k in range(len(a)):
        pivotrow = next((i for i in range(k, len(a)) if a[i][k]), None)
        if pivotrow is None:
            return Q(0)
        if pivotrow != k:
            a[k], a[pivotrow] = a[pivotrow], a[k]
            det = -det
        pivot = a[k][k]
        det *= pivot
        for i in range(k+1, len(a)):
            if a[i][k]:
                factor = a[i][k]/pivot
                for j in range(k+1, len(a)):
                    a[i][j] -= factor*a[k][j]
                a[i][k] = Q(0)
    return det


def validate_tree(n, edges):
    assert n >= 2 and len(edges) == n-1
    parent = list(range(n))
    def find(i):
        while parent[i] != i:
            i = parent[i]
        return i
    degree = [0]*n
    for i, j in edges:
        assert 0 <= i < n and 0 <= j < n and i != j
        a, b = find(i), find(j)
        assert a != b, 'Cycle or repeated edge'
        parent[a] = b
        degree[i] += 1; degree[j] += 1
    assert len({find(i) for i in range(n)}) == 1
    return degree


def build_witness(n, edges):
    degree = validate_tree(n, edges)
    leaf = min(i for i in range(n) if degree[i] == 1)
    root = next(j if i == leaf else i for i, j in edges if leaf in (i, j))
    equator = [i for i in range(n) if i != leaf]
    p = [None]*n
    tangents = [None]*n
    ts = {}
    p[leaf] = (Q(0),Q(0),Q(1))
    tangents[leaf] = ((Q(1),Q(0),Q(0)),(Q(0),Q(1),Q(0)))
    for k, i in enumerate(equator, start=1):
        t = Q(k, n+1)
        x, y = (1-t*t)/(1+t*t), 2*t/(1+t*t)
        ts[i] = t
        p[i] = (x,y,Q(0))
        tangents[i] = ((-y,x,Q(0)),(Q(0),Q(0),Q(1)))
    return leaf, root, equator, ts, p, tangents


def full_jacobian(p, tangents, edges):
    n = len(p)
    total = [sum(x[k] for x in p) for k in range(3)]
    js = [[Q(0)]*(2*n) for _ in range(n)]
    for i in range(n):
        for k in range(n):
            for a in range(2):
                val = dot(p[i],tangents[k][a])
                if i == k:
                    val += dot(tangents[i][a],total)
                js[i][2*k+a] = val
    qrows, srows = [], []
    for i,j in edges:
        row = [Q(0)]*(2*n)
        for a in range(2):
            row[2*i+a] = dot(tangents[i][a],p[j])
            row[2*j+a] = dot(p[i],tangents[j][a])
        qrows.append(row)
        srows.append([x+y for x,y in zip(js[i],js[j])])
    return qrows+srows


def certify(n, edges, name='custom'):
    leaf, root, eq, ts, points, tangents = build_witness(n, edges)
    # Oriented edges keep their supplied order. T minus a leaf is a tree.
    eq_indices = [k for k,e in enumerate(edges) if leaf not in e]
    leaf_edge = next(k for k,e in enumerate(edges) if leaf in e)
    eq_edges = [edges[k] for k in eq_indices]
    alpha_vertices = [i for i in eq if i != root]
    columns = [2*i for i in alpha_vertices]+[2*i+1 for i in eq]
    rows = eq_indices+[leaf_edge]+[n-1+k for k in eq_indices]
    df = full_jacobian(points,tangents,edges)
    minor = [[df[i][j] for j in columns] for i in rows]
    na, nb = n-2, n-1
    assert len(minor) == 2*n-3 and all(len(r)==2*n-3 for r in minor)
    # Top-right block is zero; bottom-left need not vanish.
    assert all(minor[i][j] == 0 for i in range(na) for j in range(na,na+nb))
    a = [row[:na] for row in minor[:na]]
    b = [row[na:] for row in minor[na:]]
    assert b[0] == [Q(int(i==root)) for i in eq]
    for row,(i,j) in zip(b[1:],eq_edges):
        assert row == [Q(int(k==i or k==j)) for k in eq]
    factors, normalized = [], []
    for row,(i,j) in zip(a,eq_edges):
        # sin(theta_j-theta_i) = x_i y_j-y_i x_j, nonzero by distinct t.
        sine = points[i][0]*points[j][1]-points[i][1]*points[j][0]
        assert sine
        factors.append(sine)
        norm = [x/sine for x in row]
        assert norm == [Q(int(k==i)-int(k==j)) for k in alpha_vertices]
        normalized.append(norm)
    det_incidence, det_anchor = determinant(normalized), determinant(b)
    assert abs(det_incidence) == 1 and abs(det_anchor) == 1
    cert = det_incidence*det_anchor
    for f in factors:
        cert *= f
    assert cert
    # For small cases additionally check the complete dense minor directly.
    if n <= 8:
        assert determinant(minor) == cert
    return {'name':name,'n':n,'edges':edges,'pole_leaf':leaf,'equator_root':root,
            'equatorial_parameters':{str(i):str(t) for i,t in ts.items()},
            'minor_rows':rows,'minor_columns':columns,'certified_rank':2*n-3,
            'sine_factors':[str(s) for s in factors],
            'oriented_incidence_determinant':str(det_incidence),
            'anchored_unsigned_incidence_determinant':str(det_anchor),
            'minor_determinant_rational':str(cert),
            'complete_minor_checked_directly': n <= 8,
            'exact_block_identity_checked':True}


def prufer_tree(n,rng):
    sequence = [rng.randrange(n) for _ in range(n-2)]
    degree = [1]*n
    for v in sequence:
        degree[v] += 1
    leaves = [i for i,d in enumerate(degree) if d==1]
    heapq.heapify(leaves)
    edges=[]
    for v in sequence:
        leaf = heapq.heappop(leaves)
        edges.append((leaf,v))
        degree[leaf] -= 1; degree[v] -= 1
        if degree[v] == 1:
            heapq.heappush(leaves,v)
    edges.append(tuple(leaves))
    return edges


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--edges-json',type=Path)
    args=parser.parse_args()
    if args.edges_json:
        data=json.loads(args.edges_json.read_text())
        edges=data['edges'] if isinstance(data,dict) else data
        n=data.get('n',max(max(e) for e in edges)+1) if isinstance(data,dict) else max(max(e) for e in edges)+1
        print(json.dumps(certify(n,edges),indent=2))
        return
    rng=random.Random(20261002)
    results=[]
    for n in (3,4,8,16,32):
        examples=[('path',[(i,i+1) for i in range(n-1)]),
                  ('star',[(0,i) for i in range(1,n)]),
                  ('binary', [((i-1)//2,i) for i in range(1,n)]),
                  ('random',prufer_tree(n,rng))]
        for name,edges in examples:
            result=certify(n,edges,name)
            results.append(result)
            print(f'{name:6} n={n:2} certified rank={2*n-3}; block determinant nonzero')
    payload={'statement':'Universal leaf-pole/equator witness; verification supports the all-n structural proof.',
             'seed':20261002,'tested_sizes':[3,4,8,16,32],
             'caveat':'Witness need not realize its tree as MST; generic nonzero rank extends to every nonempty open strict-MST chamber.',
             'results':results}
    (HERE/'results.json').write_text(json.dumps(payload,indent=2)+'\n')


if __name__=='__main__':
    main()
