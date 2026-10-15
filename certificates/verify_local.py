"""Exact arithmetic checks of the interval kernel and certificate formulas.

The temporal derivative is cross-checked using an independently written
Hessian/product-rule formula, not the acceleration-expression AD implementation.
"""
from fractions import Fraction as F
from pathlib import Path
import json
import random
from interval_core import I,AD,dot
from local_certificate import (fractions,chart_ad,jet,acceleration,field,
                               exact_left_inverse,attempt)

ROOT=Path(__file__).resolve().parent


def elementary_checks():
    rng=random.Random(20261002)
    for _ in range(300):
        a,b=sorted(F(rng.randint(-1000,1000),rng.randint(1,100)) for _ in range(2))
        c,d=sorted(F(rng.randint(-1000,1000),rng.randint(1,100)) for _ in range(2))
        x,y=I(a,b),I(c,d)
        assert x.contains(a) and x.contains(b)
        assert all((x+y).contains(z) for z in (a+c,b+d))
        assert all((x*y).contains(z) for z in (a*c,a*d,b*c,b*d))
        assert x.square().contains(a*a) and x.square().contains(b*b)
        if a<=0<=b:assert x.square().contains(0)
        if not c<=0<=d:
            assert all((x/y).contains(z) for z in (a/c,a/d,b/c,b/d))
        low,high=abs(a),abs(b)
        low,high=min(low,high),max(low,high)
        z=I(low,high).sqrt()
        assert z.lower()>=0 and z.lower()**2<=low and z.upper()**2>=high
    u=AD.variable(F(3,4),1,0)
    z=(1+u.square()).invsqrt()
    assert z.v.contains(F(4,5)) and z.g[0].contains(-F(48,125))
    return 300


def direct_acceleration_derivative(p,edges):
    n=len(p);dim=3*n
    S=[sum(pi[k] for pi in p) for k in range(3)]
    s=[dot(pi,S) for pi in p]
    f=[[S[k]-s[i]*p[i][k] for k in range(3)] for i in range(n)]
    FS=[sum(fi[k] for fi in f) for k in range(3)]
    rows=[[F(0)]*dim for _ in edges]
    for c in range(dim):
        v=[[F(int(3*i+k==c)) for k in range(3)] for i in range(n)]
        VS=[sum(vi[k] for vi in v) for k in range(3)]
        ds=[dot(v[i],S)+dot(p[i],VS) for i in range(n)]
        df=[[VS[k]-ds[i]*p[i][k]-s[i]*v[i][k] for k in range(3)] for i in range(n)]
        DFS=[sum(di[k] for di in df) for k in range(3)]
        for r,(i,j) in enumerate(edges):
            q=dot(p[i],p[j]);t=s[i]+s[j]
            pq=[p[i][k]+p[j][k] for k in range(3)]
            vf=[f[i][k]+f[j][k] for k in range(3)]
            vv=[v[i][k]+v[j][k] for k in range(3)]
            vdf=[df[i][k]+df[j][k] for k in range(3)]
            dqv=dot(v[i],p[j])+dot(p[i],v[j])
            dqf=dot(f[i],p[j])+dot(p[i],f[j])
            dtv=dot(vv,S)+dot(pq,VS)
            dtf=dot(vf,S)+dot(pq,FS)
            ddq=dot(f[i],v[j])+dot(v[i],f[j])
            ddt=dot(vf,VS)+dot(vv,FS)
            hess=-ddq*t-dqf*dtv-dqv*dtf+(1-q)*ddt
            dqdf=dot(df[i],p[j])+dot(p[i],df[j])
            dtdf=dot(vdf,S)+dot(pq,DFS)
            rows[r][c]=hess-dqdf*t+(1-q)*dtdf
    return rows


def main():
    elementary=elementary_checks()
    refs=json.loads((ROOT/'reference_certificates.json').read_text())['certificates']
    outputs=json.loads((ROOT/'local_certificates.json').read_text())['results']
    bylabel={x['label']:x for x in outputs}
    results=[]
    for c in refs:
        ex=c['exact'];p=fractions(ex['points']);Z=fractions(ex['tangent_blocks']);E=ex['mst_edges']
        A=fractions(ex['jet_jacobian']);L=exact_left_inverse(ex)
        q,b=jet(chart_ad(p,Z,F(0)),E)
        assert all(y.contains(x) for r,row in enumerate(q+b) for x,y in zip(A[r],row.g))
        dim=3*len(p)
        X=[[AD.variable(p[i][k],dim,3*i+k) for k in range(3)] for i in range(len(p))]
        ad=acceleration(X,E)
        direct=direct_acceleration_derivative(p,E)
        assert all(y.contains(x) for row,a in zip(direct,ad) for x,y in zip(row,a.g))
        old=bylabel[c['label']]
        fresh=attempt(c,F(old['rho']),F(old['tau']),L,details=True)
        assert fresh['accepted']
        for key,value in fresh.items():
            if isinstance(value,F):assert value==F(old[key])
        # Negative certificate controls: do not certify large
        # domains merely because the reference rank is full.
        enlarged=attempt(c,100*F(old['rho']),100*F(old['tau']),L)
        assert not enlarged['accepted']
        results.append(dict(label=c['label'],reference_jet_enclosure=True,
                            temporal_derivative_hessian_crosscheck=True,
                            final_certificate_recomputed=True,
                            enlarged_candidate_not_certified=True,
                            exact_total_error=str(fresh['total_error'])))
    report=dict(elementary_rational_interval_cases=elementary,
                inverse_sqrt_derivative_exact_check=True,results=results,
                scope='Exact software cross-checks; not a proof-assistant formalization or external peer review')
    (ROOT/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
