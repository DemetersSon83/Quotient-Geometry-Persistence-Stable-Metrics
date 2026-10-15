"""Uniform local interval certificate for temporal spherical persistence.

Reads exact reference certificates. All acceptance decisions use exact rational
arithmetic or outward-rounded integer dyadic intervals (96 fractional bits).
Floating point fields and printed numbers are display-only.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import json
import time
from interval_core import I,AD,dot,norm_upper,matmul,row_sum_norm_upper,BITS
from reference_certificate import gauge_basis,exact_jet_jacobian,mst_edges

ROOT=Path(__file__).resolve().parent


def fractions(a):
    return [fractions(x) for x in a] if isinstance(a,list) else F(a)


def exact_left_inverse(ex):
    p=fractions(ex['points']);Z=fractions(ex['tangent_blocks'])
    assert all(dot(pi,pi)==1 for pi in p)
    assert Z==gauge_basis(p)
    assert [list(e) for e in mst_edges(p)]==ex['mst_edges']
    A=fractions(ex['jet_jacobian'])
    assert A==exact_jet_jacobian(p,Z,ex['mst_edges'])
    Ginv=fractions(ex['inverse_gram'])
    d,m2=len(A[0]),len(A)
    L=[[sum((Ginv[i][k]*A[j][k] for k in range(d)),F(0))
        for j in range(m2)] for i in range(d)]
    assert all(sum((L[i][k]*A[k][j] for k in range(m2)),F(0))==int(i==j)
               for i in range(d) for j in range(d))
    return L


def chart_ad(p,Z,rho):
    n,d=len(p),len(Z[0][0])
    u=[AD.variable(I(-rho,rho),d,k) for k in range(d)]
    points=[]
    for i in range(n):
        # The reference blocks are orthonormal with disjoint coordinate support.
        ids=[] if i==0 else ([0] if i==1 else [2*i-3,2*i-2])
        rad=AD.constant(1,d)+sum((u[k].square() for k in ids),AD.constant(0,d))
        weight=rad.invsqrt()
        points.append([(AD.constant(p[i][j],d)+sum((Z[i][j][k]*u[k] for k in ids),AD.constant(0,d)))*weight
                       for j in range(3)])
    return points


def field(points):
    S=[sum(p[k] for p in points) for k in range(3)]
    s=[dot(p,S) for p in points]
    f=[[S[k]-s[i]*p[k] for k in range(3)] for i,p in enumerate(points)]
    return S,s,f


def jet(points,edges):
    S=[sum(p[k] for p in points) for k in range(3)]
    s=[dot(p,S) for p in points]
    q=[dot(points[i],points[j]) for i,j in edges]
    b=[(1-q[k])*(s[i]+s[j]) for k,(i,j) in enumerate(edges)]
    return q,b


def acceleration(points,edges):
    """L_f B, a polynomial extension valid for all ambient points."""
    S,s,f=field(points)
    Fsum=[sum(v[k] for v in f) for k in range(3)]
    ds=[dot(f[i],S)+dot(p,Fsum) for i,p in enumerate(points)]
    return [(1-dot(points[i],points[j]))*(ds[i]+ds[j]-(s[i]+s[j]).square())
            for i,j in edges]


def spatial_error(p,Z,edges,L,rho):
    q,b=jet(chart_ad(p,Z,rho),edges)
    jac=[a.g for a in q+b]
    E=matmul(L,jac)
    for i in range(len(E)):
        E[i][i]=E[i][i]-1
    return norm_upper(E),E


def tube_box(p,R):
    return [[I(x-R,x+R).intersect(I(-1,1)) for x in point] for point in p]


def speed_upper(p,R):
    _,_,f=field(tube_box(p,R))
    return max(I(sum(v.absmax()**2 for v in vec)).sqrt().upper() for vec in f)


def pick_tube(p,rho,tau):
    # A final strict first-exit inequality proves the enclosure; iteration
    # is only a search heuristic and is not itself accepted as proof.
    R=2*(rho+len(p)*tau)
    for _ in range(4):
        V=min(F(len(p)),speed_upper(p,R))
        proposal=F(101,100)*(rho+tau*V)
        if proposal<R:R=proposal
    V=min(F(len(p)),speed_upper(p,R))
    return R,V,bool(rho+tau*V<R)


def time_bound(p,edges,L1,R,tau):
    n,k=len(p),3*len(p)
    box=tube_box(p,R)
    X=[[AD.variable(box[i][j],k,3*i+j) for j in range(3)] for i in range(n)]
    C=acceleration(X,edges)
    H=norm_upper(matmul(L1,[c.g for c in C]))
    _,_,f=field(X)
    lam=min(F(3*n),norm_upper([a.g for v in f for a in v]))
    if tau*lam>=1:return None,H,lam
    beta=tau*H/(2*(1-lam*tau))
    return beta,H,lam


def tree_path(edges,source,target):
    adj={i:[] for e in edges for i in e}
    for k,(i,j) in enumerate(edges):adj[i].append((j,k));adj[j].append((i,k))
    stack=[(source,-1,[])]
    while stack:
        i,prev,path=stack.pop()
        if i==target:return path
        stack.extend((j,i,path+[k]) for j,k in adj[i] if j!=prev)
    raise ValueError('Disconnected tree')


def chamber(p,edges,R):
    n=len(p); bounds={};ref={}
    for i in range(n):
        for j in range(i+1,n):
            q=dot(p[i],p[j]);ref[(i,j)]=q
            sine=I(1-q*q).sqrt().upper()
            width=2*sine*R+2*R*R
            bounds[(i,j)]=I(q-width,q+width).intersect(I(-1,1))
    E=[tuple(sorted(e)) for e in edges]
    order_slack=min(bounds[E[k]].lower()-bounds[E[k+1]].upper() for k in range(len(E)-1))
    cycle_slack=min((min(bounds[E[k]].lower() for k in tree_path(edges,i,j))-bounds[(i,j)].upper()
                     for i,j in bounds if (i,j) not in E),default=F(1))
    endpoint_slack=min(min(1-q.upper(),1+q.lower()) for q in bounds.values())
    death_lower=(1-max(ref.values()))/2
    diagonal_slack=death_lower-5*R
    ok=min(order_slack,cycle_slack,endpoint_slack,diagonal_slack)>0
    return ok,dict(mst_order_slack=order_slack,mst_cycle_slack=cycle_slack,
                   endpoint_slack=endpoint_slack,diagonal_slack=diagonal_slack)


def attempt(c,rho,tau,context=None,details=False):
    ex=c['exact'];p=fractions(ex['points']);Z=fractions(ex['tangent_blocks']);edges=ex['mst_edges']
    L=exact_left_inverse(ex) if context is None else context
    m=len(edges);L0=[r[:m] for r in L];L1=[r[m:] for r in L]
    if rho>min(F(1,2),p[1][0]/2) or min(rho,tau)<=0:
        return dict(accepted=False,reason='gauge or positive-domain condition')
    R,V,tube_ok=pick_tube(p,rho,tau)
    ch_ok,slacks=chamber(p,edges,R)
    if not tube_ok or not ch_ok:
        return dict(accepted=False,reason='tube or persistence chamber',R=R,**slacks)
    alpha,spatial_matrix=spatial_error(p,Z,edges,L,rho)
    if alpha>=F(1,2) and not details:
        return dict(accepted=False,reason='spatial derivative',alpha=alpha)
    beta,H,lam=time_bound(p,edges,L1,R,tau)
    if beta is None:return dict(accepted=False,reason='variational denominator')
    theta=alpha+beta
    eta=rho*p[1][0]/(12*F(ex['quantitative_bounds']['N']))
    c0,c1=row_sum_norm_upper(L0),row_sum_norm_upper(L1)
    K=(2*(c0+2*c1/tau)/(1-theta)) if theta<1 else None
    ans=dict(accepted=theta<=F(1,2),rho=rho,tau=tau,eta=eta,tube_radius=R,
             speed_bound=V,tube_slack=R-rho-tau*V,spatial_error=alpha,
             temporal_error=beta,total_error=theta,preconditioned_acceleration_bound=H,
             flow_derivative_bound=lam,flow_denominator=1-lam*tau,
             left_block_norm=c0,right_block_norm=c1,K_at_tau_upper=K,**slacks)
    if details:
        ans['left_inverse']=L
        ans['spatial_matrix_intervals']=[[[x.lower(),x.upper()] for x in row] for row in spatial_matrix]
    return ans


def search(c):
    old=c['exact']['quantitative_bounds'];r0,t0=F(old['rho']),F(old['tau'])
    L=exact_left_inverse(c['exact'])
    lo,hi=F(0),F(1)
    history=[]
    for _ in range(18):
        result=attempt(c,r0*hi,t0*hi,L)
        history.append(dict(scale=str(hi),accepted=result['accepted'],reason=result.get('reason','combined derivative')))
        if not result['accepted']:break
        lo,hi=hi,hi*2
    else:raise RuntimeError('Search cap reached; no failed upper candidate')
    for _ in range(9):
        mid=(lo+hi)/2
        result=attempt(c,r0*mid,t0*mid,L)
        history.append(dict(scale=str(mid),accepted=result['accepted'],reason=result.get('reason','combined derivative')))
        if result['accepted']:lo=mid
        else:hi=mid
    # Leave modest search slack. The final chosen result is separately checked.
    scale=lo*F(19,20)
    selected_rho=r0*scale
    base_tau=t0*scale
    # With the enlarged radius fixed, use the remaining derivative budget
    # for a larger delay. Every final candidate is validated independently.
    delay_lo,delay_hi=F(1),F(2)
    delay_history=[]
    for _ in range(24):
        result=attempt(c,selected_rho,base_tau*delay_hi,L)
        delay_history.append(dict(scale=str(delay_hi),accepted=result['accepted'],reason=result.get('reason','combined derivative')))
        if not result['accepted']:break
        delay_lo,delay_hi=delay_hi,delay_hi*2
    else:raise RuntimeError('Delay search cap reached')
    for _ in range(9):
        mid=(delay_lo+delay_hi)/2
        result=attempt(c,selected_rho,base_tau*mid,L)
        delay_history.append(dict(scale=str(mid),accepted=result['accepted'],reason=result.get('reason','combined derivative')))
        if result['accepted']:delay_lo=mid
        else:delay_hi=mid
    delay_scale=max(F(1),F(19,20)*delay_lo)
    final=attempt(c,selected_rho,base_tau*delay_scale,L,details=True)
    assert final['accepted']
    final.update(scale=scale,radius_gain=scale,delay_gain=scale*delay_scale,
                 label=c['label'],n=c['n'],search_history=history,delay_search_history=delay_history,
                 original_bounds={k:F(v) for k,v in old.items()})
    return final


def encode(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:encode(v) for k,v in x.items()}
    if isinstance(x,list):return [encode(v) for v in x]
    return x


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--labels',nargs='*')
    parser.add_argument('--output',type=Path,default=ROOT/'local_certificates.json')
    args=parser.parse_args()
    references=json.loads((ROOT/'reference_certificates.json').read_text())['certificates']
    if args.labels:references=[c for c in references if c['label'] in args.labels]
    results=[];start=time.time()
    for c in references:
        r=search(c);results.append(r)
        print(json.dumps({k:(float(r[k]) if isinstance(r[k],F) else r[k])
                          for k in ['label','radius_gain','delay_gain','rho','tau','eta','spatial_error','temporal_error','K_at_tau_upper']}),flush=True)
        out=dict(scope='Uniform interval certificate; search is heuristic, final acceptance is rigorously enclosed',
                 interval_bits=BITS,results=results,elapsed_seconds=time.time()-start)
        args.output.write_text(json.dumps(encode(out),indent=2)+'\n')
    print(f'Saved {args.output}; elapsed {time.time()-start:.2f}s',flush=True)


if __name__=='__main__':main()
