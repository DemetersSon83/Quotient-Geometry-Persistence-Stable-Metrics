"""Exact reference certificates for temporal spherical H0 observability.

Only Python's standard library is used.  No floating point operation enters the
rank, chamber, orthonormality, or singular-value certificates.  Decimal fields
in JSON are explicitly display-only.  The reference chart fixes p[0]=e3 and
confines p[1] to the xz meridian; this is transverse to rotations, not horizontal.

Run: python3 reference_certificate.py --output reference_certificates.json
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from itertools import combinations
import json
from math import isqrt, sqrt
from pathlib import Path


def dot(a, b):
    return sum((x*y for x, y in zip(a, b)), F(0))


def stereographic(a, b):
    """Rational unit point; e3 corresponds to a=b=0."""
    a, b = F(a), F(b)
    den = 1+a*a+b*b
    return (2*a/den, 2*b/den, (1-a*a-b*b)/den)


def tangent_frame(p):
    """Rational orthonormal tangent frame, valid away from the south pole."""
    x, y, z = p
    if 1+z == 0:
        raise ValueError("The stereographic tangent frame excludes the south pole")
    a = (1-x*x/(1+z), -x*y/(1+z), -x)
    b = (-x*y/(1+z), 1-y*y/(1+z), -y)
    assert dot(p, a) == dot(p, b) == dot(a, b) == 0
    assert dot(a, a) == dot(b, b) == 1
    return a, b


def gauge_basis(points):
    """Return B_i, each 3 by d, with stacked B exactly orthonormal.

    The associated chart is (p_i+B_i u)/sqrt(1+||B_i u||^2).
    p0 stays fixed, p1 stays in the xz plane.  x1>0 fixes the local meridian.
    """
    n, d = len(points), 2*len(points)-3
    assert tuple(points[0]) == (0, 0, 1)
    x, y, z = points[1]
    assert y == 0 and x > 0
    blocks = [[[F(0) for _ in range(d)] for _ in range(3)] for _ in points]
    for k, a in enumerate((z, F(0), -x)):
        blocks[1][k][0] = a
    for i in range(2, n):
        a, b = tangent_frame(points[i])
        for k in range(3):
            blocks[i][k][2*i-3] = a[k]
            blocks[i][k][2*i-2] = b[k]
    cols = [[blocks[i][k][j] for i in range(n) for k in range(3)] for j in range(d)]
    assert all(dot(cols[j], cols[k]) == int(j == k) for j in range(d) for k in range(d))
    assert all(dot(points[i], [blocks[i][k][j] for k in range(3)]) == 0
               for i in range(n) for j in range(d))
    return blocks


def mst_edges(points):
    """Exact Kruskal MST: increasing geodesic distance = decreasing dot product."""
    parent = list(range(len(points)))
    def find(i):
        while parent[i] != i:
            i = parent[i]
        return i
    out = []
    for i, j in sorted(combinations(range(len(points)), 2),
                       key=lambda e: -dot(points[e[0]], points[e[1]])):
        a, b = find(i), find(j)
        if a != b:
            parent[a] = b
            out.append((i, j))
    return out


def exact_jet_jacobian(points, blocks, edges):
    """Jacobian of (q_e, (1-q_e)(s_i+s_j)) at u=0."""
    n, d = len(points), len(blocks[0][0])
    S = [sum(p[k] for p in points) for k in range(3)]
    ds = [[F(0)]*d for _ in range(n)]
    for c in range(d):
        dS = [sum(blocks[i][k][c] for i in range(n)) for k in range(3)]
        for i in range(n):
            vi = [blocks[i][k][c] for k in range(3)]
            ds[i][c] = dot(vi, S)+dot(points[i], dS)
    qrows, rrows = [], []
    for i, j in edges:
        q, s = dot(points[i], points[j]), dot(points[i], S)+dot(points[j], S)
        dq = [sum(blocks[i][k][c]*points[j][k]+points[i][k]*blocks[j][k][c]
                  for k in range(3)) for c in range(d)]
        qrows.append(dq)
        rrows.append([-dq[c]*s+(1-q)*(ds[i][c]+ds[j][c]) for c in range(d)])
    return qrows+rrows


def gram_matrix(A):
    d = len(A[0])
    return [[sum((row[i]*row[j] for row in A), F(0)) for j in range(d)] for i in range(d)]


def inverse_matrix(A):
    n = len(A)
    work = [list(row)+[F(int(i == j)) for j in range(n)] for i, row in enumerate(A)]
    for c in range(n):
        pivot = next((r for r in range(c, n) if work[r][c]), None)
        if pivot is None:
            raise ValueError("Singular matrix")
        work[c], work[pivot] = work[pivot], work[c]
        den = work[c][c]
        work[c] = [x/den for x in work[c]]
        for r in range(n):
            if r != c and work[r][c]:
                den = work[r][c]
                work[r] = [a-den*b for a, b in zip(work[r], work[c])]
    assert all(work[i][j] == int(i == j) for i in range(n) for j in range(n))
    inv = [row[n:] for row in work]
    # Verify against the input independently of the row-reduced left block.
    assert all(sum((A[i][k]*inv[k][j] for k in range(n)), F(0)) == int(i == j)
               for i in range(n) for j in range(n))
    return inv


def reciprocal_sqrt_lower(x, bits=32):
    """Dyadic lower bound for 1/sqrt(x), with an exact squared inequality."""
    if x <= 0:
        raise ValueError("Positive argument required")
    ans = F(isqrt((x.denominator << (2*bits))//x.numerator), 1 << bits)
    assert ans > 0 and ans*ans*x <= 1
    return ans


def sqrt_upper(x, bits=32):
    """Dyadic upper square-root enclosure using integer arithmetic."""
    if x < 0:
        raise ValueError("Nonnegative argument required")
    scaled_num = x.numerator << (2*bits)
    k = isqrt(scaled_num//x.denominator)
    if k*k*x.denominator < scaled_num:
        k += 1
    ans = F(k, 1 << bits)
    assert ans*ans >= x
    return ans


def ceil_sqrt_integer(n):
    k = isqrt(n)
    return k if k*k == n else k+1


def rational_rank(A):
    work = [list(row) for row in A]
    rank = 0
    for c in range(len(A[0])):
        pivot = next((r for r in range(rank, len(work)) if work[r][c]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        den = work[rank][c]
        work[rank] = [x/den for x in work[rank]]
        for r in range(rank+1, len(work)):
            den = work[r][c]
            if den:
                work[r] = [a-den*b for a, b in zip(work[r], work[rank])]
        rank += 1
    return rank


def verify_jet_via_flow(points, blocks, edges, A):
    """Independent exact differentiation of f_i=S−(p_i.S)p_i."""
    n, d = len(points), len(blocks[0][0])
    S = [sum(p[k] for p in points) for k in range(3)]
    s = [dot(p, S) for p in points]
    f = [[S[k]-s[i]*p[k] for k in range(3)] for i, p in enumerate(points)]
    for c in range(d):
        v = [[blocks[i][k][c] for k in range(3)] for i in range(n)]
        dS = [sum(v[i][k] for i in range(n)) for k in range(3)]
        ds = [dot(v[i], S)+dot(points[i], dS) for i in range(n)]
        df = [[dS[k]-ds[i]*points[i][k]-s[i]*v[i][k] for k in range(3)] for i in range(n)]
        for r, (i, j) in enumerate(edges):
            direct = dot(df[i], points[j])+dot(f[i], v[j])+dot(v[i], f[j])+dot(points[i], df[j])
            assert direct == A[len(edges)+r][c]


def quantitative_bounds(n, x, gap, collision_margin, antipodal_margin, sigma):
    """Exact conservative bounds from the accompanying quantitative theorem.

    M,N are integer upper bounds for sqrt(n−1),sqrt(n).  K is an upper bound
    at the chosen h=tau, not an estimate of the best achievable condition number.
    """
    M, N = ceil_sqrt_integer(n-1), ceil_sqrt_integer(n)
    mu = min(gap, collision_margin, antipodal_margin)/2
    rho = min(F(1, 2), x/2, mu/16, sigma/(248*n*M))
    tau = min(F(1, 3*n), mu/(16*n), sigma/(6*n*n*M*(13*N+24)))
    eta = rho*x/(12*N)
    K = 4*M/sigma*sqrt_upper(1+4/(tau*tau))
    return {"mu": mu, "M": F(M), "N": F(N), "rho": rho, "tau": tau,
            "eta": eta, "K_at_tau_upper": K}


def rational_string(x):
    return str(F(x))


def encode_matrix(A):
    return [[rational_string(x) for x in row] for row in A]


def certificate(points, label=None):
    """Compute independently checkable exact data for one rational formation.

    The bound follows because sigma_min(A)^(-2) <= tr((A^T A)^(-1)).
    Finite-time and neighborhood bounds apply the analytic constants proved
    in the accompanying theorem; no numerical Taylor fit is used.
    """
    points = [tuple(F(x) for x in p) for p in points]
    n, d = len(points), 2*len(points)-3
    assert n >= 3 and all(dot(p, p) == 1 for p in points)
    pairs = sorted(((dot(points[i], points[j]), i, j)
                    for i, j in combinations(range(n), 2)), reverse=True)
    qs = [v[0] for v in pairs]
    if len(set(qs)) != len(qs) or min(qs) <= -1 or max(qs) >= 1:
        raise ValueError("Certificate requires distinct nonantipodal points and strict edge order")
    gap = min(qs[i]-qs[i+1] for i in range(len(qs)-1))
    edges = mst_edges(points)
    B = gauge_basis(points)
    A = exact_jet_jacobian(points, B, edges)
    verify_jet_via_flow(points, B, edges, A)
    G = gram_matrix(A)
    Ginv = inverse_matrix(G)
    trace = sum((Ginv[i][i] for i in range(d)), F(0))
    sigma = reciprocal_sqrt_lower(trace)
    bounds = quantitative_bounds(n, points[1][0], gap, 1-max(qs), 1+min(qs), sigma)
    # Since G=A^T A and invertible, it is positive definite; no numerical
    # eigenvalue threshold or unproved rank decision appears in this step.
    return {
        "label": label or f"rational_n{n}", "n": n, "d": d,
        "exact": {
            "points": encode_matrix(points),
            "tangent_blocks": [encode_matrix(block) for block in B],
            "mst_edges": [list(e) for e in edges],
            "ordered_pairs": [[i, j] for _, i, j in pairs],
            "ordered_pair_dots": [rational_string(q) for q in qs],
            "dot_order_gap": rational_string(gap),
            "collision_dot_margin": rational_string(1-max(qs)),
            "antipodal_dot_margin": rational_string(1+min(qs)),
            "gauge_p1_x": rational_string(points[1][0]),
            "jet_jacobian": encode_matrix(A),
            "inverse_gram": encode_matrix(Ginv),
            "trace_inverse_gram": rational_string(trace),
            "sigma_lower": rational_string(sigma),
            "sigma_lower_squared_times_trace": rational_string(sigma*sigma*trace),
            "quantitative_bounds": {key: rational_string(value) for key, value in bounds.items()},
        },
        "display_only": {
            "dot_order_gap": float(gap),
            "collision_dot_margin": float(1-max(qs)),
            "antipodal_dot_margin": float(1+min(qs)),
            "trace_inverse_gram": float(trace),
            "sigma_lower": float(sigma),
            "trace_bound": 1/sqrt(float(trace)),
            "quantitative_bounds": {key: float(value) for key, value in bounds.items()},
        },
        "checks": {
            "unit_points_exact": True,
            "tangent_basis_orthonormal_exact": True,
            "strict_pair_order_exact": True,
            "inverse_gram_verified_exact": True,
            "sigma_squared_times_trace_at_most_one_exact": True,
            "qdot_jacobian_independently_verified_from_flow_exact": True,
            "rank": d,
        },
    }


def reference_points(n):
    parameters = [(F(0), F(0)), (F(1, 2), F(0)),
                  (F(1, 5), F(3, 5)), (F(-3, 5), F(2, 7)),
                  (F(-2, 5), F(-4, 5)), (F(4, 5), F(-3, 5)),
                  (F(5, 4), F(2, 5)), (F(-3, 2), F(-1, 4))]
    if not 3 <= n <= len(parameters):
        raise ValueError("This illustrative family supports 3 <= n <= 8")
    return [stereographic(a, b) for a, b in parameters[:n]]


def centered_points():
    """Six nonantipodal rational points summing exactly to zero.

    Reflect a rational triple across the plane perpendicular to its sum. This
    negative control has repeated pair distances and is outside the strict
    chamber theorem as well as the full-rank locus.
    """
    first = reference_points(3)
    S = [sum(v[k] for v in first) for k in range(3)]
    den = dot(S, S)
    return first+[tuple(v[k]-2*dot(v, S)*S[k]/den for k in range(3)) for v in first]


def near_centered_points(epsilon=F(1, 10000)):
    p = centered_points()
    for i in range(3, 6):
        x, y, z = p[i]
        p[i] = stereographic(x/(1+z)+(i+1)*epsilon, y/(1+z)+(i+2)*epsilon)
    return p


def centered_negative_control():
    p = centered_points()
    S = [sum(v[k] for v in p) for k in range(3)]
    assert S == [0, 0, 0]
    B = gauge_basis(p)
    E = mst_edges(p)
    A = exact_jet_jacobian(p, B, E)
    verify_jet_via_flow(p, B, E, A)
    rank = rational_rank(A)
    assert rank <= len(p)+2 < 2*len(p)-3
    return {"label": "centered_n6_negative_control", "n": 6,
            "points": encode_matrix(p), "centroid_sum": encode_matrix([S])[0],
            "chosen_mst_edges": [list(e) for e in E], "jet_rank_exact": rank,
            "quotient_dimension": 9, "jet_jacobian": encode_matrix(A),
            "no_coincident_or_antipodal_pairs_exact": all(-1 < dot(a, b) < 1 for a, b in combinations(p, 2)),
            "scope": "Centered exact rank obstruction. Repeated pair distances also violate strict chamber assumption; no quantitative certificate is claimed."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("reference_certificates.json"))
    parser.add_argument("--sizes", nargs="+", type=int, default=[3, 4, 6, 8])
    args = parser.parse_args()
    out = {
        "schema": "exact-reference-jet-certificate-v1",
        "method": "Rational meridian gauge; exact Gram inverse; reciprocal trace eigenvalue bound",
        "scope": "Exact reference data and conservative validity bounds; see accompanying theorem for hypotheses and constants",
        "certificates": [certificate(reference_points(n)) for n in args.sizes]+[certificate(near_centered_points(), "near_centered_n6_epsilon_1_over_10000")],
        "negative_controls": [centered_negative_control()],
    }
    args.output.write_text(json.dumps(out, indent=2)+"\n")
    for item in out["certificates"]:
        disp = item["display_only"]
        print(f"n={item['n']}: exact rank {item['d']}; sigma >= {disp['sigma_lower']:.9g}; dot gap {disp['dot_order_gap']:.9g}")
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
