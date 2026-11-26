#!/usr/bin/env python3
"""Independent replay of the [T2026] certificate for C_59 = K_2 < 0.302825279490988.

Written from the 59a constant page and the certificate JSON only; the claimant's
checker (t7428_c59.py) was neither read nor run.

Method (differs from the claimant's and from the #297 replay's state-space recursion):
  1. Exact Gaussian-integer check that U = Ut/den is unitary (U*U = UU* = I).
  2. f = A + B Delta (I - D Delta)^{-1} C = P/Q with Q = det(I - D Delta), a
     polynomial of bidegree (1,1). Coefficients of f come from exact power-series
     division Q f = P, kept as Gaussian integers g_jk = f_jk * q0^(j+k+1).
  3. Each |f_jk| r^(j+k) is enclosed with mpmath interval arithmetic (iv, 200 digits);
     the decision uses only the interval's lower endpoint.
  4. Canaries that must FAIL: r below the search limit, a one-unit change in Ut,
     truncation at N = 20, and an inflated claimed margin.
Usage: python3 -I replay_c59_t2026.py path/to/c59_cert.json
"""
import json, sys, hashlib, time
from mpmath import iv, mp, mpf, mpc

EXPECTED_SHA = "6d2e39e8741fb5e73593d4db007c68831243d32fe4ae298bc4fe750033ebedab"

def load(path):
    raw = open(path, "rb").read()
    d = json.loads(raw)
    den = int(d["den"])
    Ut = [[complex_int(e) for e in row] for row in d["U"]]
    R, E = (int(x) for x in d["r"])
    return hashlib.sha256(raw).hexdigest(), d, den, Ut, R, E

def complex_int(e):
    return (int(e[0]), int(e[1]))

# Gaussian integers as (re, im) tuples
def gm(a, b): return (a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def ga(a, b): return (a[0]+b[0], a[1]+b[1])
def gs(a, b): return (a[0]-b[0], a[1]-b[1])
def gk(k, a): return (k*a[0], k*a[1])
def conj(a): return (a[0], -a[1])

def unitary(Ut, den):
    n = len(Ut)
    for i in range(n):
        for j in range(n):
            s1 = (0, 0); s2 = (0, 0)
            for k in range(n):
                s1 = ga(s1, gm(conj(Ut[k][i]), Ut[k][j]))   # (U*U)_ij
                s2 = ga(s2, gm(Ut[i][k], conj(Ut[j][k])))   # (UU*)_ij
            want = (den*den, 0) if i == j else (0, 0)
            if s1 != want or s2 != want:
                return False
    return True

def coefficients(Ut, den, N):
    """Coefficients of f over the scaled integers; returns g, q0 with f_jk = g_jk / q0^(j+k+1)."""
    a = Ut[0][0]; b1, b2 = Ut[0][1], Ut[0][2]
    c1, c2 = Ut[1][0], Ut[2][0]
    d11, d12, d21, d22 = Ut[1][1], Ut[1][2], Ut[2][1], Ut[2][2]
    # All over den: A=a/den etc.  I - D Delta = [[1 - d11 z, -d12 w], [-d21 z, 1 - d22 w]] (entries /den).
    # det * den^2:   Q = den^2 - den d11 z - den d22 w + (d11 d22 - d12 d21) z w
    Q = {(0, 0): (den*den, 0), (1, 0): gk(-den, d11), (0, 1): gk(-den, d22),
         (1, 1): gs(gm(d11, d22), gm(d12, d21))}
    # adj(I - D Delta) * den = [[den - d22 w, d12 w], [d21 z, den - d11 z]]
    # B Delta adj C = (b1 z, b2 w) . adj . (c1, c2)^T ; everything over den^3 after scaling.
    # f = A + [B Delta adj C]/det  =>  f * Qd = A*Qd + B Delta adj C  with Qd = den^2 det.
    # Work with numerator over den^3:  P = a*Q + den^? ... build carefully below.
    # B Delta adj C (scaled by den^3, since B,C,adj each carry one 1/den):
    #  z b1 [ (den - d22 w) c1 + d12 w c2 ] + w b2 [ d21 z c1 + (den - d11 z) c2 ]
    BAC = {}
    def add(key, v):
        BAC[key] = ga(BAC.get(key, (0, 0)), v)
    add((1, 0), gk(den, gm(b1, c1)))
    add((1, 1), gm(b1, gs(gm(d12, c2), gm(d22, c1))))
    add((0, 1), gk(den, gm(b2, c2)))
    add((1, 1), gm(b2, gs(gm(d21, c1), gm(d11, c2))))
    # f = a/den + BAC/(den^3 det) = (a Q + BAC) / (den Q)   [Q = den^2 det]
    P = {}
    for key, v in Q.items():
        P[key] = gm(a, v)
    for key, v in BAC.items():
        P[key] = ga(P.get(key, (0, 0)), v)
    # f * (den Q) = P.  Let Qs = den*Q, q0 = Qs[0,0] = den^3.
    Qs = {k: gk(den, v) for k, v in Q.items()}
    q0 = Qs[(0, 0)][0]
    assert Qs[(0, 0)][1] == 0
    g = {}
    qpow = [1]
    for _ in range(2*N + 3):
        qpow.append(qpow[-1]*q0)
    for s in range(2*N + 1):
        for j in range(max(0, s - N), min(N, s) + 1):
            k = s - j
            acc = gk(qpow[s], P.get((j, k), (0, 0)))
            for (u, v), q in Qs.items():
                if (u, v) == (0, 0) or u > j or v > k:
                    continue
                acc = gs(acc, gk(qpow[u + v - 1], gm(q, g[(j - u, k - v)])))
            g[(j, k)] = acc   # g_jk = f_jk * q0^(j+k+1)
    return g, q0, P, Qs

def majorant(g, q0, R, E, N, dps=200):
    iv.dps = dps
    r = iv.mpf(R) / iv.mpf(E)
    q = iv.mpf(q0)
    tot = iv.mpf(0)
    for (j, k), (x, y) in g.items():
        if j > N or k > N:
            continue
        m = iv.sqrt(iv.mpf(x*x + y*y)) / q**(j + k + 1)
        tot += m * r**(j + k)
    return tot

def closed_form_check(Ut, den, g, q0, N):
    mp.dps = 120
    U = [[mpc(e[0], e[1]) / den for e in row] for row in Ut]
    for z, w in [(mpf("0.11"), mpf("-0.07")), (mpc("0.05", "0.03"), mpc("-0.02", "0.09"))]:
        M = [[1 - U[1][1]*z, -U[1][2]*w], [-U[2][1]*z, 1 - U[2][2]*w]]
        det = M[0][0]*M[1][1] - M[0][1]*M[1][0]
        inv = [[M[1][1]/det, -M[0][1]/det], [-M[1][0]/det, M[0][0]/det]]
        x1 = inv[0][0]*U[1][0] + inv[0][1]*U[2][0]
        x2 = inv[1][0]*U[1][0] + inv[1][1]*U[2][0]
        f = U[0][0] + U[0][1]*z*x1 + U[0][2]*w*x2
        s = mpc(0)
        for (j, k), (x, y) in g.items():
            s += mpc(x, y) / mpf(q0)**(j + k + 1) * z**j * w**k
        err = abs(s - f)
        print(f"  closed form vs series at {(complex(z), complex(w))}: |diff| = {mp.nstr(err, 3)}")
        if err > mpf(10)**-40:
            return False
    return True

def verdict(label, Ut, den, R, E, N, margin_num=0, margin_den=1):
    if not unitary(Ut, den):
        print(f"[{label}] unitarity FAILS -> certificate rejected")
        return False
    g, q0, _, _ = coefficients(Ut, den, N)
    tot = majorant(g, q0, R, E, N)
    thr = 1 + iv.mpf(margin_num) / iv.mpf(margin_den)
    lo = tot.a
    ok = lo > thr.b
    print(f"[{label}] N={N} r={R}/{E}: sum - 1 in [{mp.nstr(mpf((tot - 1).a.a), 12)}, {mp.nstr(mpf((tot - 1).b.b), 12)}];"
          f" need > 1 + {margin_num}/{margin_den}: {'PASS' if ok else 'FAIL'}")
    return ok

def main(path):
    t0 = time.time()
    sha, d, den, Ut, R, E = load(path)
    print("sha256", sha, "OK" if sha == EXPECTED_SHA else "MISMATCH")
    N = d["N"]
    print("den digits:", len(str(den)), " N:", N, " r:", R, "/", E)
    print("exact unitarity U*U = UU* = I:", unitary(Ut, den))
    g, q0, _, _ = coefficients(Ut, den, N)
    print("coefficients:", len(g), "(expect 1681)")
    if not closed_form_check(Ut, den, g, q0, N):
        print("closed-form cross-check FAILED"); sys.exit(1)
    main_ok = verdict("MAIN claim", Ut, den, R, E, N)
    m_ok = verdict("MAIN, stated margin 2.2e-42", Ut, den, R, E, N, 22, 10**43)
    print("--- canaries (each must FAIL) ---")
    c1 = verdict("canary r=0.30282527949097", Ut, den, 30282527949097, 10**14, N)
    bad = [row[:] for row in Ut]; bad[1][1] = (bad[1][1][0] + 1, bad[1][1][1])
    c2 = verdict("canary Ut[1][1] + 1", bad, den, R, E, N)
    c3 = verdict("canary N=20", Ut, den, R, E, 20)
    c4 = verdict("canary margin 3e-42", Ut, den, R, E, N, 3, 10**42)
    canaries_ok = not (c1 or c2 or c3 or c4)
    print("canaries all failed as required:", canaries_ok)
    print(f"runtime {time.time() - t0:.1f}s")
    print("RESULT:", "CONFIRMED K_2 < r" if (main_ok and m_ok and canaries_ok) else "NOT CONFIRMED")

if __name__ == "__main__":
    main(sys.argv[1])
