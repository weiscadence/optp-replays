#!/usr/bin/env python3
"""Independent exact replay of the C59 LOWER-bound certificate (K_2 >= 189/625), PR #403.

Own code; the author's verify_c59_lower.py was neither read nor run.
Uses only the JSON data + the argument in the PR/README.  Pure stdlib, exact arithmetic.

Run:  nice -n 19 python3 -I replay_c59_lower.py CERT.json [--r NUM/DEN] [--claim DECIMAL] [--no-tail]
      nice -n 19 python3 -I replay_c59_lower.py --canaries CERT.json

Index set: monomials z1^a z2^b with a+b <= K, degree-major, first exponent descending
  (0,0),(1,0),(0,1),(2,0),(1,1),(0,2),...   (21 for K=5).
Z (21x21 Hermitian): Z00 = -B/2, Z_{a0} = c_a u_a / 2 (c_a = r^|a|), rest = free block.
  Z0  : first column c_a/2 at a in {(1,0),(0,1)} (u=1), block = Z0_block.
  Z_a : first column e_a/2, first row 0, block = Za[a]   (non-Hermitian; enters as c_a(u Z_a + conj(u) Z_a^*)).
D_i(Z)[a,b] = Z[a+e_i, b+e_i] (0 if outside the index set) - Z[a,b]   [= J^T Z J - Z, J: e_{a-e_i} -> e_a].
Checks:
  (1) 36 cube blocks [[P, c D_i(Z_a)],[c D_i(Z_a)^*, Q]]  positive definite
  (2) 2 main LMIs  D_i(Z0) - sum_a (P_{a,i}+Q_{a,i})  positive definite
  (3) B <= 0.99315605 (rational), r = 189/625
  (4) tail  sum_{n>=6} 2 sqrt(n+1) r^n  <=  rigorous rational upper bound
  (5) B + tail < 1
PD test is exact: Hermitian M = X+iY is embedded as real symmetric [[X,-Y],[Y,X]]
(PD iff M PD), scaled to integers, and all leading principal minors are computed by
fraction-free Bareiss elimination; M is PD iff all are > 0 (Sylvester).
"""
import sys, json, math
from fractions import Fraction as F


def fr(p):
    return F(p[0], p[1])


def cx(e):
    return (fr(e[0]), fr(e[1]))


def load_mat(m, n):
    M = [[(F(0), F(0)) for _ in range(n)] for _ in range(n)]
    for k, v in m.items():
        i, j = map(int, k.split(','))
        M[i][j] = cx(v)
    return M


def cadd(a, b): return (a[0] + b[0], a[1] + b[1])
def csub(a, b): return (a[0] - b[0], a[1] - b[1])
def cscale(a, s): return (a[0] * s, a[1] * s)
def cconj(a): return (a[0], -a[1])
ZERO = (F(0), F(0))


def monos(K):
    out = []
    for d in range(K + 1):
        for a in range(d, -1, -1):
            out.append((a, d - a))
    return out


def Dop(Z, idx, pos, i):
    """D_i(Z) = J_i^T Z J_i - Z,  (J_i^T Z J_i)[a,b] = Z[a+e_i, b+e_i] if both in index set else 0."""
    n = len(idx)
    R = [[ZERO] * n for _ in range(n)]
    for a in range(n):
        for b in range(n):
            ma, mb = idx[a], idx[b]
            sa = list(ma); sa[i] += 1
            sb = list(mb); sb[i] += 1
            sa, sb = tuple(sa), tuple(sb)
            sh = Z[pos[sa]][pos[sb]] if (sa in pos and sb in pos) else ZERO
            R[a][b] = csub(sh, Z[a][b])
    return R


def is_hermitian(M):
    n = len(M)
    return all(M[i][j] == cconj(M[j][i]) for i in range(n) for j in range(n))


def lcm(a, b): return a * b // math.gcd(a, b)


def pos_def_exact(M):
    """M complex Hermitian (list of (re,im) Fractions). True iff positive definite."""
    n = len(M)
    if not is_hermitian(M):
        return False
    N = 2 * n
    R = [[F(0)] * N for _ in range(N)]
    for i in range(n):
        for j in range(n):
            x, y = M[i][j]
            R[i][j] = x; R[i + n][j + n] = x
            R[i][j + n] = -y; R[i + n][j] = y
    den = 1
    for row in R:
        for v in row:
            den = lcm(den, v.denominator)
    A = [[int(v * den) for v in row] for row in R]
    # Bareiss, no pivoting: pivots are leading principal minors
    prev = 1
    for k in range(N):
        piv = A[k][k]
        if piv <= 0:
            return False
        for i in range(k + 1, N):
            for j in range(k + 1, N):
                A[i][j] = (A[i][j] * piv - A[i][k] * A[k][j]) // prev
        prev = piv
    return True


def sqrt_upper(m, k=30):
    """rational upper bound of sqrt(m)."""
    s = math.isqrt(m * 10 ** (2 * k))
    if s * s < m * 10 ** (2 * k):
        s += 1
    return F(s, 10 ** k)


def tail_upper(r, n0=6, N0=200):
    """rigorous upper bound of sum_{n>=n0} 2 sqrt(n+1) r^n."""
    tot = F(0)
    rn = r ** n0
    for n in range(n0, N0 + 1):
        tot += 2 * sqrt_upper(n + 1) * rn
        rn *= r
    # n > N0: t_{n+1}/t_n = r sqrt((n+2)/(n+1)) <= r (1 + 1/(2(N0+2)))  for n>=N0+1
    t_next = 2 * sqrt_upper(N0 + 2) * rn  # n = N0+1
    rho = r * (1 + F(1, 2 * (N0 + 2)))
    assert rho < 1
    tot += t_next / (1 - rho)
    return tot


def check(cert, r_override=None, claim=None, do_tail=True, verbose=True):
    K = cert['K']
    r = fr(cert['r']) if r_override is None else r_override
    B = fr(cert['B'])
    idx = monos(K)
    pos = {m: i for i, m in enumerate(idx)}
    n = len(idx)
    c = {m: r ** (m[0] + m[1]) for m in idx}
    ok = True
    log = []

    def rep(name, val):
        nonlocal ok
        log.append((name, val))
        if not val:
            ok = False
        if verbose:
            print(f'  {name}: {"PASS" if val else "FAIL"}')

    # Z0
    Z0 = [[ZERO] * n for _ in range(n)]
    Z0[0][0] = (-B / 2, F(0))
    for m in [(1, 0), (0, 1)]:
        Z0[pos[m]][0] = (c[m] / 2, F(0))
        Z0[0][pos[m]] = (c[m] / 2, F(0))
    blk = load_mat(cert['Z0_block'], n - 1)
    for i in range(n - 1):
        for j in range(n - 1):
            Z0[i + 1][j + 1] = blk[i][j]
    rep('Z0 Hermitian', is_hermitian(Z0))

    robust = [m for m in idx if m[0] + m[1] >= 2]
    Za = {}
    for m in robust:
        Zm = [[ZERO] * n for _ in range(n)]
        Zm[pos[m]][0] = (F(1, 2), F(0))
        b = load_mat(cert['Za'][str(m)], n - 1)
        for i in range(n - 1):
            for j in range(n - 1):
                Zm[i + 1][j + 1] = b[i][j]
        Za[m] = Zm

    # (1) cube blocks
    nb = 0
    nb_ok = 0
    herm_ok = True
    for i in range(2):
        for m in robust:
            P = load_mat(cert['P'][f'{m}|{i}'], n)
            Q = load_mat(cert['Q'][f'{m}|{i}'], n)
            herm_ok &= is_hermitian(P) and is_hermitian(Q)
            D = Dop(Za[m], idx, pos, i)
            top = [P[a] + [cscale(D[a][b], c[m]) for b in range(n)] for a in range(n)]
            bot = [[cscale(cconj(D[b][a]), c[m]) for b in range(n)] + Q[a] for a in range(n)]
            M = [top[a] for a in range(n)] + [bot[a] for a in range(n)]
            nb += 1
            res = pos_def_exact(M)
            nb_ok += res
            if verbose and not res:
                print(f'    block {m}|{i} NOT PD')
    rep(f'P,Q Hermitian', herm_ok)
    rep(f'cube blocks PD ({nb_ok}/{nb})', nb_ok == nb and nb == 36)

    # (2) main LMIs
    main_ok = True
    for i in range(2):
        D0 = Dop(Z0, idx, pos, i)
        S = [[D0[a][b] for b in range(n)] for a in range(n)]
        for m in robust:
            P = load_mat(cert['P'][f'{m}|{i}'], n)
            Q = load_mat(cert['Q'][f'{m}|{i}'], n)
            for a in range(n):
                for b in range(n):
                    S[a][b] = csub(S[a][b], cadd(P[a][b], Q[a][b]))
        res = pos_def_exact(S)
        main_ok &= res
        if verbose:
            print(f'    main LMI i={i}: {"PD" if res else "NOT PD"}')
    rep('2 main LMIs PD', main_ok)

    # (3) B
    Bcap = F(99315605, 10 ** 8)
    Bcap_hi = F(99315606, 10 ** 8)
    rep(f'B = {float(B):.12f} <= 0.99315606 (cert states 0.99315605..)', B <= Bcap_hi)
    # (4)(5)
    if do_tail:
        tail = tail_upper(r)
        rep(f'tail <= {float(tail):.10f} (<= 0.0059699)', tail <= F(59699, 10 ** 7))
        total = B + tail
        rep(f'B + tail = {float(total):.10f} < 1', total < 1)
        if claim is not None:
            rep(f'claimed total {claim} >= B+tail', F(claim) >= total)
    return ok, log


def main():
    args = sys.argv[1:]
    if not args:
        print(__doc__); sys.exit(2)
    if args[0] == '--canaries':
        sys.exit(canaries(args[1]))
    path = args[0]
    r_over = None
    claim = None
    do_tail = True
    i = 1
    while i < len(args):
        if args[i] == '--r':
            p, q = args[i + 1].split('/'); r_over = F(int(p), int(q)); i += 2
        elif args[i] == '--claim':
            claim = args[i + 1]; i += 2
        elif args[i] == '--no-tail':
            do_tail = False; i += 1
        else:
            i += 1
    cert = json.load(open(path))
    print(f'K={cert["K"]} r={cert["r"]} B={float(fr(cert["B"])):.12f}')
    ok, _ = check(cert, r_over, claim, do_tail)
    print('RESULT:', 'CONFIRMED  sum_{n<=5} r^n|H_n|_1 <= B and B+tail<1 => K_2 >= r' if ok else 'NOT CONFIRMED')
    sys.exit(0 if ok else 1)


def canaries(path):
    import copy
    cert = json.load(open(path))
    results = []

    def run(name, c, **kw):
        print(f'[canary] {name}')
        ok, _ = check(c, verbose=True, **kw)
        results.append((name, ok))
        print(f'   -> {"ACCEPTED (BAD)" if ok else "rejected (good)"}')

    # baseline
    print('[baseline]')
    okb, _ = check(cert)
    # 1. larger r: 0.3028 (just under rho*): cert matrices were built for 189/625
    run('r = 0.3028 (cert built for 0.3024)', cert, r_override=F(3028, 10000))
    run('r = 0.31', cert, r_override=F(31, 100))
    # 2. perturb a P entry (diag) to make a cube block indefinite
    c2 = copy.deepcopy(cert)
    k = '(3, 1)|0'
    e = c2['P'][k]['5,5']
    e[0][0] = -abs(e[0][0]) * 10 ** 6 - 10 ** 12  # large negative diagonal
    run('P[(3,1)|0][5,5] made hugely negative', c2)
    # 2b. subtle perturbation in Z0 block diagonal
    c2b = copy.deepcopy(cert)
    d = c2b['Z0_block']['3,3']
    d[0][0] = d[0][0] - 10 * d[0][1]  # shifts entry by -10
    run('Z0_block[3,3] shifted by -10', c2b)
    # 3. tail matters: it consumes most of the margin.  With B inflated just enough that B+tail >= 1 the
    #    final inequality must fail (checked directly; LMIs are not the point of this control).
    r = fr(cert['r']); B = fr(cert['B'])
    tail = tail_upper(r)
    print(f'[canary] tail share: 1-B = {float(1 - B):.6f}, tail = {float(tail):.6f}, remaining margin = {float(1 - B - tail):.6f}')
    for name, Bp in [('B inflated to 0.9941', F(9941, 10000))]:
        ok_ = (Bp + tail) < 1
        print(f'[canary] {name}: B+tail = {float(Bp + tail):.6f} -> {"ACCEPTED (BAD)" if ok_ else "rejected (good)"}')
        results.append((name, ok_))
    # 4. inflated claimed bound: claim total 0.99 < true B+tail => must be refused
    run('claimed total 0.9930 (< true B+tail)', cert, claim='0.9930')
    # 5. B deflated (claim B too small): LMI no longer PD
    c5 = copy.deepcopy(cert)
    c5['B'] = [c5['B'][0] - c5['B'][1] // 100, c5['B'][1]]  # B*0.99
    run('B deflated by 1% (claims a better bound)', c5)
    bad = [n for n, ok in results if ok]
    print('baseline ok:', okb, '| canaries all rejected:', not bad)
    return 0 if (okb and not bad) else 1


if __name__ == '__main__':
    main()
