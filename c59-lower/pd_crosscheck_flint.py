"""Second, independent PD method for replay_c59_lower.py.

Replaces Sylvester/Bareiss with: exact characteristic polynomial over Q (python-flint),
then a real symmetric matrix is PD iff its charpoly coefficients strictly alternate in sign
(all roots real since symmetric; Descartes => no root <= 0). Runs the same check() path.
Run: python3 -I pd_crosscheck_flint.py CERT.json   (needs python-flint)
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flint
import replay_c59_lower as R

calls = [0]

def pos_def_charpoly(M):
    n = len(M)
    if not R.is_hermitian(M):
        return False
    N = 2 * n
    A = flint.fmpq_mat(N, N)
    for i in range(n):
        for j in range(n):
            x, y = M[i][j]
            for (a, b, v) in ((i, j, x), (i + n, j + n, x), (i, j + n, -y), (i + n, j, y)):
                A[a, b] = flint.fmpq(v.numerator, v.denominator)
    c = A.charpoly().coeffs()  # low -> high, monic
    calls[0] += 1
    # p(t) = sum c_k t^k ; PD iff (-1)^(N-k) c_k > 0 for all k
    return all((c[k] > 0) if (N - k) % 2 == 0 else (c[k] < 0) for k in range(N + 1))

R.pos_def_exact = pos_def_charpoly
sys.argv = [sys.argv[0]] + sys.argv[1:]
R.main()
print("charpoly PD checks performed:", calls[0])
