# c59: independent replay of the [T2026] certificate for K_2 < 0.302825279490988

Replays the certificate behind constant 59a of [teorth/optimizationproblems](https://github.com/teorth/optimizationproblems) (L. Tang, `gomyway/c59-bohr-bidisc-certificate`).
Written from the 59a page and the certificate JSON only; the claimant's checker was not read or run, and neither was A. Canli's replay code.

## Run

```
python3 replay_c59_t2026.py c59_cert.json
```

Python 3 and `mpmath` (tested with 1.3.0). About 20 s on 4 vCPU (longer on a busy machine). Note: do not add `-I` unless mpmath is installed system-wide; `-I` hides user-site packages.
Expected last line: `RESULT: CONFIRMED K_2 < r`. Output should match `replay_c59_t2026.log` except the `runtime` line.

## Files and sha256

| file | sha256 |
|---|---|
| `replay_c59_t2026.py` | `85713037b56ecc59db54cfa2e06e8eca4588cae3d5627a6a79adc39b4538266b` |
| `replay_c59_t2026.log` | `d28404a0794d790bfa5337bd0e8b09f4ca22afefd5ed63207cd66b90d386995d` |
| `c59_cert.json` | `6d2e39e8741fb5e73593d4db007c68831243d32fe4ae298bc4fe750033ebedab` |

`c59_cert.json` is a byte copy of `c59_cert.json` at `gomyway/c59-bohr-bidisc-certificate@a014dc9` (fetched from GitHub; its sha256 equals the one printed on the 59a page). The script prints `OK` or `MISMATCH` for the input sha on its first line.

## What it checks

1. `U = Ut/den` is exactly unitary (Gaussian-integer arithmetic, no floats).
2. `f = A + B Δ (I − DΔ)⁻¹ C = P/Q`; all 1681 coefficients `f_jk` (j,k ≤ 40) by exact series division, cross-checked against the closed form at two points (diff ~1e-70).
3. `Σ |f_jk| r^(j+k)` enclosed with mpmath interval arithmetic (200 digits) at `r = 3028252794909880/10^16`; only the lower endpoint decides. Result: sum − 1 = 2.22594666642e-42 > 0, hence `K_2 < r`.
4. Four controls that must fail: smaller r, one unit changed in `Ut`, truncation at N = 20, an inflated margin (3e-42).

## Prepared with AI assistance

Code and write-up by Cadence (Claude, Anthropic). Reviewed by: Jiggy.
No affiliation with the claimant or with A. Canli.
