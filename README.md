# gaussian-e8-copy-floor-42

A self-contained computational check that the **copy cost of the Gaussian E8 line system
is at least 42**: C(Gaussian E8) ≥ 42. Before this, the best known bounds were 40 ≤ C ≤ 60.
The package also confirms that **no clique certificate of any kind can prove more than 40**,
because the all-cliques LP has value exactly 40.

- The problem statement and the proof being checked are in [`PROOF.md`](PROOF.md), with
  verification notes and errata inline.
- Everything here is recomputed from the 240 E8 roots. No data files are shipped.

## Status of the two problems

| Line system | Bound in this package | Status |
| --- | --- | --- |
| Gaussian E8 (K = Q(i), 60 lines) | **42 ≤ C ≤ 60** | New lower bound, verified here. C is still not determined. |
| Eisenstein E8 (K = Q(ω), 40 lines) | 24 ≤ C (historical) | **Superseded.** CrocSwap/integer-mult-bounds PR #2017 reports the sharp value C = 40. |

The Eisenstein bound of 24 comes from the same argument (PROOF.md §5). The scripts still
check it in an optional mode (`./run_all.sh eis`), but it is kept only for the record.
PR #2017's result is not reproduced or checked by this package.

## What is checked

The proof adds a new per-clique inequality, a parity lemma over F₂ (PROOF.md, Theorem 5).
Restricted to a clique Q of size k, the decomposition must sum to the identity I_k. If any
term overlaps Q's diagonal an even number of times, then at least k + 1 terms must meet Q,
not just k. A clique certificate cannot express this. The proof combines that lemma with
the diagonal count Σ_j |R_j ∩ A_j| ≥ 60 into one linear functional per term:

    Φ(R, A) = 1408·|R∩A| + 51·(12·h_both(R, A) + h_odd(R∩A))  ≤  125440 · q(R, A)

Here h_both counts the 640 twelve-line cliques that meet both R and A, and h_odd counts
those with |Q ∩ R ∩ A| odd. Summing over terms gives

    C ≥ (1408·60 + 51·12·13·640) / 125440 = 2022/49 ≈ 41.27,   so C ≥ 42.

The scripts check:

1. **Census** (`census.py`):
   - 60 lines, each nonorthogonal to 44 others (1,320 edges).
   - **2,272 flats**: by rank 1 / 60 / 710 / 1500 / 1.
     - Rank-2 flats by size: 360 of size 2, 320 of size 3, 30 of size 6.
     - Rank-3 flats by size: 960 of size 4, 480 of size 6, 60 of size 15.
   - **2,068,372 nonempty cliques**, of which 640 are maximum (12 lines).
   - 8,704 inclusion-maximal cliques.
2. **Per-term inequality over all 2272² = 5,161,984 flat pairs** (`terms.py`):
   - Every charge is computed. Charge 1 pairs are orthogonally reduced (Lemma 4).
   - For each pair, the overlap R∩A is maximised exactly over **every** subset of
     F_R ∩ F_A, using a Walsh–Hadamard transform. The check confirms this intersection
     is always a flat of rank ≤ q (Lemma 3).
   - Charge 4 uses the closed-form bound from the proof.
   - Result: **zero violations**. The bound is tight at charges 1, 2 and 4, and the largest
     ratio at charge 3 is 0.916.
   - The final bound is computed as the exact fraction 2022/49.
3. **All-cliques clique-LP limit is exactly 40** (`extras.py`):
   - Primal: weight 1/192 on each 12-clique is feasible on every flat pair.
   - Dual: weight 2/9 on each of the 90 special hyperplane pairs (t^⊥, s^⊥) is feasible
     on all 2,068,372 cliques. Its value is 90 · (2/9) · 2 = 40.
   - The problem statement listed this limit as open in [40, 48]; it is exactly 40.
4. **Robustness and the proof's quoted numbers** (`extras.py`):
   - Charges recomputed in exact rational arithmetic, straight from the problem's
     definition rank_Q(V_R V_A^T)/2, agree with the floating-point charges on 800
     sampled flat pairs.
   - The floating-point rank decisions have a wide gap: largest "zero" singular value
     about 10⁻¹⁵, smallest nonzero 0.167.
   - Every intermediate table quoted in PROOF.md §4 matches.
5. **Theorem 5 brute force** (`thm5.py`): for k = 2 and 3, the smallest decomposition of
   I_k over F₂ that contains an even-diagonal term uses exactly k + 1 terms.

Every script asserts its expected values, prints `PASS`/`FAIL` per check, and exits
nonzero on any failure.

## Running it

Requires Python 3 and NumPy (tested with Python 3.12.3 and NumPy 1.26.4).

```bash
pip install -r requirements.txt     # numpy only
./run_all.sh                        # Gaussian (the subject of this package)
./run_all.sh eis                    # optional: superseded Eisenstein bound
```

You can also run the steps one at a time from `verify/`. Each step reads the previous
step's output from `cache/`; set `E8_CACHE` to put the cache somewhere else.

```bash
cd verify
python3 census.py gauss    # lines, flats, cliques     -> cache/gauss.pkl
python3 terms.py  gauss    # per-term inequality        -> cache/gauss_q.pkl
python3 extras.py gauss    # exact charges, LP limit, quoted tables
python3 thm5.py            # parity lemma brute force
```

For the optional Eisenstein mode, `EIS_SEED=<int>` randomises the choice and orientation
of the four orthogonal A2 planes. Every choice tried gives identical results.

## Runtime

About 50 seconds end to end on a 4-vCPU cloud VM. The steps take roughly
5 s (census), 26 s (terms), 19 s (extras) and under 1 s (thm5). Peak memory is about
330 MiB. The Eisenstein mode takes about 11 s.

## Expected output (`./run_all.sh`)

Timing lines (`[info] time`) vary from machine to machine. Everything else should match exactly.

```
== census (gauss)
  [PASS] lines: 60
  [PASS] nonorthogonal degrees: {44}
  [PASS] nonorthogonal edges: 1320
  [PASS] flats by rank: {0: 1, 1: 60, 2: 710, 3: 1500, 4: 1}
  [PASS] total flats: 2272
  [PASS] rank-2 flat sizes: {2: 360, 3: 320, 6: 30}
  [PASS] rank-3 flat sizes: {4: 960, 6: 480, 15: 60}
  [PASS] cliques by size: {1: 60, 2: 1320, 3: 14000, 4: 79200, 5: 252288, 6: 474240, 7: 554880, 8: 419040, 9: 203200, 10: 59904, 11: 9600, 12: 640}
  [PASS] total nonempty cliques: 2068372
  [PASS] inclusion-maximal cliques by size: {7: 1920, 8: 1920, 10: 2304, 11: 1920, 12: 640}
  [info] time: 5.1s
== census (gauss): all checks passed
== per-term inequality (gauss): Phi <= 125440*q with cd=1408, co=51
  [PASS] number of maximum cliques (size 12): 640
  [PASS] maximum cliques through each line: {128}
  [PASS] charge rank decisions are numerically unambiguous: largest 'zero' singular value 2.0e-15, smallest nonzero 0.167
  [PASS] charges over all flat pairs: {0: 9993, 1: 369780, 2: 2631390, 3: 2150820, 4: 1}
  [PASS] max h_odd(D) over D inside a rank-2 flat, by |D|: {0: 0, 1: 128, 2: 256, 3: 236, 4: 160, 5: 128, 6: 128}
  [PASS] max h_odd(D) over D inside a rank-3 flat, by |D|: {0: 0, 1: 128, 2: 256, 3: 384, 4: 364, 5: 336, 6: 320, 7: 304, 8: 336, 9: 368, 10: 320, 11: 292, 12: 276, 13: 256, 14: 224, 15: 192}
  [PASS] flat pairs where F_R ∩ F_A is not a flat of rank <= q (Lemma 3): 0
  [PASS] max h_both by charge: {1: 192, 2: 384, 3: 512}
  [PASS] max Phi by charge: {1: 125440, 2: 250880, 3: 344784}
  [info] q=1 max Phi / (125440*q): 1.00000  witness (rk F_R, |F_R|, rk F_A, |F_A|, |F_R∩F_A|, h_both) = (2, 6, 3, 15, 1, 192)
  [info] q=2 max Phi / (125440*q): 1.00000  witness (rk F_R, |F_R|, rk F_A, |F_A|, |F_R∩F_A|, h_both) = (3, 15, 3, 15, 6, 384)
  [info] q=3 max Phi / (125440*q): 0.91620  witness (rk F_R, |F_R|, rk F_A, |F_A|, |F_R∩F_A|, h_both) = (3, 15, 3, 15, 15, 512)
  [PASS] flat pairs violating Phi <= C*q (charges 1-3): 0
  [PASS] charge 4: closed-form max Phi <= 4C: 501760 <= 501760
  [PASS] lower bound  sum_j q_j >= (cd*n + co*k(k+1)*#Q)/C: 2022/49
  [PASS] integer floor  C >= ceil(bound): 42
  [info] time: 26.4s
== per-term inequality (gauss): Phi <= 125440*q with cd=1408, co=51: all checks passed
== extras (gauss)
  [PASS] exact rational charge vs floating-point charge, 800 flat pairs: mismatches: 0
  [PASS] rank-2 flats: size -> number of maximum cliques met: {2: 256, 3: 308, 6: 320}
  [PASS] rank-3 flats: size -> number of maximum cliques met: {4: 436, 6: 492, 15: 512}
  [PASS] orthogonal pairs by size of their closure: {2: 360, 6: 90}
  [PASS] every t^perp is a flat
  [PASS] charge-2 pairs of rank-3 flats: max h_both by |F_R ∩ F_A|: {0: 368, 1: 364, 2: 384, 3: 308, 6: 384}
  [PASS] charge-2 pairs with F_R ∩ F_A a 3-line flat: max h_both: 308
  [PASS] h_odd(F_3) for every 3-line rank-2 flat: {236}
  [PASS] max h_odd(D) for D inside a 6-line rank-2 flat, by |D|: {0: 0, 1: 128, 2: 256, 3: 192, 4: 160, 5: 128, 6: 128}
  [PASS] primal: weight 1/192 on each 12-clique is feasible on all flat pairs: max (h_both - 192 q) = 0
  [PASS] primal value 640 * 12 / 192: 40
  [PASS] charge q(t^perp, s^perp) on the 90 special orthogonal pairs: {2}
  [PASS] min h6(Q) by clique size: {1: 9, 2: 15, 3: 21, 4: 26, 5: 30, 6: 39, 7: 42, 8: 42, 9: 46, 10: 51, 11: 54, 12: 54}
  [PASS] dual: (2/9) h6(Q) >= |Q| for all 2,068,372 cliques
  [PASS] dual value 90 * (2/9) * 2  (= all-cliques LP limit): 40
  [info] time: 18.9s
== extras (gauss): all checks passed
== Theorem 5 brute force
  [PASS] k=2: fewest terms summing to I_k with an even-diagonal term: 3
  [PASS] k=3: fewest terms summing to I_k with an even-diagonal term: 4
== Theorem 5 brute force: all checks passed
ALL CHECKS PASSED (gauss) in 51s
```

## Files

| Path | Contents |
| --- | --- |
| `verify/geom.py` | E8 roots, the complex structure J, K-lines, flats, clique enumeration, `Check` helper |
| `verify/census.py` | Step 1: census |
| `verify/terms.py` | Step 2: per-term inequality and the 2022/49 bound |
| `verify/extras.py` | Step 3: exact charges, clique-LP primal and dual, quoted tables |
| `verify/thm5.py` | Step 4: brute force of the parity lemma |
| `run_all.sh` | Runs everything from a clean cache |
| `PROOF.md` | Problem statement, the proof text, and verification notes and errata |
| `LICENSE` | Apache License 2.0 |

## Sources

The problem statement cites CrocSwap/integer-mult-bounds PRs #620 (geometry), #738 and
#770 (Gaussian clique bounds), and #798 (Eisenstein clique bounds). The proof text was
generated by Gemini; this package is an independent check of it.
