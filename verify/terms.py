"""Step 2: the per-term inequality  Phi(R, A) <= C * q(R, A)  for every possible term.

    Phi(R, A) = cd*|R ∩ A| + co*( k*h_both(R, A) + h_odd(R ∩ A) )

where k is the maximum clique size, h_both counts maximum cliques meeting both R and A,
and h_odd counts maximum cliques Q with |Q ∩ R ∩ A| odd.

Reduction to flats.  q(R, A) = q(cl R, cl A), h_both is monotone, and R ∩ A lies in
cl R ∩ cl A, which is a flat of rank <= q (Lemma 3).  So for charge q >= 2 it suffices to
take every pair of flats (F_R, F_A) and maximise over all subsets D of F_R ∩ F_A.  For
charge 1, terms are orthogonally reduced (Lemma 4), so R ⊆ F_R minus W^perp and
A ⊆ F_A minus U^perp.  Charge 4 is handled by the closed-form bound from the proof.
The subset maximum over D uses a Walsh-Hadamard transform of the clique-trace histogram.

Summing over terms and using  sum_j |R_j ∩ A_j| >= n  together with Theorem 5 on each
maximum clique gives  sum_j q_j >= (cd*n + co*k*(k+1)*#Q_k) / C.
Writes cache/<kind>_q.pkl (charges of all flat pairs) for extras.py.
"""
import collections
import math
import sys
import time
from fractions import Fraction

import numpy as np

from geom import Check, flat_bases, flat_list, load, members, parse_kind, popcount, save

PARAMS = {"gauss": (1408, 51, 125440), "eis": (144, 1, 7200)}
EXPECTED = {
    "gauss": dict(
        per_line=128, num_q=640,
        charges={0: 9993, 1: 369780, 2: 2631390, 3: 2150820, 4: 1},
        hb_max={1: 192, 2: 384, 3: 512}, phi_max={1: 125440, 2: 250880, 3: 344784},
        hodd_rank2={0: 0, 1: 128, 2: 256, 3: 236, 4: 160, 5: 128, 6: 128},
        hodd_rank3={0: 0, 1: 128, 2: 256, 3: 384, 4: 364, 5: 336, 6: 320, 7: 304, 8: 336,
                    9: 368, 10: 320, 11: 292, 12: 276, 13: 256, 14: 224, 15: 192},
        bound=Fraction(2022, 49), floor=42,
    ),
    "eis": dict(
        per_line=504, num_q=2880,
        charges={0: 4113, 1: 87120, 2: 358110, 3: 146640, 4: 1},
        hb_max={1: 904, 2: 1872, 3: 2376}, phi_max={1: 6328, 2: 14400, 3: 19512},
        hodd_rank2={0: 0, 1: 504, 2: 1008, 3: 904, 4: 928},
        hodd_rank3={0: 0, 1: 504, 2: 1008, 3: 1512, 4: 1408, 5: 1432, 6: 1260, 7: 1284,
                    8: 1312, 9: 1440, 10: 1344, 11: 1248, 12: 1152},
        bound=Fraction(116, 5), floor=24,
    ),
}

kind = parse_kind(sys.argv)
E = EXPECTED[kind]
cd, co, C = PARAMS[kind]
t0 = time.time()
chk = Check(f"per-term inequality ({kind}): Phi <= {C}*q with cd={cd}, co={co}")
P = load(f"{kind}.pkl")
V = P["S"]["V"]
n = len(V)
kmax = max(map(popcount, P["cliques"]))
Qk = [c for c in P["cliques"] if popcount(c) == kmax]
NQ = len(Qk)
chk.eq(f"number of maximum cliques (size {kmax})", NQ, E["num_q"])

linecl = [0] * n
for qi, c in enumerate(Qk):
    for t in members(c):
        linecl[t] |= 1 << qi
per_line = {popcount(x) for x in linecl}
chk.eq("maximum cliques through each line", per_line, {E["per_line"]})
per_line = E["per_line"]


def meets(S_):
    m = 0
    for t in members(S_):
        m |= linecl[t]
    return m


flist, frank = flat_list(P["flats"])
F = len(flist)
B = flat_bases(V, flist, frank)

perp = []
for i in range(F):
    if frank[i] == 0:
        perp.append((1 << n) - 1)
    else:
        proj = np.abs(V.conj() @ B[i])
        perp.append(sum(1 << t for t in range(n) if proj[t].max() < 1e-7))

q = np.zeros((F, F), dtype=np.int8)
gap_lo, gap_hi = np.inf, 0.0
for i in range(F):
    s = np.linalg.svd(np.einsum("ab,jac->jbc", B[i].conj(), B), compute_uv=False)
    q[i] = (s > 1e-7).sum(axis=1)
    if (s > 1e-7).any():
        gap_lo = min(gap_lo, s[s > 1e-7].min())
    if (s <= 1e-7).any():
        gap_hi = max(gap_hi, s[s <= 1e-7].max())
chk.true("charge rank decisions are numerically unambiguous", gap_hi < 1e-12 and gap_lo > 1e-2,
         f"largest 'zero' singular value {gap_hi:.1e}, smallest nonzero {gap_lo:.3f}")
chk.eq("charges over all flat pairs", dict(sorted(collections.Counter(q.ravel().tolist()).items())), E["charges"])


def subset_table(X):
    """For D ⊆ X: (max of cd|D| + co*h_odd(D), {|D|: max h_odd(D)})."""
    mem = members(X)
    m = len(mem)
    hist = np.zeros(1 << m, dtype=np.int64)
    for c in Qk:
        hist[sum(1 << b for b, t in enumerate(mem) if (c >> t) & 1)] += 1
    W = hist.copy()
    h = 1
    while h < (1 << m):
        W = W.reshape(-1, 2, h)
        W = np.stack([W[:, 0] + W[:, 1], W[:, 0] - W[:, 1]], axis=1).reshape(-1)
        h *= 2
    hodd = (NQ - W) // 2
    sizes = np.array([popcount(d) for d in range(1 << m)])
    val = cd * sizes + co * hodd
    return int(val.max()), {d: int(hodd[sizes == d].max()) for d in range(m + 1)}


rank_of = dict(zip(flist, frank))
table = {f: subset_table(f) for f, r in zip(flist, frank) if r <= 3}
for r, key in ((2, "hodd_rank2"), (3, "hodd_rank3")):
    agg = {}
    for f, rr in zip(flist, frank):
        if rr == r:
            for d, h in table[f][1].items():
                agg[d] = max(agg.get(d, 0), h)
    chk.eq(f"max h_odd(D) over D inside a rank-{r} flat, by |D|", agg, E[key])

meet = [meets(f) for f in flist]
hb_max = collections.defaultdict(int)
phi_max = collections.defaultdict(int)
witness = {}
violations = 0
lemma3_failures = 0
for i in range(F):
    for j in range(F):
        qq = int(q[i, j])
        if qq in (0, 4):
            continue
        X = flist[i] & flist[j]
        if rank_of.get(X, 99) > qq:
            lemma3_failures += 1
            continue
        if qq == 1:
            hb = popcount(meets(flist[i] & ~perp[j]) & meets(flist[j] & ~perp[i]))
        else:
            hb = popcount(meet[i] & meet[j])
        phi = co * kmax * hb + table[X][0]
        hb_max[qq] = max(hb_max[qq], hb)
        if phi > phi_max[qq]:
            phi_max[qq] = phi
            witness[qq] = (frank[i], popcount(flist[i]), frank[j], popcount(flist[j]), popcount(X), hb)
        if phi > C * qq:
            violations += 1
chk.eq("flat pairs where F_R ∩ F_A is not a flat of rank <= q (Lemma 3)", lemma3_failures, 0)
chk.eq("max h_both by charge", dict(hb_max), E["hb_max"])
chk.eq("max Phi by charge", dict(phi_max), E["phi_max"])
for qq in (1, 2, 3):
    chk.info(f"q={qq} max Phi / ({C}*q)", f"{phi_max[qq] / (C * qq):.5f}  witness (rk F_R, |F_R|, rk F_A, |F_A|, |F_R∩F_A|, h_both) = {witness[qq]}")
chk.eq("flat pairs violating Phi <= C*q (charges 1-3)", violations, 0)


def hodd_cap(d):
    cap = min(NQ, per_line * d)
    if kmax % 2 == 0:
        cap = min(cap, per_line * (n - d))
    return cap


phi4 = max(cd * d + co * (kmax * NQ + hodd_cap(d)) for d in range(n + 1))
chk.true("charge 4: closed-form max Phi <= 4C", phi4 <= 4 * C, f"{phi4} <= {4 * C}")

bound = Fraction(cd * n + co * kmax * (kmax + 1) * NQ, C)
chk.eq("lower bound  sum_j q_j >= (cd*n + co*k(k+1)*#Q)/C", bound, E["bound"])
chk.eq("integer floor  C >= ceil(bound)", math.ceil(bound), E["floor"])

save(f"{kind}_q.pkl", dict(q=q, flist=flist, frank=frank, perp=perp))
chk.info("time", f"{time.time() - t0:.1f}s")
chk.done()
