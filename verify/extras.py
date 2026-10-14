"""Step 3: exact charge spot-checks, the clique-LP limit, and the tables quoted in PROOF.md.

  * charges recomputed in exact rational arithmetic straight from the problem's definition
    q(R, A) = rank_Q(V_R V_A^T) / dim_Q K, on 800 flat pairs;
  * the all-cliques clique-certificate LP has value exactly 40 (Gaussian) / 22 (Eisenstein):
    a feasible primal certificate and a feasible dual solution of the same value;
  * intermediate counts quoted in the proof text.
"""
import collections
import random
import sys
import time
from fractions import Fraction

import numpy as np

from geom import Check, closure, load, members, parse_kind, popcount

kind = parse_kind(sys.argv)
t0 = time.time()
chk = Check(f"extras ({kind})")
P = load(f"{kind}.pkl")
Qd = load(f"{kind}_q.pkl")
S = P["S"]
V, R, J, N = S["V"], S["R"], S["J"], P["N"]
n = len(V)
flist, frank, q = Qd["flist"], Qd["frank"], Qd["q"]
F = len(flist)
fidx = {f: i for i, f in enumerate(flist)}
cliques = P["cliques"]
kmax = max(map(popcount, cliques))
Qk = [c for c in cliques if popcount(c) == kmax]

# ---- exact charges -------------------------------------------------------------------
den = 2 if kind == "eis" else 1  # Eisenstein J has denominator 2
Jden = np.rint(J * den).astype(np.int64)
assert np.allclose(Jden, J * den)


def rref_rank(rows):
    M = [list(r) for r in rows]
    rk = 0
    for col in range(len(M[0]) if M else 0):
        piv = next((i for i in range(rk, len(M)) if M[i][col] != 0), None)
        if piv is None:
            continue
        M[rk], M[piv] = M[piv], M[rk]
        for i in range(len(M)):
            if i != rk and M[i][col] != 0:
                fac = M[i][col] / M[rk][col]
                M[i] = [a - fac * b for a, b in zip(M[i], M[rk])]
        rk += 1
    return rk


def rational_basis(flat):
    """Basis of X_R = span_Q{r, Jr : L_r in flat} (dim_Q K = 2, so r and Jr span each line)."""
    basis = []
    for t in members(flat):
        r = R[S["lines"][t][0]]
        for v in (r * den, Jden @ r):
            cand = basis + [[Fraction(int(x)) for x in v]]
            if rref_rank(cand) > len(basis):
                basis = cand
    return basis


rng = random.Random(1)
nonempty = [i for i in range(F) if frank[i] > 0]
mid = [i for i in range(F) if frank[i] in (2, 3)]
pairs = [(rng.choice(nonempty), rng.choice(nonempty)) for _ in range(400)]
pairs += [(rng.choice(mid), rng.choice(mid)) for _ in range(400)]
bases = {}
mismatch = 0
for i, j in pairs:
    bi = bases.setdefault(i, rational_basis(flist[i]))
    bj = bases.setdefault(j, rational_basis(flist[j]))
    gram = [[sum(a * b for a, b in zip(x, y)) for y in bj] for x in bi]
    if Fraction(rref_rank(gram), 2) != int(q[i, j]):
        mismatch += 1
chk.eq(f"exact rational charge vs floating-point charge, {len(pairs)} flat pairs: mismatches", mismatch, 0)

# ---- counts quoted in the proof --------------------------------------------------------
linecl = [0] * n
for qi, c in enumerate(Qk):
    for t in members(c):
        linecl[t] |= 1 << qi


def meets(S_):
    m = 0
    for t in members(S_):
        m |= linecl[t]
    return m


meet = [meets(f) for f in flist]
EXP_MEET = {"gauss": {2: {2: 256, 3: 308, 6: 320}, 3: {4: 436, 6: 492, 15: 512}},
            "eis": {2: {2: 1008, 4: 1408}, 3: {5: 1912, 12: 2376}}}[kind]
for r in (2, 3):
    d = collections.defaultdict(set)
    for i, f in enumerate(flist):
        if frank[i] == r:
            d[popcount(f)].add(popcount(meet[i]))
    got = {size: (v.pop() if len(v) == 1 else sorted(v)) for size, v in sorted(d.items())}
    chk.eq(f"rank-{r} flats: size -> number of maximum cliques met", got, EXP_MEET[r])

orth_pairs = [(a, b) for a in range(n) for b in range(a + 1, n) if not N[a, b]]
closure_size = {p: popcount(closure(V, p)) for p in orth_pairs}
chk.eq("orthogonal pairs by size of their closure",
       dict(sorted(collections.Counter(closure_size.values()).items())),
       {"gauss": {2: 360, 6: 90}, "eis": {2: 240}}[kind])
hyp = {t: sum(1 << s for s in range(n) if not N[t, s]) for t in range(n)}
chk.true("every t^perp is a flat", all(h in fidx for h in hyp.values()))

if kind == "gauss":
    r3 = [i for i in range(F) if frank[i] == 3]
    tab = collections.defaultdict(int)
    for i in r3:
        for j in r3:
            if q[i, j] == 2:
                x = popcount(flist[i] & flist[j])
                tab[x] = max(tab[x], popcount(meet[i] & meet[j]))
    chk.eq("charge-2 pairs of rank-3 flats: max h_both by |F_R ∩ F_A|", dict(sorted(tab.items())),
           {0: 368, 1: 364, 2: 384, 3: 308, 6: 384})
    f3 = {f for f, r in zip(flist, frank) if r == 2 and popcount(f) == 3}
    m3 = max(popcount(meet[i] & meet[j]) for i in range(F) for j in range(F)
             if q[i, j] == 2 and (flist[i] & flist[j]) in f3)
    chk.eq("charge-2 pairs with F_R ∩ F_A a 3-line flat: max h_both", m3, 308)
    hodd_f3 = {sum(1 for c in Qk if popcount(c & f) % 2) for f in f3}
    chk.eq("h_odd(F_3) for every 3-line rank-2 flat", hodd_f3, {236})

    def hodd_by_size(X):
        mem = members(X)
        out = {}
        for D in range(1 << len(mem)):
            Ds = sum(1 << mem[b] for b in range(len(mem)) if D >> b & 1)
            h = sum(1 for c in Qk if popcount(c & Ds) % 2)
            out[popcount(D)] = max(out.get(popcount(D), 0), h)
        return out

    agg = {}
    for f, r in zip(flist, frank):
        if r == 2 and popcount(f) == 6:
            for d, h in hodd_by_size(f).items():
                agg[d] = max(agg.get(d, 0), h)
    chk.eq("max h_odd(D) for D inside a 6-line rank-2 flat, by |D|", agg,
           {0: 0, 1: 128, 2: 256, 3: 192, 4: 160, 5: 128, 6: 128})
else:
    r3 = [i for i in range(F) if frank[i] == 3]
    tab = collections.defaultdict(int)
    for i in r3:
        for j in r3:
            if q[i, j] == 2:
                x = popcount(flist[i] & flist[j])
                tab[x] = max(tab[x], popcount(meet[i] & meet[j]))
    chk.eq("charge-2 pairs of rank-3 flats: max h_both by |F_R ∩ F_A|", dict(sorted(tab.items())),
           {0: 1178, 1: 1408, 2: 1872, 4: 1408})

# ---- clique-certificate LP: primal and dual of equal value ------------------------------
# Primal feasibility only needs flat pairs: h_both is monotone and q depends on closures.
sizes = np.array([popcount(c) for c in cliques])
cl_u64 = np.array(cliques, dtype=np.uint64)
meet_h = {t: (cl_u64 & np.uint64(hyp[t])) != 0 for t in range(n)}
if kind == "gauss":
    worst = max(popcount(meet[i] & meet[j]) - 192 * int(q[i, j]) for i in range(F) for j in range(F))
    chk.true("primal: weight 1/192 on each 12-clique is feasible on all flat pairs",
             worst <= 0, f"max (h_both - 192 q) = {worst}")
    chk.eq("primal value 640 * 12 / 192", Fraction(640 * 12, 192), 40)
    O6 = [p for p in orth_pairs if closure_size[p] == 6]
    chk.eq("charge q(t^perp, s^perp) on the 90 special orthogonal pairs",
           {int(q[fidx[hyp[a]], fidx[hyp[b]]]) for a, b in O6}, {2})
    h6 = np.zeros(len(cliques), dtype=np.int64)
    for a, b in O6:
        h6 += meet_h[a] & meet_h[b]
    chk.eq("min h6(Q) by clique size", {k: int(h6[sizes == k].min()) for k in range(1, kmax + 1)},
           {1: 9, 2: 15, 3: 21, 4: 26, 5: 30, 6: 39, 7: 42, 8: 42, 9: 46, 10: 51, 11: 54, 12: 54})
    chk.true("dual: (2/9) h6(Q) >= |Q| for all 2,068,372 cliques", bool((2 * h6 >= 9 * sizes).all()))
    chk.eq("dual value 90 * (2/9) * 2  (= all-cliques LP limit)", Fraction(2, 9) * 2 * len(O6), 40)
else:
    worst = max(24 * popcount(flist[i] & flist[j]) + popcount(meet[i] & meet[j]) - 960 * int(q[i, j])
                for i in range(F) for j in range(F))
    chk.true("primal: 1/40 per singleton + 1/960 per 7-clique is feasible on all flat pairs",
             worst <= 0, f"max (24|F_R∩F_A| + h_both - 960 q) = {worst}")
    chk.eq("primal value 40/40 + 2880*7/960", Fraction(40, 40) + Fraction(2880 * 7, 960), 22)
    chk.eq("charge q(t^perp, s^perp) on all orthogonal pairs",
           {int(q[fidx[hyp[a]], fidx[hyp[b]]]) for a, b in orth_pairs}, {2})
    hp = np.zeros(len(cliques), dtype=np.int64)
    for a, b in orth_pairs:
        hp += meet_h[a] & meet_h[b]
    mins = {k: int(hp[sizes == k].min()) for k in range(1, kmax + 1)}
    chk.eq("min h_perp(Q) by clique size (proof §2 quotes 12,42,66,84,108,132,156 -- erratum)",
           mins, {1: 12, 2: 48, 3: 81, 4: 106, 5: 125, 6: 145, 7: 156})
    chk.true("dual: h_perp(Q)/24 + 1/2 >= |Q| for all cliques", bool((hp + 12 >= 24 * sizes).all()))
    chk.eq("dual value 240*(1/24)*2 + (1/2)*4", Fraction(240 * 2, 24) + 2, 22)

chk.info("time", f"{time.time() - t0:.1f}s")
chk.done()
