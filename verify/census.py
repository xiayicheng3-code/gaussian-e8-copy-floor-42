"""Step 1: lines, orthogonality graph, flats, and all cliques.  Writes cache/<kind>.pkl."""
import collections
import sys
import time

from geom import Check, all_cliques, all_flats, build, nonorth_matrix, parse_kind, popcount, save

EXPECTED = {
    "gauss": dict(
        lines=60, degree={44}, edges=1320,
        flats={0: 1, 1: 60, 2: 710, 3: 1500, 4: 1}, total_flats=2272,
        rank2_sizes={2: 360, 3: 320, 6: 30}, rank3_sizes={4: 960, 6: 480, 15: 60},
        cliques={1: 60, 2: 1320, 3: 14000, 4: 79200, 5: 252288, 6: 474240, 7: 554880,
                 8: 419040, 9: 203200, 10: 59904, 11: 9600, 12: 640},
        total_cliques=2068372,
        maximal={7: 1920, 8: 1920, 10: 2304, 11: 1920, 12: 640},
    ),
    "eis": dict(
        lines=40, degree={27}, edges=540,
        flats={0: 1, 1: 40, 2: 330, 3: 400, 4: 1}, total_flats=772,
        rank2_sizes={2: 240, 4: 90}, rank3_sizes={5: 360, 12: 40},
        cliques={1: 40, 2: 540, 3: 3240, 4: 9450, 5: 13824, 6: 10080, 7: 2880},
        total_cliques=40054,
        maximal={4: 90, 7: 2880},
    ),
}

kind = parse_kind(sys.argv)
E = EXPECTED[kind]
t0 = time.time()
chk = Check(f"census ({kind})")
S = build(kind)
V = S["V"]
n = len(V)
N = nonorth_matrix(V)
chk.eq("lines", n, E["lines"])
chk.eq("nonorthogonal degrees", set((N.sum(axis=1) - 1).tolist()), E["degree"])
chk.eq("nonorthogonal edges", int((N.sum() - n) // 2), E["edges"])
adj = [sum(1 << j for j in range(n) if j != i and N[i, j]) for i in range(n)]

flats = all_flats(V)
chk.eq("flats by rank", {r: len(f) for r, f in flats.items()}, E["flats"])
chk.eq("total flats", sum(len(f) for f in flats.values()), E["total_flats"])
chk.eq("rank-2 flat sizes", dict(sorted(collections.Counter(map(popcount, flats[2])).items())), E["rank2_sizes"])
chk.eq("rank-3 flat sizes", dict(sorted(collections.Counter(map(popcount, flats[3])).items())), E["rank3_sizes"])

cliques = all_cliques(adj, n)
chk.eq("cliques by size", dict(sorted(collections.Counter(map(popcount, cliques)).items())), E["cliques"])
chk.eq("total nonempty cliques", len(cliques), E["total_cliques"])
cset = set(cliques)
maximal = [c for c in cliques if not any((c | (1 << v)) in cset for v in range(n) if not (c >> v) & 1)]
chk.eq("inclusion-maximal cliques by size", dict(sorted(collections.Counter(map(popcount, maximal)).items())), E["maximal"])

save(f"{kind}.pkl", dict(S=S, adj=adj, flats=flats, cliques=cliques, N=N))
chk.info("time", f"{time.time() - t0:.1f}s")
chk.done()
