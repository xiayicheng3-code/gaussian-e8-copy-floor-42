"""Step 4: brute-force sanity check of Theorem 5 (the F_2 parity lemma) for k = 2, 3.

If I_k = sum of rank-1 terms u v^T over F_2 and some term has v.u = 0 (even diagonal
overlap), then at least k + 1 terms are needed; k + 1 is attained.  Repeated terms cancel
in pairs over F_2, so it suffices to search sets of distinct terms.
"""
import itertools

import numpy as np

from geom import Check

chk = Check("Theorem 5 brute force")
for k in (2, 3):
    vecs = [np.array(v) for v in itertools.product((0, 1), repeat=k) if any(v)]
    terms = [(np.outer(u, v) % 2, int(u @ v) % 2) for u in vecs for v in vecs]
    eye = np.eye(k, dtype=int)
    smallest = None
    for m in range(1, k + 2):
        for combo in itertools.combinations(range(len(terms)), m):
            if all(terms[c][1] == 1 for c in combo):
                continue
            if (sum(terms[c][0] for c in combo) % 2 == eye).all():
                smallest = m
                break
        if smallest:
            break
    chk.eq(f"k={k}: fewest terms summing to I_k with an even-diagonal term", smallest, k + 1)
chk.done()
