"""Reconstruction of the complex E8 line systems from the E8 roots.

`build("gauss")` gives the 60-line Gaussian system (K = Q(i)); `build("eis")` gives the
40-line Eisenstein system (K = Q(omega)).  Each K-line is represented by a vector in C^4
via a C-linear isometry Q^8 (x) R -> C^4 under which J acts as a scalar, so the real part
of the Hermitian form is the Euclidean dot product and K-ranks equal C-ranks.
"""
import itertools
import os
import pickle
import random

import numpy as np

TOL = 1e-8
CACHE_DIR = os.environ.get(
    "E8_CACHE", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "cache"))
KINDS = ("gauss", "eis")


def cache_path(name):
    os.makedirs(CACHE_DIR, exist_ok=True)
    return os.path.join(CACHE_DIR, name)


def save(name, obj):
    with open(cache_path(name), "wb") as fh:
        pickle.dump(obj, fh)


def load(name):
    path = cache_path(name)
    if not os.path.exists(path):
        raise SystemExit(f"missing {path}: run the earlier steps first (see README)")
    with open(path, "rb") as fh:
        return pickle.load(fh)


def parse_kind(argv):
    kind = argv[1] if len(argv) > 1 else "gauss"
    if kind not in KINDS:
        raise SystemExit(f"usage: {argv[0]} [gauss|eis]")
    return kind


def e8_roots():
    """The 240 roots in doubled coordinates (squared norm 8)."""
    roots = []
    for i, j in itertools.combinations(range(8), 2):
        for si in (2, -2):
            for sj in (2, -2):
                v = [0] * 8
                v[i], v[j] = si, sj
                roots.append(v)
    for signs in itertools.product((1, -1), repeat=8):
        if signs.count(-1) % 2 == 0:
            roots.append(list(signs))
    R = np.array(roots, dtype=np.int64)
    assert R.shape == (240, 8) and (np.sum(R * R, axis=1) == 8).all()
    return R


def gaussian_J():
    """(x, y) -> (-y, x) on each consecutive coordinate pair."""
    J = np.zeros((8, 8))
    for k in range(4):
        J[2 * k + 1, 2 * k] = 1
        J[2 * k, 2 * k + 1] = -1
    return J


def reflection(a):
    a = np.asarray(a, dtype=float)
    return np.eye(8) - np.outer(a, a) / 4.0


def find_four_orthogonal_A2(R, seed=None):
    """Pairs (a, b) with a.b = -4 spanning four mutually orthogonal A2 subsystems.

    With `seed`, the search order and the orientation of each plane are randomised.
    """
    G = R @ R.T
    rng = random.Random(seed) if seed is not None else None
    order = list(range(240))
    if rng:
        rng.shuffle(order)

    def search(chosen):
        if len(chosen) == 4:
            return chosen
        used = [c for pair in chosen for c in pair]
        for a in order:
            if any(G[a, c] != 0 for c in used):
                continue
            for b in order:
                if G[a, b] == -4 and all(G[b, c] == 0 for c in used):
                    res = search(chosen + [(a, b)])
                    if res:
                        return res
            return None
        return None

    planes = search([])
    if rng:
        planes = [(b, a) if rng.random() < 0.5 else (a, b) for a, b in planes]
    return planes


def eisenstein_J(R, planes):
    J = np.eye(8)
    for a, b in planes:
        J = J @ reflection(R[a]) @ reflection(R[b])
    return J


def build(kind):
    R = e8_roots()
    if kind == "gauss":
        J, planes = gaussian_J(), None
    else:
        seed = os.environ.get("EIS_SEED")
        planes = find_four_orthogonal_A2(R, int(seed) if seed else None)
        J = eisenstein_J(R, planes)
    I = np.eye(8)
    if kind == "gauss":
        assert np.allclose(J @ J, -I)
    else:
        assert np.allclose(J @ J + J + I, 0)
    assert np.allclose(J @ J.T, I)

    def apply_J(w):
        u = J @ w
        ui = np.rint(u).astype(np.int64)
        assert np.allclose(u, ui)
        return ui

    rset = {tuple(r): i for i, r in enumerate(R)}
    assert all(tuple(apply_J(r)) in rset for r in R), "J must preserve the roots"

    # K-lines = orbits of <J, -1> on the roots
    line_of = [-1] * 240
    lines = []
    for i in range(240):
        if line_of[i] >= 0:
            continue
        orb = {tuple(R[i])}
        frontier = [R[i]]
        while frontier:
            w = frontier.pop()
            for u in (apply_J(w), -w):
                if tuple(u) not in orb:
                    orb.add(tuple(u))
                    frontier.append(u)
        idx = sorted(rset[o] for o in orb)
        for k in idx:
            line_of[k] = len(lines)
        lines.append(idx)

    F = np.zeros((4, 8), dtype=complex)
    if kind == "gauss":
        for k in range(4):
            F[k, 2 * k] = 1
            F[k, 2 * k + 1] = 1j
    else:
        for k, (a, _) in enumerate(planes):
            f1 = R[a] / np.linalg.norm(R[a])
            Ja = J @ R[a]
            g = Ja - (Ja @ f1) * f1
            F[k] = f1 + 1j * g / np.linalg.norm(g)
    Z = (F @ R.T).T
    ZJ = (F @ (J @ R.T)).T
    k0 = np.argmax(np.abs(Z[0]))
    unit = ZJ[0, k0] / Z[0, k0]
    assert np.allclose(ZJ, unit * Z), "J must act as a scalar"
    assert np.allclose((Z @ Z.conj().T).real, R @ R.T)

    V = np.array([Z[l[0]] for l in lines])
    return dict(kind=kind, R=R, J=J, lines=lines, V=V, planes=planes)


def nonorth_matrix(V):
    return np.abs(V @ V.conj().T) > TOL


def mask_of(idx):
    m = 0
    for i in idx:
        m |= 1 << int(i)
    return m


def members(m):
    out = []
    i = 0
    while m:
        if m & 1:
            out.append(i)
        m >>= 1
        i += 1
    return out


def popcount(m):
    return bin(m).count("1")


def closure(V, idx):
    """cl(S) = lines contained in span_K(S), as a bitmask."""
    idx = list(idx)
    if not idx:
        return 0
    U, s, _ = np.linalg.svd(V[idx].T, full_matrices=False)
    Q = U[:, : int((s > TOL * s[0]).sum())]
    res = V.T - Q @ (Q.conj().T @ V.T)
    return mask_of(np.nonzero(np.linalg.norm(res, axis=0) < 1e-7)[0])


def all_flats(V):
    n = len(V)
    flats = {0: {0}, 1: {1 << t for t in range(n)}, 2: set(), 3: set(), 4: {(1 << n) - 1}}
    for r in (2, 3):
        for f in flats[r - 1]:
            for t in range(n):
                if not (f >> t) & 1:
                    flats[r].add(closure(V, members(f) + [t]))
    return flats


def all_cliques(adj, n):
    """All nonempty cliques as bitmasks; adj[i] = bitmask of nonorthogonal neighbours."""
    out = []

    def rec(cur, cand):
        out.append(cur)
        while cand:
            low = cand & -cand
            v = low.bit_length() - 1
            cand ^= low
            rec(cur | low, cand & adj[v])

    for v in range(n):
        rec(1 << v, adj[v] & ~((1 << (v + 1)) - 1))
    return out


def flat_list(flats):
    flist, frank = [], []
    for r in range(5):
        for f in sorted(flats[r]):
            flist.append(f)
            frank.append(r)
    return flist, frank


def flat_bases(V, flist, frank):
    """Orthonormal C^4 bases of each flat's span, zero-padded to 4 columns."""
    B = np.zeros((len(flist), 4, 4), dtype=complex)
    for i, f in enumerate(flist):
        if frank[i]:
            U, s, _ = np.linalg.svd(V[members(f)].T, full_matrices=False)
            assert int((s > TOL * s[0]).sum()) == frank[i]
            B[i, :, : frank[i]] = U[:, : frank[i]]
    return B


class Check:
    """Collects pass/fail results and exits nonzero if any expectation fails."""

    def __init__(self, title):
        self.title = title
        self.failed = []
        print(f"== {title}")

    def eq(self, label, got, want):
        ok = got == want
        print(f"  [{'PASS' if ok else 'FAIL'}] {label}: {got}" + ("" if ok else f" (expected {want})"))
        if not ok:
            self.failed.append(label)

    def true(self, label, cond, detail=""):
        print(f"  [{'PASS' if cond else 'FAIL'}] {label}" + (f": {detail}" if detail else ""))
        if not cond:
            self.failed.append(label)

    def info(self, label, value):
        print(f"  [info] {label}: {value}")

    def done(self):
        if self.failed:
            print(f"== {self.title}: {len(self.failed)} FAILED: {self.failed}")
            raise SystemExit(1)
        print(f"== {self.title}: all checks passed")
