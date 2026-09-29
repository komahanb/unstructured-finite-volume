#!/usr/bin/env python3
"""Rank of the constraint jacobian A of the unified beam Lagrangian,
w_tt + w_xxxx = 0 with clamped ends, on a small product grid of P
beamline points and N instants. Every jet component at every point is
an unknown. The time representation is backward Euler on the chain
w -> w_t -> w_tt. The space representation is either

    family   the trapezoid rule on each segment, on the chain
             w -> w_x -> w_xx -> w_xxx -> w_xxxx (open at both ends)
    fitted   each x component from the values at the five nearest
             points by the polynomial of degree four (closed)

The variants place the initial rows, the boundary rows and the physics
rows, and the optional mixed component w_xt with its defining law.

    python3 unified-residuals-rank.py [P] [N]      defaults 6, 4
"""
import math
import sys
import numpy as np

P = int(sys.argv[1]) if len(sys.argv) > 1 else 6
N = int(sys.argv[2]) if len(sys.argv) > 2 else 4
h, dx = 0.1, 0.2
X_CHAIN = [('w', 'wx'), ('wx', 'wxx'), ('wxx', 'wxxx'), ('wxxx', 'wxxxx')]


def fitted_weights(i, order):
    """Weights of the derivative of the given order at point i from
    the five nearest points, exact on polynomials of degree four."""
    lo = min(max(i - 2, 0), P - 5)
    s = np.arange(lo, lo + 5)
    V = np.vander((s - i) * dx, 5, increasing=True).T
    e = np.zeros(5)
    e[order] = math.factorial(order)
    return s, np.linalg.solve(V, e)


def system(space, placement, mixed=None):
    C = ['w', 'wt', 'wtt', 'wx', 'wxx', 'wxxx', 'wxxxx'] + (['wxt'] if mixed else [])
    n = P * N * len(C)
    idx = lambda i, k, c: ((k * P) + i) * len(C) + C.index(c)
    R = []

    def row(entries):
        r = np.zeros(n)
        for (i, k, c), v in entries:
            r[idx(i, k, c)] += v
        R.append(r)

    ends = (0, P - 1)
    # the physics: on every point, or on the points two away from each end
    for k in range(N):
        for i in range(P):
            if placement == 'physics_interior' and k > 0 and (i < 2 or i > P - 3):
                continue
            row([((i, k, 'wtt'), 1), ((i, k, 'wxxxx'), 1)])
    # the time representation on every step
    for k in range(1, N):
        for i in range(P):
            row([((i, k, 'w'), 1), ((i, k - 1, 'w'), -1), ((i, k, 'wt'), -h)])
            row([((i, k, 'wt'), 1), ((i, k - 1, 'wt'), -1), ((i, k, 'wtt'), -h)])
    # the space representation
    x_instants = range(1, N) if placement == 'initial_jet' else range(N)
    for k in x_instants:
        if space == 'family':
            for i in range(1, P):
                for c, d in X_CHAIN:
                    row([((i, k, c), 1), ((i - 1, k, c), -1), ((i, k, d), -dx / 2), ((i - 1, k, d), -dx / 2)])
        else:
            for i in range(P):
                for o, c in enumerate(['wx', 'wxx', 'wxxx', 'wxxxx'], start=1):
                    s, wts = fitted_weights(i, o)
                    row([((i, k, c), 1)] + [((j, k, 'w'), -wj) for j, wj in zip(s, wts)])
    # the mixed component: one defining law along t, along x, or both
    if mixed in ('along_t', 'both'):
        for k in range(1, N):
            for i in range(P):
                row([((i, k, 'wx'), 1), ((i, k - 1, 'wx'), -1), ((i, k, 'wxt'), -h)])
    if mixed in ('along_x', 'both'):
        for k in range(N):
            for i in range(1, P):
                row([((i, k, 'wt'), 1), ((i - 1, k, 'wt'), -1), ((i, k, 'wxt'), -dx / 2), ((i - 1, k, 'wxt'), -dx / 2)])
    # the initial rows
    for i in range(P):
        if placement == 'drop_ic_corner' and i in ends:
            row([((i, 0, 'wt'), 1)])
        else:
            row([((i, 0, 'w'), 1)])
            row([((i, 0, 'wt'), 1)])
        if placement == 'initial_jet':
            for c in ('wx', 'wxx', 'wxxx', 'wxxxx'):
                row([((i, 0, c), 1)])
        if mixed in ('along_t', 'both'):
            row([((i, 0, 'wxt'), 1)])
    # the boundary rows, clamped: w and w_x at both ends
    b_instants = range(N) if placement in ('all', 'drop_ic_corner') else range(1, N)
    for k in b_instants:
        for i in ends:
            row([((i, k, 'w'), 1)])
            row([((i, k, 'wx'), 1)])
            if mixed == 'along_x' and i == 0:
                row([((i, k, 'wxt'), 1)])
    if mixed == 'along_x':
        row([((0, 0, 'wxt'), 1)])
    return np.array(R), n


def report(label, space, placement, mixed=None):
    A, n = system(space, placement, mixed)
    s = np.linalg.svd(A, compute_uv=False)
    rank = int(np.sum(s > s[0] * max(A.shape) * np.finfo(float).eps))
    print(f'{label:44s} unknowns {n:4d}  rows {A.shape[0]:4d}  rank {rank:4d}  '
          f'smallest singular value {s[-1]:.1e}')


report('family, initial w w_t, bc at every instant', 'family', 'all')
report('family, initial w row dropped at the ends', 'family', 'drop_ic_corner')
report('family, initial jet, x rows and bc for t > 0', 'family', 'initial_jet')
report('fitted, physics on interior for t > 0', 'fitted', 'physics_interior')
report('family + w_xt defined along t', 'family', 'initial_jet', 'along_t')
report('family + w_xt defined along x', 'family', 'initial_jet', 'along_x')
report('family + w_xt defined along t and along x', 'family', 'initial_jet', 'both')
