"""Smoothing factors, mu = rho (Q S): S a sweep's error operator for the
(2d+1)-point Laplacian on a periodic grid, Q the projector onto the high
frequencies (some |theta_a| >= pi / 2). Weighted Jacobi, lexicographic and
red-black Gauss-Seidel, in 2D and 3D. Prints JSON."""
import itertools
import json
import numpy as np


def laplacian(n, d):
    N = n ** d
    A = np.zeros((N, N))
    idx = np.arange(N).reshape((n,) * d)
    for p in itertools.product(range(n), repeat=d):
        i = idx[p]
        A[i, i] = -2 * d
        for a in range(d):
            for s in (-1, 1):
                q = list(p)
                q[a] = (q[a] + s) % n
                A[i, idx[tuple(q)]] += 1
    return A, idx


def high_projector(n, d):
    k = np.fft.fftfreq(n) * 2 * np.pi          # theta in [-pi, pi)
    grids = np.meshgrid(*([k] * d), indexing='ij')
    high = np.zeros((n,) * d, bool)
    for g in grids:
        high |= np.abs(g) >= np.pi / 2 - 1e-12
    N = n ** d
    # Q = F^-1 diag (high) F, F the d-dimensional DFT on flattened vectors
    F = np.ones((1, 1), complex)
    w = np.exp(-2j * np.pi * np.outer(np.arange(n), np.arange(n)) / n)
    for _ in range(d):
        F = np.kron(F, w)
    return np.linalg.inv(F) @ np.diag(high.flatten().astype(float)) @ F


def jacobi(A, d, omega):
    return np.eye(len(A)) - omega * A / (-2 * d)


def red_black(A, idx, d):
    N = len(A)
    colour = np.zeros(N, int)
    for p in itertools.product(range(idx.shape[0]), repeat=d):
        colour[idx[p]] = sum(p) % 2

    def half(c):                     # update the cells of colour c: e = mean of neighbours
        S = np.eye(N)
        for i in np.where(colour == c)[0]:
            S[i, :] = 0
            for j in np.nonzero(A[i])[0]:
                if j != i:
                    S[i, j] = A[i, j] / (2 * d)
        return S
    return half(1) @ half(0)


def lex_gs_symbol(d, n=256):
    k = (np.arange(n) / n) * 2 * np.pi - np.pi
    g = np.meshgrid(*([k] * d), indexing='ij')
    high = np.zeros(g[0].shape, bool)
    num = np.zeros(g[0].shape, complex)
    den = 2 * d * np.ones(g[0].shape, complex)
    for t in g:
        high |= np.abs(t) >= np.pi / 2
        num += np.exp(1j * t)
        den -= np.exp(-1j * t)
    return float(np.max(np.abs(num / den)[high]))


out = {}
for d, n in ((2, 16), (3, 8)):
    A, idx = laplacian(n, d)
    Q = high_projector(n, d)
    rho = lambda S: float(np.max(np.abs(np.linalg.eigvals(Q @ S))))
    omega = 2 * d / (2 * d + 1)
    out[d] = {
        'jacobi_omega': omega,
        'jacobi': rho(jacobi(A, d, omega)),
        'jacobi_1': rho(jacobi(A, d, 1.0)),
        'lex_gs': lex_gs_symbol(d, 128 if d == 2 else 48),
        'red_black': rho(red_black(A, idx, d)),
    }
print(json.dumps(out, indent=1))
