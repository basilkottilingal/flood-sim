"""Figures for the red-black Gauss-Seidel note. Writes fig*.png here."""
import json
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch

RED, BLACK = '#c0392b', '#2c3e50'
RED_L, BLACK_L = '#f2c4bf', '#b9c3cc'
GREEN, ORANGE, GREY = '#27ae60', '#e67e22', '#bdc3c7'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9,
                     'axes.spines.top': False, 'axes.spines.right': False})
DPI = 200
lfa = json.load(open('lfa.json'))
cost = {2: json.load(open('cost2.json')), 3: json.load(open('cost3.json'))}


def cell(ax, x, y, colour, ghost=False, edge='white', lw=0.6, alpha=1, hatch=None):
    ax.add_patch(Rectangle((x, y), 1, 1, facecolor=colour, edgecolor=edge,
                           lw=lw, alpha=alpha, hatch=hatch))


# 1. the colouring of a block, a red cell's stencil -------------------------
def fig_checkerboard():
    m, g = 8, 2
    fig, ax = plt.subplots(figsize=(5.2, 5.2))
    for j in range(-g, m + g):
        for i in range(-g, m + g):
            inner = 0 <= i < m and 0 <= j < m
            red = (i + j) % 2 == 0
            c = (RED if red else BLACK) if inner else (RED_L if red else BLACK_L)
            cell(ax, i, j, c)
    ax.add_patch(Rectangle((0, 0), m, m, fill=False, edgecolor='black', lw=1.6))
    p = (3, 3)
    cx, cy = p[0] + 0.5, p[1] + 0.5
    for d in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        ax.add_patch(FancyArrowPatch((cx, cy), (cx + 0.82 * d[0], cy + 0.82 * d[1]),
                                     arrowstyle='-|>', mutation_scale=11, color='white', lw=1.6))
    ax.add_patch(Rectangle(p, 1, 1, fill=False, edgecolor='yellow', lw=2.2))
    for q, lab in (((4, 4), 's[1,1]'), ((5, 3), 's[2]')):
        ax.add_patch(Rectangle(q, 1, 1, fill=False, edgecolor='yellow', lw=1.6, ls='--'))
        ax.text(q[0] + 0.5, q[1] + 1.15, lab, ha='center', va='bottom', color='black', fontsize=8,
                bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none'))
    ax.text(cx, p[1] - 0.25, 'the point', ha='center', va='top', fontsize=8,
            bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none'))
    ax.set_xlim(-g - 0.2, m + g + 0.2)
    ax.set_ylim(-g - 0.2, m + g + 0.2)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title('A block (m = 8 shown; m >= 16 in libb-amr), 2 ghost layers\n'
                 'interior dark, ghosts light; red: (i + j) even', fontsize=9)
    fig.savefig('fig1_checkerboard.png', dpi=DPI, bbox_inches='tight')
    plt.close(fig)


# 2. parity across blocks, and at a level jump ----------------------------
def fig_blocks():
    m = 8
    fig, ax = plt.subplots(figsize=(7.6, 4.0))
    # two blocks of level l, side by side: global i = local i + m * I
    for B in (0, 1):
        for j in range(m):
            for i in range(m):
                gi = i + m * B
                cell(ax, gi, j, RED if (gi + j) % 2 == 0 else BLACK)
        ax.add_patch(Rectangle((m * B, 0), m, m, fill=False, edgecolor='black', lw=1.6))
        ax.text(m * B + m / 2, m + 0.35, f'block I = {B} of level l\nlocal i = 0 .. {m - 1}',
                ha='center', va='bottom', fontsize=8)
    # a coarser leaf to the right: cells twice as large, its own colouring
    for j in range(m // 2):
        for i in range(m // 2):
            x, y = 2 * m + 2 * i + 0.6, 2 * j
            ax.add_patch(Rectangle((x, y), 2, 2, facecolor=RED if (i + j) % 2 == 0 else BLACK,
                                   edgecolor='white', lw=0.6, alpha=0.55))
    ax.add_patch(Rectangle((2 * m + 0.6, 0), m, m, fill=False, edgecolor='black', lw=1.6))
    ax.text(2 * m + 0.6 + m / 2, m + 0.35, 'a coarser leaf (level l - 1)\nits own colouring',
            ha='center', va='bottom', fontsize=8)
    # the ghosts of block 1 toward the coarse leaf: prolonged, no colour
    for j in range(m):
        for gx in range(2):
            ax.add_patch(Rectangle((2 * m + gx * 0.3, j), 0.3, 1, facecolor=GREY,
                                   edgecolor='white', lw=0.4, hatch='///'))
    ax.annotate('ghosts at the level jump:\nprolonged from the coarse leaf,\n'
                'no colour of level l', xy=(2 * m + 0.3, m / 2), xytext=(2 * m + 4.6, -2.9),
                fontsize=8, ha='center', arrowprops=dict(arrowstyle='->', lw=0.8))
    ax.text(0, -1.0, 'global parity = (i + m I + j + m J) mod 2\n= (i + j) mod 2, as m is even\n'
            '(i, j local; I, J the block\'s index)', ha='left', va='top', fontsize=8)
    ax.set_xlim(-0.5, 3 * m + 1.5)
    ax.set_ylim(-4.2, m + 2.0)
    ax.set_aspect('equal')
    ax.axis('off')
    fig.savefig('fig2_blocks.png', dpi=DPI, bbox_inches='tight')
    plt.close(fig)


# 3. a sweep, to scale by the measured costs ---------------------------------
def fig_sweeps():
    fig, axes = plt.subplots(2, 1, figsize=(7.6, 3.4), sharex=False)
    for ax, d in zip(axes, (2, 3)):
        R, M, G = cost[d]['R'], cost[d]['M'], cost[d]['G']
        nj = math.log(lfa[str(d)]['red_black']) / math.log(lfa[str(d)]['jacobi'])
        rows = [('weighted Jacobi\n(%.2f sweeps)' % nj,
                 [('relax: dn', R, '#7f8c8d'), ('mix: da', M, '#95a5a6'), ('ghosts', G, ORANGE)] * 3,
                 nj),
                ('red-black GS\n(1 sweep)', [('red', R / 2, RED), ('ghosts', G, ORANGE),
                                            ('black', R / 2, BLACK), ('ghosts', G, ORANGE)], 1)]
        for y, (name, parts, count) in enumerate(rows[::-1]):
            x = 0
            limit = count * (R + M + G) if 'Jacobi' in name else None
            for lab, w, c in parts:
                if limit is not None and x >= limit - 1e-9:
                    break
                ww = min(w, limit - x) if limit is not None else w
                ax.barh(y, ww, left=x, color=c, edgecolor='white', height=0.6)
                if ww > 0.06 * (2.8 * (R + M + G)):
                    ax.text(x + ww / 2, y, lab, ha='center', va='center', color='white', fontsize=7)
                x += ww
            ax.text(x + 0.01 * x, y, ' %.0f us' % x, va='center', fontsize=8)
        ax.set_yticks([0, 1])
        ax.set_yticklabels([r[0] for r in rows[::-1]], fontsize=8)
        ax.set_title('%dD, the finest grid of an adaptive tree (%d blocks): the same smoothing'
                     % (d, cost[d]['blocks']), fontsize=8.5, loc='left')
        ax.set_xlabel('time, us (measured parts: relax %.0f, mix %.0f, ghosts %.0f)' % (R, M, G),
                      fontsize=7.5)
        ax.set_xlim(0, 3.3 * (R + M + G))
    fig.tight_layout()
    fig.savefig('fig3_sweeps.png', dpi=DPI, bbox_inches='tight')
    plt.close(fig)


# 4. the parity rule on the offsets ----------------------------------------
def fig_rule():
    fig, ax = plt.subplots(figsize=(4.6, 4.6))
    R = 2
    for j in range(-R, R + 1):
        for i in range(-R, R + 1):
            s = abs(i) + abs(j)
            if s == 0:
                c, t = '#ecf0f1', 's[]\nown cell'
            elif s % 2:
                c, t = GREEN, 's[%d,%d]' % (i, j) if j else 's[%d]' % i
            else:
                c, t = ORANGE, 's[%d,%d]' % (i, j) if j else 's[%d]' % i
            ax.add_patch(Rectangle((i - 0.5, j - 0.5), 1, 1, facecolor=c, edgecolor='white', lw=2))
            ax.text(i, j, t, ha='center', va='center', fontsize=7.5,
                    color='black' if s == 0 else 'white')
    ax.set_xlim(-R - 0.6, R + 0.6)
    ax.set_ylim(-R - 1.6, R + 0.6)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.add_patch(Rectangle((-R - 0.4, -R - 1.45), 0.4, 0.4, facecolor=GREEN))
    ax.text(-R + 0.1, -R - 1.25, '|i|+|j| odd: the other colour, allowed', va='center', fontsize=7.5)
    ax.add_patch(Rectangle((-R - 0.4, -R - 0.95), 0.4, 0.4, facecolor=ORANGE))
    ax.text(-R + 0.1, -R - 0.75, '|i|+|j| even: the same colour, warned', va='center', fontsize=7.5)
    ax.set_title('Reading s at an offset, in a loop that writes s\non one colour (2D)', fontsize=9)
    fig.savefig('fig4_rule.png', dpi=DPI, bbox_inches='tight')
    plt.close(fig)


# 5. smoothing factors -----------------------------------------------------
def fig_smoothing():
    names = [('jacobi_1', 'Jacobi, w = 1'), ('jacobi', 'weighted Jacobi\n(4/5, 6/7)'),
             ('lex_gs', 'lexicographic\nGauss-Seidel'), ('red_black', 'red-black\nGauss-Seidel')]
    fig, ax = plt.subplots(figsize=(6.4, 3.0))
    x = np.arange(len(names))
    for k, (d, c) in enumerate(((2, '#2980b9'), (3, '#8e44ad'))):
        v = [lfa[str(d)][n] for n, _ in names]
        b = ax.bar(x + (k - 0.5) * 0.36, v, 0.34, color=c, label='%dD' % d)
        for xi, vi in zip(x + (k - 0.5) * 0.36, v):
            ax.text(xi, vi + 0.02, '%.3f' % vi, ha='center', va='bottom', fontsize=7.5)
    ax.set_xticks(x)
    ax.set_xticklabels([n for _, n in names], fontsize=8)
    ax.set_ylabel('smoothing factor per sweep')
    ax.set_ylim(0, 1.15)
    ax.legend(frameon=False, fontsize=8)
    ax.set_title('High-frequency error left after one sweep (local Fourier analysis, '
                 '(2d+1)-point Laplacian)', fontsize=8.5)
    fig.savefig('fig5_smoothing.png', dpi=DPI, bbox_inches='tight')
    plt.close(fig)


# 6. a random error smoothed -----------------------------------------------
def fig_error():
    n = 64
    rng = np.random.default_rng(3)
    e0 = rng.uniform(-1, 1, (n, n))

    def nb(e):
        p = np.pad(e, 1)                         # zero at the walls
        return p[2:, 1:-1] + p[:-2, 1:-1] + p[1:-1, 2:] + p[1:-1, :-2]

    def jac(e, w=0.8):
        return (1 - w) * e + w * nb(e) / 4

    I, J = np.meshgrid(np.arange(n), np.arange(n), indexing='ij')
    red = (I + J) % 2 == 0

    def rb(e):
        e = e.copy()
        e[red] = (nb(e) / 4)[red]
        e[~red] = (nb(e) / 4)[~red]
        return e

    def high_rms(e):
        F = np.fft.fft2(e)
        k = np.fft.fftfreq(n) * 2 * np.pi
        K1, K2 = np.meshgrid(k, k, indexing='ij')
        h = (np.abs(K1) >= np.pi / 2) | (np.abs(K2) >= np.pi / 2)
        return math.sqrt(np.sum(np.abs(F[h]) ** 2)) / n

    fields = [('the error at first', e0)]
    ej, er = e0, e0
    for s in range(2):
        ej, er = jac(ej), rb(er)
    fields += [('weighted Jacobi, 2 sweeps', ej), ('red-black GS, 2 sweeps', er)]
    fig, axes = plt.subplots(1, 3, figsize=(7.6, 2.8))
    h0 = high_rms(e0)
    for ax, (t, f) in zip(axes, fields):
        ax.imshow(f.T, origin='lower', cmap='RdBu_r', vmin=-1, vmax=1)
        ax.set_title('%s\nhigh frequencies: %.0f%% left' % (t, 100 * high_rms(f) / h0), fontsize=8)
        ax.axis('off')
    fig.savefig('fig6_error.png', dpi=DPI, bbox_inches='tight')
    plt.close(fig)


# 7. a row, in vector lanes -------------------------------------------------
def fig_simd():
    fig, axes = plt.subplots(2, 1, figsize=(7.6, 2.9))
    n, lanes = 12, 4
    for ax, kind in zip(axes, ('lex', 'rb')):
        for i in range(n):
            if kind == 'lex':
                c = '#7f8c8d'
            else:
                c = RED if i % 2 == 0 else BLACK
            ax.add_patch(Rectangle((i, 0), 1, 1, facecolor=c, edgecolor='white', lw=1))
            ax.text(i + 0.5, -0.25, str(i), ha='center', va='top', fontsize=7)
        for v in range(0, n, lanes):
            ax.add_patch(Rectangle((v + 0.03, -0.05), lanes - 0.06, 1.1, fill=False,
                                   edgecolor='#2980b9', lw=1.2, ls='--'))
        if kind == 'lex':
            for i in range(1, n):
                ax.add_patch(FancyArrowPatch((i - 0.35, 0.75), (i + 0.35, 0.75), arrowstyle='-|>',
                                             mutation_scale=8, color='#f1c40f', lw=1.2))
            ax.text(n + 0.3, 0.5, 'in place, lexicographic:\ncell i needs the new i - 1:\n'
                    'a chain through each vector\n(wrong under SIMD)', va='center', fontsize=7.5)
        else:
            for i in range(0, n, 2):
                for dx in (-1, 1):
                    if 0 <= i + dx < n:
                        ax.add_patch(FancyArrowPatch((i + 0.5 + dx * 0.9, 0.25), (i + 0.5 + dx * 0.15, 0.25),
                                                     arrowstyle='-|>', mutation_scale=7,
                                                     color='#f1c40f', lw=1.0))
            ax.text(n + 0.3, 0.5, 'red pass: writes even cells,\nreads odd ones, unchanged:\n'
                    'every lane independent\n(safe under SIMD)', va='center', fontsize=7.5)
        ax.set_xlim(-0.3, n + 6.5)
        ax.set_ylim(-0.8, 1.4)
        ax.set_aspect('equal')
        ax.axis('off')
    axes[0].set_title('A row of cells in vectors of 4 (dashed)', fontsize=8.5, loc='left')
    fig.savefig('fig7_simd.png', dpi=DPI, bbox_inches='tight')
    plt.close(fig)


for f in (fig_checkerboard, fig_blocks, fig_sweeps, fig_rule, fig_smoothing, fig_error, fig_simd):
    f()
print('ok')
