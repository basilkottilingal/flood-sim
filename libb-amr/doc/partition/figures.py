#!/usr/bin/env python3
"""
.. figures.py: the generated TikZ of partition-methods.tex. Each fig-*.tex
.. is computed, not drawn by hand: a quadtree refined around a circle, its
.. Morton order and cuts, and the migration of a 1D example by the two
.. methods (one pass, or moves of one rank repeated).
"""
import math

L, NPE = 4, 4

# ------------------------------------------------------------ the quadtree
def morton (i, j, l):
  s = L - l
  i, j, k = i << s, j << s, 0
  for b in range (L):
    k |= ((i >> b) & 1) << (2*b) | ((j >> b) & 1) << (2*b + 1)
  return k

def children (n):
  i, j, l = n
  return [(2*i + (c & 1), 2*j + (c >> 1), l + 1) for c in range (4)]

def parent (n):
  i, j, l = n
  return (i // 2, j // 2, l - 1)

def rect (n):
  i, j, l = n
  h = 1.0 / (1 << l)
  return (i*h, j*h, (i + 1)*h, (j + 1)*h)

def build ():
  leaves = {(i, j, 2) for i in range (4) for j in range (4)}
  def near (n):
    x0, y0, x1, y1 = rect (n)
    h = x1 - x0
    x, y = (x0 + x1)/2, (y0 + y1)/2
    d1 = abs (math.hypot (x - 0.45, y - 0.42) - 0.30)
    return d1 < 0.6*h
  changed = True
  while changed:
    changed = False
    for n in list (leaves):
      if n[2] < L and near (n):
        leaves.remove (n); leaves.update (children (n)); changed = True
  def leaf_at (x, y):
    for l in range (L, -1, -1):
      n = (int (x*(1 << l)), int (y*(1 << l)), l)
      if n in leaves:
        return n
    return None
  changed = True                    # 2:1 balance across faces
  while changed:
    changed = False
    for n in list (leaves):
      if n not in leaves:
        continue
      x0, y0, x1, y1 = rect (n)
      e = 1e-9
      for (x, y) in [(x0 - e, (y0+y1)/2), (x1 + e, (y0+y1)/2),
                     ((x0+x1)/2, y0 - e), ((x0+x1)/2, y1 + e)]:
        if 0 < x < 1 and 0 < y < 1:
          m = leaf_at (x, y)
          if m and m[2] < n[2] - 1:
            leaves.remove (m); leaves.update (children (m)); changed = True
  return sorted (leaves, key = lambda n: morton (*n))

leaves = build ()
_lset = set (leaves)
def leaf_at (x, y):
  for l in range (L, -1, -1):
    n = (int (x*(1 << l)), int (y*(1 << l)), l)
    if n in _lset:
      return n
  return None
N = len (leaves)
middle = {n: (2*NPE*q + NPE) // (2*N) for q, n in enumerate (leaves)}
family = dict (middle)
lset = set (leaves)
split = set ()
for n in leaves:
  if n[2] == 0:
    continue
  p = parent (n)
  ch = children (p)
  if all (c in lset for c in ch):
    family[n] = middle[ch[0]]
    if len ({middle[c] for c in ch}) > 1:
      split.add (p)

S = 6.0   # cm, the square's side

def r (v):
  return "%.3f" % (v*S)

def quadtree (fn, owner, mark = ()):
  out = []
  for n in leaves:
    x0, y0, x1, y1 = rect (n)
    out.append ("\\filldraw[fill=r%d!35, draw=black!35, line width=0.2pt] "
                "(%s,%s) rectangle (%s,%s);" % (owner[n], r (x0), r (y0),
                                                r (x1), r (y1)))
  pts = []
  for n in leaves:
    x0, y0, x1, y1 = rect (n)
    pts.append ("(%s,%s)" % (r ((x0+x1)/2), r ((y0+y1)/2)))
  out.append ("\\draw[black!55, line width=0.4pt, line join=round] "
              + " -- ".join (pts) + ";")
  for n in leaves:
    x0, y0, x1, y1 = rect (n)
    e = 1e-9
    for (x, y, a, b, c, d) in [(x0 - e, (y0+y1)/2, x0, y0, x0, y1),
                               (x1 + e, (y0+y1)/2, x1, y0, x1, y1),
                               ((x0+x1)/2, y0 - e, x0, y0, x1, y0),
                               ((x0+x1)/2, y1 + e, x0, y1, x1, y1)]:
      m = leaf_at (x, y)
      if m and owner[m] != owner[n]:
        out.append ("\\draw[black, line width=1.1pt, line cap=round] "
                    "(%s,%s) -- (%s,%s);" % (r (a), r (b), r (c), r (d)))
  for p in sorted (mark):
    x0, y0, x1, y1 = rect (p)
    out.append ("\\draw[red!80!black, line width=1.6pt, dash pattern=on 3pt off 1.5pt] "
                "(%s,%s) rectangle (%s,%s);" % (r (x0), r (y0), r (x1), r (y1)))
  out.append ("\\draw[black, line width=0.6pt] (0,0) rectangle (%s,%s);"
              % (r (1), r (1)))
  open (fn, "w").write ("\n".join (out) + "\n")

quadtree ("fig-middle.tex", middle, split)
quadtree ("fig-family.tex", family, split)
cnt = [sum (1 for n in leaves if family[n] == k) for k in range (NPE)]
open ("fig-counts.tex", "w").write (
  "\\newcommand{\\nleaves}{%d}\\newcommand{\\nsplit}{%d}"
  "\\newcommand{\\famcounts}{%s}\n" % (N, len (split), ", ".join (map (str, cnt))))

# ------------------------------------------------------------ the skeleton
def touch (a, b):
  ax0, ay0, ax1, ay1 = rect (a)
  bx0, by0, bx1, by1 = rect (b)
  e = 1e-12
  return ax0 <= bx1 + e and bx0 <= ax1 + e and ay0 <= by1 + e and by0 <= ay1 + e

ME = 1
own = [n for n in leaves if family[n] == ME]
halo = [n for n in leaves if family[n] != ME and any (touch (n, m) for m in own)]
skel = set (own) | set (halo)
for n in list (skel):
  while n[2] > 0:
    n = parent (n); skel.add (n)
pruned = set ()
for n in skel:
  if n in lset:
    continue
  for c in children (n):
    if c not in skel:
      pruned.add (c)
out = []
for n in sorted (pruned, key = lambda n: n[2]):
  x0, y0, x1, y1 = rect (n)
  out.append ("\\filldraw[fill=black!12, draw=black!55, line width=0.4pt] "
              "(%s,%s) rectangle (%s,%s);" % (r (x0), r (y0), r (x1), r (y1)))
for n in halo:
  x0, y0, x1, y1 = rect (n)
  out.append ("\\filldraw[fill=r%d!38, draw=r%d!80!black, line width=0.3pt, "
              "pattern color=r%d] (%s,%s) rectangle (%s,%s);"
              % (family[n], family[n], family[n], r (x0), r (y0), r (x1), r (y1)))
for n in own:
  x0, y0, x1, y1 = rect (n)
  out.append ("\\filldraw[fill=r%d!70, draw=black!50, line width=0.3pt] "
              "(%s,%s) rectangle (%s,%s);" % (ME, r (x0), r (y0), r (x1), r (y1)))
out.append ("\\draw[black, line width=0.6pt] (0,0) rectangle (%s,%s);"
            % (r (1), r (1)))
open ("fig-skeleton.tex", "w").write ("\n".join (out) + "\n")
open ("fig-skelcounts.tex", "w").write (
  "\\newcommand{\\nown}{%d}\\newcommand{\\nhalo}{%d}\\newcommand{\\npruned}{%d}"
  "\\newcommand{\\nskel}{%d}\n" % (len (own), len (halo), len (pruned),
                                   len (skel) + len (pruned)))

# ------------------------------------------------------------ edge cut
# the leaf graph G_L (leaves touching by a face or a corner), its edge cut
# under the family partition, the halos H_r and the communication graph
# G_C with V_rs = |H_r^s|, the leaves r receives from s
EL = [(a, b) for x, a in enumerate (leaves) for b in leaves[x + 1:]
      if touch (a, b)]
cut = [(a, b) for (a, b) in EL if family[a] != family[b]]
V = [[0]*NPE for _ in range (NPE)]
for r in range (NPE):
  mine = [n for n in leaves if family[n] == r]
  for n in leaves:
    if family[n] != r and any (touch (n, m) for m in mine):
      V[r][family[n]] += 1
pos = [(0.0, 2.6), (2.6, 2.6), (0.0, 0.0), (2.6, 0.0)]
o = []
for r in range (NPE):
  for s in range (r + 1, NPE):
    if V[r][s] or V[s][r]:
      (x0, y0), (x1, y1) = pos[r], pos[s]
      o.append ("\\draw[black!55, line width=%.2fpt] (%.2f,%.2f) -- (%.2f,%.2f);"
                % (0.4 + 0.08*(V[r][s] + V[s][r]), x0, y0, x1, y1))
      # V_rs written near r, V_sr near s
      for (a, b, ab) in [(r, s, V[r][s]), (s, r, V[s][r])]:
        (xa, ya), (xb, yb) = pos[a], pos[b]
        f = 0.30
        o.append ("\\node[font=\\scriptsize, fill=white, inner sep=1pt] "
                  "at (%.2f,%.2f) {%d};"
                  % (xa + f*(xb - xa), ya + f*(yb - ya), ab))
for r in range (NPE):
  x, y = pos[r]
  o.append ("\\node[circle, fill=r%d, text=white, font=\\footnotesize\\bfseries, "
            "minimum size=7mm, inner sep=0pt] at (%.2f,%.2f) {%d};" % (r, x, y, r))
open ("fig-commgraph.tex", "w").write ("\n".join (o) + "\n")
deg = [sum (1 for s in range (NPE) if V[r][s]) for r in range (NPE)]
open ("fig-cutcounts.tex", "w").write (
  "\\newcommand{\\nedges}{%d}\\newcommand{\\ncut}{%d}"
  "\\newcommand{\\halocounts}{%s}\\newcommand{\\degrees}{%s}"
  "\\newcommand{\\ncommedges}{%d}\n"
  % (len (EL), len (cut), ", ".join (str (sum (V[r])) for r in range (NPE)),
     ", ".join (map (str, deg)),
     sum (1 for r in range (NPE) for s in range (r + 1, NPE) if V[r][s])))

# ------------------------------------------------------------ 1D migration
# 4 ranks, 4 leaves each; rank 0 refines its 4 leaves twice: 64 + 12
old = [0]*64 + [1]*4 + [2]*4 + [3]*4
M = len (old)
W = 11.6 / M       # cm a leaf

def bar (y, owner, label, ticks = None):
  o = ["\\node[anchor=east, font=\\footnotesize] at (-0.15,%.3f) {%s};"
       % (y + 0.18, label)]
  for q, k in enumerate (owner):
    o.append ("\\filldraw[fill=r%d!55, draw=white, line width=0.25pt] "
              "(%.3f,%.3f) rectangle (%.3f,%.3f);" % (k, q*W, y, (q+1)*W, y + 0.36))
  for k in range (NPE):
    qs = [q for q in range (M) if owner[q] == k]
    if qs:
      o.append ("\\node[font=\\scriptsize, text=black] at (%.3f,%.3f) {%d};"
                % ((qs[0] + qs[-1] + 1)*W/2, y + 0.18, len (qs)))
  return o

def graph (x, y, moves):
  """ranks as dots, an arrow a pair (from, to) with its count"""
  o = []
  for k in range (NPE):
    o.append ("\\node[circle, fill=r%d, inner sep=1.6pt, text=white, "
              "font=\\tiny\\bfseries] (g%d) at (%.3f,%.3f) {%d};"
              % (k, k, x + 0.85*k, y, k))
  for (a, b), c in sorted (moves.items ()):
    bend = 30 + 14*(abs (b - a) - 1)
    up = a == 0            # from rank 0 above, the others below
    o.append ("\\draw[-{Stealth[length=3pt]}, line width=0.5pt, r%d] (g%d) "
              "to[bend %s=%d] node[font=\\tiny, %s, inner sep=0.5pt] {%d} (g%d);"
              % (a, a, "left" if up else "right", bend,
                 "above" if up else "below", c, b))
  o.append ("\\begin{scope}[on background layer]\\node[fit=(g0)(g%d), "
            "inner sep=0pt]{};\\end{scope}" % (NPE - 1))
  return o

def moves_of (a, b):
  m = {}
  for q in range (M):
    if a[q] != b[q]:
      m[(a[q], b[q])] = m.get ((a[q], b[q]), 0) + 1
  return m

# one pass: each leaf by its middle
new = [(2*NPE*q + NPE) // (2*M) for q in range (M)]
o = bar (0.9, old, "after adapting") + bar (0, new, "one pass")
g = ["\\begin{scope}[local bounding box=gg]"] + graph (12.1, 0.25, moves_of (old, new)) \
    + ["\\end{scope}"]
open ("fig-onepass.tex", "w").write ("\n".join (o + g) + "\n")
op = moves_of (old, new)

# moves of one rank, repeated (balance () in Basilisk): target from the
# start of the leaf (cw / (tw / npe)), clamped to pid +- 1
rows = [old]
cur = list (old)
while True:
  t = [min (NPE - 1, (q*NPE) // M) for q in range (M)]
  nxt = [max (cur[q] - 1, min (cur[q] + 1, t[q])) for q in range (M)]
  if nxt == cur:
    break
  rows.append (nxt); cur = nxt
o = []
y = 0.9*(len (rows) - 1)
for p, row in enumerate (rows):
  yy = 0.9*(len (rows) - 1 - p)
  o += bar (yy, row, "after adapting" if p == 0 else "pass %d" % p)
  if p > 0:
    o.append ("\\begin{scope}")
    o += graph (12.1, yy + 0.25, moves_of (rows[p-1], row))
    o.append ("\\end{scope}")
open ("fig-iterate.tex", "w").write ("\n".join (o) + "\n")
tot = sum (sum (moves_of (rows[p-1], rows[p]).values ()) for p in range (1, len (rows)))
open ("fig-migcounts.tex", "w").write (
  "\\newcommand{\\mleaves}{%d}\\newcommand{\\npasses}{%d}"
  "\\newcommand{\\movedonce}{%d}\\newcommand{\\moveditr}{%d}"
  "\\newcommand{\\pairsonce}{%d}\n"
  % (M, len (rows) - 1, sum (op.values ()), tot, len (op)))

# ------------------------------------------------------------ z_weights
# a binary tree (1D), leaves finer near x = 0.3; sums up, offsets down
def bin_leaves ():
  ls = [(0, 0)]
  done = False
  while not done:
    done = True
    nl = []
    for (i, l) in ls:
      x = (i + 0.5) / (1 << l)
      if l < 2 or (l < 5 and abs (x - 0.3) < 1.1 / (1 << l)):
        nl += [(2*i, l + 1), (2*i + 1, l + 1)]; done = False
      else:
        nl.append ((i, l))
    ls = nl
  return sorted (ls, key = lambda n: n[0] / (1 << n[1]))

bl = bin_leaves ()
bset = set (bl)
nodes = {}
def sw (n):
  i, l = n
  if n in bset:
    v = 1
  else:
    v = sw ((2*i, l + 1)) + sw ((2*i + 1, l + 1))
  nodes[n] = v
  return v
sw ((0, 0))
cw = {(0, 0): 0}
for n in sorted (nodes, key = lambda n: n[1]):
  if n not in bset:
    i, l = n
    cw[(2*i, l + 1)] = cw[n]
    cw[(2*i + 1, l + 1)] = cw[n] + nodes[(2*i, l + 1)]
nb = len (bl)
def zfig (fn, val, leafval, up):
  o = []
  X = 6.4
  pos = {}
  def place (n):
    i, l = n
    if n in bset:
      pos[n] = bl.index (n) * X / (len (bl) - 1)
    else:
      place ((2*i, l + 1)); place ((2*i + 1, l + 1))
      pos[n] = (pos[(2*i, l + 1)] + pos[(2*i + 1, l + 1)]) / 2
  place ((0, 0))
  for n in nodes:
    i, l = n
    x = pos[n]
    if n not in bset:
      for c in (0, 1):
        cx = pos[(2*i + c, l + 1)]
        o.append ("\\draw[black!45] (%.3f,%.3f) -- (%.3f,%.3f);"
                  % (x, -0.75*l, cx, -0.75*(l + 1)))
  for n in nodes:
    i, l = n
    x = pos[n]
    leaf = n in bset
    rk = min (NPE - 1, (cw[n]*NPE) // nb) if leaf else None
    style = ("fill=r%d!45" % rk) if leaf else "fill=white"
    o.append ("\\node[circle, draw=black!60, %s, inner sep=0pt, minimum size=4.2mm,"
              " font=\\tiny] at (%.3f,%.3f) {%s};"
              % (style, x, -0.75*l, leafval (n) if leaf else val (n)))
  a = "\\draw[-{Stealth}, line width=1pt, black!50] (%.2f,%.2f) -- (%.2f,%.2f);"
  o.append (a % ((-0.35, -3.75, -0.35, 0.0) if up else (-0.35, 0.0, -0.35, -3.75)))
  open (fn, "w").write ("\n".join (o) + "\n")
zfig ("fig-zup.tex", lambda n: str (nodes[n]), lambda n: "1", True)
zfig ("fig-zdown.tex", lambda n: str (cw[n]), lambda n: str (cw[n]), False)
open ("fig-zcounts.tex", "w").write ("\\newcommand{\\zleaves}{%d}\n" % nb)

print ("leaves %d, split families %d, per rank %s; skeleton own %d halo %d "
       "pruned %d; 1D %d leaves, %d passes, moved %d (one pass %d)"
       % (N, len (split), cnt, len (own), len (halo), len (pruned), M,
          len (rows) - 1, tot, sum (op.values ())))

# ------------------------------------------------------------ weighted cut
wts = [1, 1, 2, 1, 3, 1, 1, 1, 2, 2, 1, 1, 4, 1, 1, 2, 1, 1, 3, 1]
WT = sum (wts)
U = 13.6 / WT
o = []
before = 0
for q, w in enumerate (wts):
  k = min (NPE - 1, int (NPE*(before + w/2) / WT))
  x0, x1 = before*U, (before + w)*U
  o.append ("\\filldraw[fill=r%d!50, draw=white, line width=0.6pt] "
            "(%.3f,0) rectangle (%.3f,0.5);" % (k, x0, x1))
  o.append ("\\node[font=\\scriptsize] at (%.3f,0.25) {%d};" % ((x0 + x1)/2, w))
  o.append ("\\fill[black] (%.3f,0.5) circle (0.9pt);" % ((x0 + x1)/2))
  o.append ("\\node[font=\\tiny, text=black!60] at (%.3f,-0.18) {%d};" % (x0, before))
  before += w
for k in range (NPE + 1):
  x = k*WT/NPE*U
  o.append ("\\draw[black, dashed, line width=0.6pt] (%.3f,-0.4) -- (%.3f,0.95);" % (x, x))
  o.append ("\\node[font=\\scriptsize, above] at (%.3f,0.95) {$%s$};"
            % (x, "0" if k == 0 else ("W" if k == NPE else "\\tfrac{%d}{%d}W" % (k, NPE))))
open ("fig-weights.tex", "w").write ("\n".join (o) + "\n")
open ("fig-wcounts.tex", "w").write ("\\newcommand{\\wtotal}{%d}\n" % WT)

# ------------------------------------------------------------ the new table
# the one-pass example again: each new first key K_s given by the one old
# rank holding the first leaf with r >= s; the counts each old rank gives
# (MPI_Allgather of the last pids, then MPI_Allgatherv of the entries)
def table_entries (old, new):
  lasts = [max ([new[q] for q in range (len (old)) if old[q] == o], default = -1)
           for o in range (NPE)]
  counts, last, giver = [], -1, {}
  for o in range (NPE):
    c = lasts[o] - last if lasts[o] > last else 0
    counts.append (c)
    for s in range (last + 1, last + 1 + c):
      giver[s] = (o, min (q for q in range (len (old)) if new[q] >= s))
    last = max (last, lasts[o])
  return lasts, counts, giver

def entries_fig (old, new, y0):
  lasts, counts, giver = table_entries (old, new)
  o = []
  for q in range (M):
    o.append ("\\filldraw[fill=r%d!55, draw=white, line width=0.25pt] "
              "(%.3f,%.3f) rectangle (%.3f,%.3f);" % (new[q], q*W, y0, (q+1)*W, y0 + 0.36))
  for k in range (NPE):                  # the old pieces, outlined
    qs = [q for q in range (M) if old[q] == k]
    if qs:
      o.append ("\\draw[line width=0.9pt] (%.3f,%.3f) rectangle (%.3f,%.3f);"
                % (qs[0]*W, y0 - 0.04, (qs[-1] + 1)*W, y0 + 0.40))
      o.append ("\\node[font=\\scriptsize, anchor=north] at (%.3f,%.3f) "
                "{%d};" % ((qs[0] + qs[-1] + 1)*W/2, y0 - 0.06, counts[k]))
  o.append ("\\node[font=\\scriptsize, anchor=east, align=right] at "
            "(-0.15,%.3f) {new ranks;\\\\entries each old\\\\piece gives};"
            % (y0 + 0.05))
  for s, (g, q) in sorted (giver.items ()):
    o.append ("\\draw[-{Stealth[length=3pt]}, line width=0.6pt] (%.3f,%.3f) "
              "-- (%.3f,%.3f) node[above, font=\\scriptsize] {$K_%d$};"
              % (q*W, y0 + 0.85, q*W, y0 + 0.42, s))
  return o, counts

o1, c1 = entries_fig (old, new, 1.3)
# a small change: one pass done, then rank 2 refines one leaf (2 children)
old2 = list (new)
i2 = old2.index (2)
old2 = old2[:i2] + [2] + old2[i2:]      # one leaf more on rank 2
M2 = len (old2)
new2 = [(2*NPE*q + NPE) // (2*M2) for q in range (M2)]
Msave, M = M, M2
o2, c2 = entries_fig (old2, new2, 0)
M = Msave
open ("fig-entries.tex", "w").write ("\n".join (o1 + o2) + "\n")
open ("fig-entcounts.tex", "w").write (
  "\\newcommand{\\entbig}{%s}\\newcommand{\\entsmall}{%s}\n"
  % (", ".join (map (str, c1)), ", ".join (map (str, c2))))

# ------------------------------------------------------------ Rabenseifner
# p = 8 ranks, the vector in 8 blocks. Reduce-scatter by recursive halving
# (distance 4, 2, 1: keep the half on my side, send the other), then
# allgather by recursive doubling (distance 1, 2, 4). A cell is block b on
# rank r: empty if not held, its shade the number of ranks summed in it
P8 = 8
held = [[{r} for b in range (P8)] for r in range (P8)]   # contributors
win = [(0, P8) for r in range (P8)]
stages = [("start: all of $n$", [row[:] for row in held], "")]
d = P8 // 2
while d >= 1:
  nxt = [[None]*P8 for r in range (P8)]
  nwin = []
  for r in range (P8):
    a, b = win[r]
    mid = (a + b) // 2
    keep = (a, mid) if not (r & d) else (mid, b)
    p = r ^ d
    for x in range (P8):
      if keep[0] <= x < keep[1]:
        nxt[r][x] = held[r][x] | held[p][x]
    nwin.append (keep)
  held, win = nxt, nwin
  stages.append (("halving, $d = %d$, sends $n/%d$" % (d, P8 // d),
                  [row[:] for row in held], ""))
  d //= 2
d = 1
while d < P8:
  nxt = [row[:] for row in held]
  for r in range (P8):
    p = r ^ d
    for x in range (P8):
      if held[p][x] is not None:
        nxt[r][x] = held[p][x]
  held = nxt
  stages.append (("doubling, $d = %d$, sends $n/%d$" % (d, P8 // d),
                  [row[:] for row in held], ""))
  d *= 2
o = []
CW, RW, SH = 0.15, 1.42, 0.55
for k, (label, h, _) in enumerate (stages):
  y = -SH*k
  o.append ("\\node[anchor=east, font=\\scriptsize] at (-0.1,%.3f) {%s};"
            % (y + CW/2, label))
  for r in range (P8):
    for b in range (P8):
      x = r*RW + b*CW
      if h[r][b] is None:
        o.append ("\\draw[black!25, line width=0.2pt] (%.3f,%.3f) rectangle ++(%.3f,%.3f);"
                  % (x, y, CW, CW))
      else:
        n = len (h[r][b])
        o.append ("\\filldraw[fill=r0!%d, draw=black!40, line width=0.2pt] "
                  "(%.3f,%.3f) rectangle ++(%.3f,%.3f);"
                  % (12 + 88*int (math.log2 (n)) // 3, x, y, CW, CW))
for r in range (P8):
  o.append ("\\node[font=\\scriptsize] at (%.3f,%.3f) {rank %d};"
            % (r*RW + P8*CW/2, CW + 0.22, r))
open ("fig-rabenseifner.tex", "w").write ("\n".join (o) + "\n")

# ------------------------------------------------------------ the cost model
# the table of the repartition, old (one MPI_Allreduce of P entries of 24 B,
# a user op) and new (MPI_Allgather of 4 B, MPI_Allgatherv of 24 B a rank),
# in the alpha-beta-gamma model, for a few P: alpha 1 us, beta 0.1 ns/B
# (10 GB/s), gamma 0.5 ns/B (a user function on 24 B entries)
ALPHA, BETA, GAMMA = 1e-6, 0.1e-9, 0.5e-9
def lg (p):
  return math.ceil (math.log2 (p))
def allreduce_rd (p, n):
  return lg (p)*(ALPHA + n*BETA + n*GAMMA)
def allreduce_rab (p, n):
  return 2*lg (p)*ALPHA + 2*(p - 1)/p*n*BETA + (p - 1)/p*n*GAMMA
def allgather_bruck (p, n):            # n: the bytes gathered in all
  return lg (p)*ALPHA + (p - 1)/p*n*BETA
rows = []
for p in (16, 1024, 16384, 131072):
  n = 24*p
  old_rd, old_rab = allreduce_rd (p, n), allreduce_rab (p, n)
  new = allgather_bruck (p, 4*p) + allgather_bruck (p, 24*p)
  size = ("%d\\,B" % n if n < 1e3 else "%.0f\\,kB" % (n/1e3) if n < 1e6
          else "%.1f\\,MB" % (n/1e6))
  rows.append ("%s & %s & %.1f & %.1f & %.1f \\\\"
               % ("{:,}".format (p).replace (",", "\\,"), size,
                  old_rd*1e6, old_rab*1e6, new*1e6))
open ("fig-costs.tex", "w").write ("\\newcommand{\\costrows}{%s}\n"
                                   % " ".join (rows))
