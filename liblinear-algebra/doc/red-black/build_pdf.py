"""The design note on red-black Gauss-Seidel, as a PDF: python3 build_pdf.py OUT.pdf
Reads lfa.json, cost2.json, cost3.json and fig*.png from this directory."""
import json
import math
import sys
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm
from reportlab.lib import colors
from reportlab.lib.enums import TA_JUSTIFY, TA_CENTER
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (BaseDocTemplate, Frame, PageTemplate, Paragraph, Spacer,
                                Image, Table, TableStyle, Preformatted, KeepTogether,
                                PageBreak, ListFlowable, ListItem)

F = '/usr/share/fonts/truetype/dejavu/'
pdfmetrics.registerFont(TTFont('DV', F + 'DejaVuSans.ttf'))
pdfmetrics.registerFont(TTFont('DV-B', F + 'DejaVuSans-Bold.ttf'))
pdfmetrics.registerFont(TTFont('DV-I', F + 'DejaVuSans-Oblique.ttf'))
pdfmetrics.registerFont(TTFont('DV-BI', F + 'DejaVuSans-BoldOblique.ttf'))
pdfmetrics.registerFont(TTFont('DVM', F + 'DejaVuSansMono.ttf'))
from reportlab.pdfbase.pdfmetrics import registerFontFamily
registerFontFamily('DV', normal='DV', bold='DV-B', italic='DV-I', boldItalic='DV-BI')

lfa = json.load(open('lfa.json'))
cost = {2: json.load(open('cost2.json')), 3: json.load(open('cost3.json'))}

INK, ACCENT, SOFT = colors.HexColor('#1f2d3a'), colors.HexColor('#c0392b'), colors.HexColor('#f4f6f7')
S = {
    'title': ParagraphStyle('t', fontName='DV-B', fontSize=19, leading=24, textColor=INK, spaceAfter=4),
    'sub': ParagraphStyle('s', fontName='DV', fontSize=10.5, leading=14, textColor=colors.HexColor('#566573'),
                          spaceAfter=14),
    'h1': ParagraphStyle('h1', fontName='DV-B', fontSize=13.5, leading=18, textColor=INK, spaceBefore=14,
                         spaceAfter=6, keepWithNext=1),
    'h2': ParagraphStyle('h2', fontName='DV-B', fontSize=10.8, leading=14, textColor=INK, spaceBefore=10,
                         spaceAfter=4, keepWithNext=1),
    'body': ParagraphStyle('b', fontName='DV', fontSize=9.4, leading=13.6, alignment=TA_JUSTIFY,
                           spaceAfter=6, textColor=colors.HexColor('#222222')),
    'cap': ParagraphStyle('c', fontName='DV-I', fontSize=8.2, leading=11, textColor=colors.HexColor('#566573'),
                          alignment=TA_CENTER, spaceBefore=3, spaceAfter=10),
    'code': ParagraphStyle('code', fontName='DVM', fontSize=7.6, leading=10, textColor=colors.HexColor('#1b2631'),
                           backColor=SOFT, borderPadding=(5, 6, 5, 6), spaceBefore=4, spaceAfter=10,
                           leftIndent=4, rightIndent=4),
    'cell': ParagraphStyle('cell', fontName='DV', fontSize=8.2, leading=10.5),
    'cellb': ParagraphStyle('cellb', fontName='DV-B', fontSize=8.2, leading=10.5, textColor=colors.white),
    'box': ParagraphStyle('box', fontName='DV', fontSize=9.2, leading=13.2, alignment=TA_JUSTIFY),
}


def P(t, s='body'):
    return Paragraph(t, S[s])


def code(t):
    return Preformatted(t.strip('\n'), S['code'])


def fig(path, width_mm, caption):
    from reportlab.lib.utils import ImageReader
    iw, ih = ImageReader(path).getSize()
    w = width_mm * mm
    return KeepTogether([Image(path, width=w, height=w * ih / iw), P(caption, 'cap')])


def table(rows, widths, header=True, align_right_from=None):
    data = [[Paragraph(str(c), S['cellb'] if header and r == 0 else S['cell']) for c in row]
            for r, row in enumerate(rows)]
    t = Table(data, colWidths=[w * mm for w in widths], hAlign='LEFT')
    st = [('GRID', (0, 0), (-1, -1), 0.4, colors.HexColor('#d5d8dc')),
          ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
          ('TOPPADDING', (0, 0), (-1, -1), 3), ('BOTTOMPADDING', (0, 0), (-1, -1), 3)]
    if header:
        st += [('BACKGROUND', (0, 0), (-1, 0), INK)]
    for r in range(1, len(rows)):
        if r % 2 == 0:
            st.append(('BACKGROUND', (0, r), (-1, r), SOFT))
    t.setStyle(TableStyle(st))
    return t


def box(paragraphs):
    t = Table([[[Paragraph(p, S['box']) for p in paragraphs]]], colWidths=[170 * mm])
    t.setStyle(TableStyle([('BACKGROUND', (0, 0), (-1, -1), colors.HexColor('#fdf2f0')),
                           ('LINEBEFORE', (0, 0), (0, -1), 3, ACCENT),
                           ('LEFTPADDING', (0, 0), (-1, -1), 10), ('RIGHTPADDING', (0, 0), (-1, -1), 10),
                           ('TOPPADDING', (0, 0), (-1, -1), 7), ('BOTTOMPADDING', (0, 0), (-1, -1), 7)]))
    return t


def bullets(items):
    return ListFlowable([ListItem(P(i), leftIndent=12, value='•') for i in items],
                        bulletType='bullet', start='•', leftIndent=12, bulletFontName='DV')


def on_page(c, doc):
    c.saveState()
    c.setFont('DV', 7.5)
    c.setFillColor(colors.HexColor('#7f8c8d'))
    c.drawString(20 * mm, 12 * mm, 'Red-black Gauss-Seidel and the neighbour warning: a design note')
    c.drawRightString(190 * mm, 12 * mm, str(doc.page))
    c.restoreState()


# ------------------------------------------------------------------ numbers
mu = {d: lfa[str(d)] for d in (2, 3)}
nj = {d: math.log(mu[d]['red_black']) / math.log(mu[d]['jacobi']) for d in (2, 3)}
jac_cost = {d: nj[d] * (cost[d]['R'] + cost[d]['M'] + cost[d]['G']) for d in (2, 3)}
rb_cost = {d: cost[d]['R'] + 2 * cost[d]['G'] for d in (2, 3)}
breakeven = {d: (jac_cost[d] - 2 * cost[d]['G']) / cost[d]['R'] for d in (2, 3)}

story = []
story += [P('Red-black Gauss-Seidel for liblinear-algebra’s multigrid, and cwrap’s neighbour warning', 'title'),
          P('A design note. Proposal only: nothing here is implemented. '
            'liblinear-algebra (multigrid.c), libb-amr (block-loops.h), libc-pp (cwrap). 2026-11-12.', 'sub')]

story.append(box([
    '<b>In short.</b> The multigrid solver smooths with weighted Jacobi. Red-black Gauss-Seidel '
    'smooths much better per sweep: it leaves %.2f of the high-frequency error in 2D, against %.2f '
    '(%.3f against %.3f in 3D).' % (mu[2]['red_black'], mu[2]['jacobi'], mu[3]['red_black'], mu[3]['jacobi']),
    'On libb-amr’s blocks it needs one more ghost update per sweep. Measured on our code, it still costs '
    'about <b>%.1f× less (2D) and %.1f× less (3D)</b> for the same smoothing.'
    % (jac_cost[2] / rb_cost[2], jac_cost[3] / rb_cost[3]),
    'cwrap now warns on it, because the loop writes <font face="DVM">da[]</font> and reads '
    '<font face="DVM">da[1]</font>. The proposal: <b>no blanket opt-out</b>. Instead, a loop over one colour, '
    'in which cwrap checks that every neighbour read has the other colour. That is a rule it can prove at '
    'compile time, so the guarantee the warning exists for is kept.',
]))

# ------------------------------------------------------------------ 1
story += [P('1. Where the smoother sits', 'h1'),
          P('A multigrid cycle (liblinear-algebra/multigrid.c) solves <i>L a = b</i> for '
            '<i>L a</i> = ∇·(α∇a) + λa. It restricts the residual to every level of the tree. Then, '
            'from the coarsest level up, it prolongs a correction <font face="DVM">da</font> into each '
            'level and <b>smooths</b> it there, a few sweeps of a relaxation on the grid at that depth '
            '(<font face="DVM">foreach_level_or_leaf (l)</font>).'),
          P('A smoother does not need to solve anything. It only needs to remove the error components that '
            'oscillate from cell to cell: the <b>high frequencies</b>, which the coarser level cannot '
            'represent. The smooth part is left to the coarse levels. So a smoother is judged by its '
            '<b>smoothing factor</b> μ, the share of high-frequency error one sweep leaves at worst. '
            'μ<super>ν</super> for ν sweeps bounds what the cycle can do.'),
          P('Today: weighted Jacobi', 'h2'),
          P('Each cell gets the value that would zero its residual, given its neighbours’ old values. That '
            'value goes to a field of its own, <font face="DVM">dn</font>, and is then mixed with the old '
            'one by a weight ω:'),
          code('''
/* the faces of the cell: n = sum of alpha_f da_n, d = sum of alpha_f */
dn[] = (n / h2 - res[]) / (d / h2 - lambda[]);      /* pass 1: relax  */
da[] = (1 - omega) * da[] + omega * dn[];           /* pass 2: mix    */
'''),
          P('The weight matters. Plain Jacobi (ω = 1) has a smoothing factor of <b>1.0</b>: the '
            'checkerboard mode, +1 −1 +1 …, flips sign every sweep and is never damped. '
            'ω = 2d/(2d+1) (4/5 in 2D, 6/7 in 3D) is the optimum, giving %.2f and %.3f. Jacobi’s great '
            'merit is that every cell reads only old values, so the order of the cells does not matter, '
            'and neither do threads, vectors or a GPU. That is what cwrap’s rules protect.'
            % (mu[2]['jacobi'], mu[3]['jacobi'])),
          ]

# ------------------------------------------------------------------ 2
story += [P('2. Red-black Gauss-Seidel', 'h1'),
          P('Colour each cell by the parity of its index: <b>red</b> if <i>i + j + k</i> is even, '
            '<b>black</b> if odd. In the (2d+1)-point stencil, every neighbour of a cell has the other '
            'colour. A sweep is two passes, each updating its cells in place, with no second field:'),
          bullets(['<b>red pass</b>: every red cell takes the value that zeros its residual. It reads only '
                   'black neighbours, which this pass does not change;',
                   '<b>black pass</b>: the same for the black cells, reading the red ones just updated.']),
          Spacer(1, 4),
          fig('fig1_checkerboard.png', 98,
              'Figure 1. A block of 8 × 8 cells and its two ghost layers, coloured. The point (yellow) is '
              'red. Its four neighbours (arrows) are black, so they are unchanged during the red pass. '
              's[1,1] and s[2] are red: the same colour, possibly updated already.'),
          P('Within a pass, no cell reads a cell the pass writes. So a pass is as order-free as Jacobi: any '
            'order of blocks, rows, cells or vector lanes gives the same result. Across the two passes it is '
            'Gauss-Seidel: the black cells see the new red values. That is where the better smoothing comes '
            'from.'),
          P('How much better', 'h2'),
          P('The smoothing factors below come from a local Fourier analysis, computed for this note rather '
            'than quoted. On a periodic grid, μ is the spectral radius of <i>Q S</i>, where <i>S</i> is one '
            'sweep’s error operator and <i>Q</i> projects onto the high frequencies (some |θ<sub>a</sub>| ≥ '
            'π/2). The method reproduces the known 0.6 of weighted Jacobi in 2D.'),
          fig('fig5_smoothing.png', 140,
              'Figure 2. Smoothing factors per sweep, the (2d+1)-point Laplacian. Lower is better. Plain '
              'Jacobi does not smooth at all.'),
          table([['smoother', '2D μ', '3D μ', 'Jacobi sweeps for one of it (2D / 3D)'],
                 ['weighted Jacobi (now)', '%.3f' % mu[2]['jacobi'], '%.3f' % mu[3]['jacobi'], '1 / 1'],
                 ['lexicographic Gauss-Seidel', '%.3f' % mu[2]['lex_gs'], '%.3f' % mu[3]['lex_gs'],
                  '%.2f / %.2f' % (math.log(mu[2]['lex_gs']) / math.log(mu[2]['jacobi']),
                                   math.log(mu[3]['lex_gs']) / math.log(mu[3]['jacobi']))],
                 ['red-black Gauss-Seidel', '%.3f' % mu[2]['red_black'], '%.3f' % mu[3]['red_black'],
                  '%.2f / %.2f' % (nj[2], nj[3])]],
                [52, 22, 22, 74]),
          Spacer(1, 8),
          P('The last column is ln μ<sub>RB</sub> / ln μ<sub>Jacobi</sub>: how many Jacobi sweeps leave '
            'as little high-frequency error as one red-black sweep. Figure 3 shows the same on a field: '
            'random error on a 64 × 64 grid with zero walls, after two sweeps of each.'),
          fig('fig6_error.png', 160,
              'Figure 3. A random error, then two sweeps of weighted Jacobi (ω = 4/5) and two of red-black '
              'Gauss-Seidel. "Left" is the RMS of the high-frequency part against the start. Red-black leaves '
              'a third of what Jacobi does, and a smoother field for the coarse levels.'),
          P('The coarsest level gains too', 'h2'),
          P('The coarsest grid is a whole block, at least 16 × 16 cells, which the tree cannot coarsen. The '
            'solver gives it m² sweeps (256 for m = 16), where it must really solve, not just smooth. For '
            'the consistently ordered matrices of this stencil, red-black Gauss-Seidel’s spectral radius is '
            'Jacobi’s squared. So about half the sweeps reach the same accuracy there.'),
          ]

# ------------------------------------------------------------------ 3
story += [P('3. On libb-amr’s blocks', 'h1'),
          P('3.1 The colour of a cell, from its block’s own indices', 'h2'),
          P('A cell’s global index along an axis is <i>i + m I</i>: its index <i>i</i> in the block, '
            'plus <i>m</i> times the block’s index <i>I</i> at its level. The block size <i>m</i> is a '
            'power of two, at least 16, so <i>m I</i> is even. The colour can therefore come from the local '
            'indices alone, with nothing per block. The padding of the aligned layout (<font face="DVM">ox'
            '</font>) does not enter: the parity uses interior coordinates, <font face="DVM">ci − ox</font>, '
            '<font face="DVM">cj − g</font>.'),
          fig('fig2_blocks.png', 165,
              'Figure 4. Two blocks of level l side by side: the colouring runs on across their edge. At a '
              'level jump (right), the fine block’s ghosts are prolonged from the coarser leaf, which has '
              'its own colouring at its own level.'),
          P('At a level jump the colouring does not carry across. The ghosts there are prolonged values, '
            'neither red nor black, and the coarser leaf, a member of the same depth-l grid, is coloured at '
            'its own level. That is harmless: within a pass every ghost is a fixed value, read as such. '
            'What makes the passes order-free is that a pass never writes what it reads, and that holds '
            'everywhere.'),
          P('3.2 A ghost update between the two passes', 'h2'),
          P('The black pass reads red cells across block edges, through the ghosts, so the ghosts of the red '
            'values must be made between the passes. A red-black sweep is therefore <b>red, ghosts, black, '
            'ghosts</b>: two ghost updates, against Jacobi’s one. On adaptive levels '
            '<font face="DVM">ghosts_update_level</font> also restricts the coarser parents. So it is worth '
            'measuring whether the better smoothing pays for it.'),
          P('It does, by a wide margin. The parts of a sweep were timed on the finest grid of an adaptive tree '
            '(2D: levels %d to %d, %d blocks; 3D: levels %d to %d, %d blocks; -O2, double): one relaxation '
            'pass <i>R</i>, Jacobi’s mixing pass <i>M</i>, one ghost update <i>G</i>. Jacobi then costs '
            '<i>n</i>(R + M + G) for the smoothing of one red-black sweep, which costs R + 2G (the two half '
            'passes are taken as R together).'
            % (cost[2]['levels'][0], cost[2]['levels'][1], cost[2]['blocks'],
               cost[3]['levels'][0], cost[3]['levels'][1], cost[3]['blocks'])),
          table([['', 'R (relax)', 'M (mix)', 'G (ghosts)', 'Jacobi, n sweeps', 'red-black, 1 sweep', 'ratio'],
                 ['2D', '%.0f µs' % cost[2]['R'], '%.0f µs' % cost[2]['M'], '%.0f µs' % cost[2]['G'],
                  '%.2f × → %.0f µs' % (nj[2], jac_cost[2]), '%.0f µs' % rb_cost[2],
                  '%.1f×' % (jac_cost[2] / rb_cost[2])],
                 ['3D', '%.0f µs' % cost[3]['R'], '%.0f µs' % cost[3]['M'], '%.0f µs' % cost[3]['G'],
                  '%.2f × → %.0f µs' % (nj[3], jac_cost[3]), '%.0f µs' % rb_cost[3],
                  '%.1f×' % (jac_cost[3] / rb_cost[3])]],
                [12, 22, 20, 22, 36, 30, 16]),
          Spacer(1, 6),
          fig('fig3_sweeps.png', 165,
              'Figure 5. The same smoothing, to scale: about 2.7 sweeps of weighted Jacobi (2D; 2.4 in 3D) '
              'against one red-black sweep, from the measured parts.'),
          P('The relaxation pass dominates: G is about a third of R. So the second ghost update is cheap '
            'next to the passes it saves. The one unknown is the cost of the half passes, which visit every '
            'other cell of a row (§3.3). Red-black still breaks even if those strided passes together cost '
            '<b>%.1f R (2D) or %.1f R (3D)</b>, far above the 1–1.5 R expected.'
            % (breakeven[2], breakeven[3])),
          P('3.3 Rows, vectors and threads', 'h2'),
          P('A pass visits every other cell of a row: it starts at the first cell of its colour and steps by '
            'two. That is a strided access, which vectorises less well than a contiguous row, and the reason '
            'the half passes should be measured. It is also exactly what makes them safe to vectorise. A red '
            'pass writes even cells and reads odd ones, so no vector lane reads what another lane writes. '
            'Lexicographic Gauss-Seidel, by contrast, chains every cell to the one before it.'),
          fig('fig7_simd.png', 160,
              'Figure 6. A row in vectors of four. Above: in-place lexicographic Gauss-Seidel needs cell i−1’s '
              'new value, a dependency through each vector, so it is wrong once vectorised. Below: the red pass '
              'reads only black cells, so every lane is independent.'),
          P('Across blocks, each pass is as parallel as Jacobi. Its blocks read only ghosts, and cells of the '
            'other colour, that nothing changes during the pass. The result is therefore deterministic for '
            'any number of threads, in any order, and on a GPU.'),
          ]

# ------------------------------------------------------------------ 4
story += [P('4. cwrap’s warning, and why an opt-out is not needed', 'h1'),
          P('4.1 What the warning guards', 'h2'),
          P('cwrap warns when a stencil loop reads a neighbour of a scalar it also writes:'),
          code('''
foreach () { p[] = x; q[] = p[1,0] - p[-1,0]; }
warning: 'p[1,0]' reads a neighbour of 'p', which this foreach () writes: its ghosts
         are stale, and its neighbours old or new by the order of the cells; write it
         in a foreach () of its own first
'''),
          P('Both halves of the message are real hazards. The neighbour inside the block is old or new '
            'depending on the order the cells are visited in, which changes with the vector width, the '
            'threads or a GPU. The neighbour across a block edge is a ghost, older still. The warning keeps '
            'a stencil loop’s result independent of how it runs.'),
          P('Red-black Gauss-Seidel written today, with the loops we have, trips it, and rightly so. cwrap '
            'cannot see that the <font face="DVM">if</font> keeps the colours apart:'),
          code('''
foreach_level_or_leaf (l)
  if (((ci - ox) + (cj - g)) % 2 == colour) {           /* one colour */
    double n = 0, d = 0;
    MG_NEIGHBOURS (da, alpha, n, d);                    /* reads da[±1]: warned */
    da[] = (n / h2 - res[]) / (d / h2 - lambda[]);      /* writes da[]        */
  }
'''),
          P('It also visits every cell to update half of them.'),
          P('4.2 Option A: a blanket opt-out', 'h2'),
          P('An annotation on the loop, a <font face="DVM">_Pragma</font>, a keyword or an attribute, telling '
            'cwrap not to warn. It is simple to build, but:'),
          bullets(['it switches the check off for <b>every</b> read in the loop. A mistaken '
                   '<font face="DVM">da[1,1]</font>, the same colour and possibly updated, would pass '
                   'silently;',
                   'it adds declaration syntax to cwrap, which we have kept out so far (cwrap’s rules come '
                   'from the names loops give, not from new syntax);',
                   'it still visits every cell with a branch.']),
          P('4.3 Option B (proposed): loops over one colour, and a parity rule', 'h2'),
          P('A cell loop that visits one colour only, and says so the way libb-amr’s loops already say what '
            'they are: by a <b>name it gives</b>, <font face="DVM">checkerboard</font>. The other such names '
            'are <font face="DVM">ghosts_by_caller</font> for the level loops, <font face="DVM">cc</font> '
            'and <font face="DVM">fc</font> for the loops over two levels. In a loop whose cell loop gives '
            '<font face="DVM">checkerboard</font>, cwrap knows the point has one colour. A read of a scalar '
            'the loop writes, at an offset (i, j, k), is then judged by the parity of |i| + |j| + |k|:'),
          fig('fig4_rule.png', 82,
              'Figure 7. The parity rule. A read at an odd offset sum has the other colour, never written '
              'in this pass, and is allowed. An even sum is the same colour, so it stays warned.'),
          table([['the read', 'parity', 'cwrap'],
                 ['s[] (own cell)', '0', 'allowed, as now'],
                 ['s[1], s[0,-1], s[2,1] …', 'odd: the other colour', 'allowed: unchanged during the pass'],
                 ['s[1,1], s[2], s[0,2] …', 'even: the same colour', 'warned, as now (or an error: Q2)'],
                 ['s[i,j] with i, j expressions', 'known to gcc only',
                  'a _Static_assert that (i + j + k) is odd, as cwrap already asserts the reach'],
                 ['foreach_neighbor (r) reading s', 'both parities', 'warned: it reaches the same colour']],
                [50, 40, 80]),
          Spacer(1, 8),
          P('The rule is <b>proved, not trusted</b>. Constant offsets are decided while translating. Offsets '
            'written as expressions become a compile-time assertion, so a same-colour read is a compile '
            'error with no cost at run time. That is the same mechanism that already keeps offsets within '
            'the ghost layers:'),
          code('''
/* what cwrap would add for da[a, b], a and b constant expressions */
_Static_assert (((a) + (b)) % 2 != 0,
                "da[a, b] reads the same colour, in a loop over one colour");
'''),
          P('<b>On the leaves, the ghosts come for free.</b> With a leaf loop over one colour, '
            '<font face="DVM">foreach_colour (c)</font>, cwrap’s lazy ghosts already do the right thing. '
            'The red pass marks <font face="DVM">da</font> stale, and the black pass, reading its neighbours, '
            'has it updated first: the ghost update between the passes is written by cwrap, not by hand. '
            'On levels (<font face="DVM">ghosts_by_caller</font>) the caller makes it, as the solver does '
            'now.'),
          P('<b>What it does not do.</b> Like every cwrap rule, it reasons about scalar names. Two scalar '
            'variables holding the same id at run time are not seen, the same as today.'),
          ]

# ------------------------------------------------------------------ 5
story += [P('5. What it would take', 'h1'),
          P('libb-amr: a cell loop over one colour', 'h2'),
          code('''
/* every cell of one colour, halo included: (x + y + z) % 2 == colour in
   interior coordinates; a row from its first cell of that colour, by 2 */
@loop foreach_cell_colour (BlockShape shape, int halo, int colour)
  : CIJK, c, shape, checkerboard
{
  int checkerboard = 1;
  (void) checkerboard;
  int h = block_halo (shape, halo);
  for (int cj = shape.g - h; cj < shape.g + shape.m + h; cj++) {   /* (ck too, in 3D) */
    int y = cj - shape.g, x0 = -h;
    int first = x0 + (((x0 + y - colour) % 2) + 2) % 2;
    for (int x = first; x < shape.m + h; x += 2) {
      int ci = shape.ox + x;
      int c = ci + cj * shape.sy;
      {...}
    }
  }
}

@foreach foreach_level_or_leaf_colour (int l, int colour)
  = foreach_level_or_leaf_block (l) foreach_cell_colour (tree_block_shape (), 0, colour);
@foreach foreach_colour (int colour)              /* leaves; question 4 */
  = foreach_leaf () foreach_cell_colour (tree_block_shape (), 0, colour);
'''),
          P('cwrap', 'h2'),
          bullets(['<font face="DVM">checkerboard = loop_gives (cell, "checkerboard")</font>, beside '
                   '<font face="DVM">ghosts_by_caller</font>;',
                   'in the read-of-a-written-scalar check: skip the warning for an odd constant offset sum; '
                   'emit the parity <font face="DVM">_Static_assert</font> for an expression offset; keep '
                   'the warning for an even sum and for <font face="DVM">foreach_neighbor</font>;',
                   'tests: a mock checkerboard loop (as <font face="DVM">stencil_by_caller.c</font>), the '
                   'allowed and warned reads, an expression offset failing to compile.']),
          P('liblinear-algebra: the smoother', 'h2'),
          code('''
static void relax_red_black (scalar da, scalar res, scalar alpha, scalar lambda, int l)
{
  for (int colour = 0; colour < 2; colour++) {
    foreach_level_or_leaf_colour (l, colour) {
      double n = 0, d = 0, h2 = Delta * Delta;
      MG_NEIGHBOURS (da, alpha, n, d);          /* odd offsets: the other colour */
      da[] = (Real) ((n / h2 - res[]) / (d / h2 - lambda[]));
    }
    ghosts_update_level (da, l);                /* red, ghosts, black, ghosts */
  }
}

typedef struct {
  ...
  int smoother;        /* MG_JACOBI or MG_RED_BLACK: question 3 */
} MGParams;
'''),
          P('No <font face="DVM">dn</font> field and no mixing pass: one scalar less per solve. The '
            'Helmholtz form is unchanged; only the update is in place.'),
          P('6. How to decide by measurement', 'h1'),
          bullets(['the half passes’ cost against R: the one unknown of §3.2;',
                   'cycles to a tolerance, and time, for Jacobi (4 sweeps) against red-black (1, 2) in 1D, 2D '
                   'and 3D, on tests/test_poisson.c’s adaptive trees. Then pick nrelax and the default '
                   'smoother;',
                   'the coarsest level: the sweeps red-black needs for the same reduction (expected about '
                   'half);',
                   'the residual per cycle: it should fall from Jacobi’s 0.13–0.28 toward μ<super>ν</super>.']),
          ]

# ------------------------------------------------------------------ 7
story += [P('7. Questions', 'h1'),
          table([['', 'question', 'proposed'],
                 ['1', 'Loops over one colour with the parity rule, or a blanket opt-out?',
                  'loops over one colour: a rule cwrap proves, no new syntax'],
                 ['2', 'In such a loop, a same-colour read (even offset sum): warning or error?',
                  'an error: in a loop over one colour it is always a mistake'],
                 ['3', 'The smoother selectable in MGParams, Jacobi kept, the default by measurement?',
                  'yes; red-black the default if §6 confirms'],
                 ['4', 'Colour loops on levels only, or also foreach_colour () on the leaves?',
                  'both: the leaves get their ghosts from cwrap for free, useful to implicit schemes']],
                [8, 82, 80]),
          Spacer(1, 12),
          P('Appendix: how the numbers were made', 'h2'),
          P('<b>Smoothing factors</b> (lfa.py): the (2d+1)-point Laplacian on a periodic grid of 16<super>2'
            '</super> (2D) or 8<super>3</super> (3D) cells. S is a sweep’s error operator, written as a '
            'matrix; Q is the projector onto the Fourier modes with some |θ<sub>a</sub>| ≥ π/2; μ = ρ(Q S). '
            'Lexicographic Gauss-Seidel is computed from its symbol. The 0.600 of weighted Jacobi (2D) and '
            'the 0.250 of red-black are the textbook values.'),
          P('<b>Costs</b> (cost.c): libb-amr built with cwrap at -O2, double, m = 16. The tree is refined '
            'around (0.3, 0.3); the finest grid, foreach_level_or_leaf at the deepest level. R is one '
            'relaxation pass of multigrid.c’s operator (α a field, as in the solver), M the mixing pass, G '
            'one ghosts_update_level. Each is the best of 7 runs of 20.'),
          P('<b>Figures</b> (figures.py), from those two. The sources are in liblinear-algebra/doc/red-black/; '
            '<font face="DVM">make</font> there rebuilds this file.'),
          ]

doc = BaseDocTemplate(sys.argv[1], pagesize=A4, leftMargin=20 * mm, rightMargin=20 * mm,
                      topMargin=18 * mm, bottomMargin=20 * mm,
                      title='Red-black Gauss-Seidel and the neighbour warning',
                      author='kerala-flood-sim', subject='liblinear-algebra design note')
frame = Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height, id='f')
doc.addPageTemplates([PageTemplate(id='p', frames=[frame], onPage=on_page)])
doc.build(story)
print('written', sys.argv[1])
