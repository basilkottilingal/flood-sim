# libwebgl

The grid of libb-amr and the manifold of libhmesh in the browser: one HTML
file a scene, drawn by three.js (WebGL), opened from disk. Pan, zoom, pick
a cell or an element and read its values, choose the field and the colour
map. 2D for now.

Two ways in:

- **From a running program:** build a scene from the tree, the mesh, or any
  blocks, lines, triangles and points, and write it.
- **From a saved file:** a `.bamr` file of libb-amr's `tree_output ()`,
  read by `webgl_bamr ()` or by the command `bamr2html`, after the run and
  whatever the number of processes that wrote it.

```c
#include "webgl-tree.h"
#include "webgl-hmesh.h"

WebglScene * S = webgl_scene_new ("step 120");
webgl_tree (S, "grid", 2, (scalar []) { c, h });   /* leaves, c, h */
webgl_hmesh (S, "interface", M);                   /* points, edges */
webgl_scene_write (S, "step-120.html");
webgl_scene_free (S);
```

```sh
bamr2html step-120.bamr          # step-120.html
```

## three.js: in the page, or from the network

three.js r160 is not kept in git (700 KB). Two builds:

| | how | the pages |
|---|---|---|
| **offline** | `sh vendor/fetch.sh` (from jsDelivr, sha256 checked), then build | carry three.js inside: no network to open them |
| **online** | build with `vendor/` empty (or `WEBGL_ONLINE=1 sh gen-assets.sh`) | load three.js r160 from jsDelivr: about 1 MB smaller, the network needed |

`webgl_offline ()` says which a build does. A page that cannot start (no
WebGL2, or no three.js and no network) says so in red. The licence, MIT, is
`vendor/three.LICENSE`. The files, the same bytes as the npm package
`three@0.160.0`:

| file | sha256 |
|---|---|
| `three.module.min.js` | `3e690ac7d180b0aadf0891bea39eec643e29e2d3e75c99b18689518665f69ba6` |
| `OrbitControls.js` | `5a44a9e86a2a0fb11933eed69bc2cd33c76a496854c1aed6ed776efa87d7b064` |

## In the browser

- **Pan** by dragging, **zoom** with the wheel (at the cursor), **r** resets
  the view.
- **A panel a layer:** shown or hidden; the field drawn; the colour map
  (viridis, coolwarm, gray); the range, typed or `auto`; the blocks'
  outlines, a layer's points, a triangulation's edges.
- **Click to pick:** a cell (its block, indices, box and every field there),
  a segment or a triangle (its points, their coordinates, the fields), a
  point. Coordinates in double.
- **Categorical fields** (levels, ranks) in distinct colours (Okabe-Ito),
  the others through a colour map over their range; NaN grey.

## Layers and fields

| layer | what | fields on it |
|---|---|---|
| blocks | n blocks of m x m cells, each a box (x0, y0, width, height): a block-structured grid | `WEBGL_CELL` (n m², block by block, i fastest), `WEBGL_BLOCK` (one a block) |
| lines | points, segments between them | `WEBGL_POINT`, `WEBGL_ELEMENT` |
| triangles | points, triangles | `WEBGL_POINT`, `WEBGL_ELEMENT` |
| points | points | `WEBGL_POINT` |

A field is single precision, categorical or not. The scene copies what it
is given.

**From the tree** (`webgl_tree`, libb-amr, built with cwrap): the leaves
(`tree_caches ()->leaves`, under MPI the halo too) as blocks, a field a
scalar (named by `scalar_name ()`, else `s0`, `s1`, ...), each leaf's level
and owner's rank. Cells are read through `block_at ()`, no lookup in the
level index.

**From the mesh** (`webgl_hmesh`, libhmesh): its points and edges (m = 1)
or triangles (m = 2), the points renumbered with the holes in their ids
left out; every `HmReal` field of the points and of the top level (a
vector's components as `name.x`, `name.y`), the ids, each simplex's length
or area.

**From a .bamr file** (`webgl_bamr`, plain C, no libb-amr): its leaves as
blocks named after the file, a field a scalar (named as in the file, the
file's Reals in single precision), each leaf's level. The file is checked
first (magic, version 1, byte order, codec, size against its sections,
levels, indices) and a bad one adds nothing. A `.bamr` file has no owner
ranks: it does not depend on their number.

## How it works

- **The file** is the template `viewer/viewer.html` with four marks filled:
  the title, escaped; the import map (three.js and OrbitControls as `data:`
  URLs in an offline build, jsDelivr URLs in an online one); the data; and
  `viewer/viewer.js`, inline as a module.
- **The data** is a JSON manifest (layers, fields, their ranges, the names
  of the arrays) and one `<script type="application/octet-stream">` an
  array, base64, little-endian: f32 for boxes, positions and fields, u32
  for indices. Names have `<` escaped, so none ends its script.
- **Positions** are offsets from the scene's lower corner, in single
  precision, the corner itself in double in the manifest: far from the
  origin a float would lose the cell.
- **Compiled in:** `gen-assets.sh` turns the template, the viewer and (if
  present) three.js into C arrays (`xxd -i`, `webgl-assets.c`, not kept),
  so a program needs no file at run time.
- **The viewer** (three.js r160): an orthographic camera, OrbitControls with
  rotation off (MapControls' settings). A block is a quad with a float
  texture of its cells, coloured by a shader, nearest so cells stay sharp;
  lines, triangles and points carry a colour a vertex. three.js's colour
  management is off, so the shader's colours and the vertices' agree.
  Drawn only when something changes.

## API

| | |
|---|---|
| `webgl_scene_new (title)`, `webgl_scene_free (S)` | a scene |
| `webgl_blocks (S, name, n, m, box)` | a blocks layer, `box[4 n]`; its number |
| `webgl_lines (S, name, np, xy, ne, e)` | points `xy[2 np]`, segments `e[2 ne]` |
| `webgl_triangles (S, name, np, xy, ne, t)` | triangles `t[3 ne]` |
| `webgl_points (S, name, np, xy)` | points |
| `webgl_field (S, layer, name, where, v, categorical)` | a field on a layer |
| `webgl_scene_write (S, path)` | the HTML file |
| `webgl_offline ()` | 1 if three.js is in the pages, 0 if from jsDelivr |
| `webgl_tree (S, name, n, s)` | libb-amr's leaves and scalars (`webgl-tree.h`) |
| `webgl_hmesh (S, name, M)` | a libhmesh manifold (`webgl-hmesh.h`) |
| `webgl_bamr (S, path)` | a `.bamr` file (`webgl-bamr.h`) |
| `bamr2html FILE.bamr [OUT.html]` | the command (`tools/bamr2html.c`) |

Calls return -1 for what they cannot use, and change nothing then.

## Building

```sh
make                 # libwebgl.a: webgl.c, the adapters, webgl-assets.c
make bamr2html       # the command: plain C, no libb-amr
sh vendor/fetch.sh   # first, for offline pages
```

The adapters `webgl-tree.c` and `webgl-hmesh.c` are built with cwrap, as
libb-amr and libhmesh are; `webgl.c` and `webgl-bamr.c` are plain C.

## Tests

```sh
sh tests/run.sh      # from libwebgl
```

All under AddressSanitizer and UBSan, failing on any warning:

- **The scene** (`test_webgl.c`, plain C): every kind of layer, far from the
  origin, with fields of every kind (cells with NaN and inf, a categorical
  one); eight calls refused and nothing added by them; a title and a name
  made to break the HTML (`</script>`, quotes, `&`), escaped; as many
  `<script>` opened as closed; every array decoded by the test's own base64
  decoder and compared bit for bit; three.js and the viewer in the file
  byte for byte (56 checks); again built online: the jsDelivr URLs,
  nothing embedded (52 checks).
- **The adapters** (`test_scene.c`, cwrap, double and float): a tree refined
  over half its domain, moved and resized; a circle from libhmesh with
  points removed (holes in the ids), a scalar on its points, a vector on its
  edges. Boxes as the leaves' indices give them; each cell the function at
  its centre; every point and segment as the mesh has it (16 checks).
- **The .bamr reader** (`test_bamr.c`, cwrap, double and float): a file
  written by libb-amr's `tree_output ()`, read by `webgl_bamr ()`: boxes,
  levels, names, every cell of two scalars bit for bit against the blocks;
  a missing file and seven broken copies refused with nothing added (18
  checks). Then `bamr2html`'s page against the API's, byte for byte.
- **The JSON** (`check_json.py`): manifest and import map by Python's
  parser, every array named once.
- **viewer.js** parsed by `gjs` (SpiderMonkey, Firefox's engine) when it is
  there: syntax only. The drawing is checked by opening
  `tests/test_scene.html` in a browser.

Bugs made on purpose, each caught (checks failed):

| Bug | test | failed |
|---|---|---|
| base64: one byte left, one `=` only | test_webgl, test_scene | 3, 4 |
| boxes: the corner not subtracted | test_webgl, test_scene | 1, 1 |
| points: the x offset for y too | test_webgl, test_scene | 3, 1 |
| names: `<` not escaped | test_webgl | 2 |
| the title's mark not filled | test_webgl | 2 |
| a point field allowed on blocks | test_webgl | 12 |
| cells transposed (j fastest) | test_scene | 1 |
| a box's y from i (the tree) | test_scene | 1 |
| a segment's points by id, not renumbered | test_scene | 4 |
| a vector's components all the first | test_scene | 1 |
| .bamr: an index out of its level accepted | test_bamr | 3 |
| .bamr: a level past 30 accepted | test_bamr | 3 |
| .bamr: the size not checked | test_bamr | crash |
| .bamr: float cells read as doubles | test_bamr (float) | crash |
| .bamr: a box's y from i | test_bamr | 1 |
| .bamr: another dimension accepted | test_bamr | 2 |

## Next

- **3D:** a slice plane through the blocks (a 3D texture a leaf), a surface
  mesh whole or clipped, orbiting; `.bamr` files of 3D.
- **The mesh saved too:** a file of libhmesh's points and simplices beside
  the `.bamr`, so a saved run shows grid and interface together.
- **Time series:** several snapshots in one page, a slider.
- **Live:** the same viewer fed through a WebSocket during a run (a
  minimal server of our own, or libws: wsServer, its origin, is GPL-3.0).
- **MPI from a program:** the processes' parts gathered on rank 0
  (`webgl_tree` draws one process's leaves and halo today; a `.bamr` file
  needs none of it).
- **Many blocks:** one texture atlas instead of a texture a block, past
  about 10⁴ leaves.
- **Edge length against the cell** (libhmesh's design, (E)) as a field of
  the mesh.
