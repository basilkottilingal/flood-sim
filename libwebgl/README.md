# libwebgl

A scene for the browser: the grid of libb-amr and the manifold of libhmesh,
or any blocks, lines, triangles and points with fields on them, written as
**one HTML file** that draws them with three.js (WebGL). Nothing to install:
the viewer and the data are inside the file, and three.js too when the
library was built with it (`sh vendor/fetch.sh`, then the pages need no
network); otherwise the pages load three.js r160 from jsDelivr. Firefox
opens them from disk. 2D for now.

```c
#include "webgl-tree.h"
#include "webgl-hmesh.h"

WebglScene * S = webgl_scene_new ("step 120");
webgl_tree (S, "grid", 2, (scalar []) { c, h });   /* the leaves, c and h */
webgl_hmesh (S, "interface", M);                   /* points, edges, fields */
webgl_scene_write (S, "step-120.html");
webgl_scene_free (S);
```

## In the browser

- **Pan** by dragging, **zoom** with the wheel (at the cursor), **r** resets
  the view.
- **A panel per layer:** shown or hidden; the field drawn; the colour map
  (viridis, coolwarm, gray); the range, typed or `auto`; the blocks'
  outlines, a layer's points, a triangulation's edges.
- **Click to pick:** a cell (its block, indices, box and every field's value
  there), a segment or triangle (its points, their coordinates, the fields),
  a point. Coordinates are given back in double, the scene's origin added.
- **Categorical fields** (levels, ranks) get distinct colours (Okabe-Ito),
  the others a colour map over their range; NaN is grey.
- Errors show on the page itself, in red, and so does a viewer that did not
  start (an old browser, no WebGL2).

## How

- **The scene** (`webgl.h`, plain C) copies what it is given. Layers:
  blocks of m x m cells (each a box), lines, triangles, points. Fields live
  on cells, blocks, points or elements, single precision.
- **Positions** are stored as offsets from the scene's lower corner, in
  single precision: far from the origin a float would lose the cell.
- **The file** is the template `viewer/viewer.html` with four marks filled:
  the title (escaped); the import map, with three.js and OrbitControls as
  `data:` URLs, so the viewer's imports resolve offline; the data, a JSON
  manifest (layers, fields, ranges, blob ids) and one
  `<script type="application/octet-stream">` a typed array, base64,
  little-endian; and `viewer/viewer.js`, inline as a module. Names in the
  manifest have `<` escaped, so none ends its script.
- **Compiled in:** `gen-assets.sh` turns the template, the viewer and
  three.js into C arrays (`xxd -i`, `webgl-assets.c`, made by the Makefile
  and the tests, not kept in git), so a program needs no file at run time.
  A file is about 1 MB with no data (three.js, base64).
- **The viewer** (three.js r160): an orthographic camera and OrbitControls
  with rotation off (MapControls' settings). A block is a quad with a float
  texture of its cells, coloured by a shader (nearest, so cells stay
  sharp); lines, triangles and points are geometry with colours per vertex.
  three.js's colour management is off, so the shader's colours and the
  vertices' agree.
- **The adapters**, built with cwrap: `webgl-tree.c` reads the leaves'
  blocks through `tree_caches ()` and `block_at ()` (no lookup in the level
  index), and adds each leaf's level and owner; `webgl-hmesh.c` renumbers
  the points (holes left out), and adds every `HmReal` field of the points
  and of the top level (a vector as `name.x`, `name.y`), the ids, and each
  simplex's length or area.

## API

| | |
|---|---|
| `webgl_scene_new (title)`, `webgl_scene_free (S)` | a scene |
| `webgl_blocks (S, name, n, m, box)` | n blocks of m x m cells, `box[4 n]` = x0, y0, width, height; the layer's number |
| `webgl_lines (S, name, np, xy, ne, e)` | points `xy[2 np]`, segments `e[2 ne]` |
| `webgl_triangles (S, name, np, xy, ne, t)` | triangles `t[3 ne]` |
| `webgl_points (S, name, np, xy)` | points |
| `webgl_field (S, layer, name, where, v, categorical)` | a field: `WEBGL_CELL` (n m^2, block by block, i fastest), `WEBGL_BLOCK`, `WEBGL_POINT`, `WEBGL_ELEMENT` |
| `webgl_scene_write (S, path)` | the HTML file |
| `webgl_offline ()` | 1 if three.js is in the pages, 0 if they load it from jsDelivr |
| `webgl_tree (S, name, n, s)` | the tree's leaves, a field a scalar, level and pid (`webgl-tree.h`) |
| `webgl_hmesh (S, name, M)` | a manifold: points, lines or triangles, its fields (`webgl-hmesh.h`) |

Calls return -1 for what they cannot use, and change nothing then: a field
of the wrong kind for its layer, an index out of range, a NULL.

## Tests

```sh
sh tests/run.sh       # from libwebgl
```

- **The scene** (`test_webgl.c`, plain C): every kind of layer, far from the
  origin, with fields of every kind (cells with NaN and inf, a categorical
  one); the calls refused (eight of them) and nothing added by them; a title
  and a name made to break the HTML (`</script>`, quotes, `&`), escaped; as
  many `<script>` opened as closed; every blob decoded by the test's own
  base64 decoder and compared bit for bit; three.js, OrbitControls and the
  viewer in the file, byte for byte (56 checks); and again built without
  three.js: the import map's jsDelivr URLs, nothing embedded (52 checks).
- **The adapters** (`test_scene.c`, cwrap, 2D, double and float): a tree
  refined over half its domain, moved and resized, a scalar on it; a circle
  from libhmesh with a few points removed (holes in the ids), a scalar on
  its points, a vector on its edges. Each leaf's box as its `Index` gives
  it; each cell's value the function at the cell's centre (so the cells'
  order is right); levels and ranks; every point where `hm_x ()` has it;
  every segment between the points `hm_vertex ()` names; lengths; both
  components of the vector. 16 checks.
- **The JSON** (`check_json.py`): the manifest and the import map parsed by
  Python's parser, every blob named once.
- **viewer.js** parsed as a module by `gjs` (SpiderMonkey, Firefox's engine)
  when it is there: syntax only. The drawing itself is checked by opening
  `tests/test_scene.html` in a browser: Firefox here is a snap, which does
  not run inside Claude's sandbox.

All under AddressSanitizer and UBSan, failing on any warning. Bugs made on
purpose, each caught (checks failed):

| Bug | test_webgl | test_scene |
|---|---|---|
| base64: one byte left, one `=` only | 3 | 4 |
| boxes: the corner not subtracted | 1 | 1 |
| points: the x offset for y too | 3 | 1 |
| names: `<` not escaped | 2 | |
| the title's mark not filled | 2 | |
| a point field allowed on blocks | 12 | |
| cells transposed (j fastest) | | 1 |
| a box's y from i | | 1 |
| a segment's points by id, not renumbered | | 4 |
| a vector's components all the first | | 1 |

## three.js: embedded or online

three.js is not kept in git (700 KB). `sh vendor/fetch.sh` gets r160
(`three@0.160.0`) from jsDelivr into `vendor/` and checks the files against
the sums below; built then, the library embeds them and its pages work
offline. Built without them (or with `WEBGL_ONLINE=1 sh gen-assets.sh`),
the pages' import map points at the same files on jsDelivr, and opening a
page needs the network. `webgl_offline ()` tells which a build does; the
page says so in red if it cannot start. The licence, MIT, is kept
(`vendor/three.LICENSE`). The first copy came from the npm registry, its
SHA-512 integrity checked; jsDelivr's files are the same, byte for byte:

| file | sha256 |
|---|---|
| `three.module.min.js` | `3e690ac7d180b0aadf0891bea39eec643e29e2d3e75c99b18689518665f69ba6` |
| `OrbitControls.js` | `5a44a9e86a2a0fb11933eed69bc2cd33c76a496854c1aed6ed776efa87d7b064` |
| `three.LICENSE` | `852e0e8699169bf9f6fdc6bda3e682d078dcbc738b5d33e74df594721bff271d` |

## Next

- **3D**: a slice plane through the blocks (a 3D texture a leaf), a surface
  mesh drawn whole or clipped, orbiting.
- **MPI**: the processes' parts gathered on rank 0, or one file a rank; the
  mesh's owners and ghosts, the halves cut by the partition marked.
- **Time series**: several snapshots in one file, a slider.
- **Live**: the same viewer fed through a WebSocket (libws) during a run.
- **Many blocks**: one texture atlas instead of a texture a block, past
  about 10^4 leaves.
- **Edge length against the cell** ((E) of libhmesh's design) as a field of
  the mesh, from the tree's cell at each edge.
