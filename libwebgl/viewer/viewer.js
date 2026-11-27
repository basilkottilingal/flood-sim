/*
.. libwebgl's viewer: the scene in the page's manifest and blobs, drawn
.. with three.js, in 2D: an orthographic camera, pan and zoom, layers that
.. can be shown or hidden, a field and a colour map a layer, and picking.
..
.. Positions in the blobs are offsets from manifest.origin, in single
.. precision; picking adds the origin back in double.
..
.. Colours are written as given: three.js's colour management is off and
.. the output is not converted, so the shader's colours and the lines'
.. vertex colours agree.
*/

import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

THREE.ColorManagement.enabled = false;

const manifest = JSON.parse (document.getElementById ('manifest').textContent);
document.getElementById ('title').textContent = manifest.title || 'scene';

/* a blob as a typed array */
function decode (id, Type) {
  const s = atob (document.getElementById (id).textContent);
  const bytes = new Uint8Array (s.length);
  for (let i = 0; i < s.length; i++)
    bytes[i] = s.charCodeAt (i);
  return new Type (bytes.buffer);
}

/* -------- colours -------- */

const MAPS = {
  viridis: ['#440154', '#472d7b', '#3b528b', '#2c728e', '#21908c',
            '#27ad81', '#5dc863', '#aadc32', '#fde725'],
  coolwarm: ['#3b4cc0', '#dddddd', '#b40426'],
  gray: ['#000000', '#ffffff'],
};
/* Okabe-Ito: levels, ranks */
const PALETTE = ['#0072B2', '#E69F00', '#009E73', '#CC79A7', '#56B4E9',
                 '#D55E00', '#F0E442', '#999999'];
const NAN_RGB = [0.45, 0.45, 0.45];

function hex (h) {
  const n = parseInt (h.slice (1), 16);
  return [(n >> 16 & 255) / 255, (n >> 8 & 255) / 255, (n & 255) / 255];
}

/* 256 colours of a map, linear between its anchors */
const LUTS = {};
for (const name in MAPS) {
  const a = MAPS[name].map (hex), lut = new Float32Array (256 * 3);
  for (let i = 0; i < 256; i++) {
    const t = i / 255 * (a.length - 1), k = Math.min (Math.floor (t),
                                                      a.length - 2);
    const f = t - k;
    for (let c = 0; c < 3; c++)
      lut[3 * i + c] = a[k][c] * (1 - f) + a[k + 1][c] * f;
  }
  LUTS[name] = lut;
}
const PALETTE_RGB = PALETTE.map (hex);

function lutTexture (name) {
  const lut = LUTS[name], data = new Uint8Array (256 * 4);
  for (let i = 0; i < 256; i++) {
    for (let c = 0; c < 3; c++)
      data[4 * i + c] = Math.round (lut[3 * i + c] * 255);
    data[4 * i + 3] = 255;
  }
  const t = new THREE.DataTexture (data, 256, 1, THREE.RGBAFormat);
  t.magFilter = t.minFilter = THREE.LinearFilter;
  t.needsUpdate = true;
  return t;
}

/* the colour of value v under a layer's colouring */
function colour (v, c) {
  if (!Number.isFinite (v))
    return NAN_RGB;
  if (c.categorical)
    return PALETTE_RGB[((Math.round (v) % 8) + 8) % 8];
  const t = Math.min (1, Math.max (0, (v - c.min) / ((c.max - c.min) ||
                                                     1e-30)));
  const i = Math.round (t * 255), lut = LUTS[c.map];
  return [lut[3 * i], lut[3 * i + 1], lut[3 * i + 2]];
}

function barCanvas (c, canvas) {
  const w = canvas.width = 256, ctx = canvas.getContext ('2d');
  canvas.height = 1;
  for (let i = 0; i < w; i++) {
    const rgb = c.categorical ? PALETTE_RGB[Math.floor (i / 32) % 8]
      : LUTS[c.map].slice (3 * i, 3 * i + 3);
    ctx.fillStyle = `rgb(${rgb.map (x => Math.round (x * 255)).join (',')})`;
    ctx.fillRect (i, 0, 1, 1);
  }
}

/* -------- renderer, camera, controls -------- */

const view = document.getElementById ('view');
const renderer = new THREE.WebGLRenderer ({ antialias: true,
                                            preserveDrawingBuffer: true });
renderer.outputColorSpace = THREE.LinearSRGBColorSpace;
renderer.setPixelRatio (window.devicePixelRatio || 1);
view.appendChild (renderer.domElement);
const scene = new THREE.Scene ();
scene.background = new THREE.Color (0x1b1d22);

const W = manifest.size[0] || 1, H = manifest.size[1] || 1;
const camera = new THREE.OrthographicCamera (-1, 1, 1, -1, -100, 100);
const controls = new OrbitControls (camera, renderer.domElement);
controls.enableRotate = false;
controls.screenSpacePanning = true;
controls.zoomToCursor = true;
controls.mouseButtons = { LEFT: THREE.MOUSE.PAN, MIDDLE: THREE.MOUSE.DOLLY,
                          RIGHT: THREE.MOUSE.PAN };
controls.touches = { ONE: THREE.TOUCH.PAN, TWO: THREE.TOUCH.DOLLY_PAN };

/* the frustum for the window, the scene's height fitting at zoom 1 */
function frustum () {
  const w = view.clientWidth || 1, h = view.clientHeight || 1;
  const aspect = w / h;
  let fh = H * 1.08, fw = fh * aspect;
  if (fw < W * 1.08)
    fw = W * 1.08, fh = fw / aspect;
  camera.left = -fw / 2, camera.right = fw / 2;
  camera.top = fh / 2, camera.bottom = -fh / 2;
  camera.updateProjectionMatrix ();
  renderer.setSize (w, h);
}

function reset () {
  camera.position.set (W / 2, H / 2, 10);
  camera.zoom = 1;
  controls.target.set (W / 2, H / 2, 0);
  frustum ();
  controls.update ();
  render ();
}

function render () {
  renderer.render (scene, camera);
}
controls.addEventListener ('change', render);
window.addEventListener ('resize', () => { frustum (); render (); });
window.addEventListener ('keydown', e => { if (e.key === 'r') reset (); });

/* -------- layers -------- */

const VERTEX = `
varying vec2 vUv;
void main () {
  vUv = uv;
  gl_Position = projectionMatrix * modelViewMatrix * vec4 (position, 1.0);
}`;

const FRAGMENT = `
uniform sampler2D tex;
uniform sampler2D lut;
uniform float vmin, vmax;
uniform int categorical;
uniform vec3 palette[8];
varying vec2 vUv;
void main () {
  float v = texture2D (tex, vUv).r;
  if (!(v == v) || abs (v) > 3.0e38) {
    gl_FragColor = vec4 (0.45, 0.45, 0.45, 1.0);
    return;
  }
  if (categorical == 1) {
    int k = int (mod (floor (v + 0.5), 8.0));
    gl_FragColor = vec4 (palette[k], 1.0);
    return;
  }
  float t = clamp ((v - vmin) / max (vmax - vmin, 1e-30), 0.0, 1.0);
  gl_FragColor = vec4 (texture2D (lut, vec2 (t, 0.5)).rgb, 1.0);
}`;

const unitPlane = new THREE.PlaneGeometry (1, 1);

/* a blocks layer: a quad a block, a float texture of its cells */
function makeBlocks (L) {
  const box = decode (L.box, Float32Array), m = L.m;
  const fields = L.fields.map (f => ({ ...f, v: decode (f.data,
                                                         Float32Array) }));
  const uniforms = {
    lut: { value: lutTexture ('viridis') },
    vmin: { value: 0 }, vmax: { value: 1 }, categorical: { value: 0 },
    palette: { value: PALETTE_RGB.map (c => new THREE.Vector3 (...c)) },
  };
  const group = new THREE.Group (), textures = [], meshes = [];
  for (let b = 0; b < L.n; b++) {
    const data = new Float32Array (m * m);
    const tex = new THREE.DataTexture (data, m, m, THREE.RedFormat,
                                       THREE.FloatType);
    tex.magFilter = tex.minFilter = THREE.NearestFilter;
    tex.needsUpdate = true;
    const mat = new THREE.ShaderMaterial ({
      uniforms: { ...uniforms, tex: { value: tex } },
      vertexShader: VERTEX, fragmentShader: FRAGMENT });
    const mesh = new THREE.Mesh (unitPlane, mat);
    const [x0, y0, w, h] = box.subarray (4 * b, 4 * b + 4);
    mesh.scale.set (w, h, 1);
    mesh.position.set (x0 + w / 2, y0 + h / 2, 0);
    mesh.userData = { layer: L, block: b };
    group.add (mesh);
    textures.push (tex), meshes.push (mesh);
  }
  /* outlines */
  const p = new Float32Array (L.n * 8 * 3);
  for (let b = 0; b < L.n; b++) {
    const [x0, y0, w, h] = box.subarray (4 * b, 4 * b + 4);
    const c = [[x0, y0], [x0 + w, y0], [x0 + w, y0 + h], [x0, y0 + h]];
    for (let e = 0; e < 4; e++) {
      const a = c[e], d = c[(e + 1) % 4], o = 24 * b + 6 * e;
      p.set ([a[0], a[1], 0.0005, d[0], d[1], 0.0005], o);
    }
  }
  const og = new THREE.BufferGeometry ();
  og.setAttribute ('position', new THREE.BufferAttribute (p, 3));
  const outline = new THREE.LineSegments (og, new THREE.LineBasicMaterial (
    { color: 0x000000, transparent: true, opacity: 0.35 }));
  outline.raycast = () => {};
  group.add (outline);

  const state = { field: fields.findIndex (f => f.where === 'cell'),
                  map: 'viridis', min: 0, max: 1, categorical: 0 };
  if (state.field < 0)
    state.field = fields.length ? 0 : -1;
  function apply () {
    const f = fields[state.field];
    for (let b = 0; b < L.n; b++) {
      const data = textures[b].image.data;
      if (!f)
        data.fill (0);
      else if (f.where === 'cell')
        data.set (f.v.subarray (b * m * m, (b + 1) * m * m));
      else
        data.fill (f.v[b]);
      textures[b].needsUpdate = true;
    }
    uniforms.vmin.value = state.min, uniforms.vmax.value = state.max;
    uniforms.categorical.value = state.categorical ? 1 : 0;
    uniforms.lut.value = lutTexture (state.map);
  }
  function pick (hit) {
    const b = hit.object.userData.block, uv = hit.uv;
    const i = Math.min (m - 1, Math.floor (uv.x * m));
    const j = Math.min (m - 1, Math.floor (uv.y * m));
    const [x0, y0, w, h] = box.subarray (4 * b, 4 * b + 4);
    let s = `${L.name}: block ${b}, cell (${i}, ${j})\n` +
      `  cell box  [${fmt (manifest.origin[0] + x0 + i * w / m)}, ` +
      `${fmt (manifest.origin[0] + x0 + (i + 1) * w / m)}] x ` +
      `[${fmt (manifest.origin[1] + y0 + j * h / m)}, ` +
      `${fmt (manifest.origin[1] + y0 + (j + 1) * h / m)}]\n`;
    for (const f of fields)
      s += `  ${f.name.padEnd (10)} ${fmt (f.where === 'cell' ?
        f.v[b * m * m + j * m + i] : f.v[b])}\n`;
    return s;
  }
  return { L, group, fields, state, apply, pick, extras: [
    { name: 'outlines', object: outline, on: true }] };
}

/* lines, triangles or points: vertex colours from a point or element
   field */
function makeGeometry (L) {
  const xy = decode (L.xy, Float32Array);
  const k = L.type === 'lines' ? 2 : L.type === 'triangles' ? 3 : 1;
  const e = L.type === 'points' ? null : decode (L.e, Uint32Array);
  const fields = L.fields.map (f => ({ ...f, v: decode (f.data,
                                                         Float32Array) }));
  const nv = L.type === 'points' ? L.np : L.ne * k;
  const z = L.type === 'triangles' ? 0.001 : L.type === 'lines' ? 0.002
    : 0.003;
  const pos = new Float32Array (nv * 3), col = new Float32Array (nv * 3);
  for (let v = 0; v < nv; v++) {
    const p = e ? e[v] : v;
    pos.set ([xy[2 * p], xy[2 * p + 1], z], 3 * v);
  }
  const g = new THREE.BufferGeometry ();
  g.setAttribute ('position', new THREE.BufferAttribute (pos, 3));
  g.setAttribute ('color', new THREE.BufferAttribute (col, 3));
  const group = new THREE.Group ();
  let object;
  if (L.type === 'lines')
    object = new THREE.LineSegments (g, new THREE.LineBasicMaterial (
      { vertexColors: true }));
  else if (L.type === 'triangles')
    object = new THREE.Mesh (g, new THREE.MeshBasicMaterial (
      { vertexColors: true, side: THREE.DoubleSide }));
  else
    object = new THREE.Points (g, new THREE.PointsMaterial (
      { vertexColors: true, size: 5, sizeAttenuation: false }));
  object.userData = { layer: L };
  group.add (object);
  const extras = [];
  if (L.type !== 'points') {
    /* the points of the layer, and a triangle's edges */
    const pg = new THREE.BufferGeometry ();
    const pp = new Float32Array (L.np * 3);
    for (let p = 0; p < L.np; p++)
      pp.set ([xy[2 * p], xy[2 * p + 1], 0.003], 3 * p);
    pg.setAttribute ('position', new THREE.BufferAttribute (pp, 3));
    const pts = new THREE.Points (pg, new THREE.PointsMaterial (
      { color: 0xffffff, size: 4, sizeAttenuation: false }));
    pts.raycast = () => {};
    group.add (pts);
    extras.push ({ name: 'points', object: pts, on: L.np < 20000 });
  }
  if (L.type === 'triangles') {
    const ep = new Float32Array (L.ne * 6 * 3);
    for (let t = 0; t < L.ne; t++)
      for (let s = 0; s < 3; s++) {
        const a = e[3 * t + s], b = e[3 * t + (s + 1) % 3];
        ep.set ([xy[2 * a], xy[2 * a + 1], 0.0015,
                 xy[2 * b], xy[2 * b + 1], 0.0015], 18 * t + 6 * s);
      }
    const eg = new THREE.BufferGeometry ();
    eg.setAttribute ('position', new THREE.BufferAttribute (ep, 3));
    const edges = new THREE.LineSegments (eg, new THREE.LineBasicMaterial (
      { color: 0x000000, transparent: true, opacity: 0.5 }));
    edges.raycast = () => {};
    group.add (edges);
    extras.push ({ name: 'edges', object: edges, on: true });
  }

  const state = { field: fields.length ? 0 : -1, map: 'viridis', min: 0,
                  max: 1, categorical: 0 };
  function apply () {
    const f = fields[state.field];
    for (let v = 0; v < nv; v++) {
      let rgb = [0.95, 0.95, 0.95];
      if (f)
        rgb = colour (f.where === 'point' ? f.v[e ? e[v] : v]
                      : f.v[Math.floor (v / k)], state);
      col.set (rgb, 3 * v);
    }
    g.attributes.color.needsUpdate = true;
  }
  function pick (hit) {
    let s;
    if (L.type === 'points') {
      const p = hit.index;
      s = `${L.name}: point ${p}\n` +
        `  x ${fmt (manifest.origin[0] + xy[2 * p])}` +
        `  y ${fmt (manifest.origin[1] + xy[2 * p + 1])}\n`;
      for (const f of fields)
        s += `  ${f.name.padEnd (10)} ${fmt (f.v[p])}\n`;
      return s;
    }
    const el = L.type === 'lines' ? Math.floor (hit.index / 2)
      : hit.faceIndex;
    s = `${L.name}: ${L.type === 'lines' ? 'segment' : 'triangle'} ${el}\n`;
    for (let c = 0; c < k; c++) {
      const p = e[k * el + c];
      s += `  point ${p}: (${fmt (manifest.origin[0] + xy[2 * p])}, ` +
        `${fmt (manifest.origin[1] + xy[2 * p + 1])})\n`;
    }
    for (const f of fields)
      s += f.where === 'element'
        ? `  ${f.name.padEnd (10)} ${fmt (f.v[el])}\n`
        : `  ${f.name.padEnd (10)} ` +
          Array.from ({ length: k }, (_, c) => fmt (f.v[e[k * el + c]]))
               .join ('  ') + '\n';
    return s;
  }
  return { L, group, fields, state, apply, pick, extras };
}

function fmt (x) {
  return Number.isFinite (x) ? (Math.abs (x) >= 1e5 || (x && Math.abs (x) <
    1e-3) ? x.toExponential (6) : x.toPrecision (8).replace (/\.?0+$/, ''))
    : String (x);
}

/* -------- the panel -------- */

function element (tag, attrs = {}, ...children) {
  const n = document.createElement (tag);
  for (const a in attrs)
    if (a === 'on')
      for (const ev in attrs.on)
        n.addEventListener (ev, attrs.on[ev]);
    else
      n[a] = attrs[a];
  for (const c of children)
    n.append (c);
  return n;
}

function panel (lay) {
  const { L, fields, state } = lay;
  const minIn = element ('input', { type: 'number', step: 'any' });
  const maxIn = element ('input', { type: 'number', step: 'any' });
  const bar = element ('canvas', { className: 'bar' });
  const mapSel = element ('select', {},
    ...Object.keys (MAPS).map (k => element ('option', { value: k,
                                                         textContent: k })));
  const fieldSel = element ('select', {},
    element ('option', { value: -1, textContent: '(none)' }),
    ...fields.map ((f, i) => element ('option', { value: i,
      textContent: `${f.name} (${f.where})` })));
  function auto () {
    const f = fields[state.field];
    state.categorical = f ? f.categorical : 0;
    state.min = f && f.min !== null ? f.min : 0;
    state.max = f && f.max !== null ? f.max : 1;
    minIn.value = state.min, maxIn.value = state.max;
  }
  function update () {
    lay.apply ();
    barCanvas (state, bar);
    render ();
  }
  fieldSel.value = state.field;
  fieldSel.addEventListener ('change', () => {
    state.field = +fieldSel.value;
    auto (), update ();
  });
  mapSel.addEventListener ('change', () => { state.map = mapSel.value;
                                             update (); });
  minIn.addEventListener ('change', () => { state.min = +minIn.value;
                                            update (); });
  maxIn.addEventListener ('change', () => { state.max = +maxIn.value;
                                            update (); });
  const show = element ('input', { type: 'checkbox', checked: true,
    on: { change: e => { lay.group.visible = e.target.checked;
                         render (); } } });
  const extras = lay.extras.map (x => {
    x.object.visible = x.on;
    return element ('label', {}, element ('input', {
      type: 'checkbox', checked: x.on,
      on: { change: e => { x.object.visible = e.target.checked;
                           render (); } } }), ' ' + x.name);
  });
  const div = element ('div', { className: 'layer' },
    element ('b', {}, element ('label', {}, show, ` ${L.name}`),
             element ('span', { className: 'hint',
               textContent: L.type === 'blocks'
                 ? `  ${L.n} blocks of ${L.m}²`
                 : `  ${L.np} points, ${L.ne} ${L.type === 'points' ? ''
                    : L.type}` })),
    element ('div', { className: 'row' }, element ('span', {}, 'field'),
             fieldSel),
    element ('div', { className: 'row' }, element ('span', {}, 'colours'),
             mapSel),
    element ('div', { className: 'row' }, element ('span', {}, 'range'),
             minIn, maxIn,
             element ('button', { textContent: 'auto',
                                  on: { click: () => { auto (); update (); }
                                      } })),
    bar,
    element ('div', { className: 'row' }, ...extras));
  auto ();
  update ();
  return div;
}

/* -------- start -------- */

const layers = manifest.layers.map (L => L.type === 'blocks'
                                    ? makeBlocks (L) : makeGeometry (L));
const list = document.getElementById ('layers');
for (const lay of layers) {
  scene.add (lay.group);
  list.append (panel (lay));
}

/* picking: the nearest object under the pointer (points and lines lie
   above triangles, triangles above blocks) */
const raycaster = new THREE.Raycaster ();
let down = null;
renderer.domElement.addEventListener ('pointerdown',
  e => { down = [e.clientX, e.clientY]; });
renderer.domElement.addEventListener ('pointerup', e => {
  if (!down || Math.hypot (e.clientX - down[0], e.clientY - down[1]) > 4)
    return;
  const r = renderer.domElement.getBoundingClientRect ();
  const ndc = new THREE.Vector2 ((e.clientX - r.left) / r.width * 2 - 1,
                                 -(e.clientY - r.top) / r.height * 2 + 1);
  const pixel = (camera.right - camera.left) / camera.zoom / r.width;
  raycaster.params.Line.threshold = 4 * pixel;
  raycaster.params.Points.threshold = 4 * pixel;
  raycaster.setFromCamera (ndc, camera);
  const targets = layers.filter (l => l.group.visible).map (l => l.group);
  const hits = raycaster.intersectObjects (targets, true)
    .filter (h => h.object.visible && h.object.userData.layer);
  const info = document.getElementById ('info');
  if (!hits.length) {
    const p = raycaster.ray.origin;
    info.textContent = `(${fmt (manifest.origin[0] + p.x)}, ` +
      `${fmt (manifest.origin[1] + p.y)}): nothing here`;
    return;
  }
  hits.sort ((a, b) => b.point.z - a.point.z || a.distance - b.distance);
  const h = hits[0], lay = layers.find (l => l.L === h.object.userData.layer);
  info.textContent = lay.pick (h);
});

reset ();
window.viewerReady = true;
