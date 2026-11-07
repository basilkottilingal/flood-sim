# Physics-Informed Neural Networks (PINN) in C

(on going project)

It is not a good idea to use C for NNs. But since we are giving an attempt on C wrapper,
it might give a window to simplify an attempt. However you might have to rely on some python libraries 
for training. This is only an attempt. 

## What it is

A C interface to a network trained by PyTorch or JAX, fed with the cells of
a libb-amr quadtree: `u = u (x, y)` (or `(x, y, t)`) from the solver
reinforces the same network, fit after fit, and the network's values come
back into scalars.

    C program ── pinn-bamr.h ──> pinn.h ──> embedded CPython ──> pinn_backend.py
    (libb-amr)   cells -> batch   memoryviews,      GIL          ├─ pinn_torch.py
                 batch -> scalars  no copy                        └─ pinn_jax.py

Two layers, so that cwrap never sees `Python.h`:

| file | built by | knows |
|---|---|---|
| `pinn.h`, `pinn.c` | `cc` | Python, not the tree |
| `pinn-bamr.h`, `pinn-bamr.c` | `cwrap -dimension=D` | the tree, not Python |
| `python/pinn_backend.py` | | buffers, weights, scaling, both backends |
| `python/pinn_torch.py`, `python/pinn_jax.py` | | the network, the loss, Adam |

## Use

```c
#include "pinn-bamr.h"

PinnGrid g = pinn_grid ();                 /* the box of tree_shift_resize () */
g.x0 [0] = x0, g.x0 [1] = y0, g.h0 = h0;
g.sample = PINN_SAMPLE_GAUSS2;             /* fit cell averages, O (d^4)     */
g.weight = PINN_WEIGHT_AREA;               /* or _UNIFORM: fronts count most */

PinnConfig cfg = pinn_grid_config (&g, "torch", 3, 0, 0);   /* or "jax" */
cfg.residual = "swe:residual";             /* optional: module:function     */
Pinn * p = pinn_new (&cfg);

PinnBatch b = {0};
for (step ...) {
  solver_step ();                           /* adapt_wavelet () etc          */
  pinn_bamr_gather (&g, 3, (int []) { eta, hu, hv }, &b);
  pinn_bamr_collocation (&g, 16, step, &b); /* points crowd where refined    */
  PinnLoss l;
  if (pinn_fit (p, &b, 200, &l) < 0)        /* warm start: reinforces        */
    fprintf (stderr, "%s\n", pinn_error (p));
}
pinn_bamr_predict (p, &g, 3, (int []) { pe, phu, phv });   /* back to scalars */
pinn_bamr_residual (p, &g, 3, (int []) { r1, r2, r3 });   /* where it breaks */
pinn_save (p, "pinn.ckpt");
pinn_batch_free (&b), pinn_free (p), pinn_shutdown ();
```

Residuals are written in Python, in physical coordinates (the network maps
its inputs to [-1, 1] inside, so the chain rule is free):

```python
def residual (model, x):            # torch: x (N, nin), requires_grad
  u = model (x)                     # (N, nout)
  ...torch.autograd.grad (..., create_graph=True)...
  return r                          # (N, nout)

def residual (u, x):                # jax: one point x (nin,), u (x) -> (nout,)
  return jax.jacfwd (u) (x) [:, 0]  # (nout,); vmapped by the backend
```

See `tests/pinn_test_residual.py` for both.

## Notes

- Cell values are finite volume averages: `PINN_SAMPLE_GAUSS2` compares the
  mean of the network over the 2^D Gauss points of a cell with its value;
  `PINN_SAMPLE_CENTRE` uses the centre, `O (d^2)` off where `u` curves.
- The scaling of the outputs is set by the first fit with data, and kept,
  so that warm starts stay consistent.
- JAX: batches are padded to a power of 2 (weight 0), so the jitted step
  compiles again only when the cell count doubles.
- The module directory is compiled in (`PINN_PYDIR`); the environment
  variable `PINN_PYDIR` overrides it. torch and jax come from
  `PYTHONPATH` or the system's site-packages.
- Calls hold the GIL: safe from any thread, one at a time.
- `tree.c` keeps the box `x0, h0` static, so `PinnGrid` repeats it. A
  `tree_get_box ()` in libb-amr would remove that.

## Build and test

    make CWRAP=../libc-pp/build/cwrap DIMENSION=2   # libpinn.a
    make print-ldlibs                               # -lpython3.10 ... to link

    PYTHONPATH=/where/torch/and/jax/are sh tests/run.sh

`tests/run.sh` runs `tests/test_pinn.c` in 1D and 2D, double and float, with
each backend Python can import, and skips the others.
