"""
JAX backend of libpinn: see pinn_backend.py.

The train step is jit compiled for each shape of batch. As the tree adapts,
the number of cells changes every snapshot, so batches are padded to the
next power of 2 with weight 0: a new compile only when the size doubles.
Adam is written here (a few lines) so that optax is not needed.
"""

import math
import pickle

import numpy as np
import jax
import jax.numpy as jnp

from pinn_backend import Base


_ACT = {"tanh": jnp.tanh, "silu": jax.nn.silu, "sin": jnp.sin}


def _pow2 (n):
  return 1 if n <= 1 else 1 << (n - 1).bit_length ()


def _pad (a, n):
  """a padded along axis 0 to n rows, by repeating its first row"""
  if a.shape [0] == n:
    return a
  pad = np.repeat (a [:1], n - a.shape [0], axis=0)
  return np.concatenate ([a, pad], axis=0)


class JaxModel (Base):

  def __init__ (self, cfg):
    super ().__init__ (cfg)
    if self.dtype == "float64":
      jax.config.update ("jax_enable_x64", True)
    self.jdtype = jnp.float64 if self.dtype == "float64" else jnp.float32
    dev = (cfg.get ("device") or "cpu").split (":")
    devs = jax.devices (dev [0] if dev [0] != "cuda" else "gpu")
    self.device = devs [int (dev [1]) if len (dev) > 1 else 0]
    self.act = _ACT [cfg.get ("activation") or "tanh"]
    key = jax.random.PRNGKey (cfg.get ("seed", 0))
    n = 2 * self.B.shape [1] if self.B is not None else self.nin
    sizes = [n] + [cfg ["width"]] * cfg ["depth"] + [self.nout]
    params = []
    for a, b in zip (sizes [:-1], sizes [1:]):
      key, k = jax.random.split (key)
      lim = math.sqrt (6.0 / (a + b))                  # Glorot uniform
      W = jax.random.uniform (k, (a, b), self.jdtype, -lim, lim)
      params.append ((W, jnp.zeros ((b,), self.jdtype)))
    self.params = self._put (params)
    self.lr = cfg ["lr"]
    self.adam = self._put ({"m": jax.tree_util.tree_map (jnp.zeros_like, params),
                            "v": jax.tree_util.tree_map (jnp.zeros_like, params),
                            "t": jnp.zeros ((), self.jdtype)})
    self._build ()

  def _put (self, tree):
    return jax.device_put (jax.tree_util.tree_map (
      lambda a: jnp.asarray (a, self.jdtype), tree), self.device)

  def on_scaling (self):
    self._build ()

  def _build (self):
    """the jitted functions, closing over the constants of the model"""
    c = lambda a: jnp.asarray (a, self.jdtype)
    lo, hi = c (self.lo), c (self.hi)
    shift, scale, inv = c (self.shift), c (self.scale), c (1 / self.scale)
    B = None if self.B is None else c (self.B)
    act, nout = self.act, self.nout
    ld, lp, lr = self.lambda_data, self.lambda_phys, self.lr
    rfn = self.residual_fn

    def net (params, x):                  # one point (nin,) -> (nout,)
      z = 2 * (x - lo) / (hi - lo) - 1
      if B is not None:
        zb = 2 * jnp.pi * z @ B
        z = jnp.concatenate ([jnp.sin (zb), jnp.cos (zb)])
      for W, b in params [:-1]:
        z = act (z @ W + b)
      W, b = params [-1]
      return shift + scale * (z @ W + b)

    batched = jax.vmap (net, in_axes=(None, 0))

    def resid (params, x):                # (N, nin) -> (N, nout)
      return jax.vmap (lambda p: rfn (lambda y: net (params, y), p)) (x)

    def loss (params, x, u, w, xc, wc):
      data = phys = jnp.zeros ((), self.jdtype)
      if x is not None:
        n, nq, nin = x.shape
        pred = batched (params, x.reshape (-1, nin)).reshape (n, nq, nout)
        data = (w * (((pred.mean (axis=1) - u) * inv) ** 2).mean (axis=1)).sum ()
      if xc is not None:
        phys = (wc * ((resid (params, xc) * inv) ** 2).mean (axis=1)).sum ()
      return ld * data + lp * phys, (data, phys)

    def step (params, adam, x, u, w, xc, wc):
      (l, (d, p)), g = jax.value_and_grad (loss, has_aux=True) (
        params, x, u, w, xc, wc)
      b1, b2, eps = 0.9, 0.999, 1e-8
      t = adam ["t"] + 1
      m = jax.tree_util.tree_map (lambda m, g: b1 * m + (1 - b1) * g,
                                  adam ["m"], g)
      v = jax.tree_util.tree_map (lambda v, g: b2 * v + (1 - b2) * g * g,
                                  adam ["v"], g)
      a = lr * jnp.sqrt (1 - b2 ** t) / (1 - b1 ** t)
      params = jax.tree_util.tree_map (
        lambda q, m, v: q - a * m / (jnp.sqrt (v) + eps), params, m, v)
      return params, {"m": m, "v": v, "t": t}, l, d, p

    self._step = jax.jit (step)
    self._fwd = jax.jit (batched)
    self._res = jax.jit (resid) if rfn else None

  def train (self, x, u, w, xc, wc, steps):
    put = lambda a: None if a is None else \
      jax.device_put (jnp.asarray (a, self.jdtype), self.device)
    if x is not None:
      N = _pow2 (x.shape [0])
      x, u, w = _pad (x, N), _pad (u, N), np.pad (w, (0, N - len (w)))
    if xc is not None:
      N = _pow2 (xc.shape [0])
      xc, wc = _pad (xc, N), np.pad (wc, (0, N - len (wc)))
    x, u, w, xc, wc = map (put, (x, u, w, xc, wc))
    params, adam = self.params, self.adam
    for _ in range (steps):
      params, adam, l, d, p = self._step (params, adam, x, u, w, xc, wc)
    self.params, self.adam = params, adam
    return float (l), float (d), float (p)

  def _eval (self, f, x):
    n = x.shape [0]
    xp = jnp.asarray (_pad (x, _pow2 (n)), self.jdtype)
    return np.asarray (f (self.params, xp)) [:n]

  def forward_np (self, x):
    return self._eval (self._fwd, x)

  def residual_np (self, x):
    return self._eval (self._res, x)

  def save (self, path):
    host = lambda t: jax.tree_util.tree_map (np.asarray, t)
    with open (path, "wb") as f:
      pickle.dump ({"params": host (self.params), "adam": host (self.adam),
                    "common": self.state_common ()}, f)

  def load (self, path):
    with open (path, "rb") as f:
      s = pickle.load (f)
    self.params = self._put ([tuple (l) for l in s ["params"]])
    self.adam = self._put (s ["adam"])
    self.load_common (s ["common"])
