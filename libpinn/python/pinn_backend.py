"""
The Python side of libpinn (pinn.h): pinn.c calls create (cfg) once per
Pinn, and then the methods of the object it returns:

  fit (x, u, w, xc, wc, n, nq, nc, steps) -> (total, data, phys)
  predict (x, out, n)        writes out [n, nout]
  residual (x, out, n)       writes out [n, nout]
  save (path), load (path)

The arrays come as memoryviews of the C caller's memory (None for NULL):
x and xc of float64, the others of float32. They are only valid during the
call, so nothing here keeps them.

The network, in either backend:

  z = 2 (x - lo) / (hi - lo) - 1                   inputs into [-1, 1]
  z = [sin (2 pi z B), cos (2 pi z B)]             if fourier > 0
  y = MLP (z)                                      depth x width, act
  u = shift + scale * y                            the scaling of the
                                                   first fit, then kept

so that the network takes and gives physical values, and a residual
written in physical coordinates gets the chain rule for free.

Loss, with weights normalised to sum to 1, and errors in units of scale
so that fields of different magnitude count alike:

  data = sum_i w_i  mean_k ((mean_q u (x_iq) - u_i) / scale_k)^2
  phys = sum_j wc_j mean_k (r_k (xc_j) / scale_k)^2
  L    = lambda_data data + lambda_phys phys

Residual functions, named in the config as "module:function":

  torch   r = f (model, x)   x a (N, nin) tensor with requires_grad,
                             model (x) -> (N, nout); use
                             torch.autograd.grad (..., create_graph=True).
                             Returns (N, nout).
  jax     r = f (u, x)       x one point (nin,), u (x) -> (nout,); use
                             jax.jacfwd, jax.hessian. Returns (nout,).
                             The backend vmaps it over the points.
"""

import importlib
import math

import numpy as np


def _array (buf, dtype, shape):
  if buf is None:
    return None
  return np.frombuffer (buf, dtype=dtype).reshape (shape)


def _weights (w, n):
  """w normalised to sum to 1, or 1/n for None"""
  if n == 0:
    return np.zeros ((0,), np.float64)
  if w is None:
    return np.full ((n,), 1.0 / n)
  w = np.asarray (w, np.float64)
  s = w.sum ()
  if not s > 0:
    raise ValueError ("weights must have a positive sum")
  return w / s


def _resolve (name):
  """'module:function' -> the function"""
  if not name:
    return None
  mod, _, fn = name.partition (":")
  if not fn:
    raise ValueError ("residual must be 'module:function', not %r" % name)
  return getattr (importlib.import_module (mod), fn)


class Base:
  """what both backends share: buffers, weights, scaling"""

  def __init__ (self, cfg):
    self.cfg = cfg
    self.nin = cfg ["nin"]
    self.nout = cfg ["nout"]
    self.lo = np.asarray (cfg ["lo"], np.float64)
    self.hi = np.asarray (cfg ["hi"], np.float64)
    self.dtype = cfg.get ("dtype") or "float32"
    self.lambda_data = cfg.get ("lambda_data") or 1.0
    self.lambda_phys = cfg.get ("lambda_phys", 1.0)
    self.residual_fn = _resolve (cfg.get ("residual"))
    self.shift = np.zeros ((self.nout,))
    self.scale = np.ones ((self.nout,))
    self.scaled = False
    rng = np.random.default_rng (cfg.get ("seed", 0))
    nf = cfg.get ("fourier", 0)
    sigma = cfg.get ("fourier_sigma", 1.0)
    self.B = rng.normal (0.0, sigma, (self.nin, nf)) if nf > 0 else None

  def _set_scaling (self, u):
    """from the first batch with data; kept after, for warm starts"""
    if self.scaled or u is None or len (u) == 0:
      return
    self.shift = u.mean (axis=0).astype (np.float64)
    s = u.std (axis=0).astype (np.float64)
    self.scale = np.where (s > 1e-12 * (1 + np.abs (self.shift)), s, 1.0)
    self.scaled = True
    self.on_scaling ()

  def on_scaling (self):
    pass

  def fit (self, x, u, w, xc, wc, n, nq, nc, steps):
    x = _array (x, np.float64, (n, nq, self.nin)) if n else None
    u = _array (u, np.float32, (n, self.nout)) if n else None
    w = _weights (_array (w, np.float32, (n,)) if n else None, n)
    xc = _array (xc, np.float64, (nc, self.nin)) if nc else None
    wc = _weights (_array (wc, np.float32, (nc,)) if nc else None, nc)
    if nc and self.residual_fn is None:
      raise ValueError ("collocation points, but the config names no residual")
    self._set_scaling (u)
    if steps == 0:
      return (0.0, 0.0, 0.0)
    return tuple (float (v) for v in self.train (x, u, w, xc, wc, steps))

  def predict (self, x, out, n):
    x = _array (x, np.float64, (n, self.nin))
    np.frombuffer (out, dtype=np.float32).reshape (n, self.nout) [...] = \
      self.forward_np (x)

  def residual (self, x, out, n):
    if self.residual_fn is None:
      raise ValueError ("the config names no residual")
    x = _array (x, np.float64, (n, self.nin))
    np.frombuffer (out, dtype=np.float32).reshape (n, self.nout) [...] = \
      self.residual_np (x)

  def state_common (self):
    return {"shift": self.shift, "scale": self.scale, "scaled": self.scaled,
            "B": self.B}

  def load_common (self, s):
    self.shift, self.scale = np.asarray (s ["shift"]), np.asarray (s ["scale"])
    self.scaled, self.B = s ["scaled"], s ["B"]
    self.on_scaling ()


def create (cfg):
  backend = (cfg.get ("backend") or "").lower ()
  if backend in ("torch", "pytorch"):
    from pinn_torch import TorchModel
    return TorchModel (cfg)
  if backend == "jax":
    from pinn_jax import JaxModel
    return JaxModel (cfg)
  raise ValueError ("unknown backend %r: 'torch' or 'jax'" % backend)
