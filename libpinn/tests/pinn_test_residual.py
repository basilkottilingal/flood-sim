"""
Residuals for tests/test_pinn.c: r_k = d u_k / dx, the derivative along the
first input, in physical coordinates. The test compares it with finite
differences of predictions, which checks the autodiff path and the chain
rule through the normalisation of the inputs.
"""


def dudx_torch (model, x):
  import torch
  u = model (x)
  cols = []
  for k in range (u.shape [1]):
    g, = torch.autograd.grad (u [:, k].sum (), x, create_graph=True)
    cols.append (g [:, 0])
  return torch.stack (cols, dim=1)


def dudx_jax (u, x):
  import jax
  return jax.jacfwd (u) (x) [:, 0]
