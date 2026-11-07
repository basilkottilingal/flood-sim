"""PyTorch backend of libpinn: see pinn_backend.py"""

import math

import numpy as np
import torch

from pinn_backend import Base


class Sin (torch.nn.Module):
  def forward (self, x):
    return torch.sin (x)


_ACT = {"tanh": torch.nn.Tanh, "silu": torch.nn.SiLU, "sin": Sin}


class Net (torch.nn.Module):
  """physical x -> physical u: see pinn_backend.py"""

  def __init__ (self, cfg, B, dtype):
    super ().__init__ ()
    nin, nout = cfg ["nin"], cfg ["nout"]
    act = _ACT [cfg.get ("activation") or "tanh"]
    self.register_buffer ("lo", torch.tensor (cfg ["lo"], dtype=dtype))
    self.register_buffer ("hi", torch.tensor (cfg ["hi"], dtype=dtype))
    self.register_buffer ("shift", torch.zeros (nout, dtype=dtype))
    self.register_buffer ("scale", torch.ones (nout, dtype=dtype))
    if B is not None:
      self.register_buffer ("B", torch.tensor (B, dtype=dtype))
      n = 2 * B.shape [1]
    else:
      self.B = None
      n = nin
    layers = []
    for _ in range (cfg ["depth"]):
      layers += [torch.nn.Linear (n, cfg ["width"]), act ()]
      n = cfg ["width"]
    layers.append (torch.nn.Linear (n, nout))
    self.mlp = torch.nn.Sequential (*layers)

  def forward (self, x):
    z = 2 * (x - self.lo) / (self.hi - self.lo) - 1
    if self.B is not None:
      zb = 2 * math.pi * z @ self.B
      z = torch.cat ([torch.sin (zb), torch.cos (zb)], dim=-1)
    return self.shift + self.scale * self.mlp (z)


class TorchModel (Base):

  def __init__ (self, cfg):
    super ().__init__ (cfg)
    torch.manual_seed (cfg.get ("seed", 0))
    self.tdtype = torch.float64 if self.dtype == "float64" else torch.float32
    self.device = torch.device (cfg.get ("device") or "cpu")
    self.net = Net (cfg, self.B, self.tdtype).to (self.device)
    self.opt = torch.optim.Adam (self.net.parameters (), lr=cfg ["lr"])

  def _t (self, a):
    # a copy: the arrays are views of C memory, read only and valid only
    # during the call, which a tensor must not keep
    return torch.tensor (np.array (a), dtype=self.tdtype, device=self.device)

  def on_scaling (self):
    with torch.no_grad ():
      self.net.shift.copy_ (self._t (self.shift))
      self.net.scale.copy_ (self._t (self.scale))

  def _residual (self, x):
    x = x.detach ().requires_grad_ (True)
    return self.residual_fn (self.net, x)

  def train (self, x, u, w, xc, wc, steps):
    nout, inv = self.nout, 1 / self._t (self.scale)
    if x is not None:
      n, nq = x.shape [0], x.shape [1]
      xt, ut, wt = self._t (x.reshape (-1, self.nin)), self._t (u), self._t (w)
    if xc is not None:
      xct, wct = self._t (xc), self._t (wc)
    self.net.train ()
    for _ in range (steps):
      self.opt.zero_grad (set_to_none=True)
      data = phys = torch.zeros ((), dtype=self.tdtype, device=self.device)
      if x is not None:
        pred = self.net (xt).reshape (n, nq, nout).mean (dim=1)
        data = (wt * (((pred - ut) * inv) ** 2).mean (dim=1)).sum ()
      if xc is not None:
        r = self._residual (xct)
        phys = (wct * ((r * inv) ** 2).mean (dim=1)).sum ()
      loss = self.lambda_data * data + self.lambda_phys * phys
      loss.backward ()
      self.opt.step ()
    return loss.item (), data.item (), phys.item ()

  def forward_np (self, x):
    self.net.eval ()
    with torch.no_grad ():
      return self.net (self._t (x)).cpu ().numpy ()

  def residual_np (self, x):
    self.net.eval ()
    return self._residual (self._t (x)).detach ().cpu ().numpy ()

  def save (self, path):
    torch.save ({"net": self.net.state_dict (), "opt": self.opt.state_dict (),
                 "common": self.state_common ()}, path)

  def load (self, path):
    s = torch.load (path, map_location=self.device, weights_only=False)
    self.net.load_state_dict (s ["net"])
    self.opt.load_state_dict (s ["opt"])
    self.load_common (s ["common"])
