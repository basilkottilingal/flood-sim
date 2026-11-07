/*
.. libpinn core: the embedded interpreter, and the calls to the backend
.. object. The only file of libpinn that includes Python.h: see pinn.h.
..
.. Arrays cross without a copy: each becomes a memoryview of the caller's
.. memory, that pinn_backend.py reads with numpy.frombuffer (), and writes
.. into for predictions. The views live only during the call.
*/
#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pinn.h"

#ifndef PINN_PYDIR
  #define PINN_PYDIR "python"
#endif

#define PINN_ERRLEN 4096

struct Pinn
{
  PyObject * obj;           /* the backend object: pinn_backend.Model */
  int nin, nout;
  int has_residual;
  char err [PINN_ERRLEN];
};

static char err_global [PINN_ERRLEN];  /* pinn_new () failures         */
static int py_owned = 0;               /* we started the interpreter    */
static int py_alive = 0;               /* started, and not shut down    */
static int nlive = 0;                  /* Pinn not freed                */
static PyThreadState * py_main = NULL; /* the main thread, GIL released */
static PyObject * backend_mod = NULL;  /* pinn_backend                  */

static void set_err (char * dst, const char * msg)
{
  snprintf (dst, PINN_ERRLEN, "%s", msg);
}

/*
.. the pending Python exception into dst, traceback included, and clears
.. it. With the GIL held.
*/
static void fetch_err (char * dst, const char * what)
{
  PyObject * type, * value, * tb;
  PyErr_Fetch (&type, &value, &tb);
  if (!type)
  {
    snprintf (dst, PINN_ERRLEN, "%s: failed, no Python exception", what);
    return;
  }
  PyErr_NormalizeException (&type, &value, &tb);
  PyObject * text = NULL, * mod = PyImport_ImportModule ("traceback");
  if (mod)
  {
    PyObject * lines = PyObject_CallMethod (mod, "format_exception", "OOO",
      type, value ? value : Py_None, tb ? tb : Py_None);
    if (lines)
    {
      PyObject * sep = PyUnicode_FromString ("");
      text = sep ? PyUnicode_Join (sep, lines) : NULL;
      Py_XDECREF (sep);
      Py_DECREF (lines);
    }
    Py_DECREF (mod);
  }
  if (!text && value)
    text = PyObject_Str (value);
  const char * s = text ? PyUnicode_AsUTF8 (text) : NULL;
  snprintf (dst, PINN_ERRLEN, "%s: %s", what, s ? s : "(no message)");
  Py_XDECREF (text);
  Py_XDECREF (type);
  Py_XDECREF (value);
  Py_XDECREF (tb);
  PyErr_Clear ();
}

/*
.. starts the interpreter if no one has, puts the directory of
.. pinn_backend.py first in sys.path (PINN_PYDIR, or the environment
.. variable of the same name), imports it, and leaves the GIL released.
*/
static int python_up (void)
{
  if (backend_mod)
    return 0;
  if (!Py_IsInitialized ())
  {
    if (py_alive)
    {
      set_err (err_global, "pinn: the interpreter was shut down");
      return -1;
    }
    Py_InitializeEx (0);          /* 0: leave the program's signals alone */
    py_owned = 1;
    py_main = PyEval_SaveThread ();
  }
  py_alive = 1;
  PyGILState_STATE g = PyGILState_Ensure ();
  int ret = 0;
  const char * dir = getenv ("PINN_PYDIR");
  PyObject * path = PySys_GetObject ("path");          /* borrowed */
  PyObject * d = PyUnicode_FromString (dir && *dir ? dir : PINN_PYDIR);
  if (!path || !d || PyList_Insert (path, 0, d) < 0)
  {
    fetch_err (err_global, "pinn: sys.path");
    ret = -1;
  }
  Py_XDECREF (d);
  if (!ret && !(backend_mod = PyImport_ImportModule ("pinn_backend")))
  {
    fetch_err (err_global, "pinn: import pinn_backend");
    ret = -1;
  }
  PyGILState_Release (g);
  return ret;
}

/* a memoryview of n items of size bytes at a, or None for NULL */
static PyObject * view (const void * a, size_t n, size_t size, int writable)
{
  if (!a)
    Py_RETURN_NONE;
  return PyMemoryView_FromMemory ((char *) a, (Py_ssize_t) (n * size),
    writable ? PyBUF_WRITE : PyBUF_READ);
}

static int dict_set (PyObject * d, const char * key, PyObject * v)
{
  if (!v)
    return -1;
  int r = PyDict_SetItemString (d, key, v);
  Py_DECREF (v);
  return r;
}

static PyObject * str_or_none (const char * s)
{
  if (!s)
    Py_RETURN_NONE;
  return PyUnicode_FromString (s);
}

static PyObject * box (const double * a, int n)
{
  PyObject * l = PyList_New (n);
  for (int i = 0; l && i < n; i++)
    PyList_SET_ITEM (l, i, PyFloat_FromDouble (a [i]));
  return l;
}

static PyObject * config_dict (const PinnConfig * c)
{
  PyObject * d = PyDict_New ();
  if (!d)
    return NULL;
  int e = 0;
  e |= dict_set (d, "backend", str_or_none (c->backend));
  e |= dict_set (d, "nin", PyLong_FromLong (c->nin));
  e |= dict_set (d, "nout", PyLong_FromLong (c->nout));
  e |= dict_set (d, "width", PyLong_FromLong (c->width));
  e |= dict_set (d, "depth", PyLong_FromLong (c->depth));
  e |= dict_set (d, "activation", str_or_none (c->activation));
  e |= dict_set (d, "fourier", PyLong_FromLong (c->fourier));
  e |= dict_set (d, "fourier_sigma", PyFloat_FromDouble (c->fourier_sigma));
  e |= dict_set (d, "lr", PyFloat_FromDouble (c->lr));
  e |= dict_set (d, "lambda_data", PyFloat_FromDouble (c->lambda_data));
  e |= dict_set (d, "lambda_phys", PyFloat_FromDouble (c->lambda_phys));
  e |= dict_set (d, "lo", box (c->lo, c->nin));
  e |= dict_set (d, "hi", box (c->hi, c->nin));
  e |= dict_set (d, "residual", str_or_none (c->residual));
  e |= dict_set (d, "device", str_or_none (c->device));
  e |= dict_set (d, "dtype", str_or_none (c->dtype));
  e |= dict_set (d, "seed", PyLong_FromLong (c->seed));
  if (e)
  {
    Py_DECREF (d);
    return NULL;
  }
  return d;
}

PinnConfig pinn_config (const char * backend, int nin, int nout)
{
  PinnConfig c;
  memset (&c, 0, sizeof (c));
  c.backend = backend;
  c.nin = nin;
  c.nout = nout;
  c.width = 64;
  c.depth = 4;
  c.activation = "tanh";
  c.fourier_sigma = 1.;
  c.lr = 1e-3;
  c.lambda_data = 1.;
  c.lambda_phys = 1.;
  for (int i = 0; i < PINN_MAX_IN; i++)
    c.lo [i] = -1., c.hi [i] = 1.;
  c.device = "cpu";
  c.dtype = "float32";
  return c;
}

Pinn * pinn_new (const PinnConfig * cfg)
{
  if (!cfg || !cfg->backend || cfg->nin < 1 || cfg->nin > PINN_MAX_IN ||
      cfg->nout < 1 || cfg->nout > PINN_MAX_OUT)
  {
    set_err (err_global, "pinn_new: bad config (backend, nin or nout)");
    return NULL;
  }
  for (int i = 0; i < cfg->nin; i++)
    if (!(cfg->hi [i] > cfg->lo [i]))
    {
      set_err (err_global, "pinn_new: empty box, hi <= lo");
      return NULL;
    }
  if (python_up ())
    return NULL;
  Pinn * p = calloc (1, sizeof (Pinn));
  if (!p)
  {
    set_err (err_global, "pinn_new: out of memory");
    return NULL;
  }
  p->nin = cfg->nin;
  p->nout = cfg->nout;
  p->has_residual = cfg->residual != NULL;
  PyGILState_STATE g = PyGILState_Ensure ();
  PyObject * d = config_dict (cfg);
  if (d)
  {
    p->obj = PyObject_CallMethod (backend_mod, "create", "O", d);
    Py_DECREF (d);
  }
  if (!p->obj)
    fetch_err (err_global, "pinn_new");
  PyGILState_Release (g);
  if (!p->obj)
  {
    free (p);
    return NULL;
  }
  nlive++;
  return p;
}

void pinn_free (Pinn * p)
{
  if (!p)
    return;
  if (p->obj && Py_IsInitialized ())
  {
    PyGILState_STATE g = PyGILState_Ensure ();
    Py_DECREF (p->obj);
    PyGILState_Release (g);
  }
  free (p);
  nlive--;
}

/* calls the method, keeps a failure in p->err; returns its result or NULL */
static PyObject * call (Pinn * p, const char * method, PyObject * args)
{
  PyObject * f = args ? PyObject_GetAttrString (p->obj, method) : NULL;
  PyObject * r = f ? PyObject_CallObject (f, args) : NULL;
  if (!r)
  {
    char what [64];
    snprintf (what, sizeof (what), "pinn %s", method);
    fetch_err (p->err, what);
  }
  Py_XDECREF (f);
  Py_XDECREF (args);
  return r;
}

int pinn_fit (Pinn * p, const PinnBatch * b, int steps, PinnLoss * loss)
{
  if (!p || !b || b->n < 0 || b->nc < 0 || (b->n > 0 && (!b->x || !b->u)) ||
      (b->nc > 0 && !b->xc) || steps < 0)
  {
    if (p)
      set_err (p->err, "pinn_fit: bad batch");
    return -1;
  }
  int nq = b->nq > 0 ? b->nq : 1;
  PyGILState_STATE g = PyGILState_Ensure ();
  size_t n = (size_t) b->n, nc = (size_t) b->nc;
  PyObject * args = Py_BuildValue ("(NNNNNiiii)",
    view (b->n ? b->x : NULL, n * nq * p->nin, sizeof (double), 0),
    view (b->n ? b->u : NULL, n * p->nout, sizeof (float), 0),
    view (b->n ? b->w : NULL, n, sizeof (float), 0),
    view (b->nc ? b->xc : NULL, nc * p->nin, sizeof (double), 0),
    view (b->nc ? b->wc : NULL, nc, sizeof (float), 0),
    b->n, nq, b->nc, steps);
  PyObject * r = call (p, "fit", args);
  int ret = -1;
  if (r)
  {
    double t = 0, d = 0, ph = 0;
    if (PyArg_ParseTuple (r, "ddd", &t, &d, &ph))
    {
      if (loss)
        loss->total = t, loss->data = d, loss->phys = ph;
      ret = 0;
    }
    else
      fetch_err (p->err, "pinn_fit: fit () must return 3 floats");
    Py_DECREF (r);
  }
  PyGILState_Release (g);
  return ret;
}

/* method (x, out, n): out [n * nout] from x [n * nin] */
static int eval (Pinn * p, const char * method, int n, const double * x,
                 float * out)
{
  if (!p || n < 0 || (n > 0 && (!x || !out)))
  {
    if (p)
      set_err (p->err, "pinn: bad arguments");
    return -1;
  }
  if (n == 0)
    return 0;
  PyGILState_STATE g = PyGILState_Ensure ();
  PyObject * args = Py_BuildValue ("(NNi)",
    view (x, (size_t) n * p->nin, sizeof (double), 0),
    view (out, (size_t) n * p->nout, sizeof (float), 1), n);
  PyObject * r = call (p, method, args);
  Py_XDECREF (r);
  PyGILState_Release (g);
  return r ? 0 : -1;
}

int pinn_predict (Pinn * p, int n, const double * x, float * u)
{
  return eval (p, "predict", n, x, u);
}

int pinn_residual (Pinn * p, int n, const double * x, float * r)
{
  if (p && !p->has_residual)
  {
    set_err (p->err, "pinn_residual: the config names no residual");
    return -1;
  }
  return eval (p, "residual", n, x, r);
}

static int file_call (Pinn * p, const char * method, const char * path)
{
  if (!p || !path)
    return -1;
  PyGILState_STATE g = PyGILState_Ensure ();
  PyObject * r = call (p, method, Py_BuildValue ("(s)", path));
  Py_XDECREF (r);
  PyGILState_Release (g);
  return r ? 0 : -1;
}

int pinn_save (Pinn * p, const char * path)
{
  return file_call (p, "save", path);
}

int pinn_load (Pinn * p, const char * path)
{
  return file_call (p, "load", path);
}

const char * pinn_error (const Pinn * p)
{
  return p ? p->err : err_global;
}

void pinn_shutdown (void)
{
  if (!py_alive || nlive > 0)
    return;
  if (backend_mod)
  {
    PyGILState_STATE g = PyGILState_Ensure ();
    Py_CLEAR (backend_mod);
    PyGILState_Release (g);
  }
  if (py_owned)
  {
    PyEval_RestoreThread (py_main);
    Py_FinalizeEx ();
    py_owned = 0;
    py_main = NULL;
  }
}
