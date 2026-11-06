#ifndef _QUADTREE_MORTON_CURVE_H_
#define _QUADTREE_MORTON_CURVE_H_

  #include <assert.h>
  #include "tree.h"
  #include "liballoc/array.h"

  /*
  .. morton_start (IJKL) { ... } morton_end ();
  .. runs the body once for every node of the subtree at IJKL (TREE_ROOT for
  .. the whole tree), in any dimension. A leaf is visited when it is reached,
  .. a parent after all its children. In the body, IJKL, node and sblocks are
  .. those of the node visited. 'break' leaves the body, not the traversal.
  ..
  .. Each step down pushes the parent and all but the first child, so the
  .. stack never holds more than TREE_NCHILDREN per level.
  */
  #define morton_start(...) do {                                   \
    Index _cur_ = {__VA_ARGS__};                                   \
    struct _mstack_ { Index p; int visited; }                      \
      _stack_ [TREE_NCHILDREN * BAMR_MAX_DEPTH + 1];               \
    int _depth_ = 0, _visited_ = 0;                                \
    if (!index_is_allocated (_cur_))                               \
      break;                                                       \
    for (;;)                                                       \
    {                                                              \
      while (!_visited_)                                           \
      {                                                            \
        BAMRNode * _n_ = index_node (_cur_);                       \
        if (!tree_has_children (_n_))                              \
          break;                                                   \
        assert (_cur_.level < BAMR_MAX_DEPTH);                     \
        _stack_ [_depth_++] = (struct _mstack_) {_cur_, 1};        \
        for (int _c_ = TREE_NCHILDREN - 1; _c_ > 0; _c_--)         \
          _stack_ [_depth_++] =                                    \
            (struct _mstack_) {index_child (_cur_, _c_), 0};       \
        _cur_ = index_child (_cur_, 0);                            \
      }                                                            \
      do {                                                         \
        INDEX_UNPACK (_cur_);                                      \
        BAMR_NODE_ENV ();                                          \
        /* do something here */
  #define morton_end()                                             \
      } while (0);                                                 \
      if (!_depth_)                                                \
        break;                                                     \
      _depth_--;                                                   \
      _cur_     = _stack_ [_depth_].p;                             \
      _visited_ = _stack_ [_depth_].visited;                       \
    }                                                              \
  } while (0)

  /* API */
  Array * tree_block_indices (void);

#endif
