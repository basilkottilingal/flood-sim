#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>

#include "quadtree.h"

#define G(x) grid_##x

static int G(nb) = 1;         /* number of blocks                 */
static int G(nlb) = 1;        /* number of leaf blocks            */
static double G(x0) = 0.;     /* left corner of domain            */
static double G(y0) = 0.;     /* bottom corner of domain          */
static double G(h) = 1.;      /* domain [x0, x0+h] x [y0, y0+h]   */
static int G(depth) = 1;      /* depth of the current bAMR        */
static int G(depthmin) = 1;   /* min {level of leaves}            */
static int G(depthmax) = 1;   /* max depth forced by user         */
static int G(blockdepth) = 1; /* number of cells in block 2^(2*blockdepth) */

static int grid_init (int depth, int depthmax, int blockdepth)
{
  G(nb) = 
}

void quadtree (int depth, int depthmax, int blockdepth)
{
  depthmax = depthmax > QAMR_MAX_DEPTH ? QAMR_MAX_DEPTH :
    depthmax <= 0 ? 1 : depthmax;
  depth = depth > depthmax ? depthmax : depth < 0 : 0;
  blockdepth = blockdepth < 4 ? 4 : blockdepth > 10 ? 10 : blockdepth;
  grid_init (depth, depthmax, blockdepth);
}

void quadtree_offset (double x0, double y0)
{
  G(x0) = x0, G(y0) = y0;
}
