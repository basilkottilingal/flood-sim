# liblinear-algebra: solvers on libb-amr's tree, built by cwrap
#   make                  liblinear-algebra.a, for DIMENSION (default 2)
#   sh tests/run.sh       the tests, 1-3D, double and float
CWRAP      ?= ../libc-pp/build/cwrap
DIMENSION  ?= 2
CFLAGS     += -Wall -Wextra -O2 -D_GNU_SOURCE
LIB        := liblinear-algebra.a

$(LIB): multigrid.o
	$(AR) rcs $@ $^

multigrid.o: multigrid.c multigrid.h Makefile
	$(CWRAP) -dimension=$(DIMENSION) $(CFLAGS) -I. -I.. -I../libb-amr -c $< -o $@

clean:
	$(RM) multigrid.o $(LIB)

.PHONY: clean
