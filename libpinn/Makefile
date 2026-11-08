# libpinn: two layers, two compilers (see README.md)
#   pinn.o       the Python bridge, by $(CC): the only file with Python.h
#   pinn-bamr.o  the libb-amr adapter, by cwrap -dimension=$(DIMENSION)
# Programs link libpinn.a, libb-amr, and $(PINN_LDLIBS):
#   make print-ldlibs
CC         ?= cc
CWRAP      ?= ../libc-pp/build/cwrap
DIMENSION  ?= 2
PYTHON     ?= python3
PYDIR      := $(abspath python)

PY_CFLAGS  := $(shell $(PYTHON)-config --includes)
PINN_LDLIBS := $(shell $(PYTHON)-config --ldflags --embed) -lm

CFLAGS     += -Wall -Wextra -O2 -D_GNU_SOURCE
LIB        := libpinn.a

$(LIB): pinn.o pinn-bamr.o
	$(AR) rcs $@ $^

pinn.o: pinn.c pinn.h Makefile
	$(CC) $(CFLAGS) $(PY_CFLAGS) -DPINN_PYDIR='"$(PYDIR)"' -c $< -o $@

pinn-bamr.o: pinn-bamr.c pinn-bamr.h pinn.h Makefile
	$(CWRAP) -dimension=$(DIMENSION) $(CFLAGS) -I. -I.. -I../libb-amr -c $< -o $@

print-ldlibs:
	@echo $(PINN_LDLIBS)

clean:
	$(RM) pinn.o pinn-bamr.o $(LIB)

.PHONY: clean print-ldlibs
