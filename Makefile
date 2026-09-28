# dependents first, providers last (link order)
LINK_ORDER := flow-network geotiff view alloc
LIBS     := $(addprefix lib,$(LINK_ORDER))

#LIBS     := $(patsubst %/,%,$(dir $(wildcard lib*/Makefile))) 

NAMES    := $(patsubst lib%,%,$(LIBS))
ARCHIVES := $(foreach l,$(LIBS),$(l)/$(l).a)

all: $(ARCHIVES)

# Phony so we always descend; the sub-make decides what is out of date.
$(ARCHIVES):
	$(MAKE) -C $(@D) $(@F)

# Inter-library dependencies, if any
libview/libview.a: libgeotiff/libgeotiff.a

test: all
	$(MAKE) -C tests LIBS="$(NAMES)"

clean:
	for d in $(LIBS) tests; do $(MAKE) -C $$d clean; done

.PHONY: all test clean $(ARCHIVES)

