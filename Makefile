LIBS     := $(patsubst %/,%,$(dir $(wildcard lib*/Makefile))) 
NAMES    := $(patsubst lib%,%,$(LIBS))
ARCHIVES := $(foreach l,$(LIBS),$(l)/$(l).a)

all: $(ARCHIVES)

# Phony so we always descend; the sub-make decides what is out of date.
$(ARCHIVES):
	$(MAKE) -C $(@D) $(@F)

# Inter-library dependencies, if any (libbar uses libfoo):
# libbar/libbar.a: libfoo/libfoo.a

test: all
	$(MAKE) -C tests LIBS="$(NAMES)"

clean:
	for d in $(LIBS) tests; do $(MAKE) -C $$d clean; done

.PHONY: all test clean $(ARCHIVES)

