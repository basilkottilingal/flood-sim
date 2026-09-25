CC        = gcc
LIBS     ?= geotiff      # names without the "lib" prefix; top level overrides this
CPPFLAGS += $(addprefix -I../lib,$(LIBS))
LDFLAGS  += $(addprefix -L../lib,$(LIBS))
LDLIBS   += $(addprefix -l,$(LIBS))
CFLAGS   ?= -Wall -Wextra -g

TESTS    := $(wildcard *.c)
SERIAL   := $(TESTS:.c=.s)
ARCHIVES := $(foreach l,$(LIBS),../lib$(l)/lib$(l).a)

#all: $(SERIAL)

# Relink tests when a library changes
$(SERIAL): $(ARCHIVES)

%.s       : %.c Makefile
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -o $@ $< $(LDLIBS) -lm

all:
	make read-geotiff.s && ./read-geotiff.s ../topography_ESA_Copernicus_30m_resolution.tif
	make geotiff-crs-paris.s && ./geotiff-crs-paris.s

clean:
	$(RM) *.s

.PHONY: all clean
