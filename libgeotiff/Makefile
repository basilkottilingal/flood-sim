CC       = gcc
CCFLAGS += -D_XOPEN_SOURCE=700 -O2

AR       = ar
ARFLAGS  = cr

SRC      = filemap.c parse.c tiff.c lzw.c coordinates.c geodetic.c geotiff.c
OBJ      = $(SRC:.c=.o)

%.o: %.c %.h Makefile
	$(CC) $< $(CCFLAGS) -c -o $@ -lm

libgeotiff.a : $(OBJ)
	$(AR) $(ARFLAGS) $@ $^

%.tst       : %.c libgeotiff.a
	$(CC) $(CCFLAGS) -o $* $< -L./ -I./ -lgeotiff -lm

all:
	make test.tst && ./test ../topography_ESA_Copernicus_30m_resolution.tif
	make test-crs.tst && ./test-crs

clean        :
	rm *.o
	rm libgeotiff.a

