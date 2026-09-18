CC       = gcc
CCFLAGS += -D_XOPEN_SOURCE=700 -O2

AR       = ar
ARFLAGS  = cr

SRC      = filemap.c parse.c tiff.c lzw.c coordinates.c geotiff.c
OBJ      = $(SRC:.c=.o)

%.o: %.c %.h Makefile
	$(CC) $< $(CCFLAGS) -c -o $@ -lm

libgeotiff.a : $(OBJ)
	$(AR) $(ARFLAGS) $@ $^

test         : test.c libgeotiff.a
	$(CC) $(CCFLAGS) -o $@ $< -L./ -I./ -lgeotiff

clean        :
	rm *.o
	rm libgeotiff.a

