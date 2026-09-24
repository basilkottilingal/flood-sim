CC       = gcc
CCFLAGS += -D_XOPEN_SOURCE=700 -O2

AR       = ar
ARFLAGS  = cr

SRC      =
OBJ      = $(SRC:.c=.o)

%.o: %.c %.h Makefile
	$(CC) $< $(CCFLAGS) -c -o $@ -lm

libflow-network.a : $(OBJ)
	$(AR) $(ARFLAGS) $@ $^

%.tst       : %.c libflow-network.a
	$(CC) $(CCFLAGS) -o $* $< -L./ -I./ -lflow-network -lm

all:

clean        :
	rm *.o
	rm libflow-network.a

