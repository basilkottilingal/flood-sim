NAME     := $(notdir $(CURDIR))
LIB      := $(NAME).a
SRCS     := $(wildcard *.c)
OBJS     := $(SRCS:.c=.o)
CFLAGS   ?= -Wall -Wextra -O2 -D_XOPEN_SOURCE=700
CPPFLAGS += -I.

$(LIB): $(OBJS)
	$(AR) rcs $@ $^

%.o: %.c %.h Makefile
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

clean:
	$(RM) $(OBJS) $(LIB)

.PHONY: clean
