NAME     := $(notdir $(CURDIR))
LIB      := $(NAME).a
SRCS     := $(wildcard *.c)
OBJS     := $(SRCS:.c=.o)
CFLAGS   ?= -Wall -Wextra -O2 -D_XOPEN_SOURCE=700

DEPS     ?=
CPPFLAGS += -I. $(foreach d,$(DEPS),-I../lib$(d))
LDFLAGS  += $(foreach d,$(DEPS),-L../lib$(d))
LDLIBS   += $(foreach d,$(DEPS),-l$(d))
DEPLIBS  := $(foreach d,$(DEPS),../lib$(d)/lib$(d).a)


$(LIB): $(OBJS)
	$(AR) rcs $@ $^

%.o: %.c %.h Makefile
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

clean:
	$(RM) $(OBJS) $(LIB)

.PHONY: clean
