BUILD  := build
LEXER  ?= lexer
CFLAGS ?= -g -O1
CFLAGS += -std=gnu11 -Wall -Wextra -D_GNU_SOURCE -Isrc -I$(BUILD)

OBJS := $(BUILD)/driver.o $(BUILD)/translate.o $(BUILD)/ast.o \
        $(BUILD)/allocator.o $(BUILD)/scope.o $(BUILD)/diag.o \
        $(BUILD)/parser.tab.o $(BUILD)/lexer.yy.o

all: $(BUILD)/cwrap $(BUILD)/scalar.o

$(BUILD)/cwrap: $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^

# The runtime that cwrap adds to what it compiles: see src/scalar.h.
$(BUILD)/driver.o: CFLAGS += -DCWRAP_INCLUDE_DIR='"$(abspath src)"' \
                             -DCWRAP_LIB_DIR='"$(abspath $(BUILD))"' \
                             -DCWRAP_HEADER_DIR='"$(abspath include)"'

$(BUILD)/scalar.o: src/scalar.c src/scalar.h | $(BUILD)
	$(CC) $(CFLAGS) -c -o $@ $<

$(BUILD)/parser.tab.c $(BUILD)/parser.tab.h: src/parser.y | $(BUILD)
	bison -Wall -d -o $(BUILD)/parser.tab.c $<

$(BUILD)/lexer.yy.c: src/lexer.l | $(BUILD)
	flex -o $@ $<

# Table-based lexer, -d for debug. Not linked into cwrap until it implements
# lexer.h, which needs a reentrant mode: it keeps its state in static variables.
$(BUILD)/lexer.c: src/lexer.l | $(BUILD)
	$(LEXER) -d -o $@ $<

lexer: $(BUILD)/lexer.c

# Placeholder lexer: generated code is not held to our warning flags.
$(BUILD)/lexer.yy.o: $(BUILD)/lexer.yy.c $(BUILD)/parser.tab.h src/ast.h src/lexer.h src/scope.h src/allocator.h
	$(CC) $(CFLAGS) -Wno-unused-function -Wno-sign-compare -c -o $@ $<

$(BUILD)/parser.tab.o: $(BUILD)/parser.tab.c src/ast.h src/lexer.h src/scope.h src/allocator.h src/translate.h src/diag.h
	$(CC) $(CFLAGS) -c -o $@ $<

$(BUILD)/%.o: src/%.c src/ast.h src/lexer.h src/scope.h src/allocator.h src/translate.h src/diag.h $(BUILD)/parser.tab.h | $(BUILD)
	$(CC) $(CFLAGS) -c -o $@ $<

$(BUILD):
	mkdir -p $@

# tests/scope_test.c includes scope.c, to see its free lists.
$(BUILD)/scope_test: tests/scope_test.c src/scope.c src/scope.h $(BUILD)/allocator.o
	$(CC) $(CFLAGS) -o $@ tests/scope_test.c $(BUILD)/allocator.o

test: $(BUILD)/cwrap $(BUILD)/scalar.o $(BUILD)/scope_test
	sh tests/run.sh

clean:
	rm -rf $(BUILD)

.PHONY: all test clean lexer
