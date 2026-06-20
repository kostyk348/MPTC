# Makefile for MPTC library — used for PC build/test without CMake.
# This is a fallback when cmake is not available. The CMake build is the
# canonical one for distribution.

CC      ?= cc
CFLAGS  ?= -std=c99 -O2 -Wall -Wextra -Wno-unused-parameter -Wno-unused-function
CFLAGS  += -D_POSIX_C_SOURCE=199309L -DMPTC_TEST_HOST=1 -DMPTC_USE_FLOAT=1 -Iinclude
LDFLAGS ?= -lm

SRC     = src/mptc.c src/sdft.c src/integral_img.c
OBJ     = $(SRC:.c=.o)

LIB     = libmptc.a

.PHONY: all clean test bench

all: $(LIB)

$(LIB): $(OBJ)
	ar rcs $@ $(OBJ)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

test: tests/test_basic.c $(LIB)
	$(CC) $(CFLAGS) tests/test_basic.c $(LIB) $(LDFLAGS) -o test_basic
	./test_basic

bench: examples/pc/pc_benchmark.c $(LIB)
	$(CC) $(CFLAGS) examples/pc/pc_benchmark.c $(LIB) $(LDFLAGS) -o pc_benchmark
	./pc_benchmark

clean:
	rm -f $(OBJ) $(LIB) test_basic pc_benchmark
