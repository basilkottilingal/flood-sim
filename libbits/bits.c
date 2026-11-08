#include "bits.h"

/*
.. The de Bruijn tables: for p = 2^k, k is
..
..   bits_debruijn_table_32 [(p * 0x077CB531) >> 27]                (32 bits)
..   bits_debruijn_table_64 [(p * 0x03F79D71B4CB0A89) >> 58]        (64 bits)
..
.. Each sequence holds every 5 (6) bit pattern once in its windows, so the
.. top 5 (6) bits of p * B, which is B shifted left by k, differ for every
.. k. Generated from that definition, and checked for collisions.
*/
const unsigned char bits_debruijn_table_32 [32] =
{
   0,  1, 28,  2, 29, 14, 24,  3, 30, 22, 20, 15, 25, 17,  4,  8,
  31, 27, 13, 23, 21, 19, 16,  7, 26, 12, 18,  6, 11,  5, 10,  9
};

const unsigned char bits_debruijn_table_64 [64] =
{
   0,  1, 48,  2, 57, 49, 28,  3, 61, 58, 50, 42, 38, 29, 17,  4,
  62, 55, 59, 36, 53, 51, 43, 22, 45, 39, 33, 30, 24, 18, 12,  5,
  63, 47, 56, 27, 60, 41, 37, 16, 54, 35, 52, 21, 44, 32, 23, 11,
  46, 26, 40, 15, 34, 20, 31, 10, 25, 14, 19,  9, 13,  8,  7,  6
};
