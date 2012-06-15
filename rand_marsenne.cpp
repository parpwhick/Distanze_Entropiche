// This is the ``Mersenne Twister'' random number generator MT19937, which
// generates pseudorandom integers uniformly distributed in 0..(2^32 - 1)
// starting from any odd seed in 0..(2^32 - 1).  This version is a recode
// by Shawn Cokus (Cokus@math.washington.edu) on March 8, 1998 of a version by
// Takuji Nishimura (who had suggestions from Topher Cooper and Marc Rieffel in
// July-August 1997).
//
// According to the URL <http://www.math.keio.ac.jp/~matumoto/emt.html>
// (and paraphrasing a bit in places), the Mersenne Twister is ``designed
// with consideration of the flaws of various existing generators,'' has
// a period of 2^19937 - 1, gives a sequence that is 623-dimensionally
// equidistributed, and ``has passed many stringent tests, including the
// die-hard test of G. Marsaglia and the load test of P. Hellekalek and
// S. Wegenkittl.''  It is efficient in memory usage (typically using 2506
// to 5012 bytes of static data, depending on data type sizes, and the code
// is quite short as well).  It generates random numbers in batches of 624
// at a time, so the caching and pipelining of modern systems is exploited.
// It is also divide- and mod-free.
//
// This library is free software; you can redistribute it and/or modify it
// under the terms of the GNU Library General Public License as published by
// the Free Software Foundation (either version 2 of the License or, at your
// option, any later version).  This library is distributed in the hope that
// it will be useful, but WITHOUT ANY WARRANTY, without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See
// the GNU Library General Public License for more details.  You should have
// received a copy of the GNU Library General Public License along with this
// library; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307, USA.

#include "rand_marsenne.h"
#include "time.h"

RandMT::RandMT() {
    uint32 seed;
#ifdef __linux__
    FILE *out = fopen("/dev/urandom", "r");
    int bytes_read = fread(&seed, sizeof (long), 1, out);
    if (!bytes_read) {
        printf("Failed reading the seed, which wasn't provided\n");
        exit(1);
    }
    fclose(out);
#else
    seed = time(0);
#endif
    seedMT(seed);
}

RandMT::RandMT(uint32 seed) {
    seedMT(seed);
}

void RandMT::seedMT(uint32 seed) {
    //
    // We initialize state[0..(N-1)] via the generator
    //
    //   x_new = (69069 * x_old) mod 2^32
    //
    // from Line 15 of Table 1, p. 106, Sec. 3.3.4 of Knuth's
    // _The Art of Computer Programming_, Volume 2, 3rd ed.
    //
    // Notes (SJC): I do not know what the initial state requirements
    // of the Mersenne Twister are, but it seems this seeding generator
    // could be better.  It achieves the maximum period for its modulus
    // (2^30) iff x_initial is odd (p. 20-21, Sec. 3.2.1.2, Knuth); if
    // x_initial can be even, you have sequences like 0, 0, 0, ...;
    // 2^31, 2^31, 2^31, ...; 2^30, 2^30, 2^30, ...; 2^29, 2^29 + 2^31,
    // 2^29, 2^29 + 2^31, ..., etc. so I force seed to be odd below.
    //
    // Even if x_initial is odd, if x_initial is 1 mod 4 then
    //
    //   the          lowest bit of x is always 1,
    //   the  next-to-lowest bit of x is always 0,
    //   the 2nd-from-lowest bit of x alternates      ... 0 1 0 1 0 1 0 1 ... ,
    //   the 3rd-from-lowest bit of x 4-cycles        ... 0 1 1 0 0 1 1 0 ... ,
    //   the 4th-from-lowest bit of x has the 8-cycle ... 0 0 0 1 1 1 1 0 ... ,
    //    ...
    //
    // and if x_initial is 3 mod 4 then
    //
    //   the          lowest bit of x is always 1,
    //   the  next-to-lowest bit of x is always 1,
    //   the 2nd-from-lowest bit of x alternates      ... 0 1 0 1 0 1 0 1 ... ,
    //   the 3rd-from-lowest bit of x 4-cycles        ... 0 0 1 1 0 0 1 1 ... ,
    //   the 4th-from-lowest bit of x has the 8-cycle ... 0 0 1 1 1 1 0 0 ... ,
    //    ...
    //
    // The generator's potency (min. s>=0 with (69069-1)^s = 0 mod 2^32) is
    // 16, which seems to be alright by p. 25, Sec. 3.2.1.3 of Knuth.  It
    // also does well in the dimension 2..5 spectral tests, but it could be
    // better in dimension 6 (Line 15, Table 1, p. 106, Sec. 3.3.4, Knuth).
    //
    // Note that the random number user does not see the values generated
    // here directly since reloadMT() will always munge them first, so maybe
    // none of all of this matters.  In fact, the seed values made here could
    // even be extra-special desirable if the Mersenne Twister theory says
    // so-- that's why the only change I made is to restrict to odd seeds.
    //
    initseed = seed;

    register uint32 x = (seed | 1U) & 0xFFFFFFFFU, *s = state;
    register int j;
    left = 0;
    for (*s++ = x, j = N; --j; *s++ = (x *= 69069U) & 0xFFFFFFFFU);
}

uint32 RandMT::reloadMT(void) {
    register uint32 *p0 = state, *p2 = state + 2, *pM = state + M, s0, s1;
    register int j;

    if (left < -1)
        seedMT(initseed);

    left = N - 1, next = state + 1;

    for (s0 = state[0], s1 = state[1], j = N - M + 1; --j; s0 = s1, s1 = *p2++)
        *p0++ = *pM++ ^ (mixBits(s0, s1) >> 1) ^ (loBit(s1) ? K : 0U);

    for (pM = state, j = M; --j; s0 = s1, s1 = *p2++)
        *p0++ = *pM++ ^ (mixBits(s0, s1) >> 1) ^ (loBit(s1) ? K : 0U);

    s1 = state[0], *p0 = *pM ^ (mixBits(s0, s1) >> 1) ^ (loBit(s1) ? K : 0U);
    s1 ^= (s1 >> 11);
    s1 ^= (s1 << 7) & 0x9D2C5680U;
    s1 ^= (s1 << 15) & 0xEFC60000U;
    return (s1 ^ (s1 >> 18));
}

uint32 RandMT::rand_int32(void) {
    uint32 y;

    if (--left < 0)
        return (reloadMT());

    y = *next++;
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9D2C5680U;
    y ^= (y << 15) & 0xEFC60000U;
    return (y ^ (y >> 18));
}

#ifdef RANDOM_TEST
#include <cstdio>
// A simple test

int main(void) {
    int j;
    unsigned long l = 0;
    double l2 = 0;

    RandMT r(4357U);

    // you can seed with any uint32, but the best are odds in 0..(2^32 - 1)

    // Run this 400 million times
    for (j = 400000000; j > 0; --j) {
        l += r.get_int();
    }
    printf("control value: %lu\n", l);
    //printf("control value: %g\n",l2);
    return (0);
}

#endif
