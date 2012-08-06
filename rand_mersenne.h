/**
 * @file rand_mersenne.h
 * @author Shawn Cokus (Cokus@math.washington.edu)
 * @date March 8, 1998
 * @brief Classe per generare numeri tramite Mersenne Prime Twister, @ref RandMT
 */
#ifndef _RANDMT_H_
#define _RANDMT_H_

#include <stdint.h>

///Mersenne Twister random number generator
class RandMT {
    /// length of state vector
    static const int N = 624;
    /// a period parameter
    static const int M = 397;
    /// a magic constant
    static const uint32_t K = 0x9908B0DFU; 

    // If you want a single generator, consider using a singleton class 
    // instead of trying to make these static.

    /// state vector + 1 extra to not violate ANSI C
    uint32_t state[N + 1];
    /// next random value is computed from here
    uint32_t *next;
    /// initial seed
    uint32_t initseed;
    /// can *next++ this many times before reloading - the number of random nums in the buffer
    int left; 

    /// mask all but highest   bit of u
    static inline uint32_t hiBit(uint32_t u) {
        return u & 0x80000000U; 
    }
    /// mask all but lowest    bit of u
    static inline uint32_t loBit(uint32_t u) {
        return u & 0x00000001U; 
    }
    /// mask     the highest   bit of u
    static inline uint32_t loBits(uint32_t u) {
        return u & 0x7FFFFFFFU; 
    }
    /// move hi bit of u to hi bit of v
    static inline uint32_t mixBits(uint32_t u, uint32_t v) {
        return hiBit(u) | loBits(v); 
    }
    ///Reload the generator from seed
    uint32_t reloadMT(void);

public:
    ///Constructor with auto-generated random seed by system interface
    RandMT();
    ///Costructor with fixed (repeatable) seed
    RandMT(uint32_t seed);

    ///return 32 bit unsigned integer (32 random bits)
    inline uint32_t rand_int32(void) {
        uint32_t y;

        if (--left < 0)
            return (reloadMT());

        y = *next++;
        y ^= (y >> 11);
        y ^= (y << 7) & 0x9D2C5680U;
        y ^= (y << 15) & 0xEFC60000U;
        return (y ^ (y >> 18));
    }

    ///an alias to rand_int32
    inline uint32_t get_int(void) {
        return rand_int32();
    }

    ///random 53 bit double in the interval [0,1)
    inline double get_double53() {
        return (static_cast<double> (rand_int32() >> 5) * 67108864. +
                static_cast<double> (rand_int32() >> 6)) * (1. / 9007199254740992.);
    }

    ///random double in the interval [0,1]
    inline double get_double() {
        return static_cast<double> (rand_int32()) * (1. / 4294967295.);
    } // divided by 2^32 - 1

    ///random double in the interval (0,1)
    inline double get_open_double() {
        return (static_cast<double> (rand_int32()) + .5) * (1. / 4294967296.);
    } // divided by 2^32

    ///seed by provided number
    void seedMT(uint32_t seed);
};

#endif // _RANDMT_H_
