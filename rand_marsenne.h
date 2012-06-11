/* 
 * File:   rand_marsenne.h
 * Author: XP
 *
 * Created on 9 czerwiec 2012, 22:56
 */

#ifndef _RANDMT_H_
#define _RANDMT_H_

typedef unsigned long uint32;

class RandMT {

  static const int N =          624;                // length of state vector
  static const int M =          397;                // a period parameter
  static const uint32 K =       0x9908B0DFU;        // a magic constant

  // If you want a single generator, consider using a singleton class 
  // instead of trying to make these static.
  uint32   state[N+1];  // state vector + 1 extra to not violate ANSI C
  uint32   *next;       // next random value is computed from here
  uint32   initseed;    //
  int      left;        // can *next++ this many times before reloading

  inline uint32 hiBit(uint32 u) { 
    return u & 0x80000000U;    // mask all but highest   bit of u
  }

  inline uint32 loBit(uint32 u) { 
    return u & 0x00000001U;    // mask all but lowest    bit of u
  }

  inline uint32 loBits(uint32 u) { 
    return u & 0x7FFFFFFFU;   // mask     the highest   bit of u
  }

  inline uint32 mixBits(uint32 u, uint32 v) {
    return hiBit(u)|loBits(v);  // move hi bit of u to hi bit of v
  }
  
  uint32 reloadMT(void) ;

public:
  RandMT() ;
  RandMT(uint32 seed) ;
  uint32 rand_int32();
  inline uint32 get_int(void) {return rand_int32() ;}
   //53 bit double in the interval [0,1)
 inline double get_double53() {
    return (static_cast<double>(rand_int32() >> 5) * 67108864. + 
      static_cast<double>(rand_int32() >> 6)) * (1. / 9007199254740992.); }
  //return in the interval [0,1]
 inline double get_float() {
    return static_cast<double>(rand_int32()) * (1. / 4294967295.); } // divided by 2^32 - 1
  //return in the interval (0,1)
 inline double get_open_float() {
    return (static_cast<double>(rand_int32()) + .5) * (1. / 4294967296.); } // divided by 2^32
  void seedMT(uint32 seed) ;
};

#endif // _RANDMT_H_
