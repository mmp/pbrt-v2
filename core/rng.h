
/*
    pbrt source code Copyright(c) 1998-2009 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    pbrt is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.  Note that the text contents of
    the book "Physically Based Rendering" are *not* licensed under the
    GNU GPL.

    pbrt is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */

#ifndef PBRT_CORE_RNG_H
#define PBRT_CORE_RNG_H

// core/rng.h*
#include "pbrt.h"

// Random Number Declarations
#define MTRNG RNG
//#define TURNG RNG

class MTRNG {
public:
    MTRNG(u_long seed = 5489UL) {
        mti = N+1; /* mti==N+1 means mt[N] is not initialized */
        Seed(seed);
    }

    void Seed(u_long seed) const;
    float RandomFloat() const;
    unsigned long RandomUInt() const;

private:
    static const int N = 624;
    mutable unsigned long mt[N]; /* the array for the state vector  */
    mutable int mti;
};



// Thatcher Ulrich (tu(at)tulrich(.)com
class TURNG {
public:
    // PRNG code adapted from the complimentary-multiply-with-carry
    // code in the article: George Marsaglia, "Seeds for Random Number
    // Generators", Communications of the ACM, May 2003, Vol 46 No 5,
    // pp90-93.
    //
    // The article says:
    //
    // "Any one of the choices for seed table size and multiplier will
    // provide a RNG that has passed extensive tests of randomness,
    // particularly those in [3], yet is simple and fast --
    // approximately 30 million random 32-bit integers per second on a
    // 850MHz PC.  The period is a*b^n, where a is the multiplier, n
    // the size of the seed table and b=2^32-1.  (a is chosen so that
    // b is a primitive root of the prime a*b^n + 1.)"
    //
    // [3] Marsaglia, G., Zaman, A., and Tsang, W.  Toward a universal
    // random number generator.  _Statistics and Probability Letters
    // 8_ (1990), 35-39.

//    const Uint64    a = 18782;    // for SEED_COUNT=4096, period approx 2^131104 (from Marsaglia usenet post 2003-05-13)
//    const Uint64    a = 123471786;    // for SEED_COUNT=1024, period approx 2^32794
//    const Uint64    a = 123554632;    // for SEED_COUNT=512, period approx 2^16410
//    const Uint64    a = 8001634;    // for SEED_COUNT=256, period approx 2^8182
//    const Uint64    a = 8007626;    // for SEED_COUNT=128, period approx 2^4118
//    const Uint64    a = 647535442;    // for SEED_COUNT=64, period approx 2^2077
//    const Uint64    a = 547416522;    // for SEED_COUNT=32, period approx 2^1053
//    const Uint64    a = 487198574;    // for SEED_COUNT=16, period approx  2^540

    static const int SEED_COUNT = 8;

    TURNG() {
        Seed(987654321);
    }


    void    Seed(u_int seed) {
        if (seed == 0) {
            // 0 is a terrible seed (probably the only bad
            // choice), substitute something else:
            seed = 12345;
        }

        // Simple pseudo-random to reseed the seeds.
        // Suggested by the above article.
        uint32_t j = seed;
        for (int i = 0; i < SEED_COUNT; i++) {
            j = j ^ (j << 13);
            j = j ^ (j >> 17);
            j = j ^ (j << 5);
            m_Q[i] = j;
        }

        m_c = 362436;
        m_i = SEED_COUNT - 1;
    }


    unsigned long RandomUInt() const {
        uint64_t t;
        uint32_t x;

        const uint64_t    a = 716514398;    // for SEED_COUNT=8, period approx 2^285

        //static uint32_t    c = 362436;
        //static uint32_t    i = SEED_COUNT - 1;
        const uint32_t    r = 0xFFFFFFFE;

        m_i = (m_i + 1) & (SEED_COUNT - 1);
        t = a * m_Q[m_i] + m_c;
        m_c = (uint32_t) (t >> 32);
        x = (uint32_t) (t + m_c);
        if (x < m_c) {
            x++;
            m_c++;
        }

        uint32_t    val = r - x;
        m_Q[m_i] = val;
        return val;
    }

    float RandomFloat() {
        uint32_t    r = RandomUInt();
        // 24 bits of precision.
        return float(r >> 8) / (16777216.0f - 1.0f);
    }

private:
    mutable uint32_t m_Q[SEED_COUNT];
    mutable uint32_t m_c;
    mutable uint32_t m_i;
};



#endif // PBRT_CORE_RNG_H
