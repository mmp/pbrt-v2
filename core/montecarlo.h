
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

#ifndef PBRT_CORE_MONTECARLO_H
#define PBRT_CORE_MONTECARLO_H

// core/montecarlo.h*
#include "pbrt.h"
#include "geometry.h"
#include "rng.h"

// MC Utility Declarations
void ComputeStep1dCDF(const float *f, int nValues, float *c, float *cdf);
float SampleStep1d(const float *f, const float *cdf, float c,
    int nSteps, float u, float *weight);
void RejectionSampleDisk(float u1, float u2, float *x, float *y);
Vector UniformSampleHemisphere(float u1, float u2);
float  UniformHemispherePdf();
Vector UniformSampleSphere(float u1, float u2);
float  UniformSpherePdf();
Vector UniformSampleCone(float u1, float u2, float thetamax);
Vector UniformSampleCone(float u1, float u2, float thetamax,
    const Vector &x, const Vector &y, const Vector &z);
float  UniformConePdf(float thetamax);
void UniformSampleDisk(float u1, float u2, float *x, float *y);
void ConcentricSampleDisk(float u1, float u2, float *dx, float *dy);
inline Vector CosineSampleHemisphere(float u1, float u2) {
    Vector ret;
    ConcentricSampleDisk(u1, u2, &ret.x, &ret.y);
    ret.z = sqrtf(max(0.f, 1.f - ret.x*ret.x - ret.y*ret.y));
    return ret;
}


inline float CosineHemispherePdf(float costheta, float phi) {
    return costheta * INV_PI;
}


void UniformSampleTriangle(float ud1, float ud2, float *u, float *v);
void StratifiedSample1D(float *samples, int nsamples, RNG &rng,
                        bool jitter = true);
void StratifiedSample2D(float *samples, int nx, int ny, RNG &rng,
                        bool jitter = true);
void Shuffle(float *samp, u_int count, u_int dims, RNG &rng);
void LatinHypercube(float *samples, u_int nSamples, u_int nDim, RNG &rng);
inline double RadicalInverse(int n, int base) {
    double val = 0;
    double invBase = 1. / base, invBi = invBase;
    while (n > 0) {
        // Compute next digit of radical inverse
        int d_i = (n % base);
        val += d_i * invBi;
        n *= invBase;
        invBi *= invBase;
    }
    return val;
}


inline double FoldedRadicalInverse(int n, int base) {
    double val = 0;
    double invBase = 1.f/base, invBi = invBase;
    int modOffset = 0;
    while (val + base * invBi != val) {
        // Compute next digit of folded radical inverse
        int digit = ((n+modOffset) % base);
        val += digit * invBi;
        n *= invBase;
        invBi *= invBase;
        ++modOffset;
    }
    return val;
}


inline void GeneratePermutation(u_int *buf, u_int n, RNG &rng) {
    for (u_int i = 0; i < n; ++i)
        buf[i] = i;
    for (u_int i = 0; i < n; ++i) {
        u_int other = i + (rng.RandomUInt() % (n - i));
        swap(buf[i], buf[other]);
    }
}


inline double PermutedRadicalInverse(u_int n, u_int base, const u_int *permute) {
    double val = 0;
    double invBase = 1. / base, invBi = invBase;
    while (n > 0) {
        u_int d_i = (n % base);
        d_i = permute[d_i];
        val += d_i * invBi;
        n *= invBase;
        invBi *= invBase;
    }
    return val;
}


class PermutedHalton {
public:
    // PermutedHalton Public Methods
    PermutedHalton(u_int d, RNG &rng);
    ~PermutedHalton() {
        delete[] bases;
        delete[] permute;
    }
    void Sample(u_int n, float *out) const {
        u_int *p = permute;
        for (u_int i = 0; i < dims; ++i) {
            out[i] = PermutedRadicalInverse(n, bases[i], p);
            p += bases[i];
        }
    }
private:
    // PermutedHalton Private Data
    u_int dims;
    u_int *bases, *permute;
    PermutedHalton(const PermutedHalton &);
    PermutedHalton &operator=(const PermutedHalton &);
};


inline float VanDerCorput(u_int n, u_int scramble = 0);
inline float Sobol2(u_int n, u_int scramble = 0);
inline float LarcherPillichshammer2(u_int n, u_int scramble = 0);
inline void Sample02(u_int n, u_int scramble[2], float sample[2]);
int LDPixelSampleFloatsNeeded(const Sample *sample, int pixelSamples);
void LDPixelSample(int xPos, int yPos, float ShutterOpen,
    float ShutterClose, int pixelSamples, Sample *samples, float *buf);
void SampleBlinn(const Vector &wo, Vector *wi, float u1, float u2, float *pdf, float exponent);
float BlinnPdf(const Vector &wo, const Vector &wi, float exponent);
struct Distribution1D {
    // Distribution1D Methods
    Distribution1D(float *f, int n) {
        func = new float[n];
        cdf = new float[n+1];
        count = n;
        memcpy(func, f, n*sizeof(float));
        ComputeStep1dCDF(func, n, &funcInt, cdf);
        invFuncInt = 1.f / funcInt;
        invCount = 1.f / count;
    }
    ~Distribution1D() {
        delete[] func;
        delete[] cdf;
    }
    float Sample(float u, float *pdf) {
        // Find surrounding cdf segments
        float *ptr = std::lower_bound(cdf, cdf+count+1, u);
        int offset = max(0, (int) (ptr-cdf-1));
        // Return offset along current cdf segment
        u = (u - cdf[offset]) / (cdf[offset+1] - cdf[offset]);
        *pdf = func[offset] * invFuncInt;
        return offset + u;
    }

    // Distribution1D Data
    float *func, *cdf;
    float funcInt, invFuncInt, invCount;
    int count;
};


Vector SampleHG(const Vector &w, float g, float u1, float u2);
float HGPdf(const Vector &w, const Vector &wp, float g);

// MC Inline Functions
inline float BalanceHeuristic(int nf, float fPdf, int ng, float gPdf) {
    return (nf * fPdf) / (nf * fPdf + ng * gPdf);
}


inline float PowerHeuristic(int nf, float fPdf, int ng, float gPdf) {
    float f = nf * fPdf, g = ng * gPdf;
    return (f*f) / (f*f + g*g);
}



// Sampling Inline Functions
inline void Sample02(u_int n, u_int scramble[2], float sample[2]) {
    sample[0] = VanDerCorput(n, scramble[0]);
    sample[1] = Sobol2(n, scramble[1]);
}


inline float VanDerCorput(u_int n, u_int scramble) {
    n = (n << 16) | (n >> 16);
    n = ((n & 0x00ff00ff) << 8) | ((n & 0xff00ff00) >> 8);
    n = ((n & 0x0f0f0f0f) << 4) | ((n & 0xf0f0f0f0) >> 4);
    n = ((n & 0x33333333) << 2) | ((n & 0xcccccccc) >> 2);
    n = ((n & 0x55555555) << 1) | ((n & 0xaaaaaaaa) >> 1);
    n ^= scramble;
    return (float)n / (float)0x100000000LL;
}


inline float Sobol2(u_int n, u_int scramble) {
    for (u_int v = 1 << 31; n != 0; n >>= 1, v ^= v >> 1)
        if (n & 0x1) scramble ^= v;
    return (float)scramble / (float)0x100000000LL;
}


inline float
LarcherPillichshammer2(u_int n, u_int scramble) {
    for (u_int v = 1 << 31; n != 0; n >>= 1, v |= v >> 1)
        if (n & 0x1) scramble ^= v;
    return (float)scramble / (float)0x100000000LL;
}


inline void LDShuffleScrambled1D(int nSamples,
        int nPixel, float *samples, RNG &rng) {
    u_int scramble = rng.RandomUInt();
    for (int i = 0; i < nSamples * nPixel; ++i)
        samples[i] = VanDerCorput(i, scramble);
    for (int i = 0; i < nPixel; ++i)
        Shuffle(samples + i * nSamples, nSamples, 1, rng);
    Shuffle(samples, nPixel, nSamples, rng);
}


inline void LDShuffleScrambled2D(int nSamples,
        int nPixel, float *samples, RNG &rng) {
    u_int scramble[2] = { rng.RandomUInt(), rng.RandomUInt() };
    for (int i = 0; i < nSamples * nPixel; ++i)
        Sample02(i, scramble, &samples[2*i]);
    for (int i = 0; i < nPixel; ++i)
        Shuffle(samples + 2 * i * nSamples, nSamples, 2, rng);
    Shuffle(samples, nPixel, 2 * nSamples, rng);
}



#endif // PBRT_CORE_MONTECARLO_H
