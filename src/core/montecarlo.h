
/*
    pbrt source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.

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

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CORE_MONTECARLO_H
#define PBRT_CORE_MONTECARLO_H

// core/montecarlo.h*
#include "pbrt.h"
#include "geometry.h"
#include "rng.h"

// smallest floating point value less than one; all canonical random samples
// should be <= this.
#ifdef PBRT_IS_WINDOWS
// sadly, MSVC2008 (at least) doesn't support hexidecimal fp constants...
static const float OneMinusEpsilon=0.9999999403953552f;
#else
static const float OneMinusEpsilon=0x1.fffffep-1;
#endif

// Monte Carlo Utility Declarations
struct Distribution1D {
    // Distribution1D Public Methods
    Distribution1D(const float *f, int n) {
        count = n;
        func = new float[n];
        memcpy(func, f, n*sizeof(float));
        cdf = new float[n+1];
        // Compute integral of step function at $x_i$
        cdf[0] = 0.;
        for (int i = 1; i < count+1; ++i)
            cdf[i] = cdf[i-1] + func[i-1] / n;

        // Transform step function integral into CDF
        funcInt = cdf[count];
        if (funcInt == 0.f) {
            for (int i = 1; i < n+1; ++i)
                cdf[i] = float(i) / float(n);
        }
        else {
            for (int i = 1; i < n+1; ++i)
                cdf[i] /= funcInt;
        }
    }
    ~Distribution1D() {
        delete[] func;
        delete[] cdf;
    }
    float SampleContinuous(float u, float *pdf, int *off = NULL) const {
        // Find surrounding CDF segments and _offset_
        float *ptr = std::upper_bound(cdf, cdf+count+1, u);
        int offset = max(0, int(ptr-cdf-1));
        if (off) *off = offset;
        Assert(offset < count);
        Assert(u >= cdf[offset] && u < cdf[offset+1]);

        // Compute offset along CDF segment
        float du = (u - cdf[offset]) / (cdf[offset+1] - cdf[offset]);
        Assert(!isnan(du));

        // Compute PDF for sampled offset
        if (pdf) *pdf = func[offset] / funcInt;

        // Return $x\in{}[0,1)$ corresponding to sample
        return (offset + du) / count;
    }
    int SampleDiscrete(float u, float *pdf) const {
        // Find surrounding CDF segments and _offset_
        float *ptr = std::upper_bound(cdf, cdf+count+1, u);
        int offset = max(0, int(ptr-cdf-1));
        Assert(offset < count);
        Assert(u >= cdf[offset] && u < cdf[offset+1]);
        if (pdf) *pdf = func[offset] / (funcInt * count);
        return offset;
    }
private:
    friend struct Distribution2D;
    // Distribution1D Private Data
    float *func, *cdf;
    float funcInt;
    int count;
};


void RejectionSampleDisk(float *x, float *y, RNG &rng);
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
struct Distribution2D {
    // Distribution2D Public Methods
    Distribution2D(const float *data, int nu, int nv);
    ~Distribution2D();
    void SampleContinuous(float u0, float u1, float uv[2],
                          float *pdf) const {
        float pdfs[2];
        int v;
        uv[1] = pMarginal->SampleContinuous(u1, &pdfs[1], &v);
        uv[0] = pConditionalV[v]->SampleContinuous(u0, &pdfs[0]);
        *pdf = pdfs[0] * pdfs[1];
    }
    float Pdf(float u, float v) const {
        int iu = Clamp(Float2Int(u * pConditionalV[0]->count), 0,
                       pConditionalV[0]->count-1);
        int iv = Clamp(Float2Int(v * pMarginal->count), 0,
                       pMarginal->count-1);
        if (pConditionalV[iv]->funcInt * pMarginal->funcInt == 0.f) return 0.f;
        return (pConditionalV[iv]->func[iu] * pMarginal->func[iv]) /
               (pConditionalV[iv]->funcInt * pMarginal->funcInt);
    }
private:
    // Distribution2D Private Data
    vector<Distribution1D *> pConditionalV;
    Distribution1D *pMarginal;
};


void StratifiedSample1D(float *samples, int nsamples, RNG &rng,
                        bool jitter = true);
void StratifiedSample2D(float *samples, int nx, int ny, RNG &rng,
                        bool jitter = true);
template <typename T>
void Shuffle(T *samp, uint32_t count, uint32_t dims, RNG &rng) {
    for (uint32_t i = 0; i < count; ++i) {
        uint32_t other = i + (rng.RandomUInt() % (count - i));
        for (uint32_t j = 0; j < dims; ++j)
            swap(samp[dims*i + j], samp[dims*other + j]);
    }
}


void LatinHypercube(float *samples, uint32_t nSamples, uint32_t nDim, RNG &rng);
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


inline void GeneratePermutation(uint32_t *buf, uint32_t b, RNG &rng) {
    for (uint32_t i = 0; i < b; ++i)
        buf[i] = i;
    Shuffle(buf, b, 1, rng);
}


inline double PermutedRadicalInverse(uint32_t n, uint32_t base,
                                     const uint32_t *p) {
    double val = 0;
    double invBase = 1. / base, invBi = invBase;

    while (n > 0) {
        uint32_t d_i = p[n % base];
        val += d_i * invBi;
        n *= invBase;
        invBi *= invBase;
    }
    return val;
}


class PermutedHalton {
public:
    // PermutedHalton Public Methods
    PermutedHalton(uint32_t d, RNG &rng);
    ~PermutedHalton() {
        delete[] b;
        delete[] permute;
    }
    void Sample(uint32_t n, float *out) const {
        uint32_t *p = permute;
        for (uint32_t i = 0; i < dims; ++i) {
            out[i] = min(float(PermutedRadicalInverse(n, b[i], p)), 
                         OneMinusEpsilon);
            p += b[i];
        }
    }
private:
    // PermutedHalton Private Data
    uint32_t dims;
    uint32_t *b, *permute;
    PermutedHalton(const PermutedHalton &);
    PermutedHalton &operator=(const PermutedHalton &);
};


inline float VanDerCorput(uint32_t n, uint32_t scramble = 0);
inline float Sobol2(uint32_t n, uint32_t scramble = 0);
inline float LarcherPillichshammer2(uint32_t n, uint32_t scramble = 0);
inline void Sample02(uint32_t n, const uint32_t scramble[2], float sample[2]);
int LDPixelSampleFloatsNeeded(const Sample *sample, int nPixelSamples);
void LDPixelSample(int xPos, int yPos, float shutterOpen,
    float shutterClose, int nPixelSamples, Sample *samples, float *buf, RNG &rng);
Vector SampleHG(const Vector &w, float g, float u1, float u2);
float HGPdf(const Vector &w, const Vector &wp, float g);

// Monte Carlo Inline Functions
inline float BalanceHeuristic(int nf, float fPdf, int ng, float gPdf) {
    return (nf * fPdf) / (nf * fPdf + ng * gPdf);
}


inline float PowerHeuristic(int nf, float fPdf, int ng, float gPdf) {
    float f = nf * fPdf, g = ng * gPdf;
    return (f*f) / (f*f + g*g);
}



// Sampling Inline Functions
inline void Sample02(uint32_t n, const uint32_t scramble[2],
                     float sample[2]) {
    sample[0] = VanDerCorput(n, scramble[0]);
    sample[1] = Sobol2(n, scramble[1]);
}


inline float VanDerCorput(uint32_t n, uint32_t scramble) {
    // Reverse bits of _n_
    n = (n << 16) | (n >> 16);
    n = ((n & 0x00ff00ff) << 8) | ((n & 0xff00ff00) >> 8);
    n = ((n & 0x0f0f0f0f) << 4) | ((n & 0xf0f0f0f0) >> 4);
    n = ((n & 0x33333333) << 2) | ((n & 0xcccccccc) >> 2);
    n = ((n & 0x55555555) << 1) | ((n & 0xaaaaaaaa) >> 1);
    n ^= scramble;
    return min(((n>>8) & 0xffffff) / float(1 << 24), OneMinusEpsilon);
}


inline float Sobol2(uint32_t n, uint32_t scramble) {
    for (uint32_t v = 1 << 31; n != 0; n >>= 1, v ^= v >> 1)
        if (n & 0x1) scramble ^= v;
    return min(((scramble>>8) & 0xffffff) / float(1 << 24), OneMinusEpsilon);
}


inline float
LarcherPillichshammer2(uint32_t n, uint32_t scramble) {
    for (uint32_t v = 1 << 31; n != 0; n >>= 1, v |= v >> 1)
        if (n & 0x1) scramble ^= v;
    return min(((scramble>>8) & 0xffffff) / float(1 << 24), OneMinusEpsilon);
}


inline void LDShuffleScrambled1D(int nSamples, int nPixel,
                                 float *samples, RNG &rng) {
    uint32_t scramble = rng.RandomUInt();
    for (int i = 0; i < nSamples * nPixel; ++i)
        samples[i] = VanDerCorput(i, scramble);
    for (int i = 0; i < nPixel; ++i)
        Shuffle(samples + i * nSamples, nSamples, 1, rng);
    Shuffle(samples, nPixel, nSamples, rng);
}


inline void LDShuffleScrambled2D(int nSamples, int nPixel,
                                 float *samples, RNG &rng) {
    uint32_t scramble[2] = { rng.RandomUInt(), rng.RandomUInt() };
    for (int i = 0; i < nSamples * nPixel; ++i)
        Sample02(i, scramble, &samples[2*i]);
    for (int i = 0; i < nPixel; ++i)
        Shuffle(samples + 2 * i * nSamples, nSamples, 2, rng);
    Shuffle(samples, nPixel, 2 * nSamples, rng);
}



#endif // PBRT_CORE_MONTECARLO_H
