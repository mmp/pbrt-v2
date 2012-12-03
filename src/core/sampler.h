
/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CORE_SAMPLER_H
#define PBRT_CORE_SAMPLER_H

// core/sampler.h*
#include "pbrt.h"
#include "geometry.h"
#include "rng.h"
#include "memory.h"

// Sampling Declarations
class Sampler {
public:
    // Sampler Interface
    virtual ~Sampler();
    Sampler(int xstart, int xend, int ystart, int yend,
            int spp, float sopen, float sclose);
    virtual int GetMoreSamples(Sample *sample, RNG &rng) = 0;
    virtual int MaximumSampleCount() = 0;
    virtual bool ReportResults(Sample *samples, const RayDifferential *rays,
        const Spectrum *Ls, const Intersection *isects, int count);
    virtual Sampler *GetSubSampler(int num, int count) = 0;
    virtual int RoundSize(int size) const = 0;

    // Sampler Public Data
    const int xPixelStart, xPixelEnd, yPixelStart, yPixelEnd;
    const int samplesPerPixel;
    const float shutterOpen, shutterClose;
protected:
    // Sampler Protected Methods
    void ComputeSubWindow(int num, int count, int *xstart, int *xend, int *ystart, int *yend) const;
};


struct CameraSample {
    float imageX, imageY;
    float lensU, lensV;
    float time;
};


struct Sample : public CameraSample {
    // Sample Public Methods
    Sample(Sampler *sampler, SurfaceIntegrator *surf, VolumeIntegrator *vol,
        const Scene *scene);
    uint32_t Add1D(uint32_t num) {
        n1D.push_back(num);
        return n1D.size()-1;
    }
    uint32_t Add2D(uint32_t num) {
        n2D.push_back(num);
        return n2D.size()-1;
    }
    ~Sample() {
        if (oneD != NULL) {
            FreeAligned(oneD[0]);
            FreeAligned(oneD);
        }
    }
    Sample *Duplicate(int count) const;

    // Sample Public Data
    vector<uint32_t> n1D, n2D;
    float **oneD, **twoD;
private:
    // Sample Private Methods
    void AllocateSampleMemory();
    Sample() { oneD = twoD = NULL; }
};



#endif // PBRT_CORE_SAMPLER_H
