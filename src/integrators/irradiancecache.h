
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

#ifndef PBRT_INTEGRATORS_IRRADIANCECACHE_H
#define PBRT_INTEGRATORS_IRRADIANCECACHE_H

// integrators/irradiancecache.h*
#include "pbrt.h"
#include "integrator.h"
#include "octree.h"
#include "parallel.h"

// IrradianceCacheIntegrator Forward Declarations
struct IrradianceSample;

// IrradianceCacheIntegrator Declarations
class IrradianceCacheIntegrator : public SurfaceIntegrator {
public:
    // IrradianceCacheIntegrator Public Methods
    IrradianceCacheIntegrator(float minwt, float minsp, float maxsp,
                              float maxang, int maxspec, int maxind, int ns) {
        minWeight = minwt;
        minSamplePixelSpacing = minsp;
        maxSamplePixelSpacing = maxsp;
        cosMaxSampleAngleDifference = cosf(Degrees(maxang));
        nSamples = ns;
        maxSpecularDepth = maxspec;
        maxIndirectDepth = maxind;
        mutex = RWMutex::Create();
        lightSampleOffsets = NULL;
        bsdfSampleOffsets = NULL;
    }
    ~IrradianceCacheIntegrator();
    Spectrum Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential &ray, const Intersection &isect,
        const Sample *sample, RNG &rng, MemoryArena &arena) const;
    void RequestSamples(Sampler *sampler, Sample *sample, const Scene *scene);
    void Preprocess(const Scene *, const Camera *, const Renderer *);
private:
    // IrradianceCacheIntegrator Private Data
    float minSamplePixelSpacing, maxSamplePixelSpacing;
    float minWeight, cosMaxSampleAngleDifference;
    int nSamples, maxSpecularDepth, maxIndirectDepth;
    mutable RWMutex *mutex;

    // Declare sample parameters for light source sampling
    LightSampleOffsets *lightSampleOffsets;
    BSDFSampleOffsets *bsdfSampleOffsets;
    mutable Octree<IrradianceSample *> *octree;

    // IrradianceCacheIntegrator Private Methods
    Spectrum indirectLo(const Point &p, const Normal &ng, float pixelSpacing,
        const Vector &wo, float rayEpsilon,BSDF *bsdf, BxDFType flags, RNG &rng,
        const Scene *scene, const Renderer *renderer, MemoryArena &arena) const;
    bool interpolateE(const Scene *scene,
            const Point &p, const Normal &n, Spectrum *E, Vector *wi) const;
    Spectrum pathL(Ray &r, const Scene *scene, const Renderer *renderer,
        RNG &rng, MemoryArena &arena) const;
};


IrradianceCacheIntegrator *CreateIrradianceCacheIntegrator(const ParamSet &params);

#endif // PBRT_INTEGRATORS_IRRADIANCECACHE_H
