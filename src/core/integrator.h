
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

#ifndef PBRT_CORE_INTEGRATOR_H
#define PBRT_CORE_INTEGRATOR_H

// core/integrator.h*
#include "pbrt.h"
#include "primitive.h"
#include "spectrum.h"
#include "light.h"
#include "reflection.h"
#include "sampler.h"
#include "material.h"
#include "probes.h"
#include "renderer.h"

// Integrator Declarations
class Integrator {
public:
    // Integrator Interface
    virtual ~Integrator();
    virtual void Preprocess(const Scene *scene, const Camera *camera,
                            const Renderer *renderer) {
    }
    virtual void RequestSamples(Sampler *sampler, Sample *sample,
                                const Scene *scene) {
    }
};


class SurfaceIntegrator : public Integrator {
public:
    // SurfaceIntegrator Interface
    virtual Spectrum Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential &ray, const Intersection &isect,
        const Sample *sample, RNG &rng, MemoryArena &arena) const = 0;
};


Spectrum UniformSampleAllLights(const Scene *scene, const Renderer *renderer,
    MemoryArena &arena, const Point &p, const Normal &n, const Vector &wo,
    float rayEpsilon, float time, BSDF *bsdf, const Sample *sample, RNG &rng,
    const LightSampleOffsets *lightSampleOffsets,
    const BSDFSampleOffsets *bsdfSampleOffsets);
Spectrum UniformSampleOneLight(const Scene *scene, const Renderer *renderer,
    MemoryArena &arena, const Point &p, const Normal &n, const Vector &wo,
    float rayEpsilon, float time, BSDF *bsdf,
    const Sample *sample, RNG &rng, int lightNumOffset = -1,
    const LightSampleOffsets *lightSampleOffset = NULL,
    const BSDFSampleOffsets *bsdfSampleOffset = NULL);
Spectrum EstimateDirect(const Scene *scene, const Renderer *renderer,
    MemoryArena &arena, const Light *light, const Point &p,
    const Normal &n, const Vector &wo, float rayEpsilon, float time, const BSDF *bsdf,
    RNG &rng, const LightSample &lightSample, const BSDFSample &bsdfSample,
    BxDFType flags);
Spectrum SpecularReflect(const RayDifferential &ray, BSDF *bsdf, RNG &rng,
    const Intersection &isect, const Renderer *renderer, const Scene *scene,
    const Sample *sample, MemoryArena &arena);
Spectrum SpecularTransmit(const RayDifferential &ray, BSDF *bsdf, RNG &rng,
    const Intersection &isect, const Renderer *renderer, const Scene *scene,
    const Sample *sample, MemoryArena &arena);
Distribution1D *ComputeLightSamplingCDF(const Scene *scene);

#endif // PBRT_CORE_INTEGRATOR_H
