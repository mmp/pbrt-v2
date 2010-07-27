
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
