
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
