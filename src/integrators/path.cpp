
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


// integrators/path.cpp*
#include "stdafx.h"
#include "integrators/path.h"
#include "scene.h"
#include "intersection.h"
#include "paramset.h"

// PathIntegrator Method Definitions
void PathIntegrator::RequestSamples(Sampler *sampler, Sample *sample,
                                    const Scene *scene) {
    for (int i = 0; i < SAMPLE_DEPTH; ++i) {
        lightSampleOffsets[i] = LightSampleOffsets(1, sample);
        lightNumOffset[i] = sample->Add1D(1);
        bsdfSampleOffsets[i] = BSDFSampleOffsets(1, sample);
        pathSampleOffsets[i] = BSDFSampleOffsets(1, sample);
    }
}


Spectrum PathIntegrator::Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential &r, const Intersection &isect,
        const Sample *sample, RNG &rng, MemoryArena &arena) const {
    // Declare common path integration variables
    Spectrum pathThroughput = 1., L = 0.;
    RayDifferential ray(r);
    bool specularBounce = false;
    Intersection localIsect;
    const Intersection *isectp = &isect;
    for (int bounces = 0; ; ++bounces) {
        // Possibly add emitted light at path vertex
        if (bounces == 0 || specularBounce)
            L += pathThroughput * isectp->Le(-ray.d);

        // Sample illumination from lights to find path contribution
        BSDF *bsdf = isectp->GetBSDF(ray, arena);
        const Point &p = bsdf->dgShading.p;
        const Normal &n = bsdf->dgShading.nn;
        Vector wo = -ray.d;
        if (bounces < SAMPLE_DEPTH)
            L += pathThroughput *
                 UniformSampleOneLight(scene, renderer, arena, p, n, wo,
                     isectp->rayEpsilon, ray.time, bsdf, sample, rng,
                     lightNumOffset[bounces], &lightSampleOffsets[bounces],
                     &bsdfSampleOffsets[bounces]);
        else
            L += pathThroughput *
                 UniformSampleOneLight(scene, renderer, arena, p, n, wo,
                     isectp->rayEpsilon, ray.time, bsdf, sample, rng);

        // Sample BSDF to get new path direction

        // Get _outgoingBSDFSample_ for sampling new path direction
        BSDFSample outgoingBSDFSample;
        if (bounces < SAMPLE_DEPTH)
            outgoingBSDFSample = BSDFSample(sample, pathSampleOffsets[bounces],
                                            0);
        else
            outgoingBSDFSample = BSDFSample(rng);
        Vector wi;
        float pdf;
        BxDFType flags;
        Spectrum f = bsdf->Sample_f(wo, &wi, outgoingBSDFSample, &pdf,
                                    BSDF_ALL, &flags);
        if (f.IsBlack() || pdf == 0.)
            break;
        specularBounce = (flags & BSDF_SPECULAR) != 0;
        pathThroughput *= f * AbsDot(wi, n) / pdf;
        ray = RayDifferential(p, wi, ray, isectp->rayEpsilon);

        // Possibly terminate the path
        if (bounces > 3) {
            float continueProbability = min(.5f, pathThroughput.y());
            if (rng.RandomFloat() > continueProbability)
                break;
            pathThroughput /= continueProbability;
        }
        if (bounces == maxDepth)
            break;

        // Find next vertex of path
        if (!scene->Intersect(ray, &localIsect)) {
            if (specularBounce)
                for (uint32_t i = 0; i < scene->lights.size(); ++i)
                   L += pathThroughput * scene->lights[i]->Le(ray);
            break;
        }
        pathThroughput *= renderer->Transmittance(scene, ray, NULL, rng, arena);
        isectp = &localIsect;
    }
    return L;
}


PathIntegrator *CreatePathSurfaceIntegrator(const ParamSet &params) {
    int maxDepth = params.FindOneInt("maxdepth", 5);
    return new PathIntegrator(maxDepth);
}


