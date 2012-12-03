
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


// integrators/directlighting.cpp*
#include "stdafx.h"
#include "integrators/directlighting.h"
#include "intersection.h"
#include "paramset.h"

// DirectLightingIntegrator Method Definitions
DirectLightingIntegrator::DirectLightingIntegrator(LightStrategy st, int md) {
    maxDepth = md;
    strategy = st;
    lightSampleOffsets = NULL;
    bsdfSampleOffsets = NULL;
}


DirectLightingIntegrator::~DirectLightingIntegrator() {
    delete[] lightSampleOffsets;
    delete[] bsdfSampleOffsets;
}


void DirectLightingIntegrator::RequestSamples(Sampler *sampler,
        Sample *sample, const Scene *scene) {
    if (strategy == SAMPLE_ALL_UNIFORM) {
        // Allocate and request samples for sampling all lights
        uint32_t nLights = scene->lights.size();
        lightSampleOffsets = new LightSampleOffsets[nLights];
        bsdfSampleOffsets = new BSDFSampleOffsets[nLights];
        for (uint32_t i = 0; i < nLights; ++i) {
            const Light *light = scene->lights[i];
            int nSamples = light->nSamples;
            if (sampler) nSamples = sampler->RoundSize(nSamples);
            lightSampleOffsets[i] = LightSampleOffsets(nSamples, sample);
            bsdfSampleOffsets[i] = BSDFSampleOffsets(nSamples, sample);
        }
        lightNumOffset = -1;
    }
    else {
        // Allocate and request samples for sampling one light
        lightSampleOffsets = new LightSampleOffsets[1];
        lightSampleOffsets[0] = LightSampleOffsets(1, sample);
        lightNumOffset = sample->Add1D(1);
        bsdfSampleOffsets = new BSDFSampleOffsets[1];
        bsdfSampleOffsets[0] = BSDFSampleOffsets(1, sample);
    }
}


Spectrum DirectLightingIntegrator::Li(const Scene *scene,
        const Renderer *renderer, const RayDifferential &ray,
        const Intersection &isect, const Sample *sample, RNG &rng, MemoryArena &arena) const {
    Spectrum L(0.f);
    // Evaluate BSDF at hit point
    BSDF *bsdf = isect.GetBSDF(ray, arena);
    Vector wo = -ray.d;
    const Point &p = bsdf->dgShading.p;
    const Normal &n = bsdf->dgShading.nn;
    // Compute emitted light if ray hit an area light source
    L += isect.Le(wo);

    // Compute direct lighting for _DirectLightingIntegrator_ integrator
    if (scene->lights.size() > 0) {
        // Apply direct lighting strategy
        switch (strategy) {
            case SAMPLE_ALL_UNIFORM:
                L += UniformSampleAllLights(scene, renderer, arena, p, n, wo,
                    isect.rayEpsilon, ray.time, bsdf, sample, rng,
                    lightSampleOffsets, bsdfSampleOffsets);
                break;
            case SAMPLE_ONE_UNIFORM:
                L += UniformSampleOneLight(scene, renderer, arena, p, n, wo,
                    isect.rayEpsilon, ray.time, bsdf, sample, rng,
                    lightNumOffset, lightSampleOffsets, bsdfSampleOffsets);
                break;
        }
    }
    if (ray.depth + 1 < maxDepth) {
        Vector wi;
        // Trace rays for specular reflection and refraction
        L += SpecularReflect(ray, bsdf, rng, isect, renderer, scene, sample,
                             arena);
        L += SpecularTransmit(ray, bsdf, rng, isect, renderer, scene, sample,
                              arena);
    }
    return L;
}


DirectLightingIntegrator *CreateDirectLightingIntegrator(const ParamSet &params) {
    int maxDepth = params.FindOneInt("maxdepth", 5);
    LightStrategy strategy;
    string st = params.FindOneString("strategy", "all");
    if (st == "one") strategy = SAMPLE_ONE_UNIFORM;
    else if (st == "all") strategy = SAMPLE_ALL_UNIFORM;
    else {
        Warning("Strategy \"%s\" for direct lighting unknown. "
            "Using \"all\".", st.c_str());
        strategy = SAMPLE_ALL_UNIFORM;
    }
    return new DirectLightingIntegrator(strategy, maxDepth);
}


