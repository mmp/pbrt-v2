
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


