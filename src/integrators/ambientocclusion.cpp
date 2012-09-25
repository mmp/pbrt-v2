
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


// integrators/ambientocclusion.cpp*
#include "stdafx.h"
#include "integrators/ambientocclusion.h"
#include "paramset.h"
#include "montecarlo.h"
#include "scene.h"
#include "intersection.h"

// AmbientOcclusionIntegrator Method Definitions
Spectrum AmbientOcclusionIntegrator::Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential &ray, const Intersection &isect,
        const Sample *sample, RNG &rng, MemoryArena &arena) const {

    BSDF *bsdf = isect.GetBSDF(ray, arena);
    const Point &p = bsdf->dgShading.p;
    Normal n = Faceforward(isect.dg.nn, -ray.d);

    uint32_t scramble[2] = { rng.RandomUInt(), rng.RandomUInt() };
    float u[2];
    int nClear = 0;
    for (int i = 0; i < nSamples; ++i) {
        Sample02(i, scramble, u);
        Vector w = UniformSampleSphere(u[0], u[1]);
        if (Dot(w, n) < 0.) w = -w;
        Ray r(p, w, .01f, maxDist);
        if (!scene->IntersectP(r)) ++nClear;
    }
    return Spectrum(float(nClear) / float(nSamples));
}


AmbientOcclusionIntegrator *CreateAmbientOcclusionIntegrator(const ParamSet &params) {
    int nSamples = params.FindOneInt("nsamples", 2048);
    float maxDist = params.FindOneFloat("maxdist", INFINITY);
    if (PbrtOptions.quickRender) { nSamples = max(1, nSamples / 4); }
    return new AmbientOcclusionIntegrator(nSamples, maxDist);
}


