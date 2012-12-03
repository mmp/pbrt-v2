
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


// integrators/diffuseprt.cpp*
#include "stdafx.h"
#include "integrators/diffuseprt.h"
#include "sh.h"
#include "light.h"
#include "scene.h"
#include "camera.h"
#include "intersection.h"
#include "paramset.h"
#include "montecarlo.h"

// DiffusePRTIntegrator Method Definitions
DiffusePRTIntegrator::DiffusePRTIntegrator(int lm, int ns)
    : lmax(lm), nSamples(RoundUpPow2(ns)) {
    c_in = new Spectrum[SHTerms(lmax)];
}


DiffusePRTIntegrator::~DiffusePRTIntegrator() {
    delete[] c_in;
}


void DiffusePRTIntegrator::Preprocess(const Scene *scene,
        const Camera *camera, const Renderer *renderer) {
    BBox bbox = scene->WorldBound();
    Point p = .5f * bbox.pMin + .5f * bbox.pMax;
    RNG rng;
    MemoryArena arena;
    SHProjectIncidentDirectRadiance(p, 0.f, camera->shutterOpen, arena,
                                    scene, false, lmax, rng, c_in);
}


void DiffusePRTIntegrator::RequestSamples(Sampler *sampler, Sample *sample, const Scene *scene) {
}


Spectrum DiffusePRTIntegrator::Li(const Scene *scene, const Renderer *,
            const RayDifferential &ray, const Intersection &isect,
            const Sample *sample, RNG &rng, MemoryArena &arena) const {
    Spectrum L = 0.f;
    Vector wo = -ray.d;
    // Compute emitted light if ray hit an area light source
    L += isect.Le(wo);

    // Evaluate BSDF at hit point
    BSDF *bsdf = isect.GetBSDF(ray, arena);
    const Point &p = bsdf->dgShading.p;
    const Normal &n = bsdf->dgShading.nn;
    // Compute reflected radiance using diffuse PRT

    // Project diffuse transfer function at point to SH
    Spectrum *c_transfer = arena.Alloc<Spectrum>(SHTerms(lmax));
    SHComputeDiffuseTransfer(p, Faceforward(n, wo), isect.rayEpsilon,
                             scene, rng, nSamples, lmax, c_transfer);

    // Compute integral of product of incident radiance and transfer function
    Spectrum Kd = bsdf->rho(wo, rng, BSDF_ALL_REFLECTION) * INV_PI;
    Spectrum Lo = 0.f;
    for (int i = 0; i < SHTerms(lmax); ++i)
        Lo += c_in[i] * c_transfer[i];
    return L + Kd * Lo.Clamp();
}


DiffusePRTIntegrator *CreateDiffusePRTIntegratorSurfaceIntegrator(const ParamSet &params) {
    int lmax = params.FindOneInt("lmax", 4);
    int ns = params.FindOneInt("nsamples", 4096);
    return new DiffusePRTIntegrator(lmax, ns);
}


