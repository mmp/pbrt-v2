
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


