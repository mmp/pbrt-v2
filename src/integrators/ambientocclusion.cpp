
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


