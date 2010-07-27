
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


// integrators/whitted.cpp*
#include "stdafx.h"
#include "integrators/whitted.h"
#include "intersection.h"
#include "paramset.h"

// WhittedIntegrator Method Definitions
Spectrum WhittedIntegrator::Li(const Scene *scene,
        const Renderer *renderer, const RayDifferential &ray,
        const Intersection &isect, const Sample *sample, RNG &rng,
        MemoryArena &arena) const {
    Spectrum L(0.);
    // Compute emitted and reflected light at ray intersection point

    // Evaluate BSDF at hit point
    BSDF *bsdf = isect.GetBSDF(ray, arena);

    // Initialize common variables for Whitted integrator
    const Point &p = bsdf->dgShading.p;
    const Normal &n = bsdf->dgShading.nn;
    Vector wo = -ray.d;

    // Compute emitted light if ray hit an area light source
    L += isect.Le(wo);

    // Add contribution of each light source
    for (uint32_t i = 0; i < scene->lights.size(); ++i) {
        Vector wi;
        float pdf;
        VisibilityTester visibility;
        Spectrum Li = scene->lights[i]->Sample_L(p, isect.rayEpsilon,
            LightSample(rng), ray.time, &wi, &pdf, &visibility);
        if (Li.IsBlack() || pdf == 0.f) continue;
        Spectrum f = bsdf->f(wo, wi);
        if (!f.IsBlack() && visibility.Unoccluded(scene))
            L += f * Li * AbsDot(wi, n) *
                 visibility.Transmittance(scene, renderer,
                                          sample, rng, arena) / pdf;
    }
    if (ray.depth + 1 < maxDepth) {
        // Trace rays for specular reflection and refraction
        L += SpecularReflect(ray, bsdf, rng, isect, renderer, scene, sample,
                             arena);
        L += SpecularTransmit(ray, bsdf, rng, isect, renderer, scene, sample,
                              arena);
    }
    return L;
}


WhittedIntegrator *CreateWhittedSurfaceIntegrator(const ParamSet &params)
{
    int maxDepth = params.FindOneInt("maxdepth", 5);
    return new WhittedIntegrator(maxDepth);
}


