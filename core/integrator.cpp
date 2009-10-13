
/*
    pbrt source code Copyright(c) 1998-2009 Matt Pharr and Greg Humphreys.

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


// core/integrator.cpp*
#include "integrator.h"
#include "scene.h"
#include "intersection.h"
#include "montecarlo.h"

// Integrator Method Definitions
Integrator::~Integrator() {
}



// Integrator Utility Functions
Spectrum UniformSampleAllLights(const Scene *scene,
        const Renderer *renderer, MemoryArena &arena, const Point &p,
        const Normal &n, const Vector &wo, float rayEpsilon,
        BSDF *bsdf, const Sample *sample,
        const LightSampleOffsets *lightSampleOffsets,
        const BSDFSampleOffsets *bsdfSampleOffsets) {
    Spectrum L(0.);
    for (u_int i = 0; i < scene->lights.size(); ++i) {
        Light *light = scene->lights[i];
        int nSamples = lightSampleOffsets ?
                       lightSampleOffsets[i].nSamples : 1;
        // Estimate direct lighting from _light_ samples
        Spectrum Ld(0.);
        for (int j = 0; j < nSamples; ++j) {
            // Find light and BSDF sample values for direct lighting estimate
            LightSample lightSample;
            BSDFSample bsdfSample;
            if (lightSampleOffsets != NULL && bsdfSampleOffsets != NULL) {
                lightSample = LightSample(sample, lightSampleOffsets[i], j);
                bsdfSample = BSDFSample(sample, bsdfSampleOffsets[i], j);
            }
            else {
                lightSample = LightSample(*sample->rng);
                bsdfSample = BSDFSample(*sample->rng);
            }
            Ld += EstimateDirect(scene, renderer, arena, light, p, n, wo,
                rayEpsilon, sample->Time, bsdf, sample->rng, lightSample, bsdfSample);
        }
        L += Ld / nSamples;
    }
    return L;
}


Spectrum UniformSampleOneLight(const Scene *scene,
        const Renderer *renderer, MemoryArena &arena, const Point &p,
        const Normal &n, const Vector &wo, float rayEpsilon,
        BSDF *bsdf, const Sample *sample, int lightNumOffset,
        const LightSampleOffsets *lightSampleOffset,
        const BSDFSampleOffsets *bsdfSampleOffset) {
    // Randomly choose a single light to sample, _light_
    int nLights = int(scene->lights.size());
    if (nLights == 0.) return Spectrum(0.);
    int lightNum;
    if (lightNumOffset != -1)
        lightNum = Floor2Int(sample->oneD[lightNumOffset][0] * nLights);
    else
        lightNum = Floor2Int(sample->rng->RandomFloat() * nLights);
    lightNum = min(lightNum, nLights-1);
    Light *light = scene->lights[lightNum];

    // Initialize light and bsdf samples for single light sample
    LightSample lightSample;
    BSDFSample bsdfSample;
    if (lightSampleOffset != NULL && bsdfSampleOffset != NULL) {
        lightSample = LightSample(sample, *lightSampleOffset, 0);
        bsdfSample = BSDFSample(sample, *bsdfSampleOffset, 0);
    }
    else {
        lightSample = LightSample(*sample->rng);
        bsdfSample = BSDFSample(*sample->rng);
    }
    return (float)nLights *
        EstimateDirect(scene, renderer, arena, light, p, n, wo,
                       rayEpsilon, sample->Time, bsdf, sample->rng, lightSample, bsdfSample);
}


Spectrum EstimateDirect(const Scene *scene, const Renderer *renderer,
        MemoryArena &arena, const Light *light, const Point &p,
        const Normal &n, const Vector &wo, float rayEpsilon, float time,
        BSDF *bsdf, RNG *rng, const LightSample &lightSample,
        const BSDFSample &bsdfSample) {
    Spectrum Ld(0.);
    // Sample light source with multiple importance sampling
    Vector wi;
    float lightPdf, bsdfPdf;
    VisibilityTester visibility;
    Spectrum Li = light->Sample_L(p, rayEpsilon, lightSample,
                                  &wi, &lightPdf, &visibility);
    if (lightPdf > 0. && !Li.IsBlack()) {
        Spectrum f = bsdf->f(wo, wi);
        if (!f.IsBlack() && visibility.Unoccluded(scene, time)) {
            // Add light's contribution to reflected radiance
            Li *= visibility.Transmittance(scene, renderer, time, NULL, rng, arena);
            if (light->IsDeltaLight())
                Ld += f * Li * AbsDot(wi, n) / lightPdf;
            else {
                bsdfPdf = bsdf->Pdf(wo, wi);
                float weight = PowerHeuristic(1, lightPdf, 1, bsdfPdf);
                Ld += f * Li * AbsDot(wi, n) * weight / lightPdf;
            }
        }
    }

    // Sample BSDF with multiple importance sampling
    if (!light->IsDeltaLight()) {
        BxDFType flags = BxDFType(BSDF_ALL & ~BSDF_SPECULAR);
        Spectrum f = bsdf->Sample_f(wo, &wi, bsdfSample, &bsdfPdf, flags);
        if (!f.IsBlack() && bsdfPdf > 0.) {
            lightPdf = light->Pdf(p, wi);
            if (lightPdf > 0.) {
                // Add light contribution from BSDF sampling
                float weight = PowerHeuristic(1, bsdfPdf, 1, lightPdf);
                Intersection lightIsect;
                Spectrum Li(0.f);
                RayDifferential ray(p, wi, rayEpsilon, INFINITY, time);
                if (scene->Intersect(ray, &lightIsect)) {
                    if (lightIsect.primitive->GetAreaLight() == light)
                        Li = lightIsect.Le(-wi);
                }
                else
                    Li = light->Le(ray);
                if (!Li.IsBlack()) {
                    Li *= renderer->Transmittance(scene, ray, NULL, arena, rng);
                    Ld += f * Li * AbsDot(wi, n) * weight / bsdfPdf;
                }
            }
        }
    }
    return Ld;
}


Spectrum SpecularReflect(const RayDifferential &ray, BSDF *bsdf,
        RNG &rng, const Intersection &isect, const Renderer *renderer,
        const Scene *scene, const Sample *sample, MemoryArena &arena) {
    Vector wo = -ray.d, wi;
    const Point &p = bsdf->dgShading.p;
    const Normal &n = bsdf->dgShading.nn;
    Spectrum f = bsdf->Sample_f(wo, &wi, rng,
        BxDFType(BSDF_REFLECTION | BSDF_SPECULAR));
    Spectrum L = 0.f;
    if (!f.IsBlack() && AbsDot(wi, n) != 0.f) {
        // Compute ray differential _rd_ for specular reflection
        RayDifferential rd(p, wi, isect.RayEpsilon);
        rd.hasDifferentials = true;
        rd.depth = ray.depth + 1;
        rd.time = ray.time;
        rd.rxOrigin = p + isect.dg.dpdx;
        rd.ryOrigin = p + isect.dg.dpdy;

        // Compute differential reflected directions
        Normal dndx = bsdf->dgShading.dndu * bsdf->dgShading.dudx +
                      bsdf->dgShading.dndv * bsdf->dgShading.dvdx;
        Normal dndy = bsdf->dgShading.dndu * bsdf->dgShading.dudy +
                      bsdf->dgShading.dndv * bsdf->dgShading.dvdy;
        Vector dwodx = -ray.rxDirection - wo, dwody = -ray.ryDirection - wo;
        float dDNdx = Dot(dwodx, n) + Dot(wo, dndx);
        float dDNdy = Dot(dwody, n) + Dot(wo, dndy);
        rd.rxDirection = wi - dwodx + 2 * Vector(Dot(wo, n) * dndx +
                                          dDNdx * n);
        rd.ryDirection = wi - dwody + 2 * Vector(Dot(wo, n) * dndy +
                                          dDNdy * n);
        PBRT_STARTED_SPECULAR_REFLECTION_RAY(const_cast<RayDifferential *>(&rd));
        L = renderer->Li(scene, rd, sample, arena) * f * AbsDot(wi, n);
        PBRT_FINISHED_SPECULAR_REFLECTION_RAY(const_cast<RayDifferential *>(&rd));
    }
    return L;
}


Spectrum SpecularTransmit(const RayDifferential &ray, BSDF *bsdf,
        RNG &rng, const Intersection &isect, const Renderer *renderer,
        const Scene *scene, const Sample *sample, MemoryArena &arena) {
    Vector wo = -ray.d, wi;
    const Point &p = bsdf->dgShading.p;
    const Normal &n = bsdf->dgShading.nn;
    Spectrum f = bsdf->Sample_f(wo, &wi, rng,
        BxDFType(BSDF_TRANSMISSION | BSDF_SPECULAR));
    Spectrum L = 0.f;
    if (!f.IsBlack() && AbsDot(wi, n) != 0.f) {
        // Compute ray differential _rd_ for specular transmission
        RayDifferential rd(p, wi, isect.RayEpsilon);
        rd.hasDifferentials = true;
        rd.depth = ray.depth + 1;
        rd.time = ray.time;
        rd.rxOrigin = p + isect.dg.dpdx;
        rd.ryOrigin = p + isect.dg.dpdy;
        
        float eta = bsdf->eta;
        Vector w = -wo;
        if (Dot(wo, n) < 0) eta = 1.f / eta;
        
        Normal dndx = bsdf->dgShading.dndu * bsdf->dgShading.dudx + bsdf->dgShading.dndv * bsdf->dgShading.dvdx;
        Normal dndy = bsdf->dgShading.dndu * bsdf->dgShading.dudy + bsdf->dgShading.dndv * bsdf->dgShading.dvdy;
        
        Vector dwodx = -ray.rxDirection - wo, dwody = -ray.ryDirection - wo;
        float dDNdx = Dot(dwodx, n) + Dot(wo, dndx);
        float dDNdy = Dot(dwody, n) + Dot(wo, dndy);
        
        float mu = eta * Dot(w, n) - Dot(wi, n);
        float dmudx = (eta - (eta*eta*Dot(w,n))/Dot(wi, n)) * dDNdx;
        float dmudy = (eta - (eta*eta*Dot(w,n))/Dot(wi, n)) * dDNdy;
        
        rd.rxDirection = wi + eta * dwodx - Vector(mu * dndx + dmudx * n);
        rd.ryDirection = wi + eta * dwody - Vector(mu * dndy + dmudy * n);
        PBRT_STARTED_SPECULAR_REFRACTION_RAY(const_cast<RayDifferential *>(&rd));
        L = renderer->Li(scene, rd, sample, arena) * f * AbsDot(wi, n);
        PBRT_FINISHED_SPECULAR_REFRACTION_RAY(const_cast<RayDifferential *>(&rd));
    }
    return L;
}


Distribution1D *ComputeLightSamplingCDF(const Scene *scene) {
    u_int nLights = int(scene->lights.size());
    Assert(nLights > 0);
    vector<float>lightPower(nLights, 0.f);
    for (u_int i = 0; i < nLights; ++i)
        lightPower[i] = scene->lights[i]->Power(scene).y();
    return new Distribution1D(&lightPower[0], nLights);
}


