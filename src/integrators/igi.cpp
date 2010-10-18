
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


// integrators/igi.cpp*
#include "stdafx.h"
#include "integrators/igi.h"
#include "scene.h"
#include "montecarlo.h"
#include "progressreporter.h"
#include "sampler.h"
#include "intersection.h"
#include "paramset.h"
#include "camera.h"

// IGIIntegrator Method Definitions
IGIIntegrator::~IGIIntegrator() {
    delete[] lightSampleOffsets;
    delete[] bsdfSampleOffsets;
}


void IGIIntegrator::RequestSamples(Sampler *sampler, Sample *sample,
                                   const Scene *scene) {
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
    vlSetOffset = sample->Add1D(1);

    if (sampler) nGatherSamples = sampler->RoundSize(nGatherSamples);
    gatherSampleOffset = BSDFSampleOffsets(nGatherSamples, sample);
}


void IGIIntegrator::Preprocess(const Scene *scene, const Camera *camera,
                               const Renderer *renderer) {
    if (scene->lights.size() == 0) return;
    MemoryArena arena;
    RNG rng;
    // Compute samples for emitted rays from lights
    vector<float> lightNum(nLightPaths * nLightSets);
    vector<float> lightSampPos(2 * nLightPaths * nLightSets, 0.f);
    vector<float> lightSampComp(nLightPaths * nLightSets, 0.f);
    vector<float> lightSampDir(2 * nLightPaths * nLightSets, 0.f);
    LDShuffleScrambled1D(nLightPaths, nLightSets, &lightNum[0], rng);
    LDShuffleScrambled2D(nLightPaths, nLightSets, &lightSampPos[0], rng);
    LDShuffleScrambled1D(nLightPaths, nLightSets, &lightSampComp[0], rng);
    LDShuffleScrambled2D(nLightPaths, nLightSets, &lightSampDir[0], rng);

    // Precompute information for light sampling densities
    Distribution1D *lightDistribution = ComputeLightSamplingCDF(scene);
    for (uint32_t s = 0; s < nLightSets; ++s) {
        for (uint32_t i = 0; i < nLightPaths; ++i) {
            // Follow path _i_ from light to create virtual lights
            int sampOffset = s*nLightPaths + i;

            // Choose light source to trace virtual light path from
            float lightPdf;
            int ln = lightDistribution->SampleDiscrete(lightNum[sampOffset],
                                                       &lightPdf);
            Light *light = scene->lights[ln];

            // Sample ray leaving light source for virtual light path
            RayDifferential ray;
            float pdf;
            LightSample ls(lightSampPos[2*sampOffset], lightSampPos[2*sampOffset+1],
                           lightSampComp[sampOffset]);
            Normal Nl;
            Spectrum alpha = light->Sample_L(scene, ls, lightSampDir[2*sampOffset],
                                             lightSampDir[2*sampOffset+1],
                                             camera->shutterOpen, &ray, &Nl, &pdf);
            if (pdf == 0.f || alpha.IsBlack()) continue;
            alpha /= pdf * lightPdf;
            Intersection isect;
            while (scene->Intersect(ray, &isect) && !alpha.IsBlack()) {
                // Create virtual light and sample new ray for path
                alpha *= renderer->Transmittance(scene, RayDifferential(ray), NULL,
                                                 rng, arena);
                Vector wo = -ray.d;
                BSDF *bsdf = isect.GetBSDF(ray, arena);

                // Create virtual light at ray intersection point
                Spectrum contrib = alpha * bsdf->rho(wo, rng) / M_PI;
                virtualLights[s].push_back(VirtualLight(isect.dg.p, isect.dg.nn, contrib,
                                                        isect.rayEpsilon));

                // Sample new ray direction and update weight for virtual light path
                Vector wi;
                float pdf;
                BSDFSample bsdfSample(rng);
                Spectrum fr = bsdf->Sample_f(wo, &wi, bsdfSample, &pdf);
                if (fr.IsBlack() || pdf == 0.f)
                    break;
                Spectrum contribScale = fr * AbsDot(wi, bsdf->dgShading.nn) / pdf;

                // Possibly terminate virtual light path with Russian roulette
                float rrProb = min(1.f, contribScale.y());
                if (rng.RandomFloat() > rrProb)
                    break;
                alpha *= contribScale / rrProb;
                ray = RayDifferential(isect.dg.p, wi, ray, isect.rayEpsilon);
            }
            arena.FreeAll();
        }
    }
    delete lightDistribution;
}


Spectrum IGIIntegrator::Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential &ray, const Intersection &isect,
        const Sample *sample, RNG &rng, MemoryArena &arena) const {
    Spectrum L(0.);
    Vector wo = -ray.d;
    // Compute emitted light if ray hit an area light source
    L += isect.Le(wo);

    // Evaluate BSDF at hit point
    BSDF *bsdf = isect.GetBSDF(ray, arena);
    const Point &p = bsdf->dgShading.p;
    const Normal &n = bsdf->dgShading.nn;
    L += UniformSampleAllLights(scene, renderer, arena, p, n,
                    wo, isect.rayEpsilon, ray.time, bsdf, sample, rng,
                    lightSampleOffsets, bsdfSampleOffsets);
    // Compute indirect illumination with virtual lights
    uint32_t lSet = min(uint32_t(sample->oneD[vlSetOffset][0] * nLightSets),
                        nLightSets-1);
    for (uint32_t i = 0; i < virtualLights[lSet].size(); ++i) {
        const VirtualLight &vl = virtualLights[lSet][i];
        // Compute virtual light's tentative contribution _Llight_
        float d2 = DistanceSquared(p, vl.p);
        Vector wi = Normalize(vl.p - p);
        float G = AbsDot(wi, n) * AbsDot(wi, vl.n) / d2;
        G = min(G, gLimit);
        Spectrum f = bsdf->f(wo, wi);
        if (G == 0.f || f.IsBlack()) continue;
        Spectrum Llight = f * G * vl.pathContrib / nLightPaths;
        RayDifferential connectRay(p, wi, ray, isect.rayEpsilon,
                                   sqrtf(d2) * (1.f - vl.rayEpsilon));
        Llight *= renderer->Transmittance(scene, connectRay, NULL, rng, arena);

        // Possibly skip virtual light shadow ray with Russian roulette
        if (Llight.y() < rrThreshold) {
            float continueProbability = .1f;
            if (rng.RandomFloat() > continueProbability)
                continue;
            Llight /= continueProbability;
        }

        // Add contribution from _VirtualLight_ _vl_
        if (!scene->IntersectP(connectRay))
            L += Llight;
    }
    if (ray.depth < maxSpecularDepth) {
        // Do bias compensation for bounding geometry term
        int nSamples = (ray.depth == 0) ? nGatherSamples : 1;
        for (int i = 0; i < nSamples; ++i) {
            Vector wi;
            float pdf;
            BSDFSample bsdfSample = (ray.depth == 0) ?
                BSDFSample(sample, gatherSampleOffset, i) : BSDFSample(rng);
            Spectrum f = bsdf->Sample_f(wo, &wi, bsdfSample,
                                        &pdf, BxDFType(BSDF_ALL & ~BSDF_SPECULAR));
            if (!f.IsBlack() && pdf > 0.f) {
                // Trace ray for bias compensation gather sample
                float maxDist = sqrtf(AbsDot(wi, n) / gLimit);
                RayDifferential gatherRay(p, wi, ray, isect.rayEpsilon, maxDist);
                Intersection gatherIsect;
                Spectrum Li = renderer->Li(scene, gatherRay, sample, rng, arena,
                                           &gatherIsect);
                if (Li.IsBlack()) continue;

                // Add bias compensation ray contribution to radiance sum
                float Ggather = AbsDot(wi, n) * AbsDot(-wi, gatherIsect.dg.nn) /
                    DistanceSquared(p, gatherIsect.dg.p);
                if (Ggather - gLimit > 0.f && !isinf(Ggather)) {
                    float gs = (Ggather - gLimit) / Ggather;
                    L += f * Li * (AbsDot(wi, n) * gs / (nSamples * pdf));
                }
            }
        }
    }
    if (ray.depth + 1 < maxSpecularDepth) {
        Vector wi;
        // Trace rays for specular reflection and refraction
        L += SpecularReflect(ray, bsdf, rng, isect, renderer, scene, sample,
                             arena);
        L += SpecularTransmit(ray, bsdf, rng, isect, renderer, scene, sample,
                              arena);
    }
    return L;
}


IGIIntegrator *CreateIGISurfaceIntegrator(const ParamSet &params) {
    int nLightPaths = params.FindOneInt("nlights", 64);
    if (PbrtOptions.quickRender) nLightPaths = max(1, nLightPaths / 4);
    int nLightSets = params.FindOneInt("nsets", 4);
    float rrThresh = params.FindOneFloat("rrthreshold", .0001f);
    int maxDepth = params.FindOneInt("maxdepth", 5);
    float glimit = params.FindOneFloat("glimit", 10.f);
    int gatherSamples = params.FindOneInt("gathersamples", 16);
    return new IGIIntegrator(nLightPaths, nLightSets, rrThresh,
                             maxDepth, glimit, gatherSamples);
}


