
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


// integrators/irradiancecache.cpp*
#include "stdafx.h"
#include "integrators/irradiancecache.h"
#include "camera.h"
#include "progressreporter.h"
#include "scene.h"
#include "montecarlo.h"
#include "film.h"
#include "samplers/halton.h"
#include "intersection.h"
#include "paramset.h"

// IrradianceCacheIntegrator Local Declarations
struct IrradiancePrimeTask : public Task {
    IrradiancePrimeTask(const Scene *sc, const Renderer *sr, const Camera *c, Sampler *samp,
                        Sample *s, IrradianceCacheIntegrator *ic, ProgressReporter &pr,
                        int tn, int nt) : progress(pr) {
        scene = sc;
        renderer = sr;
        camera = c;
        origSample = s;
        sampler = samp->GetSubSampler(tn, nt);
        irradianceCache = ic;
        taskNum = tn;
        numTasks = nt;
    }
    void Run();

    const Scene *scene;
    const Camera *camera;
    const Renderer *renderer;
    Sampler *sampler;
    Sample *origSample;
    IrradianceCacheIntegrator *irradianceCache;
    ProgressReporter &progress;
    int taskNum, numTasks;
};


struct IrradProcess {
    // IrradProcess Public Methods
    IrradProcess(const Point &P, const Normal &N, float mw, float cmsad) {
        p = P;
        n = N;
        minWeight = mw;
        cosMaxSampleAngleDifference = cmsad;
        nFound = 0;
        sumWt = 0.;
        E = 0.;
        wAvg = Vector(0,0,0);
    }
    bool operator()(const IrradianceSample *sample);
    bool Successful() {
        return sumWt >= minWeight;
    }
    Spectrum GetIrradiance() const { return E / sumWt; }
    Vector GetAverageDirection() const { return wAvg; }

    // IrradProcess Data
    Point p;
    Normal n;
    float minWeight, cosMaxSampleAngleDifference, sumWt;
    int nFound;
    Spectrum E;
    Vector wAvg;
};


struct IrradianceSample {
    // IrradianceSample Constructor
    IrradianceSample() { }
    IrradianceSample(const Spectrum &e, const Point &P, const Normal &N,
        const Vector &pd, float md) : E(e), n(N), p(P), wAvg(pd) {
        maxDist = md;
    }
    Spectrum E;
    Normal n;
    Point p;
    Vector wAvg;
    float maxDist;
};



// IrradianceCacheIntegrator Method Definitions
void IrradianceCacheIntegrator::RequestSamples(Sampler *sampler,
        Sample *sample, const Scene *scene) {
//    if (lightSampleOffsets != NULL) return;
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
}


void IrradianceCacheIntegrator::Preprocess(const Scene *scene,
        const Camera *camera, const Renderer *renderer) {
    BBox wb = scene->WorldBound();
    Vector delta = .01f * (wb.pMax - wb.pMin);
    wb.pMin -= delta;
    wb.pMax += delta;
    octree = new Octree<IrradianceSample *>(wb);
    // Prime irradiance cache
    minWeight *= 1.5f;
    int xstart, xend, ystart, yend;
    camera->film->GetSampleExtent(&xstart, &xend, &ystart, &yend);
    HaltonSampler sampler(xstart, xend, ystart, yend, 1,
                          camera->shutterOpen, camera->shutterClose);
    Sample *sample = new Sample(&sampler, this, NULL, scene);
    const int nTasks = 64;
    ProgressReporter progress(nTasks, "Priming irradiance cache");
    vector<Task *> tasks;
    for (int i = 0; i < nTasks; ++i)
        tasks.push_back(new IrradiancePrimeTask(scene, renderer, camera,
                                                &sampler, sample, this,
                                                progress, i, nTasks));
    EnqueueTasks(tasks);
    WaitForAllTasks();
    for (uint32_t i = 0; i < tasks.size(); ++i)
        delete tasks[i];
    progress.Done();
    delete sample;
    minWeight /= 1.5f;
}


IrradianceCacheIntegrator::~IrradianceCacheIntegrator() {
    delete octree;
    RWMutex::Destroy(mutex);
    delete[] lightSampleOffsets;
    delete[] bsdfSampleOffsets;
}


void IrradiancePrimeTask::Run() {
    if (!sampler) { progress.Update(); return; }
    MemoryArena arena;
    int sampleCount;
    RNG rng(29 * taskNum);
    int maxSamples = sampler->MaximumSampleCount();
    Sample *samples = origSample->Duplicate(maxSamples);
    while ((sampleCount = sampler->GetMoreSamples(samples, rng)) > 0) {
        for (int i = 0; i < sampleCount; ++i) {
            RayDifferential ray;
            camera->GenerateRayDifferential(samples[i], &ray);
            Intersection isect;
            if (scene->Intersect(ray, &isect))
                (void)irradianceCache->Li(scene, renderer, ray, isect, &samples[i], rng, arena);
        }
        arena.FreeAll();
    }
    delete[] samples;
    delete sampler;
    progress.Update();
}


Spectrum IrradianceCacheIntegrator::Li(const Scene *scene,
        const Renderer *renderer, const RayDifferential &ray, const Intersection &isect,
        const Sample *sample, RNG &rng, MemoryArena &arena) const {
    Spectrum L(0.);
    // Evaluate BSDF at hit point
    BSDF *bsdf = isect.GetBSDF(ray, arena);
    Vector wo = -ray.d;
    const Point &p = bsdf->dgShading.p;
    const Normal &n = bsdf->dgShading.nn;
    L += isect.Le(wo);
    // Compute direct lighting for irradiance cache
    L += UniformSampleAllLights(scene, renderer, arena, p, n, wo,
             isect.rayEpsilon, ray.time, bsdf, sample, rng,
             lightSampleOffsets, bsdfSampleOffsets);

    // Compute indirect lighting for irradiance cache
    if (ray.depth + 1 < maxSpecularDepth) {
        Vector wi;
        // Trace rays for specular reflection and refraction
        L += SpecularReflect(ray, bsdf, rng, isect, renderer, scene, sample,
                             arena);
        L += SpecularTransmit(ray, bsdf, rng, isect, renderer, scene, sample,
                              arena);
    }

    // Estimate indirect lighting with irradiance cache
    Normal ng = isect.dg.nn;
    ng = Faceforward(ng, wo);

    // Compute pixel spacing in world space at intersection point
    float pixelSpacing = sqrtf(Cross(isect.dg.dpdx, isect.dg.dpdy).Length());
    BxDFType flags = BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE | BSDF_GLOSSY);
    L += indirectLo(p, ng, pixelSpacing, wo, isect.rayEpsilon,
                    bsdf, flags, rng, scene, renderer, arena);
    flags = BxDFType(BSDF_TRANSMISSION | BSDF_DIFFUSE | BSDF_GLOSSY);
    L += indirectLo(p, -ng, pixelSpacing, wo, isect.rayEpsilon,
                    bsdf, flags, rng, scene, renderer, arena);
    return L;
}


Spectrum IrradianceCacheIntegrator::indirectLo(const Point &p,
        const Normal &ng, float pixelSpacing, const Vector &wo,
        float rayEpsilon, BSDF *bsdf, BxDFType flags, RNG &rng,
        const Scene *scene, const Renderer *renderer,
        MemoryArena &arena) const {
    if (bsdf->NumComponents(flags) == 0)
        return Spectrum(0.);
    Spectrum E;
    Vector wi;
    // Get irradiance _E_ and average incident direction _wi_ at point _p_
    if (!interpolateE(scene, p, ng, &E, &wi)) {
        // Compute irradiance at current point
        PBRT_IRRADIANCE_CACHE_STARTED_COMPUTING_IRRADIANCE(const_cast<Point *>(&p), const_cast<Normal *>(&ng));
        uint32_t scramble[2] = { rng.RandomUInt(), rng.RandomUInt() };
        float minHitDistance = INFINITY;
        Vector wAvg(0,0,0);
        Spectrum LiSum = 0.f;
        for (int i = 0; i < nSamples; ++i) {
            // Sample direction for irradiance estimate ray
            float u[2];
            Sample02(i, scramble, u);
            Vector w = CosineSampleHemisphere(u[0], u[1]);
            RayDifferential r(p, bsdf->LocalToWorld(w), rayEpsilon);
            r.d = Faceforward(r.d, ng);

            // Trace ray to sample radiance for irradiance estimate
            PBRT_IRRADIANCE_CACHE_STARTED_RAY(&r);
            Spectrum L = pathL(r, scene, renderer, rng, arena);
            LiSum += L;
            wAvg += r.d * L.y();
            minHitDistance = min(minHitDistance, r.maxt);
            PBRT_IRRADIANCE_CACHE_FINISHED_RAY(&r, r.maxt, &L);
        }
        E = (M_PI / float(nSamples)) * LiSum;
        PBRT_IRRADIANCE_CACHE_FINISHED_COMPUTING_IRRADIANCE(const_cast<Point *>(&p), const_cast<Normal *>(&ng));

        // Add computed irradiance value to cache

        // Compute irradiance sample's contribution extent and bounding box
        float maxDist = maxSamplePixelSpacing * pixelSpacing;
        float minDist = minSamplePixelSpacing * pixelSpacing;
        float contribExtent = Clamp(minHitDistance / 2.f, minDist, maxDist);
        BBox sampleExtent(p);
        sampleExtent.Expand(contribExtent);
        PBRT_IRRADIANCE_CACHE_ADDED_NEW_SAMPLE(const_cast<Point *>(&p), const_cast<Normal *>(&ng), contribExtent, &E, &wAvg, pixelSpacing);

        // Allocate _IrradianceSample_, get write lock, add to octree
        IrradianceSample *sample = new IrradianceSample(E, p, ng, wAvg,
                                                        contribExtent);
        RWMutexLock lock(*mutex, WRITE);
        octree->Add(sample, sampleExtent);
        wi = wAvg;
    }

    // Compute reflected radiance due to irradiance and BSDF
    if (wi.LengthSquared() == 0.f) return Spectrum(0.);
    return bsdf->f(wo, Normalize(wi), flags) * E;
}


bool IrradianceCacheIntegrator::interpolateE(const Scene *scene,
        const Point &p, const Normal &n, Spectrum *E,
        Vector *wi) const {
    if (!octree) return false;
    PBRT_IRRADIANCE_CACHE_STARTED_INTERPOLATION(const_cast<Point *>(&p), const_cast<Normal *>(&n));
    IrradProcess proc(p, n, minWeight, cosMaxSampleAngleDifference);
    RWMutexLock lock(*mutex, READ);
    octree->Lookup(p, proc);
    PBRT_IRRADIANCE_CACHE_FINISHED_INTERPOLATION(const_cast<Point *>(&p), const_cast<Normal *>(&n),
        proc.Successful() ? 1 : 0, proc.nFound);
    if (!proc.Successful()) return false;
    *E = proc.GetIrradiance();
    *wi = proc.GetAverageDirection();
    return true;
}


bool IrradProcess::operator()(const IrradianceSample *sample) {
    // Compute estimate error term and possibly use sample
    float perr = Distance(p, sample->p) / sample->maxDist;
    float nerr = sqrtf((1.f - Dot(n, sample->n)) /
                       (1.f - cosMaxSampleAngleDifference));
    float err = max(perr, nerr);
    PBRT_IRRADIANCE_CACHE_CHECKED_SAMPLE(const_cast<IrradianceSample *>(sample), perr, nerr);
    if (err < 1.) {
        ++nFound;
        float wt = 1.f - err;
        E += wt * sample->E;
        wAvg += wt * sample->wAvg;
        sumWt += wt;
    }
    return true;
}


Spectrum IrradianceCacheIntegrator::pathL(Ray &r, const Scene *scene,
        const Renderer *renderer, RNG &rng, MemoryArena &arena) const {
    Spectrum L(0.f);
    Spectrum pathThroughput = 1.;
    RayDifferential ray(r);
    bool specularBounce = false;
    for (int pathLength = 0; ; ++pathLength) {
        // Find next vertex of path
        Intersection isect;
        if (!scene->Intersect(ray, &isect))
            break;
        if (pathLength == 0)
            r.maxt = ray.maxt;
        pathThroughput *= renderer->Transmittance(scene, ray, NULL, rng, arena);
        // Possibly add emitted light at path vertex
        if (specularBounce)
            L += pathThroughput * isect.Le(-ray.d);
        // Evaluate BSDF at hit point
        BSDF *bsdf = isect.GetBSDF(ray, arena);
        // Sample illumination from lights to find path contribution
        const Point &p = bsdf->dgShading.p;
        const Normal &n = bsdf->dgShading.nn;
        Vector wo = -ray.d;
        L += pathThroughput *
            UniformSampleOneLight(scene, renderer, arena, p, n, wo, isect.rayEpsilon,
                                  ray.time, bsdf, NULL, rng);
        if (pathLength+1 == maxIndirectDepth) break;
        // Sample BSDF to get new path direction
        // Get random numbers for sampling new direction, \mono{bs1}, \mono{bs2}, and \mono{bcs}
        Vector wi;
        float pdf;
        BxDFType flags;
        Spectrum f = bsdf->Sample_f(wo, &wi, BSDFSample(rng),
            &pdf, BSDF_ALL, &flags);
        if (f.IsBlack() || pdf == 0.)
            break;
        specularBounce = (flags & BSDF_SPECULAR) != 0;
        pathThroughput *= f * AbsDot(wi, n) / pdf;
        ray = RayDifferential(p, wi, ray, isect.rayEpsilon);
        // Possibly terminate the path
        if (pathLength > 2) {
            float rrProb = min(1.f, pathThroughput.y());
            if (rng.RandomFloat() > rrProb)
                break;
            pathThroughput /= rrProb;
        }
    }
    return L;
}


IrradianceCacheIntegrator *CreateIrradianceCacheIntegrator(const ParamSet &params) {
    float minWeight = params.FindOneFloat("minweight", 0.5f);
    float minSpacing = params.FindOneFloat("minpixelspacing", 2.5f);
    float maxSpacing = params.FindOneFloat("maxpixelspacing", 15.f);
    float maxAngle = params.FindOneFloat("maxangledifference", 10.f);
    int maxSpecularDepth = params.FindOneInt("maxspeculardepth", 5);
    int maxIndirectDepth = params.FindOneInt("maxindirectdepth", 3);
    int nSamples = params.FindOneInt("nsamples", 4096);
    if (PbrtOptions.quickRender) nSamples = max(1, nSamples / 16);
    return new IrradianceCacheIntegrator(minWeight, minSpacing, maxSpacing, maxAngle,
        maxSpecularDepth, maxIndirectDepth, nSamples);
}


