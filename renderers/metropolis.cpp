
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


// renderers/metropolis.cpp*
#include "renderers/metropolis.h"
#include "renderers/samplerrenderer.h"
#include "scene.h"
#include "imageio.h"
#include "spectrum.h"
#include "camera.h"
#include "film.h"
#include "rng.h"
#include "progressreporter.h"
#include "paramset.h"
#include "parallel.h"
#include "probes.h"
#include "intersection.h"
#include "montecarlo.h"
#include "samplers/lowdiscrepancy.h"
#include "integrators/directlighting.h"

// Metropolis Local Declarations
struct PathSample {
    BSDFSample bsdfSample;
    float rrSample;
};


struct LightingSample {
    BSDFSample bsdfSample;
    float lightNum;
    LightSample lightSample;
};


struct MLTSample {
    MLTSample(int maxLength) {
        cameraPathSamples.resize(maxLength);
        lightPathSamples.resize(maxLength);
        lightingSamples.resize(maxLength);
    }
    CameraSample cameraSample;
    float lightNumSample, lightRaySamples[5];
    vector<PathSample> cameraPathSamples, lightPathSamples;
    vector<LightingSample> lightingSamples;
};


static void LargeStep(RNG &rng, MLTSample *sample, int maxDepth,
        float x, float y, float t0, float t1, bool bidirectional) {
    // Do large step mutation of _cameraSample_
    sample->cameraSample.imageX = x;
    sample->cameraSample.imageY = y;
    sample->cameraSample.time = Lerp(rng.RandomFloat(), t0, t1);
    sample->cameraSample.lensU = rng.RandomFloat();
    sample->cameraSample.lensV = rng.RandomFloat();
    for (int i = 0; i < maxDepth; ++i) {
        // Apply large step to $i$th camera _PathSample_
        PathSample &eps = sample->cameraPathSamples[i];
        eps.bsdfSample.uComponent = rng.RandomFloat();
        eps.bsdfSample.uDir[0] = rng.RandomFloat();
        eps.bsdfSample.uDir[1] = rng.RandomFloat();
        eps.rrSample = rng.RandomFloat();

        // Apply large step to $i$th _LightingSample_
        LightingSample &ls = sample->lightingSamples[i];
        ls.bsdfSample.uComponent = rng.RandomFloat();
        ls.bsdfSample.uDir[0] = rng.RandomFloat();
        ls.bsdfSample.uDir[1] = rng.RandomFloat();
        ls.lightNum = rng.RandomFloat();
        ls.lightSample.uComponent = rng.RandomFloat();
        ls.lightSample.uPos[0] = rng.RandomFloat();
        ls.lightSample.uPos[1] = rng.RandomFloat();
    }
    if (bidirectional) {
        // Mutate bidirectional light ray and light path samples
        sample->lightNumSample = rng.RandomFloat();
        for (int i = 0; i < 5; ++i)
            sample->lightRaySamples[i] = rng.RandomFloat();
        for (int i = 0; i < maxDepth; ++i) {
            // Apply large step to $i$th light _PathSample_
            PathSample &lps = sample->lightPathSamples[i];
            lps.bsdfSample.uComponent = rng.RandomFloat();
            lps.bsdfSample.uDir[0] = rng.RandomFloat();
            lps.bsdfSample.uDir[1] = rng.RandomFloat();
            lps.rrSample = rng.RandomFloat();
        }
    }
}


static inline void mutate(RNG &rng, float *v, float a = 0.f,
                          float b = 1.f);
static inline void mutate(RNG &rng, float *v, float min, float max) {
    if (min == max) { *v = min; return; }
    Assert(min < max);
    float s1 = 1.f / 1024.f, s2 = 1.f / 64.f;
    static const float negLogRatio = -logf(s2/s1);
    float delta = (max - min) * s2 * expf(negLogRatio * rng.RandomFloat());
    if (rng.RandomFloat() < 0.5f) {
        *v += delta;
        if (*v > max) *v = min + (*v - max);
    }
    else {
        *v -= delta;
        if (*v < min) *v = max - (min - *v);
    }
    Assert(*v >= min && *v <= max);
}


static void SmallStep(RNG &rng, MLTSample *sample, int maxDepth,
        int x0, int x1, int y0, int y1, float t0, float t1,
        bool bidirectional) {
    mutate(rng, &sample->cameraSample.imageX, x0, x1);
    mutate(rng, &sample->cameraSample.imageY, y0, y1);
    mutate(rng, &sample->cameraSample.time, t0, t1);
    mutate(rng, &sample->cameraSample.lensU);
    mutate(rng, &sample->cameraSample.lensV);
    // Apply small step mutation to camera, lighting, and light samples
    for (int i = 0; i < maxDepth; ++i) {
        // Apply small step to $i$th camera _PathSample_
        PathSample &eps = sample->cameraPathSamples[i];
        mutate(rng, &eps.bsdfSample.uComponent);
        mutate(rng, &eps.bsdfSample.uDir[0]);
        mutate(rng, &eps.bsdfSample.uDir[1]);
        mutate(rng, &eps.rrSample);

        // Apply small step to $i$th _LightingSample_
        LightingSample &ls = sample->lightingSamples[i];
        mutate(rng, &ls.bsdfSample.uComponent);
        mutate(rng, &ls.bsdfSample.uDir[0]);
        mutate(rng, &ls.bsdfSample.uDir[1]);
        mutate(rng, &ls.lightNum);
        mutate(rng, &ls.lightSample.uComponent);
        mutate(rng, &ls.lightSample.uPos[0]);
        mutate(rng, &ls.lightSample.uPos[1]);
    }
    
    if (bidirectional) {
        mutate(rng, &sample->lightNumSample);
        for (int i = 0; i < 5; ++i)
            mutate(rng, &sample->lightRaySamples[i]);
        for (int i = 0; i < maxDepth; ++i) {
            // Apply small step to $i$th light _PathSample_
            PathSample &lps = sample->lightPathSamples[i];
            mutate(rng, &lps.bsdfSample.uComponent);
            mutate(rng, &lps.bsdfSample.uDir[0]);
            mutate(rng, &lps.bsdfSample.uDir[1]);
            mutate(rng, &lps.rrSample);
        }
    }
}


struct PathVertex {
    Intersection isect;
    Vector wPrev, wNext;
    BSDF *bsdf;
    bool specularBounce;
    Spectrum alpha;
};


static uint32_t GeneratePath(const RayDifferential &r, const Spectrum &alpha,
    const Scene *scene, MemoryArena &arena, const vector<PathSample> &samples,
    PathVertex *path, RayDifferential *escapedRay,
    Spectrum *escapedAlpha);
inline float I(const Spectrum &L);
class MLTTask : public Task {
public:
    MLTTask(ProgressReporter &prog, uint32_t pfreq, uint32_t taskNum,
        float dx, float dy, int xx0, int xx1, int yy0, int yy1, float tt0, float tt1,
        float bb, const MLTSample &is, bool id,
        const Scene *sc, const Camera *c, MetropolisRenderer *renderer,
        Mutex *filmMutex, Distribution1D *lightDistribution);
    void Run();

private:
    ProgressReporter &progress;
    uint32_t progressUpdateFrequency, taskNum;
    float dx, dy;
    int currentPixelSample;
    int x0, x1, y0, y1;
    float t0, t1;
    float b;
    bool directIsSeparate;
    const MLTSample &initialSample;
    const Scene *scene;
    const Camera *camera;
    MetropolisRenderer *renderer;
    Mutex *filmMutex;
    Distribution1D *lightDistribution;
};



// Metropolis Method Definitions
static uint32_t GeneratePath(const RayDifferential &r,
        const Spectrum &a, const Scene *scene, MemoryArena &arena,
        const vector<PathSample> &samples, PathVertex *path,
        RayDifferential *escapedRay, Spectrum *escapedAlpha) {
    RayDifferential ray = r;
    Spectrum alpha = a;
    if (escapedAlpha) *escapedAlpha = 0.f;
    uint32_t length = 0;
    for (; length < samples.size(); ++length) {
        // Try to generate next vertex of ray path
        if (!scene->Intersect(ray, &path[length].isect)) {
            // Handle ray that leaves the scene during path generation
            if (escapedAlpha) *escapedAlpha = alpha;
            if (escapedRay)   *escapedRay = ray;
            break;
        }

        // Record information for current path vertex
        PathVertex &v = path[length];
        v.alpha = alpha;
        BSDF *bsdf = v.isect.GetBSDF(ray, arena);
        v.bsdf = bsdf;
        v.wPrev = -ray.d;

        // Sample direction for outgoing Metropolis path direction
        Vector wo = -ray.d;
        float pdf;
        BxDFType flags;
        Spectrum f = bsdf->Sample_f(wo, &v.wNext, samples[length].bsdfSample,
                                    &pdf, BSDF_ALL, &flags);
        if (f.IsBlack() || pdf == 0.f) break;
        v.specularBounce = (flags & BSDF_SPECULAR) != 0;

        // Terminate path with RR or prepare for finding next vertex
        const Point &p = bsdf->dgShading.p;
        const Normal &n = bsdf->dgShading.nn;
        Spectrum pathScale = f * AbsDot(v.wNext, n) / pdf;
        float rrSurviveProb = min(1.f, pathScale.y());
        if (samples[length].rrSample > rrSurviveProb)
            return length+1;
        alpha *= pathScale / rrSurviveProb;
        //alpha *= renderer->Transmittance(scene, ray, NULL, rng, arena);
        ray = RayDifferential(p, v.wNext, ray, v.isect.rayEpsilon);
    }
    return length;
}


Spectrum MetropolisRenderer::PathL(const MLTSample &sample,
        const Scene *scene, MemoryArena &arena, const Camera *camera,
        const Distribution1D *lightDistribution,
        PathVertex *cameraPath, PathVertex *lightPath,
        RNG &rng) const {
    // Generate camera path from camera path samples
    RayDifferential cameraRay;
    float cameraWt = camera->GenerateRayDifferential(sample.cameraSample,
                                                     &cameraRay);
    cameraRay.ScaleDifferentials(1.f / sqrtf(nPixelSamples));
    RayDifferential escapedRay;
    Spectrum escapedAlpha;
    uint32_t cameraLength = GeneratePath(cameraRay, cameraWt, scene, arena,
        sample.cameraPathSamples, cameraPath, &escapedRay, &escapedAlpha);
    if (!bidirectional) {
        // Compute radiance along paths using path tracing
        return Lpath(scene, cameraPath, cameraLength, arena, rng,
                     sample.lightingSamples, sample.cameraSample.time,
                     lightDistribution, escapedRay, escapedAlpha);
    }
    else {
        // Sample light ray and apply bidirectional path tracing

        // Choose light and sample ray to start light path
        float lightPdf, lightRayPdf;
        uint32_t lightNum = lightDistribution->SampleDiscrete(sample.lightNumSample,
                                                              &lightPdf);
        const Light *light = scene->lights[lightNum];
        Ray lightRay;
        Normal Nl;
        LightSample lrs(sample.lightRaySamples[0], sample.lightRaySamples[1],
                        sample.lightRaySamples[2]);
        Spectrum lightWt = light->Sample_L(scene, lrs, sample.lightRaySamples[3],
            sample.lightRaySamples[4], sample.cameraSample.time,  &lightRay,
            &Nl, &lightRayPdf);
        if (lightWt.IsBlack() || lightRayPdf == 0.f) {
            // Compute radiance along paths using path tracing
            return Lpath(scene, cameraPath, cameraLength, arena, rng,
                         sample.lightingSamples, sample.cameraSample.time,
                         lightDistribution, escapedRay, escapedAlpha);
        }
        else {
            // Compute radiance along paths using bidirectional path tracing
            lightWt *= AbsDot(Normalize(Nl), lightRay.d) / (lightPdf * lightRayPdf);
            uint32_t lightLength = GeneratePath(RayDifferential(lightRay), lightWt,
                scene, arena, sample.lightPathSamples, lightPath, NULL, NULL);
            return Lbidir(scene, cameraPath, cameraLength, lightPath, lightLength,
                          arena, rng, sample.lightingSamples, sample.cameraSample.time,
                          lightDistribution, escapedRay, escapedAlpha);
        }
    }
}


Spectrum MetropolisRenderer::Lpath(const Scene *scene,
        const PathVertex *cameraPath, int cameraPathLength,
        MemoryArena &arena, RNG &rng, const vector<LightingSample> &samples,
        float time, const Distribution1D *lightDistribution,
        const RayDifferential &eRay, const Spectrum &eAlpha) const {
    Spectrum L = 0.;
    bool previousSpecular = true, allSpecular = true;
    for (int i = 0; i < cameraPathLength; ++i) {
        const PathVertex &vc = cameraPath[i];
        const Point &pc = vc.bsdf->dgShading.p;
        const Normal &nc = vc.bsdf->dgShading.nn;
        // Add emitted light from vertex if appropriate
        bool directIsSeparate = (directLighting != NULL);
        if (previousSpecular && (!directIsSeparate || !allSpecular))
            L += vc.alpha * vc.isect.Le(vc.wPrev);

        // Compute direct illumination for Metropolis path vertex
        Spectrum Ld(0.f);
        if (!directIsSeparate || !allSpecular) {
            // Choose light and call _EstimateDirect()_ for Metropolis vertex
            const LightingSample &ls = samples[i];
            float lightPdf;
            uint32_t lightNum = lightDistribution->SampleDiscrete(ls.lightNum,
                                                                  &lightPdf);
            const Light *light = scene->lights[lightNum];
            Ld = vc.alpha *
                 EstimateDirect(scene, this, arena, light, pc, nc, vc.wPrev,
                                vc.isect.rayEpsilon, time, vc.bsdf, rng,
                                ls.lightSample, ls.bsdfSample) / lightPdf;
        }
        previousSpecular = vc.specularBounce;
        allSpecular &= previousSpecular;
        L += Ld;
    }
    // Add contribution of escaped ray, if any
    if (!eAlpha.IsBlack() && previousSpecular &&
        (directLighting == NULL || !allSpecular))
        for (uint32_t i = 0; i < scene->lights.size(); ++i)
           L += eAlpha * scene->lights[i]->Le(eRay);
    return L;
}


Spectrum MetropolisRenderer::Lbidir(const Scene *scene,
        const PathVertex *cameraPath, int cameraPathLength,
        const PathVertex *lightPath, int lightPathLength,
        MemoryArena &arena, RNG &rng, const vector<LightingSample> &samples,
        float time, const Distribution1D *lightDistribution,
        const RayDifferential &eRay, const Spectrum &eAlpha) const {
    Spectrum L = 0.;
    bool previousSpecular = true, allSpecular = true;
    int nCameraSpecular = 0;
    for (int i = 0; i < cameraPathLength; ++i) {
        // Compute reflected light at camera path vertex
        const PathVertex &vc = cameraPath[i];
        const Point &pc = vc.bsdf->dgShading.p;
        const Normal &nc = vc.bsdf->dgShading.nn;

        // Add emitted light from vertex if appropriate
        bool directIsSeparate = (directLighting != NULL);
        if (previousSpecular && (!directIsSeparate || !allSpecular))
            L += vc.alpha * vc.isect.Le(vc.wPrev);

        // Compute direct illumination for Metropolis path vertex
        Spectrum Ld(0.f);
        if (!directIsSeparate || !allSpecular) {
            // Choose light and call _EstimateDirect()_ for Metropolis vertex
            const LightingSample &ls = samples[i];
            float lightPdf;
            uint32_t lightNum = lightDistribution->SampleDiscrete(ls.lightNum,
                                                                  &lightPdf);
            const Light *light = scene->lights[lightNum];
            Ld = vc.alpha *
                 EstimateDirect(scene, this, arena, light, pc, nc, vc.wPrev,
                                vc.isect.rayEpsilon, time, vc.bsdf, rng,
                                ls.lightSample, ls.bsdfSample) / lightPdf;
        }
        previousSpecular = vc.specularBounce;
        allSpecular &= previousSpecular;
        L += Ld / (i + 1);
        if (vc.specularBounce) ++nCameraSpecular;
        else {
            // Loop over light path vertices and connect to camera vertex
            int nLightSpecular = 0;
            for (int j = 0; j < lightPathLength; ++j) {
                const PathVertex &vl = lightPath[j];
                const Point &pl = vl.bsdf->dgShading.p;
                const Normal &nl = vl.bsdf->dgShading.nn;
                if (vl.specularBounce) ++nLightSpecular;
                else {
                    // Compute contribution between camera and light vertices
                    Vector w = Normalize(pl - pc);
                    Spectrum fc = vc.bsdf->f(vc.wPrev, w);
                    Spectrum fl = vl.bsdf->f(-w, vl.wPrev);
                    if (fc.IsBlack() || fl.IsBlack()) continue;
                    float G = AbsDot(nc, w) * AbsDot(nl, w) / DistanceSquared(pl, pc);
                    Ray r(pc, pl - pc, 1e-3f, .999f, time);
                    if (!scene->IntersectP(r)) {
                        // Compute weight for bidirectional path, _pathWt_
                        float pathWt = i + j + 2 - nCameraSpecular - nLightSpecular;
                        L += (vc.alpha * fc * G * fl * vl.alpha) / pathWt;
                    }
                }
            }
        }
    }
    // Add contribution of escaped ray, if any
    if (!eAlpha.IsBlack() && previousSpecular &&
        (directLighting == NULL || !allSpecular))
        for (uint32_t i = 0; i < scene->lights.size(); ++i)
           L += eAlpha * scene->lights[i]->Le(eRay);
    return L;
}


MetropolisRenderer::MetropolisRenderer(int perPixelSamples,
        int nboot, int dps, float lsp, bool dds, int mr, int md,
        Camera *c, bool db) {
    camera = c;

    nPixelSamples = perPixelSamples;
    float largeStepProbability = lsp;
    largeStepsPerPixel = max(1u, RoundUpPow2(largeStepProbability * nPixelSamples));
    if (largeStepsPerPixel >= nPixelSamples) largeStepsPerPixel /= 2;
    Assert(largeStepsPerPixel >= 1 && largeStepsPerPixel < nPixelSamples);
    Assert(IsPowerOf2(largeStepsPerPixel));
    if ((nPixelSamples % largeStepsPerPixel) != 0) {
        int origPixelSamples = nPixelSamples;
        nPixelSamples += largeStepsPerPixel - (nPixelSamples % largeStepsPerPixel);
        Warning("Rounding up to %d Metropolis samples per pixel (from %d)",
                nPixelSamples, origPixelSamples);
    }

    nBootstrap = nboot;
    nDirectPixelSamples = dps;

    maxDepth = md;
    maxConsecutiveRejects = mr;
    nTasksFinished  = 0;
    directLighting = dds ? new DirectLightingIntegrator(SAMPLE_ALL_UNIFORM, maxDepth) : NULL;
    bidirectional = db;
}


MetropolisRenderer::~MetropolisRenderer() {
    delete camera;
    delete directLighting;
}


MetropolisRenderer *CreateMetropolisRenderer(const ParamSet &params,
        Camera *camera) {
    float largeStepProbability = params.FindOneFloat("largestepprobability", .25f);
    int perPixelSamples = params.FindOneInt("samplesperpixel", 100);
    int nBootstrap = params.FindOneInt("bootstrapsamples", 100000);
    int nDirectPixelSamples = params.FindOneInt("directsamples", 4);
    bool doDirectSeparately = params.FindOneBool("dodirectseparately", true);
    int mr = params.FindOneInt("maxconsecutiverejects", 512);
    int md = params.FindOneInt("maxdepth", 7);
    bool doBidirectional = params.FindOneBool("bidirectional", true);

    if (PbrtOptions.quickRender) {
        perPixelSamples = max(1, perPixelSamples / 4);
        nBootstrap = max(1, nBootstrap / 4);
        nDirectPixelSamples = max(1, nDirectPixelSamples / 4);
    }

    return new MetropolisRenderer(perPixelSamples, nBootstrap,
        nDirectPixelSamples, largeStepProbability, doDirectSeparately,
        mr, md, camera, doBidirectional);
}


void MetropolisRenderer::Render(const Scene *scene) {
    if (scene->lights.size() > 0) {
        int x0, x1, y0, y1;
        camera->film->GetPixelExtent(&x0, &x1, &y0, &y1);
        float t0 = camera->shutterOpen, t1 = camera->shutterClose;
        Distribution1D *lightDistribution = ComputeLightSamplingCDF(scene);

        if (directLighting != NULL) {
            // Compute direct lighting before Metropolis light transport
            LDSampler sampler(x0, x1, y0, y1, nDirectPixelSamples, t0, t1);
            Sample *sample = new Sample(&sampler, directLighting, NULL, scene);
            vector<Task *> directTasks;
            int nDirectTasks = max(32 * NumSystemCores(),
                             (camera->film->xResolution * camera->film->yResolution) / (16*16));
            nDirectTasks = RoundUpPow2(nDirectTasks);
            ProgressReporter directProgress(nDirectTasks, "Direct Lighting");
            for (int i = 0; i < nDirectTasks; ++i)
                directTasks.push_back(new SamplerRendererTask(scene, this, camera, directProgress,
                    &sampler, sample, i, nDirectTasks));
            std::reverse(directTasks.begin(), directTasks.end());
            EnqueueTasks(directTasks);
            WaitForAllTasks();
            for (uint32_t i = 0; i < directTasks.size(); ++i)
                delete directTasks[i];
            delete sample;
            directProgress.Done();
        }
        // Take initial set of samples to compute $b$
        RNG rng(0);
        MemoryArena arena;
        vector<float> bootstrapSamples;
        vector<PathVertex> cameraPath(maxDepth, PathVertex());
        vector<PathVertex> lightPath(maxDepth, PathVertex());
        float sumContrib = 0.f;
        bootstrapSamples.reserve(nBootstrap);
        MLTSample sample(maxDepth);
        for (uint32_t i = 0; i < nBootstrap; ++i) {
            // Generate random sample and path radiance for MLT bootstrapping
            float x = Lerp(rng.RandomFloat(), x0, x1);
            float y = Lerp(rng.RandomFloat(), y0, y1);
            LargeStep(rng, &sample, maxDepth, x, y, t0, t1, bidirectional);
            Spectrum L = PathL(sample, scene, arena, camera, lightDistribution,
                               &cameraPath[0], &lightPath[0], rng);

            // Compute contribution for random sample for MLT bootstrapping
            float I = ::I(L);
            sumContrib += I;
            bootstrapSamples.push_back(I);
            arena.FreeAll();
        }
        float b = sumContrib / nBootstrap;

        // Select initial sample from bootstrap samples
        float contribOffset = rng.RandomFloat() * sumContrib;
        rng.Seed(0);
        sumContrib = 0.f;
        MLTSample initialSample(maxDepth);
        for (uint32_t i = 0; i < nBootstrap; ++i) {
            float x = Lerp(rng.RandomFloat(), x0, x1);
            float y = Lerp(rng.RandomFloat(), y0, y1);
            LargeStep(rng, &initialSample, maxDepth, x, y, t0, t1,
                      bidirectional);
            sumContrib += bootstrapSamples[i];
            if (contribOffset < sumContrib)
                break;
        }

        // Launch tasks to generate Metropolis samples
        uint32_t nTasks = largeStepsPerPixel;
        uint32_t largeStepRate = nPixelSamples / largeStepsPerPixel;
        fprintf(stderr, "n tasks %d, large step rate %d - %d\n", nTasks,
                                 largeStepRate, nTasks * largeStepRate);
        ProgressReporter progress(nTasks * largeStepRate, "Metropolis");
        vector<Task *> tasks;
        Mutex *filmMutex = Mutex::Create();
        Assert(IsPowerOf2(nTasks));
        uint32_t scramble[2] = { rng.RandomUInt(), rng.RandomUInt() };
        uint32_t pfreq = (x1-x0) * (y1-y0);
        for (uint32_t i = 0; i < nTasks; ++i) {
            float d[2];
            Sample02(i, scramble, d);
            tasks.push_back(new MLTTask(progress, pfreq, i,
                d[0], d[1], x0, x1, y0, y1, t0, t1, b, initialSample,
                directLighting != NULL, scene, camera, this,
                filmMutex, lightDistribution));
        }
        EnqueueTasks(tasks);
        WaitForAllTasks();
        for (uint32_t i = 0; i < tasks.size(); ++i)
            delete tasks[i];
        progress.Done();
        Mutex::Destroy(filmMutex);
    }
    camera->film->WriteImage();
}


inline float I(const Spectrum &L) {
    return L.y();
}


MLTTask::MLTTask(ProgressReporter &prog, uint32_t pfreq, uint32_t tn,
        float ddx, float ddy, int xx0, int xx1, int yy0, int yy1, float tt0, float tt1,
        float bb, const MLTSample &is, bool id, const Scene *sc, const Camera *c,
        MetropolisRenderer *ren, Mutex *fm, Distribution1D *ld)
    : progress(prog), initialSample(is) {
    progressUpdateFrequency = pfreq;
    taskNum = tn;
    dx = ddx;
    dy = ddy;
    x0 = xx0;
    x1 = xx1;
    y0 = yy0;
    y1 = yy1;
    t0 = tt0;
    t1 = tt1;
    currentPixelSample = 0;
    b = bb;
    directIsSeparate = id;
    scene = sc;
    camera = c;
    renderer = ren;
    filmMutex = fm;
    lightDistribution = ld;
}


void MLTTask::Run() {
    PBRT_MLT_STARTED_MLT_TASK(this);
    // Declare basic _MLTTask_ variables and prepare for sampling
    uint32_t nPixels = (x1-x0) * (y1-y0);
    uint32_t nPixelSamples = renderer->nPixelSamples;
    uint32_t largeStepRate = nPixelSamples / renderer->largeStepsPerPixel;
    Assert(largeStepRate > 1);
    uint64_t nTaskSamples = uint64_t(nPixels) * uint64_t(largeStepRate);
    RNG rng(taskNum);
    uint32_t pixelNumOffset = 0;
    
    vector<PathVertex> cameraPath(renderer->maxDepth, PathVertex());
    vector<PathVertex> lightPath(renderer->maxDepth, PathVertex());
    
    uint32_t consecutiveRejects = 0;
    float largeStepProbability = float(renderer->largeStepsPerPixel) /
                                 float(nPixelSamples);
    uint32_t progressCounter = progressUpdateFrequency;
    

    // Declare variables to store current and proposed MLT samples
    MemoryArena arena;
    vector<MLTSample> samples(2, MLTSample(renderer->maxDepth));
    Spectrum L[2];
    float I[2];
    uint32_t current = 0, proposed = 1;

    // Compute _L[current]_ for initial sample
    samples[current] = initialSample;
    L[current] = renderer->PathL(initialSample, scene, arena, camera,
                     lightDistribution, &cameraPath[0], &lightPath[0], rng);
    I[current] = ::I(L[current]);
    arena.FreeAll();

    // Compute randomly permuted table of pixel indices for large steps
    vector<int> largeStepPixelNum;
    for (uint32_t i = 0; i < nPixels; ++i) largeStepPixelNum.push_back(i);
    Shuffle(&largeStepPixelNum[0], nPixels, 1, rng);
    for (uint64_t s = 0; s < nTaskSamples; ++s) {
        // Compute proposed mutation to current sample
        samples[proposed] = samples[current];
        bool largeStep = ((s % largeStepRate) == 0);
        if (largeStep) {
            int x = x0 + largeStepPixelNum[pixelNumOffset] % (x1 - x0);
            int y = y0 + largeStepPixelNum[pixelNumOffset] / (x1 - x0);
            LargeStep(rng, &samples[proposed], renderer->maxDepth,
                      x + dx, y + dy, t0, t1, renderer->bidirectional);
            ++pixelNumOffset;
        }
        else
            SmallStep(rng, &samples[proposed], renderer->maxDepth,
                      x0, x1, y0, y1, t0, t1, renderer->bidirectional);

        // Compute contribution of proposed sample
        L[proposed] = renderer->PathL(samples[proposed], scene, arena, camera,
                         lightDistribution, &cameraPath[0], &lightPath[0], rng);
        I[proposed] = ::I(L[proposed]);
        arena.FreeAll();

        // Compute acceptance probability for proposed sample
        float a = min(1.f, I[proposed] / I[current]);
        float currentWeight = (1.f - a) /
                              (I[current] / b + largeStepProbability) /
                              float(nPixelSamples);
        float proposedWeight = (a + (largeStep ? 1.f : 0.f)) /
                               (I[proposed] / b + largeStepProbability) /
                               float(nPixelSamples);

        // Splat current and proposed samples to _Film_
        if (currentWeight > 0.f && I[current] > 0.f)
            camera->film->Splat(samples[current].cameraSample,
                                L[current] * currentWeight);
        if (proposedWeight > 0.f && I[proposed] > 0.f)
            camera->film->Splat(samples[proposed].cameraSample,
                                L[proposed] * proposedWeight);

        // Randomly accept proposed path mutation (or not)
        if (consecutiveRejects >= renderer->maxConsecutiveRejects ||
            rng.RandomFloat() < a) {
            PBRT_MLT_ACCEPTED_MUTATION(a, &samples[current], &samples[proposed]);
            current ^= 1;
            proposed ^= 1;
            consecutiveRejects = 0;
        }
        else
        {
            PBRT_MLT_REJECTED_MUTATION(a, &samples[current], &samples[proposed]);
            ++consecutiveRejects;
        }
        if (--progressCounter == 0) {
            progress.Update();
            progressCounter = progressUpdateFrequency;
        }
    }
    Assert(pixelNumOffset == nPixels);
    // Update display for recently computed Metropolis samples
    int ntf = AtomicAdd(&renderer->nTasksFinished, 1);
    int64_t totalSamples = int64_t(nPixels) * int64_t(nPixelSamples);
    float splatScale = float(double(totalSamples) / double(ntf * nTaskSamples));
    camera->film->UpdateDisplay(x0, y0, x1, y1, splatScale);
    if ((taskNum % 8) == 0) {
        MutexLock lock(*filmMutex);
        camera->film->WriteImage(splatScale);
    }
    PBRT_MLT_FINISHED_MLT_TASK(this);
}


Spectrum MetropolisRenderer::Li(const Scene *scene, const RayDifferential &ray,
        const Sample *sample, RNG &rng, MemoryArena &arena, Intersection *isect,
        Spectrum *T) const {
    Spectrum localT;
    if (!T) T = &localT;
    Intersection localIsect;
    if (!isect) isect = &localIsect;
    Spectrum Lo = 0.f;
    if (scene->Intersect(ray, isect))
        Lo = directLighting->Li(scene, this, ray, *isect, sample,
                                rng, arena);
    else {
        // Handle ray that doesn't intersect any geometry
        for (uint32_t i = 0; i < scene->lights.size(); ++i)
           Lo += scene->lights[i]->Le(ray);
    }
    return Lo;
}


Spectrum MetropolisRenderer::Transmittance(const Scene *scene, const RayDifferential &ray,
    const Sample *sample, RNG &rng, MemoryArena &arena) const {
    return 1.f;
}


