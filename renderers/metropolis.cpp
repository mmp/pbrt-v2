
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
struct PathVertex;
inline float I(const Spectrum &L, const MLTSample &sample,
               const PathVertex *path, int pathLength);
class MLTTask : public Task {
public:
    MLTTask(ProgressReporter &prog, uint32_t taskNum, uint32_t nps,
        float dx, float dy, int xx0, int xx1, int yy0, int yy1, float tt0, float tt1,
        float bb, float lsp, uint32_t lspp, const MLTSample &is, bool id, int mr,
        int md, const Scene *sc, const Camera *c, const Renderer *renderer,
        AtomicInt32 *nTasksFinished, Mutex *filmMutex,
        Distribution1D *lightDistribution);
    void Run();

private:
    ProgressReporter &progress;
    uint32_t taskNum;
    uint32_t nPixelSamples;
    float dx, dy;
    int currentPixelSample;
    int x0, x1, y0, y1;
    float t0, t1;
    float b, largeStepProbability;
    uint32_t largeStepsPerPixel;
    bool ignoreDirect;
    const MLTSample &initialSample;
    const Scene *scene;
    const Camera *camera;
    const Renderer *renderer;
    int maxConsecutiveRejects, maxDepth;
    AtomicInt32 *nTasksFinished;
    Mutex *filmMutex;
    Distribution1D *lightDistribution;
};


struct PathSample {
    float bsdfComponent, bsdfDir0, bsdfDir1;
    float rrSample;
    float bsdfLightComponent;
    float bsdfLightDir0, bsdfLightDir1;
    float lightNum0, lightNum1;
    float lightDir0, lightDir1;
};



struct MLTSample {
    MLTSample(int maxDepth) { pathSamples.resize(maxDepth); }
    CameraSample cameraSample;
    vector<PathSample> pathSamples;
};


static void LargeStep(RNG &rng, MLTSample *sample, int maxDepth,
        float x, float y, float t0, float t1) {
    sample->cameraSample.imageX = x;
    sample->cameraSample.imageY = y;
    sample->cameraSample.time = Lerp(rng.RandomFloat(), t0, t1);
    sample->cameraSample.lensU = rng.RandomFloat();
    sample->cameraSample.lensV = rng.RandomFloat();
    for (int i = 0; i < maxDepth; ++i) {
        // Apply large step to $i$th _PathSample_
        PathSample &ps = sample->pathSamples[i];
        ps.bsdfComponent = rng.RandomFloat();
        ps.bsdfDir0 = rng.RandomFloat();
        ps.bsdfDir1 = rng.RandomFloat();
        ps.rrSample = rng.RandomFloat();
        ps.bsdfLightComponent = rng.RandomFloat();
        ps.bsdfLightDir0 = rng.RandomFloat();
        ps.bsdfLightDir1 = rng.RandomFloat();
        ps.lightNum0 = rng.RandomFloat();
        ps.lightNum1 = rng.RandomFloat();
        ps.lightDir0 = rng.RandomFloat();
        ps.lightDir1 = rng.RandomFloat();
    }
}


static inline void mutate(RNG &rng, float *v, float min = 0.f, float max = 1.f) {
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
        int x0, int x1, int y0, int y1, float t0, float t1) {
    mutate(rng, &sample->cameraSample.imageX, x0, x1);
    mutate(rng, &sample->cameraSample.imageY, y0, y1);
    mutate(rng, &sample->cameraSample.time, t0, t1);
    mutate(rng, &sample->cameraSample.lensU);
    mutate(rng, &sample->cameraSample.lensV);
    for (int i = 0; i < maxDepth; ++i) {
        // Apply small step to $i$th _PathSample_
        PathSample &ps = sample->pathSamples[i];
        mutate(rng, &ps.bsdfComponent);
        mutate(rng, &ps.bsdfDir0);
        mutate(rng, &ps.bsdfDir1);
        mutate(rng, &ps.rrSample);
        mutate(rng, &ps.bsdfLightComponent);
        mutate(rng, &ps.bsdfLightDir0);
        mutate(rng, &ps.bsdfLightDir1);
        mutate(rng, &ps.lightNum0);
        mutate(rng, &ps.lightNum1);
        mutate(rng, &ps.lightDir0);
        mutate(rng, &ps.lightDir1);
    }
}


struct PathVertex {
    Intersection isect;
    BSDF *bsdf;
    Spectrum alpha;
    Vector w;
    bool specularBounce;
};


static int GeneratePath(const RayDifferential &r, const Spectrum &alpha, const Scene *scene,
    MemoryArena &arena, int maxDepth, const MLTSample &sample, PathVertex *path,
    RayDifferential *escapedRay, Spectrum *escapedAlpha, bool ignoreDirect);
static Spectrum L(const Scene *scene, const Renderer *renderer,
    const Vector &wo, const PathVertex *path, int pathLength,
    MemoryArena &arena, RNG &rng, bool ignoreDirect, const MLTSample &sample,
    const Distribution1D *lightDistribution, const RayDifferential &escapedRay,
    const Spectrum &escapedAlpha);

// Metropolis Method Definitions
MetropolisRenderer::MetropolisRenderer(int perPixelSamples,
        int nboot, int dps, float lsp, bool dds, int mr, int md,
        Camera *c) {
    camera = c;

    nPixelSamples = perPixelSamples;
    largeStepProbability = lsp;
    largeStepsPerPixel = max(1u, RoundUpPow2(largeStepProbability * nPixelSamples));
    if (largeStepsPerPixel >= nPixelSamples) largeStepsPerPixel /= 2;
    Assert(largeStepsPerPixel >= 1 && largeStepsPerPixel < nPixelSamples);
    Assert(IsPowerOf2(largeStepsPerPixel));
printf("orig n pixel samples %d, large per (from %f) %d ", nPixelSamples, largeStepProbability, largeStepsPerPixel);
    if ((nPixelSamples % largeStepsPerPixel) != 0) {
        int origPixelSamples = nPixelSamples;
        nPixelSamples += largeStepsPerPixel - (nPixelSamples % largeStepsPerPixel);
        Warning("Rounding up to %d Metropolis samples per pixel (from %d)",
                nPixelSamples, origPixelSamples);
    }

    int x0, x1, y0, y1;
    camera->film->GetPixelExtent(&x0, &x1, &y0, &y1);
    nSamples = uint64_t(perPixelSamples) * uint64_t((x1-x0) * (y1-y0));

    nBootstrap = nboot;
    nDirectPixelSamples = dps;

    maxDepth = md;
    maxConsecutiveRejects = mr;
    nTasksFinished  = 0;
    directLighting = dds ? new DirectLightingIntegrator(SAMPLE_ALL_UNIFORM, maxDepth) : NULL;
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
    int md = params.FindOneInt("maxdepth", 4);

    if (PbrtOptions.quickRender) {
        perPixelSamples = max(1, perPixelSamples / 4);
        nBootstrap = max(1, nBootstrap / 4);
        nDirectPixelSamples = max(1, nDirectPixelSamples / 4);
    }

    return new MetropolisRenderer(perPixelSamples, nBootstrap,
        nDirectPixelSamples, largeStepProbability, doDirectSeparately,
        mr, md, camera);
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
        vector<PathVertex> path(maxDepth, PathVertex());
        float sumContrib = 0.f;
        bootstrapSamples.reserve(nBootstrap);
        MLTSample sample(maxDepth);
        for (uint32_t i = 0; i < nBootstrap; ++i) {
            // Generate random sample and path for MLT bootstrapping
            float x = Lerp(rng.RandomFloat(), x0, x1);
            float y = Lerp(rng.RandomFloat(), y0, y1);
            LargeStep(rng, &sample, maxDepth, x, y, t0, t1);
            RayDifferential ray;
            float cameraWeight = camera->GenerateRayDifferential(sample.cameraSample,
                                                                 &ray);
            RayDifferential escapedRay;
            Spectrum escapedAlpha;
            int pathLength = GeneratePath(ray, cameraWeight, scene, arena, maxDepth,
                                          sample, &path[0], &escapedRay, &escapedAlpha,
                                          directLighting != NULL);

            // Compute contribution for random sample for MLT bootstrapping
            Spectrum pathL = L(scene, this, -ray.d, &path[0], pathLength,
                               arena, rng, directLighting != NULL, sample,
                               lightDistribution, escapedRay, escapedAlpha);
            float contrib = I(pathL, sample, &path[0], pathLength);
            sumContrib += contrib;
            bootstrapSamples.push_back(contrib);
            arena.FreeAll();
        }
        float b = sumContrib / nBootstrap;

        // Select initial sample from bootstrap samples
        rng.Seed(0);
        float contribOffset = rng.RandomFloat() * sumContrib;
        sumContrib = 0.f;
        MLTSample initialSample(maxDepth);
        for (uint32_t i = 0; i < nBootstrap; ++i) {
            float x = Lerp(rng.RandomFloat(), x0, x1);
            float y = Lerp(rng.RandomFloat(), y0, y1);
            LargeStep(rng, &initialSample, maxDepth, x, y, t0, t1);
            sumContrib += bootstrapSamples[i];
            if (contribOffset < sumContrib)
                break;
        }

        // Launch tasks to generate Metropolis samples
        uint32_t nTasks = largeStepsPerPixel;
        ProgressReporter progress(nTasks, "Metropolis");
        vector<Task *> tasks;
        Mutex *filmMutex = Mutex::Create();
        Assert(IsPowerOf2(nTasks));
        uint32_t scramble[2] = { rng.RandomUInt(), rng.RandomUInt() };
        for (uint32_t i = 0; i < nTasks; ++i) {
            float d[2];
            Sample02(i, scramble, d);
            tasks.push_back(new MLTTask(progress, i, nPixelSamples,
                d[0], d[1], x0, x1, y0, y1, t0, t1, b, largeStepProbability,
                largeStepsPerPixel, initialSample,
                directLighting != NULL, maxConsecutiveRejects, maxDepth, scene, camera, this,
                &nTasksFinished, filmMutex, lightDistribution));
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


inline float I(const Spectrum &L, const MLTSample &sample,
               const PathVertex *path, int pathLength) {
    return L.y();
}


MLTTask::MLTTask(ProgressReporter &prog, uint32_t tn, uint32_t nps,
        float ddx, float ddy, int xx0, int xx1, int yy0, int yy1, float tt0, float tt1,
        float bb, float lsp, uint32_t lspp, const MLTSample &is, bool id, int mr,
        int md, const Scene *sc, const Camera *c, const Renderer *ren,
        AtomicInt32 *nfinished, Mutex *fm, Distribution1D *ld)
    : progress(prog), initialSample(is) {
    taskNum = tn;
    nPixelSamples = nps;
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
    largeStepProbability = lsp;
    largeStepsPerPixel = lspp;
    ignoreDirect = id;
    maxConsecutiveRejects = mr;
    maxDepth = md;
    scene = sc;
    camera = c;
    renderer = ren;
    nTasksFinished = nfinished;
    filmMutex = fm;
    lightDistribution = ld;
}


void MLTTask::Run() {
    PBRT_MLT_STARTED_MLT_TASK(this);
    // Declare basic _MLTTask_ variables and prepare for sampling
    uint32_t nPixels = (x1-x0) * (y1-y0);
    int64_t totalSamples = int64_t(nPixels) * int64_t(nPixelSamples);
    RNG rng(taskNum);
    uint32_t pixelNumOffset = 0;
    
    vector<PathVertex> path[2];
    path[0].resize(maxDepth);
    path[1].resize(maxDepth);
    int pathLength[2];
    
    MemoryArena arena[2];
    vector<MLTSample> samples(2, MLTSample(maxDepth));
    Spectrum sampleLs[2];
    uint32_t current = 0, proposed = 1;
    samples[current] = initialSample;
    
    int consecutiveRejects = 0;
    uint32_t largeStepRate = nPixelSamples / largeStepsPerPixel;
    Assert(largeStepRate > 1);

    // Compute _sampleLs[current]_ for initial sample
    // xxxx initial sample vs mltSmaples[current]
    RayDifferential ray;
    float cameraWeight = camera->GenerateRayDifferential(initialSample.cameraSample, &ray);
    RayDifferential escapedRay;
    Spectrum escapedAlpha;
    pathLength[current] = GeneratePath(ray, cameraWeight, scene, arena[current], maxDepth,
        initialSample, &path[current][0], &escapedRay, &escapedAlpha, ignoreDirect);
    sampleLs[current] = L(scene, renderer, -ray.d, &path[current][0],
        pathLength[current], arena[current], rng, ignoreDirect, initialSample,
        lightDistribution, escapedRay, escapedAlpha);

    // Compute randomly permuted table of pixel indices for large steps
    vector<int> largeStepPixelNum;
    for (uint32_t i = 0; i < nPixels; ++i) largeStepPixelNum.push_back(i);
    Shuffle(&largeStepPixelNum[0], nPixels, 1, rng);
    uint64_t nTaskSamples = uint64_t(nPixels) * uint64_t(largeStepRate);
    for (uint64_t sampleNum = 0; sampleNum < nTaskSamples; ++sampleNum) {
        // Compute proposed mutation to current sample
        samples[proposed] = samples[current];
        bool largeStep = ((sampleNum % largeStepRate) == 0);
        if (largeStep) {
            int x = x0 + largeStepPixelNum[pixelNumOffset] % (x1 - x0);
            int y = y0 + largeStepPixelNum[pixelNumOffset] / (x1 - x0);
            LargeStep(rng, &samples[proposed], maxDepth,
                      x + dx, y + dy, t0, t1);
            ++pixelNumOffset;
        }
        else
            SmallStep(rng, &samples[proposed], maxDepth,
                      x0, x1, y0, y1, t0, t1);

        // Compute contribution of proposed sample
        cameraWeight = camera->GenerateRayDifferential(samples[proposed].cameraSample, &ray);
        RayDifferential escapedRay;
        Spectrum escapedAlpha;
        pathLength[proposed] = GeneratePath(ray, cameraWeight, scene, arena[proposed], maxDepth,
            samples[proposed], &path[proposed][0], &escapedRay, &escapedAlpha, ignoreDirect);
        sampleLs[proposed] = L(scene, renderer, -ray.d, &path[proposed][0],
            pathLength[proposed], arena[proposed], rng, ignoreDirect,
            samples[proposed], lightDistribution, escapedRay, escapedAlpha);

        // Compute acceptance probability for proposed sample
        float currentI = I(sampleLs[current], samples[current],
                           &path[current][0], pathLength[current]);
        float proposedI = I(sampleLs[proposed], samples[proposed],
                            &path[proposed][0], pathLength[proposed]);
        float a = min(1.f, proposedI / currentI);
        float currentWeight = (1.f - a) /
                              (currentI / b + largeStepProbability) /
                              float(nPixelSamples);
        float proposedWeight = (a + (largeStep ? 1.f : 0.f)) /
                               (proposedI / b + largeStepProbability) /
                               float(nPixelSamples);

        // Splat current and proposed samples to _Film_
        if (currentWeight > 0.f && currentI > 0.f)
            camera->film->Splat(samples[current].cameraSample,
                                sampleLs[current] * currentWeight);
        if (proposedWeight > 0.f && proposedI > 0.f)
            camera->film->Splat(samples[proposed].cameraSample,
                                sampleLs[proposed] * proposedWeight);

        // Randomly accept proposed path mutation (or not)
        if (consecutiveRejects >= maxConsecutiveRejects ||
            rng.RandomFloat() < a) {
            PBRT_MLT_ACCEPTED_MUTATION(a, &samples[current], &samples[proposed]);
            current ^= 1;
            proposed ^= 1;
            consecutiveRejects = 0;
        }
        else {
            PBRT_MLT_REJECTED_MUTATION(a, &samples[current], &samples[proposed]);
            ++consecutiveRejects;
        }
        arena[proposed].FreeAll();
    }
    Assert(pixelNumOffset == nPixels);
    // Update display for recently computed Metropolis samples
    int ntf = AtomicAdd(nTasksFinished, 1);
    float splatScale = float(double(totalSamples) / double(ntf * nTaskSamples));
    camera->film->UpdateDisplay(x0, y0, x1, y1, splatScale);
    if ((taskNum % 8) == 0) {
        MutexLock lock(*filmMutex);
        camera->film->WriteImage(splatScale);
    }
    progress.Update();
    PBRT_MLT_FINISHED_MLT_TASK(this);
}


static int GeneratePath(const RayDifferential &r, const Spectrum &a, const Scene *scene,
        MemoryArena &arena, int maxDepth, const MLTSample &sample,
        PathVertex *path, RayDifferential *escapedRay, Spectrum *escapedAlpha,
        bool ignoreDirect) {
    RayDifferential ray = r;
    Spectrum alpha = a;
    bool allSpecular = true;
    *escapedAlpha = 0.f;
    int length = 0;
    for (; length < maxDepth; ++length) {
        if (!scene->Intersect(ray, &path[length].isect)) {
            // Handle ray that leaves the scene during path generation
            bool includeLe = ignoreDirect ? (!allSpecular && path[length-1].specularBounce) :
                                    (length == 0 || path[length-1].specularBounce);
            if (includeLe) {
                *escapedAlpha = alpha;
                *escapedRay = ray;
            }
            break;
        }
        // Record information for current path vertex
        PathVertex &v = path[length];
        v.alpha = alpha;
        BSDF *bsdf = v.isect.GetBSDF(ray, arena);
        v.bsdf = bsdf;

        // Sample direction for outgoing Metropolis path direction
        const Point &p = bsdf->dgShading.p;
        const Normal &n = bsdf->dgShading.nn;
        Vector wo = -ray.d;
        const PathSample &ps = sample.pathSamples[length];
        BSDFSample outgoingBSDFSample(ps.bsdfDir0, ps.bsdfDir1,
                                      ps.bsdfComponent);
        float pdf;
        BxDFType flags;
        Spectrum f = bsdf->Sample_f(wo, &v.w, outgoingBSDFSample,
                                    &pdf, BSDF_ALL, &flags);
        if (f.IsBlack() || pdf == 0.f) break;
        v.specularBounce = (flags & BSDF_SPECULAR) != 0;
        allSpecular &= v.specularBounce;

        // Terminate path with RR or prepare for finding next vertex
        Spectrum pathScale = f * AbsDot(v.w, n) / pdf;
        float rrSurviveProb = min(1.f, pathScale.y());
        if (ps.rrSample > rrSurviveProb)
            return length+1;
        alpha *= pathScale / rrSurviveProb;
        ray = RayDifferential(p, v.w, ray, v.isect.rayEpsilon);
        //alpha *= renderer->Transmittance(scene, ray, NULL, rng, arena);
    }
    return length;
}



static Spectrum L(const Scene *scene, const Renderer *renderer,
        const Vector &w, const PathVertex *path, int length,
        MemoryArena &arena, RNG &rng, bool ignoreDirect, const MLTSample &sample,
        const Distribution1D *lightDistribution, const RayDifferential &escapedRay,
        const Spectrum &escapedAlpha) {
    Spectrum L = 0.;
    Vector wo = w;
    bool allSpecular = true;
    for (int pathLength = 0; pathLength < length; ++pathLength) {
        const PathVertex &v = path[pathLength];
        // Compute direct illumination for Metropolis path vertex
        const BSDF *bsdf = v.bsdf;
        const Point &p = bsdf->dgShading.p;
        const Normal &n = bsdf->dgShading.nn;
        const PathSample &ps = sample.pathSamples[pathLength];
        if (ignoreDirect ? (v.specularBounce && !allSpecular) :
                           (v.specularBounce || pathLength == 0))
            L += v.alpha * v.isect.Le(wo);
        if (!ignoreDirect || !allSpecular) {
            LightSample lightSample(ps.lightDir0, ps.lightDir1, ps.lightNum0);
            BSDFSample bsdfSample(ps.bsdfLightDir0, ps.bsdfLightDir1,
                                  ps.bsdfLightComponent);
            float lightPdf;
            uint32_t lightNum = lightDistribution->SampleDiscrete(ps.lightNum1, &lightPdf);
            const Light *light = scene->lights[lightNum];
            L += v.alpha *
                 EstimateDirect(scene, renderer, arena, light, p, n, wo,
                     v.isect.rayEpsilon, sample.cameraSample.time, bsdf, rng,
                     lightSample, bsdfSample) / lightPdf;
        }
        allSpecular &= v.specularBounce;
        wo = -v.w;
    }
    if (!escapedAlpha.IsBlack())
        for (uint32_t i = 0; i < scene->lights.size(); ++i)
           L += escapedAlpha * scene->lights[i]->Le(escapedRay);
    return L;
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


