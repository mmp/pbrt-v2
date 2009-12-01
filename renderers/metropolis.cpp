
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
inline float I(const Spectrum &L, const MLTSample &sample) {
    return L.y();
}


class MLTTask : public Task {
public:
    MLTTask(ProgressReporter &prog, int tn, int ns, int64_t ts, int np,
        int xx0, int xx1, int yy0, int yy1, float tt0, float tt1,
        float bb, float lsp, const MLTSample &is, bool id, int mr,
        int md, const Scene *sc, const Camera *c, const Renderer *renderer,
        volatile float *nSamplesFinished, Mutex *filmMutex,
        Distribution1D *lightDistribution);
    void Run();

private:
    ProgressReporter &progress;
    int taskNum, nSamples, nPixels;
    int64_t totalSamples;
    int x0, x1, y0, y1;
    float t0, t1;
    float b, largeStepProbability;
    bool ignoreDirect;
    const MLTSample &initialSample;
    const Scene *scene;
    const Camera *camera;
    const Renderer *renderer;
    int maxConsecutiveRejects;
    int maxDepth;
    volatile float *nSamplesFinished;
    Mutex *filmMutex;
    Distribution1D *lightDistribution;
};


struct PathSample {
    float bsdfComponent, bsdfDir0, bsdfDir1;
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
        int x0, int x1, int y0, int y1, float t0, float t1) {
    sample->cameraSample.imageX = Lerp(rng.RandomFloat(), x0, x1);
    sample->cameraSample.imageY = Lerp(rng.RandomFloat(), y0, y1);
    sample->cameraSample.time = Lerp(rng.RandomFloat(), t0, t1);
    sample->cameraSample.lensU = rng.RandomFloat();
    sample->cameraSample.lensV = rng.RandomFloat();
    for (int i = 0; i < maxDepth; ++i) {
        // Apply large step to $i$th _PathSample_
        PathSample &ps = sample->pathSamples[i];
        ps.bsdfComponent = rng.RandomFloat();
        ps.bsdfDir0 = rng.RandomFloat();
        ps.bsdfDir1 = rng.RandomFloat();
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
        mutate(rng, &ps.bsdfLightComponent);
        mutate(rng, &ps.bsdfLightDir0);
        mutate(rng, &ps.bsdfLightDir1);
        mutate(rng, &ps.lightNum0);
        mutate(rng, &ps.lightNum1);
        mutate(rng, &ps.lightDir0);
        mutate(rng, &ps.lightDir1);
    }
}


static Spectrum L(const Scene *scene, const Renderer *renderer,
    const Camera *camera, MemoryArena &arena,
    RNG &rng, int maxDepth, bool ignoreDirect, const MLTSample &sample,
    const Distribution1D *lightDistribution);

// Metropolis Method Definitions
MetropolisRenderer::MetropolisRenderer(int ts, int perPixelSamples,
    int nboot, int dps, float lsp, bool dds, int mr, int md,
    Camera *c) {
    camera = c;
    if (perPixelSamples > 0) {
        int x0, x1, y0, y1;
        camera->film->GetPixelExtent(&x0, &x1, &y0, &y1);
        int64_t nPixels = (x1-x0) * (y1-y0);
        nSamples = int64_t(perPixelSamples) * nPixels;
        if (ts != 0) {
            nSamples = max(uint64_t(ts), nSamples);
            Warning("Both a total sample count (%d) and a "
                    "per-pixel sample count (%d) specified for "
                    "MetropolisRenderer.  Using maximum as total (%lld).",
                    ts, perPixelSamples, nSamples);
        }
    }
    else
        nSamples = ts;
    nBootstrap = nboot;
    directPixelSamples = dps;
    largeStepProbability = lsp;
    doDirectSeparately = dds;
    maxDepth = md;
    maxConsecutiveRejects = mr;
    nSamplesFinished = 0.f;
    directLighting = new DirectLightingIntegrator(SAMPLE_ALL_UNIFORM, maxDepth);
}


MetropolisRenderer::~MetropolisRenderer() {
    delete camera;
    delete directLighting;
}


MetropolisRenderer *CreateMetropolisRenderer(const ParamSet &params,
        Camera *camera) {
    float largeStepProbability = params.FindOneFloat("largestepprobability", .5f);
    int nSamples = params.FindOneInt("nsamples", 0);
    int perPixelSamples = params.FindOneInt("samplesperpixel", 100);
    int nBootstrap = params.FindOneInt("bootstrapsamples", 100000);
    int directPixelSamples = params.FindOneInt("directsamples", 4);
    bool doDirectSeparately = params.FindOneBool("dodirectseparately", true);
    int mr = params.FindOneInt("maxconsecutiverejects", 512);
    int md = params.FindOneInt("maxdepth", 4);

    if (PbrtOptions.quickRender) {
        if (nSamples > 0)
            nSamples = max(1, nSamples / 4);
        if (perPixelSamples > 0)
            perPixelSamples = max(1, perPixelSamples / 4);
        nBootstrap = max(1, nBootstrap / 4);
        directPixelSamples = max(1, directPixelSamples / 4);
    }

    return new MetropolisRenderer(nSamples, perPixelSamples, nBootstrap,
        directPixelSamples, largeStepProbability, doDirectSeparately,
        mr, md, camera);
}


void MetropolisRenderer::Render(const Scene *scene) {
    int x0, x1, y0, y1;
    camera->film->GetPixelExtent(&x0, &x1, &y0, &y1);
    int nPixels = (x1-x0) * (y1-y0);
    float t0 = camera->shutterOpen;
    float t1 = camera->shutterClose;
    Distribution1D *lightDistribution = ComputeLightSamplingCDF(scene);

    if (doDirectSeparately) {
        // Compute direct lighting before Metropolis light transport
        LDSampler sampler(x0, x1, y0, y1, directPixelSamples, t0, t1);
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
    float sumContrib = 0.f;
    bootstrapSamples.reserve(nBootstrap);
    MLTSample sample(maxDepth);
    for (int i = 0; i < nBootstrap; ++i) {
        // Compute contribution for random sample for MLT bootstrapping
        LargeStep(rng, &sample, maxDepth, x0, x1, y0, y1, t0, t1);
        float contrib = I(L(scene, this, camera, arena, rng, maxDepth,
                            doDirectSeparately, sample, lightDistribution), sample);
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
    for (int i = 0; i < nBootstrap; ++i) {
        LargeStep(rng, &initialSample, maxDepth, x0, x1, y0, y1, t0, t1);
        sumContrib += bootstrapSamples[i];
        if (contribOffset < sumContrib)
            break;
    }

    // Launch tasks to generate Metropolis samples
    if (scene->lights.size() > 0) {
        int nTasks = int(nSamples / 50000);
        nTasks = max(nTasks, 32 * NumSystemCores());
        nTasks = min(nTasks, 32768);
        nSamples = (nSamples / nTasks) * nTasks;
        ProgressReporter progress(nTasks, "Metropolis");
        vector<Task *> tasks;
        Mutex *filmMutex = Mutex::Create();
        for (int i = 0; i < nTasks; ++i)
            tasks.push_back(new MLTTask(progress, i, int(nSamples/nTasks), nSamples,
                nPixels, x0, x1, y0, y1, t0, t1, b, largeStepProbability, initialSample,
                doDirectSeparately, maxConsecutiveRejects, maxDepth, scene, camera, this,
                &nSamplesFinished, filmMutex, lightDistribution));
        EnqueueTasks(tasks);
        WaitForAllTasks();
        for (uint32_t i = 0; i < tasks.size(); ++i)
            delete tasks[i];
        progress.Done();
    }
    camera->film->WriteImage();
}


MLTTask::MLTTask(ProgressReporter &prog, int tn, int ns, int64_t ts, int np,
    int xx0, int xx1, int yy0, int yy1, float tt0, float tt1,
    float bb, float lsp, const MLTSample &is, bool id, int mr,
    int md, const Scene *sc, const Camera *c, const Renderer *ren,
    volatile float *nfinished, Mutex *fm, Distribution1D *ld)
    : progress(prog), initialSample(is) {
    taskNum = tn;
    nSamples = ns;
    totalSamples = ts;
    nPixels = np;
    x0 = xx0;
    x1 = xx1;
    y0 = yy0;
    y1 = yy1;
    t0 = tt0;
    t1 = tt1;
    b = bb;
    largeStepProbability = lsp;
    ignoreDirect = id;
    maxConsecutiveRejects = mr;
    maxDepth = md;
    scene = sc;
    camera = c;
    renderer = ren;
    nSamplesFinished = nfinished;
    filmMutex = fm;
    lightDistribution = ld;
}


void MLTTask::Run() {
    PBRT_MLT_STARTED_MLT_TASK(this);
    // Declare basic _MLTTask_ variables and prepare for sampling
    RNG rng(taskNum);
    MemoryArena arena;
    vector<MLTSample> mltSamples(2, MLTSample(maxDepth));
    Spectrum sampleLs[2];
    uint32_t currentSample = 0, proposedSample = 1;
    mltSamples[currentSample] = initialSample;
    sampleLs[currentSample] = L(scene, renderer, camera, arena, rng, maxDepth,
                                ignoreDirect, mltSamples[currentSample], lightDistribution);
    int consecutiveRejects = 0;
    for (int sampleNum = 0; sampleNum < nSamples; ++sampleNum) {
        // Compute proposed mutation to current sample
        bool largeStep = rng.RandomFloat() < largeStepProbability;
        mltSamples[proposedSample] = mltSamples[currentSample];
        if (largeStep)
            LargeStep(rng, &mltSamples[proposedSample], maxDepth,
                      x0, x1, y0, y1, t0, t1);
        else
            SmallStep(rng, &mltSamples[proposedSample], maxDepth,
                      x0, x1, y0, y1, t0, t1);

        // Compute contribution of proposed sample and acceptance probability
        sampleLs[proposedSample] = L(scene, renderer, camera, arena, rng, maxDepth,
                                     ignoreDirect, mltSamples[proposedSample], lightDistribution);
        float currentI = I(sampleLs[currentSample], mltSamples[currentSample]);
        float proposedI = I(sampleLs[proposedSample], mltSamples[proposedSample]);
        float a = min(1.f, proposedI / currentI);
        float currentWeight = (1.f - a) /
                              (currentI / b + largeStepProbability) *
                              float(nPixels) / float(totalSamples);
        float proposedWeight = (a + (largeStep ? 1.f : 0.f)) /
                               (proposedI / b + largeStepProbability) *
                               float(nPixels) / float(totalSamples);

        // Splat current and proposed samples to _Film_
        if (currentWeight > 0.f && currentI > 0.f)
            camera->film->Splat(mltSamples[currentSample].cameraSample,
                sampleLs[currentSample] * currentWeight);
        if (proposedWeight > 0.f && proposedI > 0.f)
            camera->film->Splat(mltSamples[proposedSample].cameraSample,
               sampleLs[proposedSample] * proposedWeight);

        // Randomly accept proposed path mutation (or not)
        if (consecutiveRejects >= maxConsecutiveRejects ||
            rng.RandomFloat() < a) {
            PBRT_MLT_ACCEPTED_MUTATION(a, &mltSamples[currentSample], &mltSamples[proposedSample]);
            currentSample ^= 1;
            proposedSample ^= 1;
            consecutiveRejects = 0;
        }
        else {
            PBRT_MLT_REJECTED_MUTATION(a, &mltSamples[currentSample], &mltSamples[proposedSample]);
            ++consecutiveRejects;
        }
        arena.FreeAll();
    }
    // Update display for recently computed Metropolis samples
    float nf = AtomicAdd(nSamplesFinished, nSamples);
    float splatScale = float(totalSamples)/nf;
    camera->film->UpdateDisplay(x0, y0, x1, y1, splatScale);
    if ((taskNum % 32) == 0) {
        MutexLock lock(*filmMutex);
        camera->film->WriteImage(splatScale);
    }
    progress.Update();
    PBRT_MLT_FINISHED_MLT_TASK(this);
}


static Spectrum L(const Scene *scene, const Renderer *renderer,
        const Camera *camera, MemoryArena &arena, RNG &rng, int maxDepth,
        bool ignoreDirect, const MLTSample &sample, const Distribution1D *lightDistribution) {
    // Generate camera ray from Metropolis sample
    RayDifferential ray;
    float cameraWeight = camera->GenerateRayDifferential(sample.cameraSample,
                                                         &ray);
    Spectrum pathThroughput = cameraWeight, L = 0.;
    bool specularBounce = false, allSpecular = true;
    for (int pathLength = 0; pathLength < maxDepth; ++pathLength) {
        // Find next intersection in Metropolis light path
        Intersection isect;
        if (!scene->Intersect(ray, &isect)) {
            bool includeLe = ignoreDirect ? (specularBounce && !allSpecular) :
                                            (pathLength == 0 || specularBounce);
            if (includeLe)
                for (uint32_t i = 0; i < scene->lights.size(); ++i)
                   L += pathThroughput * scene->lights[i]->Le(ray);
            break;
        }
        if (ignoreDirect ? (specularBounce && !allSpecular) :
                           (specularBounce || pathLength == 0))
            L += pathThroughput * isect.Le(-ray.d);
        BSDF *bsdf = isect.GetBSDF(ray, arena);
        const Point &p = bsdf->dgShading.p;
        const Normal &n = bsdf->dgShading.nn;
        Vector wo = -ray.d;
        const PathSample &ps = sample.pathSamples[pathLength];
        // Sample direct illumination for Metropolis path vertex
        if ((!ignoreDirect || !allSpecular) && scene->lights.size() > 0) {
            LightSample lightSample(ps.lightDir0, ps.lightDir1, ps.lightNum0);
            BSDFSample bsdfSample(ps.bsdfLightDir0, ps.bsdfLightDir1,
                                  ps.bsdfLightComponent);
            float lightPdf;
            uint32_t lightNum = lightDistribution->SampleDiscrete(ps.lightNum1, &lightPdf);
            const Light *light = scene->lights[lightNum];
            L += pathThroughput *
                 EstimateDirect(scene, renderer, arena, light, p, n, wo,
                     isect.rayEpsilon, sample.cameraSample.time, bsdf, rng,
                     lightSample, bsdfSample) / lightPdf;
        }

        // Sample direction for outgoing Metropolis path direction
        BSDFSample outgoingBSDFSample(ps.bsdfDir0, ps.bsdfDir1,
                                      ps.bsdfComponent);
        Vector wi;
        float pdf;
        BxDFType flags;
        Spectrum f = bsdf->Sample_f(wo, &wi, outgoingBSDFSample,
                                    &pdf, BSDF_ALL, &flags);
        if (f.IsBlack() || pdf == 0.)
            break;
        specularBounce = (flags & BSDF_SPECULAR) != 0;
        allSpecular &= specularBounce;
        pathThroughput *= f * AbsDot(wi, n) / pdf;
        ray = RayDifferential(p, wi, ray, isect.rayEpsilon);
        //pathThroughput *= renderer->Transmittance(scene, ray, NULL, rng, arena);
    }
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


