
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
               const vector<PathVertex> *eyePath,
               const vector<PathVertex> *lightPath);
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
    bool ignoreDirect;
    const MLTSample &initialSample;
    const Scene *scene;
    const Camera *camera;
    MetropolisRenderer *renderer;
    Mutex *filmMutex;
    Distribution1D *lightDistribution;
};


struct PathSample {
    float bsdfComponent, bsdfDir0, bsdfDir1;
    float rrSample;
};



struct LightingSample {
    float bsdfLightComponent;
    float bsdfLightDir0, bsdfLightDir1;
    float lightNum0, lightNum1;
    float lightDir0, lightDir1;
};



struct MLTSample {
    MLTSample(int maxDepth) {
        eyePathSamples.resize(maxDepth);
        lightPathSamples.resize(maxDepth);
        lightingSamples.resize(maxDepth);
    }
    CameraSample cameraSample;
    float lightNumSample, lightRaySamples[5];
    vector<PathSample> eyePathSamples, lightPathSamples;
    vector<LightingSample> lightingSamples;
};


static void LargeStep(RNG &rng, MLTSample *sample, int maxDepth,
        float x, float y, float t0, float t1, bool bidirectional) {
    sample->cameraSample.imageX = x;
    sample->cameraSample.imageY = y;
    sample->cameraSample.time = Lerp(rng.RandomFloat(), t0, t1);
    sample->cameraSample.lensU = rng.RandomFloat();
    sample->cameraSample.lensV = rng.RandomFloat();

    for (int i = 0; i < maxDepth; ++i) {
        // Apply large step to $i$th eye _PathSample_
        PathSample &eps = sample->eyePathSamples[i];
        eps.bsdfComponent = rng.RandomFloat();
        eps.bsdfDir0 = rng.RandomFloat();
        eps.bsdfDir1 = rng.RandomFloat();
        eps.rrSample = rng.RandomFloat();

        // Apply large step to $i$th _LightingSample_
        LightingSample &ls = sample->lightingSamples[i];
        ls.bsdfLightComponent = rng.RandomFloat();
        ls.bsdfLightDir0 = rng.RandomFloat();
        ls.bsdfLightDir1 = rng.RandomFloat();
        ls.lightNum0 = rng.RandomFloat();
        ls.lightNum1 = rng.RandomFloat();
        ls.lightDir0 = rng.RandomFloat();
        ls.lightDir1 = rng.RandomFloat();
    }

    if (bidirectional) {
        sample->lightNumSample = rng.RandomFloat();
        for (int i = 0; i < 5; ++i)
            sample->lightRaySamples[i] = rng.RandomFloat();
        for (int i = 0; i < maxDepth; ++i) {
            // Apply large step to $i$th light _PathSample_
            PathSample &lps = sample->lightPathSamples[i];
            lps.bsdfComponent = rng.RandomFloat();
            lps.bsdfDir0 = rng.RandomFloat();
            lps.bsdfDir1 = rng.RandomFloat();
            lps.rrSample = rng.RandomFloat();
        }
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
        int x0, int x1, int y0, int y1, float t0, float t1, bool bidirectional) {
    mutate(rng, &sample->cameraSample.imageX, x0, x1);
    mutate(rng, &sample->cameraSample.imageY, y0, y1);
    mutate(rng, &sample->cameraSample.time, t0, t1);
    mutate(rng, &sample->cameraSample.lensU);
    mutate(rng, &sample->cameraSample.lensV);

    for (int i = 0; i < maxDepth; ++i) {
        // Apply small step to $i$th eye _PathSample_
        PathSample &eps = sample->eyePathSamples[i];
        mutate(rng, &eps.bsdfComponent);
        mutate(rng, &eps.bsdfDir0);
        mutate(rng, &eps.bsdfDir1);
        mutate(rng, &eps.rrSample);

        // Apply small step to $i$th _LightingSample_
        LightingSample &ls = sample->lightingSamples[i];
        mutate(rng, &ls.bsdfLightComponent);
        mutate(rng, &ls.bsdfLightDir0);
        mutate(rng, &ls.bsdfLightDir1);
        mutate(rng, &ls.lightNum0);
        mutate(rng, &ls.lightNum1);
        mutate(rng, &ls.lightDir0);
        mutate(rng, &ls.lightDir1);
    }

    if (bidirectional) {
        mutate(rng, &sample->lightNumSample);
        for (int i = 0; i < 5; ++i)
            mutate(rng, &sample->lightRaySamples[i]);
        for (int i = 0; i < maxDepth; ++i) {
            // Apply small step to $i$th light _PathSample_
            PathSample &lps = sample->lightPathSamples[i];
            mutate(rng, &lps.bsdfComponent);
            mutate(rng, &lps.bsdfDir0);
            mutate(rng, &lps.bsdfDir1);
            mutate(rng, &lps.rrSample);
        }
    }
}


struct PathVertex {
    Intersection isect;
    BSDF *bsdf;
    Spectrum alpha;
    Vector wPrev, wNext;
    bool specularBounce;
};


static void SampleIL(const MLTSample &sample, const Scene *scene,
    MemoryArena &arena, const Camera *camera, const Renderer *renderer,
    const Distribution1D *lightDistribution, vector<PathVertex> &eyePath,
    vector<PathVertex> &lightPath, RNG &rng, bool ignoreDirect,
    bool bidirectional, Spectrum *L, float *I);
static uint32_t GeneratePath(const RayDifferential &r, const Spectrum &alpha,
    const Scene *scene, MemoryArena &arena, const vector<PathSample> &sample,
    vector<PathVertex> &path, RayDifferential *escapedRay,
    Spectrum *escapedAlpha, bool ignoreDirect);
static Spectrum L(const Scene *scene, const Renderer *renderer,
    const vector<PathVertex> &path, int pathLength,
    MemoryArena &arena, RNG &rng, bool ignoreDirect, const vector<LightingSample> &samples,
    float time, const Distribution1D *lightDistribution, const RayDifferential &escapedRay,
    const Spectrum &escapedAlpha);
static Spectrum L(const Scene *scene, const Renderer *renderer,
    const vector<PathVertex> &eyePath, int eyePathLength,
    const vector<PathVertex> &lightPath, int lightPathLength,
    MemoryArena &arena, RNG &rng, bool ignoreDirect, const vector<LightingSample> &samples,
    float time, const Distribution1D *lightDistribution, const RayDifferential &escapedRay,
    const Spectrum &escapedAlpha);

// Metropolis Method Definitions
MetropolisRenderer::MetropolisRenderer(int perPixelSamples,
        int nboot, int dps, float lsp, bool dds, int mr, int md,
        Camera *c, bool db) {
    camera = c;

    nPixelSamples = perPixelSamples;
    largeStepProbability = lsp;
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

    int x0, x1, y0, y1;
    camera->film->GetPixelExtent(&x0, &x1, &y0, &y1);

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
        vector<PathVertex> eyePath(maxDepth, PathVertex());
        vector<PathVertex> lightPath(maxDepth, PathVertex());
        float sumContrib = 0.f;
        bootstrapSamples.reserve(nBootstrap);
        MLTSample sample(maxDepth);
        for (uint32_t i = 0; i < nBootstrap; ++i) {
            // Generate random sample and path for MLT bootstrapping
            float x = Lerp(rng.RandomFloat(), x0, x1);
            float y = Lerp(rng.RandomFloat(), y0, y1);
            LargeStep(rng, &sample, maxDepth, x, y, t0, t1, bidirectional);
            Spectrum L;
            float I;
            SampleIL(sample, scene, arena, camera, this, lightDistribution, eyePath,
                     lightPath, rng, directLighting != NULL, bidirectional,
                     &L, &I);

            // Compute contribution for random sample for MLT bootstrapping
            sumContrib += I;
            bootstrapSamples.push_back(I);
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
            LargeStep(rng, &initialSample, maxDepth, x, y, t0, t1, bidirectional);
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


inline float I(const Spectrum &L, const MLTSample &sample,
               const vector<PathVertex> *eyePath,
               const vector<PathVertex> *lightPath) {
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
    ignoreDirect = id;
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
    int64_t totalSamples = int64_t(nPixels) * int64_t(nPixelSamples);
    RNG rng(taskNum);
    uint32_t pixelNumOffset = 0;
    
    vector<PathVertex> eyePath, lightPath;
    eyePath.reserve(renderer->maxDepth);
    lightPath.reserve(renderer->maxDepth);
    
    MemoryArena arena[2];
    vector<MLTSample> samples(2, MLTSample(renderer->maxDepth));
    Spectrum L[2];
    float I[2];
    uint32_t current = 0, proposed = 1;
    samples[current] = initialSample;
    
    uint32_t consecutiveRejects = 0;
    uint32_t largeStepRate = nPixelSamples / renderer->largeStepsPerPixel;
    Assert(largeStepRate > 1);
    
    uint32_t progressCounter = progressUpdateFrequency;

    // Compute _L[current]_ for initial sample
    SampleIL(initialSample, scene, arena[current], camera, renderer,
             lightDistribution, eyePath, lightPath,
             rng, ignoreDirect, renderer->bidirectional, &L[current], &I[current]);

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
            LargeStep(rng, &samples[proposed], renderer->maxDepth,
                      x + dx, y + dy, t0, t1, renderer->bidirectional);
            ++pixelNumOffset;
        }
        else
            SmallStep(rng, &samples[proposed], renderer->maxDepth,
                      x0, x1, y0, y1, t0, t1, renderer->bidirectional);

        // Compute contribution of proposed sample
        SampleIL(samples[proposed], scene, arena[proposed], camera, renderer,
                 lightDistribution, eyePath, lightPath,
                 rng, ignoreDirect, renderer->bidirectional, &L[proposed], &I[proposed]);

        // Compute acceptance probability for proposed sample
        float a = min(1.f, I[proposed] / I[current]);
        float currentWeight = (1.f - a) /
                              (I[current] / b + renderer->largeStepProbability) /
                              float(nPixelSamples);
        float proposedWeight = (a + (largeStep ? 1.f : 0.f)) /
                               (I[proposed] / b + renderer->largeStepProbability) /
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
        else {
            PBRT_MLT_REJECTED_MUTATION(a, &samples[current], &samples[proposed]);
            ++consecutiveRejects;
        }
        arena[proposed].FreeAll();
        if (--progressCounter == 0) {
            progress.Update();
              progressCounter = progressUpdateFrequency;
        }
    }
    Assert(pixelNumOffset == nPixels);
    // Update display for recently computed Metropolis samples
    int ntf = AtomicAdd(&renderer->nTasksFinished, 1);
    float splatScale = float(double(totalSamples) / double(ntf * nTaskSamples));
    camera->film->UpdateDisplay(x0, y0, x1, y1, splatScale);
    if ((taskNum % 8) == 0) {
        MutexLock lock(*filmMutex);
        camera->film->WriteImage(splatScale);
    }
    PBRT_MLT_FINISHED_MLT_TASK(this);
}


static void SampleIL(const MLTSample &sample, const Scene *scene,
        MemoryArena &arena, const Camera *camera, const Renderer *renderer,
        const Distribution1D *lightDistribution, vector<PathVertex> &eyePath,
        vector<PathVertex> &lightPath, RNG &rng, bool ignoreDirect,
        bool bidirectional, Spectrum *L, float *I) {
    RayDifferential eyeRay;
    float eyeWt = camera->GenerateRayDifferential(sample.cameraSample, &eyeRay);
    RayDifferential escapedRay;
    Spectrum escapedAlpha;
    uint32_t eyeLength = GeneratePath(eyeRay, eyeWt, scene, arena,
        sample.eyePathSamples, eyePath, &escapedRay, &escapedAlpha,
        ignoreDirect);

    if (bidirectional) {
        // Compute radiance along paths using bidirectional path tracing
        float lightPdf, lightRayPdf;
        uint32_t lightNum = lightDistribution->SampleDiscrete(sample.lightNumSample, &lightPdf);
        const Light *light = scene->lights[lightNum];
        Ray lightRay;
        Normal Nl;
        Spectrum lightWt = light->Sample_L(scene,
            LightSample(sample.lightRaySamples[0],
                        sample.lightRaySamples[1],
                        sample.lightRaySamples[2]),
            sample.lightRaySamples[3], sample.lightRaySamples[4], sample.cameraSample.time,
            &lightRay, &Nl, &lightRayPdf);
        uint32_t lightLength = 0;
        if (lightWt.IsBlack() || lightRayPdf == 0.f) {
            *L = 0.f;
            *I = 0.f;
        }
        else {
            lightWt *= AbsDot(Normalize(Nl), lightRay.d) / (lightPdf * lightRayPdf);
            lightLength = GeneratePath(RayDifferential(lightRay), lightWt, scene, arena,
               sample.lightPathSamples, lightPath, NULL, NULL, ignoreDirect);
            *L = ::L(scene, renderer, eyePath, eyeLength, lightPath, lightLength,
                     arena, rng, ignoreDirect, sample.lightingSamples,
                     sample.cameraSample.time,
                     lightDistribution, escapedRay, escapedAlpha);
            *I = ::I(*L, sample, &eyePath, &lightPath);
        }
    }
    else {
        *L = ::L(scene, renderer, eyePath, eyeLength, arena, rng, ignoreDirect,
                sample.lightingSamples, sample.cameraSample.time,
                lightDistribution, escapedRay, escapedAlpha);
        *I = ::I(*L, sample, &eyePath, NULL);
    }
}


static uint32_t GeneratePath(const RayDifferential &r, const Spectrum &a, const Scene *scene,
        MemoryArena &arena, const vector<PathSample> &sample,
        vector<PathVertex> &path, RayDifferential *escapedRay, Spectrum *escapedAlpha,
        bool ignoreDirect) {
    RayDifferential ray = r;
    Spectrum alpha = a;
    bool allSpecular = true;
    if (escapedAlpha) *escapedAlpha = 0.f;
    uint32_t length = 0;
    for (; length < sample.size(); ++length) {
        if (!scene->Intersect(ray, &path[length].isect)) {
            // Handle ray that leaves the scene during path generation
            bool includeLe = ignoreDirect ? (!allSpecular && path[length-1].specularBounce) :
                                    (length == 0 || path[length-1].specularBounce);
            if (includeLe) {
                if (escapedAlpha) *escapedAlpha = alpha;
                if (escapedRay) *escapedRay = ray;
            }
            break;
        }
        // Record information for current path vertex
        PathVertex &v = path[length];
        v.alpha = alpha;
        BSDF *bsdf = v.isect.GetBSDF(ray, arena);
        v.bsdf = bsdf;
        v.wPrev = -ray.d;

        // Sample direction for outgoing Metropolis path direction
        const Point &p = bsdf->dgShading.p;
        const Normal &n = bsdf->dgShading.nn;
        Vector wo = -ray.d;
        const PathSample &ps = sample[length];
        BSDFSample outgoingBSDFSample(ps.bsdfDir0, ps.bsdfDir1,
                                      ps.bsdfComponent);
        float pdf;
        BxDFType flags;
        Spectrum f = bsdf->Sample_f(wo, &v.wNext, outgoingBSDFSample,
                                    &pdf, BSDF_ALL, &flags);
        if (f.IsBlack() || pdf == 0.f) break;
        v.specularBounce = (flags & BSDF_SPECULAR) != 0;
        allSpecular &= v.specularBounce;

        // Terminate path with RR or prepare for finding next vertex
        Spectrum pathScale = f * AbsDot(v.wNext, n) / pdf;
        float rrSurviveProb = min(1.f, pathScale.y());
        if (ps.rrSample > rrSurviveProb)
            return length+1;
        alpha *= pathScale / rrSurviveProb;
        ray = RayDifferential(p, v.wNext, ray, v.isect.rayEpsilon);
        //alpha *= renderer->Transmittance(scene, ray, NULL, rng, arena);
    }
    return length;
}



static Spectrum L(const Scene *scene, const Renderer *renderer,
        const vector<PathVertex> &path, int length,
        MemoryArena &arena, RNG &rng, bool ignoreDirect, const vector<LightingSample> &samples,
        float time, const Distribution1D *lightDistribution, const RayDifferential &escapedRay,
        const Spectrum &escapedAlpha) {
    Spectrum L = 0.;
    bool allSpecular = true;
    for (int e = 0; e < length; ++e) {
        const PathVertex &ve = path[e];
        const Point &pe = ve.bsdf->dgShading.p;
        const Normal &ne = ve.bsdf->dgShading.nn;
        // Compute direct illumination for Metropolis path vertex
        const LightingSample &ls = samples[e];
        if (ignoreDirect ? (ve.specularBounce && !allSpecular) :
                           (ve.specularBounce || e == 0))
            L += ve.alpha * ve.isect.Le(ve.wPrev);
        Spectrum Ld(0.f);
        if (!ignoreDirect || !allSpecular) {
            LightSample lightSample(ls.lightDir0, ls.lightDir1, ls.lightNum0);
            BSDFSample bsdfSample(ls.bsdfLightDir0, ls.bsdfLightDir1,
                                  ls.bsdfLightComponent);
            float lightPdf;
            uint32_t lightNum = lightDistribution->SampleDiscrete(ls.lightNum1, &lightPdf);
            const Light *light = scene->lights[lightNum];
            Ld = ve.alpha *
                 EstimateDirect(scene, renderer, arena, light, pe, ne, ve.wPrev,
                     ve.isect.rayEpsilon, time, ve.bsdf, rng,
                     lightSample, bsdfSample) / lightPdf;
        }
        allSpecular &= ve.specularBounce;
        L += Ld;
    }
    // Add contribution of escaped ray, if any
    if (!escapedAlpha.IsBlack())
        for (uint32_t i = 0; i < scene->lights.size(); ++i)
           L += escapedAlpha * scene->lights[i]->Le(escapedRay);
    return L;
}


static Spectrum L(const Scene *scene, const Renderer *renderer,
        const vector<PathVertex> &eyePath, int eyePathLength,
        const vector<PathVertex> &lightPath, int lightPathLength,
        MemoryArena &arena, RNG &rng, bool ignoreDirect, const vector<LightingSample> &samples,
        float time, const Distribution1D *lightDistribution, const RayDifferential &escapedRay,
        const Spectrum &escapedAlpha) {
    Spectrum L = 0.;
    bool allSpecular = true;
    for (int e = 0; e < eyePathLength; ++e) {
        const PathVertex &ve = eyePath[e];
        const Point &pe = ve.bsdf->dgShading.p;
        const Normal &ne = ve.bsdf->dgShading.nn;
        // Compute direct illumination for Metropolis path vertex
        const LightingSample &ls = samples[e];
        if (ignoreDirect ? (ve.specularBounce && !allSpecular) :
                           (ve.specularBounce || e == 0))
            L += ve.alpha * ve.isect.Le(ve.wPrev);
        Spectrum Ld(0.f);
        if (!ignoreDirect || !allSpecular) {
            LightSample lightSample(ls.lightDir0, ls.lightDir1, ls.lightNum0);
            BSDFSample bsdfSample(ls.bsdfLightDir0, ls.bsdfLightDir1,
                                  ls.bsdfLightComponent);
            float lightPdf;
            uint32_t lightNum = lightDistribution->SampleDiscrete(ls.lightNum1, &lightPdf);
            const Light *light = scene->lights[lightNum];
            Ld = ve.alpha *
                 EstimateDirect(scene, renderer, arena, light, pe, ne, ve.wPrev,
                     ve.isect.rayEpsilon, time, ve.bsdf, rng,
                     lightSample, bsdfSample) / lightPdf;
        }
        allSpecular &= ve.specularBounce;
        L += Ld / (e + 1);
        for (int l = 0; l < lightPathLength; ++l) {
            const PathVertex &vl = lightPath[l];
            const Point &pl = vl.bsdf->dgShading.p;
            const Normal &nl = vl.bsdf->dgShading.nn;
            // Compute bidirectional contribution between eye and light vertices
            Vector w = Normalize(pl - pe);
            Spectrum fe = ve.bsdf->f(ve.wPrev, w);
            Spectrum fl = vl.bsdf->f(-w, vl.wPrev);
            if (fe.IsBlack() || fl.IsBlack()) continue;
            float G = AbsDot(ne, w) * AbsDot(nl, w) / DistanceSquared(pl, pe);
            Ray r(pe, pl - pe, 1e-3f, .999f, time);
            if (!scene->IntersectP(r))
                L += (ve.alpha * fe * G * fl * vl.alpha) / (e + l + 2);
        }
    }
    // Add contribution of escaped ray, if any
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


