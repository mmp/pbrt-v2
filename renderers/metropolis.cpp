
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
#include "paramset.h"
#include "samplers/lowdiscrepancy.h"

// Metropolis Local Declarations
class MLTDirectIntegrator : public SurfaceIntegrator {
public:
    MLTDirectIntegrator() {
        lightSampleOffsets = NULL;
        bsdfSampleOffsets = NULL;
    }
    ~MLTDirectIntegrator() {
        delete[] lightSampleOffsets;
        delete[] bsdfSampleOffsets;
    }
    Spectrum Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential &ray, const Intersection &isect,
        const Sample *sample, RNG &rng, MemoryArena &arena) const;
    void RequestSamples(Sampler *sampler, Sample *sample, const Scene *scene);
private:
    LightSampleOffsets *lightSampleOffsets;
    BSDFSampleOffsets *bsdfSampleOffsets;
};


class MLTDirectTask : public Task {
public:
    MLTDirectTask(const Scene *sc, const Camera *c, const Renderer *ren,
        const MLTDirectIntegrator *di, Sampler *s, Sample *sa,
        int n, int nt, ProgressReporter &prog) : progress(prog) {
        scene = sc;
        camera = c;
        renderer = ren;
        direct = di;
        mainSampler = s;
        origSample = sa;
        taskNum = n;
        taskCount = nt;
    }
    void Run();

    const Scene *scene;
    const Camera *camera;
    const Renderer *renderer;
    const MLTDirectIntegrator *direct;
    Sampler *mainSampler;
    Sample *origSample;
    int taskNum, taskCount;
    ProgressReporter &progress;
};


inline float I(const Spectrum &L, const MLTSample &sample) {
    return L.y();
}


class MLTTask : public Task {
public:
    MLTTask(ProgressReporter &prog, int tn, int ns, int64_t ts, int np,
        int xx0, int xx1, int yy0, int yy1, float tt0, float tt1,
        float bb, float lsp, const MLTSample &is, bool id, int mr,
        int md, const Scene *sc, const Camera *c, const Renderer *renderer,
        volatile float *nSamplesFinished);
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
    RNG &rng, int maxDepth, bool ignoreDirect, const MLTSample &sample);

// Metropolis Method Definitions
MetropolisRenderer::MetropolisRenderer(int ts, int perPixelSamples,
    int nboot, int dps, float lsp, bool dds, int mr, int md,
    bool io, Camera *c) {
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
    indirectOnly = io;
    nSamplesFinished = 0.f;
}


MetropolisRenderer::~MetropolisRenderer() {
    delete camera;
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
    bool io = params.FindOneBool("indirectonly", false);

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
        mr, md, io, camera);
}


void MetropolisRenderer::Render(const Scene *scene) {
    int x0, x1, y0, y1;
    camera->film->GetPixelExtent(&x0, &x1, &y0, &y1);
    int nPixels = (x1-x0) * (y1-y0);
    float t0 = camera->shutterOpen;
    float t1 = camera->shutterClose;

    if (!indirectOnly && doDirectSeparately) {
        // Compute direct lighting before Metropolis light transport
        LDSampler sampler(x0, x1, y0, y1, directPixelSamples, t0, t1);
        MLTDirectIntegrator directLighting;
        Sample *sample = new Sample(&sampler, &directLighting, NULL, scene);
        vector<Task *> directTasks;
        int nDirectTasks = max(32 * NumSystemCores(),
                         (camera->film->xResolution * camera->film->yResolution) / (16*16));
        nDirectTasks = RoundUpPow2(nDirectTasks);
        ProgressReporter directProgress(nDirectTasks, "Direct Lighting");
        for (int i = 0; i < nDirectTasks; ++i)
            directTasks.push_back(new MLTDirectTask(scene, camera, this, &directLighting,
                &sampler, sample, i, nDirectTasks, directProgress));
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
    bool ignoreDirect = doDirectSeparately || indirectOnly;
    MLTSample sample(maxDepth);
    for (int i = 0; i < nBootstrap; ++i) {
        // Compute contribution for random sample for MLT bootstrapping
        LargeStep(rng, &sample, maxDepth, x0, x1, y0, y1, t0, t1);
        float contrib = I(L(scene, this, camera, arena, rng, maxDepth, ignoreDirect, sample), sample);
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
    if (scene->lights.size() > 0) {
        // Launch tasks to generate Metropolis samples
        int nTasks = int(nSamples / 50000);
        nTasks = max(nTasks, 32 * NumSystemCores());
        nTasks = min(nTasks, 32768);
        nSamples = (nSamples / nTasks) * nTasks;
        ProgressReporter progress(nTasks, "Metropolis");
        vector<Task *> tasks;
        for (int i = 0; i < nTasks; ++i)
            tasks.push_back(new MLTTask(progress, i, int(nSamples/nTasks), nSamples,
                nPixels, x0, x1, y0, y1, t0, t1, b, largeStepProbability, initialSample,
                ignoreDirect, maxConsecutiveRejects, maxDepth, scene, camera, this,
                &nSamplesFinished));
        EnqueueTasks(tasks);
        WaitForAllTasks();
        for (uint32_t i = 0; i < tasks.size(); ++i)
            delete tasks[i];
        progress.Done();
    }
    camera->film->WriteImage();
}


void MLTDirectIntegrator::RequestSamples(Sampler *sampler, Sample *sample,
        const Scene *scene) {
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


Spectrum MLTDirectIntegrator::Li(const Scene *scene,
        const Renderer *renderer, const RayDifferential &ray, const Intersection &isect,
        const Sample *sample, RNG &rng, MemoryArena &arena) const {
    Spectrum L(0.f);
    // Evaluate BSDF at hit point
    BSDF *bsdf = isect.GetBSDF(ray, arena);
    Vector wo = -ray.d;
    const Point &p = bsdf->dgShading.p;
    const Normal &n = bsdf->dgShading.nn;
    // Compute emitted light if ray hit an area light source
    L += isect.Le(wo);
    L += UniformSampleAllLights(scene, renderer, arena, p, n, wo,
            isect.rayEpsilon, ray.time, bsdf, sample, rng, lightSampleOffsets,
            bsdfSampleOffsets);
    return L;
}


void MLTDirectTask::Run() {
    Sampler *sampler = mainSampler->GetSubSampler(taskNum, taskCount);
    if (!sampler) {
        progress.Update();
        return;
    }
    MemoryArena arena;
    RNG rng(taskNum);
    // Allocate space for samples and intersections
    int maxSamples = sampler->MaximumSampleCount();
    Sample *samples = origSample->Duplicate(maxSamples);
    RayDifferential *rays = new RayDifferential[maxSamples];
    Spectrum *Ls = new Spectrum[maxSamples];
    Spectrum *Ts = new Spectrum[maxSamples];
    Intersection *isects = new Intersection[maxSamples];
    int sampleCount;
    while ((sampleCount = sampler->GetMoreSamples(samples, rng)) > 0) {
        for (int i = 0; i < sampleCount; ++i) {
            // Find camera ray for _sample[i]_
            PBRT_STARTED_GENERATING_CAMERA_RAY(&samples[i]);
            float rayWeight = camera->GenerateRayDifferential(samples[i], &rays[i]);
            PBRT_FINISHED_GENERATING_CAMERA_RAY(&samples[i], &rays[i], rayWeight);
            Ls[i] = 0.f;
            if (rayWeight > 0.f) {
                if (scene->Intersect(rays[i], &isects[i]))
                    Ls[i] = rayWeight * direct->Li(scene, renderer,
                        rays[i], isects[i], &samples[i], rng, arena);
                else {
                    for (uint32_t j = 0; j < scene->lights.size(); ++j)
                        Ls[i] += rayWeight * scene->lights[j]->Le(rays[i]);
                }
            }
            Ts[i] = 1.f;
            // Issue warning if unexpected radiance value returned
            if (Ls[i].HasNaNs()) {
                Error("Not-a-number radiance value returned "
                      "for image sample.  Setting to black.");
                Ls[i] = Spectrum(0.f);
            }
            else if (Ls[i].y() < -1e-5) {
                Error("Negative luminance value, %f, returned"
                      "for image sample.  Setting to black.", Ls[i].y());
                Ls[i] = Spectrum(0.f);
            }
            else if (isinf(Ls[i].y())) {
                Error("Infinite luminance value returned"
                      "for image sample.  Setting to black.");
                Ls[i] = Spectrum(0.f);
            }
            PBRT_FINISHED_CAMERA_RAY_INTEGRATION(&rays[i], &samples[i], &Ls[i]);
        }
        // Report sample results to _Sampler_, add contributions to image
        if (sampler->ReportResults(samples, rays, Ls, isects, sampleCount))
        {
            for (int i = 0; i < sampleCount; ++i)
            {
                PBRT_STARTED_ADDING_IMAGE_SAMPLE(&samples[i], &rays[i], &Ls[i], &Ts[i]);
                camera->film->AddSample(samples[i], Ls[i]);
                PBRT_FINISHED_ADDING_IMAGE_SAMPLE();
            }
        }

        // Free \use{MemoryArena} memory from computing image sample values
        arena.FreeAll();
    }
    camera->film->UpdateDisplay(sampler->xPixelStart,
        sampler->yPixelStart, sampler->xPixelEnd+1, sampler->yPixelEnd+1);
    delete sampler;
    delete[] samples;
    delete[] rays;
    delete[] Ls;
    delete[] Ts;
    delete[] isects;
    progress.Update();
    PBRT_FINISHED_RENDERTASK(taskNum);
}


MLTTask::MLTTask(ProgressReporter &prog, int tn, int ns, int64_t ts, int np,
    int xx0, int xx1, int yy0, int yy1, float tt0, float tt1,
    float bb, float lsp, const MLTSample &is, bool id, int mr,
    int md, const Scene *sc, const Camera *c, const Renderer *ren,
    volatile float *nfinished)
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
    sampleLs[currentSample] = L(scene, renderer, camera, arena, rng, maxDepth, ignoreDirect,
        mltSamples[currentSample]);
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
        sampleLs[proposedSample] = L(scene, renderer, camera, arena, rng, maxDepth, ignoreDirect,
            mltSamples[proposedSample]);
        float currentI = I(sampleLs[currentSample], mltSamples[currentSample]);
        float proposedI = I(sampleLs[proposedSample], mltSamples[proposedSample]);
        float a = min(1.f, proposedI / currentI);
        float currentWeight = (1.f - a) / (currentI / b + largeStepProbability) *
            float(nPixels) / float(totalSamples);
        float proposedWeight = (a + (largeStep ? 1.f : 0.f)) /
            (proposedI / b + largeStepProbability) * float(nPixels) / float(totalSamples);

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
    camera->film->UpdateDisplay(x0, y0, x1, y1, float(totalSamples)/nf);
    progress.Update();
    PBRT_MLT_FINISHED_MLT_TASK(this);
}


static Spectrum L(const Scene *scene, const Renderer *renderer,
        const Camera *camera, MemoryArena &arena,
        RNG &rng, int maxDepth, bool ignoreDirect, const MLTSample &sample) {
    // Generate camera ray from Metropolis sample
    RayDifferential ray;
    float cameraWeight = camera->GenerateRayDifferential(sample.cameraSample, &ray);
    Spectrum pathThroughput = cameraWeight, L = 0.;
    bool specularBounce = false;
    for (int pathLength = 0; pathLength < maxDepth; ++pathLength) {
        // Find next intersection in Metropolis light path
        Intersection isect;
        if (!scene->Intersect(ray, &isect)) {
            bool includeLe = (specularBounce && pathLength >= 1) ||
                             (!ignoreDirect && pathLength == 0);
            if (includeLe)
                for (uint32_t i = 0; i < scene->lights.size(); ++i)
                   L += pathThroughput * scene->lights[i]->Le(ray);
            break;
        }
        if (specularBounce)
            L += pathThroughput * isect.Le(-ray.d);
        BSDF *bsdf = isect.GetBSDF(ray, arena);
        const Point &p = bsdf->dgShading.p;
        const Normal &n = bsdf->dgShading.nn;
        Vector wo = -ray.d;
        const PathSample &ps = sample.pathSamples[pathLength];
        // Sample direct illumination for Metropolis path vertex
        if (!ignoreDirect || pathLength > 0) {
            LightSample lightSample(ps.lightDir0, ps.lightDir1, ps.lightNum0);
            BSDFSample bsdfSample(ps.bsdfLightDir0, ps.bsdfLightDir1, ps.bsdfLightComponent);
            uint32_t lightNum = Floor2Int(ps.lightNum1 * scene->lights.size());
            lightNum = min(lightNum, (uint32_t)(scene->lights.size()-1));
            const Light *light = scene->lights[lightNum];
            L += pathThroughput *
                 EstimateDirect(scene, renderer, arena, light, p, n, wo,
                     isect.rayEpsilon, sample.cameraSample.time, bsdf, rng,
                     lightSample, bsdfSample);
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
        pathThroughput *= f * AbsDot(wi, n) / pdf;
        ray = RayDifferential(p, wi, ray, isect.rayEpsilon);
        
        //pathThroughput *= renderer->Transmittance(scene, ray, NULL, rng, arena);
    }
    return L;
}


Spectrum MetropolisRenderer::Li(const Scene *scene, const RayDifferential &ray,
    const Sample *sample, RNG &rng, MemoryArena &arena, Intersection *isect,
    Spectrum *T) const {
Severe("WHA MLT::Li()");
return 0.f;
}


Spectrum MetropolisRenderer::Transmittance(const Scene *scene, const RayDifferential &ray,
    const Sample *sample, RNG &rng, MemoryArena &arena) const {
// FIXME
    return 1.f;
}


