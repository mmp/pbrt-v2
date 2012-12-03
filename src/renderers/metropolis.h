
/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_RENDERERS_METROPOLIS_H
#define PBRT_RENDERERS_METROPOLIS_H

// renderers/metropolis.h*
#include "pbrt.h"
#include "renderer.h"
#include "parallel.h"
struct MLTSample;
class DirectLightingIntegrator;
struct LightingSample;

// Metropolis Declarations
struct PathVertex;

class MetropolisRenderer : public Renderer {
public:
    // MetropolisRenderer Public Methods
    MetropolisRenderer(int perPixelSamples, int nBootstrap,
        int directPixelSamples, float largeStepProbability,
        bool doDirectSeparately, int maxConsecutiveRejects, int maxDepth,
        Camera *camera, bool doBidirectional);
    ~MetropolisRenderer();
    void Render(const Scene *scene);
    Spectrum Li(const Scene *scene, const RayDifferential &ray,
        const Sample *sample, RNG &rng, MemoryArena &arena,
        Intersection *isect = NULL, Spectrum *T = NULL) const;
    Spectrum Transmittance(const Scene *scene, const RayDifferential &ray,
        const Sample *sample, RNG &rng, MemoryArena &arena) const;
private:
    // MetropolisRenderer Private Methods
    Spectrum PathL(const MLTSample &sample, const Scene *scene,
        MemoryArena &arena, const Camera *camera,
        const Distribution1D *lightDistribution, PathVertex *cameraPath,
        PathVertex *lightPath, RNG &rng) const;
    Spectrum Lpath(const Scene *scene, const PathVertex *path, int pathLength,
        MemoryArena &arena, const vector<LightingSample> &samples,
        RNG &rng, float time, const Distribution1D *lightDistribution,
        const RayDifferential &escapedRay, const Spectrum &escapedAlpha) const;
    Spectrum Lbidir(const Scene *scene,
        const PathVertex *cameraPath, int cameraPathLength,
        const PathVertex *lightPath, int lightPathLength,
        MemoryArena &arena, const vector<LightingSample> &samples,
        RNG &rng, float time, const Distribution1D *lightDistribution,
        const RayDifferential &escapedRay, const Spectrum &escapedAlpha) const;

    // MetropolisRenderer Private Data
    Camera *camera;
    bool bidirectional;
    uint32_t nDirectPixelSamples, nPixelSamples, maxDepth;
    uint32_t largeStepsPerPixel, nBootstrap, maxConsecutiveRejects;
    DirectLightingIntegrator *directLighting;
    AtomicInt32 nTasksFinished;
    friend class MLTTask;
};


MetropolisRenderer *CreateMetropolisRenderer(const ParamSet &params,
    Camera *camera);

#endif // PBRT_RENDERERS_METROPOLIS_H
