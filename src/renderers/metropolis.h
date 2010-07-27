
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
