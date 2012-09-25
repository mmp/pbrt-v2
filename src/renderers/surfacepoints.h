
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

#ifndef PBRT_RENDERERS_SURFACEPOINTS_H
#define PBRT_RENDERERS_SURFACEPOINTS_H

// renderers/surfacepoints.h*
#include "pbrt.h"
#include "geometry.h"
#include "renderer.h"

// SurfacePointsRenderer Declarations
struct SurfacePoint {
    SurfacePoint() { }
    SurfacePoint(const Point &pp, const Normal &nn, float a, float eps)
        : p(pp), n(nn), area(a), rayEpsilon(eps) { }
    // SurfacePoint Data
    Point p;
    Normal n;
    float area, rayEpsilon;
};


class SurfacePointsRenderer : public Renderer {
public:
    // SurfacePointsRenderer Public Methods
    SurfacePointsRenderer(float md, const Point &pc, float t,
                          const string &fn)
        : minDist(md), time(t), pCamera(pc), filename(fn) { }
    void Render(const Scene *scene);
    Spectrum Li(const Scene *scene, const RayDifferential &ray,
        const Sample *sample, RNG &rng, MemoryArena &arena,
        Intersection *isect, Spectrum *T) const;
    Spectrum Transmittance(const Scene *scene, const RayDifferential &ray,
        const Sample *sample, RNG &rng, MemoryArena &arena) const;
private:
    // SurfacePointsRenderer Private Data
    float minDist, time;
    Point pCamera;
    string filename;
    friend void FindPoissonPointDistribution(const Point &pCamera, float time,
        float minDist, const Scene *scene, vector<SurfacePoint> *points);
    vector<SurfacePoint> points;
};


void FindPoissonPointDistribution(const Point &pCamera, float time, float minDist,
    const Scene *scene, vector<SurfacePoint> *points);
SurfacePointsRenderer *CreateSurfacePointsRenderer(const ParamSet &params,
    const Point &pCamera, float time);

#endif // PBRT_RENDERERS_SURFACEPOINTS_H
