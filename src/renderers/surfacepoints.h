
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
