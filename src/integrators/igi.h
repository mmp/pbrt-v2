
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

#ifndef PBRT_INTEGRATORS_IGI_H
#define PBRT_INTEGRATORS_IGI_H

// integrators/igi.h*
#include "pbrt.h"
#include "integrator.h"

// IGIIntegrator Local Structures
struct VirtualLight {
    VirtualLight() { }
    VirtualLight(const Point &pp, const Normal &nn, const Spectrum &c,
                 float reps)
        : p(pp), n(nn), pathContrib(c), rayEpsilon(reps) { }
    Point p;
    Normal n;
    Spectrum pathContrib;
    float rayEpsilon;
};



// IGIIntegrator Declarations
class IGIIntegrator : public SurfaceIntegrator {
public:
    // IGIIntegrator Public Methods
    ~IGIIntegrator();
    Spectrum Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential &ray, const Intersection &isect,
        const Sample *sample, RNG &rng, MemoryArena &arena) const;
    void RequestSamples(Sampler *sampler, Sample *sample, const Scene *scene);
    void Preprocess(const Scene *, const Camera *, const Renderer *);
    IGIIntegrator(uint32_t nl, uint32_t ns, float rrt, int maxd, float gl, int ng) {
        nLightPaths = RoundUpPow2(nl);
        nLightSets = RoundUpPow2(ns);
        rrThreshold = rrt;
        maxSpecularDepth = maxd;
        virtualLights.resize(nLightSets);
        gLimit = gl;
        nGatherSamples = ng;
        lightSampleOffsets = NULL;
        bsdfSampleOffsets = NULL;
    }
private:
    // IGIIntegrator Private Data

    // Declare sample parameters for light source sampling
    LightSampleOffsets *lightSampleOffsets;
    BSDFSampleOffsets *bsdfSampleOffsets;
    uint32_t nLightPaths, nLightSets;
    float gLimit;
    int nGatherSamples;
    float rrThreshold;
    int maxSpecularDepth;
    int vlSetOffset;
    BSDFSampleOffsets gatherSampleOffset;
    vector<vector<VirtualLight> > virtualLights;
};


IGIIntegrator *CreateIGISurfaceIntegrator(const ParamSet &params);

#endif // PBRT_INTEGRATORS_IGI_H
