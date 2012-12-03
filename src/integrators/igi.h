
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
