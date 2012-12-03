
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

#ifndef PBRT_INTEGRATORS_USEPROBES_H
#define PBRT_INTEGRATORS_USEPROBES_H

// integrators/useprobes.h*
#include "pbrt.h"
#include "integrator.h"
#include "sh.h"

// UseRadianceProbes Declarations
class UseRadianceProbes : public SurfaceIntegrator {
public:
    // UseRadianceProbes Public Methods
    UseRadianceProbes(const string &filename);
    ~UseRadianceProbes();
    void RequestSamples(Sampler *sampler, Sample *sample, const Scene *scene);
    Spectrum Li(const Scene *scene, const Renderer *,
                const RayDifferential &ray, const Intersection &isect,
                const Sample *sample, RNG &rng, MemoryArena &arena) const;
private:
    // UseRadianceProbes Private Methods
    const Spectrum *c_inXYZ(int lmax, int vx, int vy, int vz) const {
        vx = Clamp(vx, 0, nProbes[0]-1);
        vy = Clamp(vy, 0, nProbes[1]-1);
        vz = Clamp(vz, 0, nProbes[2]-1);
        int offset = vx + vy * nProbes[0] + vz * nProbes[0] * nProbes[1];
        return &c_in[SHTerms(lmax) * offset];
    }

    // UseRadianceProbes Private Data
    BBox bbox;
    int lmax, includeDirectInProbes, includeIndirectInProbes;
    int nProbes[3];
    Spectrum *c_in;

    // Declare sample parameters for light source sampling
    LightSampleOffsets *lightSampleOffsets;
    BSDFSampleOffsets *bsdfSampleOffsets;
};


extern UseRadianceProbes *CreateRadianceProbesSurfaceIntegrator(const ParamSet &paramSet);

#endif // PBRT_INTEGRATORS_USEPROBES_H
