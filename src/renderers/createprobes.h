
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

#ifndef PBRT_RENDERERS_CREATEPROBES_H
#define PBRT_RENDERERS_CREATEPROBES_H

// renderers/createprobes.h*
#include "pbrt.h"
#include "renderer.h"
#include "geometry.h"

// CreateRadianceProbes Declarations
class CreateRadianceProbes : public Renderer {
public:
    // CreateRadianceProbes Public Methods
    CreateRadianceProbes(SurfaceIntegrator *surf, VolumeIntegrator *vol,
        const Camera *camera, int lmax, float probeSpacing, const BBox &bbox,
        int nIndirSamples, bool includeDirect, bool includeIndirect,
        float time, const string &filename);
    ~CreateRadianceProbes();
    void Render(const Scene *scene);
    Spectrum Li(const Scene *scene, const RayDifferential &ray,
        const Sample *sample, RNG &rng, MemoryArena &arena, Intersection *isect,
        Spectrum *T) const;
    Spectrum Transmittance(const Scene *scene, const RayDifferential &ray,
        const Sample *sample, RNG &rng, MemoryArena &arena) const;
private:
    // CreateRadianceProbes Private Data
    SurfaceIntegrator *surfaceIntegrator;
    VolumeIntegrator *volumeIntegrator;
    const Camera *camera;
    int lmax, nIndirSamples;
    BBox bbox;
    bool includeDirectInProbes, includeIndirectInProbes;
    float time, probeSpacing;
    string filename;
};


CreateRadianceProbes *CreateRadianceProbesRenderer(const Camera *camera,
    SurfaceIntegrator *surf, VolumeIntegrator *vol, const ParamSet &params);

#endif // PBRT_RENDERERS_CREATEPROBES_H
