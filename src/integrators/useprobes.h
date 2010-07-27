
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
