
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

#ifndef PBRT_INTEGRATORS_EMISSION_H
#define PBRT_INTEGRATORS_EMISSION_H

// integrators/emission.h*
#include "volume.h"
#include "integrator.h"
#include "scene.h"

// EmissionIntegrator Declarations
class EmissionIntegrator : public VolumeIntegrator {
public:
    // EmissionIntegrator Public Methods
    EmissionIntegrator(float ss) { stepSize = ss; }
    void RequestSamples(Sampler *sampler, Sample *sample, const Scene *scene);
    Spectrum Li(const Scene *scene, const Renderer *renderer,
            const RayDifferential &ray, const Sample *sample, RNG &rng,
            Spectrum *transmittance, MemoryArena &arena) const;
    Spectrum Transmittance(const Scene *scene, const Renderer *,
            const RayDifferential &ray, const Sample *sample, RNG &rng,
            MemoryArena &arena) const;
private:
    // EmissionIntegrator Private Data
    float stepSize;
    int tauSampleOffset, scatterSampleOffset;
};


EmissionIntegrator *CreateEmissionVolumeIntegrator(const ParamSet &params);

#endif // PBRT_INTEGRATORS_EMISSION_H
