
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

#ifndef PBRT_RENDERERS_AGGREGATETEST_H
#define PBRT_RENDERERS_AGGREGATETEST_H

// renderers/aggregatetest.h*
#include "pbrt.h"
#include "renderer.h"
#include "memory.h"

// AggregateTest Declarations
class AggregateTest : public Renderer {
public:
    // AggregateTest Public Methods
    AggregateTest(int nIters, const vector<Reference<Primitive> > &primitives);
    void Render(const Scene *scene);
    Spectrum Li(const Scene *scene, const RayDifferential &ray,
        const Sample *sample, RNG &rng, MemoryArena &arena, Intersection *isect = NULL,
        Spectrum *T = NULL) const;
    Spectrum Transmittance(const Scene *scene, const RayDifferential &ray,
            const Sample *sample, RNG &rng, MemoryArena &arena) const;
private:
    // AggregateTest Private Data
    int nIterations;
    vector<Reference<Primitive> > primitives;
    vector<BBox> bboxes;
};


AggregateTest *CreateAggregateTestRenderer(const ParamSet &params,
    const vector<Reference<Primitive> > &primitives);

#endif // PBRT_RENDERERS_AGGREGATETEST_H
