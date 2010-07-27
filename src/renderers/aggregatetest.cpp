
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


// renderers/aggregatetest.cpp*
#include "stdafx.h"
#include "renderers/aggregatetest.h"
#include "progressreporter.h"
#include "scene.h"
#include "paramset.h"
#include "montecarlo.h"
#include "primitive.h"
#include "intersection.h"

// AggregateTest Method Definitions
AggregateTest::AggregateTest(int niters,
        const vector<Reference<Primitive> > &p) {
    nIterations = niters;
    for (uint32_t i = 0; i < p.size(); ++i)
        p[i]->FullyRefine(primitives);
    for (uint32_t i = 0; i < primitives.size(); ++i)
        bboxes.push_back(primitives[i]->WorldBound());
}


AggregateTest *CreateAggregateTestRenderer(const ParamSet &params,
    const vector<Reference<Primitive> > &primitives) {
    int niters = params.FindOneInt("niters", 100000);
    return new AggregateTest(niters, primitives);
}


void AggregateTest::Render(const Scene *scene) {
    RNG rng;
    ProgressReporter prog(nIterations, "Aggregate Test");
    // Compute bounding box of region used to generate random rays
    BBox bbox = scene->WorldBound();
    bbox.Expand(bbox.pMax[bbox.MaximumExtent()] -
                bbox.pMin[bbox.MaximumExtent()]);
    Point lastHit;
    float lastEps = 0.f;
    for (int i = 0; i < nIterations; ++i) {
        // Choose random rays, _rayAccel_ and _rayAll_ for testing

        // Choose ray origin for testing accelerator
        Point org(Lerp(rng.RandomFloat(), bbox.pMin.x, bbox.pMax.x),
                  Lerp(rng.RandomFloat(), bbox.pMin.y, bbox.pMax.y),
                  Lerp(rng.RandomFloat(), bbox.pMin.z, bbox.pMax.z));
        if ((rng.RandomUInt() % 4) == 0) org = lastHit;

        // Choose ray direction for testing accelerator
        Vector dir = UniformSampleSphere(rng.RandomFloat(), rng.RandomFloat());
        if ((rng.RandomUInt() % 32) == 0)      dir.x = dir.y = 0.f;
        else if ((rng.RandomUInt() % 32) == 0) dir.x = dir.z = 0.f;
        else if ((rng.RandomUInt() % 32) == 0) dir.y = dir.z = 0.f;

        // Choose ray epsilon for testing accelerator
        float eps = 0.f;
        if (rng.RandomFloat() < .25) eps = lastEps;
        else if (rng.RandomFloat() < .25) eps = 1e-3f;
        Ray rayAccel(org, dir, eps);
        Ray rayAll = rayAccel;

        // Compute intersections using accelerator and exhaustive testing
        Intersection isectAccel, isectAll;
        bool hitAccel = scene->Intersect(rayAccel, &isectAccel);
        bool hitAll = false;
        bool inconsistentBounds = false;
        for (uint32_t j = 0; j < primitives.size(); ++j) {
            if (bboxes[j].IntersectP(rayAll))
                hitAll |= primitives[j]->Intersect(rayAll, &isectAll);
            else if (primitives[j]->Intersect(rayAll, &isectAll))
                inconsistentBounds = true;
        }

        // Report any inconsistencies between intersections
        if (!inconsistentBounds &&
            ((hitAccel != hitAll) || (rayAccel.maxt != rayAll.maxt)))
            Warning("Disagreement: t accel %.16g [%a] t exhaustive %.16g [%a]\n"
                    "Ray: org [%a, %a, %a], dir [%a, %a, %a], mint = %a",
                    rayAccel.maxt, rayAll.maxt, rayAccel.maxt, rayAll.maxt,
                    rayAll.o.x, rayAll.o.y, rayAll.o.z,
                    rayAll.d.x, rayAll.d.y, rayAll.d.z, rayAll.mint);
        if (hitAll) {
            lastHit = rayAll(rayAll.maxt);
            lastEps = isectAll.rayEpsilon;
        }
        prog.Update();
    }
    prog.Done();
}


Spectrum AggregateTest::Li(const Scene *scene, const RayDifferential &ray,
        const Sample *sample, RNG &rng, MemoryArena &arena, Intersection *isect,
        Spectrum *T) const {
    return 0.f;
}


Spectrum AggregateTest::Transmittance(const Scene *scene,
        const RayDifferential &ray, const Sample *sample, RNG &rng,
        MemoryArena &arena) const {
    return 0.f;
}


