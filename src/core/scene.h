
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

#ifndef PBRT_CORE_SCENE_H
#define PBRT_CORE_SCENE_H

// core/scene.h*
#include "pbrt.h"
#include "primitive.h"
#include "integrator.h"

// Scene Declarations
class Scene {
public:
    // Scene Public Methods
    Scene(Primitive *accel, const vector<Light *> &lts, VolumeRegion *vr);
    ~Scene();
    bool Intersect(const Ray &ray, Intersection *isect) const {
        PBRT_STARTED_RAY_INTERSECTION(const_cast<Ray *>(&ray));
        bool hit = aggregate->Intersect(ray, isect);
        PBRT_FINISHED_RAY_INTERSECTION(const_cast<Ray *>(&ray), isect, int(hit));
        return hit;
    }
    bool IntersectP(const Ray &ray) const {
        PBRT_STARTED_RAY_INTERSECTIONP(const_cast<Ray *>(&ray));
        bool hit = aggregate->IntersectP(ray);
        PBRT_FINISHED_RAY_INTERSECTIONP(const_cast<Ray *>(&ray), int(hit));
        return hit;
    }
    const BBox &WorldBound() const;

    // Scene Public Data
    Primitive *aggregate;
    vector<Light *> lights;
    VolumeRegion *volumeRegion;
    BBox bound;
};



#endif // PBRT_CORE_SCENE_H
