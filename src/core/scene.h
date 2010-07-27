
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
