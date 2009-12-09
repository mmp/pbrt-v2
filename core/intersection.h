
/*
    pbrt source code Copyright(c) 1998-2009 Matt Pharr and Greg Humphreys.

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

#ifndef PBRT_CORE_INTERSECTION_H
#define PBRT_CORE_INTERSECTION_H

// core/intersection.h*
#include "pbrt.h"
#include "diffgeom.h"
#include "transform.h"

// Intersection Declarations
struct Intersection {
    // Intersection Public Methods
    Intersection() {
        primitive = NULL;
        shapeId = primitiveId = 0;
        rayEpsilon = 0.f;
    }
    BSDF *GetBSDF(const RayDifferential &ray, MemoryArena &arena) const;
    BSSRDF *GetBSSRDF(const RayDifferential &ray, MemoryArena &arena) const;
    Spectrum Le(const Vector &wo) const;

    // Intersection Public Data
    DifferentialGeometry dg;
    const Primitive *primitive;
    Transform WorldToObject, ObjectToWorld;
    uint32_t shapeId, primitiveId;
    float rayEpsilon;
};



#endif // PBRT_CORE_INTERSECTION_H
