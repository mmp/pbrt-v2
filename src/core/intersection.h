
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
