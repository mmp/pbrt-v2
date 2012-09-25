
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


// core/intersection.cpp*
#include "stdafx.h"
#include "intersection.h"
#include "shape.h"
#include "primitive.h"
#include "light.h"

// Intersection Method Definitions
BSDF *Intersection::GetBSDF(const RayDifferential &ray,
                            MemoryArena &arena) const {
    PBRT_STARTED_BSDF_SHADING(const_cast<RayDifferential *>(&ray));
    dg.ComputeDifferentials(ray);
    BSDF *bsdf = primitive->GetBSDF(dg, ObjectToWorld, arena);
    PBRT_FINISHED_BSDF_SHADING(const_cast<RayDifferential *>(&ray), bsdf);
    return bsdf;
}


BSSRDF *Intersection::GetBSSRDF(const RayDifferential &ray,
          MemoryArena &arena) const {
    PBRT_STARTED_BSSRDF_SHADING(const_cast<RayDifferential *>(&ray));
    dg.ComputeDifferentials(ray);
    BSSRDF *bssrdf = primitive->GetBSSRDF(dg, ObjectToWorld, arena);
    PBRT_FINISHED_BSSRDF_SHADING(const_cast<RayDifferential *>(&ray), bssrdf);
    return bssrdf;
}


Spectrum Intersection::Le(const Vector &w) const {
    const AreaLight *area = primitive->GetAreaLight();
    return area ? area->L(dg.p, dg.nn, w) : Spectrum(0.);
}


