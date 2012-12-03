
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

#ifndef PBRT_CORE_SHAPE_H
#define PBRT_CORE_SHAPE_H

// core/shape.h*
#include "pbrt.h"
#include "geometry.h"
#include "transform.h"
#include "diffgeom.h"
#include "memory.h"

// Shape Declarations
class Shape : public ReferenceCounted {
public:
    // Shape Interface
    Shape(const Transform *o2w, const Transform *w2o, bool ro);
    virtual ~Shape();
    virtual BBox ObjectBound() const = 0;
    virtual BBox WorldBound() const;
    virtual bool CanIntersect() const;
    virtual void Refine(vector<Reference<Shape> > &refined) const;
    virtual bool Intersect(const Ray &ray, float *tHit,
                           float *rayEpsilon, DifferentialGeometry *dg) const;
    virtual bool IntersectP(const Ray &ray) const;
    virtual void GetShadingGeometry(const Transform &obj2world,
            const DifferentialGeometry &dg,
            DifferentialGeometry *dgShading) const {
        *dgShading = dg;
    }
    virtual float Area() const;
    virtual Point Sample(float u1, float u2, Normal *Ns) const {
        Severe("Unimplemented Shape::Sample() method called");
        return Point();
    }
    virtual float Pdf(const Point &Pshape) const {
        return 1.f / Area();
    }
    virtual Point Sample(const Point &P, float u1, float u2,
                         Normal *Ns) const {
        return Sample(u1, u2, Ns);
    }
    virtual float Pdf(const Point &p, const Vector &wi) const;

    // Shape Public Data
    const Transform *ObjectToWorld, *WorldToObject;
    const bool ReverseOrientation, TransformSwapsHandedness;
    const uint32_t shapeId;
    static uint32_t nextshapeId;
};



#endif // PBRT_CORE_SHAPE_H
