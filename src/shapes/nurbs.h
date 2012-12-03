
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

#ifndef PBRT_SHAPES_NURBS_H
#define PBRT_SHAPES_NURBS_H

// shapes/nurbs.h*
#include "pbrt.h"
#include "shape.h"
#include "geometry.h"


// NURBS Declarations
class NURBS : public Shape {
public:
    // NURBS Methods
    NURBS(const Transform *o2w, const Transform *w2o,
        bool ReverseOrientation, int nu, int uorder,
        const float *uknot, float umin, float umax,
        int nv, int vorder, const float *vknot, float vmin, float vmax,
        const float *P, bool isHomogeneous);
    ~NURBS();
    BBox ObjectBound() const;
    BBox WorldBound() const;
    bool CanIntersect() const { return false; }
    void Refine(vector<Reference<Shape> > &refined) const;
private:
    // NURBS Data
    int nu, uorder, nv, vorder;
    float umin, umax, vmin, vmax;
    float *uknot, *vknot;
    bool isHomogeneous;
    float *P;
};




extern NURBS *CreateNURBSShape(const Transform *o2w, const Transform *w2o,
    bool ReverseOrientation, const ParamSet &params);


#endif // PBRT_SHAPES_NURBS_H
