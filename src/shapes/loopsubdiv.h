
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

#ifndef PBRT_SHAPES_LOOPSUBDIV_H
#define PBRT_SHAPES_LOOPSUBDIV_H

// shapes/loopsubdiv.h*
#include "shape.h"
struct SDVertex;
struct SDFace;

// LoopSubdiv Declarations
class LoopSubdiv : public Shape {
public:
    // LoopSubdiv Public Methods
    LoopSubdiv(const Transform *o2w, const Transform *w2o, bool ro,
               int nt, int nv, const int *vi,
               const Point *P, int nlevels);
    ~LoopSubdiv();
    bool CanIntersect() const;
    void Refine(vector<Reference<Shape> > &refined) const;
    BBox ObjectBound() const;
    BBox WorldBound() const;
private:
    // LoopSubdiv Private Methods
    static float beta(int valence) {
        if (valence == 3) return 3.f/16.f;
        else return 3.f / (8.f * valence);
    }
    static Point weightOneRing(SDVertex *vert, float beta);
    static Point weightBoundary(SDVertex *vert, float beta);
    static float gamma(int valence) {
        return 1.f / (valence + 3.f / (8.f * beta(valence)));
    }

    // LoopSubdiv Private Data
    int nLevels;
    vector<SDVertex *> vertices;
    vector<SDFace *> faces;
};


LoopSubdiv *CreateLoopSubdivShape(const Transform *o2w, const Transform *w2o,
        bool reverseOrientation, const ParamSet &params);

#endif // PBRT_SHAPES_LOOPSUBDIV_H
