
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
