
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
