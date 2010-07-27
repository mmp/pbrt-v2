
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

#ifndef PBRT_SHAPES_PARABOLOID_H
#define PBRT_SHAPES_PARABOLOID_H

// shapes/paraboloid.h*
#include "shape.h"

// Paraboloid Declarations
class Paraboloid : public Shape {
public:
    // Paraboloid Public Methods
    Paraboloid(const Transform *o2w, const Transform *w2o, bool ro, float rad,
               float z0, float z1, float tm );
    BBox ObjectBound() const;
    bool Intersect(const Ray &ray, float *tHit, float *rayEpsilon,
                      DifferentialGeometry *dg) const;
    bool IntersectP(const Ray &ray) const;
    float Area() const;
protected:
    // Paraboloid Private Data
    float radius;
    float zmin, zmax;
    float phiMax;
};


Paraboloid *CreateParaboloidShape(const Transform *o2w, const Transform *w2o,
        bool reverseOrientation, const ParamSet &params);

#endif // PBRT_SHAPES_PARABOLOID_H
