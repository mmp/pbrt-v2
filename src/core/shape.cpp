
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


// core/shape.cpp*
#include "stdafx.h"
#include "shape.h"

// Shape Method Definitions
Shape::~Shape() {
}


Shape::Shape(const Transform *o2w, const Transform *w2o, bool ro)
    : ObjectToWorld(o2w), WorldToObject(w2o), ReverseOrientation(ro),
      TransformSwapsHandedness(o2w->SwapsHandedness()),
      shapeId(nextshapeId++) {
    // Update shape creation statistics
    PBRT_CREATED_SHAPE(this);
}


uint32_t Shape::nextshapeId = 1;
BBox Shape::WorldBound() const {
    return (*ObjectToWorld)(ObjectBound());
}


bool Shape::CanIntersect() const {
    return true;
}


void Shape::Refine(vector<Reference<Shape> > &refined) const {
    Severe("Unimplemented Shape::Refine() method called");
}


bool Shape::Intersect(const Ray &ray, float *tHit, float *rayEpsilon,
                      DifferentialGeometry *dg) const {
    Severe("Unimplemented Shape::Intersect() method called");
    return false;
}


bool Shape::IntersectP(const Ray &ray) const {
    Severe("Unimplemented Shape::IntersectP() method called");
    return false;
}


float Shape::Area() const {
    Severe("Unimplemented Shape::Area() method called");
    return 0.;
}


float Shape::Pdf(const Point &p, const Vector &wi) const {
    // Intersect sample ray with area light geometry
    DifferentialGeometry dgLight;
    Ray ray(p, wi, 1e-3f);
    ray.depth = -1; // temporary hack to ignore alpha mask
    float thit, rayEpsilon;
    if (!Intersect(ray, &thit, &rayEpsilon, &dgLight)) return 0.;

    // Convert light sample weight to solid angle measure
    float pdf = DistanceSquared(p, ray(thit)) /
                (AbsDot(dgLight.nn, -wi) * Area());
    if (isinf(pdf)) pdf = 0.f;
    return pdf;
}


