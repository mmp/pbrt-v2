
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


