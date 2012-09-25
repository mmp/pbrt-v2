
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


// core/primitive.cpp*
#include "stdafx.h"
#include "primitive.h"
#include "light.h"
#include "intersection.h"

// Primitive Method Definitions
uint32_t Primitive::nextprimitiveId = 1;
Primitive::~Primitive() { }

bool Primitive::CanIntersect() const {
    return true;
}



void Primitive::Refine(vector<Reference<Primitive> > &refined) const {
    Severe("Unimplemented Primitive::Refine() method called!");
}


void
Primitive::FullyRefine(vector<Reference<Primitive> > &refined) const {
    vector<Reference<Primitive> > todo;
    todo.push_back(const_cast<Primitive *>(this));
    while (todo.size()) {
        // Refine last primitive in todo list
        Reference<Primitive> prim = todo.back();
        todo.pop_back();
        if (prim->CanIntersect())
            refined.push_back(prim);
        else
            prim->Refine(todo);
    }
}


const AreaLight *Aggregate::GetAreaLight() const {
    Severe("Aggregate::GetAreaLight() method"
         "called; should have gone to GeometricPrimitive");
    return NULL;
}


BSDF *Aggregate::GetBSDF(const DifferentialGeometry &,
        const Transform &, MemoryArena &) const {
    Severe("Aggregate::GetBSDF() method"
        "called; should have gone to GeometricPrimitive");
    return NULL;
}


BSSRDF *Aggregate::GetBSSRDF(const DifferentialGeometry &,
        const Transform &, MemoryArena &) const {
    Severe("Aggregate::GetBSSRDF() method"
        "called; should have gone to GeometricPrimitive");
    return NULL;
}



// TransformedPrimitive Method Definitions
bool TransformedPrimitive::Intersect(const Ray &r,
                                     Intersection *isect) const {
    Transform w2p;
    WorldToPrimitive.Interpolate(r.time, &w2p);
    Ray ray = w2p(r);
    if (!primitive->Intersect(ray, isect))
        return false;
    r.maxt = ray.maxt;
    isect->primitiveId = primitiveId;
    if (!w2p.IsIdentity()) {
        // Compute world-to-object transformation for instance
        isect->WorldToObject = isect->WorldToObject * w2p;
        isect->ObjectToWorld = Inverse(isect->WorldToObject);

        // Transform instance's differential geometry to world space
        Transform PrimitiveToWorld = Inverse(w2p);
        isect->dg.p = PrimitiveToWorld(isect->dg.p);
        isect->dg.nn = Normalize(PrimitiveToWorld(isect->dg.nn));
        isect->dg.dpdu = PrimitiveToWorld(isect->dg.dpdu);
        isect->dg.dpdv = PrimitiveToWorld(isect->dg.dpdv);
        isect->dg.dndu = PrimitiveToWorld(isect->dg.dndu);
        isect->dg.dndv = PrimitiveToWorld(isect->dg.dndv);
    }
    return true;
}


bool TransformedPrimitive::IntersectP(const Ray &r) const {
    return primitive->IntersectP(WorldToPrimitive(r));
}



// GeometricPrimitive Method Definitions
BBox GeometricPrimitive::WorldBound() const {
    return shape->WorldBound();
}


bool GeometricPrimitive::IntersectP(const Ray &r) const {
    return shape->IntersectP(r);
}


bool GeometricPrimitive::CanIntersect() const {
    return shape->CanIntersect();
}


void GeometricPrimitive::
        Refine(vector<Reference<Primitive> > &refined)
        const {
    vector<Reference<Shape> > r;
    shape->Refine(r);
    for (uint32_t i = 0; i < r.size(); ++i) {
        GeometricPrimitive *gp = new GeometricPrimitive(r[i],
               material, areaLight);
        refined.push_back(gp);
    }
}


GeometricPrimitive::GeometricPrimitive(const Reference<Shape> &s,
        const Reference<Material> &m, AreaLight *a)
    : shape(s), material(m), areaLight(a) {
}


bool GeometricPrimitive::Intersect(const Ray &r,
                                   Intersection *isect) const {
    float thit, rayEpsilon;
    if (!shape->Intersect(r, &thit, &rayEpsilon, &isect->dg))
        return false;
    isect->primitive = this;
    isect->WorldToObject = *shape->WorldToObject;
    isect->ObjectToWorld = *shape->ObjectToWorld;
    isect->shapeId = shape->shapeId;
    isect->primitiveId = primitiveId;
    isect->rayEpsilon = rayEpsilon;
    r.maxt = thit;
    return true;
}


const AreaLight *GeometricPrimitive::GetAreaLight() const {
    return areaLight;
}


BSDF *GeometricPrimitive::GetBSDF(const DifferentialGeometry &dg,
                                  const Transform &ObjectToWorld,
                                  MemoryArena &arena) const {
    DifferentialGeometry dgs;
    shape->GetShadingGeometry(ObjectToWorld, dg, &dgs);
    return material->GetBSDF(dg, dgs, arena);
}


BSSRDF *GeometricPrimitive::GetBSSRDF(const DifferentialGeometry &dg,
                                  const Transform &ObjectToWorld,
                                  MemoryArena &arena) const {
    DifferentialGeometry dgs;
    shape->GetShadingGeometry(ObjectToWorld, dg, &dgs);
    return material->GetBSSRDF(dg, dgs, arena);
}


