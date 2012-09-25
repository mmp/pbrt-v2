
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


// shapes/disk.cpp*
#include "stdafx.h"
#include "shapes/disk.h"
#include "paramset.h"
#include "montecarlo.h"

// Disk Method Definitions
Disk::Disk(const Transform *o2w, const Transform *w2o, bool ro,
           float ht, float r, float ri, float tmax)
    : Shape(o2w, w2o, ro) {
    height = ht;
    radius = r;
    innerRadius = ri;
    phiMax = Radians(Clamp(tmax, 0.0f, 360.0f));
}


BBox Disk::ObjectBound() const {
    return BBox(Point(-radius, -radius, height),
                Point( radius,  radius, height));
}


bool Disk::Intersect(const Ray &r, float *tHit, float *rayEpsilon,
                     DifferentialGeometry *dg) const {
    // Transform _Ray_ to object space
    Ray ray;
    (*WorldToObject)(r, &ray);

    // Compute plane intersection for disk
    if (fabsf(ray.d.z) < 1e-7) return false;
    float thit = (height - ray.o.z) / ray.d.z;
    if (thit < ray.mint || thit > ray.maxt)
        return false;

    // See if hit point is inside disk radii and $\phimax$
    Point phit = ray(thit);
    float dist2 = phit.x * phit.x + phit.y * phit.y;
    if (dist2 > radius * radius || dist2 < innerRadius * innerRadius)
        return false;

    // Test disk $\phi$ value against $\phimax$
    float phi = atan2f(phit.y, phit.x);
    if (phi < 0) phi += 2. * M_PI;
    if (phi > phiMax)
        return false;

    // Find parametric representation of disk hit
    float u = phi / phiMax;
    float oneMinusV = ((sqrtf(dist2)-innerRadius) /
                       (radius-innerRadius));
    float invOneMinusV = (oneMinusV > 0.f) ? (1.f / oneMinusV) : 0.f;
    float v = 1.f - oneMinusV;
    Vector dpdu(-phiMax * phit.y, phiMax * phit.x, 0.);
    Vector dpdv(-phit.x * invOneMinusV, -phit.y * invOneMinusV, 0.);
    dpdu *= phiMax * INV_TWOPI;
    dpdv *= (radius - innerRadius) / radius;
    Normal dndu(0,0,0), dndv(0,0,0);

    // Initialize _DifferentialGeometry_ from parametric information
    const Transform &o2w = *ObjectToWorld;
    *dg = DifferentialGeometry(o2w(phit), o2w(dpdu), o2w(dpdv),
                               o2w(dndu), o2w(dndv), u, v, this);

    // Update _tHit_ for quadric intersection
    *tHit = thit;

    // Compute _rayEpsilon_ for quadric intersection
    *rayEpsilon = 5e-4f * *tHit;
    return true;
}


bool Disk::IntersectP(const Ray &r) const {
    // Transform _Ray_ to object space
    Ray ray;
    (*WorldToObject)(r, &ray);

    // Compute plane intersection for disk
    if (fabsf(ray.d.z) < 1e-7) return false;
    float thit = (height - ray.o.z) / ray.d.z;
    if (thit < ray.mint || thit > ray.maxt)
        return false;

    // See if hit point is inside disk radii and $\phimax$
    Point phit = ray(thit);
    float dist2 = phit.x * phit.x + phit.y * phit.y;
    if (dist2 > radius * radius || dist2 < innerRadius * innerRadius)
        return false;

    // Test disk $\phi$ value against $\phimax$
    float phi = atan2f(phit.y, phit.x);
    if (phi < 0) phi += 2. * M_PI;
    if (phi > phiMax)
        return false;
    return true;
}


float Disk::Area() const {
    return phiMax * 0.5f *
       (radius * radius - innerRadius * innerRadius);
}


Disk *CreateDiskShape(const Transform *o2w, const Transform *w2o,
        bool reverseOrientation, const ParamSet &params) {
    float height = params.FindOneFloat("height", 0.);
    float radius = params.FindOneFloat("radius", 1);
    float inner_radius = params.FindOneFloat("innerradius", 0);
    float phimax = params.FindOneFloat("phimax", 360);
    return new Disk(o2w, w2o, reverseOrientation, height, radius, inner_radius, phimax);
}


Point Disk::Sample(float u1, float u2, Normal *Ns) const {
    Point p;
    ConcentricSampleDisk(u1, u2, &p.x, &p.y);
    p.x *= radius;
    p.y *= radius;
    p.z = height;
    *Ns = Normalize((*ObjectToWorld)(Normal(0,0,1)));
    if (ReverseOrientation) *Ns *= -1.f;
    return (*ObjectToWorld)(p);
}


