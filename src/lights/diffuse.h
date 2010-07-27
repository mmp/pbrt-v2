
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

#ifndef PBRT_LIGHTS_DIFFUSE_H
#define PBRT_LIGHTS_DIFFUSE_H

// lights/diffuse.h*
#include "pbrt.h"
#include "light.h"
#include "primitive.h"

// DiffuseAreaLight Declarations
class DiffuseAreaLight : public AreaLight {
public:
    // DiffuseAreaLight Public Methods
    DiffuseAreaLight(const Transform &light2world,
        const Spectrum &Le, int ns, const Reference<Shape> &shape);
    ~DiffuseAreaLight();
    Spectrum L(const Point &p, const Normal &n, const Vector &w) const {
        return Dot(n, w) > 0.f ? Lemit : 0.f;
    }
    Spectrum Power(const Scene *) const;
    bool IsDeltaLight() const { return false; }
    float Pdf(const Point &, const Vector &) const;
    Spectrum Sample_L(const Point &P, float pEpsilon, const LightSample &ls, float time,
        Vector *wo, float *pdf, VisibilityTester *visibility) const;
    Spectrum Sample_L(const Scene *scene, const LightSample &ls, float u1, float u2,
        float time, Ray *ray, Normal *Ns, float *pdf) const;
protected:
    // DiffuseAreaLight Protected Data
    Spectrum Lemit;
    ShapeSet *shapeSet;
    float area;
};


AreaLight *CreateDiffuseAreaLight(const Transform &light2world, const ParamSet &paramSet,
        const Reference<Shape> &shape);

#endif // PBRT_LIGHTS_DIFFUSE_H
