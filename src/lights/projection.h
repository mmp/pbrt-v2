
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

#ifndef PBRT_LIGHTS_PROJECTION_H
#define PBRT_LIGHTS_PROJECTION_H

// lights/projection.h*
#include "pbrt.h"
#include "light.h"
#include "shape.h"
#include "mipmap.h"

// ProjectionLight Declarations
class ProjectionLight : public Light {
public:
    // ProjectionLight Public Methods
    ProjectionLight(const Transform &light2world, const Spectrum &intensity,
        const string &texname, float fov);
    ~ProjectionLight();
    Spectrum Sample_L(const Point &p, float pEpsilon, const LightSample &ls, float time,
        Vector *wi, float *pdf, VisibilityTester *vis) const;
    bool IsDeltaLight() const { return true; }
    Spectrum Projection(const Vector &w) const;
    Spectrum Power(const Scene *) const;
    Spectrum Sample_L(const Scene *scene, const LightSample &ls, float u1, float u2,
            float time, Ray *ray, Normal *Ns, float *pdf) const;
    float Pdf(const Point &, const Vector &) const;
private:
    // ProjectionLight Private Data
    MIPMap<RGBSpectrum> *projectionMap;
    Point lightPos;
    Spectrum Intensity;
    Transform lightProjection;
    float hither, yon;
    float screenX0, screenX1, screenY0, screenY1;
    float cosTotalWidth;
};


ProjectionLight *CreateProjectionLight(const Transform &light2world,
        const ParamSet &paramSet);

#endif // PBRT_LIGHTS_PROJECTION_H
