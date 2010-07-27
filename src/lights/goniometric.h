
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

#ifndef PBRT_LIGHTS_GONIOMETRIC_H
#define PBRT_LIGHTS_GONIOMETRIC_H

// lights/goniometric.h*
#include "pbrt.h"
#include "light.h"
#include "shape.h"
#include "scene.h"
#include "mipmap.h"

// GonioPhotometricLight Declarations
class GonioPhotometricLight : public Light {
public:
    // GonioPhotometricLight Public Methods
    GonioPhotometricLight(const Transform &light2world, const Spectrum &, const
    string &texname);
    Spectrum Sample_L(const Point &p, float pEpsilon, const LightSample &ls,
        float time, Vector *wi, float *pdf, VisibilityTester *vis) const;
    ~GonioPhotometricLight() { delete mipmap; }
    bool IsDeltaLight() const { return true; }
    Spectrum Scale(const Vector &w) const {
        Vector wp = Normalize(WorldToLight(w));
        swap(wp.y, wp.z);
        float theta = SphericalTheta(wp);
        float phi   = SphericalPhi(wp);
        float s = phi * INV_TWOPI, t = theta * INV_PI;
        return (mipmap == NULL) ? 1.f :
               Spectrum(mipmap->Lookup(s, t, SPECTRUM_ILLUMINANT));
    }
    Spectrum Power(const Scene *) const;
    Spectrum Sample_L(const Scene *scene, const LightSample &ls, float u1, float u2,
        float time, Ray *ray, Normal *Ns, float *pdf) const;
    float Pdf(const Point &, const Vector &) const;
private:
    // GonioPhotometricLight Private Data
    Point lightPos;
    Spectrum Intensity;
    MIPMap<RGBSpectrum> *mipmap;
};


GonioPhotometricLight *CreateGoniometricLight(const Transform &light2world,
        const ParamSet &paramSet);

#endif // PBRT_LIGHTS_GONIOMETRIC_H
