
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
