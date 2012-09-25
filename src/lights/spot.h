
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

#ifndef PBRT_LIGHTS_SPOT_H
#define PBRT_LIGHTS_SPOT_H

// lights/spot.h*
#include "pbrt.h"
#include "light.h"
#include "shape.h"

// SpotLight Declarations
class SpotLight : public Light {
public:
    // SpotLight Public Methods
    SpotLight(const Transform &light2world, const Spectrum &, float width, float fall);
    Spectrum Sample_L(const Point &p, float pEpsilon, const LightSample &ls, float time,
        Vector *wi, float *pdf, VisibilityTester *vis) const;
    bool IsDeltaLight() const { return true; }
    float Falloff(const Vector &w) const;
    Spectrum Power(const Scene *) const;
    Spectrum Sample_L(const Scene *scene, const LightSample &ls,
        float u1, float u2, float time, Ray *ray, Normal *Ns, float *pdf) const;
    float Pdf(const Point &, const Vector &) const;
private:
    // SpotLight Private Data
    Point lightPos;
    Spectrum Intensity;
    float cosTotalWidth, cosFalloffStart;
};


SpotLight *CreateSpotLight(const Transform &l2w, const ParamSet &paramSet);

#endif // PBRT_LIGHTS_SPOT_H
