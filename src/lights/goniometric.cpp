
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


// lights/goniometric.cpp*
#include "stdafx.h"
#include "lights/goniometric.h"
#include "paramset.h"
#include "montecarlo.h"
#include "imageio.h"

// GonioPhotometricLight Method Definitions
Spectrum GonioPhotometricLight::Sample_L(const Point &p, float pEpsilon,
        const LightSample &ls, float time, Vector *wi, float *pdf, VisibilityTester *visibility) const {
    *wi = Normalize(lightPos - p);
    *pdf = 1.f;
    visibility->SetSegment(p, pEpsilon, lightPos, 0., time);
    return Intensity * Scale(-*wi) / DistanceSquared(lightPos, p);
}


GonioPhotometricLight::GonioPhotometricLight(const Transform &light2world,
        const Spectrum &intensity, const string &texname)
    : Light(light2world) {
    lightPos = LightToWorld(Point(0,0,0));
    Intensity = intensity;
    // Create _mipmap_ for _GonioPhotometricLight_
    int width, height;
    RGBSpectrum *texels = ReadImage(texname, &width, &height);
    if (texels) {
        mipmap = new MIPMap<RGBSpectrum>(width, height, texels);
        delete[] texels;
    }
    else mipmap = NULL;
}


Spectrum GonioPhotometricLight::Power(const Scene *) const {
    return 4.f * M_PI * Intensity *
        Spectrum(mipmap ? mipmap->Lookup(.5f, .5f, .5f) : 1.f, SPECTRUM_ILLUMINANT);
}


GonioPhotometricLight *CreateGoniometricLight(const Transform &light2world,
        const ParamSet &paramSet) {
    Spectrum I = paramSet.FindOneSpectrum("I", Spectrum(1.0));
    Spectrum sc = paramSet.FindOneSpectrum("scale", Spectrum(1.0));
    string texname = paramSet.FindOneFilename("mapname", "");
    return new GonioPhotometricLight(light2world, I * sc, texname);
}


Spectrum GonioPhotometricLight::Sample_L(const Scene *scene, const LightSample &ls,
        float u1, float u2, float time, Ray *ray, Normal *Ns, float *pdf) const {
    *ray = Ray(lightPos, UniformSampleSphere(ls.uPos[0], ls.uPos[1]), 0.f, INFINITY, time);
    *Ns = (Normal)ray->d;
    *pdf = UniformSpherePdf();
    return Intensity * Scale(ray->d);
}


float GonioPhotometricLight::Pdf(const Point &, const Vector &) const {
    return 0.;
}


