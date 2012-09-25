
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


// lights/projection.cpp*
#include "stdafx.h"
#include "lights/projection.h"
#include "montecarlo.h"
#include "paramset.h"
#include "imageio.h"

// ProjectionLight Method Definitions
ProjectionLight::ProjectionLight(const Transform &light2world,
        const Spectrum &intensity, const string &texname,
        float fov)
    : Light(light2world) {
    lightPos = LightToWorld(Point(0,0,0));
    Intensity = intensity;
    // Create _ProjectionLight_ MIP-map
    int width, height;
    RGBSpectrum *texels = ReadImage(texname, &width, &height);
    if (texels)
        projectionMap = new MIPMap<RGBSpectrum>(width, height, texels);
    else
        projectionMap = NULL;
    delete[] texels;

    // Initialize _ProjectionLight_ projection matrix
    float aspect = projectionMap ? float(width) / float(height) : 1.f;
    if (aspect > 1.f)  {
        screenX0 = -aspect; screenX1 = aspect;
        screenY0 = -1.f;    screenY1 = 1.f;
    }
    else {
        screenX0 = -1.f;            screenX1 = 1.f;
        screenY0 = -1.f / aspect;   screenY1 = 1.f / aspect;
    }
    hither = 1e-3f;
    yon = 1e30f;
    lightProjection = Perspective(fov, hither, yon);

    // Compute cosine of cone surrounding projection directions
    float opposite = tanf(Radians(fov) / 2.f);
    float tanDiag = opposite * sqrtf(1.f + 1.f/(aspect*aspect));
    cosTotalWidth = cosf(atanf(tanDiag));
}


ProjectionLight::~ProjectionLight() { delete projectionMap; }
Spectrum ProjectionLight::Sample_L(const Point &p, float pEpsilon,
         const LightSample &ls, float time, Vector *wi,
         float *pdf, VisibilityTester *visibility) const {
    *wi = Normalize(lightPos - p);
    *pdf = 1.f;
    visibility->SetSegment(p, pEpsilon, lightPos, 0., time);
    return Intensity * Projection(-*wi) /
        DistanceSquared(lightPos, p);
}


Spectrum ProjectionLight::Projection(const Vector &w) const {
    Vector wl = WorldToLight(w);
    // Discard directions behind projection light
    if (wl.z < hither) return 0.;

    // Project point onto projection plane and compute light
    Point Pl = lightProjection(Point(wl.x, wl.y, wl.z));
    if (Pl.x < screenX0 || Pl.x > screenX1 ||
        Pl.y < screenY0 || Pl.y > screenY1) return 0.;
    if (!projectionMap) return 1;
    float s = (Pl.x - screenX0) / (screenX1 - screenX0);
    float t = (Pl.y - screenY0) / (screenY1 - screenY0);
    return Spectrum(projectionMap->Lookup(s, t), SPECTRUM_ILLUMINANT);
}


Spectrum ProjectionLight::Power(const Scene *) const {
    return (projectionMap ? Spectrum(projectionMap->Lookup(.5f, .5f, .5f),
                                     SPECTRUM_ILLUMINANT) : Spectrum(1.f)) *
        Intensity * 2.f * M_PI * (1.f - cosTotalWidth);
}


ProjectionLight *CreateProjectionLight(const Transform &light2world,
        const ParamSet &paramSet) {
    Spectrum I = paramSet.FindOneSpectrum("I", Spectrum(1.0));
    Spectrum sc = paramSet.FindOneSpectrum("scale", Spectrum(1.0));
    float fov = paramSet.FindOneFloat("fov", 45.);
    string texname = paramSet.FindOneFilename("mapname", "");
    return new ProjectionLight(light2world, I * sc, texname, fov);
}


Spectrum ProjectionLight::Sample_L(const Scene *scene, const LightSample &ls,
        float u1, float u2, float time, Ray *ray, Normal *Ns, float *pdf) const {
    Vector v = UniformSampleCone(ls.uPos[0], ls.uPos[1], cosTotalWidth);
    *ray = Ray(lightPos, LightToWorld(v), 0.f, INFINITY, time);
    *Ns = (Normal)ray->d;
    *pdf = UniformConePdf(cosTotalWidth);
    return Intensity * Projection(ray->d);
}


float ProjectionLight::Pdf(const Point &, const Vector &) const {
    return 0.;
}


