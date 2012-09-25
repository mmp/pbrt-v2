
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


// lights/point.cpp*
#include "stdafx.h"
#include "lights/point.h"
#include "sh.h"
#include "scene.h"
#include "paramset.h"
#include "montecarlo.h"

// PointLight Method Definitions
PointLight::PointLight(const Transform &light2world,
                       const Spectrum &intensity)
    : Light(light2world) {
    lightPos = LightToWorld(Point(0,0,0));
    Intensity = intensity;
}


Spectrum PointLight::Sample_L(const Point &p, float pEpsilon,
         const LightSample &ls, float time, Vector *wi, float *pdf,
         VisibilityTester *visibility) const {
    *wi = Normalize(lightPos - p);
    *pdf = 1.f;
    visibility->SetSegment(p, pEpsilon, lightPos, 0., time);
    return Intensity / DistanceSquared(lightPos, p);
}


Spectrum PointLight::Power(const Scene *) const {
    return 4.f * M_PI * Intensity;
}


PointLight *CreatePointLight(const Transform &light2world,
        const ParamSet &paramSet) {
    Spectrum I = paramSet.FindOneSpectrum("I", Spectrum(1.0));
    Spectrum sc = paramSet.FindOneSpectrum("scale", Spectrum(1.0));
    Point P = paramSet.FindOnePoint("from", Point(0,0,0));
    Transform l2w = Translate(Vector(P.x, P.y, P.z)) * light2world;
    return new PointLight(l2w, I * sc);
}


float PointLight::Pdf(const Point &, const Vector &) const {
    return 0.;
}


Spectrum PointLight::Sample_L(const Scene *scene, const LightSample &ls,
        float u1, float u2, float time, Ray *ray, Normal *Ns,
        float *pdf) const {
    *ray = Ray(lightPos, UniformSampleSphere(ls.uPos[0], ls.uPos[1]),
               0.f, INFINITY, time);
    *Ns = (Normal)ray->d;
    *pdf = UniformSpherePdf();
    return Intensity;
}


void PointLight::SHProject(const Point &p, float pEpsilon, int lmax,
        const Scene *scene, bool computeLightVisibility, float time,
        RNG &rng, Spectrum *coeffs) const {
    for (int i = 0; i < SHTerms(lmax); ++i)
        coeffs[i] = 0.f;
    if (computeLightVisibility &&
        scene->IntersectP(Ray(p, Normalize(lightPos - p), pEpsilon,
                              Distance(lightPos, p), time)))
        return;
    // Project point light source to SH
    float *Ylm = ALLOCA(float, SHTerms(lmax));
    Vector wi = Normalize(lightPos - p);
    SHEvaluate(wi, lmax, Ylm);
    Spectrum Li = Intensity / DistanceSquared(lightPos, p);
    for (int i = 0; i < SHTerms(lmax); ++i)
        coeffs[i] = Li * Ylm[i];
}


