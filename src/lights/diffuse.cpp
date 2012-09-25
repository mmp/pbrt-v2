
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


// lights/diffuse.cpp*
#include "stdafx.h"
#include "lights/diffuse.h"
#include "paramset.h"
#include "montecarlo.h"

// DiffuseAreaLight Method Definitions
DiffuseAreaLight::~DiffuseAreaLight() {
    delete shapeSet;
}


DiffuseAreaLight::DiffuseAreaLight(const Transform &light2world,
        const Spectrum &le, int ns, const Reference<Shape> &s)
    : AreaLight(light2world, ns) {
    Lemit = le;
    shapeSet = new ShapeSet(s);
    area = shapeSet->Area();
}


Spectrum DiffuseAreaLight::Power(const Scene *) const {
    return Lemit * area * M_PI;
}


AreaLight *CreateDiffuseAreaLight(const Transform &light2world, const ParamSet &paramSet,
        const Reference<Shape> &shape) {
    Spectrum L = paramSet.FindOneSpectrum("L", Spectrum(1.0));
    Spectrum sc = paramSet.FindOneSpectrum("scale", Spectrum(1.0));
    int nSamples = paramSet.FindOneInt("nsamples", 1);
    if (PbrtOptions.quickRender) nSamples = max(1, nSamples / 4);
    return new DiffuseAreaLight(light2world, L * sc, nSamples, shape);
}


Spectrum DiffuseAreaLight::Sample_L(const Point &p, float pEpsilon,
        const LightSample &ls, float time, Vector *wi, float *pdf,
        VisibilityTester *visibility) const {
    PBRT_AREA_LIGHT_STARTED_SAMPLE();
    Normal ns;
    Point ps = shapeSet->Sample(p, ls, &ns);
    *wi = Normalize(ps - p);
    *pdf = shapeSet->Pdf(p, *wi);
    visibility->SetSegment(p, pEpsilon, ps, 1e-3f, time);
    Spectrum Ls = L(ps, ns, -*wi);
    PBRT_AREA_LIGHT_FINISHED_SAMPLE();
    return Ls;
}


float DiffuseAreaLight::Pdf(const Point &p, const Vector &wi) const {
    return shapeSet->Pdf(p, wi);
}


Spectrum DiffuseAreaLight::Sample_L(const Scene *scene,
        const LightSample &ls, float u1, float u2, float time,
        Ray *ray, Normal *Ns, float *pdf) const {
    PBRT_AREA_LIGHT_STARTED_SAMPLE();
    Point org = shapeSet->Sample(ls, Ns);
    Vector dir = UniformSampleSphere(u1, u2);
    if (Dot(dir, *Ns) < 0.) dir *= -1.f;
    *ray = Ray(org, dir, 1e-3f, INFINITY, time);
    *pdf = shapeSet->Pdf(org) * INV_TWOPI;
    Spectrum Ls = L(org, *Ns, dir);
    PBRT_AREA_LIGHT_FINISHED_SAMPLE();
    return Ls;
}


