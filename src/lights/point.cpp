
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


