
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


// lights/spot.cpp*
#include "stdafx.h"
#include "lights/spot.h"
#include "paramset.h"
#include "montecarlo.h"

// SpotLight Method Definitions
SpotLight::SpotLight(const Transform &light2world,
                     const Spectrum &intensity, float width, float fall)
    : Light(light2world) {
    lightPos = LightToWorld(Point(0,0,0));
    Intensity = intensity;
    cosTotalWidth = cosf(Radians(width));
    cosFalloffStart = cosf(Radians(fall));
}


Spectrum SpotLight::Sample_L(const Point &p, float pEpsilon,
        const LightSample &ls, float time, Vector *wi,
        float *pdf, VisibilityTester *visibility) const {
    *wi = Normalize(lightPos - p);
    *pdf = 1.f;
    visibility->SetSegment(p, pEpsilon, lightPos, 0., time);
    return Intensity * Falloff(-*wi) / DistanceSquared(lightPos, p);
}


float SpotLight::Falloff(const Vector &w) const {
    Vector wl = Normalize(WorldToLight(w));
    float costheta = wl.z;
    if (costheta < cosTotalWidth)     return 0.;
    if (costheta > cosFalloffStart)   return 1.;
    // Compute falloff inside spotlight cone
    float delta = (costheta - cosTotalWidth) /
                  (cosFalloffStart - cosTotalWidth);
    return delta*delta*delta*delta;
}


Spectrum SpotLight::Power(const Scene *) const {
    return Intensity * 2.f * M_PI *
           (1.f - .5f * (cosFalloffStart + cosTotalWidth));
}


SpotLight *CreateSpotLight(const Transform &l2w, const ParamSet &paramSet) {
    Spectrum I = paramSet.FindOneSpectrum("I", Spectrum(1.0));
    Spectrum sc = paramSet.FindOneSpectrum("scale", Spectrum(1.0));
    float coneangle = paramSet.FindOneFloat("coneangle", 30.);
    float conedelta = paramSet.FindOneFloat("conedeltaangle", 5.);
    // Compute spotlight world to light transformation
    Point from = paramSet.FindOnePoint("from", Point(0,0,0));
    Point to = paramSet.FindOnePoint("to", Point(0,0,1));
    Vector dir = Normalize(to - from);
    Vector du, dv;
    CoordinateSystem(dir, &du, &dv);
    Transform dirToZ =
        Transform(Matrix4x4( du.x,  du.y,  du.z, 0.,
                             dv.x,  dv.y,  dv.z, 0.,
                            dir.x, dir.y, dir.z, 0.,
                                0,     0,     0, 1.));
    Transform light2world = l2w *
        Translate(Vector(from.x, from.y, from.z)) * Inverse(dirToZ);
    return new SpotLight(light2world, I * sc, coneangle,
        coneangle-conedelta);
}


float SpotLight::Pdf(const Point &, const Vector &) const {
    return 0.;
}


Spectrum SpotLight::Sample_L(const Scene *scene, const LightSample &ls,
        float u1, float u2, float time, Ray *ray, Normal *Ns,
        float *pdf) const {
    Vector v = UniformSampleCone(ls.uPos[0], ls.uPos[1], cosTotalWidth);
    *ray = Ray(lightPos, LightToWorld(v), 0.f, INFINITY, time);
    *Ns = (Normal)ray->d;
    *pdf = UniformConePdf(cosTotalWidth);
    return Intensity * Falloff(ray->d);
}


