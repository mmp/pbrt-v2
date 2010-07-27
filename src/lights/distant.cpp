
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


// lights/distant.cpp*
#include "stdafx.h"
#include "lights/distant.h"
#include "paramset.h"
#include "montecarlo.h"

// DistantLight Method Definitions
DistantLight::DistantLight(const Transform &light2world,
        const Spectrum &radiance, const Vector &dir)
    : Light(light2world) {
    lightDir = Normalize(LightToWorld(dir));
    L = radiance;
}


Spectrum DistantLight::Sample_L(const Point &p, float pEpsilon,
        const LightSample &ls, float time, Vector *wi, float *pdf,
        VisibilityTester *visibility) const {
    *wi = lightDir;
    *pdf = 1.f;
    visibility->SetRay(p, pEpsilon, *wi, time);
    return L;
}


Spectrum DistantLight::Power(const Scene *scene) const {
    Point worldCenter;
    float worldRadius;
    scene->WorldBound().BoundingSphere(&worldCenter, &worldRadius);
    return L * M_PI * worldRadius * worldRadius;
}


DistantLight *CreateDistantLight(const Transform &light2world,
        const ParamSet &paramSet) {
    Spectrum L = paramSet.FindOneSpectrum("L", Spectrum(1.0));
    Spectrum sc = paramSet.FindOneSpectrum("scale", Spectrum(1.0));
    Point from = paramSet.FindOnePoint("from", Point(0,0,0));
    Point to = paramSet.FindOnePoint("to", Point(0,0,1));
    Vector dir = from-to;
    return new DistantLight(light2world, L * sc, dir);
}


float DistantLight::Pdf(const Point &, const Vector &) const {
    return 0.;
}


Spectrum DistantLight::Sample_L(const Scene *scene,
        const LightSample &ls, float u1, float u2, float time,
        Ray *ray, Normal *Ns, float *pdf) const {
    // Choose point on disk oriented toward infinite light direction
    Point worldCenter;
    float worldRadius;
    scene->WorldBound().BoundingSphere(&worldCenter, &worldRadius);
    Vector v1, v2;
    CoordinateSystem(lightDir, &v1, &v2);
    float d1, d2;
    ConcentricSampleDisk(ls.uPos[0], ls.uPos[1], &d1, &d2);
    Point Pdisk = worldCenter + worldRadius * (d1 * v1 + d2 * v2);

    // Set ray origin and direction for infinite light ray
    *ray = Ray(Pdisk + worldRadius * lightDir, -lightDir, 0.f, INFINITY,
               time);
    *Ns = (Normal)ray->d;

    *pdf = 1.f / (M_PI * worldRadius * worldRadius);
    return L;
}


