
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


