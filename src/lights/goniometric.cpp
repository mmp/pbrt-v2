
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


