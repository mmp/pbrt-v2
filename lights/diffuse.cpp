
/*
    pbrt source code Copyright(c) 1998-2009 Matt Pharr and Greg Humphreys.

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


// lights/diffuse.cpp*
#include "lights/diffuse.h"
#include "paramset.h"
#include "montecarlo.h"

// DiffuseAreaLight Method Definitions
DiffuseAreaLight::~DiffuseAreaLight() {
    delete shapeSet;
}


DiffuseAreaLight::DiffuseAreaLight(const Transform &light2world,
                                   const Spectrum &le, int ns,
                                   const Reference<Shape> &s)
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
    if (getenv("PBRT_QUICK_RENDER")) nSamples = max(1, nSamples / 4);
    return new DiffuseAreaLight(light2world, L * sc, nSamples, shape);
}


Spectrum DiffuseAreaLight::Sample_L(const Point &p,
        float pEpsilon, const LightSample &ls, Vector *wi, float *pdf,
        VisibilityTester *visibility) const {
    Normal ns;
    Point ps = shapeSet->Sample(p, ls, &ns);
    *wi = Normalize(ps - p);
    *pdf = shapeSet->Pdf(p, *wi);
    visibility->SetSegment(p, pEpsilon, ps, 1e-3f);
    return L(ps, ns, -*wi);
}


float DiffuseAreaLight::Pdf(const Point &p,
                     const Vector &wi) const {
    return shapeSet->Pdf(p, wi);
}


Spectrum DiffuseAreaLight::Sample_L(const Scene *scene, const LightSample &ls,
                             float u1, float u2,
                             Ray *ray, Normal *Ns, float *pdf) const {
    ray->o = shapeSet->Sample(ls, Ns);
    ray->d = UniformSampleSphere(u1, u2);
    ray->mint = 1e-3f;
    if (Dot(ray->d, *Ns) < 0.) ray->d *= -1;
    *pdf = shapeSet->Pdf(ray->o) * INV_TWOPI;
    return L(ray->o, *Ns, ray->d);
}


