
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


// lights/infinite.cpp*
#include "lights/infinite.h"
#include "sh.h"
#include "montecarlo.h"
#include "paramset.h"
#include "imageio.h"

// InfiniteAreaLight Utility Classes
struct InfiniteAreaCube {
    // InfiniteAreaCube Public Methods
    InfiniteAreaCube(const InfiniteAreaLight *l, const Scene *s,
                     float t, bool cv, float pe) {
        light = l; scene = s; time = t;
        computeVis = cv; pEpsilon = pe;
    }
    Spectrum operator()(int, int, const Point &p, const Vector &w) {
        Ray ray(p, w, pEpsilon, INFINITY, time);
        if (!computeVis || !scene->IntersectP(ray))
            return light->Le(RayDifferential(ray));
        return 0.f;
    }
    const InfiniteAreaLight *light;
    const Scene *scene;
    float time, pEpsilon;
    bool computeVis;
};



// InfiniteAreaLight Method Definitions
InfiniteAreaLight::~InfiniteAreaLight() {
    delete uDistrib;
    for (u_int i = 0; i < vDistribs.size(); ++i)
        delete vDistribs[i];
    delete radianceMap;
}


InfiniteAreaLight::InfiniteAreaLight(const Transform &light2world,
                                     const Spectrum &L, int ns,
                                     const string &texmap)
    : Light(light2world, ns) {
    radianceMap = NULL;
    int width = 0, height = 0;
    RGBSpectrum *texels = NULL;
    if (texmap != "") {
        texels = ReadImage(texmap, &width, &height);
        for (int i = 0; i < width * height; ++i)
            texels[i] *= L.ToRGBSpectrum();
    }
    if (!texels) {
        width = height = 1;
        texels = new RGBSpectrum[1];
        texels[0] = L.ToRGBSpectrum();
    }
    radianceMap = new MIPMap<RGBSpectrum>(width, height, texels);
    delete[] texels;
    // Initialize sampling PDFs for infinite area light

    // Compute scalar-valued image from environment map
    float filter = 1.f / max(width, height);
    int nu = width, nv = height;
    float *img = new float[width*height];
    for (int u = 0; u < nu; ++u) {
        float up = (float)u / (float)nu;
        for (int v = 0; v < nv; ++v) {
            float vp = (float)v / (float)nv;
            img[v+u*nv] = radianceMap->Lookup(up, vp, filter).y();
        }
    }
    vector<float> func(max(nu, nv), 0.f);
    vector<float> sinVals(nv, 0.f);
    for (int i = 0; i < nv; ++i)
        sinVals[i] = sin(M_PI * float(i+.5)/float(nv));
    vDistribs.reserve(nu);
    for (int u = 0; u < nu; ++u) {
        // Compute sampling distribution for column _u_
        for (int v = 0; v < nv; ++v)
            func[v] = img[u*nv+v] * sinVals[v];
        vDistribs.push_back(new Distribution1D(&func[0], nv));
    }

    // Compute sampling distribution for columns of image
    for (int u = 0; u < nu; ++u)
        func[u] = vDistribs[u]->funcInt;
    uDistrib = new Distribution1D(&func[0], nu);
    delete[] img;
}


Spectrum InfiniteAreaLight::Power(const Scene *scene) const {
    Point worldCenter;
    float worldRadius;
    scene->WorldBound().BoundingSphere(&worldCenter, &worldRadius);
    return Spectrum(radianceMap->Lookup(.5f, .5f, .5f)) *
           M_PI * worldRadius * worldRadius;
}


Spectrum InfiniteAreaLight::Le(const RayDifferential &r) const {
    Vector w = r.d;
    // Compute infinite light radiance for direction
    Vector wh = Normalize(WorldToLight(w));
    float s = SphericalPhi(wh) * INV_TWOPI;
    float t = SphericalTheta(wh) * INV_PI;
    Spectrum L = radianceMap->Lookup(s, t);
    return L;
}


void InfiniteAreaLight::SHProject(const Point &p, float pEpsilon,
        int lmax, const Scene *scene, bool computeLightVis,
        float time, RNG &rng, Spectrum *coeffs) const {
    if (computeLightVis) {
        Light::SHProject(p, pEpsilon, lmax, scene, computeLightVis,
                         time, rng, coeffs);
        return;
    }
    for (int i = 0; i < SHTerms(lmax); ++i)
        coeffs[i] = 0.f;
    int ntheta = radianceMap->Height(), nphi = radianceMap->Width();
    if (min(ntheta, nphi) > 50) {
        // Project _InfiniteAreaLight_ to SH from lat-long representation

        // Precompute $\theta$ and $\phi$ values for lat-long map projection
        float *buf = new float[2*ntheta + 2*nphi];
        float *bufp = buf;
        float *sintheta = bufp;  bufp += ntheta;
        float *costheta = bufp;  bufp += ntheta;
        float *sinphi = bufp;    bufp += nphi;
        float *cosphi = bufp;
        for (int theta = 0; theta < ntheta; ++theta) {
            sintheta[theta] = sinf((theta + .5f)/ntheta * M_PI);
            costheta[theta] = cosf((theta + .5f)/ntheta * M_PI);
        }
        for (int phi = 0; phi < nphi; ++phi) {
            sinphi[phi] = sinf((phi + .5f)/nphi * 2.f * M_PI);
            cosphi[phi] = cosf((phi + .5f)/nphi * 2.f * M_PI);
        }
        float *Ylm = ALLOCA(float, SHTerms(lmax));
        for (int theta = 0; theta < ntheta; ++theta) {
            for (int phi = 0; phi < nphi; ++phi) {
                // Add _InfiniteAreaLight_ texel's contribution to SH coefficients
                Vector w = Vector(sintheta[theta] * cosphi[phi],
                                  sintheta[theta] * sinphi[phi],
                                  costheta[theta]);
                w = Normalize(LightToWorld(w));
                if (!computeLightVis || !scene->IntersectP(Ray(p, w, pEpsilon))) {
                    Spectrum Le = radianceMap->Texel(0, phi, theta);
                    SHEvaluate(w, lmax, Ylm);
                    for (int i = 0; i < SHTerms(lmax); ++i)
                        coeffs[i] += Le * Ylm[i] * sintheta[theta] *
                            (M_PI / ntheta) * (2.f * M_PI / nphi);
                }
            }
        }

        // Free memory used for lat-long theta and phi values
        delete[] buf;
    }
    else {
        // Project _InfiniteAreaLight_ to SH from cube map sampling
        SHProjectCube(InfiniteAreaCube(this, scene, time, computeLightVis,
                                       pEpsilon),
                      p, 200, lmax, coeffs);
    }
}


InfiniteAreaLight *CreateInfiniteLight(const Transform &light2world,
        const ParamSet &paramSet) {
    Spectrum L = paramSet.FindOneSpectrum("L", Spectrum(1.0));
    Spectrum sc = paramSet.FindOneSpectrum("scale", Spectrum(1.0));
    string texmap = paramSet.FindOneString("mapname", "");
    int nSamples = paramSet.FindOneInt("nsamples", 1);
    if (getenv("PBRT_QUICK_RENDER")) nSamples = max(1, nSamples / 4);
    return new InfiniteAreaLight(light2world, L * sc, nSamples, texmap);
}


Spectrum InfiniteAreaLight::Sample_L(const Point &p,
        float pEpsilon, const LightSample &ls, Vector *wi, float *pdf,
        VisibilityTester *visibility) const {
    // Find $(u,v)$ sample coordinates in infinite light texture
    float pdfs[2];
    float fu = uDistrib->Sample(ls.uPos[0], &pdfs[0]);
    int u = Clamp(Float2Int(fu), 0, uDistrib->count-1);
    float fv = vDistribs[u]->Sample(ls.uPos[1], &pdfs[1]);

    // Convert infinite light sample point to direction
    float theta = fv * vDistribs[u]->invCount * M_PI;
    float phi = fu * uDistrib->invCount * 2.f * M_PI;
    float costheta = cos(theta), sintheta = sin(theta);
    float sinphi = sin(phi), cosphi = cos(phi);
    *wi = LightToWorld(Vector(sintheta * cosphi, sintheta * sinphi,
                              costheta));

    // Compute PDF for sampled infinite light direction
    *pdf = (pdfs[0] * pdfs[1]) / (2.f * M_PI * M_PI * sintheta);

    // Return radiance value for infinite light direction
    visibility->SetRay(p, pEpsilon, *wi);
    return radianceMap->Lookup(fu * uDistrib->invCount,
                               fv * vDistribs[u]->invCount);
}


float InfiniteAreaLight::Pdf(const Point &,
        const Vector &w) const {
    Vector wi = WorldToLight(w);
    float theta = SphericalTheta(wi), phi = SphericalPhi(wi);
    int u = Clamp(Float2Int(phi * INV_TWOPI * uDistrib->count),
                  0, uDistrib->count-1);
    int v = Clamp(Float2Int(theta * INV_PI * vDistribs[u]->count),
                  0, vDistribs[u]->count-1);
    return (uDistrib->func[u] * vDistribs[u]->func[v]) /
           (uDistrib->funcInt * vDistribs[u]->funcInt) *
           1.f / (2.f * M_PI * M_PI * sinf(theta));
}


Spectrum InfiniteAreaLight::Sample_L(const Scene *scene,
        const LightSample &ls, float u1, float u2,
        Ray *ray, Normal *Ns, float *pdf) const {
    // Compute direction for infinite light sample ray

    // Find $(u,v)$ sample coordinates in infinite light texture
    float pdfs[2];
    float fu = uDistrib->Sample(ls.uPos[0], &pdfs[0]);
    int u = Clamp(Float2Int(fu), 0, uDistrib->count-1);
    float fv = vDistribs[u]->Sample(ls.uPos[1], &pdfs[1]);
    float theta = fv * vDistribs[u]->invCount * M_PI;
    float phi = fu * uDistrib->invCount * 2.f * M_PI;
    float costheta = cos(theta), sintheta = sin(theta);
    float sinphi = sin(phi), cosphi = cos(phi);
    ray->d = -LightToWorld(Vector(sintheta * cosphi, sintheta * sinphi,
                                  costheta));
    *Ns = (Normal)ray->d;

    // Compute origin for infinite light sample ray
    Point worldCenter;
    float worldRadius;
    scene->WorldBound().BoundingSphere(&worldCenter, &worldRadius);
    Vector v1, v2;
    CoordinateSystem(-ray->d, &v1, &v2);
    float d1, d2;
    ConcentricSampleDisk(u1, u2, &d1, &d2);
    Point Pdisk = worldCenter + worldRadius * (d1 * v1 + d2 * v2);
    ray->o = Pdisk + worldRadius * -ray->d;

    // Compute _InfiniteAreaLight_ ray PDF
    float directionPdf = (pdfs[0] * pdfs[1]) / (2.f * M_PI * M_PI * sintheta);
    float areaPdf = 1.f / (M_PI * worldRadius * worldRadius);
    *pdf = directionPdf * areaPdf;

    // Return radiance value for sampled infinite light ray
    return radianceMap->Lookup(fu * uDistrib->invCount,
                               fv * vDistribs[u]->invCount);
}


