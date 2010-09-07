
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


// integrators/useprobes.cpp*
#include "stdafx.h"
#include "integrators/useprobes.h"
#include "integrators/photonmap.h"
#include "integrators/directlighting.h"
#include "integrators/emission.h"
#include "parallel.h"
#include "scene.h"
#include "progressreporter.h"
#include "camera.h"
#include "intersection.h"
#include "paramset.h"
#include "montecarlo.h"

// UseRadianceProbes Method Definitions
UseRadianceProbes *CreateRadianceProbesSurfaceIntegrator(const ParamSet &paramSet) {
    string filename = paramSet.FindOneFilename("filename", "probes.out");
    return new UseRadianceProbes(filename);
}


UseRadianceProbes::UseRadianceProbes(const string &filename) {
    lightSampleOffsets = NULL;
    bsdfSampleOffsets = NULL;
    // Read precomputed radiance probe values from file
    FILE *f = fopen(filename.c_str(), "r");
    if (f) {
        if (fscanf(f, "%d %d %d", &lmax, &includeDirectInProbes,
                   &includeIndirectInProbes) != 3 ||
            fscanf(f, "%d %d %d", &nProbes[0], &nProbes[1], &nProbes[2]) != 3 ||
            fscanf(f, "%f %f %f %f %f %f", &bbox.pMin.x, &bbox.pMin.y, &bbox.pMin.z,
                &bbox.pMax.x, &bbox.pMax.y, &bbox.pMax.z) != 6)
            Severe("Error reading data from radiance probe file \"%s\"", filename.c_str());
    
        c_in = new Spectrum[SHTerms(lmax) * nProbes[0] * nProbes[1] * nProbes[2]];
        int offset = 0;
        for (int i = 0; i < nProbes[0] * nProbes[1] * nProbes[2]; ++i) {
            for (int j = 0; j < SHTerms(lmax); ++j)
                if (!c_in[offset++].Read(f))
                    Severe("Error reading data from radiance probe file \"%s\"",
                        filename.c_str());
        }
        fclose(f);
    }
    else
        Error("Unable to read saved radiance volume values from file \"%s\"",
              filename.c_str());
}


UseRadianceProbes::~UseRadianceProbes() {
    delete[] lightSampleOffsets;
    delete[] bsdfSampleOffsets;
    delete[] c_in;
}


void UseRadianceProbes::RequestSamples(Sampler *sampler, Sample *sample, const Scene *scene) {
    // Allocate and request samples for sampling all lights
    uint32_t nLights = scene->lights.size();
    lightSampleOffsets = new LightSampleOffsets[nLights];
    bsdfSampleOffsets = new BSDFSampleOffsets[nLights];
    for (uint32_t i = 0; i < nLights; ++i) {
        const Light *light = scene->lights[i];
        int nSamples = light->nSamples;
        if (sampler) nSamples = sampler->RoundSize(nSamples);
        lightSampleOffsets[i] = LightSampleOffsets(nSamples, sample);
        bsdfSampleOffsets[i] = BSDFSampleOffsets(nSamples, sample);
    }
}


Spectrum UseRadianceProbes::Li(const Scene *scene, const Renderer *renderer,
            const RayDifferential &ray, const Intersection &isect,
            const Sample *sample, RNG &rng, MemoryArena &arena) const {
    Spectrum L(0.);
    Vector wo = -ray.d;
    // Compute emitted light if ray hit an area light source
    L += isect.Le(wo);

    // Evaluate BSDF at hit point
    BSDF *bsdf = isect.GetBSDF(ray, arena);
    const Point &p = bsdf->dgShading.p;
    const Normal &n = bsdf->dgShading.nn;
    // Compute reflection for radiance probes integrator
    if (!includeDirectInProbes)
        L += UniformSampleAllLights(scene, renderer, arena, p, n,
                wo, isect.rayEpsilon, ray.time, bsdf, sample, rng,
                lightSampleOffsets, bsdfSampleOffsets);

    // Compute reflected lighting using radiance probes

    // Compute probe coordinates and offsets for lookup point
    Vector offset = bbox.Offset(p);
    float voxx = (offset.x * nProbes[0]) - 0.5f;
    float voxy = (offset.y * nProbes[1]) - 0.5f;
    float voxz = (offset.z * nProbes[2]) - 0.5f;
    int vx = Floor2Int(voxx), vy = Floor2Int(voxy), vz = Floor2Int(voxz);
    float dx = voxx - vx, dy = voxy - vy, dz = voxz - vz;

    // Get radiance probe coefficients around lookup point
    const Spectrum *b000 = c_inXYZ(lmax, vx,   vy,   vz);
    const Spectrum *b100 = c_inXYZ(lmax, vx+1, vy,   vz);
    const Spectrum *b010 = c_inXYZ(lmax, vx,   vy+1, vz);
    const Spectrum *b110 = c_inXYZ(lmax, vx+1, vy+1, vz);
    const Spectrum *b001 = c_inXYZ(lmax, vx,   vy,   vz+1);
    const Spectrum *b101 = c_inXYZ(lmax, vx+1, vy,   vz+1);
    const Spectrum *b011 = c_inXYZ(lmax, vx,   vy+1, vz+1);
    const Spectrum *b111 = c_inXYZ(lmax, vx+1, vy+1, vz+1);

    // Compute incident radiance from radiance probe coefficients
    Spectrum *c_inp = arena.Alloc<Spectrum>(SHTerms(lmax));
    for (int i = 0; i < SHTerms(lmax); ++i) {
        // Do trilinear interpolation to compute SH coefficients at point
        Spectrum c00 = Lerp(dx, b000[i], b100[i]);
        Spectrum c10 = Lerp(dx, b010[i], b110[i]);
        Spectrum c01 = Lerp(dx, b001[i], b101[i]);
        Spectrum c11 = Lerp(dx, b011[i], b111[i]);
        Spectrum c0 = Lerp(dy, c00, c10);
        Spectrum c1 = Lerp(dy, c01, c11);
        c_inp[i] = Lerp(dz, c0, c1);
    }

    // Convolve incident radiance to compute irradiance function
    Spectrum *c_E = arena.Alloc<Spectrum>(SHTerms(lmax));
    SHConvolveCosTheta(lmax, c_inp, c_E);

    // Evaluate irradiance function and accumulate reflection
    Spectrum rho = bsdf->rho(wo, rng, BSDF_ALL_REFLECTION);
    float *Ylm = ALLOCA(float, SHTerms(lmax));
    SHEvaluate(Vector(Faceforward(n, wo)), lmax, Ylm);
    Spectrum E = 0.f;
    for (int i = 0; i < SHTerms(lmax); ++i)
        E += c_E[i] * Ylm[i];
    L += rho * INV_PI * E.Clamp();
    return L;
}


