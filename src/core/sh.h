
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

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CORE_SH_H
#define PBRT_CORE_SH_H

// core/sh.h*
#include "pbrt.h"
#include "geometry.h"
#include "spectrum.h"

// Spherical Harmonics Declarations
inline int SHTerms(int lmax) {
    return (lmax + 1) * (lmax + 1);
}


inline int SHIndex(int l, int m) {
    return l*l + l + m;
}


void SHEvaluate(const Vector &v, int lmax, float *out);
void SHWriteImage(const char *filename, const Spectrum *c, int lmax, int yres);
template <typename Func>
void SHProjectCube(Func func, const Point &p, int res, int lmax,
                   Spectrum *coeffs) {
    float *Ylm = ALLOCA(float, SHTerms(lmax));
    for (int u = 0; u < res; ++u) {
        float fu = -1.f + 2.f * (float(u) + 0.5f) / float(res);
        for (int v = 0; v < res; ++v) {
            float fv = -1.f + 2.f * (float(v) + 0.5f) / float(res);
            // Incorporate results from $+z$ face to coefficients
            Vector w(fu, fv, 1);
            SHEvaluate(Normalize(w), lmax, Ylm);
            Spectrum f = func(u, v, p, w);
            float dA = 1.f / powf(Dot(w, w), 3.f/2.f);
            for (int k = 0; k < SHTerms(lmax); ++k)
                coeffs[k] += f * Ylm[k] * dA * (4.f / (res * res));

            // Incorporate results from other faces to coefficients
            w = Vector(fu, fv, -1);
            SHEvaluate(Normalize(w), lmax, Ylm);
            f = func(u, v, p, w);
            for (int k = 0; k < SHTerms(lmax); ++k)
                coeffs[k] += f * Ylm[k] * dA * (4.f / (res * res));
            w = Vector(fu, 1, fv);
            SHEvaluate(Normalize(w), lmax, Ylm);
            f = func(u, v, p, w);
            for (int k = 0; k < SHTerms(lmax); ++k)
                coeffs[k] += f * Ylm[k] * dA * (4.f / (res * res));
            w = Vector(fu, -1, fv);
            SHEvaluate(Normalize(w), lmax, Ylm);
            f = func(u, v, p, w);
            for (int k = 0; k < SHTerms(lmax); ++k)
                coeffs[k] += f * Ylm[k] * dA * (4.f / (res * res));
            w = Vector(1, fu, fv);
            SHEvaluate(Normalize(w), lmax, Ylm);
            f = func(u, v, p, w);
            for (int k = 0; k < SHTerms(lmax); ++k)
                coeffs[k] += f * Ylm[k] * dA * (4.f / (res * res));
            w = Vector(-1, fu, fv);
            SHEvaluate(Normalize(w), lmax, Ylm);
            f = func(u, v, p, w);
            for (int k = 0; k < SHTerms(lmax); ++k)
                coeffs[k] += f * Ylm[k] * dA * (4.f / (res * res));
        }
    }
}


void SHProjectIncidentDirectRadiance(const Point &p, float pEpsilon, float time,
    MemoryArena &arena, const Scene *scene, bool computeLightVisibility,
    int lmax, RNG &rng, Spectrum *c_d);
void SHProjectIncidentIndirectRadiance(const Point &p, float pEpsilon,
    float time, const Renderer *renderer, Sample *origSample,
    const Scene *scene, int lmax, RNG &rng, int nSamples, Spectrum *c_i);
void SHReduceRinging(Spectrum *c, int lmax, float lambda = .005f);
void SHRotate(const Spectrum *c_in, Spectrum *c_out, const Matrix4x4 &m,
              int lmax, MemoryArena &arena);
void SHRotateZ(const Spectrum *c_in, Spectrum *c_out, float alpha, int lmax);
void SHRotateXMinus(const Spectrum *c_in, Spectrum *c_out, int lmax);
void SHRotateXPlus(const Spectrum *c_in, Spectrum *c_out, int lmax);
//void SHSwapYZ(const Spectrum *c_in, Spectrum *c_out, int lmax);
void SHConvolveCosTheta(int lmax, const Spectrum *c_in, Spectrum *c_out);
void SHConvolvePhong(int lmax, float n, const Spectrum *c_in, Spectrum *c_out);
void SHComputeDiffuseTransfer(const Point &p, const Normal &n, float rayEpsilon,
    const Scene *scene, RNG &rng, int nSamples, int lmax, Spectrum *c_transfer);
void SHComputeTransferMatrix(const Point &p, float rayEpsilon,
    const Scene *scene, RNG &rng, int nSamples, int lmax, Spectrum *T);
void SHComputeBSDFMatrix(const Spectrum &Kd, const Spectrum &Ks,
    float roughness, RNG &rng, int nSamples, int lmax, Spectrum *B);
void SHMatrixVectorMultiply(const Spectrum *M, const Spectrum *v,
                            Spectrum *vout, int lmax);

#endif // PBRT_CORE_SH_H
