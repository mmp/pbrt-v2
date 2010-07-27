
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


// core/sh.cpp*
#include "stdafx.h"
#include "sh.h"
#include "scene.h"
#include "integrator.h"
#include "intersection.h"
#include "montecarlo.h"
#include "imageio.h"
#include <float.h>

// Spherical Harmonics Local Definitions
static void legendrep(float x, int lmax, float *out) {
#define P(l,m) out[SHIndex(l,m)]
    // Compute $m=0$ Legendre values using recurrence
    P(0,0) = 1.f;
    P(1,0) = x;
    for (int l = 2; l <= lmax; ++l)
    {
        P(l, 0) = ((2*l-1)*x*P(l-1,0) - (l-1)*P(l-2,0)) / l;
        Assert(!isnan(P(l, 0)));
        Assert(!isinf(P(l, 0)));
    }

    // Compute $m=l$ edge using Legendre recurrence
    float neg = -1.f;
    float dfact = 1.f;
    float xroot = sqrtf(max(0.f, 1.f - x*x));
    float xpow = xroot;
    for (int l = 1; l <= lmax; ++l) {
        P(l, l) = neg * dfact * xpow;
        Assert(!isnan(P(l, l)));
        Assert(!isinf(P(l, l)));
        neg *= -1.f;      // neg = (-1)^l
        dfact *= 2*l + 1; // dfact = (2*l-1)!!
        xpow *= xroot;    // xpow = powf(1.f - x*x, float(l) * 0.5f);
    }

    // Compute $m=l-1$ edge using Legendre recurrence
    for (int l = 2; l <= lmax; ++l)
    {
        P(l, l-1) = x * (2*l-1) * P(l-1, l-1);
        Assert(!isnan(P(l, l-1)));
        Assert(!isinf(P(l, l-1)));
    }

    // Compute $m=1, \ldots, l-2$ values using Legendre recurrence
    for (int l = 3; l <= lmax; ++l)
        for (int m = 1; m <= l-2; ++m)
        {
            P(l, m) = ((2 * (l-1) + 1) * x * P(l-1,m) -
                       (l-1+m) * P(l-2,m)) / (l - m);
            Assert(!isnan(P(l, m)));
            Assert(!isinf(P(l, m)));
        }
    #if 0
        // wrap up with the negative m ones now
        // P(l,-m)(x) = -1^m (l-m)!/(l+m)! P(l,m)(x)
        for (int l = 1; l <= lmax; ++l) {
            float fa = 1.f, fb = fact(2*l);
            // fa = fact(l+m), fb = fact(l-m)
            for (int m = -l; m < 0; ++m) {
                float neg = ((-m) & 0x1) ? -1.f : 1.f;
                P(l,m) = neg * fa/fb * P(l,-m);
                fb /= l-m;
                fa *= (l+m+1) > 1 ? (l+m+1) : 1.;
            }
        }
    #endif
#undef P
}


static inline float fact(float v);
static inline float divfact(int a, int b);
static inline float K(int l, int m) {
    return sqrtf((2.f * l + 1.f) * INV_FOURPI * divfact(l, m));
}


static inline float divfact(int a, int b) {
    if (b == 0) return 1.f;
    float fa = a, fb = fabsf(b);
    float v = 1.f;
    for (float x = fa-fb+1.f; x <= fa+fb; x += 1.f)
        v *= x;
    return 1.f / v;
}


// n!! = 1 if n==0 or 1, otherwise n * (n-2)!!
static float dfact(float v) {
    if (v <= 1.f) return 1.f;
    return v * dfact(v - 2.f);
}


static inline float fact(float v) {
    if (v <= 1.f) return 1.f;
    return v * fact(v - 1.f);
}


static void sinCosIndexed(float s, float c, int n,
                          float *sout, float *cout) {
    float si = 0, ci = 1;
    for (int i = 0; i < n; ++i) {
        // Compute $\sin{}i\phi$ and $\cos{}i\phi$ using recurrence
        *sout++ = si;
        *cout++ = ci;
        float oldsi = si;
        si = si * c + ci * s;
        ci = ci * c - oldsi * s;
    }
}


static void toZYZ(const Matrix4x4 &m, float *alpha, float *beta, float *gamma) {
#define M(a, b) (m.m[a][b])

    float sy = sqrtf(M(2,1)*M(2,1) + M(2,0)*M(2,0));
    if (sy > 16*FLT_EPSILON) {
        *gamma = -atan2f(M(1,2), -M(0,2));
        *beta  = -atan2f(sy, M(2,2));
        *alpha = -atan2f(M(2,1), M(2,0));
    } else {
        *gamma =  0;
        *beta  = -atan2f(sy, M(2,2));
        *alpha = -atan2f(-M(1,0), M(1,1));
    }
#undef M
}


static inline float lambda(float l) {
    return sqrtf((4.f * M_PI) / (2.f * l + 1.));
}



// Spherical Harmonics Definitions
void SHEvaluate(const Vector &w, int lmax, float *out) {
    if (lmax > 28)
        Severe("SHEvaluate() runs out of numerical precision for lmax > 28. "
               "If you need more bands, try recompiling using doubles.");
    // Compute Legendre polynomial values for $\cos\theta$
    Assert(w.Length() > .995f && w.Length() < 1.005f);
    legendrep(w.z, lmax, out);

    // Compute $K_l^m$ coefficients
    float *Klm = ALLOCA(float, SHTerms(lmax));
    for (int l = 0; l <= lmax; ++l)
        for (int m = -l; m <= l; ++m)
            Klm[SHIndex(l, m)] = K(l, m);

    // Compute $\sin\phi$ and $\cos\phi$ values
    float *sins = ALLOCA(float, lmax+1), *coss = ALLOCA(float, lmax+1);
    float xyLen = sqrtf(max(0.f, 1.f - w.z*w.z));
    if (xyLen == 0.f) {
        for (int i = 0; i <= lmax; ++i) sins[i] = 0.f;
        for (int i = 0; i <= lmax; ++i) coss[i] = 1.f;
    }
    else
        sinCosIndexed(w.y / xyLen, w.x / xyLen, lmax+1, sins, coss);

    // Apply SH definitions to compute final $(l,m)$ values
    static const float sqrt2 = sqrtf(2.f);
    for (int l = 0; l <= lmax; ++l) {
        for (int m = -l; m < 0; ++m)
        {
            out[SHIndex(l, m)] = sqrt2 * Klm[SHIndex(l, m)] *
                out[SHIndex(l, -m)] * sins[-m];
            Assert(!isnan(out[SHIndex(l,m)]));
            Assert(!isinf(out[SHIndex(l,m)]));
        }
        out[SHIndex(l, 0)] *= Klm[SHIndex(l, 0)];
        for (int m = 1; m <= l; ++m)
        {
            out[SHIndex(l, m)] *= sqrt2 * Klm[SHIndex(l, m)] * coss[m];
            Assert(!isnan(out[SHIndex(l,m)]));
            Assert(!isinf(out[SHIndex(l,m)]));
        }
    }
}


#if 0
// Believe this is correct, but not well tested
void SHEvaluate(float costheta, float cosphi, float sinphi, int lmax, float *out) {
    legendrep(costheta, lmax, out);

    float *Klm = ALLOCA(float, SHTerms(lmax));
    klm(lmax, Klm);

    float sqrt2 = sqrtf(2.f);
    float sins[(lmax+1)], coss[(lmax+1)];
    sinCosIndexed(sinphi, cosphi, lmax+1, sins, coss);

    for (int l = 0; l <= lmax; ++l) {
        for (int m = -l; m < 0; ++m)
            // sin(-x) = -sin(x)
            out[SHIndex(l, m)] = sqrt2 * Klm[SHIndex(l, m)] *
                out[SHIndex(l, -m)] * sins[-m];

        out[SHIndex(l, 0)] *= Klm[SHIndex(l, 0)];

        for (int m = 1; m <= l; ++m)
            out[SHIndex(l, m)] *= sqrt2 * Klm[SHIndex(l, m)] * coss[m];
    }
}


#endif
void SHWriteImage(const char *filename, const Spectrum *c, int lmax, int yres) {
    int xres = 2 * yres;
    float *rgb = new float[xres * yres * 3];
    float *Ylm = ALLOCA(float, SHTerms(lmax));
    for (int y = 0; y < yres; ++y) {
        float theta = (float(y) + 0.5f) / float(yres) * M_PI;
        for (int x = 0; x < xres; ++x) {
            float phi = (float(x) + 0.5f) / float(xres) * 2.f * M_PI;
            // Compute RGB color for direction $(\theta,\phi)$ from SH coefficients
            Vector w = SphericalDirection(sinf(theta), cosf(theta), phi);
            SHEvaluate(w, lmax, Ylm);
            Spectrum val = 0.f;
            for (int i = 0; i < SHTerms(lmax); ++i)
                val += Ylm[i] * c[i];
            int offset = xres * y + x;
            val.ToRGB(&rgb[3*offset]);
        }
    }

    WriteImage(filename, rgb, NULL, xres, yres, xres, yres, 0, 0);
    delete[] rgb;
}


void SHProjectIncidentDirectRadiance(const Point &p, float pEpsilon,
        float time, MemoryArena &arena, const Scene *scene,
        bool computeLightVis, int lmax, RNG &rng, Spectrum *c_d) {
    // Loop over light sources and sum their SH coefficients
    Spectrum *c = arena.Alloc<Spectrum>(SHTerms(lmax));
    for (uint32_t i = 0; i < scene->lights.size(); ++i) {
        Light *light = scene->lights[i];
        light->SHProject(p, pEpsilon, lmax, scene, computeLightVis, time,
                         rng, c);
        for (int j = 0; j < SHTerms(lmax); ++j)
            c_d[j] += c[j];
    }
    SHReduceRinging(c_d, lmax);
}


void SHProjectIncidentIndirectRadiance(const Point &p, float pEpsilon,
        float time, const Renderer *renderer, Sample *origSample,
        const Scene *scene, int lmax, RNG &rng, int ns, Spectrum *c_i) {
    Sample *sample = origSample->Duplicate(1);
    MemoryArena arena;
    uint32_t scramble[2] = { rng.RandomUInt(), rng.RandomUInt() };
    int nSamples = RoundUpPow2(ns);
    float *Ylm = ALLOCA(float, SHTerms(lmax));
    for (int i = 0; i < nSamples; ++i) {
        // Sample incident direction for radiance probe
        float u[2];
        Sample02(i, scramble, u);
        Vector wi = UniformSampleSphere(u[0], u[1]);
        float pdf = UniformSpherePdf();

        // Compute incident radiance along direction for probe
        Spectrum Li = 0.f;
        RayDifferential ray(p, wi, pEpsilon, INFINITY, time);

        // Fill in values in _sample_ for radiance probe ray
        sample->time = time;
        for (uint32_t j = 0; j < sample->n1D.size(); ++j)
            for (uint32_t k = 0; k < sample->n1D[j]; ++k)
                sample->oneD[j][k] = rng.RandomFloat();
        for (uint32_t j = 0; j < sample->n2D.size(); ++j)
            for (uint32_t k = 0; k < 2 * sample->n2D[j]; ++k)
                sample->twoD[j][k] = rng.RandomFloat();
        Li = renderer->Li(scene, ray, sample, rng, arena);

        // Update SH coefficients for probe sample point
        SHEvaluate(wi, lmax, Ylm);
        for (int j = 0; j < SHTerms(lmax); ++j)
            c_i[j] += Ylm[j] * Li / (pdf * nSamples);
        arena.FreeAll();
    }
    delete[] sample;
}


void SHReduceRinging(Spectrum *c, int lmax, float lambda) {
    for (int l = 0; l <= lmax; ++l) {
        float scale = 1.f / (1.f + lambda * l * l * (l + 1) * (l + 1));
        for (int m = -l; m <= l; ++m)
            c[SHIndex(l, m)] *= scale;
    }
}


void SHRotate(const Spectrum *c_in, Spectrum *c_out, const Matrix4x4 &m,
              int lmax, MemoryArena &arena) {
    float alpha, beta, gamma;
    toZYZ(m, &alpha, &beta, &gamma);
    Spectrum *work = arena.Alloc<Spectrum>(SHTerms(lmax));
    SHRotateZ(c_in, c_out, gamma, lmax);
    SHRotateXPlus(c_out, work, lmax);
    SHRotateZ(work, c_out, beta, lmax);
    SHRotateXMinus(c_out, work, lmax);
    SHRotateZ(work, c_out, alpha, lmax);
}


void SHRotateZ(const Spectrum *c_in, Spectrum *c_out, float alpha,
               int lmax) {
    Assert(c_in != c_out);
    c_out[0] = c_in[0];
    if (lmax == 0) return;
    // Precompute sine and cosine terms for $z$-axis SH rotation
    float *ct = ALLOCA(float, lmax+1);
    float *st = ALLOCA(float, lmax+1);
    sinCosIndexed(sinf(alpha), cosf(alpha), lmax+1, st, ct);
    for (int l = 1; l <= lmax; ++l) {
        // Rotate coefficients for band _l_ about $z$
        for (int m = -l; m < 0; ++m)
            c_out[SHIndex(l, m)] =
                ( ct[-m] * c_in[SHIndex(l,  m)] +
                 -st[-m] * c_in[SHIndex(l, -m)]);
        c_out[SHIndex(l, 0)] = c_in[SHIndex(l, 0)];
        for (int m = 1; m <= l; ++m)
            c_out[SHIndex(l, m)] =
                (ct[m] * c_in[SHIndex(l,  m)] +
                 st[m] * c_in[SHIndex(l, -m)]);
    }
}


void SHConvolveCosTheta(int lmax, const Spectrum *c_in,
                        Spectrum *c_out) {
    static const float c_costheta[18] = { 0.8862268925, 1.0233267546,
        0.4954159260, 0.0000000000, -0.1107783690, 0.0000000000,
        0.0499271341, 0.0000000000, -0.0285469331, 0.0000000000,
        0.0185080823, 0.0000000000, -0.0129818395, 0.0000000000,
        0.0096125342, 0.0000000000, -0.0074057109, 0.0000000000 };
    for (int l = 0; l <= lmax; ++l)
        for (int m = -l; m <= l; ++m) {
            int o = SHIndex(l, m);
            if (l < 18) c_out[o] = lambda(l) * c_in[o] * c_costheta[l];
            else        c_out[o] = 0.f;
        }
}


void SHConvolvePhong(int lmax, float n, const Spectrum *c_in,
        Spectrum *c_out) {
    for (int l = 0; l <= lmax; ++l) {
        float c_phong = expf(-(l*l) / (2.f * n));
        for (int m = -l; m <= l; ++m) {
            int o = SHIndex(l, m);
            c_out[o] = lambda(l) * c_in[o] * c_phong;
        }
    }
}


void SHComputeDiffuseTransfer(const Point &p, const Normal &n,
        float rayEpsilon, const Scene *scene, RNG &rng, int nSamples,
        int lmax, Spectrum *c_transfer) {
    for (int i = 0; i < SHTerms(lmax); ++i)
        c_transfer[i] = 0.f;
    uint32_t scramble[2] = { rng.RandomUInt(), rng.RandomUInt() };
    float *Ylm = ALLOCA(float, SHTerms(lmax));
    for (int i = 0; i < nSamples; ++i) {
        // Sample _i_th direction and compute estimate for transfer coefficients
        float u[2];
        Sample02(i, scramble, u);
        Vector w = UniformSampleSphere(u[0], u[1]);
        float pdf = UniformSpherePdf();
        if (Dot(w, n) > 0.f && !scene->IntersectP(Ray(p, w, rayEpsilon))) {
            // Accumulate contribution of direction $\w{}$ to transfer coefficients
            SHEvaluate(w, lmax, Ylm);
            for (int j = 0; j < SHTerms(lmax); ++j)
                c_transfer[j] += (Ylm[j] * AbsDot(w, n)) / (pdf * nSamples);
        }
    }
}


void SHComputeTransferMatrix(const Point &p, float rayEpsilon,
        const Scene *scene, RNG &rng, int nSamples, int lmax,
        Spectrum *T) {
    for (int i = 0; i < SHTerms(lmax)*SHTerms(lmax); ++i)
        T[i] = 0.f;
    uint32_t scramble[2] = { rng.RandomUInt(), rng.RandomUInt() };
    float *Ylm = ALLOCA(float, SHTerms(lmax));
    for (int i = 0; i < nSamples; ++i) {
        // Compute Monte Carlo estimate of $i$th sample for transfer matrix
        float u[2];
        Sample02(i, scramble, u);
        Vector w = UniformSampleSphere(u[0], u[1]);
        float pdf = UniformSpherePdf();
        if (!scene->IntersectP(Ray(p, w, rayEpsilon))) {
            // Update transfer matrix for unoccluded direction
            SHEvaluate(w, lmax, Ylm);
            for (int j = 0; j < SHTerms(lmax); ++j)
                for (int k = 0; k < SHTerms(lmax); ++k)
                    T[j*SHTerms(lmax)+k] += (Ylm[j] * Ylm[k]) / (pdf * nSamples);
        }
    }
}


void SHComputeBSDFMatrix(const Spectrum &Kd, const Spectrum &Ks,
        float roughness, RNG &rng, int nSamples, int lmax, Spectrum *B) {
    for (int i = 0; i < SHTerms(lmax)*SHTerms(lmax); ++i)
        B[i] = 0.f;
    // Create _BSDF_ for computing BSDF transfer matrix
    MemoryArena arena;
    DifferentialGeometry dg(Point(0,0,0), Vector(1,0,0), Vector(0,1,0),
         Normal(0,0,0), Normal(0,0,0), 0, 0, NULL);
    BSDF *bsdf = BSDF_ALLOC(arena, BSDF)(dg, Normal(0,0,1));
    bsdf->Add(BSDF_ALLOC(arena, Lambertian)(Spectrum(Kd)));
    Fresnel *fresnel = BSDF_ALLOC(arena, FresnelDielectric)(1.5f, 1.f);
    bsdf->Add(BSDF_ALLOC(arena, Microfacet)(Ks, fresnel,
                                            BSDF_ALLOC(arena, Blinn)(1.f / roughness)));

    // Precompute directions $\w{}$ and SH values for directions
    float *Ylm = new float[SHTerms(lmax) * nSamples];
    Vector *w = new Vector[nSamples];
    uint32_t scramble[2] = { rng.RandomUInt(), rng.RandomUInt() };
    for (int i = 0; i < nSamples; ++i) {
        float u[2];
        Sample02(i, scramble, u);
        w[i] = UniformSampleSphere(u[0], u[1]);
        SHEvaluate(w[i], lmax, &Ylm[SHTerms(lmax)*i]);
    }

    // Compute double spherical integral for BSDF matrix
    for (int osamp = 0; osamp < nSamples; ++osamp) {
        const Vector &wo = w[osamp];
        for (int isamp = 0; isamp < nSamples; ++isamp) {
            const Vector &wi = w[isamp];
            // Update BSDF matrix elements for sampled directions
            Spectrum f = bsdf->f(wo, wi);
            if (!f.IsBlack()) {
                float pdf = UniformSpherePdf() * UniformSpherePdf();
                f *= fabsf(CosTheta(wi)) / (pdf * nSamples * nSamples);
                for (int i = 0; i < SHTerms(lmax); ++i)
                    for (int j = 0; j < SHTerms(lmax); ++j)
                        B[i*SHTerms(lmax)+j] += f * Ylm[isamp*SHTerms(lmax)+j] *
                                                Ylm[osamp*SHTerms(lmax)+i];
            }
        }
    }

    // Free memory allocated for SH matrix computation
    delete[] w;
    delete[] Ylm;
}


void SHMatrixVectorMultiply(const Spectrum *M, const Spectrum *v,
        Spectrum *vout, int lmax) {
    for (int i = 0; i < SHTerms(lmax); ++i) {
        vout[i] = 0.f;
        for (int j = 0; j < SHTerms(lmax); ++j)
            vout[i] += M[SHTerms(lmax) * i + j] * v[j];
    }
}


