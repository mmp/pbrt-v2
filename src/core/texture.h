
/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CORE_TEXTURE_H
#define PBRT_CORE_TEXTURE_H

// core/texture.h*
#include "pbrt.h"
#include "spectrum.h"
#include "geometry.h"
#include "transform.h"
#include "memory.h"

// Texture Declarations
class TextureMapping2D {
public:
    // TextureMapping2D Interface
    virtual ~TextureMapping2D() { }
    virtual void Map(const DifferentialGeometry &dg,
                     float *s, float *t, float *dsdx, float *dtdx,
                     float *dsdy, float *dtdy) const = 0;
};


class UVMapping2D : public TextureMapping2D {
public:
    // UVMapping2D Public Methods
    UVMapping2D(float su = 1, float sv = 1, float du = 0, float dv = 0);
    void Map(const DifferentialGeometry &dg, float *s, float *t,
        float *dsdx, float *dtdx, float *dsdy, float *dtdy) const;
private:
    float su, sv, du, dv;
};


class SphericalMapping2D : public TextureMapping2D {
public:
    // SphericalMapping2D Public Methods
    SphericalMapping2D(const Transform &toSph)
        : WorldToTexture(toSph) {
    }
    void Map(const DifferentialGeometry &dg, float *s, float *t,
        float *dsdx, float *dtdx,
        float *dsdy, float *dtdy) const;
private:
    void sphere(const Point &P, float *s, float *t) const;
    Transform WorldToTexture;
};


class CylindricalMapping2D : public TextureMapping2D {
public:
    // CylindricalMapping2D Public Methods
    CylindricalMapping2D(const Transform &toCyl)
        : WorldToTexture(toCyl) {
    }
    void Map(const DifferentialGeometry &dg, float *s, float *t,
        float *dsdx, float *dtdx, float *dsdy, float *dtdy) const;
private:
    // CylindricalMapping2D Private Methods
    void cylinder(const Point &p, float *s, float *t) const {
        Vector vec = Normalize(WorldToTexture(p) - Point(0,0,0));
        *s = (M_PI + atan2f(vec.y, vec.x)) / (2.f * M_PI);
        *t = vec.z;
    }
    Transform WorldToTexture;
};


class PlanarMapping2D : public TextureMapping2D {
public:
    // PlanarMapping2D Public Methods
    void Map(const DifferentialGeometry &dg, float *s, float *t,
             float *dsdx, float *dtdx, float *dsdy, float *dtdy) const;
    PlanarMapping2D(const Vector &vv1, const Vector &vv2,
                    float dds = 0, float ddt = 0)
        : vs(vv1), vt(vv2), ds(dds), dt(ddt) { }
private:
    const Vector vs, vt;
    const float ds, dt;
};


class TextureMapping3D {
public:
    // TextureMapping3D Interface
    virtual ~TextureMapping3D() { }
    virtual Point Map(const DifferentialGeometry &dg,
                      Vector *dpdx, Vector *dpdy) const = 0;
};


class IdentityMapping3D : public TextureMapping3D {
public:
    IdentityMapping3D(const Transform &x)
        : WorldToTexture(x) { }
    Point Map(const DifferentialGeometry &dg, Vector *dpdx,
              Vector *dpdy) const;
private:
    Transform WorldToTexture;
};


template <typename T> class Texture : public ReferenceCounted {
public:
    // Texture Interface
    virtual T Evaluate(const DifferentialGeometry &) const = 0;
    virtual ~Texture() { }
};


float Lanczos(float, float tau=2);
float Noise(float x, float y = .5f, float z = .5f);
float Noise(const Point &P);
float FBm(const Point &P, const Vector &dpdx, const Vector &dpdy,
    float omega, int octaves);
float Turbulence(const Point &P, const Vector &dpdx, const Vector &dpdy,
    float omega, int octaves);

#endif // PBRT_CORE_TEXTURE_H
