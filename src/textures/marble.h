
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

#ifndef PBRT_TEXTURES_MARBLE_H
#define PBRT_TEXTURES_MARBLE_H

// textures/marble.h*
#include "pbrt.h"
#include "texture.h"
#include "paramset.h"

// MarbleTexture Declarations
class MarbleTexture : public Texture<Spectrum> {
public:
    // MarbleTexture Public Methods
    ~MarbleTexture() {
        delete mapping;
    }
    MarbleTexture(int oct, float roughness, float sc, float var,
                  TextureMapping3D *map)
        : octaves(oct), omega(roughness), scale(sc), variation(var),
          mapping(map) { }
    Spectrum Evaluate(const DifferentialGeometry &dg) const {
        Vector dpdx, dpdy;
        Point P = mapping->Map(dg, &dpdx, &dpdy);
        P *= scale;
        float marble = P.y + variation *
                       FBm(P, scale * dpdx, scale * dpdy, omega, octaves);
        float t = .5f + .5f * sinf(marble);
        // Evaluate marble spline at _t_
        static float c[][3] = { { .58f, .58f, .6f }, { .58f, .58f, .6f }, { .58f, .58f, .6f },
            { .5f, .5f, .5f }, { .6f, .59f, .58f }, { .58f, .58f, .6f },
            { .58f, .58f, .6f }, {.2f, .2f, .33f }, { .58f, .58f, .6f }, };
#define NC  sizeof(c) / sizeof(c[0])
#define NSEG (NC-3)
        int first = Floor2Int(t * NSEG);
        t = (t * NSEG - first);
        Spectrum c0 = Spectrum::FromRGB(c[first]);
        Spectrum c1 = Spectrum::FromRGB(c[first+1]);
        Spectrum c2 = Spectrum::FromRGB(c[first+2]);
        Spectrum c3 = Spectrum::FromRGB(c[first+3]);
        // Bezier spline evaluated with de Castilejau's algorithm
        Spectrum s0 = (1.f - t) * c0 + t * c1;
        Spectrum s1 = (1.f - t) * c1 + t * c2;
        Spectrum s2 = (1.f - t) * c2 + t * c3;
        s0 = (1.f - t) * s0 + t * s1;
        s1 = (1.f - t) * s1 + t * s2;
        // Extra scale of 1.5 to increase variation among colors
        return 1.5f * ((1.f - t) * s0 + t * s1);
    }
private:
    // MarbleTexture Private Data
    int octaves;
    float omega, scale, variation;
    TextureMapping3D *mapping;
};


Texture<float> *CreateMarbleFloatTexture(const Transform &tex2world,
        const TextureParams &tp);
MarbleTexture *CreateMarbleSpectrumTexture(const Transform &tex2world,
        const TextureParams &tp);

#endif // PBRT_TEXTURES_MARBLE_H
