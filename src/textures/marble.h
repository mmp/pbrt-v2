
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
