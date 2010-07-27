
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

#ifndef PBRT_TEXTURES_FBM_H
#define PBRT_TEXTURES_FBM_H

// textures/fbm.h*
#include "pbrt.h"
#include "texture.h"
#include "paramset.h"

// FBmTexture Declarations
template <typename T> class FBmTexture : public Texture<T> {
public:
    // FBmTexture Public Methods
    ~FBmTexture() {
        delete mapping;
    }
    FBmTexture(int oct, float roughness, TextureMapping3D *map)
        : omega(roughness), octaves(oct), mapping(map) { }
    T Evaluate(const DifferentialGeometry &dg) const {
        Vector dpdx, dpdy;
        Point P = mapping->Map(dg, &dpdx, &dpdy);
        return FBm(P, dpdx, dpdy, omega, octaves);
    }
private:
    // FBmTexture Private Data
    float omega;
    int octaves;
    TextureMapping3D *mapping;
};


FBmTexture<float> *CreateFBmFloatTexture(const Transform &tex2world,
        const TextureParams &tp);
FBmTexture<Spectrum> *CreateFBmSpectrumTexture(const Transform &tex2world,
        const TextureParams &tp);

#endif // PBRT_TEXTURES_FBM_H
