
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

#ifndef PBRT_TEXTURES_WRINKLED_H
#define PBRT_TEXTURES_WRINKLED_H

// textures/wrinkled.h*
#include "pbrt.h"
#include "texture.h"
#include "paramset.h"

// WrinkledTexture Declarations
template <typename T> class WrinkledTexture : public Texture<T> {
public:
    // WrinkledTexture Public Methods
    ~WrinkledTexture() {
        delete mapping;
    }
    WrinkledTexture(int oct, float roughness, TextureMapping3D *map) {
        omega = roughness;
        octaves = oct;
        mapping = map;
    }
    T Evaluate(const DifferentialGeometry &dg) const {
        Vector dpdx, dpdy;
        Point P = mapping->Map(dg, &dpdx, &dpdy);
        return Turbulence(P, dpdx, dpdy, omega, octaves);
    }
private:
    // WrinkledTexture Private Data
    int octaves;
    float omega;
    TextureMapping3D *mapping;
};


WrinkledTexture<float> *CreateWrinkledFloatTexture(const Transform &tex2world,
        const TextureParams &tp);
WrinkledTexture<Spectrum> *CreateWrinkledSpectrumTexture(const Transform &tex2world,
        const TextureParams &tp);

#endif // PBRT_TEXTURES_WRINKLED_H
