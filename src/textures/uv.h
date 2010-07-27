
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

#ifndef PBRT_TEXTURES_UV_H
#define PBRT_TEXTURES_UV_H

// textures/uv.h*
#include "pbrt.h"
#include "texture.h"
#include "paramset.h"

// UVTexture Declarations
class UVTexture : public Texture<Spectrum> {
public:
    // UVTexture Public Methods
    UVTexture(TextureMapping2D *m) {
        mapping = m;
    }
    ~UVTexture() {
        delete mapping;
    }
    Spectrum Evaluate(const DifferentialGeometry &dg) const {
        float s, t, dsdx, dtdx, dsdy, dtdy;
        mapping->Map(dg, &s, &t, &dsdx, &dtdx, &dsdy, &dtdy);
        float rgb[3] = { s - Floor2Int(s), t - Floor2Int(t), 0.f };
        return Spectrum::FromRGB(rgb);
    }
private:
    TextureMapping2D *mapping;
};


Texture<float> *CreateUVFloatTexture(const Transform &tex2world,
        const TextureParams &tp);
UVTexture *CreateUVSpectrumTexture(const Transform &tex2world,
        const TextureParams &tp);

#endif // PBRT_TEXTURES_UV_H
