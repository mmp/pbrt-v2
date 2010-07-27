
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

#ifndef PBRT_TEXTURES_BILERP_H
#define PBRT_TEXTURES_BILERP_H

// textures/bilerp.h*
#include "pbrt.h"
#include "texture.h"
#include "paramset.h"

// BilerpTexture Declarations
template <typename T> class BilerpTexture : public Texture<T> {
public:
    // BilerpTexture Public Methods
    BilerpTexture(TextureMapping2D *m, const T &t00, const T &t01,
                  const T &t10, const T &t11)
        : mapping(m), v00(t00), v01(t01), v10(t10), v11(t11) {
    }
    ~BilerpTexture() {
        delete mapping;
    }
    T Evaluate(const DifferentialGeometry &dg) const {
        float s, t, dsdx, dtdx, dsdy, dtdy;
        mapping->Map(dg, &s, &t, &dsdx, &dtdx, &dsdy, &dtdy);
        return (1-s)*(1-t) * v00 + (1-s)*(  t) * v01 +
               (  s)*(1-t) * v10 + (  s)*(  t) * v11;
    }
private:
    // BilerpTexture Private Data
    TextureMapping2D *mapping;
    T v00, v01, v10, v11;
};


BilerpTexture<float> *CreateBilerpFloatTexture(const Transform &tex2world,
        const TextureParams &tp);
BilerpTexture<Spectrum> *CreateBilerpSpectrumTexture(const Transform &tex2world,
        const TextureParams &tp);

#endif // PBRT_TEXTURES_BILERP_H
