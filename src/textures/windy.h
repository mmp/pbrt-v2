
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

#ifndef PBRT_TEXTURES_WINDY_H
#define PBRT_TEXTURES_WINDY_H

// textures/windy.h*
#include "pbrt.h"
#include "texture.h"
#include "paramset.h"

// WindyTexture Declarations
template <typename T> class WindyTexture : public Texture<T> {
public:
    // WindyTexture Public Methods
    ~WindyTexture() {
        delete mapping;
    }
    WindyTexture(TextureMapping3D *map) : mapping(map) { }
    T Evaluate(const DifferentialGeometry &dg) const {
        Vector dpdx, dpdy;
        Point P = mapping->Map(dg, &dpdx, &dpdy);
        float windStrength = FBm(.1f * P, .1f * dpdx, .1f * dpdy, .5f, 3);
        float waveHeight = FBm(P, dpdx, dpdy, .5f, 6);
        return fabsf(windStrength) * waveHeight;
    }
private:
    // WindyTexture Private Data
    TextureMapping3D *mapping;
};


WindyTexture<float> *CreateWindyFloatTexture(const Transform &tex2world,
        const TextureParams &tp);
WindyTexture<Spectrum> *CreateWindySpectrumTexture(const Transform &tex2world,
        const TextureParams &tp);

#endif // PBRT_TEXTURES_WINDY_H
