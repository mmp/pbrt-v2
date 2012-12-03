
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
