
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

#ifndef PBRT_TEXTURES_MIX_H
#define PBRT_TEXTURES_MIX_H

// textures/mix.h*
#include "pbrt.h"
#include "texture.h"
#include "paramset.h"

// MixTexture Declarations
template <typename T> class MixTexture : public Texture<T> {
public:
    // MixTexture Public Methods
    MixTexture(Reference<Texture<T> > t1, Reference<Texture<T> > t2,
               Reference<Texture<float> > amt)
        : tex1(t1), tex2(t2), amount(amt) { }
    T Evaluate(const DifferentialGeometry &dg) const {
        T t1 = tex1->Evaluate(dg), t2 = tex2->Evaluate(dg);
        float amt = amount->Evaluate(dg);
        return (1.f - amt) * t1 + amt * t2;
    }
private:
    Reference<Texture<T> > tex1, tex2;
    Reference<Texture<float> > amount;
};


MixTexture<float> *CreateMixFloatTexture(const Transform &tex2world,
        const TextureParams &tp);
MixTexture<Spectrum> *CreateMixSpectrumTexture(const Transform &tex2world,
        const TextureParams &tp);

#endif // PBRT_TEXTURES_MIX_H
