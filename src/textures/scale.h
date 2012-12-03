
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

#ifndef PBRT_TEXTURES_SCALE_H
#define PBRT_TEXTURES_SCALE_H

// textures/scale.h*
#include "pbrt.h"
#include "texture.h"
#include "paramset.h"

// ScaleTexture Declarations
template <typename T1, typename T2>
class ScaleTexture : public Texture<T2> {
public:
    // ScaleTexture Public Methods
    ScaleTexture(Reference<Texture<T1> > t1, Reference<Texture<T2> > t2)
        : tex1(t1), tex2(t2) { }
    T2 Evaluate(const DifferentialGeometry &dg) const {
        return tex1->Evaluate(dg) * tex2->Evaluate(dg);
    }
private:
    Reference<Texture<T1> > tex1;
    Reference<Texture<T2> > tex2;
};


ScaleTexture<float, float> *CreateScaleFloatTexture(const Transform &tex2world,
        const TextureParams &tp);
ScaleTexture<Spectrum, Spectrum> *CreateScaleSpectrumTexture(const Transform &tex2world,
        const TextureParams &tp);

#endif // PBRT_TEXTURES_SCALE_H
