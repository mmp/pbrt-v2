
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
