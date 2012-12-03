
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

#ifndef PBRT_TEXTURES_DOTS_H
#define PBRT_TEXTURES_DOTS_H

// textures/dots.h*
#include "pbrt.h"
#include "texture.h"
#include "paramset.h"

// DotsTexture Declarations
template <typename T> class DotsTexture : public Texture<T> {
public:
    // DotsTexture Public Methods
    ~DotsTexture() {
        delete mapping;
    }
    DotsTexture(TextureMapping2D *m, Reference<Texture<T> > t1,
                Reference<Texture<T> > t2)
        : mapping(m), outsideDot(t1), insideDot(t2) {
    }
    T Evaluate(const DifferentialGeometry &dg) const {
        // Compute cell indices for dots
        float s, t, dsdx, dtdx, dsdy, dtdy;
        mapping->Map(dg, &s, &t, &dsdx, &dtdx, &dsdy, &dtdy);
        int sCell = Floor2Int(s + .5f), tCell = Floor2Int(t + .5f);

        // Return _insideDot_ result if point is inside dot
        if (Noise(sCell+.5f, tCell+.5f) > 0) {
            float radius = .35f;
            float maxShift = 0.5f - radius;
            float sCenter = sCell + maxShift *
                Noise(sCell + 1.5f, tCell + 2.8f);
            float tCenter = tCell + maxShift *
                Noise(sCell + 4.5f, tCell + 9.8f);
            float ds = s - sCenter, dt = t - tCenter;
            if (ds*ds + dt*dt < radius*radius)
                return insideDot->Evaluate(dg);
        }
        return outsideDot->Evaluate(dg);
    }
private:
    // DotsTexture Private Data
    TextureMapping2D *mapping;
    Reference<Texture<T> > outsideDot, insideDot;
};


DotsTexture<float> *CreateDotsFloatTexture(const Transform &tex2world,
        const TextureParams &tp);
DotsTexture<Spectrum> *CreateDotsSpectrumTexture(const Transform &tex2world,
        const TextureParams &tp);

#endif // PBRT_TEXTURES_DOTS_H
