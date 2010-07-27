
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
