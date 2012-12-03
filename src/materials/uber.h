
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

#ifndef PBRT_MATERIALS_UBER_H
#define PBRT_MATERIALS_UBER_H

// materials/uber.h*
#include "pbrt.h"
#include "material.h"

// UberMaterial Declarations
class UberMaterial : public Material {
public:
    UberMaterial(Reference<Texture<Spectrum> > kd,
        Reference<Texture<Spectrum> > ks,
        Reference<Texture<Spectrum> > kr,
        Reference<Texture<Spectrum> > kt,
        Reference<Texture<float> > rough,
        Reference<Texture<Spectrum> > op,
        Reference<Texture<float> > e,
        Reference<Texture<float> > bump) {
        Kd = kd;
        Ks = ks;
        Kr = kr;
        Kt = kt;
        roughness = rough;
        opacity = op;
        eta = e;
        bumpMap = bump;
    }
    BSDF *GetBSDF(const DifferentialGeometry &dgGeom, const DifferentialGeometry &dgShading, MemoryArena &arena) const;
private:
    // UberMaterial Private Data
    Reference<Texture<Spectrum> > Kd, Ks, Kr, Kt, opacity;
    Reference<Texture<float> > roughness, eta, bumpMap;
};


UberMaterial *CreateUberMaterial(const Transform &xform,
        const TextureParams &mp);

#endif // PBRT_MATERIALS_UBER_H
