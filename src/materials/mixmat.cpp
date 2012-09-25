
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


// materials/mixmat.cpp*
#include "stdafx.h"
#include "materials/mixmat.h"
#include "materials/matte.h"
#include "spectrum.h"
#include "reflection.h"
#include "paramset.h"
#include "texture.h"

// MixMaterial Method Definitions
BSDF *MixMaterial::GetBSDF(const DifferentialGeometry &dgGeom,
                           const DifferentialGeometry &dgShading,
                           MemoryArena &arena) const {
    BSDF *b1 = m1->GetBSDF(dgGeom, dgShading, arena);
    BSDF *b2 = m2->GetBSDF(dgGeom, dgShading, arena);
    Spectrum s1 = scale->Evaluate(dgShading).Clamp();
    Spectrum s2 = (Spectrum(1.f) - s1).Clamp();
    int n1 = b1->NumComponents(), n2 = b2->NumComponents();
    for (int i = 0; i < n1; ++i)
        b1->bxdfs[i] = BSDF_ALLOC(arena, ScaledBxDF)(b1->bxdfs[i], s1);
    for (int i = 0; i < n2; ++i)
        b1->Add(BSDF_ALLOC(arena, ScaledBxDF)(b2->bxdfs[i], s2));
    return b1;
}


MixMaterial *CreateMixMaterial(const Transform &xform,
        const TextureParams &mp, const Reference<Material> &m1,
        const Reference<Material> &m2) {
    Reference<Texture<Spectrum> > scale = mp.GetSpectrumTexture("amount",
        Spectrum(0.5f));
    return new MixMaterial(m1, m2, scale);
}


