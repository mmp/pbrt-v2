
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


