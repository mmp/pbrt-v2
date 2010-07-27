
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


// materials/substrate.cpp*
#include "stdafx.h"
#include "materials/substrate.h"
#include "spectrum.h"
#include "reflection.h"
#include "paramset.h"
#include "texture.h"

// SubstrateMaterial Method Definitions
BSDF *SubstrateMaterial::GetBSDF(const DifferentialGeometry &dgGeom, const DifferentialGeometry &dgShading, MemoryArena &arena) const {
    // Allocate _BSDF_, possibly doing bump mapping with _bumpMap_
    DifferentialGeometry dgs;
    if (bumpMap)
        Bump(bumpMap, dgGeom, dgShading, &dgs);
    else
        dgs = dgShading;
    BSDF *bsdf = BSDF_ALLOC(arena, BSDF)(dgs, dgGeom.nn);
    Spectrum d = Kd->Evaluate(dgs).Clamp();
    Spectrum s = Ks->Evaluate(dgs).Clamp();
    float u = nu->Evaluate(dgs);
    float v = nv->Evaluate(dgs);

    bsdf->Add(BSDF_ALLOC(arena, FresnelBlend)(d, s, BSDF_ALLOC(arena, Anisotropic)(1.f/u, 1.f/v)));
    return bsdf;
}


SubstrateMaterial *CreateSubstrateMaterial(const Transform &xform,
        const TextureParams &mp) {
    Reference<Texture<Spectrum> > Kd = mp.GetSpectrumTexture("Kd", Spectrum(.5f));
    Reference<Texture<Spectrum> > Ks = mp.GetSpectrumTexture("Ks", Spectrum(.5f));
    Reference<Texture<float> > uroughness = mp.GetFloatTexture("uroughness", .1f);
    Reference<Texture<float> > vroughness = mp.GetFloatTexture("vroughness", .1f);
    Reference<Texture<float> > bumpMap = mp.GetFloatTexture("bumpmap", 0.f);
    return new SubstrateMaterial(Kd, Ks, uroughness, vroughness, bumpMap);
}


