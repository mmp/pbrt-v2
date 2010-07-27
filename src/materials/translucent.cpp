
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


// materials/translucent.cpp*
#include "stdafx.h"
#include "materials/translucent.h"
#include "spectrum.h"
#include "reflection.h"
#include "paramset.h"
#include "texture.h"

// TranslucentMaterial Method Definitions
BSDF *TranslucentMaterial::GetBSDF(const DifferentialGeometry &dgGeom, const DifferentialGeometry &dgShading, MemoryArena &arena) const {
    float ior = 1.5f;
    DifferentialGeometry dgs;
    if (bumpMap)
        Bump(bumpMap, dgGeom, dgShading, &dgs);
    else
        dgs = dgShading;
    BSDF *bsdf = BSDF_ALLOC(arena, BSDF)(dgs, dgGeom.nn, ior);

    Spectrum r = reflect->Evaluate(dgs).Clamp();
    Spectrum t = transmit->Evaluate(dgs).Clamp();
    if (r.IsBlack() && t.IsBlack()) return bsdf;

    Spectrum kd = Kd->Evaluate(dgs).Clamp();
    if (!kd.IsBlack()) {
        if (!r.IsBlack()) bsdf->Add(BSDF_ALLOC(arena, Lambertian)(r * kd));
        if (!t.IsBlack()) bsdf->Add(BSDF_ALLOC(arena, BRDFToBTDF)(BSDF_ALLOC(arena, Lambertian)(t * kd)));
    }
    Spectrum ks = Ks->Evaluate(dgs).Clamp();
    if (!ks.IsBlack()) {
        float rough = roughness->Evaluate(dgs);
        if (!r.IsBlack()) {
            Fresnel *fresnel = BSDF_ALLOC(arena, FresnelDielectric)(ior, 1.f);
            bsdf->Add(BSDF_ALLOC(arena, Microfacet)(r * ks, fresnel,
                BSDF_ALLOC(arena, Blinn)(1.f / rough)));
        }
        if (!t.IsBlack()) {
            Fresnel *fresnel = BSDF_ALLOC(arena, FresnelDielectric)(ior, 1.f);
            bsdf->Add(BSDF_ALLOC(arena, BRDFToBTDF)(BSDF_ALLOC(arena, Microfacet)(t * ks, fresnel,
                BSDF_ALLOC(arena, Blinn)(1.f / rough))));
        }
    }
    return bsdf;
}


TranslucentMaterial *CreateTranslucentMaterial(const Transform &xform,
        const TextureParams &mp) {
    Reference<Texture<Spectrum> > Kd = mp.GetSpectrumTexture("Kd", Spectrum(0.25f));
    Reference<Texture<Spectrum> > Ks = mp.GetSpectrumTexture("Ks", Spectrum(0.25f));
    Reference<Texture<Spectrum> > reflect = mp.GetSpectrumTexture("reflect", Spectrum(0.5f));
    Reference<Texture<Spectrum> > transmit = mp.GetSpectrumTexture("transmit", Spectrum(0.5f));
    Reference<Texture<float> > roughness = mp.GetFloatTexture("roughness", .1f);
    Reference<Texture<float> > bumpMap = mp.GetFloatTexture("bumpmap", 0.f);
    return new TranslucentMaterial(Kd, Ks, roughness, reflect, transmit, bumpMap);
}


