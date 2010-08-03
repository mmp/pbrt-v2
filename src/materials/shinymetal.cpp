
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

// shinymetal.cpp*
#include "stdafx.h"
#include "materials/shinymetal.h"
#include "spectrum.h"
#include "reflection.h"
#include "paramset.h"
#include "texture.h"

static inline Spectrum FresnelApproxEta(const Spectrum &Fr) {
    Spectrum reflectance = Fr.Clamp(0.f, .999f);
    return (Spectrum(1.) + Sqrt(reflectance)) /
        (Spectrum(1.) - Sqrt(reflectance));
}

static inline Spectrum FresnelApproxK(const Spectrum &Fr) {
    Spectrum reflectance = Fr.Clamp(0.f, .999f);
    return 2.f * Sqrt(reflectance / (Spectrum(1.) - reflectance));
}


// ShinyMetalMaterial Method Definitions
BSDF *ShinyMetalMaterial::GetBSDF(const DifferentialGeometry &dgGeom, 
                                  const DifferentialGeometry &dgShading,
                                  MemoryArena &arena) const {
    // Allocate _BSDF_, possibly doing bump-mapping with _bumpMap_
    DifferentialGeometry dgs;
    if (bumpMap)
        Bump(bumpMap, dgGeom, dgShading, &dgs);
    else
        dgs = dgShading;
    BSDF *bsdf = BSDF_ALLOC(arena, BSDF)(dgs, dgGeom.nn);
    Spectrum spec = Ks->Evaluate(dgs).Clamp();
    float rough = roughness->Evaluate(dgs);
    Spectrum R = Kr->Evaluate(dgs).Clamp();

    MicrofacetDistribution *md = BSDF_ALLOC(arena, Blinn)(1.f / rough);
    Spectrum k = 0.;
    Fresnel *frMf = BSDF_ALLOC(arena, FresnelConductor)(FresnelApproxEta(spec), k);
    Fresnel *frSr = BSDF_ALLOC(arena, FresnelConductor)(FresnelApproxEta(R), k);
    bsdf->Add(BSDF_ALLOC(arena, Microfacet)(1., frMf, md));
    bsdf->Add(BSDF_ALLOC(arena, SpecularReflection)(1., frSr));
    return bsdf;
}

ShinyMetalMaterial *CreateShinyMetalMaterial(const Transform &xform,
                                             const TextureParams &mp) {
    Reference<Texture<Spectrum> > Kr = mp.GetSpectrumTexture("Kr", Spectrum(1.f));
    Reference<Texture<Spectrum> > Ks = mp.GetSpectrumTexture("Ks", Spectrum(1.f));
    Reference<Texture<float> > roughness = mp.GetFloatTexture("roughness", .1f);
    Reference<Texture<float> > bumpMap = mp.GetFloatTexture("bumpmap", 0.f);
    return new ShinyMetalMaterial(Ks, roughness, Kr, bumpMap);
}
