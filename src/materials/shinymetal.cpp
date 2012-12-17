
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
    Reference<Texture<float> > bumpMap = mp.GetFloatTextureOrNull("bumpmap");
    return new ShinyMetalMaterial(Ks, roughness, Kr, bumpMap);
}
