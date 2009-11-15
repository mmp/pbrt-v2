
/*
    pbrt source code Copyright(c) 1998-2009 Matt Pharr and Greg Humphreys.

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


// materials/subsurface.cpp*
#include "materials/subsurface.h"
#include "textures/constant.h"
#include "volume.h"
#include "spectrum.h"
#include "reflection.h"
#include "texture.h"
#include "paramset.h"

// SubsurfaceMaterial Method Definitions
BSDF *SubsurfaceMaterial::GetBSDF(const DifferentialGeometry &dgGeom,
        const DifferentialGeometry &dgShading, MemoryArena &arena) const {
    // Allocate _BSDF_, possibly doing bump mapping with _bumpMap_
    DifferentialGeometry dgs;
    if (bumpMap)
        Bump(bumpMap, dgGeom, dgShading, &dgs);
    else
        dgs = dgShading;
    BSDF *bsdf = BSDF_ALLOC(arena, BSDF)(dgs, dgGeom.nn);
    Spectrum R = Kr->Evaluate(dgs).Clamp();
    float e = eta->Evaluate(dgs);
    if (!R.IsBlack())
        bsdf->Add(BSDF_ALLOC(arena, SpecularReflection)(R,
            BSDF_ALLOC(arena, FresnelDielectric)(1., e)));
    return bsdf;
}


BSSRDF *SubsurfaceMaterial::GetBSSRDF(const DifferentialGeometry &dgGeom,
        const DifferentialGeometry &dgShading, MemoryArena &arena) const {
    float e = eta->Evaluate(dgShading);
    return BSDF_ALLOC(arena, BSSRDF)(scale * sigma_a->Evaluate(dgShading),
        scale * sigma_prime_s->Evaluate(dgShading), e);
}


SubsurfaceMaterial *CreateSubsurfaceMaterial(const Transform &xform,
        const TextureParams &mp) {
    float sa_rgb[3] = { .0011f, .0024f, .014f }, sps_rgb[3] = { 2.55f, 3.21f, 3.77f };
    Spectrum sa = Spectrum::FromRGB(sa_rgb), sps = Spectrum::FromRGB(sps_rgb);
    string name = mp.FindString("name");
    bool found = GetVolumeScatteringProperties(name, &sa, &sps);
    if (name != "" && !found)
        Warning("Named material \"%s\" not found.  Using defaults.", name.c_str());
    float scale = mp.FindFloat("scale", 1.f);

    Reference<Texture<Spectrum> > sigma_a, sigma_prime_s;
    sigma_a = mp.GetSpectrumTexture("sigma_a", sa);
    sigma_prime_s = mp.GetSpectrumTexture("sigma_prime_s", sps);
    Reference<Texture<float> > ior = mp.GetFloatTexture("index", 1.3f);
    Reference<Texture<Spectrum> > Kr = mp.GetSpectrumTexture("Kr", Spectrum(1.f));
    Reference<Texture<float> > bumpMap = mp.GetFloatTexture("bumpmap", 0.f);
    return new SubsurfaceMaterial(scale, Kr, sigma_a, sigma_prime_s, ior, bumpMap);
}


