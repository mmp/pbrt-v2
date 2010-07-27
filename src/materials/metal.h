
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

#ifndef PBRT_MATERIALS_METAL_H
#define PBRT_MATERIALS_METAL_H

// materials/metal.h*
#include "pbrt.h"
#include "material.h"
#include "spectrum.h"

// MetalMaterial Declarations
class MetalMaterial : public Material {
public:
    // MetalMaterial Public Methods
    MetalMaterial(Reference<Texture<Spectrum> > eta,
        Reference<Texture<Spectrum> > k, Reference<Texture<float> > rough,
        Reference<Texture<float> > bump);
    BSDF *GetBSDF(const DifferentialGeometry &dgGeom,
        const DifferentialGeometry &dgShading, MemoryArena &arena) const;
private:
    // MetalMaterial Private Data
    Reference<Texture<Spectrum> > eta, k;
    Reference<Texture<float> > roughness;
    Reference<Texture<float> > bumpMap;
};


MetalMaterial *CreateMetalMaterial(const Transform &xform,
        const TextureParams &mp);

#endif // PBRT_MATERIALS_METAL_H
