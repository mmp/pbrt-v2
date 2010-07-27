
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

#ifndef PBRT_MATERIALS_GLASS_H
#define PBRT_MATERIALS_GLASS_H

// materials/glass.h*
#include "pbrt.h"
#include "material.h"

// GlassMaterial Declarations
class GlassMaterial : public Material {
public:
    // GlassMaterial Public Methods
    GlassMaterial(Reference<Texture<Spectrum> > r, Reference<Texture<Spectrum> > t,
            Reference<Texture<float> > i, Reference<Texture<float> > bump) {
        Kr = r;
        Kt = t;
        index = i;
        bumpMap = bump;
    }
    BSDF *GetBSDF(const DifferentialGeometry &dgGeom, const DifferentialGeometry &dgShading, MemoryArena &arena) const;
private:
    // GlassMaterial Private Data
    Reference<Texture<Spectrum> > Kr, Kt;
    Reference<Texture<float> > index;
    Reference<Texture<float> > bumpMap;
};


GlassMaterial *CreateGlassMaterial(const Transform &xform,
        const TextureParams &mp);

#endif // PBRT_MATERIALS_GLASS_H
