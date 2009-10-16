
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

#ifndef PBRT_MATERIALS_SUBSURFACE_H
#define PBRT_MATERIALS_SUBSURFACE_H

// materials/subsurface.h*
#include "pbrt.h"
#include "material.h"

// SubsurfaceMaterial Declarations
class SubsurfaceMaterial : public Material {
public:
    // SubsurfaceMaterial Public Methods
    SubsurfaceMaterial(float sc, Reference<Texture<Spectrum> > kr,
            Reference<Texture<Spectrum> > sa,
            Reference<Texture<Spectrum> > sps,
            Reference<Texture<float> > e,
            Reference<Texture<float> > bump) {
        scale = sc;
        Kr = kr;
        sigma_a = sa;
        sigma_prime_s = sps;
        eta = e;
        bumpMap = bump;
    }
    BSDF *GetBSDF(const DifferentialGeometry &dgGeom,
                  const DifferentialGeometry &dgShading,
                  MemoryArena &arena) const;
    BSSRDF *GetBSSRDF(const DifferentialGeometry &dgGeom,
                      const DifferentialGeometry &dgShading,
                      MemoryArena &arena) const;
private:
    // SubsurfaceMaterial Private Data
    float scale;
    Reference<Texture<Spectrum> > Kr, sigma_a, sigma_prime_s;
    Reference<Texture<float> > eta, bumpMap;
};


SubsurfaceMaterial *CreateSubsurfaceMaterial(const Transform &xform,
        const TextureParams &mp);

#endif // PBRT_MATERIALS_SUBSURFACE_H
