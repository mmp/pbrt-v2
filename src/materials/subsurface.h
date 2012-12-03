
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

#if defined(_MSC_VER)
#pragma once
#endif

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
