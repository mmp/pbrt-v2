
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

#ifndef PBRT_VOLUMES_EXPONENTIAL_H
#define PBRT_VOLUMES_EXPONENTIAL_H

// volumes/exponential.h*
#include "volume.h"

// ExponentialDensity Declarations
class ExponentialDensity : public DensityRegion {
public:
    // ExponentialDensity Public Methods
    ExponentialDensity(const Spectrum &sa, const Spectrum &ss,
                       float gg, const Spectrum &emit, const BBox &e,
                       const Transform &v2w, float aa, float bb,
                       const Vector &up)
        : DensityRegion(sa, ss, gg, emit, v2w), extent(e), a(aa), b(bb) {
        upDir = Normalize(up);
    }
    BBox WorldBound() const { return Inverse(WorldToVolume)(extent); }
    bool IntersectP(const Ray &r, float *t0, float *t1) const {
        Ray ray = WorldToVolume(r);
        return extent.IntersectP(ray, t0, t1);
    }
    float Density(const Point &Pobj) const {
        if (!extent.Inside(Pobj)) return 0;
        float height = Dot(Pobj - extent.pMin, upDir);
        return a * expf(-b * height);
    }
private:
    // ExponentialDensity Private Data
    BBox extent;
    float a, b;
    Vector upDir;
};


ExponentialDensity *CreateExponentialVolumeRegion(const Transform &volume2world,
        const ParamSet &params);

#endif // PBRT_VOLUMES_EXPONENTIAL_H
