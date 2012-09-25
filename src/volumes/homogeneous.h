
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

#ifndef PBRT_VOLUMES_HOMOGENEOUS_H
#define PBRT_VOLUMES_HOMOGENEOUS_H

// volumes/homogeneous.h*
#include "volume.h"

// HomogeneousVolumeDensity Declarations
class HomogeneousVolumeDensity : public VolumeRegion {
public:
    // HomogeneousVolumeDensity Public Methods
    HomogeneousVolumeDensity(const Spectrum &sa, const Spectrum &ss, float gg,
            const Spectrum &emit, const BBox &e, const Transform &v2w) {
        WorldToVolume = Inverse(v2w);
        sig_a = sa;
        sig_s = ss;
        g = gg;
        le = emit;
        extent = e;
    }
    BBox WorldBound() const {
        return Inverse(WorldToVolume)(extent);
    }
    bool IntersectP(const Ray &r, float *t0, float *t1) const {
        Ray ray = WorldToVolume(r);
        return extent.IntersectP(ray, t0, t1);
    }
    Spectrum sigma_a(const Point &p, const Vector &, float) const {
        return extent.Inside(WorldToVolume(p)) ? sig_a : 0.;
    }
    Spectrum sigma_s(const Point &p, const Vector &, float) const {
        return extent.Inside(WorldToVolume(p)) ? sig_s : 0.;
    }
    Spectrum sigma_t(const Point &p, const Vector &, float) const {
        return extent.Inside(WorldToVolume(p)) ? (sig_a + sig_s) : 0.;
    }
    Spectrum Lve(const Point &p, const Vector &, float) const {
        return extent.Inside(WorldToVolume(p)) ? le : 0.;
    }
    float p(const Point &p, const Vector &wi, const Vector &wo, float) const {
        if (!extent.Inside(WorldToVolume(p))) return 0.;
        return PhaseHG(wi, wo, g);
    }
    Spectrum tau(const Ray &ray, float, float) const {
        float t0, t1;
        if (!IntersectP(ray, &t0, &t1)) return 0.;
        return Distance(ray(t0), ray(t1)) * (sig_a + sig_s);
    }
private:
    // HomogeneousVolumeDensity Private Data
    Spectrum sig_a, sig_s, le;
    float g;
    BBox extent;
    Transform WorldToVolume;
};


HomogeneousVolumeDensity *CreateHomogeneousVolumeDensityRegion(const Transform &volume2world,
        const ParamSet &params);

#endif // PBRT_VOLUMES_HOMOGENEOUS_H
