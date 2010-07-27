
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
