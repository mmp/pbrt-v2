
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
