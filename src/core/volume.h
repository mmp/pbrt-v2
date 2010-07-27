
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

#ifndef PBRT_CORE_VOLUME_H
#define PBRT_CORE_VOLUME_H

// core/volume.h*
#include "pbrt.h"
#include "spectrum.h"
#include "geometry.h"
#include "transform.h"
#include "integrator.h"

// Volume Scattering Declarations
float PhaseIsotropic(const Vector &w, const Vector &wp);
float PhaseRayleigh(const Vector &w, const Vector &wp);
float PhaseMieHazy(const Vector &w, const Vector &wp);
float PhaseMieMurky(const Vector &w, const Vector &wp);
float PhaseHG(const Vector &w, const Vector &wp, float g);
float PhaseSchlick(const Vector &w, const Vector &wp, float g);
class VolumeRegion {
public:
    // VolumeRegion Interface
    virtual ~VolumeRegion();
    virtual BBox WorldBound() const = 0;
    virtual bool IntersectP(const Ray &ray, float *t0, float *t1) const = 0;
    virtual Spectrum sigma_a(const Point &, const Vector &,
                             float time) const = 0;
    virtual Spectrum sigma_s(const Point &, const Vector &,
                             float time) const = 0;
    virtual Spectrum Lve(const Point &, const Vector &,
                         float time) const = 0;
    virtual float p(const Point &, const Vector &,
                    const Vector &, float time) const = 0;
    virtual Spectrum sigma_t(const Point &p, const Vector &wo, float time) const;
    virtual Spectrum tau(const Ray &ray, float step = 1.f,
                         float offset = 0.5) const = 0;
};


class DensityRegion : public VolumeRegion {
public:
    // DensityRegion Public Methods
    DensityRegion(const Spectrum &sa, const Spectrum &ss, float gg,
                  const Spectrum &emit, const Transform &VolumeToWorld)
        : sig_a(sa), sig_s(ss), le(emit), g(gg),
          WorldToVolume(Inverse(VolumeToWorld)) { }
    virtual float Density(const Point &Pobj) const = 0;
    Spectrum sigma_a(const Point &p, const Vector &, float) const {
        return Density(WorldToVolume(p)) * sig_a;
    }
    Spectrum sigma_s(const Point &p, const Vector &, float) const {
        return Density(WorldToVolume(p)) * sig_s;
    }
    Spectrum sigma_t(const Point &p, const Vector &, float) const {
        return Density(WorldToVolume(p)) * (sig_a + sig_s);
    }
    Spectrum Lve(const Point &p, const Vector &, float) const {
        return Density(WorldToVolume(p)) * le;
    }
    float p(const Point &p, const Vector &w, const Vector &wp, float) const {
        return PhaseHG(w, wp, g);
    }
    Spectrum tau(const Ray &r, float stepSize, float offset) const;
protected:
    // DensityRegion Protected Data
    Spectrum sig_a, sig_s, le;
    float g;
    Transform WorldToVolume;
};


class AggregateVolume : public VolumeRegion {
public:
    // AggregateVolume Public Methods
    AggregateVolume(const vector<VolumeRegion *> &r);
    ~AggregateVolume();
    BBox WorldBound() const;
    bool IntersectP(const Ray &ray, float *t0, float *t1) const;
    Spectrum sigma_a(const Point &, const Vector &, float) const;
    Spectrum sigma_s(const Point &, const Vector &, float) const;
    Spectrum Lve(const Point &, const Vector &, float) const;
    float p(const Point &, const Vector &, const Vector &, float) const;
    Spectrum sigma_t(const Point &, const Vector &, float) const;
    Spectrum tau(const Ray &ray, float, float) const;
private:
    // AggregateVolume Private Data
    vector<VolumeRegion *> regions;
    BBox bound;
};


bool GetVolumeScatteringProperties(const string &name, Spectrum *sigma_a,
                                   Spectrum *sigma_prime_s);
class VolumeIntegrator : public Integrator {
public:
    // VolumeIntegrator Interface
    virtual Spectrum Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential &ray, const Sample *sample, RNG &rng,
        Spectrum *transmittance, MemoryArena &arena) const = 0;
    virtual Spectrum Transmittance(const Scene *scene,
        const Renderer *renderer, const RayDifferential &ray,
        const Sample *sample, RNG &rng, MemoryArena &arena) const = 0;
};


void SubsurfaceFromDiffuse(const Spectrum &Kd, float meanPathLength, float eta,
        Spectrum *sigma_a, Spectrum *sigma_prime_s);

#endif // PBRT_CORE_VOLUME_H
