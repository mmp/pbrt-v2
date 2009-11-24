
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

#ifndef PBRT_CORE_LIGHT_H
#define PBRT_CORE_LIGHT_H

// core/light.h*
#include "pbrt.h"
#include "geometry.h"
#include "transform.h"
#include "spectrum.h"
#include "rng.h"
#include "memory.h"

// Light Declarations
class Light {
public:
    // Light Interface
    virtual ~Light();
    Light(const Transform &l2w, int ns = 1)
        : nSamples(max(1, ns)), LightToWorld(l2w),
          WorldToLight(Inverse(l2w)) {
        // Warn if light has transformation with scale
        if (WorldToLight.HasScale())
            Warning("Scaling detected in world to light transformation!\n"
                    "The system has numerous assumptions, implicit and explicit,\n"
                    "that this transform will have no scale factors in it.\n"
                    "Proceed at your own risk; your image may have errors or\n"
                    "the system may crash as a result of this.");
    }
    virtual Spectrum Sample_L(const Point &p, float pEpsilon,
        const LightSample &ls, float time, Vector *wi, float *pdf,
        VisibilityTester *vis) const = 0;
    virtual Spectrum Power(const Scene *) const = 0;
    virtual bool IsDeltaLight() const = 0;
    virtual Spectrum Le(const RayDifferential &r) const;
    virtual float Pdf(const Point &p, const Vector &wi) const = 0;
    virtual Spectrum Sample_L(const Scene *scene, const LightSample &ls,
                              float u1, float u2, float time, Ray *ray,
                              Normal *Ns, float *pdf) const = 0;
    virtual void SHProject(const Point &p, float pEpsilon, int lmax,
        const Scene *scene, bool computeLightVisibility, float time,
        RNG &rng, Spectrum *coeffs) const;

    // Light Public Data
    const int nSamples;
protected:
    // Light Protected Data
    const Transform LightToWorld, WorldToLight;
};


struct VisibilityTester {
    // VisibilityTester Public Methods
    void SetSegment(const Point &p1, float eps1,
                    const Point &p2, float eps2, float time) {
        float dist = Distance(p1, p2);
        r = Ray(p1, (p2-p1) / dist, eps1, dist * (1.f - eps2), time);
        Assert(!r.HasNaNs());
    }
    void SetRay(const Point &p, float eps, const Vector &w, float time) {
        r = Ray(p, w, eps, INFINITY, time);
        Assert(!r.HasNaNs());
    }
    bool Unoccluded(const Scene *scene) const;
    Spectrum Transmittance(const Scene *scene, const Renderer *renderer,
        const Sample *sample, RNG &rng, MemoryArena &arena) const;
    Ray r;
};


class AreaLight : public Light {
public:
    // AreaLight Interface
    AreaLight(const Transform &l2w, int ns) : Light(l2w, ns) { }
    virtual Spectrum L(const Point &p, const Normal &n,
                       const Vector &w) const = 0;
};


struct LightSample {
   // LightSample Public Methods
   LightSample() { }
   LightSample(const Sample *sample, const LightSampleOffsets &offsets, uint32_t num);
   LightSample(RNG &rng) {
       uPos[0] = rng.RandomFloat();
       uPos[1] = rng.RandomFloat();
       uComponent = rng.RandomFloat();
   }
   LightSample(float up0, float up1, float ucomp) {
       uPos[0] = up0; uPos[1] = up1;
       uComponent = ucomp;
   }
   float uPos[2], uComponent;
};


struct LightSampleOffsets {
    LightSampleOffsets() { }
    LightSampleOffsets(int count, Sample *sample);
    int nSamples, componentOffset, posOffset;
};



// ShapeSet Declarations
class ShapeSet {
public:
    // ShapeSet Public Methods
    ShapeSet(const Reference<Shape> &s);
    float Area() const { return sumArea; }
    ~ShapeSet();
    Point Sample(const Point &p, const LightSample &ls, Normal *Ns) const;
    Point Sample(const LightSample &ls, Normal *Ns) const;
    float Pdf(const Point &p, const Vector &wi) const;
    float Pdf(const Point &p) const;
private:
    // ShapeSet Private Data
    vector<Reference<Shape> > shapes;
    float sumArea;
    vector<float> areas;
    Distribution1D *areaDistribution;
};



#endif // PBRT_CORE_LIGHT_H
