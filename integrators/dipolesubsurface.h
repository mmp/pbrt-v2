
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

#ifndef PBRT_INTEGRATORS_DIPOLESUBSURFACE_H
#define PBRT_INTEGRATORS_DIPOLESUBSURFACE_H

// integrators/dipolesubsurface.h*
#include "pbrt.h"
#include "integrator.h"
#include "kdtree.h"
struct SubsurfaceOctreeNode;

// DipoleSubsurfaceIntegrator Helper Declarations
struct IrradiancePoint {
    Point p;
    Normal n;
    Spectrum E;
    float area, rayEpsilon;
};



// DipoleSubsurfaceIntegrator Declarations
class DipoleSubsurfaceIntegrator : public SurfaceIntegrator {
public:
    // DipoleSubsurfaceIntegrator Public Methods
    ~DipoleSubsurfaceIntegrator();
    Spectrum Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential &ray, const Intersection &isect, const Sample *sample,
        RNG &rng, MemoryArena &arena) const;
    void RequestSamples(Sampler *sampler, Sample *sample, const Scene *scene);
    void Preprocess(const Scene *, const Camera *, const Renderer *);
    DipoleSubsurfaceIntegrator(int mdepth, float merror, float mindist) {
        maxSpecularDepth = mdepth;
        maxError = merror;
        minSampleDist = mindist;
        octree = NULL;
    }
private:
    // DipoleSubsurfaceIntegrator Private Data
    int maxSpecularDepth;
    float maxError, minSampleDist;
    vector<IrradiancePoint> irradiancePoints;
    BBox octreeBounds;
    SubsurfaceOctreeNode *octree;
    MemoryArena octreeArena;

    // Declare sample parameters for light source sampling
    LightSampleOffsets *lightSampleOffsets;
    BSDFSampleOffsets *bsdfSampleOffsets;
};


DipoleSubsurfaceIntegrator *CreateDipoleSubsurfaceIntegrator(const ParamSet &params);

#endif // PBRT_INTEGRATORS_DIPOLESUBSURFACE_H
