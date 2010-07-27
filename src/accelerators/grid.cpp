
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


// accelerators/grid.cpp*
#include "stdafx.h"
#include "accelerators/grid.h"
#include "probes.h"
#include "paramset.h"

// GridAccel Method Definitions
GridAccel::GridAccel(const vector<Reference<Primitive> > &p,
                     bool refineImmediately) {
    PBRT_GRID_STARTED_CONSTRUCTION(this, p.size());
    // Initialize _primitives_ with primitives for grid
    if (refineImmediately)
        for (uint32_t i = 0; i < p.size(); ++i)
            p[i]->FullyRefine(primitives);
    else
        primitives = p;

    // Compute bounds and choose grid resolution
    for (uint32_t i = 0; i < primitives.size(); ++i)
        bounds = Union(bounds, primitives[i]->WorldBound());
    Vector delta = bounds.pMax - bounds.pMin;

    // Find _voxelsPerUnitDist_ for grid
    int maxAxis = bounds.MaximumExtent();
    float invMaxWidth = 1.f / delta[maxAxis];
    Assert(invMaxWidth > 0.f);
    float cubeRoot = 3.f * powf(float(primitives.size()), 1.f/3.f);
    float voxelsPerUnitDist = cubeRoot * invMaxWidth;
    for (int axis = 0; axis < 3; ++axis) {
        nVoxels[axis] = Round2Int(delta[axis] * voxelsPerUnitDist);
        nVoxels[axis] = Clamp(nVoxels[axis], 1, 64);
    }
    PBRT_GRID_BOUNDS_AND_RESOLUTION(&bounds, nVoxels);

    // Compute voxel widths and allocate voxels
    for (int axis = 0; axis < 3; ++axis) {
        width[axis] = delta[axis] / nVoxels[axis];
        invWidth[axis] = (width[axis] == 0.f) ? 0.f : 1.f / width[axis];
    }
    int nv = nVoxels[0] * nVoxels[1] * nVoxels[2];
    voxels = AllocAligned<Voxel *>(nv);
    memset(voxels, 0, nv * sizeof(Voxel *));

    // Add primitives to grid voxels
    for (uint32_t i = 0; i < primitives.size(); ++i) {
        // Find voxel extent of primitive
        BBox pb = primitives[i]->WorldBound();
        int vmin[3], vmax[3];
        for (int axis = 0; axis < 3; ++axis) {
            vmin[axis] = posToVoxel(pb.pMin, axis);
            vmax[axis] = posToVoxel(pb.pMax, axis);
        }

        // Add primitive to overlapping voxels
        PBRT_GRID_VOXELIZED_PRIMITIVE(vmin, vmax);
        for (int z = vmin[2]; z <= vmax[2]; ++z)
            for (int y = vmin[1]; y <= vmax[1]; ++y)
                for (int x = vmin[0]; x <= vmax[0]; ++x) {
                    int o = offset(x, y, z);
                    if (!voxels[o]) {
                        // Allocate new voxel and store primitive in it
                        voxels[o] = voxelArena.Alloc<Voxel>();
                        *voxels[o] = Voxel(primitives[i]);
                    }
                    else {
                        // Add primitive to already-allocated voxel
                        voxels[o]->AddPrimitive(primitives[i]);
                    }
                }
    }

    // Create reader-writer mutex for grid
    rwMutex = RWMutex::Create();
    PBRT_GRID_FINISHED_CONSTRUCTION(this);
}


BBox GridAccel::WorldBound() const {
    return bounds;
}


GridAccel::~GridAccel() {
    for (int i = 0; i < nVoxels[0]*nVoxels[1]*nVoxels[2]; ++i)
        if (voxels[i]) voxels[i]->~Voxel();
    FreeAligned(voxels);
    RWMutex::Destroy(rwMutex);
}


bool GridAccel::Intersect(const Ray &ray, Intersection *isect) const {
    PBRT_GRID_INTERSECTION_TEST(const_cast<GridAccel *>(this), const_cast<Ray *>(&ray));
    // Check ray against overall grid bounds
    float rayT;
    if (bounds.Inside(ray(ray.mint)))
        rayT = ray.mint;
    else if (!bounds.IntersectP(ray, &rayT))
    {
        PBRT_GRID_RAY_MISSED_BOUNDS();
        return false;
    }
    Point gridIntersect = ray(rayT);

    // Set up 3D DDA for ray
    float NextCrossingT[3], DeltaT[3];
    int Step[3], Out[3], Pos[3];
    for (int axis = 0; axis < 3; ++axis) {
        // Compute current voxel for axis
        Pos[axis] = posToVoxel(gridIntersect, axis);
        if (ray.d[axis] >= 0) {
            // Handle ray with positive direction for voxel stepping
            NextCrossingT[axis] = rayT +
                (voxelToPos(Pos[axis]+1, axis) - gridIntersect[axis]) / ray.d[axis];
            DeltaT[axis] = width[axis] / ray.d[axis];
            Step[axis] = 1;
            Out[axis] = nVoxels[axis];
        }
        else {
            // Handle ray with negative direction for voxel stepping
            NextCrossingT[axis] = rayT +
                (voxelToPos(Pos[axis], axis) - gridIntersect[axis]) / ray.d[axis];
            DeltaT[axis] = -width[axis] / ray.d[axis];
            Step[axis] = -1;
            Out[axis] = -1;
        }
    }

    // Walk ray through voxel grid
    RWMutexLock lock(*rwMutex, READ);
    bool hitSomething = false;
    for (;;) {
        // Check for intersection in current voxel and advance to next
        Voxel *voxel = voxels[offset(Pos[0], Pos[1], Pos[2])];
        PBRT_GRID_RAY_TRAVERSED_VOXEL(Pos, voxel ? voxel->size() : 0);
        if (voxel != NULL)
            hitSomething |= voxel->Intersect(ray, isect, lock);

        // Advance to next voxel

        // Find _stepAxis_ for stepping to next voxel
        int bits = ((NextCrossingT[0] < NextCrossingT[1]) << 2) +
                   ((NextCrossingT[0] < NextCrossingT[2]) << 1) +
                   ((NextCrossingT[1] < NextCrossingT[2]));
        const int cmpToAxis[8] = { 2, 1, 2, 1, 2, 2, 0, 0 };
        int stepAxis = cmpToAxis[bits];
        if (ray.maxt < NextCrossingT[stepAxis])
            break;
        Pos[stepAxis] += Step[stepAxis];
        if (Pos[stepAxis] == Out[stepAxis])
            break;
        NextCrossingT[stepAxis] += DeltaT[stepAxis];
    }
    return hitSomething;
}


bool Voxel::Intersect(const Ray &ray, Intersection *isect,
                      RWMutexLock &lock) {
    // Refine primitives in voxel if needed
    if (!allCanIntersect) {
        lock.UpgradeToWrite();
        for (uint32_t i = 0; i < primitives.size(); ++i) {
            Reference<Primitive> &prim = primitives[i];
            // Refine primitive _prim_ if it's not intersectable
            if (!prim->CanIntersect()) {
                vector<Reference<Primitive> > p;
                prim->FullyRefine(p);
                Assert(p.size() > 0);
                if (p.size() == 1)
                    primitives[i] = p[0];
                else
                    primitives[i] = new GridAccel(p, false);
            }
        }
        allCanIntersect = true;
        lock.DowngradeToRead();
    }

    // Loop over primitives in voxel and find intersections
    bool hitSomething = false;
    for (uint32_t i = 0; i < primitives.size(); ++i) {
        Reference<Primitive> &prim = primitives[i];
        PBRT_GRID_RAY_PRIMITIVE_INTERSECTION_TEST(const_cast<Primitive *>(prim.GetPtr()));
        if (prim->Intersect(ray, isect))
        {
        PBRT_GRID_RAY_PRIMITIVE_HIT(const_cast<Primitive *>(prim.GetPtr()));
            hitSomething = true;
        }
    }
    return hitSomething;
}


bool GridAccel::IntersectP(const Ray &ray) const {
    PBRT_GRID_INTERSECTIONP_TEST(const_cast<GridAccel *>(this), const_cast<Ray *>(&ray));
    RWMutexLock lock(*rwMutex, READ);
    // Check ray against overall grid bounds
    float rayT;
    if (bounds.Inside(ray(ray.mint)))
        rayT = ray.mint;
    else if (!bounds.IntersectP(ray, &rayT))
    {
        PBRT_GRID_RAY_MISSED_BOUNDS();
        return false;
    }
    Point gridIntersect = ray(rayT);

    // Set up 3D DDA for ray
    float NextCrossingT[3], DeltaT[3];
    int Step[3], Out[3], Pos[3];
    for (int axis = 0; axis < 3; ++axis) {
        // Compute current voxel for axis
        Pos[axis] = posToVoxel(gridIntersect, axis);
        if (ray.d[axis] >= 0) {
            // Handle ray with positive direction for voxel stepping
            NextCrossingT[axis] = rayT +
                (voxelToPos(Pos[axis]+1, axis) - gridIntersect[axis]) / ray.d[axis];
            DeltaT[axis] = width[axis] / ray.d[axis];
            Step[axis] = 1;
            Out[axis] = nVoxels[axis];
        }
        else {
            // Handle ray with negative direction for voxel stepping
            NextCrossingT[axis] = rayT +
                (voxelToPos(Pos[axis], axis) - gridIntersect[axis]) / ray.d[axis];
            DeltaT[axis] = -width[axis] / ray.d[axis];
            Step[axis] = -1;
            Out[axis] = -1;
        }
    }

    // Walk grid for shadow ray
    for (;;) {
        int o = offset(Pos[0], Pos[1], Pos[2]);
        Voxel *voxel = voxels[o];
        PBRT_GRID_RAY_TRAVERSED_VOXEL(Pos, voxel ? voxel->size() : 0);
        if (voxel && voxel->IntersectP(ray, lock))
            return true;
        // Advance to next voxel

        // Find _stepAxis_ for stepping to next voxel
        int bits = ((NextCrossingT[0] < NextCrossingT[1]) << 2) +
                   ((NextCrossingT[0] < NextCrossingT[2]) << 1) +
                   ((NextCrossingT[1] < NextCrossingT[2]));
        const int cmpToAxis[8] = { 2, 1, 2, 1, 2, 2, 0, 0 };
        int stepAxis = cmpToAxis[bits];
        if (ray.maxt < NextCrossingT[stepAxis])
            break;
        Pos[stepAxis] += Step[stepAxis];
        if (Pos[stepAxis] == Out[stepAxis])
            break;
        NextCrossingT[stepAxis] += DeltaT[stepAxis];
    }
    return false;
}


bool Voxel::IntersectP(const Ray &ray, RWMutexLock &lock) {
    // Refine primitives in voxel if needed
    if (!allCanIntersect) {
        lock.UpgradeToWrite();
        for (uint32_t i = 0; i < primitives.size(); ++i) {
            Reference<Primitive> &prim = primitives[i];
            // Refine primitive _prim_ if it's not intersectable
            if (!prim->CanIntersect()) {
                vector<Reference<Primitive> > p;
                prim->FullyRefine(p);
                Assert(p.size() > 0);
                if (p.size() == 1)
                    primitives[i] = p[0];
                else
                    primitives[i] = new GridAccel(p, false);
            }
        }
        allCanIntersect = true;
        lock.DowngradeToRead();
    }
    for (uint32_t i = 0; i < primitives.size(); ++i) {
        Reference<Primitive> &prim = primitives[i];
        PBRT_GRID_RAY_PRIMITIVE_INTERSECTIONP_TEST(const_cast<Primitive *>(prim.GetPtr()));
        if (prim->IntersectP(ray)) {
            PBRT_GRID_RAY_PRIMITIVE_HIT(const_cast<Primitive *>(prim.GetPtr()));
            return true;
        }
    }
    return false;
}


GridAccel *CreateGridAccelerator(const vector<Reference<Primitive> > &prims,
        const ParamSet &ps) {
    bool refineImmediately = ps.FindOneBool("refineimmediately", false);
    return new GridAccel(prims, refineImmediately);
}


