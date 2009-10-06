
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

#ifndef PBRT_ACCELERATORS_GRID_H
#define PBRT_ACCELERATORS_GRID_H

// accelerators/grid.h*
#include "pbrt.h"
#include "primitive.h"

// GridAccel Forward Declarations
struct Voxel;

// Voxel Declarations
struct Voxel {
    // Voxel Public Methods
    u_int size() const { return primitives.size(); }
    Voxel() { }
    Voxel(Reference<Primitive> op) {
        allCanIntersect = false;
        primitives.push_back(op);
    }
    void AddPrimitive(Reference<Primitive> prim) {
        primitives.push_back(prim);
    }
    bool Intersect(const Ray &ray, Intersection *isect, RWMutexLock &lock);
    bool IntersectP(const Ray &ray, RWMutexLock &lock);
private:
    vector<Reference<Primitive> > primitives;
    bool allCanIntersect;
};



// GridAccel Declarations
class  GridAccel : public Aggregate {
public:
    // GridAccel Public Methods
    GridAccel(const vector<Reference<Primitive> > &p,
              bool forRefined, bool refineImmediately);
    BBox WorldBound() const;
    bool CanIntersect() const { return true; }
    ~GridAccel();
    bool Intersect(const Ray &ray, Intersection *isect) const;
    bool IntersectP(const Ray &ray) const;
private:
    // GridAccel Private Methods
    int posToVoxel(const Point &P, int axis) const {
        int v = Float2Int((P[axis] - bounds.pMin[axis]) *
                          InvWidth[axis]);
        return Clamp(v, 0, NVoxels[axis]-1);
    }
    float voxelToPos(int p, int axis) const {
        return bounds.pMin[axis] + p * Width[axis];
    }
    Point voxelToPos(int x, int y, int z) const {
        return bounds.pMin +
            Vector(x * Width[0], y * Width[1], z * Width[2]);
    }
    inline int offset(int x, int y, int z) const {
        return z*NVoxels[0]*NVoxels[1] + y*NVoxels[0] + x;
    }

    // GridAccel Private Data
    bool gridForRefined;
    vector<Reference<Primitive> > primitives;
    int NVoxels[3];
    BBox bounds;
    Vector Width, InvWidth;
    Voxel **voxels;
    MemoryArena voxelArena;
    mutable RWMutex *rwMutex;
};


GridAccel *CreateGridAccelerator(const vector<Reference<Primitive> > &prims,
        const ParamSet &ps);

#endif // PBRT_ACCELERATORS_GRID_H
