
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
    uint32_t size() const { return primitives.size(); }
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
    GridAccel(const vector<Reference<Primitive> > &p, bool refineImmediately);
    BBox WorldBound() const;
    bool CanIntersect() const { return true; }
    ~GridAccel();
    bool Intersect(const Ray &ray, Intersection *isect) const;
    bool IntersectP(const Ray &ray) const;
private:
    // GridAccel Private Methods
    int posToVoxel(const Point &P, int axis) const {
        int v = Float2Int((P[axis] - bounds.pMin[axis]) *
                          invWidth[axis]);
        return Clamp(v, 0, nVoxels[axis]-1);
    }
    float voxelToPos(int p, int axis) const {
        return bounds.pMin[axis] + p * width[axis];
    }
    inline int offset(int x, int y, int z) const {
        return z*nVoxels[0]*nVoxels[1] + y*nVoxels[0] + x;
    }

    // GridAccel Private Data
    vector<Reference<Primitive> > primitives;
    int nVoxels[3];
    BBox bounds;
    Vector width, invWidth;
    Voxel **voxels;
    MemoryArena voxelArena;
    mutable RWMutex *rwMutex;
};


GridAccel *CreateGridAccelerator(const vector<Reference<Primitive> > &prims,
        const ParamSet &ps);

#endif // PBRT_ACCELERATORS_GRID_H
