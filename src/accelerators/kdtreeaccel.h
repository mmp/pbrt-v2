
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

#ifndef PBRT_ACCELERATORS_KDTREEACCEL_H
#define PBRT_ACCELERATORS_KDTREEACCEL_H

// accelerators/kdtreeaccel.h*
#include "pbrt.h"
#include "primitive.h"

// KdTreeAccel Declarations
struct KdAccelNode;
struct BoundEdge;
class KdTreeAccel : public Aggregate {
public:
    // KdTreeAccel Public Methods
    KdTreeAccel(const vector<Reference<Primitive> > &p,
                int icost = 80, int scost = 1,  float ebonus = 0.5f, int maxp = 1,
                int maxDepth = -1);
    BBox WorldBound() const { return bounds; }
    bool CanIntersect() const { return true; }
    ~KdTreeAccel();
    bool Intersect(const Ray &ray, Intersection *isect) const;
    bool IntersectP(const Ray &ray) const;
private:
    // KdTreeAccel Private Methods
    void buildTree(int nodeNum, const BBox &bounds,
        const vector<BBox> &primBounds, uint32_t *primNums, int nprims, int depth,
        BoundEdge *edges[3], uint32_t *prims0, uint32_t *prims1, int badRefines = 0);

    // KdTreeAccel Private Data
    int isectCost, traversalCost, maxPrims, maxDepth;
    float emptyBonus;
    vector<Reference<Primitive> > primitives;
    KdAccelNode *nodes;
    int nAllocedNodes, nextFreeNode;
    BBox bounds;
    MemoryArena arena;
};


struct KdToDo {
    const KdAccelNode *node;
    float tmin, tmax;
};


KdTreeAccel *CreateKdTreeAccelerator(const vector<Reference<Primitive> > &prims,
        const ParamSet &ps);

#endif // PBRT_ACCELERATORS_KDTREEACCEL_H
