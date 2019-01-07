
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
    void buildTree(int nodeNum, const BBox &bounds, int nprims, int depth,
        BoundEdge *edges[3], const uint32_t* edge_count );

	// pick best sah split
	float	pickBestSplit(BoundEdge* edges[3], const uint32_t* edge_count, int tri_num, const BBox& box,
		unsigned& splitAxis, float& split_pos, bool& left);

	// evaluate sah value of a specific split
	float _sah(unsigned l, unsigned r, unsigned f, unsigned axis, float split, const BBox& box, bool& left);

    // KdTreeAccel Private Data
    int isectCost, traversalCost, maxPrims, maxDepth;
    float emptyBonus;
    vector<Reference<Primitive> > primitives;
    KdAccelNode *nodes;
    int nAllocedNodes, nextFreeNode;
    BBox bounds;
    MemoryArena arena;

	// temporary buffer for marking triangles
	unsigned char*	m_temp;
};


struct KdToDo {
    const KdAccelNode *node;
    float tmin, tmax;
};


KdTreeAccel *CreateKdTreeAccelerator(const vector<Reference<Primitive> > &prims,
        const ParamSet &ps);

#endif // PBRT_ACCELERATORS_KDTREEACCEL_H
