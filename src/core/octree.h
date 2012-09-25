
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

#ifndef PBRT_CORE_OCTREE_H
#define PBRT_CORE_OCTREE_H

// core/octree.h*
#include "pbrt.h"
#include "geometry.h"

// Octree Declarations
template <typename NodeData> struct OctNode {
    OctNode() {
        for (int i = 0; i < 8; ++i)
            children[i] = NULL;
    }
    ~OctNode() {
        for (int i = 0; i < 8; ++i)
            delete children[i];
    }
    OctNode *children[8];
    vector<NodeData> data;
};


template <typename NodeData> class Octree {
public:
    // Octree Public Methods
    Octree(const BBox &b, int md = 16)
        : maxDepth(md), bound(b) { }
    void Add(const NodeData &dataItem, const BBox &dataBound) {
        addPrivate(&root, bound, dataItem, dataBound,
                   DistanceSquared(dataBound.pMin, dataBound.pMax));
    }
    template <typename LookupProc> void Lookup(const Point &p,
                                               LookupProc &process) {
        if (!bound.Inside(p)) return;
        lookupPrivate(&root, bound, p, process);
    }
private:
    // Octree Private Methods
    void addPrivate(OctNode<NodeData> *node, const BBox &nodeBound,
        const NodeData &dataItem, const BBox &dataBound, float diag2,
        int depth = 0);
    template <typename LookupProc> bool lookupPrivate(OctNode<NodeData> *node,
            const BBox &nodeBound, const Point &P, LookupProc &process);

    // Octree Private Data
    int maxDepth;
    BBox bound;
    OctNode<NodeData> root;
};


inline BBox octreeChildBound(int child, const BBox &nodeBound,
                             const Point &pMid) {
    BBox childBound;
    childBound.pMin.x = (child & 4) ? pMid.x : nodeBound.pMin.x;
    childBound.pMax.x = (child & 4) ? nodeBound.pMax.x : pMid.x;
    childBound.pMin.y = (child & 2) ? pMid.y : nodeBound.pMin.y;
    childBound.pMax.y = (child & 2) ? nodeBound.pMax.y : pMid.y;
    childBound.pMin.z = (child & 1) ? pMid.z : nodeBound.pMin.z;
    childBound.pMax.z = (child & 1) ? nodeBound.pMax.z : pMid.z;
    return childBound;
}



// Octree Method Definitions
template <typename NodeData>
void Octree<NodeData>::addPrivate(
        OctNode<NodeData> *node, const BBox &nodeBound,
        const NodeData &dataItem, const BBox &dataBound,
        float diag2, int depth) {
    // Possibly add data item to current octree node
    if (depth == maxDepth ||
        DistanceSquared(nodeBound.pMin, nodeBound.pMax) < diag2) {
        node->data.push_back(dataItem);
        return;
    }

    // Otherwise add data item to octree children
    Point pMid = .5 * nodeBound.pMin + .5 * nodeBound.pMax;

    // Determine which children the item overlaps
    bool x[2] = { dataBound.pMin.x <= pMid.x, dataBound.pMax.x > pMid.x };
    bool y[2] = { dataBound.pMin.y <= pMid.y, dataBound.pMax.y > pMid.y };
    bool z[2] = { dataBound.pMin.z <= pMid.z, dataBound.pMax.z > pMid.z };
    bool over[8] = { x[0] & y[0] & z[0], x[0] & y[0] & z[1],
                     x[0] & y[1] & z[0], x[0] & y[1] & z[1],
                     x[1] & y[0] & z[0], x[1] & y[0] & z[1],
                     x[1] & y[1] & z[0], x[1] & y[1] & z[1] };
    for (int child = 0; child < 8; ++child) {
        if (!over[child]) continue;
        // Allocate octree node if needed and continue recursive traversal
        if (!node->children[child])
            node->children[child] = new OctNode<NodeData>;
        BBox childBound = octreeChildBound(child, nodeBound, pMid);
        addPrivate(node->children[child], childBound,
                   dataItem, dataBound, diag2, depth+1);
    }
}


template <typename NodeData> template <typename LookupProc>
bool Octree<NodeData>::lookupPrivate(OctNode<NodeData> *node,
        const BBox &nodeBound, const Point &p, LookupProc &process) {
    for (uint32_t i = 0; i < node->data.size(); ++i)
        if (!process(node->data[i]))
            return false;
    // Determine which octree child node _p_ is inside
    Point pMid = .5f * nodeBound.pMin + .5f * nodeBound.pMax;
    int child = (p.x > pMid.x ? 4 : 0) + (p.y > pMid.y ? 2 : 0) +
                (p.z > pMid.z ? 1 : 0);
    if (!node->children[child])
        return true;
    BBox childBound = octreeChildBound(child, nodeBound, pMid);
    return lookupPrivate(node->children[child], childBound, p, process);
}



#endif // PBRT_CORE_OCTREE_H
