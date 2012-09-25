
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

#ifndef PBRT_CORE_KDTREE_H
#define PBRT_CORE_KDTREE_H

// core/kdtree.h*
#include "pbrt.h"
#include "geometry.h"

// KdTree Declarations
struct KdNode {
    void init(float p, uint32_t a) {
        splitPos = p;
        splitAxis = a;
        rightChild = (1<<29)-1;
        hasLeftChild = 0;
    }
    void initLeaf() {
        splitAxis = 3;
        rightChild = (1<<29)-1;
        hasLeftChild = 0;
    }
    // KdNode Data
    float splitPos;
    uint32_t splitAxis:2;
    uint32_t hasLeftChild:1, rightChild:29;
};


template <typename NodeData> class KdTree {
public:
    // KdTree Public Methods
    KdTree(const vector<NodeData> &data);
    ~KdTree() {
        FreeAligned(nodes);
        FreeAligned(nodeData);
    }
    template <typename LookupProc> void Lookup(const Point &p,
            LookupProc &process, float &maxDistSquared) const;
private:
    // KdTree Private Methods
    void recursiveBuild(uint32_t nodeNum, int start, int end,
        const NodeData **buildNodes);
    template <typename LookupProc> void privateLookup(uint32_t nodeNum,
        const Point &p, LookupProc &process, float &maxDistSquared) const;

    // KdTree Private Data
    KdNode *nodes;
    NodeData *nodeData;
    uint32_t nNodes, nextFreeNode;
};


template <typename NodeData> struct CompareNode {
    CompareNode(int a) { axis = a; }
    int axis;
    bool operator()(const NodeData *d1, const NodeData *d2) const {
        return d1->p[axis] == d2->p[axis] ? (d1 < d2) :
                                            d1->p[axis] < d2->p[axis];
    }
};



// KdTree Method Definitions
template <typename NodeData>
KdTree<NodeData>::KdTree(const vector<NodeData> &d) {
    nNodes = d.size();
    nextFreeNode = 1;
    nodes = AllocAligned<KdNode>(nNodes);
    nodeData = AllocAligned<NodeData>(nNodes);
    vector<const NodeData *> buildNodes(nNodes, NULL);
    for (uint32_t i = 0; i < nNodes; ++i)
        buildNodes[i] = &d[i];
    // Begin the KdTree building process
    recursiveBuild(0, 0, nNodes, &buildNodes[0]);
}


template <typename NodeData> void
KdTree<NodeData>::recursiveBuild(uint32_t nodeNum, int start, int end,
        const NodeData **buildNodes) {
    // Create leaf node of kd-tree if we've reached the bottom
    if (start + 1 == end) {
        nodes[nodeNum].initLeaf();
        nodeData[nodeNum] = *buildNodes[start];
        return;
    }

    // Choose split direction and partition data

    // Compute bounds of data from _start_ to _end_
    BBox bound;
    for (int i = start; i < end; ++i)
        bound = Union(bound, buildNodes[i]->p);
    int splitAxis = bound.MaximumExtent();
    int splitPos = (start+end)/2;
    std::nth_element(&buildNodes[start], &buildNodes[splitPos],
                     &buildNodes[end], CompareNode<NodeData>(splitAxis));

    // Allocate kd-tree node and continue recursively
    nodes[nodeNum].init(buildNodes[splitPos]->p[splitAxis], splitAxis);
    nodeData[nodeNum] = *buildNodes[splitPos];
    if (start < splitPos) {
        nodes[nodeNum].hasLeftChild = 1;
        uint32_t childNum = nextFreeNode++;
        recursiveBuild(childNum, start, splitPos, buildNodes);
    }
    if (splitPos+1 < end) {
        nodes[nodeNum].rightChild = nextFreeNode++;
        recursiveBuild(nodes[nodeNum].rightChild, splitPos+1,
                       end, buildNodes);
    }
}


template <typename NodeData> template <typename LookupProc>
void KdTree<NodeData>::Lookup(const Point &p, LookupProc &proc,
                              float &maxDistSquared) const {
    privateLookup(0, p, proc, maxDistSquared);
}


template <typename NodeData> template <typename LookupProc>
void KdTree<NodeData>::privateLookup(uint32_t nodeNum, const Point &p,
        LookupProc &process, float &maxDistSquared) const {
    KdNode *node = &nodes[nodeNum];
    // Process kd-tree node's children
    int axis = node->splitAxis;
    if (axis != 3) {
        float dist2 = (p[axis] - node->splitPos) * (p[axis] - node->splitPos);
        if (p[axis] <= node->splitPos) {
            if (node->hasLeftChild)
                privateLookup(nodeNum+1, p, process, maxDistSquared);
            if (dist2 < maxDistSquared && node->rightChild < nNodes)
                privateLookup(node->rightChild, p, process, maxDistSquared);
        }
        else {
            if (node->rightChild < nNodes)
                privateLookup(node->rightChild, p, process, maxDistSquared);
            if (dist2 < maxDistSquared && node->hasLeftChild)
                privateLookup(nodeNum+1, p, process, maxDistSquared);
        }
    }

    // Hand kd-tree node to processing function
    float dist2 = DistanceSquared(nodeData[nodeNum].p, p);
    if (dist2 < maxDistSquared)
        process(p, nodeData[nodeNum], dist2, maxDistSquared);
}



#endif // PBRT_CORE_KDTREE_H
