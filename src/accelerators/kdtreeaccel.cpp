
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


// accelerators/kdtreeaccel.cpp*
#include "stdafx.h"
#include "accelerators/kdtreeaccel.h"
#include "paramset.h"

// KdTreeAccel Local Declarations
struct KdAccelNode {
    // KdAccelNode Methods
    void initLeaf(uint32_t *primNums, int np, MemoryArena &arena);
    void initInterior(uint32_t axis, uint32_t ac, float s) {
        split = s;
        flags = axis;
        aboveChild |= (ac << 2);
    }
    float SplitPos() const { return split; }
    uint32_t nPrimitives() const { return nPrims >> 2; }
    uint32_t SplitAxis() const { return flags & 3; }
    bool IsLeaf() const { return (flags & 3) == 3; }
    uint32_t AboveChild() const { return aboveChild >> 2; }
    union {
        float split;            // Interior
        uint32_t onePrimitive;  // Leaf
        uint32_t *primitives;   // Leaf
    };

private:
    union {
        uint32_t flags;         // Both
        uint32_t nPrims;        // Leaf
        uint32_t aboveChild;    // Interior
    };
};


struct BoundEdge {
    // BoundEdge Public Methods
    BoundEdge() { }
    BoundEdge(float tt, int pn, bool starting) {
        t = tt;
        primNum = pn;
        type = starting ? START : END;
    }
    bool operator<(const BoundEdge &e) const {
        if (t == e.t)
            return (int)type < (int)e.type;
        else return t < e.t;
    }
    float t;
    int primNum;
    enum { START, END } type;
};



// KdTreeAccel Method Definitions
KdTreeAccel::KdTreeAccel(const vector<Reference<Primitive> > &p,
                         int icost, int tcost, float ebonus, int maxp,
                         int md)
    : isectCost(icost), traversalCost(tcost), maxPrims(maxp), maxDepth(md),
      emptyBonus(ebonus) {
    PBRT_KDTREE_STARTED_CONSTRUCTION(this, p.size());
    for (uint32_t i = 0; i < p.size(); ++i)
        p[i]->FullyRefine(primitives);
    // Build kd-tree for accelerator
    nextFreeNode = nAllocedNodes = 0;
    if (maxDepth <= 0)
        maxDepth = Round2Int(8 + 1.3f * Log2Int(float(primitives.size())));

    // Compute bounds for kd-tree construction
    vector<BBox> primBounds;
    primBounds.reserve(primitives.size());
    for (uint32_t i = 0; i < primitives.size(); ++i) {
        BBox b = primitives[i]->WorldBound();
        bounds = Union(bounds, b);
        primBounds.push_back(b);
    }

    // Allocate working memory for kd-tree construction
    BoundEdge *edges[3];
    for (int i = 0; i < 3; ++i)
        edges[i] = new BoundEdge[2*primitives.size()];
    uint32_t *prims0 = new uint32_t[primitives.size()];
    uint32_t *prims1 = new uint32_t[(maxDepth+1) * primitives.size()];

    // Initialize _primNums_ for kd-tree construction
    uint32_t *primNums = new uint32_t[primitives.size()];
    for (uint32_t i = 0; i < primitives.size(); ++i)
        primNums[i] = i;

    // Start recursive construction of kd-tree
    buildTree(0, bounds, primBounds, primNums, primitives.size(),
              maxDepth, edges, prims0, prims1);

    // Free working memory for kd-tree construction
    delete[] primNums;
    for (int i = 0; i < 3; ++i)
        delete[] edges[i];
    delete[] prims0;
    delete[] prims1;
    PBRT_KDTREE_FINISHED_CONSTRUCTION(this);
}


void KdAccelNode::initLeaf(uint32_t *primNums, int np,
                           MemoryArena &arena) {
    flags = 3;
    nPrims |= (np << 2);
    // Store primitive ids for leaf node
    if (np == 0)
        onePrimitive = 0;
    else if (np == 1)
        onePrimitive = primNums[0];
    else {
        primitives = arena.Alloc<uint32_t>(np);
        for (int i = 0; i < np; ++i)
            primitives[i] = primNums[i];
    }
}


KdTreeAccel::~KdTreeAccel() {
    FreeAligned(nodes);
}


void KdTreeAccel::buildTree(int nodeNum, const BBox &nodeBounds,
        const vector<BBox> &allPrimBounds, uint32_t *primNums,
        int nPrimitives, int depth, BoundEdge *edges[3],
        uint32_t *prims0, uint32_t *prims1, int badRefines) {
    Assert(nodeNum == nextFreeNode);
    // Get next free node from _nodes_ array
    if (nextFreeNode == nAllocedNodes) {
        int nAlloc = max(2 * nAllocedNodes, 512);
        KdAccelNode *n = AllocAligned<KdAccelNode>(nAlloc);
        if (nAllocedNodes > 0) {
            memcpy(n, nodes, nAllocedNodes * sizeof(KdAccelNode));
            FreeAligned(nodes);
        }
        nodes = n;
        nAllocedNodes = nAlloc;
    }
    ++nextFreeNode;

    // Initialize leaf node if termination criteria met
    if (nPrimitives <= maxPrims || depth == 0) {
        PBRT_KDTREE_CREATED_LEAF(nPrimitives, maxDepth-depth);
        nodes[nodeNum].initLeaf(primNums, nPrimitives, arena);
        return;
    }

    // Initialize interior node and continue recursion

    // Choose split axis position for interior node
    int bestAxis = -1, bestOffset = -1;
    float bestCost = INFINITY;
    float oldCost = isectCost * float(nPrimitives);
    float totalSA = nodeBounds.SurfaceArea();
    float invTotalSA = 1.f / totalSA;
    Vector d = nodeBounds.pMax - nodeBounds.pMin;

    // Choose which axis to split along
    uint32_t axis = nodeBounds.MaximumExtent();
    int retries = 0;
    retrySplit:

    // Initialize edges for _axis_
    for (int i = 0; i < nPrimitives; ++i) {
        int pn = primNums[i];
        const BBox &bbox = allPrimBounds[pn];
        edges[axis][2*i] =   BoundEdge(bbox.pMin[axis], pn, true);
        edges[axis][2*i+1] = BoundEdge(bbox.pMax[axis], pn, false);
    }
    sort(&edges[axis][0], &edges[axis][2*nPrimitives]);

    // Compute cost of all splits for _axis_ to find best
    int nBelow = 0, nAbove = nPrimitives;
    for (int i = 0; i < 2*nPrimitives; ++i) {
        if (edges[axis][i].type == BoundEdge::END) --nAbove;
        float edget = edges[axis][i].t;
        if (edget > nodeBounds.pMin[axis] &&
            edget < nodeBounds.pMax[axis]) {
            // Compute cost for split at _i_th edge
            uint32_t otherAxis0 = (axis + 1) % 3, otherAxis1 = (axis + 2) % 3;
            float belowSA = 2 * (d[otherAxis0] * d[otherAxis1] +
                                 (edget - nodeBounds.pMin[axis]) *
                                 (d[otherAxis0] + d[otherAxis1]));
            float aboveSA = 2 * (d[otherAxis0] * d[otherAxis1] +
                                 (nodeBounds.pMax[axis] - edget) *
                                 (d[otherAxis0] + d[otherAxis1]));
            float pBelow = belowSA * invTotalSA;
            float pAbove = aboveSA * invTotalSA;
            float eb = (nAbove == 0 || nBelow == 0) ? emptyBonus : 0.f;
            float cost = traversalCost +
                         isectCost * (1.f - eb) * (pBelow * nBelow + pAbove * nAbove);

            // Update best split if this is lowest cost so far
            if (cost < bestCost)  {
                bestCost = cost;
                bestAxis = axis;
                bestOffset = i;
            }
        }
        if (edges[axis][i].type == BoundEdge::START) ++nBelow;
    }
    Assert(nBelow == nPrimitives && nAbove == 0);

    // Create leaf if no good splits were found
    if (bestAxis == -1 && retries < 2) {
        ++retries;
        axis = (axis+1) % 3;
        goto retrySplit;
    }
    if (bestCost > oldCost) ++badRefines;
    if ((bestCost > 4.f * oldCost && nPrimitives < 16) ||
        bestAxis == -1 || badRefines == 3) {
        PBRT_KDTREE_CREATED_LEAF(nPrimitives, maxDepth-depth);
        nodes[nodeNum].initLeaf(primNums, nPrimitives, arena);
        return;
    }

    // Classify primitives with respect to split
    int n0 = 0, n1 = 0;
    for (int i = 0; i < bestOffset; ++i)
        if (edges[bestAxis][i].type == BoundEdge::START)
            prims0[n0++] = edges[bestAxis][i].primNum;
    for (int i = bestOffset+1; i < 2*nPrimitives; ++i)
        if (edges[bestAxis][i].type == BoundEdge::END)
            prims1[n1++] = edges[bestAxis][i].primNum;

    // Recursively initialize children nodes
    float tsplit = edges[bestAxis][bestOffset].t;
    PBRT_KDTREE_CREATED_INTERIOR_NODE(bestAxis, tsplit);
    BBox bounds0 = nodeBounds, bounds1 = nodeBounds;
    bounds0.pMax[bestAxis] = bounds1.pMin[bestAxis] = tsplit;
    buildTree(nodeNum+1, bounds0,
              allPrimBounds, prims0, n0, depth-1, edges,
              prims0, prims1 + nPrimitives, badRefines);
    uint32_t aboveChild = nextFreeNode;
    nodes[nodeNum].initInterior(bestAxis, aboveChild, tsplit);
    buildTree(aboveChild, bounds1, allPrimBounds, prims1, n1,
              depth-1, edges, prims0, prims1 + nPrimitives, badRefines);
}


bool KdTreeAccel::Intersect(const Ray &ray,
                            Intersection *isect) const {
    PBRT_KDTREE_INTERSECTION_TEST(const_cast<KdTreeAccel *>(this), const_cast<Ray *>(&ray));
    // Compute initial parametric range of ray inside kd-tree extent
    float tmin, tmax;
    if (!bounds.IntersectP(ray, &tmin, &tmax))
    {
        PBRT_KDTREE_RAY_MISSED_BOUNDS();
        return false;
    }

    // Prepare to traverse kd-tree for ray
    Vector invDir(1.f/ray.d.x, 1.f/ray.d.y, 1.f/ray.d.z);
#define MAX_TODO 64
    KdToDo todo[MAX_TODO];
    int todoPos = 0;

    // Traverse kd-tree nodes in order for ray
    bool hit = false;
    const KdAccelNode *node = &nodes[0];
    while (node != NULL) {
        // Bail out if we found a hit closer than the current node
        if (ray.maxt < tmin) break;
        if (!node->IsLeaf()) {
            PBRT_KDTREE_INTERSECTION_TRAVERSED_INTERIOR_NODE(const_cast<KdAccelNode *>(node));
            // Process kd-tree interior node

            // Compute parametric distance along ray to split plane
            int axis = node->SplitAxis();
            float tplane = (node->SplitPos() - ray.o[axis]) * invDir[axis];

            // Get node children pointers for ray
            const KdAccelNode *firstChild, *secondChild;
            int belowFirst = (ray.o[axis] <  node->SplitPos()) ||
                             (ray.o[axis] == node->SplitPos() && ray.d[axis] <= 0);
            if (belowFirst) {
                firstChild = node + 1;
                secondChild = &nodes[node->AboveChild()];
            }
            else {
                firstChild = &nodes[node->AboveChild()];
                secondChild = node + 1;
            }

            // Advance to next child node, possibly enqueue other child
            if (tplane > tmax || tplane <= 0)
                node = firstChild;
            else if (tplane < tmin)
                node = secondChild;
            else {
                // Enqueue _secondChild_ in todo list
                todo[todoPos].node = secondChild;
                todo[todoPos].tmin = tplane;
                todo[todoPos].tmax = tmax;
                ++todoPos;
                node = firstChild;
                tmax = tplane;
            }
        }
        else {
            PBRT_KDTREE_INTERSECTION_TRAVERSED_LEAF_NODE(const_cast<KdAccelNode *>(node), node->nPrimitives());
            // Check for intersections inside leaf node
            uint32_t nPrimitives = node->nPrimitives();
            if (nPrimitives == 1) {
                const Reference<Primitive> &prim = primitives[node->onePrimitive];
                // Check one primitive inside leaf node
                PBRT_KDTREE_INTERSECTION_PRIMITIVE_TEST(const_cast<Primitive *>(prim.GetPtr()));
                if (prim->Intersect(ray, isect))
                {
                    PBRT_KDTREE_INTERSECTION_HIT(const_cast<Primitive *>(prim.GetPtr()));
                    hit = true;
                }
            }
            else {
                uint32_t *prims = node->primitives;
                for (uint32_t i = 0; i < nPrimitives; ++i) {
                    const Reference<Primitive> &prim = primitives[prims[i]];
                    // Check one primitive inside leaf node
                    PBRT_KDTREE_INTERSECTION_PRIMITIVE_TEST(const_cast<Primitive *>(prim.GetPtr()));
                    if (prim->Intersect(ray, isect))
                    {
                        PBRT_KDTREE_INTERSECTION_HIT(const_cast<Primitive *>(prim.GetPtr()));
                        hit = true;
                    }
                }
            }

            // Grab next node to process from todo list
            if (todoPos > 0) {
                --todoPos;
                node = todo[todoPos].node;
                tmin = todo[todoPos].tmin;
                tmax = todo[todoPos].tmax;
            }
            else
                break;
        }
    }
    PBRT_KDTREE_INTERSECTION_FINISHED();
    return hit;
}


bool KdTreeAccel::IntersectP(const Ray &ray) const {
    PBRT_KDTREE_INTERSECTIONP_TEST(const_cast<KdTreeAccel *>(this), const_cast<Ray *>(&ray));
    // Compute initial parametric range of ray inside kd-tree extent
    float tmin, tmax;
    if (!bounds.IntersectP(ray, &tmin, &tmax))
    {
        PBRT_KDTREE_RAY_MISSED_BOUNDS();
        return false;
    }

    // Prepare to traverse kd-tree for ray
    Vector invDir(1.f/ray.d.x, 1.f/ray.d.y, 1.f/ray.d.z);
#define MAX_TODO 64
    KdToDo todo[MAX_TODO];
    int todoPos = 0;
    const KdAccelNode *node = &nodes[0];
    while (node != NULL) {
        if (node->IsLeaf()) {
            PBRT_KDTREE_INTERSECTIONP_TRAVERSED_LEAF_NODE(const_cast<KdAccelNode *>(node), node->nPrimitives());
            // Check for shadow ray intersections inside leaf node
            uint32_t nPrimitives = node->nPrimitives();
            if (nPrimitives == 1) {
                const Reference<Primitive> &prim = primitives[node->onePrimitive];
                PBRT_KDTREE_INTERSECTIONP_PRIMITIVE_TEST(const_cast<Primitive *>(prim.GetPtr()));
                if (prim->IntersectP(ray)) {
                    PBRT_KDTREE_INTERSECTIONP_HIT(const_cast<Primitive *>(prim.GetPtr()));
                    return true;
                }
            }
            else {
                uint32_t *prims = node->primitives;
                for (uint32_t i = 0; i < nPrimitives; ++i) {
                    const Reference<Primitive> &prim = primitives[prims[i]];
                    PBRT_KDTREE_INTERSECTIONP_PRIMITIVE_TEST(const_cast<Primitive *>(prim.GetPtr()));
                    if (prim->IntersectP(ray)) {
                        PBRT_KDTREE_INTERSECTIONP_HIT(const_cast<Primitive *>(prim.GetPtr()));
                        return true;
                    }
                }
            }

            // Grab next node to process from todo list
            if (todoPos > 0) {
                --todoPos;
                node = todo[todoPos].node;
                tmin = todo[todoPos].tmin;
                tmax = todo[todoPos].tmax;
            }
            else
                break;
        }
        else {
            PBRT_KDTREE_INTERSECTIONP_TRAVERSED_INTERIOR_NODE(const_cast<KdAccelNode *>(node));
            // Process kd-tree interior node

            // Compute parametric distance along ray to split plane
            int axis = node->SplitAxis();
            float tplane = (node->SplitPos() - ray.o[axis]) * invDir[axis];

            // Get node children pointers for ray
            const KdAccelNode *firstChild, *secondChild;
            int belowFirst = (ray.o[axis] <  node->SplitPos()) ||
                             (ray.o[axis] == node->SplitPos() && ray.d[axis] <= 0);
            if (belowFirst) {
                firstChild = node + 1;
                secondChild = &nodes[node->AboveChild()];
            }
            else {
                firstChild = &nodes[node->AboveChild()];
                secondChild = node + 1;
            }

            // Advance to next child node, possibly enqueue other child
            if (tplane > tmax || tplane <= 0)
                node = firstChild;
            else if (tplane < tmin)
                node = secondChild;
            else {
                // Enqueue _secondChild_ in todo list
                todo[todoPos].node = secondChild;
                todo[todoPos].tmin = tplane;
                todo[todoPos].tmax = tmax;
                ++todoPos;
                node = firstChild;
                tmax = tplane;
            }
        }
    }
    PBRT_KDTREE_INTERSECTIONP_MISSED();
    return false;
}


KdTreeAccel *CreateKdTreeAccelerator(const vector<Reference<Primitive> > &prims,
        const ParamSet &ps) {
    int isectCost = ps.FindOneInt("intersectcost", 80);
    int travCost = ps.FindOneInt("traversalcost", 1);
    float emptyBonus = ps.FindOneFloat("emptybonus", 0.5f);
    int maxPrims = ps.FindOneInt("maxprims", 1);
    int maxDepth = ps.FindOneInt("maxdepth", -1);
    return new KdTreeAccel(prims, isectCost, travCost,
        emptyBonus, maxPrims, maxDepth);
}


