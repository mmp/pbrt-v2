
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


// accelerators/bvh.cpp*
#include "accelerators/bvh.h"
#include "probes.h"
#include "paramset.h"

// BVHAccel Local Declarations
struct BVHPrimitiveInfo {
    BVHPrimitiveInfo() { }
    BVHPrimitiveInfo(int pn, const BBox &b)
        : primitiveNumber(pn), bounds(b) {
        centroid = .5f * b.pMin + .5f * b.pMax;
    }
    int primitiveNumber;
    Point centroid;
    BBox bounds;
};


struct BVHBuildNode {
    // BVHBuildNode Public Methods
    BVHBuildNode() { children[0] = children[1] = NULL; }
    void InitLeaf(u_int first, u_int n,
            const BBox &b) {
        firstPrimOffset = first;
        nPrimitives = n;
        bounds = b;
    }
    void InitInterior(u_int axis, BVHBuildNode *c0,
            BVHBuildNode *c1) {
        children[0] = c0;
        children[1] = c1;
        bounds = Union(c0->bounds, c1->bounds);
        splitAxis = axis;
        nPrimitives = 0;
    }
    BBox bounds;
    BVHBuildNode *children[2];
    u_int splitAxis;
    u_int firstPrimOffset;
    u_int nPrimitives;
};


struct CompareToMid {
    CompareToMid(int d, float m) { dim = d; mid = m; }
    int dim;
    float mid;
    bool operator()(const BVHPrimitiveInfo &a) const {
        return a.centroid[dim] < mid;
    }
};


struct ComparePoints {
    ComparePoints(int d) { dim = d; }
    int dim;
    bool operator()(const BVHPrimitiveInfo &a,
                    const BVHPrimitiveInfo &b) const {
        return a.centroid[dim] < b.centroid[dim];
    }
};


struct CompareToBucket {
    CompareToBucket(int split, int num, int d, const BBox &b)
        : centroidBounds(b)
    { splitBucket = split; nBuckets = num; dim = d; }
    bool operator()(const BVHPrimitiveInfo &p) const;

    int splitBucket;
    int nBuckets, dim;
    const BBox &centroidBounds;
};


bool CompareToBucket::operator()(const BVHPrimitiveInfo &p) const {
    int b = nBuckets * ((p.centroid[dim] - centroidBounds.pMin[dim]) /
        (centroidBounds.pMax[dim] - centroidBounds.pMin[dim]));
    if (b == nBuckets) b = nBuckets-1;
    Assert(b >= 0 && b < nBuckets);
    return b <= splitBucket;
};


struct LinearBVHNode {
    BBox bounds;
    union {
        uint32_t primitivesOffset;
        uint32_t secondChildOffset;
    };
    uint8_t nPrimitives;  // 0 -> interior node
    uint8_t axis;         // interior node: xyz
    uint8_t pad[2];       // ensure 32 byte total size
};


static inline bool IntersectP(const BBox &bounds, const Ray &ray,
                       const Vector &invDir, const u_int dirIsNeg[3]) {
    // Check for ray intersection against $x$ and $y$ slabs
    float tmin = (bounds[  dirIsNeg[0]].x - ray.o.x) * invDir.x;
    float tmax = (bounds[1-dirIsNeg[0]].x - ray.o.x) * invDir.x;
    float tymin = (bounds[  dirIsNeg[1]].y - ray.o.y) * invDir.y;
    float tymax = (bounds[1-dirIsNeg[1]].y - ray.o.y) * invDir.y;
    if ((tmin > tymax) || (tymin > tmax))
        return false;
    if (tymin > tmin) tmin = tymin;
    if (tymax < tmax) tmax = tymax;

    // Check for ray intersection against $z$ slab
    float tzmin = (bounds[  dirIsNeg[2]].z - ray.o.z) * invDir.z;
    float tzmax = (bounds[1-dirIsNeg[2]].z - ray.o.z) * invDir.z;
    if ((tmin > tzmax) || (tzmin > tmax))
        return false;
    if (tzmin > tmin)
        tmin = tzmin;
    if (tzmax < tmax)
        tmax = tzmax;
    return (tmin < ray.maxt) && (tmax > ray.mint);
}



// BVHAccel Method Definitions
BVHAccel::BVHAccel(const vector<Reference<Primitive> > &p,
                   u_int mp, const string &sm) {
    maxPrimsInNode = mp;
    for (u_int i = 0; i < p.size(); ++i)
        p[i]->FullyRefine(primitives);
    if (sm == "sah")         splitMethod = SPLIT_SAH;
    else if (sm == "middle") splitMethod = SPLIT_MIDDLE;
    else if (sm == "equal")  splitMethod = SPLIT_EQUAL_COUNTS;
    else {
        Warning("BVH split method \"%s\" unknown.  Using \"sah\".",
                sm.c_str());
        splitMethod = SPLIT_SAH;
    }
    if (primitives.size() == 0) {
        nodes = NULL;
        return;
    }
    // Build BVH from _primitives_
    PBRT_BVH_STARTED_CONSTRUCTION(this, primitives.size());

    // Initialize _buildData_ array for primitives
    vector<BVHPrimitiveInfo> buildData;
    buildData.reserve(primitives.size());
    for (u_int i = 0; i < primitives.size(); ++i) {
        BBox bbox = primitives[i]->WorldBound();
        buildData.push_back(BVHPrimitiveInfo(i, bbox));
    }

    // Recursively build BVH tree for primitives
    MemoryArena buildArena;
    u_int totalNodes = 0;
    vector<Reference<Primitive> > orderedPrims;
    orderedPrims.reserve(primitives.size());
    BVHBuildNode *root = recursiveBuild(buildArena, buildData,
                                        0, primitives.size(), &totalNodes, orderedPrims);
    primitives.swap(orderedPrims);

    // Compute representation of depth-first traversal of BVH tree
    nodes = AllocAligned<LinearBVHNode>(totalNodes);
    for (u_int i = 0; i < totalNodes; ++i)
        new (&nodes[i]) LinearBVHNode;
    u_int offset = 0;
    flattenBVHTree(root, &offset);
    Assert(offset == totalNodes);
    PBRT_BVH_FINISHED_CONSTRUCTION(this);
}


BBox BVHAccel::WorldBound() const {
    return nodes ? nodes[0].bounds : BBox();
}


BVHBuildNode *
BVHAccel::recursiveBuild(MemoryArena &buildArena,
        vector<BVHPrimitiveInfo> &buildData, u_int start, u_int end,
        u_int *totalNodes, vector<Reference<Primitive> > &orderedPrims) {
    Assert(start != end);
    (*totalNodes)++;
    BVHBuildNode *node = buildArena.Alloc<BVHBuildNode>();
    u_int nPrimitives = end - start;
    if (nPrimitives <= maxPrimsInNode) {
        // Create leaf _BVHBuildNode_
        u_int firstPrimOffset = orderedPrims.size();
        BBox bbox;
        for (u_int i = start; i < end; ++i) {
            u_int primNum = buildData[i].primitiveNumber;
            orderedPrims.push_back(primitives[primNum]);
            bbox = Union(bbox, primitives[primNum]->WorldBound());
        }
        node->InitLeaf(firstPrimOffset, nPrimitives, bbox);
    }
    else {
        // Compute bound of primitive centroids, choose split dimension _dim_
        BBox centroidBounds;
        for (u_int i = start; i < end; ++i)
            centroidBounds = Union(centroidBounds, buildData[i].centroid);
        int dim = centroidBounds.MaximumExtent();

        // Partition primitives into two sets and build children
        u_int mid = (start + end) / 2;
        if (centroidBounds.pMax[dim] == centroidBounds.pMin[dim]) {
            // Create leaf _BVHBuildNode_
            u_int firstPrimOffset = orderedPrims.size();
            BBox bbox;
            for (u_int i = start; i < end; ++i) {
                u_int primNum = buildData[i].primitiveNumber;
                orderedPrims.push_back(primitives[primNum]);
                bbox = Union(bbox, primitives[primNum]->WorldBound());
            }
            node->InitLeaf(firstPrimOffset, nPrimitives, bbox);
            return node;
        }
        else {
            // Partition primitives based on _splitMethod_
            switch (splitMethod) {
            case SPLIT_MIDDLE: {
                // Partition primitives through node's midpoint
                float pmid = .5f * (centroidBounds.pMin[dim] + centroidBounds.pMax[dim]);
                BVHPrimitiveInfo *midPtr = std::partition(&buildData[start],
                                                          &buildData[end-1]+1,
                                                          CompareToMid(dim, pmid));
                mid = midPtr - &buildData[0];
                break;
            }
            case SPLIT_EQUAL_COUNTS: {
                // Partition primitives into equally-sized subsets
                mid = (start + end) / 2;
                std::nth_element(&buildData[start], &buildData[mid],
                                 &buildData[end-1]+1, ComparePoints(dim));
                break;
            }
            case SPLIT_SAH: default: {
                // Partition primitives using approximate SAH
                if (end-start <= 4) {
                    // Partition primitives into equally-sized subsets
                    mid = (start + end) / 2;
                    std::nth_element(&buildData[start], &buildData[mid],
                                     &buildData[end-1]+1, ComparePoints(dim));
                }
                else {
                    // Allocate _BucketInfo_ for SAH partition buckets
                    const int nBuckets = 12;
                    struct BucketInfo {
                        BucketInfo() { count = 0; }
                        int count;
                        BBox bounds;
                    };
                    BucketInfo buckets[nBuckets];

                    // Initialize _BucketInfo_ for SAH partition buckets
                    for (u_int i = start; i < end; ++i) {
                        int b = nBuckets *
                            ((buildData[i].centroid[dim] - centroidBounds.pMin[dim]) /
                             (centroidBounds.pMax[dim] - centroidBounds.pMin[dim]));
                        if (b == nBuckets) b = nBuckets-1;
                        Assert(b >= 0 && b < nBuckets);
                        buckets[b].count++;
                        buckets[b].bounds = Union(buckets[b].bounds, buildData[i].bounds);
                    }

                    // Compute costs for splitting after each bucket
                    float cost[nBuckets];
                    for (int i = 0; i < nBuckets; ++i) {
                        BBox b0, b1;
                        int count0 = 0, count1 = 0;
                        for (int j = 0; j <= i; ++j) {
                            b0 = Union(b0, buckets[j].bounds);
                            count0 += buckets[j].count;
                        }
                        for (int j = i+1; j < nBuckets; ++j) {
                            b1 = Union(b1, buckets[j].bounds);
                            count1 += buckets[j].count;
                        }
                        cost[i] = count0 * b0.SurfaceArea() + count1 * b1.SurfaceArea();
                    }

                    // Find bucket to split at that minimizes SAH metric
                    float minCost = cost[0];
                    u_int minCostSplit = 0;
                    for (int i = 1; i < nBuckets; ++i) {
                        if (cost[i] < minCost) {
                            minCost = cost[i];
                            minCostSplit = i;
                        }
                    }

                    // Split primitives at selected SAH bucket
                    BVHPrimitiveInfo *pmid = std::partition(&buildData[start],
                        &buildData[end-1]+1,
                        CompareToBucket(minCostSplit, nBuckets, dim, centroidBounds));
                    mid = pmid - &buildData[0];
                }
                break;
            }
            }
        }
        node->InitInterior(dim,
                           recursiveBuild(buildArena, buildData, start, mid,
                                          totalNodes, orderedPrims),
                           recursiveBuild(buildArena, buildData, mid, end,
                                          totalNodes, orderedPrims));
    }
    return node;
}


u_int BVHAccel::flattenBVHTree(BVHBuildNode *node, u_int *offset) {
    LinearBVHNode *linearNode = &nodes[*offset];
    linearNode->bounds = node->bounds;
    u_int myOffset = (*offset)++;
    if (node->nPrimitives > 0) {
        Assert(!node->children[0] && !node->children[1]);
        linearNode->primitivesOffset = node->firstPrimOffset;
        linearNode->nPrimitives = node->nPrimitives;
    }
    else {
        // Creater interior flattened BVH node
        linearNode->axis = node->splitAxis;
        linearNode->nPrimitives = 0;
        flattenBVHTree(node->children[0], offset);
        linearNode->secondChildOffset = flattenBVHTree(node->children[1],
                                                       offset);
    }
    return myOffset;
}


BVHAccel::~BVHAccel() {
    FreeAligned(nodes);
}


bool BVHAccel::Intersect(const Ray &ray,
                         Intersection *isect) const {
    if (!nodes) return false;
    PBRT_BVH_INTERSECTION_STARTED(const_cast<BVHAccel *>(this), const_cast<Ray *>(&ray));
    bool hit = false;
    Point origin = ray(ray.mint);
    Vector invDir(1.f / ray.d.x, 1.f / ray.d.y, 1.f / ray.d.z);
    u_int dirIsNeg[3] = { ray.d.x < 0, ray.d.y < 0, ray.d.z < 0 };
    // Follow ray through BVH nodes to find primitive intersections
    uint32_t todoOffset = 0, nodeNum = 0;
    uint32_t todo[64];
    while (true) {
        const LinearBVHNode *node = &nodes[nodeNum];
        // Process ray with BVH node
        if (::IntersectP(node->bounds, ray, invDir, dirIsNeg)) {
            if (node->nPrimitives > 0) {
                // Intersect ray with primitives in leaf BVH node
                PBRT_BVH_INTERSECTION_TRAVERSED_LEAF_NODE(const_cast<LinearBVHNode *>(node));
                for (u_int i = 0; i < node->nPrimitives; ++i)
                {
                    PBRT_BVH_INTERSECTION_PRIMITIVE_TEST(const_cast<Primitive *>(primitives[node->primitivesOffset+i].GetPtr()));
                    if (primitives[node->primitivesOffset+i]->Intersect(ray, isect))
                    {
                        PBRT_BVH_INTERSECTION_PRIMITIVE_HIT(const_cast<Primitive *>(primitives[node->primitivesOffset+i].GetPtr()));
                        hit = true;
                    }
                    else {
                        PBRT_BVH_INTERSECTION_PRIMITIVE_MISSED(const_cast<Primitive *>(primitives[node->primitivesOffset+i].GetPtr()));
                   }
                }
                if (todoOffset == 0) break;
                nodeNum = todo[--todoOffset];
            }
            else {
                // Put far BVH node on _todo_ stack, advance to near node
                PBRT_BVH_INTERSECTION_TRAVERSED_INTERIOR_NODE(const_cast<LinearBVHNode *>(node));
                if (dirIsNeg[node->axis]) {
                   todo[todoOffset++] = nodeNum + 1;
                   nodeNum = node->secondChildOffset;
                }
                else {
                   todo[todoOffset++] = node->secondChildOffset;
                   nodeNum = nodeNum + 1;
                }
            }
        }
        else {
            if (todoOffset == 0) break;
            nodeNum = todo[--todoOffset];
        }
    }
    PBRT_BVH_INTERSECTION_FINISHED();
    return hit;
}


bool BVHAccel::IntersectP(const Ray &ray) const {
    if (!nodes) return false;
    PBRT_BVH_INTERSECTIONP_STARTED(const_cast<BVHAccel *>(this), const_cast<Ray *>(&ray));
    Vector invDir(1.f / ray.d.x, 1.f / ray.d.y, 1.f / ray.d.z);
    u_int dirIsNeg[3] = { ray.d.x < 0, ray.d.y < 0, ray.d.z < 0 };
    uint32_t todo[64];
    uint32_t todoOffset = 0, nodeNum = 0;
    while (true) {
        const LinearBVHNode *node = &nodes[nodeNum];
        if (::IntersectP(node->bounds, ray, invDir, dirIsNeg)) {
            // Process BVH node _node_ for traversal
            if (node->nPrimitives > 0) {
                PBRT_BVH_INTERSECTIONP_TRAVERSED_LEAF_NODE(const_cast<LinearBVHNode *>(node));
                  for (u_int i = 0; i < node->nPrimitives; ++i) {
                    PBRT_BVH_INTERSECTIONP_PRIMITIVE_TEST(const_cast<Primitive *>(primitives[node->primitivesOffset + i].GetPtr()));
                    if (primitives[node->primitivesOffset+i]->IntersectP(ray)) {
                        PBRT_BVH_INTERSECTIONP_PRIMITIVE_HIT(const_cast<Primitive *>(primitives[node->primitivesOffset+i].GetPtr()));
                        return true;
                    }
                else {
                        PBRT_BVH_INTERSECTIONP_PRIMITIVE_MISSED(const_cast<Primitive *>(primitives[node->primitivesOffset + i].GetPtr()));
                    }
                }
                if (todoOffset == 0) break;
                nodeNum = todo[--todoOffset];
            }
            else {
                PBRT_BVH_INTERSECTIONP_TRAVERSED_INTERIOR_NODE(const_cast<LinearBVHNode *>(node));
                if (dirIsNeg[node->axis]) {
                   /// second child first
                   todo[todoOffset++] = nodeNum + 1;
                   nodeNum = node->secondChildOffset;
                }
                else {
                   todo[todoOffset++] = node->secondChildOffset;
                   nodeNum = nodeNum + 1;
                }
            }
        }
        else {
            if (todoOffset == 0) break;
            nodeNum = todo[--todoOffset];
        }
    }
    PBRT_BVH_INTERSECTIONP_FINISHED();
    return false;
}


BVHAccel *CreateBVHAccelerator(const vector<Reference<Primitive> > &prims,
        const ParamSet &ps) {
    string splitMethod = ps.FindOneString("splitmethod", "sah");
    u_int maxPrimsInNode = ps.FindOneInt("maxnodeprims", 1);
    return new BVHAccel(prims, maxPrimsInNode, splitMethod);
}


