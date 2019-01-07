
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


// accelerators/kdtreeaccel.cpp*
#include "stdafx.h"
#include "accelerators/kdtreeaccel.h"
#include "paramset.h"

// KdTreeAccel Local Declarations
struct KdAccelNode {
    // KdAccelNode Methods
	void initLeaf(BoundEdge *edges[3], const uint32_t* edge_count, int np, MemoryArena &arena);
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

    union {
        uint32_t flags;         // Both
        uint32_t nPrims;        // Leaf
        uint32_t aboveChild;    // Interior
    };
};


struct BoundEdge {
    // BoundEdge Public Methods
    BoundEdge() { }
    BoundEdge(float tt, int pn, int in_type) {
        t = tt;
        primNum = pn;
		switch (in_type)
		{
		case 0:
			type = START;
			break;
		case 1:
			type = END;
			break;
		default:
			type = FLAT;
			break;
		}
    }
    bool operator<(const BoundEdge &e) const {
        if (t == e.t)
            return (int)type < (int)e.type;
        else return t < e.t;
    }
    float t;
    int primNum;
    enum { START, END, FLAT } type;
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
	for (uint32_t i = 0; i < primitives.size(); ++i) {
		BBox b = primitives[i]->WorldBound();
		bounds = Union(bounds, b);
	}

	// Allocate working memory for kd-tree construction
	BoundEdge *edges[3];
	for (int i = 0; i < 3; ++i)
		edges[i] = new BoundEdge[2 * primitives.size()];
	uint32_t edges_count[3] = { 0, 0, 0 };
	m_temp = new unsigned char[primitives.size()];

	// create all the split candidates
	for (uint32_t k = 0; k < 3; ++k)
	{
		unsigned int& ec = edges_count[k];
		for (unsigned i = 0; i < primitives.size(); ++i)
		{
			const Reference<Primitive> &prim = primitives[i];
			const BBox& b = prim->WorldBound();
			if (b.pMax[k] != b.pMin[k])
			{
				edges[k][ec] = BoundEdge(b.pMin[k], (int)i, 0);
				edges[k][ec + 1] = BoundEdge(b.pMax[k], (int)i, 1);
				ec += 2;
			}
			else
			{
				edges[k][ec] = BoundEdge(b.pMin[k], (int)i, 2);
				++ec;
			}
		}
	}
	for (int k = 0; k < 3; ++k)
		sort(edges[k], edges[k] + edges_count[k]);

    // Start recursive construction of kd-tree
    buildTree(0, bounds, primitives.size(), maxDepth, edges, edges_count);

	// delete temporary memory
	delete[] m_temp;

    // Free working memory for kd-tree construction
    PBRT_KDTREE_FINISHED_CONSTRUCTION(this);
}


void KdAccelNode::initLeaf(BoundEdge *edges[3], const uint32_t* edge_count, int np, MemoryArena &arena) {
    flags = 3;
    nPrims |= (np << 2);
    // Store primitive ids for leaf node
    if (np == 0)
        onePrimitive = 0;
	else if (np == 1)
	{
		for (unsigned i = 0; i < edge_count[0]; i++)
			if (edges[0][i].type == BoundEdge::START || edges[0][i].type == BoundEdge::FLAT)
			{
				onePrimitive = edges[0][i].primNum;
				break;
			}
	}else
	{
		primitives = arena.Alloc<uint32_t>(np);
		for (unsigned i = 0; i < edge_count[0]; i++)
			if (edges[0][i].type == BoundEdge::START || edges[0][i].type == BoundEdge::FLAT)
				primitives[i] = edges[0][i].primNum;
		for (int i = 0; i < 3; ++i)
			delete[] edges[i];
    }
}


KdTreeAccel::~KdTreeAccel() {
    FreeAligned(nodes);
}


void KdTreeAccel::buildTree(int nodeNum, const BBox &nodeBounds,int nPrimitives, int depth, BoundEdge *edges[3], const uint32_t* edge_count) {
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
		PBRT_KDTREE_CREATED_LEAF(nPrimitives, maxDepth - depth);
		nodes[nodeNum].initLeaf(edges, edge_count, nPrimitives, arena);
		return;
	}

	// Step 1
	// pick best sah split
	float		split_pos;
	unsigned	split_Axis;
	bool		left = false;
	float sah = pickBestSplit(edges, edge_count, nPrimitives, nodeBounds, split_Axis, split_pos, left);
	nodes[nodeNum].split = split_pos;
	nodes[nodeNum].flags = split_Axis;
	if (sah >= nPrimitives)
	{
		PBRT_KDTREE_CREATED_LEAF(nPrimitives, maxDepth - depth);
		nodes[nodeNum].initLeaf(edges, edge_count, nPrimitives, arena);
		return;
	}

	// Step 2
	// mark triangles
	BoundEdge* _edges = edges[split_Axis];
	unsigned l_num = 0, r_num = 0, b_num = nPrimitives;
	for (unsigned i = 0; i < edge_count[split_Axis]; i++)
		m_temp[_edges[i].primNum] = 2;
	for (unsigned i = 0; i < edge_count[split_Axis]; i++)
	{
		if (_edges[i].type == BoundEdge::END && _edges[i].t <= split_pos)
		{
			m_temp[_edges[i].primNum] = 0;
			b_num--;
			l_num++;
		}
		else if (_edges[i].type == BoundEdge::START && _edges[i].t >= split_pos)
		{
			m_temp[_edges[i].primNum] = 1;
			b_num--;
			r_num++;
		}
		else if (_edges[i].type == BoundEdge::FLAT)
		{
			if (_edges[i].t < split_pos || (_edges[i].t == split_pos && left))
			{
				m_temp[_edges[i].primNum] = 0;
				b_num--;
				l_num++;
			}
			else if (_edges[i].t > split_pos || (_edges[i].t == split_pos && !left))
			{
				m_temp[_edges[i].primNum] = 1;
				b_num--;
				r_num++;
			}
		}
	}

	// Step 3
	// Generate new triangle lists
	BoundEdge *l_splits[3];
	BoundEdge *r_splits[3];

	uint32_t ledges_count[3] = { 0, 0, 0 };
	uint32_t redges_count[3] = { 0, 0, 0 };
	for (int i = 0; i < 3; ++i)
	{
		l_splits[i] = new BoundEdge[2 * (l_num + b_num)];
		r_splits[i] = new BoundEdge[2 * (r_num + b_num)];
	}

	for (int k = 0; k < 3; k++)
	{
		unsigned split_count = edge_count[k];
		for (unsigned i = 0; i < split_count; i++)
		{
			const BoundEdge& old = edges[k][i];
			unsigned id = old.primNum;
			if (m_temp[id] == 0 || m_temp[id] == 2)
			{
				l_splits[k][ledges_count[k]] = old;
				ledges_count[k]++;
			}
			if (m_temp[id] == 1 || m_temp[id] == 2)
			{
				r_splits[k][redges_count[k]] = old;
				redges_count[k]++;
			}
		}
	}
	for (int i = 0; i < 3; ++i)
		delete[] edges[i];

	// Recursively initialize children nodes
	PBRT_KDTREE_CREATED_INTERIOR_NODE(split_Axis, split_pos);
	BBox bounds0 = nodeBounds, bounds1 = nodeBounds;
	bounds0.pMax[split_Axis] = bounds1.pMin[split_Axis] = split_pos;
	buildTree(nodeNum + 1, bounds0, l_num + b_num, depth - 1, l_splits, ledges_count);
	uint32_t aboveChild = nextFreeNode;
	nodes[nodeNum].initInterior(split_Axis, aboveChild, split_pos);
	buildTree(aboveChild, bounds1, r_num + b_num, depth - 1, r_splits, redges_count); 
}

// pick best sah split
float KdTreeAccel::pickBestSplit(BoundEdge* edges[3], const uint32_t* edge_count, 
		int tri_num, const BBox& box, unsigned& splitAxis, float& split_pos, bool& left)
{
	float min_sah = INFINITY;
	for (int k = 0; k < 3; k++)
	{
		unsigned n_l = 0;
		unsigned n_r = tri_num;
		unsigned split_count = edge_count[k];
		unsigned i = 0;
		while (i < split_count)
		{
			unsigned pe = 0, pf = 0, ps = 0;
			float split = edges[k][i].t;

			while (i < split_count && edges[k][i].t == split && edges[k][i].type == BoundEdge::END)
			{
				i++; pe++;
			}
			while (i < split_count && edges[k][i].t == split && edges[k][i].type == BoundEdge::FLAT)
			{
				i++; pf++;
			}
			while (i < split_count && edges[k][i].t == split && edges[k][i].type == BoundEdge::START)
			{
				i++; ps++;
			}

			n_r -= pe;
			n_r -= pf;

			// get the sah value
			bool _left = false;
			float sah = _sah(n_l, n_r, pf, k, split, box, _left);
			if (sah < min_sah)
			{
				min_sah = sah;
				splitAxis = k;
				split_pos = split;
				left = _left;
			}

			n_l += ps;
			n_l += pf;
		}
	}
	
	return min_sah;
}

// evaluate sah value of a specific split
float KdTreeAccel::_sah(unsigned l, unsigned r, unsigned f, unsigned axis, float split, const BBox& box, bool& left)
{
	float inv_sarea = 2.0f / box.SurfaceArea();

	Vector delta = box.pMax - box.pMin;
	delta[axis] = split - box.pMin[axis];
	float l_sarea = delta.x * delta.y + delta.y * delta.z + delta.z * delta.x;
	delta[axis] = box.pMax[axis] - split;
	float r_sarea = delta.x * delta.y + delta.y * delta.z + delta.z * delta.x;

	if (f != 0)
	{
		l_sarea *= inv_sarea;
		r_sarea *= inv_sarea;
		float sah0 = l_sarea * (l + f) + r_sarea * r;
		float sah1 = l_sarea * l + r_sarea * (r + f);
		left = (sah0 <= sah1);
		return min(sah0, sah1);
	}

	return (l * l_sarea + r * r_sarea) * inv_sarea;
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


