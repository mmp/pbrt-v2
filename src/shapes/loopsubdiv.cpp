
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


// shapes/loopsubdiv.cpp*
#include "stdafx.h"
#include "shapes/loopsubdiv.h"
#include "shapes/trianglemesh.h"
#include "paramset.h"
#include <set>
#include <map>
using std::set;
using std::map;

// LoopSubdiv Macros
#define NEXT(i) (((i)+1)%3)
#define PREV(i) (((i)+2)%3)

// LoopSubdiv Local Structures
struct SDFace;
struct SDFace;
struct SDVertex {
    // SDVertex Constructor
    SDVertex(Point pt = Point(0,0,0))
        : P(pt), startFace(NULL), child(NULL),
          regular(false), boundary(false) { }

    // SDVertex Methods
    int valence();
    void oneRing(Point *P);
    Point P;
    SDFace *startFace;
    SDVertex *child;
    bool regular, boundary;
};


struct SDFace {
    // SDFace Constructor
    SDFace() {
        int i;
        for (i = 0; i < 3; ++i) {
            v[i] = NULL;
            f[i] = NULL;
        }
        for (i = 0; i < 4; ++i)
            children[i] = NULL;
    }

    // SDFace Methods
    int vnum(SDVertex *vert) const {
        for (int i = 0; i < 3; ++i)
            if (v[i] == vert) return i;
        Severe("Basic logic error in SDFace::vnum()");
        return -1;
    }
    SDFace *nextFace(SDVertex *vert) {
        return f[vnum(vert)];
    }
    SDFace *prevFace(SDVertex *vert) {
        return f[PREV(vnum(vert))];
    }
    SDVertex *nextVert(SDVertex *vert) {
        return v[NEXT(vnum(vert))];
    }
    SDVertex *prevVert(SDVertex *vert) {
        return v[PREV(vnum(vert))];
    }
    SDVertex *otherVert(SDVertex *v0, SDVertex *v1) {
        for (int i = 0; i < 3; ++i)
            if (v[i] != v0 && v[i] != v1)
                return v[i];
        Severe("Basic logic error in SDVertex::otherVert()");
        return NULL;
    }
    SDVertex *v[3];
    SDFace *f[3];
    SDFace *children[4];
};


struct SDEdge {
    // SDEdge Constructor
    SDEdge(SDVertex *v0 = NULL, SDVertex *v1 = NULL) {
        v[0] = min(v0, v1);
        v[1] = max(v0, v1);
        f[0] = f[1] = NULL;
        f0edgeNum = -1;
    }

    // SDEdge Comparison Function
    bool operator<(const SDEdge &e2) const {
        if (v[0] == e2.v[0]) return v[1] < e2.v[1];
        return v[0] < e2.v[0];
    }
    SDVertex *v[2];
    SDFace *f[2];
    int f0edgeNum;
};



// LoopSubdiv Inline Functions
inline int SDVertex::valence() {
    SDFace *f = startFace;
    if (!boundary) {
        // Compute valence of interior vertex
        int nf = 1;
        while ((f = f->nextFace(this)) != startFace)
            ++nf;
        return nf;
    }
    else {
        // Compute valence of boundary vertex
        int nf = 1;
        while ((f = f->nextFace(this)) != NULL)
            ++nf;
        f = startFace;
        while ((f = f->prevFace(this)) != NULL)
            ++nf;
        return nf+1;
    }
}



// LoopSubdiv Method Definitions
LoopSubdiv::LoopSubdiv(const Transform *o2w, const Transform *w2o,
                       bool ro, int nfaces, int nvertices,
                       const int *vertexIndices, const Point *P, int nl)
    : Shape(o2w, w2o, ro) {
    nLevels = nl;
    // Allocate _LoopSubdiv_ vertices and faces
    int i;
    SDVertex *verts = new SDVertex[nvertices];
    for (i = 0; i < nvertices; ++i) {
        verts[i] = SDVertex(P[i]);
        vertices.push_back(&verts[i]);
    }
    SDFace *fs = new SDFace[nfaces];
    for (i = 0; i < nfaces; ++i)
        faces.push_back(&fs[i]);

    // Set face to vertex pointers
    const int *vp = vertexIndices;
    for (i = 0; i < nfaces; ++i) {
        SDFace *f = faces[i];
        for (int j = 0; j < 3; ++j) {
            SDVertex *v = vertices[vp[j]];
            f->v[j] = v;
            v->startFace = f;
        }
        vp += 3;
    }

    // Set neighbor pointers in _faces_
    set<SDEdge> edges;
    for (i = 0; i < nfaces; ++i) {
        SDFace *f = faces[i];
        for (int edgeNum = 0; edgeNum < 3; ++edgeNum) {
            // Update neighbor pointer for _edgeNum_
            int v0 = edgeNum, v1 = NEXT(edgeNum);
            SDEdge e(f->v[v0], f->v[v1]);
            if (edges.find(e) == edges.end()) {
                // Handle new edge
                e.f[0] = f;
                e.f0edgeNum = edgeNum;
                edges.insert(e);
            }
            else {
                // Handle previously seen edge
                e = *edges.find(e);
                e.f[0]->f[e.f0edgeNum] = f;
                f->f[edgeNum] = e.f[0];
                edges.erase(e);
            }
        }
    }

    // Finish vertex initialization
    for (i = 0; i < nvertices; ++i) {
        SDVertex *v = vertices[i];
        SDFace *f = v->startFace;
        do {
            f = f->nextFace(v);
        } while (f && f != v->startFace);
        v->boundary = (f == NULL);
        if (!v->boundary && v->valence() == 6)
            v->regular = true;
        else if (v->boundary && v->valence() == 4)
            v->regular = true;
        else
            v->regular = false;
    }
}


LoopSubdiv::~LoopSubdiv() {
    delete[] vertices[0];
    delete[] faces[0];
}


BBox LoopSubdiv::ObjectBound() const {
    BBox b;
    for (uint32_t i = 0; i < vertices.size(); i++)
        b = Union(b, vertices[i]->P);
    return b;
}


BBox LoopSubdiv::WorldBound() const {
    BBox b;
    for (uint32_t i = 0; i < vertices.size(); i++)
        b = Union(b, (*ObjectToWorld)(vertices[i]->P));
    return b;
}


bool LoopSubdiv::CanIntersect() const {
    return false;
}


void LoopSubdiv::Refine(vector<Reference<Shape> > &refined) const {
    vector<SDFace *> f = faces;
    vector<SDVertex *> v = vertices;
    MemoryArena arena;
    for (int i = 0; i < nLevels; ++i) {
        // Update _f_ and _v_ for next level of subdivision
        vector<SDFace *> newFaces;
        vector<SDVertex *> newVertices;

        // Allocate next level of children in mesh tree
        for (uint32_t j = 0; j < v.size(); ++j) {
            v[j]->child = arena.Alloc<SDVertex>();
            v[j]->child->regular = v[j]->regular;
            v[j]->child->boundary = v[j]->boundary;
            newVertices.push_back(v[j]->child);
        }
        for (uint32_t j = 0; j < f.size(); ++j)
            for (int k = 0; k < 4; ++k) {
                f[j]->children[k] = arena.Alloc<SDFace>();
                newFaces.push_back(f[j]->children[k]);
            }

        // Update vertex positions and create new edge vertices

        // Update vertex positions for even vertices
        for (uint32_t j = 0; j < v.size(); ++j) {
            if (!v[j]->boundary) {
                // Apply one-ring rule for even vertex
                if (v[j]->regular)
                    v[j]->child->P = weightOneRing(v[j], 1.f/16.f);
                else
                    v[j]->child->P = weightOneRing(v[j], beta(v[j]->valence()));
            }
            else {
                // Apply boundary rule for even vertex
                v[j]->child->P = weightBoundary(v[j], 1.f/8.f);
            }
        }

        // Compute new odd edge vertices
        map<SDEdge, SDVertex *> edgeVerts;
        for (uint32_t j = 0; j < f.size(); ++j) {
            SDFace *face = f[j];
            for (int k = 0; k < 3; ++k) {
                // Compute odd vertex on _k_th edge
                SDEdge edge(face->v[k], face->v[NEXT(k)]);
                SDVertex *vert = edgeVerts[edge];
                if (!vert) {
                    // Create and initialize new odd vertex
                    vert = arena.Alloc<SDVertex>();
                    newVertices.push_back(vert);
                    vert->regular = true;
                    vert->boundary = (face->f[k] == NULL);
                    vert->startFace = face->children[3];

                    // Apply edge rules to compute new vertex position
                    if (vert->boundary) {
                        vert->P =  0.5f * edge.v[0]->P;
                        vert->P += 0.5f * edge.v[1]->P;
                    }
                    else {
                        vert->P =  3.f/8.f * edge.v[0]->P;
                        vert->P += 3.f/8.f * edge.v[1]->P;
                        vert->P += 1.f/8.f * face->otherVert(edge.v[0], edge.v[1])->P;
                        vert->P += 1.f/8.f *
                            face->f[k]->otherVert(edge.v[0], edge.v[1])->P;
                    }
                    edgeVerts[edge] = vert;
                }
            }
        }

        // Update new mesh topology

        // Update even vertex face pointers
        for (uint32_t j = 0; j < v.size(); ++j) {
            SDVertex *vert = v[j];
            int vertNum = vert->startFace->vnum(vert);
            vert->child->startFace =
                vert->startFace->children[vertNum];
        }

        // Update face neighbor pointers
        for (uint32_t j = 0; j < f.size(); ++j) {
            SDFace *face = f[j];
            for (int k = 0; k < 3; ++k) {
                // Update children _f_ pointers for siblings
                face->children[3]->f[k] = face->children[NEXT(k)];
                face->children[k]->f[NEXT(k)] = face->children[3];

                // Update children _f_ pointers for neighbor children
                SDFace *f2 = face->f[k];
                face->children[k]->f[k] =
                    f2 ? f2->children[f2->vnum(face->v[k])] : NULL;
                f2 = face->f[PREV(k)];
                face->children[k]->f[PREV(k)] =
                    f2 ? f2->children[f2->vnum(face->v[k])] : NULL;
            }
        }

        // Update face vertex pointers
        for (uint32_t j = 0; j < f.size(); ++j) {
            SDFace *face = f[j];
            for (int k = 0; k < 3; ++k) {
                // Update child vertex pointer to new even vertex
                face->children[k]->v[k] = face->v[k]->child;

                // Update child vertex pointer to new odd vertex
                SDVertex *vert = edgeVerts[SDEdge(face->v[k], face->v[NEXT(k)])];
                face->children[k]->v[NEXT(k)] = vert;
                face->children[NEXT(k)]->v[k] = vert;
                face->children[3]->v[k] = vert;
            }
        }

        // Prepare for next level of subdivision
        f = newFaces;
        v = newVertices;
    }
    // Push vertices to limit surface
    Point *Plimit = new Point[v.size()];
    for (uint32_t i = 0; i < v.size(); ++i) {
        if (v[i]->boundary)
            Plimit[i] =  weightBoundary(v[i], 1.f/5.f);
        else
            Plimit[i] =  weightOneRing(v[i], gamma(v[i]->valence()));
    }
    for (uint32_t i = 0; i < v.size(); ++i)
        v[i]->P = Plimit[i];

    // Compute vertex tangents on limit surface
    vector<Normal> Ns;
    Ns.reserve(v.size());
    vector<Point> Pring(16, Point());
    for (uint32_t i = 0; i < v.size(); ++i) {
        SDVertex *vert = v[i];
        Vector S(0,0,0), T(0,0,0);
        int valence = vert->valence();
        if (valence > (int)Pring.size())
            Pring.resize(valence);
        vert->oneRing(&Pring[0]);
        if (!vert->boundary) {
            // Compute tangents of interior face
            for (int k = 0; k < valence; ++k) {
                S += cosf(2.f*M_PI*k/valence) * Vector(Pring[k]);
                T += sinf(2.f*M_PI*k/valence) * Vector(Pring[k]);
            }
        } else {
            // Compute tangents of boundary face
            S = Pring[valence-1] - Pring[0];
            if (valence == 2)
                T = Vector(Pring[0] + Pring[1] - 2 * vert->P);
            else if (valence == 3)
                T = Pring[1] - vert->P;
            else if (valence == 4) // regular
                T = Vector(-1*Pring[0] + 2*Pring[1] + 2*Pring[2] +
                           -1*Pring[3] + -2*vert->P);
            else {
                float theta = M_PI / float(valence-1);
                T = Vector(sinf(theta) * (Pring[0] + Pring[valence-1]));
                for (int k = 1; k < valence-1; ++k) {
                    float wt = (2*cosf(theta) - 2) * sinf((k) * theta);
                    T += Vector(wt * Pring[k]);
                }
                T = -T;
            }
        }
        Ns.push_back(Normal(Cross(S, T)));
    }

    // Create _TriangleMesh_ from subdivision mesh
    uint32_t ntris = uint32_t(f.size());
    int *verts = new int[3*ntris];
    int *vp = verts;
    uint32_t totVerts = uint32_t(v.size());
    map<SDVertex *, int> usedVerts;
    for (uint32_t i = 0; i < totVerts; ++i)
        usedVerts[v[i]] = i;
    for (uint32_t i = 0; i < ntris; ++i) {
        for (int j = 0; j < 3; ++j) {
            *vp = usedVerts[f[i]->v[j]];
            ++vp;
        }
    }
    ParamSet paramSet;
    paramSet.AddInt("indices", verts, 3*ntris);
    paramSet.AddPoint("P", Plimit, totVerts);
    paramSet.AddNormal("N", &Ns[0], int(Ns.size()));
    refined.push_back(CreateTriangleMeshShape(ObjectToWorld,
            WorldToObject, ReverseOrientation, paramSet));
    delete[] verts;
    delete[] Plimit;
}


Point LoopSubdiv::weightOneRing(SDVertex *vert, float beta) {
    // Put _vert_ one-ring in _Pring_
    int valence = vert->valence();
    Point *Pring = ALLOCA(Point, valence);
    vert->oneRing(Pring);
    Point P = (1 - valence * beta) * vert->P;
    for (int i = 0; i < valence; ++i)
        P += beta * Pring[i];
    return P;
}


void SDVertex::oneRing(Point *P) {
    if (!boundary) {
        // Get one-ring vertices for interior vertex
        SDFace *face = startFace;
        do {
            *P++ = face->nextVert(this)->P;
            face = face->nextFace(this);
        } while (face != startFace);
    }
    else {
        // Get one-ring vertices for boundary vertex
        SDFace *face = startFace, *f2;
        while ((f2 = face->nextFace(this)) != NULL)
            face = f2;
        *P++ = face->nextVert(this)->P;
        do {
            *P++ = face->prevVert(this)->P;
            face = face->prevFace(this);
        } while (face != NULL);
    }
}


Point LoopSubdiv::weightBoundary(SDVertex *vert, float beta) {
    // Put _vert_ one-ring in _Pring_
    int valence = vert->valence();
    Point *Pring = ALLOCA(Point, valence);
    vert->oneRing(Pring);
    Point P = (1-2*beta) * vert->P;
    P += beta * Pring[0];
    P += beta * Pring[valence-1];
    return P;
}


LoopSubdiv *CreateLoopSubdivShape(const Transform *o2w, const Transform *w2o,
        bool reverseOrientation, const ParamSet &params) {
    int nlevels = params.FindOneInt("nlevels", 3);
    int nps, nIndices;
    const int *vi = params.FindInt("indices", &nIndices);
    const Point *P = params.FindPoint("P", &nps);
    if (!vi || !P) return NULL;

    // don't actually use this for now...
    string scheme = params.FindOneString("scheme", "loop");

    return new LoopSubdiv(o2w, w2o, reverseOrientation, nIndices/3, nps,
        vi, P, nlevels);
}


