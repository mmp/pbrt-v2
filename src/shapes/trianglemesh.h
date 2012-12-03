
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

#ifndef PBRT_SHAPES_TRIANGLEMESH_H
#define PBRT_SHAPES_TRIANGLEMESH_H

// shapes/trianglemesh.h*
#include "shape.h"
#include <map>
using std::map;

// TriangleMesh Declarations
class TriangleMesh : public Shape {
public:
    // TriangleMesh Public Methods
    TriangleMesh(const Transform *o2w, const Transform *w2o, bool ro,
                 int ntris, int nverts, const int *vptr,
                 const Point *P, const Normal *N, const Vector *S,
                 const float *uv, const Reference<Texture<float> > &atex);
    ~TriangleMesh();
    BBox ObjectBound() const;
    BBox WorldBound() const;
    bool CanIntersect() const { return false; }
    void Refine(vector<Reference<Shape> > &refined) const;
    friend class Triangle;
    template <typename T> friend class VertexTexture;
protected:
    // TriangleMesh Protected Data
    int ntris, nverts;
    int *vertexIndex;
    Point *p;
    Normal *n;
    Vector *s;
    float *uvs;
    Reference<Texture<float> > alphaTexture;
};


class Triangle : public Shape {
public:
    // Triangle Public Methods
    Triangle(const Transform *o2w, const Transform *w2o, bool ro,
             TriangleMesh *m, int n)
        : Shape(o2w, w2o, ro) {
        mesh = m;
        v = &mesh->vertexIndex[3*n];
        PBRT_CREATED_TRIANGLE(this);
    }
    BBox ObjectBound() const;
    BBox WorldBound() const;
    bool Intersect(const Ray &ray, float *tHit, float *rayEpsilon,
                   DifferentialGeometry *dg) const;
    bool IntersectP(const Ray &ray) const;
    void GetUVs(float uv[3][2]) const {
        if (mesh->uvs) {
            uv[0][0] = mesh->uvs[2*v[0]];
            uv[0][1] = mesh->uvs[2*v[0]+1];
            uv[1][0] = mesh->uvs[2*v[1]];
            uv[1][1] = mesh->uvs[2*v[1]+1];
            uv[2][0] = mesh->uvs[2*v[2]];
            uv[2][1] = mesh->uvs[2*v[2]+1];
        }
        else {
            uv[0][0] = 0.; uv[0][1] = 0.;
            uv[1][0] = 1.; uv[1][1] = 0.;
            uv[2][0] = 1.; uv[2][1] = 1.;
        }
    }
    float Area() const;
    virtual void GetShadingGeometry(const Transform &obj2world,
            const DifferentialGeometry &dg,
            DifferentialGeometry *dgShading) const;
    Point Sample(float u1, float u2, Normal *Ns) const;
private:
    // Triangle Private Data
    Reference<TriangleMesh> mesh;
    int *v;
};


TriangleMesh *CreateTriangleMeshShape(const Transform *o2w, const Transform *w2o,
    bool reverseOrientation, const ParamSet &params,
    map<string, Reference<Texture<float> > > *floatTextures = NULL);

#endif // PBRT_SHAPES_TRIANGLEMESH_H
