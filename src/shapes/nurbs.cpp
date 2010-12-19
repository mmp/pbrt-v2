
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


// shapes/nurbs.cpp*
#include "stdafx.h"
#include "shapes/nurbs.h"
#include "shapes/trianglemesh.h"
#include "paramset.h"
#include "texture.h"

// NURBS Evaluation Functions
static int KnotOffset(const float *knot, int order, int np, float t) {
    int firstKnot = order - 1;

    int knotOffset = firstKnot;
    while (t > knot[knotOffset+1])
        ++knotOffset;
    Assert(knotOffset < np); // np == lastKnot
    Assert(t >= knot[knotOffset] && t <= knot[knotOffset + 1]);
    return knotOffset;
}



// doesn't handle flat out discontinuities in the curve...

struct Homogeneous3 {
Homogeneous3() { x = y = z = w = 0.; }
Homogeneous3(float xx, float yy, float zz, float ww) {
    x = xx; y = yy; z = zz; w = ww;
}


float x, y, z, w;
};


static Homogeneous3
NURBSEvaluate(int order, const float *knot, const Homogeneous3 *cp, int np,
          int cpStride, float t, Vector *deriv = NULL) {
//    int nKnots = np + order;
    float alpha;

    int knotOffset = KnotOffset(knot, order, np, t);
    knot += knotOffset;

    int cpOffset = knotOffset - order + 1;
    Assert(cpOffset >= 0 && cpOffset < np);

    Homogeneous3 *cpWork = ALLOCA(Homogeneous3, order);
    for (int i = 0; i < order; ++i)
    cpWork[i] = cp[(cpOffset+i) * cpStride];

    for (int i = 0; i < order - 2; ++i)
    for (int j = 0; j < order - 1 - i; ++j) {
        alpha = (knot[1 + j] - t) /
        (knot[1 + j] - knot[j + 2 - order + i]);
        Assert(alpha >= 0. && alpha <= 1.);

        cpWork[j].x = cpWork[j].x * alpha + cpWork[j+1].x * (1 - alpha);
        cpWork[j].y = cpWork[j].y * alpha + cpWork[j+1].y * (1 - alpha);
        cpWork[j].z = cpWork[j].z * alpha + cpWork[j+1].z * (1 - alpha);
        cpWork[j].w = cpWork[j].w * alpha + cpWork[j+1].w * (1 - alpha);
    }

    alpha = (knot[1] - t) / (knot[1] - knot[0]);
    Assert(alpha >= 0. && alpha <= 1.);

    Homogeneous3 val(cpWork[0].x * alpha + cpWork[1].x * (1 - alpha),
             cpWork[0].y * alpha + cpWork[1].y * (1 - alpha),
             cpWork[0].z * alpha + cpWork[1].z * (1 - alpha),
             cpWork[0].w * alpha + cpWork[1].w * (1 - alpha));

    if (deriv) {
        float factor = (order - 1) / (knot[1] - knot[0]);
        Homogeneous3 delta((cpWork[1].x - cpWork[0].x) * factor,
               (cpWork[1].y - cpWork[0].y) * factor,
               (cpWork[1].z - cpWork[0].z) * factor,
               (cpWork[1].w - cpWork[0].w) * factor);

        deriv->x = delta.x / val.w - (val.x * delta.w / (val.w * val.w));
        deriv->y = delta.y / val.w - (val.y * delta.w / (val.w * val.w));
        deriv->z = delta.z / val.w - (val.z * delta.w / (val.w * val.w));
    }

    return val;
}



static Point
NURBSEvaluateSurface(int uOrder, const float *uKnot, int ucp, float u,
             int vOrder, const float *vKnot, int vcp, float v,
             const Homogeneous3 *cp, Vector *dPdu, Vector *dPdv) {
    Homogeneous3 *iso = ALLOCA(Homogeneous3, max(uOrder, vOrder));

    int uOffset = KnotOffset(uKnot, uOrder, ucp, u);
    int uFirstCp = uOffset - uOrder + 1;
    Assert(uFirstCp >= 0 && uFirstCp + uOrder - 1 < ucp);

    for (int i = 0; i < uOrder; ++i)
        iso[i] = NURBSEvaluate(vOrder, vKnot, &cp[uFirstCp + i], vcp,
                   ucp, v);

    int vOffset = KnotOffset(vKnot, vOrder, vcp, v);
    int vFirstCp = vOffset - vOrder + 1;
    Assert(vFirstCp >= 0 && vFirstCp + vOrder - 1 < vcp);

    Homogeneous3 P = NURBSEvaluate(uOrder, uKnot, iso - uFirstCp, ucp,
                   1, u, dPdu);

    if (dPdv) {
        for (int i = 0; i < vOrder; ++i)
            iso[i] = NURBSEvaluate(uOrder, uKnot, &cp[(vFirstCp+i)*ucp], ucp,
                   1, u);
        (void)NURBSEvaluate(vOrder, vKnot, iso - vFirstCp, vcp, 1, v, dPdv);
    }
    return Point(P.x/P.w, P.y/P.w, P.z/P.w);;
}




// NURBS Method Definitions
NURBS::NURBS(const Transform *o2w, const Transform *w2o,
        bool ro, int numu, int uo, const float *uk,
        float u0, float u1, int numv, int vo, const float *vk,
        float v0, float v1, const float *p, bool homogeneous)
    : Shape(o2w, w2o, ro) {
    nu = numu;    uorder = uo;
    umin = u0;    umax = u1;
    nv = numv;    vorder = vo;
    vmin = v0;    vmax = v1;
    isHomogeneous = homogeneous;
    if (isHomogeneous) {
        P = new float[4*nu*nv];
        memcpy(P, p, 4*nu*nv*sizeof(float));
    } else {
        P = new float[3*nu*nv];
        memcpy(P, p, 3*nu*nv*sizeof(float));
    }
    uknot = new float[nu + uorder];
    memcpy(uknot, uk, (nu + uorder) * sizeof(float));
    vknot = new float[nv + vorder];
    memcpy(vknot, vk, (nv + vorder) * sizeof(float));
}


NURBS::~NURBS() {
    delete[] P;
    delete[] uknot;
    delete[] vknot;
}


BBox NURBS::ObjectBound() const {
    if (!isHomogeneous) {
        // Compute object-space bound of non-homogeneous NURBS
        float *pp = P;
        BBox bound;
        for (int i = 0; i < nu*nv; ++i, pp += 3)
            bound = Union(bound, Point(pp[0], pp[1], pp[2]));
        return bound;
    } else {
        // Compute object-space bound of homogeneous NURBS
        float *pp = P;
        BBox bound;
        for (int i = 0; i < nu*nv; ++i, pp += 4)
            bound = Union(bound, Point(pp[0] / pp[3], pp[1] / pp[3], pp[2] / pp[3]));
        return bound;
    }
}


BBox NURBS::WorldBound() const {
    if (!isHomogeneous) {
        // Compute world-space bound of non-homogeneous NURBS
        float *pp = P;
        BBox bound;
        for (int i = 0; i < nu*nv; ++i, pp += 3) {
            Point pt = (*ObjectToWorld)(Point(pp[0], pp[1], pp[2]));
            bound = Union(bound, pt);
        }
        return bound;
    } else {
        // Compute world-space bound of homogeneous NURBS
        float *pp = P;
        BBox bound;
        for (int i = 0; i < nu*nv; ++i, pp += 4) {
            Point pt = (*ObjectToWorld)(Point(pp[0]/pp[3],
                pp[1]/pp[3], pp[2]/pp[3]));
            bound = Union(bound, pt);
        }
        return bound;
    }
}



void NURBS::Refine(vector<Reference<Shape> > &refined) const {
    // Compute NURBS dicing rates
    int diceu = 30, dicev = 30;
    float *ueval = new float[diceu];
    float *veval = new float[dicev];
    Point *evalPs = new Point[diceu*dicev];
    Normal *evalNs = new Normal[diceu*dicev];
    int i;
    for (i = 0; i < diceu; ++i)
        ueval[i] = Lerp((float)i / (float)(diceu-1), umin, umax);
    for (i = 0; i < dicev; ++i)
        veval[i] = Lerp((float)i / (float)(dicev-1), vmin, vmax);
    // Evaluate NURBS over grid of points
    memset(evalPs, 0, diceu*dicev*sizeof(Point));
    memset(evalNs, 0, diceu*dicev*sizeof(Point));
    float *uvs = new float[2*diceu*dicev];
    // Turn NURBS into triangles
    Homogeneous3 *Pw = (Homogeneous3 *)P;
    if (!isHomogeneous) {
        Pw = new Homogeneous3[nu * nv];
        for (int i = 0; i < nu*nv; ++i) {
            Pw[i].x = P[3*i];
            Pw[i].y = P[3*i+1];
            Pw[i].z = P[3*i+2];
            Pw[i].w = 1.;
        }
    }
    for (int v = 0; v < dicev; ++v) {
        for (int u = 0; u < diceu; ++u) {
            uvs[2*(v*diceu+u)]   = ueval[u];
            uvs[2*(v*diceu+u)+1] = veval[v];

            Vector dPdu, dPdv;
            Point pt = NURBSEvaluateSurface(uorder, uknot, nu, ueval[u],
                vorder, vknot, nv, veval[v], Pw, &dPdu, &dPdv);
            evalPs[v*diceu + u].x = pt.x;
            evalPs[v*diceu + u].y = pt.y;
            evalPs[v*diceu + u].z = pt.z;
            evalNs[v*diceu + u] = Normal(Normalize(Cross(dPdu, dPdv)));
        }
    }
    // Generate points-polygons mesh
    int nTris = 2*(diceu-1)*(dicev-1);
    int *vertices = new int[3 * nTris];
    int *vertp = vertices;
    // Compute the vertex offset numbers for the triangles
    for (int v = 0; v < dicev-1; ++v) {
        for (int u = 0; u < diceu-1; ++u) {
    #define VN(u,v) ((v)*diceu+(u))
            *vertp++ = VN(u,   v);
            *vertp++ = VN(u+1, v);
            *vertp++ = VN(u+1, v+1);

            *vertp++ = VN(u,   v);
            *vertp++ = VN(u+1, v+1);
            *vertp++ = VN(u,   v+1);
    #undef VN
        }
    }
    int nVerts = diceu*dicev;
    ParamSet paramSet;
    paramSet.AddInt("indices", vertices, 3*nTris);
    paramSet.AddPoint("P", evalPs, nVerts);
    paramSet.AddFloat("uv", uvs, 2 * nVerts);
    paramSet.AddNormal("N", evalNs, nVerts);
    refined.push_back(CreateTriangleMeshShape(ObjectToWorld, WorldToObject,
            ReverseOrientation, paramSet));
    // Cleanup from NURBS refinement
    if (Pw != (Homogeneous3 *)P) delete[] Pw;
    delete[] uvs;
    delete[] ueval;
    delete[] veval;
    delete[] evalPs;
    delete[] evalNs;
    delete[] vertices;
}



NURBS *CreateNURBSShape(const Transform *o2w, const Transform *w2o,
        bool ReverseOrientation, const ParamSet &params) {
    int nu = params.FindOneInt("nu", -1);
    int uorder = params.FindOneInt("uorder", -1);
    int nuknots, nvknots;
    const float *uknots = params.FindFloat("uknots", &nuknots);
    Assert(nu != -1 && uorder != -1 && uknots != NULL);
    Assert(nuknots == nu + uorder);
    float u0 = params.FindOneFloat("u0", uknots[uorder-1]);
    float u1 = params.FindOneFloat("u1", uknots[nu]);

    int nv = params.FindOneInt("nv", -1);
    int vorder = params.FindOneInt("vorder", -1);
    const float *vknots = params.FindFloat("vknots", &nvknots);
    Assert(nv != -1 && vorder != -1 && vknots != NULL);
    Assert(nvknots == nv + vorder);
    float v0 = params.FindOneFloat("v0", vknots[vorder-1]);
    float v1 = params.FindOneFloat("v1", vknots[nv]);

    bool isHomogeneous = false;
    int npts;
    const float *P = (const float *)params.FindPoint("P", &npts);
    if (!P) {
        P = params.FindFloat("Pw", &npts);
        if (!P) {
            Error("Must provide control points via \"P\" or \"Pw\" parameter to "
                  "NURBS shape.");
            return NULL;
        }
        if ((npts % 4) != 0) {
            Error("Number of \"Pw\" control points provided to NURBS shape must be "
                  "multiple of four");
            return NULL;
        }

        npts /= 4;
        isHomogeneous = true;
    }
    if (npts != nu*nv) {
        Error("NURBS shape was expecting %dx%d=%d control points, was given %d",
              nu, nv, nu*nv, npts);
        return NULL;
    }

    return new NURBS(o2w, w2o, ReverseOrientation, nu, uorder, uknots, u0, u1,
                         nv, vorder, vknots, v0, v1, (float *)P,
                         isHomogeneous);
}


