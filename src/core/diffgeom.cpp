
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


// core/diffgeom.cpp*
#include "stdafx.h"
#include "diffgeom.h"
#include "transform.h"
#include "shape.h"

// DifferentialGeometry Method Definitions
DifferentialGeometry::DifferentialGeometry(const Point &P,
        const Vector &DPDU, const Vector &DPDV,
        const Normal &DNDU, const Normal &DNDV,
        float uu, float vv, const Shape *sh)
    : p(P), dpdu(DPDU), dpdv(DPDV), dndu(DNDU), dndv(DNDV) {
    // Initialize _DifferentialGeometry_ from parameters
    nn = Normal(Normalize(Cross(dpdu, dpdv)));
    u = uu;
    v = vv;
    shape = sh;
    dudx = dvdx = dudy = dvdy = 0;

    // Adjust normal based on orientation and handedness
    if (shape && (shape->ReverseOrientation ^ shape->TransformSwapsHandedness))
        nn *= -1.f;
}


void DifferentialGeometry::ComputeDifferentials(
        const RayDifferential &ray) const {
    if (ray.hasDifferentials) {
        // Estimate screen space change in $\pt{}$ and $(u,v)$

        // Compute auxiliary intersection points with plane
        float d = -Dot(nn, Vector(p.x, p.y, p.z));
        Vector rxv(ray.rxOrigin.x, ray.rxOrigin.y, ray.rxOrigin.z);
        float tx = -(Dot(nn, rxv) + d) / Dot(nn, ray.rxDirection);
        if (isnan(tx)) goto fail;
        Point px = ray.rxOrigin + tx * ray.rxDirection;
        Vector ryv(ray.ryOrigin.x, ray.ryOrigin.y, ray.ryOrigin.z);
        float ty = -(Dot(nn, ryv) + d) / Dot(nn, ray.ryDirection);
        if (isnan(ty)) goto fail;
        Point py = ray.ryOrigin + ty * ray.ryDirection;
        dpdx = px - p;
        dpdy = py - p;

        // Compute $(u,v)$ offsets at auxiliary points

        // Initialize _A_, _Bx_, and _By_ matrices for offset computation
        float A[2][2], Bx[2], By[2];
        int axes[2];
        if (fabsf(nn.x) > fabsf(nn.y) && fabsf(nn.x) > fabsf(nn.z)) {
            axes[0] = 1; axes[1] = 2;
        }
        else if (fabsf(nn.y) > fabsf(nn.z)) {
            axes[0] = 0; axes[1] = 2;
        }
        else {
            axes[0] = 0; axes[1] = 1;
        }

        // Initialize matrices for chosen projection plane
        A[0][0] = dpdu[axes[0]];
        A[0][1] = dpdv[axes[0]];
        A[1][0] = dpdu[axes[1]];
        A[1][1] = dpdv[axes[1]];
        Bx[0] = px[axes[0]] - p[axes[0]];
        Bx[1] = px[axes[1]] - p[axes[1]];
        By[0] = py[axes[0]] - p[axes[0]];
        By[1] = py[axes[1]] - p[axes[1]];
        if (!SolveLinearSystem2x2(A, Bx, &dudx, &dvdx)) {
            dudx = 0.; dvdx = 0.;
        }
        if (!SolveLinearSystem2x2(A, By, &dudy, &dvdy)) {
            dudy = 0.; dvdy = 0.;
        }
    }
    else {
fail:
        dudx = dvdx = 0.;
        dudy = dvdy = 0.;
        dpdx = dpdy = Vector(0,0,0);
    }
}


