
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


// core/material.cpp*
#include "stdafx.h"
#include "material.h"
#include "primitive.h"
#include "texture.h"
#include "spectrum.h"
#include "reflection.h"

// Material Method Definitions
Material::~Material() {
}


void Material::Bump(const Reference<Texture<float> > &d,
                    const DifferentialGeometry &dgGeom,
                    const DifferentialGeometry &dgs,
                    DifferentialGeometry *dgBump) {
    // Compute offset positions and evaluate displacement texture
    DifferentialGeometry dgEval = dgs;

    // Shift _dgEval_ _du_ in the $u$ direction
    float du = .5f * (fabsf(dgs.dudx) + fabsf(dgs.dudy));
    if (du == 0.f) du = .01f;
    dgEval.p = dgs.p + du * dgs.dpdu;
    dgEval.u = dgs.u + du;
    dgEval.nn = Normalize((Normal)Cross(dgs.dpdu, dgs.dpdv) +
                          du * dgs.dndu);
    float uDisplace = d->Evaluate(dgEval);

    // Shift _dgEval_ _dv_ in the $v$ direction
    float dv = .5f * (fabsf(dgs.dvdx) + fabsf(dgs.dvdy));
    if (dv == 0.f) dv = .01f;
    dgEval.p = dgs.p + dv * dgs.dpdv;
    dgEval.u = dgs.u;
    dgEval.v = dgs.v + dv;
    dgEval.nn = Normalize((Normal)Cross(dgs.dpdu, dgs.dpdv) +
                          dv * dgs.dndv);
    float vDisplace = d->Evaluate(dgEval);
    float displace = d->Evaluate(dgs);

    // Compute bump-mapped differential geometry
    *dgBump = dgs;
    dgBump->dpdu = dgs.dpdu + (uDisplace - displace) / du * Vector(dgs.nn) +
                   displace * Vector(dgs.dndu);
    dgBump->dpdv = dgs.dpdv + (vDisplace - displace) / dv * Vector(dgs.nn) +
                   displace * Vector(dgs.dndv);
    dgBump->nn = Normal(Normalize(Cross(dgBump->dpdu, dgBump->dpdv)));
    if (dgs.shape->ReverseOrientation ^ dgs.shape->TransformSwapsHandedness)
        dgBump->nn *= -1.f;

    // Orient shading normal to match geometric normal
    dgBump->nn = Faceforward(dgBump->nn, dgGeom.nn);
}


