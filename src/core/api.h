
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

#ifndef PBRT_CORE_API_H
#define PBRT_CORE_API_H

// core/api.h*
#include "pbrt.h"

// API Function Declarations
void pbrtInit(const Options &opt);
void pbrtCleanup();
void pbrtIdentity();
void pbrtTranslate(float dx, float dy, float dz);
void pbrtRotate(float angle, float ax, float ay, float az);
void pbrtScale(float sx, float sy, float sz);
void pbrtLookAt(float ex, float ey, float ez,
                float lx, float ly, float lz,
                float ux, float uy, float uz);
void pbrtConcatTransform(float transform[16]);
void pbrtTransform(float transform[16]);
void pbrtCoordinateSystem(const string &);
void pbrtCoordSysTransform(const string &);
void pbrtActiveTransformAll();
void pbrtActiveTransformEndTime();
void pbrtActiveTransformStartTime();
void pbrtTransformTimes(float start, float end);
void pbrtPixelFilter(const string &name, const ParamSet &params);
void pbrtFilm(const string &type, const ParamSet &params);
void pbrtSampler(const string &name, const ParamSet &params);
void pbrtAccelerator(const string &name, const ParamSet &params);
void pbrtSurfaceIntegrator(const string &name, const ParamSet &params);
void pbrtVolumeIntegrator(const string &name, const ParamSet &params);
void pbrtRenderer(const string &name, const ParamSet &params);
void pbrtCamera(const string &, const ParamSet &cameraParams);
void pbrtWorldBegin();
void pbrtAttributeBegin();
void pbrtAttributeEnd();
void pbrtTransformBegin();
void pbrtTransformEnd();
void pbrtTexture(const string &name, const string &type,
    const string &texname, const ParamSet &params);
void pbrtMaterial(const string &name, const ParamSet &params);
void pbrtMakeNamedMaterial(const string &name, const ParamSet &params);
void pbrtNamedMaterial(const string &name);
void pbrtLightSource(const string &name, const ParamSet &params);
void pbrtAreaLightSource(const string &name, const ParamSet &params);
void pbrtShape(const string &name, const ParamSet &params);
void pbrtReverseOrientation();
void pbrtVolume(const string &name, const ParamSet &params);
void pbrtObjectBegin(const string &name);
void pbrtObjectEnd();
void pbrtObjectInstance(const string &name);
void pbrtWorldEnd();

#endif // PBRT_CORE_API_H
