
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

#ifndef PBRT_CORE_CAMERA_H
#define PBRT_CORE_CAMERA_H

// core/camera.h*
#include "pbrt.h"
#include "geometry.h"
#include "transform.h"

// Camera Declarations
class Camera {
public:
    // Camera Interface
    Camera(const AnimatedTransform &cam2world, float sopen, float sclose,
           Film *film);
    virtual ~Camera();
    virtual float GenerateRay(const CameraSample &sample,
                              Ray *ray) const = 0;
    virtual float GenerateRayDifferential(const CameraSample &sample, RayDifferential *rd) const;

    // Camera Public Data
    AnimatedTransform CameraToWorld;
    const float shutterOpen, shutterClose;
    Film *film;
};


class ProjectiveCamera : public Camera {
public:
    // ProjectiveCamera Public Methods
    ProjectiveCamera(const AnimatedTransform &cam2world,
        const Transform &proj, const float screenWindow[4],
        float sopen, float sclose, float lensr, float focald, Film *film);
protected:
    // ProjectiveCamera Protected Data
    Transform CameraToScreen, RasterToCamera;
    Transform ScreenToRaster, RasterToScreen;
    float lensRadius, focalDistance;
};



#endif // PBRT_CORE_CAMERA_H
