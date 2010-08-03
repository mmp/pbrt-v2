
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


// cameras/orthographic.cpp*
#include "stdafx.h"
#include "cameras/orthographic.h"
#include "paramset.h"
#include "sampler.h"
#include "montecarlo.h"

// OrthographicCamera Definitions
OrthoCamera::OrthoCamera(const AnimatedTransform &cam2world,
        const float screenWindow[4], float sopen, float sclose,
        float lensr, float focald, Film *f)
    : ProjectiveCamera(cam2world, Orthographic(0., 1.), screenWindow,
                       sopen, sclose, lensr, focald, f) {
    // Compute differential changes in origin for ortho camera rays
    dxCamera = RasterToCamera(Vector(1, 0, 0));
    dyCamera = RasterToCamera(Vector(0, 1, 0));
}


float OrthoCamera::GenerateRay(const CameraSample &sample, Ray *ray) const {
    // Generate raster and camera samples
    Point Pras(sample.imageX, sample.imageY, 0);
    Point Pcamera;
    RasterToCamera(Pras, &Pcamera);
    *ray = Ray(Pcamera, Vector(0,0,1), 0.f, INFINITY);
    // Modify ray for depth of field
    if (lensRadius > 0.) {
        // Sample point on lens
        float lensU, lensV;
        ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);
        lensU *= lensRadius;
        lensV *= lensRadius;

        // Compute point on plane of focus
        float ft = focalDistance / ray->d.z;
        Point Pfocus = (*ray)(ft);

        // Update ray for effect of lens
        ray->o = Point(lensU, lensV, 0.f);
        ray->d = Normalize(Pfocus - ray->o);
    }
    ray->time = Lerp(sample.time, shutterOpen, shutterClose);
    CameraToWorld(*ray, ray);
    return 1.f;
}


float OrthoCamera::GenerateRayDifferential(const CameraSample &sample,
        RayDifferential *ray) const {
    // Compute main orthographic viewing ray

    // Generate raster and camera samples
    Point Pras(sample.imageX, sample.imageY, 0);
    Point Pcamera;
    RasterToCamera(Pras, &Pcamera);
    *ray = RayDifferential(Pcamera, Vector(0,0,1), 0., INFINITY);

    // Modify ray for depth of field
    if (lensRadius > 0.) {
        // Sample point on lens
        float lensU, lensV;
        ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);
        lensU *= lensRadius;
        lensV *= lensRadius;

        // Compute point on plane of focus
        float ft = focalDistance / ray->d.z;
        Point Pfocus = (*ray)(ft);

        // Update ray for effect of lens
        ray->o = Point(lensU, lensV, 0.f);
        ray->d = Normalize(Pfocus - ray->o);
    }
    ray->time = Lerp(sample.time, shutterOpen, shutterClose);
    ray->rxOrigin = ray->o + dxCamera;
    ray->ryOrigin = ray->o + dyCamera;
    ray->rxDirection = ray->ryDirection = ray->d;
    ray->hasDifferentials = true;
    CameraToWorld(*ray, ray);
    return 1.f;
}


OrthoCamera *CreateOrthographicCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film) {
    // Extract common camera parameters from _ParamSet_
    float shutteropen = params.FindOneFloat("shutteropen", 0.f);
    float shutterclose = params.FindOneFloat("shutterclose", 1.f);
    if (shutterclose < shutteropen) {
        Warning("Shutter close time [%f] < shutter open [%f].  Swapping them.",
                shutterclose, shutteropen);
        swap(shutterclose, shutteropen);
    }
    float lensradius = params.FindOneFloat("lensradius", 0.f);
    float focaldistance = params.FindOneFloat("focaldistance", 1e30f);
    float frame = params.FindOneFloat("frameaspectratio",
        float(film->xResolution)/float(film->yResolution));
    float screen[4];
    if (frame > 1.f) {
        screen[0] = -frame;
        screen[1] =  frame;
        screen[2] = -1.f;
        screen[3] =  1.f;
    }
    else {
        screen[0] = -1.f;
        screen[1] =  1.f;
        screen[2] = -1.f / frame;
        screen[3] =  1.f / frame;
    }
    int swi;
    const float *sw = params.FindFloat("screenwindow", &swi);
    if (sw && swi == 4)
        memcpy(screen, sw, 4*sizeof(float));
    return new OrthoCamera(cam2world, screen, shutteropen, shutterclose,
        lensradius, focaldistance, film);
}


