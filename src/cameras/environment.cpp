
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


// cameras/environment.cpp*
#include "stdafx.h"
#include "cameras/environment.h"
#include "paramset.h"
#include "sampler.h"

// EnvironmentCamera Method Definitions
float EnvironmentCamera::GenerateRay(const CameraSample &sample,
                                     Ray *ray) const {
    float time = Lerp(sample.time, shutterOpen, shutterClose);
    // Compute environment camera ray direction
    float theta = M_PI * sample.imageY / film->yResolution;
    float phi = 2 * M_PI * sample.imageX / film->xResolution;
    Vector dir(sinf(theta) * cosf(phi), cosf(theta),
               sinf(theta) * sinf(phi));
    *ray = Ray(Point(0,0,0), dir, 0.f, INFINITY, time);
    CameraToWorld(*ray, ray);
    return 1.f;
}


EnvironmentCamera *CreateEnvironmentCamera(const ParamSet &params,
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
    (void) lensradius; // don't need this
    (void) focaldistance; // don't need this
    return new EnvironmentCamera(cam2world, shutteropen, shutterclose, film);
}


