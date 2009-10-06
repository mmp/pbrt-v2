
/*
    pbrt source code Copyright(c) 1998-2009 Matt Pharr and Greg Humphreys.

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

#ifndef PBRT_CAMERAS_REALISTIC_H
#define PBRT_CAMERAS_REALISTIC_H

// cameras/realistic.h*
#include "pbrt.h"
#include "camera.h"
#include "geometry.h"

// RealisticCamera Declarations
class LensSystem;
class RealisticCamera : public Camera {
public:
    // RealisticCamera Public Methods
    RealisticCamera(const AnimatedTransform &cam2world,
        float sopen, float sclose, Film *film, float focaldistance,
        float fstop, float filmdiagonal, float aspect, float scale,
        const string &specfile, const string &mode, bool constantWeight);
    ~RealisticCamera();
    float GenerateRay(const CameraSample &sample, Ray *) const;
private:
    // RealisticCamera Private Data
    LensSystem *lensSystem;
    float wfilm, hfilm;
    bool accurate, constantWeightRays;
};


RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film);

#endif // PBRT_CAMERAS_REALISTIC_H
