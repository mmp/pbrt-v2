
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


// samplers/halton.cpp*
#include "stdafx.h"
#include "samplers/halton.h"
#include "paramset.h"
#include "camera.h"
#include "montecarlo.h"

// HaltonSampler Method Definitions
Sampler *HaltonSampler::GetSubSampler(int num, int count) {
    int x0, x1, y0, y1;
    ComputeSubWindow(num, count, &x0, &x1, &y0, &y1);
    if (x0 == x1 || y0 == y1) return NULL;
    return new HaltonSampler(x0, x1, y0, y1, samplesPerPixel, shutterOpen,
        shutterClose);
}


HaltonSampler::HaltonSampler(int xs, int xe, int ys, int ye, int ps,
        float sopen, float sclose)
    : Sampler(xs, xe, ys, ye, ps, sopen, sclose) {
    int delta = max(xPixelEnd - xPixelStart,
                    yPixelEnd - yPixelStart);
    wantedSamples = samplesPerPixel * delta * delta;
    currentSample = 0;
}


int HaltonSampler::GetMoreSamples(Sample *samples, RNG &rng) {
retry:
    if (currentSample >= wantedSamples) return 0;
    // Generate sample with Halton sequence and reject if outside image extent
    float u = (float)RadicalInverse(currentSample, 3);
    float v = (float)RadicalInverse(currentSample, 2);
    float lerpDelta = float(max(xPixelEnd - xPixelStart,
                                yPixelEnd - yPixelStart));
    samples->imageX = Lerp(u, xPixelStart, xPixelStart + lerpDelta);
    samples->imageY = Lerp(v, yPixelStart, yPixelStart + lerpDelta);
    ++currentSample;
    if (samples->imageX >= xPixelEnd || samples->imageY >= yPixelEnd)
        goto retry;

    // Generate lens, time, and integrator samples for _HaltonSampler_
    samples->lensU = (float)RadicalInverse(currentSample, 5);
    samples->lensV = (float)RadicalInverse(currentSample, 7);
    samples->time = Lerp((float)RadicalInverse(currentSample, 11),
                         shutterOpen, shutterClose);
    for (uint32_t i = 0; i < samples->n1D.size(); ++i)
        LatinHypercube(samples->oneD[i], samples->n1D[i], 1, rng);
    for (uint32_t i = 0; i < samples->n2D.size(); ++i)
        LatinHypercube(samples->twoD[i], samples->n2D[i], 2, rng);
    return 1;
}


HaltonSampler *CreateHaltonSampler(const ParamSet &params, const Film *film,
         const Camera *camera) {
    // Initialize common sampler parameters
    int xstart, xend, ystart, yend;
    film->GetSampleExtent(&xstart, &xend, &ystart, &yend);
    int nsamp = params.FindOneInt("pixelsamples", 4);
    if (PbrtOptions.quickRender) nsamp = 1;
    return new HaltonSampler(xstart, xend, ystart, yend, nsamp,
         camera->shutterOpen, camera->shutterClose);
}


