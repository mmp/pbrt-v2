
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


// samplers/lowdiscrepancy.cpp*
#include "stdafx.h"
#include "samplers/lowdiscrepancy.h"
#include "camera.h"
#include "montecarlo.h"

// LDSampler Method Definitions
LDSampler::LDSampler(int xstart, int xend, int ystart, int yend, int ps,
                     float sopen, float sclose)
    : Sampler(xstart, xend, ystart, yend, RoundUpPow2(ps), sopen, sclose) {
    xPos = xPixelStart;
    yPos = yPixelStart;
    if (!IsPowerOf2(ps)) {
        Warning("Pixel samples being rounded up to power of 2");
        nPixelSamples = RoundUpPow2(ps);
    } else
        nPixelSamples = ps;
    sampleBuf = NULL;
}


LDSampler::~LDSampler() {
    delete[] sampleBuf;
}


Sampler *LDSampler::GetSubSampler(int num, int count) {
    int x0, x1, y0, y1;
    ComputeSubWindow(num, count, &x0, &x1, &y0, &y1);
    if (x0 == x1 || y0 == y1) return NULL;
    return new LDSampler(x0, x1, y0, y1, nPixelSamples, shutterOpen, shutterClose);
}


int LDSampler::GetMoreSamples(Sample *samples, RNG &rng) {
    if (yPos == yPixelEnd) return 0;
    if (sampleBuf == NULL)
        sampleBuf = new float[LDPixelSampleFloatsNeeded(samples,
                                                        nPixelSamples)];
    LDPixelSample(xPos, yPos, shutterOpen, shutterClose,
                  nPixelSamples, samples, sampleBuf, rng);
    if (++xPos == xPixelEnd) {
        xPos = xPixelStart;
        ++yPos;
    }
    return nPixelSamples;
}


LDSampler *CreateLowDiscrepancySampler(const ParamSet &params, const Film *film,
        const Camera *camera) {
    // Initialize common sampler parameters
    int xstart, xend, ystart, yend;
    film->GetSampleExtent(&xstart, &xend, &ystart, &yend);
    int nsamp = params.FindOneInt("pixelsamples", 4);
    if (PbrtOptions.quickRender) nsamp = 1;
    return new LDSampler(xstart, xend, ystart, yend, nsamp,
        camera->shutterOpen, camera->shutterClose);
}


