
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


