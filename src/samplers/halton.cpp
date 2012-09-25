
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


