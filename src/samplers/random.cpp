
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


// samplers/random.cpp*
#include "stdafx.h"
#include "samplers/random.h"
#include "montecarlo.h"
#include "camera.h"

RandomSampler::RandomSampler(int xstart, int xend,
        int ystart, int yend, int ns, float sopen, float sclose)
    : Sampler(xstart, xend, ystart, yend, ns, sopen, sclose) {
    xPos = xPixelStart;
    yPos = yPixelStart;
    nSamples = ns;
    // Get storage for a pixel's worth of stratified samples
    imageSamples = AllocAligned<float>(5 * nSamples);
    lensSamples = imageSamples + 2 * nSamples;
    timeSamples = lensSamples + 2 * nSamples;

    RNG rng(xstart + ystart * (xend-xstart));
    for (int i = 0; i < 5 * nSamples; ++i)
        imageSamples[i] = rng.RandomFloat();

    // Shift image samples to pixel coordinates
    for (int o = 0; o < 2 * nSamples; o += 2) {
        imageSamples[o]   += xPos;
        imageSamples[o+1] += yPos;
    }
    samplePos = 0;
}



Sampler *RandomSampler::GetSubSampler(int num, int count) {
    int x0, x1, y0, y1;
    ComputeSubWindow(num, count, &x0, &x1, &y0, &y1);
    if (x0 == x1 || y0 == y1) return NULL;
    return new RandomSampler(x0, x1, y0, y1, nSamples,
       shutterOpen, shutterClose);
}



int RandomSampler::GetMoreSamples(Sample *sample, RNG &rng) {
    if (samplePos == nSamples) {
        if (xPixelStart == xPixelEnd || yPixelStart == yPixelEnd)
            return 0;
        if (++xPos == xPixelEnd) {
            xPos = xPixelStart;
            ++yPos;
        }
        if (yPos == yPixelEnd)
            return 0;

        for (int i = 0; i < 5 * nSamples; ++i)
            imageSamples[i] = rng.RandomFloat();

        // Shift image samples to pixel coordinates
        for (int o = 0; o < 2 * nSamples; o += 2) {
            imageSamples[o]   += xPos;
            imageSamples[o+1] += yPos;
        }
        samplePos = 0;
    }
    // Return next \mono{RandomSampler} sample point
    sample->imageX = imageSamples[2*samplePos];
    sample->imageY = imageSamples[2*samplePos+1];
    sample->lensU = lensSamples[2*samplePos];
    sample->lensV = lensSamples[2*samplePos+1];
    sample->time = Lerp(timeSamples[samplePos], shutterOpen, shutterClose);
    // Generate stratified samples for integrators
    for (uint32_t i = 0; i < sample->n1D.size(); ++i)
        for (uint32_t j = 0; j < sample->n1D[i]; ++j)
            sample->oneD[i][j] = rng.RandomFloat();
    for (uint32_t i = 0; i < sample->n2D.size(); ++i)
        for (uint32_t j = 0; j < 2*sample->n2D[i]; ++j)
            sample->twoD[i][j] = rng.RandomFloat();
    ++samplePos;
    return 1;
}



Sampler *CreateRandomSampler(const ParamSet &params,
                       const Film *film, const Camera *camera) {
    int ns = params.FindOneInt("pixelsamples", 4);
    int xstart, xend, ystart, yend;
    film->GetSampleExtent(&xstart, &xend, &ystart, &yend);
    return new RandomSampler(xstart, xend, ystart, yend, ns,
                             camera->shutterOpen, camera->shutterClose);
}


