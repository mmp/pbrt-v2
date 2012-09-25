
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


// samplers/bestcandidate.cpp*
#include "stdafx.h"
#include "samplers/bestcandidate.h"
#include "camera.h"
#include "montecarlo.h"

// BestCandidateSampler Method Definitions
#include "samplers/bestcandidate.out"
Sampler *BestCandidateSampler::GetSubSampler(int num, int count) {
    int x0, x1, y0, y1;
    ComputeSubWindow(num, count, &x0, &x1, &y0, &y1);
    if (x0 == x1 || y0 == y1) return NULL;
    return new BestCandidateSampler(x0, x1, y0, y1, samplesPerPixel,
        shutterOpen, shutterClose);
}


int BestCandidateSampler::GetMoreSamples(Sample *sample, RNG &rng) {
again:
    if (tableOffset == SAMPLE_TABLE_SIZE) {
        // Advance to next best-candidate sample table position
        tableOffset = 0;
        if (++xTile > xTileEnd) {
            xTile = xTileStart;
            if (++yTile > yTileEnd)
                return 0;
        }

        // Update sample shifts
        RNG tileRng(xTile + (yTile<<8));
        for (int i = 0; i < 3; ++i)
            sampleOffsets[i] = tileRng.RandomFloat();
    }
    // Compute raster sample from table
#define WRAP(x) ((x) > 1 ? ((x)-1) : (x))
    sample->imageX = (xTile + sampleTable[tableOffset][0]) * tableWidth;
    sample->imageY = (yTile + sampleTable[tableOffset][1]) * tableWidth;
    sample->time  = Lerp(WRAP(sampleOffsets[0] + sampleTable[tableOffset][2]),
                              shutterOpen, shutterClose);
    sample->lensU = WRAP(sampleOffsets[1] +
                         sampleTable[tableOffset][3]);
    sample->lensV = WRAP(sampleOffsets[2] +
                         sampleTable[tableOffset][4]);

    // Check sample against crop window, goto _again_ if outside
    if (sample->imageX < xPixelStart || sample->imageX >= xPixelEnd ||
        sample->imageY < yPixelStart || sample->imageY >= yPixelEnd) {
        ++tableOffset;
        goto again;
    }

    // Compute integrator samples for best-candidate sample
    for (uint32_t i = 0; i < sample->n1D.size(); ++i)
         LDShuffleScrambled1D(sample->n1D[i], 1, sample->oneD[i], rng);
    for (uint32_t i = 0; i < sample->n2D.size(); ++i)
         LDShuffleScrambled2D(sample->n2D[i], 1, sample->twoD[i], rng);
    ++tableOffset;
    return 1;
}


BestCandidateSampler *CreateBestCandidateSampler(const ParamSet &params, const Film *film,
        const Camera *camera) {
    // Initialize common sampler parameters
    int xstart, xend, ystart, yend;
    film->GetSampleExtent(&xstart, &xend, &ystart, &yend);
    int nsamp = params.FindOneInt("pixelsamples", 4);
    if (PbrtOptions.quickRender) nsamp = 1;
    return new BestCandidateSampler(xstart, xend, ystart, yend, nsamp,
         camera->shutterOpen, camera->shutterClose);
}


