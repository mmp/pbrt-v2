
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


