
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


// samplers/bestcandidate.cpp*
#include "samplers/bestcandidate.h"
#include "camera.h"
#include "montecarlo.h"

// BestCandidateSampler Method Definitions
BestCandidateSampler::BestCandidateSampler(int xstart, int xend,
                                           int ystart, int yend,
                                           int pixelSamples, float sopen,
                                           float sclose, u_long rngSeed)
    : Sampler(xstart, xend, ystart, yend, pixelSamples, sopen, sclose),
      rng(rngSeed) {
    tableWidth = (float)SQRT_SAMPLE_TABLE_SIZE / (float)sqrtf(pixelSamples);
    xTableCorner = float(xPixelStart) - tableWidth;
    yTableCorner = float(yPixelStart);
    tableOffset = SAMPLE_TABLE_SIZE;

    // _BestCandidateSampler_ constructor implementation
    sampleFloats = NULL;
    oneDSamples = twoDSamples = NULL;
    strat2D = NULL;
}


#include "samplers/sampledata.out"
Sampler *BestCandidateSampler::GetSubSampler(int num, int count) {
    int x0, x1, y0, y1;
    ComputeSubWindow(num, count, &x0, &x1, &y0, &y1);
    if (x0 == x1 || y0 == y1) return NULL;
    return new BestCandidateSampler(x0, x1, y0, y1, samplesPerPixel,
        shutterOpen, shutterClose, 1024*num);
}


int BestCandidateSampler::GetMoreSamples(Sample *sample) {
again:
    if (tableOffset == SAMPLE_TABLE_SIZE) {
        // Advance to next best-candidate sample table position
        tableOffset = 0;
        xTableCorner += tableWidth;
        if (xTableCorner >= xPixelEnd) {
            xTableCorner = float(xPixelStart);
            yTableCorner += tableWidth;
            if (yTableCorner >= yPixelEnd)
                return 0;
        }
        if (!oneDSamples) {
            // Initialize sample tables and precompute _strat2D_ values
            int totalFloats = 0;
            for (u_int i = 0; i < sample->n1D.size(); ++i)
                totalFloats += (sample->n1D[i] == 1) ? SAMPLE_TABLE_SIZE : 0;
            for (u_int i = 0; i < sample->n2D.size(); ++i)
                totalFloats += (sample->n2D[i] == 1) ? 2 * SAMPLE_TABLE_SIZE : 0;
            sampleFloats = new float[totalFloats];
            float *sfp = sampleFloats;
            oneDSamples = new float *[sample->n1D.size()];
            for (u_int i = 0; i < sample->n1D.size(); ++i) {
                if (sample->n1D[i] == 1) {
                    oneDSamples[i] = sfp;
                    sfp += SAMPLE_TABLE_SIZE;
                }
                else
                    oneDSamples[i] = NULL;
            }
            twoDSamples = new float *[sample->n2D.size()];
            strat2D = new int[sample->n2D.size()];
            for (u_int i = 0; i < sample->n2D.size(); ++i) {
                if (sample->n2D[i] == 1) {
                    twoDSamples[i] = sfp;
                    sfp += 2 * SAMPLE_TABLE_SIZE;
                }
                else
                    twoDSamples[i] = NULL;
                strat2D[i] = Ceil2Int(sqrtf((float)sample->n2D[i] - .5f));
            }
        }

        // Update sample shifts
        for (int i = 0; i < 3; ++i)
            sampleOffsets[i] = sample->rng->RandomFloat();

        // Generate _SAMPLE\_TABLE\_SIZE_-sized tables for single samples
        for (u_int i = 0; i < sample->n1D.size(); ++i)
            if (sample->n1D[i] == 1)
                LDShuffleScrambled1D(SAMPLE_TABLE_SIZE, 1, oneDSamples[i], rng);
        for (u_int i = 0; i < sample->n2D.size(); ++i)
            if (sample->n2D[i] == 1)
                LDShuffleScrambled2D(SAMPLE_TABLE_SIZE, 1, twoDSamples[i], rng);
    }
    // Compute raster sample from table
#define WRAP(x) ((x) > 1 ? ((x)-1) : (x))
    sample->imageX = xTableCorner + tableWidth *
                     sampleTable[tableOffset][0];
    sample->imageY = yTableCorner + tableWidth *
                     sampleTable[tableOffset][1];
    sample->time  = Lerp(WRAP(sampleOffsets[0] +
                         sampleTable[tableOffset][2]), shutterOpen, shutterClose);
    sample->lensU = WRAP(sampleOffsets[1] +
                         sampleTable[tableOffset][3]);
    sample->lensV = WRAP(sampleOffsets[2] +
                         sampleTable[tableOffset][4]);

    // Check sample against crop window, goto _again_ if outside
    if (sample->imageX <  xPixelStart ||
        sample->imageX >= xPixelEnd   ||
        sample->imageY <  yPixelStart ||
        sample->imageY >= yPixelEnd) {
        ++tableOffset;
        goto again;
    }

    // Compute integrator samples for best-candidate sample
    for (u_int i = 0; i < sample->n1D.size(); ++i) {
        if (sample->n1D[i] == 1)
            sample->oneD[i][0] = oneDSamples[i][tableOffset];
        else
            StratifiedSample1D(sample->oneD[i], sample->n1D[i], rng);
    }
    for (u_int i = 0; i < sample->n2D.size(); ++i) {
        if (sample->n2D[i] == 1) {
           sample->twoD[i][0] = twoDSamples[i][2*tableOffset];
           sample->twoD[i][1] = twoDSamples[i][2*tableOffset+1];
        }
        else
            StratifiedSample2D(sample->twoD[i], strat2D[i], strat2D[i], rng);
    }
    ++tableOffset;
    return 1;
}


BestCandidateSampler *CreateBestCandidateSampler(const ParamSet &params, const Film *film,
        const Camera *camera) {
    // Initialize common sampler parameters
    int xstart, xend, ystart, yend;
    film->GetSampleExtent(&xstart, &xend, &ystart, &yend);
    int nsamp = params.FindOneInt("pixelsamples", 4);
    if (getenv("PBRT_QUICK_RENDER")) nsamp = 1;
    return new BestCandidateSampler(xstart, xend, ystart, yend, nsamp,
         camera->shutterOpen, camera->shutterClose, 0);
}


