
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

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_SAMPLERS_BESTCANDIDATE_H
#define PBRT_SAMPLERS_BESTCANDIDATE_H

// samplers/bestcandidate.h*
#include "sampler.h"
#include "paramset.h"
#include "film.h"

// BestCandidate Sampling Constants
#define SQRT_SAMPLE_TABLE_SIZE 64
#define SAMPLE_TABLE_SIZE (SQRT_SAMPLE_TABLE_SIZE * \
                           SQRT_SAMPLE_TABLE_SIZE)

// BestCandidateSampler Declarations
class BestCandidateSampler : public Sampler {
public:
    // BestCandidateSampler Public Methods
    BestCandidateSampler(int xstart, int xend, int ystart, int yend,
                         int nPixelSamples, float sopen, float sclose)
        : Sampler(xstart, xend, ystart, yend, nPixelSamples, sopen, sclose) {
        tableWidth = (float)SQRT_SAMPLE_TABLE_SIZE /
                     (float)sqrtf(nPixelSamples);
        xTileStart = Floor2Int(xstart / tableWidth);
        xTileEnd = Floor2Int(xend / tableWidth);
        yTileStart = Floor2Int(ystart / tableWidth);
        yTileEnd = Floor2Int(yend / tableWidth);
        xTile = xTileStart;
        yTile = yTileStart;
        tableOffset = 0;
      // Update sample shifts
      RNG tileRng(xTile + (yTile<<8));
      for (int i = 0; i < 3; ++i)
          sampleOffsets[i] = tileRng.RandomFloat();
    }
    Sampler *GetSubSampler(int num, int count);
    int RoundSize(int size) const {
        return RoundUpPow2(size);
    }
    int MaximumSampleCount() { return 1; }
    int GetMoreSamples(Sample *sample, RNG &rng);
private:
    // BestCandidateSampler Private Data
    float tableWidth;
    int tableOffset;
    int xTileStart, xTileEnd, yTileStart, yTileEnd;
    int xTile, yTile;
    static const float sampleTable[SAMPLE_TABLE_SIZE][5];
    float sampleOffsets[3];
};


BestCandidateSampler *CreateBestCandidateSampler(const ParamSet &params, const Film *film,
    const Camera *camera);

#endif // PBRT_SAMPLERS_BESTCANDIDATE_H
