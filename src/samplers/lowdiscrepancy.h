
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

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_SAMPLERS_LOWDISCREPANCY_H
#define PBRT_SAMPLERS_LOWDISCREPANCY_H

// samplers/lowdiscrepancy.h*
#include "sampler.h"
#include "paramset.h"
#include "film.h"
#include "parallel.h"

// LDSampler Declarations
class LDSampler : public Sampler {
public:
    // LDSampler Public Methods
    LDSampler(int xstart, int xend, int ystart, int yend,
              int nsamp, float sopen, float sclose);
    ~LDSampler();
    Sampler *GetSubSampler(int num, int count);
    int RoundSize(int size) const { return RoundUpPow2(size); }
    int GetMoreSamples(Sample *sample, RNG &rng);
    int MaximumSampleCount() { return nPixelSamples; }
private:
    // LDSampler Private Data
    int xPos, yPos, nPixelSamples;
    float *sampleBuf;
};


LDSampler *CreateLowDiscrepancySampler(const ParamSet &params, const Film *film,
        const Camera *camera);

#endif // PBRT_SAMPLERS_LOWDISCREPANCY_H
