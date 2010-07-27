
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

#ifndef PBRT_SAMPLERS_RANDOM_H
#define PBRT_SAMPLERS_RANDOM_H

// samplers/random.h*
#include "sampler.h"
#include "paramset.h"
#include "film.h"
class RandomSampler : public Sampler {
public:
    RandomSampler(int xstart, int xend, int ystart,
        int yend, int ns, float sopen, float sclose);
    ~RandomSampler() {
        FreeAligned(imageSamples);
    }
    int MaximumSampleCount() { return 1; }
    int GetMoreSamples(Sample *sample, RNG &rng);
    int RoundSize(int sz) const { return sz; }
    Sampler *GetSubSampler(int num, int count);
private:
    // RandomSampler Private Data
    bool jitterSamples;
    int xPos, yPos, nSamples;
    float *imageSamples, *lensSamples, *timeSamples;
    int samplePos;
};


Sampler *CreateRandomSampler(const ParamSet &params, const Film *film,
    const Camera *camera);

#endif // PBRT_SAMPLERS_RANDOM_H
