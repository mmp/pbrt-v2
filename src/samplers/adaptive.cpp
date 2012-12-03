
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


// samplers/adaptive.cpp*
#include "stdafx.h"
#include "samplers/adaptive.h"
#include "paramset.h"
#include "film.h"
#include "primitive.h"
#include "intersection.h"
#include "camera.h"
#include "montecarlo.h"

// AdaptiveSampler Method Definitions
AdaptiveSampler::AdaptiveSampler(int xstart, int xend,
                     int ystart, int yend, int mins, int maxs, const string &m,
                     float sopen, float sclose)
    : Sampler(xstart, xend, ystart, yend, RoundUpPow2(max(mins, maxs)),
              sopen, sclose) {
    xPos = xPixelStart;
    yPos = yPixelStart;
    supersamplePixel = false;
    if (mins > maxs) std::swap(mins, maxs);

    if (!IsPowerOf2(mins)) {
        Warning("Minimum pixel samples being rounded up to power of 2");
        minSamples = RoundUpPow2(mins);
    }
    else
        minSamples = mins;
    if (!IsPowerOf2(maxs)) {
        Warning("Maximum pixel samples being rounded up to power of 2");
        maxSamples = RoundUpPow2(maxs);
    }
    else
        maxSamples = maxs;

    if (minSamples < 2) {
        Warning("Adaptive sampler needs at least two initial pixel samples.  Using two.");
        minSamples = 2;
    }
    if (minSamples == maxSamples) {
        maxSamples *= 2;
        Warning("Adaptive sampler must have more maximum samples than minimum.  Using %d - %d",
                minSamples, maxSamples);
    }
    if (m == "contrast") method = ADAPTIVE_CONTRAST_THRESHOLD;
    else if (m == "shapeid") method = ADAPTIVE_COMPARE_SHAPE_ID;
    else {
        Warning("Adaptive sampling metric \"%s\" unknown.  Using \"contrast\".",
                m.c_str());
        method = ADAPTIVE_CONTRAST_THRESHOLD;
    }
    sampleBuf = NULL;
}


AdaptiveSampler::~AdaptiveSampler() {
    delete[] sampleBuf;
}


Sampler *AdaptiveSampler::GetSubSampler(int num, int count) {
    int x0, x1, y0, y1;
    ComputeSubWindow(num, count, &x0, &x1, &y0, &y1);
    if (x0 == x1 || y0 == y1) return NULL;
    return new AdaptiveSampler(x0, x1, y0, y1, minSamples, maxSamples,
        method == ADAPTIVE_CONTRAST_THRESHOLD ? "contrast" : "shapeid",
        shutterOpen, shutterClose);
}


int AdaptiveSampler::GetMoreSamples(Sample *samples, RNG &rng) {
    if (!sampleBuf)
        sampleBuf = new float[LDPixelSampleFloatsNeeded(samples,
                                                        maxSamples)];
    if (supersamplePixel) {
        LDPixelSample(xPos, yPos, shutterOpen, shutterClose, maxSamples,
                      samples, sampleBuf, rng);
        return maxSamples;
    }
    else {
        if (yPos == yPixelEnd) return 0;
        LDPixelSample(xPos, yPos, shutterOpen, shutterClose, minSamples,
                      samples, sampleBuf, rng);
        return minSamples;
    }
}


bool AdaptiveSampler::ReportResults(Sample *samples,
        const RayDifferential *rays, const Spectrum *Ls,
        const Intersection *isects, int count) {
    if (supersamplePixel) {
        supersamplePixel = false;
        // Advance to next pixel for sampling for _AdaptiveSampler_
        if (++xPos == xPixelEnd) {
            xPos = xPixelStart;
            ++yPos;
        }
        return true;
    }
    else if (needsSupersampling(samples, rays, Ls, isects, count)) {
        PBRT_SUPERSAMPLE_PIXEL_YES(xPos, yPos);
        supersamplePixel = true;
        return false;
    }

    else {
        PBRT_SUPERSAMPLE_PIXEL_NO(xPos, yPos);
        // Advance to next pixel for sampling for _AdaptiveSampler_
        if (++xPos == xPixelEnd) {
            xPos = xPixelStart;
            ++yPos;
        }
        return true;
    }
}


bool AdaptiveSampler::needsSupersampling(Sample *samples,
        const RayDifferential *rays, const Spectrum *Ls,
        const Intersection *isects, int count) {
    switch (method) {
    case ADAPTIVE_COMPARE_SHAPE_ID:
        // See if any shape ids differ within samples
        for (int i = 0; i < count-1; ++i)
            if (isects[i].shapeId != isects[i+1].shapeId ||
                isects[i].primitiveId != isects[i+1].primitiveId)
                return true;
        return false;
    case ADAPTIVE_CONTRAST_THRESHOLD:
        // Compare contrast of sample differences to threshold
        float Lavg = 0.f;
        for (int i = 0; i < count; ++i)
            Lavg += Ls[i].y();
        Lavg /= count;
        const float maxContrast = 0.5f;
        for (int i = 0; i < count; ++i)
            if (fabsf(Ls[i].y() - Lavg) / Lavg > maxContrast)
                return true;
        return false;
    }
    return false;
}


AdaptiveSampler *CreateAdaptiveSampler(const ParamSet &params, const Film *film,
        const Camera *camera) {
    // Initialize common sampler parameters
    int xstart, xend, ystart, yend;
    film->GetSampleExtent(&xstart, &xend, &ystart, &yend);
    int minsamp = params.FindOneInt("minsamples", 4);
    int maxsamp = params.FindOneInt("maxsamples", 32);
    if (PbrtOptions.quickRender) { minsamp = 2; maxsamp = 4; }
    string method = params.FindOneString("method", "contrast");
    return new AdaptiveSampler(xstart, xend, ystart, yend, minsamp, maxsamp, method,
         camera->shutterOpen, camera->shutterClose);
}


