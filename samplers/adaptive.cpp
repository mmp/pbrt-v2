
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


// samplers/adaptive.cpp*
#include "samplers/adaptive.h"
#include "paramset.h"
#include "film.h"
#include "primitive.h"
#include "camera.h"
#include "montecarlo.h"

// AdaptiveSampler Method Definitions
AdaptiveSampler::AdaptiveSampler(int xstart, int xend,
                     int ystart, int yend, int mins, int maxs,
                     float sopen, float sclose)
    : Sampler(xstart, xend, ystart, yend, RoundUpPow2(maxSamples),
              sopen, sclose) {
    xPos = xPixelStart;
    yPos = yPixelStart;
    superSamplePixel = false;
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
        ShutterOpen, ShutterClose);
}


int AdaptiveSampler::GetMoreSamples(Sample *samples) {
    if (!sampleBuf)
        sampleBuf = new float[LDPixelSampleFloatsNeeded(samples, maxSamples)];
    if (superSamplePixel) {
        LDPixelSample(xPos, yPos, ShutterOpen, ShutterClose, maxSamples, samples,
                      sampleBuf);
        return maxSamples;
    }
    else {
        if (yPos == yPixelEnd) return 0;
        LDPixelSample(xPos, yPos, ShutterOpen, ShutterClose, minSamples, samples,
                      sampleBuf);
        return minSamples;
    }
}


bool AdaptiveSampler::ReportResults(Sample *samples, const RayDifferential *rays,
        const Spectrum *Ls, const Intersection *isects,
        int count) {
    if (superSamplePixel) {
        superSamplePixel = false;
        if (++xPos == xPixelEnd) {
            xPos = xPixelStart;
            ++yPos;
        }
        return true;
    }
    else if (needsSupersampling(samples, rays, Ls, isects, count)) {
        PBRT_SUPERSAMPLE_PIXEL_YES(xPos, yPos);
        superSamplePixel = true;
        return false;
    }
    else {
        PBRT_SUPERSAMPLE_PIXEL_NO(xPos, yPos);
        if (++xPos == xPixelEnd) {
            xPos = xPixelStart;
            ++yPos;
        }
        return true;
    }
}


bool AdaptiveSampler::needsSupersampling(Sample *samples, const RayDifferential *rays,
        const Spectrum *Ls, const Intersection *isects,
        int count) {
#if 0
// diff primitives
    for (int i = 0; i < count; ++i)
        for (int j = i+1; j < count; ++j)
            if (isects[i].ShapeId != isects[j].ShapeId) return true;
    return false;
#else
// mitchell's contrast metric (I think)
    float Lavg = 0.f;
    for (int i = 0; i < count; ++i)
        Lavg += Ls[i].y();
    Lavg /= count;
    float maxc = 0.f;
    for (int i = 0; i < count; ++i)
        maxc = max(maxc, fabsf(Ls[i].y() - Lavg) / Lavg);
    return maxc > .5f;
#endif
}


AdaptiveSampler *CreateAdaptiveSampler(const ParamSet &params, const Film *film,
        const Camera *camera) {
    // Initialize common sampler parameters
    int xstart, xend, ystart, yend;
    film->GetSampleExtent(&xstart, &xend, &ystart, &yend);
    int minsamp = params.FindOneInt("minsamples", 4);
    int maxsamp = params.FindOneInt("maxsamples", 32);
    if (getenv("PBRT_QUICK_RENDER")) { minsamp = 2; maxsamp = 4; }
    return new AdaptiveSampler(xstart, xend, ystart, yend, minsamp, maxsamp,
         camera->ShutterOpen, camera->ShutterClose);
}


