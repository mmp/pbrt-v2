
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


// core/sampler.cpp*
#include "stdafx.h"
#include "sampler.h"
#include "integrator.h"
#include "volume.h"

// Sampler Method Definitions
Sampler::~Sampler() {
}


Sampler::Sampler(int xstart, int xend, int ystart, int yend, int spp,
                 float sopen, float sclose)
    : xPixelStart(xstart), xPixelEnd(xend), yPixelStart(ystart),
      yPixelEnd(yend), samplesPerPixel(spp), shutterOpen(sopen),
      shutterClose(sclose) { }
bool Sampler::ReportResults(Sample *samples, const RayDifferential *rays,
        const Spectrum *Ls, const Intersection *isects, int count) {
    return true;
}


void Sampler::ComputeSubWindow(int num, int count, int *newXStart,
        int *newXEnd, int *newYStart, int *newYEnd) const {
    // Determine how many tiles to use in each dimension, _nx_ and _ny_
    int dx = xPixelEnd - xPixelStart, dy = yPixelEnd - yPixelStart;
    int nx = count, ny = 1;
    while ((nx & 0x1) == 0 && 2 * dx * ny < dy * nx) {
        nx >>= 1;
        ny <<= 1;
    }
    Assert(nx * ny == count);

    // Compute $x$ and $y$ pixel sample range for sub-window
    int xo = num % nx, yo = num / nx;
    float tx0 = float(xo) / float(nx), tx1 = float(xo+1) / float(nx);
    float ty0 = float(yo) / float(ny), ty1 = float(yo+1) / float(ny);
    *newXStart = Floor2Int(Lerp(tx0, xPixelStart, xPixelEnd));
    *newXEnd   = Floor2Int(Lerp(tx1, xPixelStart, xPixelEnd));
    *newYStart = Floor2Int(Lerp(ty0, yPixelStart, yPixelEnd));
    *newYEnd   = Floor2Int(Lerp(ty1, yPixelStart, yPixelEnd));
}



// Sample Method Definitions
Sample::Sample(Sampler *sampler, SurfaceIntegrator *surf,
               VolumeIntegrator *vol, const Scene *scene) {
    if (surf) surf->RequestSamples(sampler, this, scene);
    if (vol)  vol->RequestSamples(sampler, this, scene);
    AllocateSampleMemory();
}


void Sample::AllocateSampleMemory() {
    // Allocate storage for sample pointers
    int nPtrs = n1D.size() + n2D.size();
    if (!nPtrs) {
        oneD = twoD = NULL;
        return;
    }
    oneD = AllocAligned<float *>(nPtrs);
    twoD = oneD + n1D.size();

    // Compute total number of sample values needed
    int totSamples = 0;
    for (uint32_t i = 0; i < n1D.size(); ++i)
        totSamples += n1D[i];
    for (uint32_t i = 0; i < n2D.size(); ++i)
        totSamples += 2 * n2D[i];

    // Allocate storage for sample values
    float *mem = AllocAligned<float>(totSamples);
    for (uint32_t i = 0; i < n1D.size(); ++i) {
        oneD[i] = mem;
        mem += n1D[i];
    }
    for (uint32_t i = 0; i < n2D.size(); ++i) {
        twoD[i] = mem;
        mem += 2 * n2D[i];
    }
}


Sample *Sample::Duplicate(int count) const {
    Sample *ret = new Sample[count];
    for (int i = 0; i < count; ++i) {
        ret[i].n1D = n1D;
        ret[i].n2D = n2D;
        ret[i].AllocateSampleMemory();
    }
    return ret;
}


