
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

#ifndef PBRT_FILM_IMAGE_H
#define PBRT_FILM_IMAGE_H

// film/image.h*
#include "pbrt.h"
#include "film.h"
#include "sampler.h"
#include "filter.h"
#include "paramset.h"

// ImageFilm Declarations
class ImageFilm : public Film {
public:
    // ImageFilm Public Methods
    ImageFilm(int xres, int yres, Filter *filt, const float crop[4],
              const string &filename, bool openWindow);
    ~ImageFilm() {
        delete pixels;
        delete filter;
        delete[] filterTable;
    }
    void AddSample(const CameraSample &sample, const Spectrum &L);
    void Splat(const CameraSample &sample, const Spectrum &L);
    void GetSampleExtent(int *xstart, int *xend, int *ystart, int *yend) const;
    void GetPixelExtent(int *xstart, int *xend, int *ystart, int *yend) const;
    void WriteImage(float splatScale);
    void UpdateDisplay(int x0, int y0, int x1, int y1, float splatScale);
private:
    // ImageFilm Private Data
    Filter *filter;
    float cropWindow[4];
    string filename;
    int xPixelStart, yPixelStart, xPixelCount, yPixelCount;
    struct Pixel {
        Pixel() {
            for (int i = 0; i < 3; ++i) Lxyz[i] = splatXYZ[i] = 0.f;
            weightSum = 0.f;
        }
        float Lxyz[3];
        float weightSum;
        float splatXYZ[3];
        float pad;
    };
    BlockedArray<Pixel> *pixels;
    float *filterTable;
};


ImageFilm *CreateImageFilm(const ParamSet &params, Filter *filter);

#endif // PBRT_FILM_IMAGE_H
