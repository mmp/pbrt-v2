
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


// film/image.cpp*
#include "film/image.h"
#include "spectrum.h"
#include "parallel.h"
#include "imageio.h"
#ifdef PBRT_HAS_LIBSDL
#include "SDL.h"
#endif // PBRT_HAS_LIBSDL

// ImageFilm Method Definitions
ImageFilm::ImageFilm(int xres, int yres, Filter *filt, const float crop[4],
                     const string &fn, bool openWindow)
    : Film(xres, yres) {
    filter = filt;
    memcpy(cropWindow, crop, 4 * sizeof(float));
    filename = fn;
    // Compute film image extent
    xPixelStart = Ceil2Int(xResolution * cropWindow[0]);
    xPixelCount = max(1, Ceil2Int(xResolution * cropWindow[1]) - xPixelStart);
    yPixelStart = Ceil2Int(yResolution * cropWindow[2]);
    yPixelCount = max(1, Ceil2Int(yResolution * cropWindow[3]) - yPixelStart);

    // Allocate film image storage
    pixels = new BlockedArray<Pixel>(xPixelCount, yPixelCount);

    // Precompute filter weight table
#define FILTER_TABLE_SIZE 16
    filterTable = new float[FILTER_TABLE_SIZE * FILTER_TABLE_SIZE];
    float *ftp = filterTable;
    for (int y = 0; y < FILTER_TABLE_SIZE; ++y) {
        float fy = ((float)y + .5f) *
                   filter->yWidth / FILTER_TABLE_SIZE;
        for (int x = 0; x < FILTER_TABLE_SIZE; ++x) {
            float fx = ((float)x + .5f) *
                       filter->xWidth / FILTER_TABLE_SIZE;
            *ftp++ = filter->Evaluate(fx, fy);
        }
    }

    // Possibly open window for image display
    if (openWindow || PbrtOptions.openWindow) {
#ifdef PBRT_HAS_LIBSDL
        SDL_Init(SDL_INIT_VIDEO);
        SDL_WM_SetCaption(filename.c_str(), filename.c_str());
        sdlWindow = SDL_SetVideoMode(xPixelCount, yPixelCount, 32, SDL_SWSURFACE);
        if (sdlWindow == NULL)
            Error("Unable to create window to display image");
        if (SDL_LockSurface(sdlWindow) == -1)
            Error("Unable to lock surface for image display");
        else {
            for (int y = 0; y < yPixelCount; ++y) {
                for (int x = 0; x < xPixelCount; ++x) {
                    uint32_t *bufp = (uint32_t *)sdlWindow->pixels + y*sdlWindow->pitch/4 + x;
                    *bufp = SDL_MapRGB(sdlWindow->format, 64, 64, 64);
                }
            }
            SDL_UnlockSurface(sdlWindow);
            SDL_UpdateRect(sdlWindow, 0, 0, xPixelCount, yPixelCount);
        }
#else
        Warning("Support for opening image display window not available in this build.");
#endif // PBRT_HAS_LIBSDL
    }
#ifdef PBRT_HAS_LIBSDL
    else
        sdlWindow = NULL;
#endif // PBRT_HAS_LIBSDL
}


void ImageFilm::AddSample(const CameraSample &sample,
                          const Spectrum &L) {
    // Compute sample's raster extent
    float dimageX = sample.imageX - 0.5f;
    float dimageY = sample.imageY - 0.5f;
    int x0 = Ceil2Int (dimageX - filter->xWidth);
    int x1 = Floor2Int(dimageX + filter->xWidth);
    int y0 = Ceil2Int (dimageY - filter->yWidth);
    int y1 = Floor2Int(dimageY + filter->yWidth);
    x0 = max(x0, xPixelStart);
    x1 = min(x1, xPixelStart + xPixelCount - 1);
    y0 = max(y0, yPixelStart);
    y1 = min(y1, yPixelStart + yPixelCount - 1);
    if ((x1-x0) < 0 || (y1-y0) < 0)
    {
        PBRT_SAMPLE_OUTSIDE_IMAGE_EXTENT(const_cast<CameraSample *>(&sample));
        return;
    }

    // Loop over filter support and add sample to pixel arrays
    float xyz[3];
    L.ToXYZ(xyz);

    // Precompute $x$ and $y$ filter table offsets
    int *ifx = ALLOCA(int, x1 - x0 + 1);
    for (int x = x0; x <= x1; ++x) {
        float fx = fabsf((x - dimageX) *
                         filter->invXWidth * FILTER_TABLE_SIZE);
        ifx[x-x0] = min(Floor2Int(fx), FILTER_TABLE_SIZE-1);
    }
    int *ify = ALLOCA(int, y1 - y0 + 1);
    for (int y = y0; y <= y1; ++y) {
        float fy = fabsf((y - dimageY) *
                         filter->invYWidth * FILTER_TABLE_SIZE);
        ify[y-y0] = min(Floor2Int(fy), FILTER_TABLE_SIZE-1);
    }
    bool syncNeeded = (filter->xWidth > 0.5f || filter->yWidth > 0.5f);
    for (int y = y0; y <= y1; ++y) {
        for (int x = x0; x <= x1; ++x) {
            // Evaluate filter value at $(x,y)$ pixel
            int offset = ify[y-y0]*FILTER_TABLE_SIZE + ifx[x-x0];
            float filterWt = filterTable[offset];

            // Update pixel values with filtered sample contribution
            Pixel &pixel = (*pixels)(x - xPixelStart, y - yPixelStart);
            if (!syncNeeded) {
                pixel.Lxyz[0] += filterWt * xyz[0];
                pixel.Lxyz[1] += filterWt * xyz[1];
                pixel.Lxyz[2] += filterWt * xyz[2];
                pixel.weightSum += filterWt;
            }
            else {
                // Safely update _Lxyz_ and _weightSum_ even with concurrency
                AtomicAdd(&pixel.Lxyz[0], filterWt * xyz[0]);
                AtomicAdd(&pixel.Lxyz[1], filterWt * xyz[1]);
                AtomicAdd(&pixel.Lxyz[2], filterWt * xyz[2]);
                AtomicAdd(&pixel.weightSum, filterWt);
            }
        }
    }
}


void ImageFilm::Splat(const CameraSample &sample, const Spectrum &L) {
    float xyz[3];
    L.ToXYZ(xyz);
    int x = Floor2Int(sample.imageX), y = Floor2Int(sample.imageY);
    if (x < xPixelStart || x - xPixelStart >= xPixelCount ||
        y < yPixelStart || y - yPixelStart >= yPixelCount) return;
    Pixel &pixel = (*pixels)(x - xPixelStart, y - yPixelStart);
    AtomicAdd(&pixel.splatXYZ[0], xyz[0]);
    AtomicAdd(&pixel.splatXYZ[1], xyz[1]);
    AtomicAdd(&pixel.splatXYZ[2], xyz[2]);
}


void ImageFilm::GetSampleExtent(int *xstart, int *xend,
                                int *ystart, int *yend) const {
    *xstart = Floor2Int(xPixelStart + 0.5f - filter->xWidth);
    *xend   = Floor2Int(xPixelStart + 0.5f + xPixelCount  +
                        filter->xWidth);
    *ystart = Floor2Int(yPixelStart + 0.5f - filter->yWidth);
    *yend   = Floor2Int(yPixelStart + 0.5f + yPixelCount +
                        filter->yWidth);
}


void ImageFilm::GetPixelExtent(int *xstart, int *xend,
                               int *ystart, int *yend) const {
    *xstart = xPixelStart;
    *xend   = xPixelStart + xPixelCount;
    *ystart = yPixelStart;
    *yend   = yPixelStart + yPixelCount;
}


void ImageFilm::WriteImage(float splatScale) {
    // Convert image to RGB and compute final pixel values
    int nPix = xPixelCount * yPixelCount;
    float *rgb = new float[3*nPix];
    int offset = 0;
    for (int y = 0; y < yPixelCount; ++y) {
        for (int x = 0; x < xPixelCount; ++x) {
            // Convert pixel XYZ color to RGB
            XYZToRGB((*pixels)(x, y).Lxyz, &rgb[3*offset]);

            // Normalize pixel with weight sum
            float weightSum = (*pixels)(x, y).weightSum;
            if (weightSum != 0.f) {
                float invWt = 1.f / weightSum;
                rgb[3*offset  ] = max(0.f, rgb[3*offset  ] * invWt);
                rgb[3*offset+1] = max(0.f, rgb[3*offset+1] * invWt);
                rgb[3*offset+2] = max(0.f, rgb[3*offset+2] * invWt);
            }

            // Add splat value at pixel
            float splatRGB[3];
            XYZToRGB((*pixels)(x,y).splatXYZ, splatRGB);
            rgb[3*offset  ] += splatScale * splatRGB[0];
            rgb[3*offset+1] += splatScale * splatRGB[1];
            rgb[3*offset+2] += splatScale * splatRGB[2];
            ++offset;
        }
    }

    // Write RGB image
    ::WriteImage(filename, rgb, NULL, xPixelCount, yPixelCount,
                 xResolution, yResolution, xPixelStart, yPixelStart);

    // Release temporary image memory
    delete[] rgb;
}


void ImageFilm::UpdateDisplay(int x0, int y0, int x1, int y1,
    float splatScale) {
#ifdef PBRT_HAS_LIBSDL
    if (!sdlWindow) return;
    // Compute window coordinates for pixels to update
    x0 -= xPixelStart;
    x1 -= xPixelStart;
    y0 -= yPixelStart;
    y1 -= yPixelStart;
    x0 = Clamp(x0, 0, xPixelCount);
    x1 = Clamp(x1, 0, xPixelCount);
    y0 = Clamp(y0, 0, yPixelCount);
    y1 = Clamp(y1, 0, yPixelCount);
    uint32_t *pix = new uint32_t[(x1-x0)*(y1-y0)];
    uint32_t *pp = pix;
    for (int y = y0; y < y1; ++y) {
        for (int x = x0; x < x1; ++x) {
            // Compute weighted pixel value and update window pixel
            Pixel &pixel = (*pixels)(x, y);
            float rgb[3];
            XYZToRGB(pixel.Lxyz, rgb);
            float weightSum = pixel.weightSum;
            if (weightSum != 0.f) {
                float invWt = 1.f / weightSum;
                rgb[0] *= invWt;
                rgb[1] *= invWt;
                rgb[2] *= invWt;
            }
            float splatRGB[3];
            XYZToRGB(pixel.splatXYZ, splatRGB);
            rgb[0] += splatScale * splatRGB[0];
            rgb[1] += splatScale * splatRGB[1];
            rgb[2] += splatScale * splatRGB[2];
            
            *pp++ = SDL_MapRGB(sdlWindow->format,
                uint8_t(Clamp(powf(rgb[0], 1./1.8), 0.f, 1.f) * 255),
                uint8_t(Clamp(powf(rgb[1], 1./1.8), 0.f, 1.f) * 255),
                uint8_t(Clamp(powf(rgb[2], 1./1.8), 0.f, 1.f) * 255));
        }
    }
    // Update window pixels, redraw
    if (SDL_LockSurface(sdlWindow) == -1) { }
    pp=pix;
    for (int y = y0; y < y1; ++y) {
        for (int x = x0; x < x1; ++x) {
            uint32_t *bufp = (uint32_t *)sdlWindow->pixels + y*sdlWindow->pitch/4 + x;
            *bufp = *pp++;
        }
    }
    SDL_UnlockSurface(sdlWindow);
    SDL_UpdateRect(sdlWindow, x0, y0, x1-x0, y1-y0);
    delete[] pix;
#if 0
    SDL_Event event;
    while (SDL_PollEvent(&event)) {
        if (event.type == SDL_QUIT)
            SDL_Quit();
    }
#endif
#endif // PBRT_HAS_LIBSDL
}


ImageFilm *CreateImageFilm(const ParamSet &params, Filter *filter) {
    string filename = params.FindOneString("filename",
#ifdef PBRT_NO_OPENEXR
        "pbrt.tga"
#else
        "pbrt.exr"
#endif
                  );

    int xres = params.FindOneInt("xresolution", 640);
    int yres = params.FindOneInt("yresolution", 480);
    if (PbrtOptions.quickRender) xres = max(1, xres / 4);
    if (PbrtOptions.quickRender) yres = max(1, yres / 4);
    bool openwin = params.FindOneBool("display", false);
    float crop[4] = { 0, 1, 0, 1 };
    int cwi;
    const float *cr = params.FindFloat("cropwindow", &cwi);
    if (cr && cwi == 4) {
        crop[0] = Clamp(min(cr[0], cr[1]), 0., 1.);
        crop[1] = Clamp(max(cr[0], cr[1]), 0., 1.);
        crop[2] = Clamp(min(cr[2], cr[3]), 0., 1.);
        crop[3] = Clamp(max(cr[2], cr[3]), 0., 1.);
    }

    return new ImageFilm(xres, yres, filter, crop,
        filename, openwin);
}


