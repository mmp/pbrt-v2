
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


// core/imageio.cpp*
#include "imageio.h"
#include "spectrum.h"

// ImageIO Local Declarations
static RGBSpectrum *ReadImageEXR(const string &name, int *width, int *height);
static void WriteImageEXR(const string &name, float *pixels,
        float *alpha, int xRes, int yRes,
        int totalXRes, int totalYRes,
        int xOffset, int yOffset);
static void WriteImageTGA(const string &name, float *pixels,
        float *alpha, int xRes, int yRes,
        int totalXRes, int totalYRes,
        int xOffset, int yOffset);
static RGBSpectrum *ReadImageTGA(const string &name, int *w, int *h);

// ImageIO Function Definitions
RGBSpectrum *ReadImage(const string &name, int *width, int *height) {
    if (name.size() >= 5) {
        uint32_t suffixOffset = name.size() - 4;
#ifdef PBRT_HAS_OPENEXR
        if (!strcmp(name.c_str() + suffixOffset, ".exr") ||
            !strcmp(name.c_str() + suffixOffset, ".EXR"))
             return ReadImageEXR(name, width, height);
        else
#endif // PBRT_HAS_OPENEXR
              if (!strcmp(name.c_str() + suffixOffset, ".tga") ||
                 !strcmp(name.c_str() + suffixOffset, ".TGA"))
             return ReadImageTGA(name, width, height);
    }
    Error("Can't determine image file type from suffix of filename \"%s\"",
          name.c_str());
    RGBSpectrum *ret = new RGBSpectrum[1];
    ret[0] = 0.5f;
    *width = *height = 1;
    return ret;
}


void WriteImage(const string &name, float *pixels,
        float *alpha, int xRes, int yRes,
        int totalXRes, int totalYRes,
        int xOffset, int yOffset) {
    if (name.size() >= 5) {
        uint32_t suffixOffset = name.size() - 4;
#ifdef PBRT_HAS_OPENEXR
        if (!strcmp(name.c_str() + suffixOffset, ".exr") ||
            !strcmp(name.c_str() + suffixOffset, ".EXR")) {
             WriteImageEXR(name, pixels, alpha, xRes, yRes, totalXRes,
                           totalYRes, xOffset, yOffset);
             return;
        }
        else
#endif // PBRT_HAS_OPENEXR
             if (!strcmp(name.c_str() + suffixOffset, ".tga") ||
                 !strcmp(name.c_str() + suffixOffset, ".TGA")) {
             WriteImageTGA(name, pixels, alpha, xRes, yRes, totalXRes,
                           totalYRes, xOffset, yOffset);
             return;
        }
    }
    Error("Can't determine image file type from suffix of filename \"%s\"",
          name.c_str());
}


#ifdef PBRT_HAS_OPENEXR
#ifdef WIN32
#define hypotf hypot // For the OpenEXR headers
#endif
#include <ImfInputFile.h>
#include <ImfRgbaFile.h>
#include <ImfChannelList.h>
#include <ImfFrameBuffer.h>
#include <half.h>
using namespace Imf;
using namespace Imath;

// EXR Function Definitions
static RGBSpectrum *ReadImageEXR(const string &name, int *width, int *height) {
    try {
    InputFile file(name.c_str());
    Box2i dw = file.header().dataWindow();
    *width  = dw.max.x - dw.min.x + 1;
    *height = dw.max.y - dw.min.y + 1;

    half *rgb = new half[3 * *width * *height];

    FrameBuffer frameBuffer;
    frameBuffer.insert("R", Slice(HALF, (char *)rgb,
        3*sizeof(half), *width * 3 * sizeof(half), 1, 1, 0.0));
    frameBuffer.insert("G", Slice(HALF, (char *)rgb+sizeof(half),
        3*sizeof(half), *width * 3 * sizeof(half), 1, 1, 0.0));
    frameBuffer.insert("B", Slice(HALF, (char *)rgb+2*sizeof(half),
        3*sizeof(half), *width * 3 * sizeof(half), 1, 1, 0.0));

    file.setFrameBuffer(frameBuffer);
    file.readPixels(dw.min.y, dw.max.y);

    RGBSpectrum *ret = new RGBSpectrum[*width * *height];
    for (int i = 0; i < *width * *height; ++i) {
        float frgb[3] = { rgb[3*i], rgb[3*i+1], rgb[3*i+2] };
        ret[i] = RGBSpectrum::FromRGB(frgb);
    }
    delete[] rgb;
    return ret;
    } catch (const std::exception &e) {
        Error("Unable to read image file \"%s\": %s", name.c_str(),
            e.what());
        return NULL;
    }
}


static void WriteImageEXR(const string &name, float *pixels,
        float *alpha, int xRes, int yRes,
        int totalXRes, int totalYRes,
        int xOffset, int yOffset) {
    Rgba *hrgba = new Rgba[xRes * yRes];
    for (int i = 0; i < xRes * yRes; ++i)
        hrgba[i] = Rgba(pixels[3*i], pixels[3*i+1], pixels[3*i+2],
                        alpha ? alpha[i]: 1.f);

    Box2i displayWindow(V2i(0,0), V2i(totalXRes-1, totalYRes-1));
    Box2i dataWindow(V2i(xOffset, yOffset), V2i(xOffset + xRes - 1, yOffset + yRes - 1));

    RgbaOutputFile file(name.c_str(), displayWindow, dataWindow, WRITE_RGBA);
    file.setFrameBuffer(hrgba - xOffset - yOffset * xRes, 1, xRes);
    try {
        file.writePixels(yRes);
    }
    catch (const std::exception &e) {
        Error("Unable to write image file \"%s\": %s", name.c_str(),
            e.what());
    }

    delete[] hrgba;
}


#endif // PBRT_HAS_OPENEXR

// TGA Function Definitions
/**\file
 *\section License
 * License: GPL
 * Online License Link: http://www.gnu.org/licenses/gpl.html
 *
 *\author Copyright (c) 2003-2009 Jaakko Keranen <jaakko.keranen@iki.fi>
 *\author Copyright (c) 2009 Daniel Swanson <danij@dengine.net>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

/**
 * gl_tga.c: TGA file format (TARGA) reader/writer.
 */

// HEADER FILES ------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef WIN32
typedef unsigned char byte;
#else
#include <stdint.h>
typedef uint8_t byte;
#endif

typedef unsigned char uchar;

// MACROS ------------------------------------------------------------------

#undef SHORT
#ifdef __BIG_ENDIAN__
#define SHORT(x)            shortSwap(x)
# else // Little-endian.
#define SHORT(x)            (x)
#endif

// TYPES -------------------------------------------------------------------

typedef struct {
    uchar           idLength; // Identification field size in bytes.
    uchar           colorMapType; // Type of the color map.
    uchar           imageType; // Image type code.
} tga_header_t;

// Color map specification.
typedef struct {
    int16_t         index; // Index of first color map entry.
    int16_t         length; // Number of color map entries.
    uchar           entrySize; // Number of bits in a color map entry (16/24/32).
} tga_colormapspec_t;

// Image specification flags:
#define ISF_SCREEN_ORIGIN_UPPER 0x1 // Upper left-hand corner screen origin.
// Data interleaving:
#define ISF_INTERLEAVE_TWOWAY   0x2 // Two-way (even/odd) interleaving.
#define ISF_INTERLEAVE_FOURWAY  0x4 // Four-way interleaving.

// Image specification.
typedef struct {
    uchar           flags;
    int16_t         xOrigin; // X coordinate of lower left corner.
    int16_t         yOrigin; // Y coordinate of lower left corner.
    int16_t         width; // Width of the image in pixels.
    int16_t         height; // Height of the image in pixels.
    uchar           pixelDepth; // Number of bits in a pixel (16/24/32).
    uchar           attributeBits;
} tga_imagespec_t;


#ifdef __BIG_ENDIAN__
static int16_t shortSwap(int16_t n)
{
    return ((n & 0xff) << 8) | ((n & 0xff00) >> 8);
}


#endif

static bool writeByte(FILE* f, uchar b)
{
    return (fwrite(&b, 1, 1, f) == 1);
}



static bool writeShort(FILE* f, int16_t s)
{
    int16_t             v = SHORT(s);
    return (fwrite(&v, sizeof(v), 1, f) == 1);
}



static uchar readByte(FILE* f)
{
    uchar               v;
    fread(&v, sizeof(v), 1, f);
    return v;
}



static int16_t readShort(FILE* f)
{
    int16_t             v;
    fread(&v, sizeof(v), 1, f);
    return v;
}



/**
 * @param idLength      Identification field size in bytes (max 255).
 *                      @c 0 indicates no identification field.
 * @param colorMapType  Type of the color map, @c 0 or @c 1:
 *                      @c 0 = color map data is not present.
 *                      @c 1 = color map data IS present.
 * @param imageType     Image data type code, one of:
 *                      @c 0 = no image data is present.
 *                      @c 1 = uncompressed, color mapped image.
 *                      @c 2 = uncompressed, true-color image.
 *                      @c 3 = uncompressed, grayscale image.
 *                      @c 9 = run-length encoded, color mapped image.
 *                      @c 10 = run-length encoded, true-color image.
 *                      @c 11 = run-length encoded, grayscale image.
 * @param file          Handle to the file to be written to.
 */
static void writeHeader(uchar idLength, uchar colorMapType, uchar imageType,
                        FILE* file)
{
    writeByte(file, idLength);
    writeByte(file, colorMapType? 1 : 0);
    writeByte(file, imageType);
}



static void readHeader(tga_header_t* dst, FILE* file)
{
    dst->idLength = readByte(file);
    dst->colorMapType = readByte(file);
    dst->imageType = readByte(file);
}



/**
 * @param index         Index of first color map entry.
 * @param length        Total number of color map entries.
 * @param entrySize     Number of bits in a color map entry; 15/16/24/32.
 * @param file          Handle to the file to be written to.
 */
static void writeColorMapSpec(int16_t index, int16_t length,
                              uchar entrySize, FILE* file)
{
    writeShort(file, index);
    writeShort(file, length);
    writeByte(file, entrySize);
}



static void readColorMapSpec(tga_colormapspec_t* dst, FILE* file)
{
    dst->index = readShort(file);
    dst->length = readShort(file);
    dst->entrySize = readByte(file);
}



/**
 * @param xOrigin       X coordinate of lower left corner.
 * @param yOrigin       Y coordinate of lower left corner.
 * @param width         Width of the image in pixels.
 * @param height        Height of the image in pixels.
 * @param pixDepth      Number of bits per pixel, one of; 16/24/32.
 * @param file          Handle to the file to be written to.
 */
static void writeImageSpec(int16_t xOrigin, int16_t yOrigin,
                           int16_t width, int16_t height, uchar pixDepth,
                           FILE* file)
{
    writeShort(file, xOrigin);
    writeShort(file, yOrigin);
    writeShort(file, width);
    writeShort(file, height);
    writeByte(file, pixDepth);

    /**
     * attributeBits:4; // Attribute bits associated with each pixel.
     * reserved:1; // A reserved bit; must be 0.
     * screenOrigin:1; // Location of screen origin; must be 0.
     * dataInterleave:2; // TGA_INTERLEAVE_*
     */
    writeByte(file, 0);
}



static void readImageSpec(tga_imagespec_t* dst, FILE* file)
{
    uchar               bits;

    dst->xOrigin = readShort(file);
    dst->yOrigin = readShort(file);
    dst->width = readShort(file);
    dst->height = readShort(file);
    dst->pixelDepth = readByte(file);
    bits = readByte(file);

    /**
     * attributeBits:4; // Attribute bits associated with each pixel.
     * reserved:1; // A reserved bit; must be 0.
     * screenOrigin:1; // Location of screen origin; must be 0.
     * dataInterleave:2; // TGA_INTERLEAVE_*
     */
    dst->flags = 0 | ((bits & 6)? ISF_SCREEN_ORIGIN_UPPER : 0) |
        ((bits & 7)? ISF_INTERLEAVE_TWOWAY : (bits & 8)? ISF_INTERLEAVE_FOURWAY : 0);
    dst->attributeBits = (bits & 0xf);
}


/**
 * Save the rgb8888 buffer as Targa 24.
 *
 * @param filename      Path to the file to be written to (need not exist).
 * @param w             Width of the image in pixels.
 * @param h             Height of the image in pixels.
 * @param buf           Ptr to the image data to be written.
 *
 * @return              Non-zero iff successful.
 */
void WriteImageTGA(const string &name, float *pixels,
        float *alpha, int xRes, int yRes,
        int totalXRes, int totalYRes,
        int xOffset, int yOffset)
{
    FILE*               file;
    uchar*              outBuf;

    if ((file = fopen(name.c_str(), "wb")) == NULL) {
        Error("Unable to open output filename \"%s\"", name.c_str());
        return;
    }

    // No identification field, no color map, Targa type 2 (unmapped RGB).
    writeHeader(0, 0, 2, file);
    writeColorMapSpec(0, 0, 0, file);
    writeImageSpec(0, 0, xRes, yRes, 24, file);

    // The save format is BGR.
    outBuf = (uchar *)malloc(xRes * yRes * 3);
    uchar *dst = outBuf;
    for (int y = yRes-1; y >= 0; --y) {
        for (int x = 0; x < xRes; ++x) {
#define TO_BYTE(v) (uint8_t(Clamp(255.f * powf((v), 1.f/2.3f), 0.f, 255.f)))
            dst[0] = TO_BYTE(pixels[3*(y*xRes+x)+2]);
            dst[1] = TO_BYTE(pixels[3*(y*xRes+x)+1]);
            dst[2] = TO_BYTE(pixels[3*(y*xRes+x)+0]);
            dst += 3;
        }
    }
    if (fwrite(outBuf, 1, 3 * xRes * yRes, file) != uint32_t(3*xRes*yRes))
        Error("Error writing TGA image file \"%s\"", name.c_str());
    free(outBuf);
    fclose(file);
}




/**
 * Loads a TGA image (not the RLE types though)
 */
static RGBSpectrum *ReadImageTGA(const string &name, int *width, int *height)
{
    int                 x, y, pixbytes;
    tga_header_t        header;
    tga_colormapspec_t  colorMapSpec;
    tga_imagespec_t     imageSpec;
    uchar*              srcBuf;
    const uchar*        src;

    FILE *file = fopen(name.c_str(), "rb");
    if (!file) {
        Error("Unable to open TGA file \"%s\"", name.c_str());
        return NULL;
    }

    // Read and check the header.
    readHeader(&header, file);
    readColorMapSpec(&colorMapSpec, file);
    readImageSpec(&imageSpec, file);

    if ((imageSpec.attributeBits != 8 && imageSpec.attributeBits != 0) ||
        (imageSpec.flags & ISF_SCREEN_ORIGIN_UPPER) ||
        (header.imageType == 2 &&
          (imageSpec.pixelDepth != 32 && imageSpec.pixelDepth != 24)) ||
        (header.imageType == 3 &&
          (imageSpec.pixelDepth != 8)) ||
        (header.imageType != 2 && header.imageType != 3)) {
        Error("ReadImageTGA: I don't know this format "
              "(type=%i pxsize=%i abits=%i)", header.imageType,
              imageSpec.pixelDepth,
              imageSpec.attributeBits);
        return NULL;
    }

    *width = imageSpec.width;
    *height = imageSpec.height;

    // Determine format.
    if (imageSpec.pixelDepth == 32)
        pixbytes = 4;
    else if (imageSpec.pixelDepth == 24)
        pixbytes = 3;
    else {
        Assert(imageSpec.pixelDepth == 8);
        pixbytes = 1;
    }

    // Read the pixel data.
    int size = *width * *height * pixbytes;
    srcBuf = (uchar *)malloc(size);
    if (fread(srcBuf, 1, size, file) != (uint32_t)size) {
        Error("Premature end-of-file when reading TGA image \"%s\"", name.c_str());
        free(srcBuf);
        return NULL;
    }

    // "Unpack" the pixels (origin in the lower left corner).
    // TGA pixels are in BGRA format.
    src = srcBuf;
    RGBSpectrum *ret = new RGBSpectrum[*width * *height];
    RGBSpectrum *dst = ret;
    for(y = *height - 1; y >= 0; y--)
        for(x = 0; x < *width; x++) {
            if (pixbytes == 1)
                *dst++ = RGBSpectrum((*src++) / 255.f);
            else {
                float c[3];
                c[2] = (*src++) / 255.f;
                c[1] = (*src++) / 255.f;
                c[0] = (*src++) / 255.f;
                *dst++ = RGBSpectrum::FromRGB(c);
                if (pixbytes == 4)
                    ++src;
            }
    }
    free(srcBuf);
    return ret;
}


