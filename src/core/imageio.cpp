
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


// core/imageio.cpp*
#include "stdafx.h"
#include "imageio.h"
#include "spectrum.h"
#include "targa.h"
#include <string.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

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
static bool WriteImagePFM(const string &filename, const float *rgb, int xres, int yres);
static RGBSpectrum *ReadImagePFM(const string &filename, int *xres, int *yres);

// ImageIO Function Definitions
RGBSpectrum *ReadImage(const string &name, int *width, int *height) {
    if (name.size() >= 5) {
        uint32_t suffixOffset = name.size() - 4;
#ifdef PBRT_HAS_OPENEXR
        if (!strcmp(name.c_str() + suffixOffset, ".exr") ||
            !strcmp(name.c_str() + suffixOffset, ".EXR"))
             return ReadImageEXR(name, width, height);
#endif // PBRT_HAS_OPENEXR
        if (!strcmp(name.c_str() + suffixOffset, ".tga") ||
            !strcmp(name.c_str() + suffixOffset, ".TGA"))
            return ReadImageTGA(name, width, height);
        if (!strcmp(name.c_str() + suffixOffset, ".pfm") ||
            !strcmp(name.c_str() + suffixOffset, ".PFM"))
            return ReadImagePFM(name, width, height);
    }
    Error("Unable to load image stored in format \"%s\" for filename \"%s\". "
          "Returning a constant grey image instead.",
          strrchr(name.c_str(), '.') ? (strrchr(name.c_str(), '.') + 1) : "(unknown)",
          name.c_str());
    RGBSpectrum *ret = new RGBSpectrum[1];
    ret[0] = 0.5f;
    *width = *height = 1;
    return ret;
}


void WriteImage(const string &name, float *pixels, float *alpha, int xRes,
                int yRes, int totalXRes, int totalYRes, int xOffset, int yOffset) {
    if (name.size() >= 5) {
        uint32_t suffixOffset = name.size() - 4;
#ifdef PBRT_HAS_OPENEXR
        if (!strcmp(name.c_str() + suffixOffset, ".exr") ||
            !strcmp(name.c_str() + suffixOffset, ".EXR")) {
             WriteImageEXR(name, pixels, alpha, xRes, yRes, totalXRes,
                           totalYRes, xOffset, yOffset);
             return;
        }
#endif // PBRT_HAS_OPENEXR
        if (!strcmp(name.c_str() + suffixOffset, ".tga") ||
            !strcmp(name.c_str() + suffixOffset, ".TGA")) {
            WriteImageTGA(name, pixels, alpha, xRes, yRes, totalXRes,
                          totalYRes, xOffset, yOffset);
            return;
        }
        if (!strcmp(name.c_str() + suffixOffset, ".pfm") ||
            !strcmp(name.c_str() + suffixOffset, ".PFM")) {
            WriteImagePFM(name, pixels, xRes, yRes);
            return;
        }
        if (!strcmp(name.c_str() + suffixOffset, ".png") ||
            !strcmp(name.c_str() + suffixOffset, ".PNG")) {
            uint8_t *rgb8 = new uint8_t[3 * xRes * yRes];
            uint8_t *dst = rgb8;
            for (int y = 0; y < yRes; ++y) {
                for (int x = 0; x < xRes; ++x) {
#define TO_BYTE(v) (uint8_t(Clamp(255.f * powf((v), 1.f/2.2f), 0.f, 255.f)))
                    dst[0] = TO_BYTE(pixels[3*(y*xRes+x)+2]);
                    dst[1] = TO_BYTE(pixels[3*(y*xRes+x)+1]);
                    dst[2] = TO_BYTE(pixels[3*(y*xRes+x)+0]);
#undef TO_BYTE
                    dst += 3;
                }
            }
            if (stbi_write_png(name.c_str(), xRes, yRes, 3, rgb8,
                               3 * xRes) == 0)
                Error("Error writing PNG \"%s\"", name.c_str());
            delete[] rgb8;
            return;
        }
    }
    Error("Can't determine image file type from suffix of filename \"%s\"",
          name.c_str());
}


#ifdef PBRT_HAS_OPENEXR
#if defined(PBRT_IS_WINDOWS)
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
    RgbaInputFile file (name.c_str());
    Box2i dw = file.dataWindow();
    *width = dw.max.x - dw.min.x + 1;
    *height = dw.max.y - dw.min.y + 1;
    std::vector<Rgba> pixels(*width * *height);
    file.setFrameBuffer(&pixels[0] - dw.min.x - dw.min.y * *width, 1, *width);
    file.readPixels(dw.min.y, dw.max.y);

    RGBSpectrum *ret = new RGBSpectrum[*width * *height];
    for (int i = 0; i < *width * *height; ++i) {
        float frgb[3] = { pixels[i].r, pixels[i].g, pixels[i].b };
        ret[i] = RGBSpectrum::FromRGB(frgb);
    }
    Info("Read EXR image %s (%d x %d)", name.c_str(), *width, *height);
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

    try {
        RgbaOutputFile file(name.c_str(), displayWindow, dataWindow, WRITE_RGBA);
        file.setFrameBuffer(hrgba - xOffset - yOffset * xRes, 1, xRes);
        file.writePixels(yRes);
    }
    catch (const std::exception &e) {
        Error("Unable to write image file \"%s\": %s", name.c_str(),
            e.what());
    }

    delete[] hrgba;
}


#endif // PBRT_HAS_OPENEXR


void WriteImageTGA(const string &name, float *pixels,
                   float *alpha, int xRes, int yRes,
                   int totalXRes, int totalYRes,
                   int xOffset, int yOffset)
{
    // Reformat to BGR layout.
    uint8_t *outBuf = new uint8_t[3 * xRes * yRes];
    uint8_t *dst = outBuf;
    for (int y = 0; y < yRes; ++y) {
        for (int x = 0; x < xRes; ++x) {
#define TO_BYTE(v) (uint8_t(Clamp(255.f * powf((v), 1.f/2.2f), 0.f, 255.f)))
            dst[0] = TO_BYTE(pixels[3*(y*xRes+x)+2]);
            dst[1] = TO_BYTE(pixels[3*(y*xRes+x)+1]);
            dst[2] = TO_BYTE(pixels[3*(y*xRes+x)+0]);
            dst += 3;
        }
    }

    tga_result result;
    if ((result = tga_write_bgr(name.c_str(), outBuf, xRes, yRes, 24)) != TGA_NOERR)
        Error("Unable to write output file \"%s\" (%s)", name.c_str(),
              tga_error(result));

    delete[] outBuf;
}




static RGBSpectrum *ReadImageTGA(const string &name, int *width, int *height)
{
    tga_image img;
    tga_result result;
    if ((result = tga_read(&img, name.c_str())) != TGA_NOERR) {
        Error("Unable to read from TGA file \"%s\" (%s)", name.c_str(),
              tga_error(result));
        return NULL;
    }

    if (tga_is_right_to_left(&img))
        tga_flip_horiz(&img);
    if (!tga_is_top_to_bottom(&img))
        tga_flip_vert(&img);
    if (tga_is_colormapped(&img))
        tga_color_unmap(&img);

    *width = img.width;
    *height = img.height;

    // "Unpack" the pixels (origin in the lower left corner).
    // TGA pixels are in BGRA format.
    RGBSpectrum *ret = new RGBSpectrum[*width * *height];
    RGBSpectrum *dst = ret;
    for (int y = *height - 1; y >= 0; y--)
        for (int x = 0; x < *width; x++) {
            uint8_t *src = tga_find_pixel(&img, x, y);
            if (tga_is_mono(&img))
                *dst++ = RGBSpectrum(*src / 255.f);
            else {
                float c[3];
                c[2] = src[0] / 255.f;
                c[1] = src[1] / 255.f;
                c[0] = src[2] / 255.f;
                *dst++ = RGBSpectrum::FromRGB(c);
            }
    }

    tga_free_buffers(&img);
    Info("Read TGA image %s (%d x %d)", name.c_str(), *width, *height);

    return ret;
}



// PFM Function Definitions
/*
 * PFM reader/writer code courtesy Jiawen "Kevin" Chen (http://people.csail.mit.edu/jiawen/)
 */

static bool hostLittleEndian =
#if defined(__LITTLE_ENDIAN__) || defined(__i386__) || defined(__x86_64__) || defined(PBRT_IS_WINDOWS)
true
#elif defined(__BIG_ENDIAN__)
false
#else
#error "Can't detect machine endian-ness at compile-time."
#endif
    ;

#define BUFFER_SIZE 80

static inline int isWhitespace( char c )
{
    return c == ' ' || c == '\n' || c == '\t';
}



// reads a "word" from the fp and puts it into buffer
// and adds a null terminator
// i.e. it keeps reading until a whitespace is reached
// returns the number of characters read
// *not* including the whitespace
// return -1 on an error
static int readWord(FILE* fp, char* buffer, int bufferLength) {
    int n;
    char c;

    if (bufferLength < 1)
        return -1;

    n = 0;
    c = fgetc( fp );
    while( c != EOF && !isWhitespace( c ) && n < bufferLength ) {
        buffer[ n ] = c;
        ++n;
        c = fgetc( fp );
    }

    if( n < bufferLength ) {
        buffer[ n ] = '\0';
        return n;
    }

    return -1;
}



static RGBSpectrum *ReadImagePFM(const string &filename, int *xres, int *yres) {
    float *data = NULL;
    RGBSpectrum *rgb = NULL;
    char buffer[ BUFFER_SIZE ];
    unsigned int nFloats;
    int nChannels, width, height;
    float scale;
    bool fileLittleEndian;

    FILE *fp = fopen(filename.c_str(), "rb");
    if (!fp)
        goto fail;

    // read either "Pf" or "PF"
    if (readWord( fp, buffer, BUFFER_SIZE ) == -1)
        goto fail;

    if (strcmp( buffer, "Pf" ) == 0)
        nChannels = 1;
    else if (strcmp( buffer, "PF" ) == 0)
        nChannels = 3;
    else
        goto fail;

    // read the rest of the header
    // read width
    if (readWord( fp, buffer, BUFFER_SIZE ) == -1)
        goto fail;
    width = atoi( buffer );
    *xres = width;

    // read height
    if (readWord( fp, buffer, BUFFER_SIZE ) == -1)
        goto fail;
    height = atoi( buffer );
    *yres = height;

    // read scale
    if (readWord( fp, buffer, BUFFER_SIZE ) == -1)
        goto fail;
    sscanf( buffer, "%f", &scale );

    // read the data
    nFloats = nChannels * width * height;
    data = new float[nFloats];
    if (fread(data, sizeof( float ), nFloats, fp ) != nFloats)
        goto fail;

    // apply endian conversian and scale if appropriate
    fileLittleEndian = (scale < 0.f);
    if (hostLittleEndian ^ fileLittleEndian) {
        uint8_t bytes[4];
        for (unsigned int i = 0; i < nFloats; ++i) {
            memcpy(bytes, &data[i], 4);
            swap(bytes[0], bytes[3]);
            swap(bytes[1], bytes[2]);
            memcpy(&data[i], bytes, 4);
        }
    }
    if (fabsf(scale) != 1.f)
        for (unsigned int i = 0; i < nFloats; ++i)
            data[i] *= fabsf(scale);

    // create RGBs...
    rgb = new RGBSpectrum[width*height];
    if (nChannels == 1) {
        for (int i = 0; i < width * height; ++i)
            rgb[i] = data[i];
    }
    else {
        for (int i = 0; i < width * height; ++i)
            rgb[i] = RGBSpectrum::FromRGB(&data[3*i]);
    }

    delete[] data;
    fclose(fp);
    return rgb;

 fail:
    Error("Error reading PFM file \"%s\"", filename.c_str());
    fclose(fp);
    delete[] data;
    delete[] rgb;
    return NULL;
}




static bool WriteImagePFM(const string &filename, const float *rgb,
                          int width, int height) {
    FILE* fp;
    unsigned int nFloats;
    float scale;

    fp = fopen(filename.c_str(), "wb");
    if (!fp) {
        Error("Unable to open output PFM file \"%s\"", filename.c_str());
        return false;
    }

    // only write 3 channel PFMs here...
    if (fprintf(fp, "PF\n") < 0)
        goto fail;

    // write the width and height, which must be positive
    if (fprintf(fp, "%d %d\n", width, height) < 0)
        goto fail;

    // write the scale, which encodes endianness
    scale = hostLittleEndian ? -1.f : 1.f;
    if (fprintf(fp, "%f\n", scale) < 0)
        goto fail;

    // write the data from bottom left to upper right as specified by 
    // http://netpbm.sourceforge.net/doc/pfm.html
    // The raster is a sequence of pixels, packed one after another, with no
    // delimiters of any kind. They are grouped by row, with the pixels in each
    // row ordered left to right and the rows ordered bottom to top.
    nFloats = 3 * width * height;
    for (int j=height-1; j>=0; j--) {
        if (fwrite(rgb + j*width*3, sizeof(float), width*3, fp) < (size_t)(width*3))
            goto fail;
    }

    fclose(fp);
    return true;

 fail:
    Error("Error writing PFM file \"%s\"", filename.c_str());
    fclose(fp);
    return false;
}


