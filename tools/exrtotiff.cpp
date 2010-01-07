/*
 * exrtoiff.cpp
 *
 * $Id: exrtotiff.cpp 1063 2004-05-14 01:23:50Z mmp $
 */

#include <stdio.h>
#include <stdlib.h>
#include <tiffio.h>
#include <ImfInputFile.h>
#include <ImfChannelList.h>
#include <ImfFrameBuffer.h>
#include <half.h>
#include <algorithm>

using namespace Imf;
using namespace Imath;

static bool ReadEXR(const char *name, float *&rgba, int &xRes, int &yRes, bool &hasAlpha);
static void WriteTIFF(const char *name, float *rgba, int xRes, int yRes, bool hasAlpha);

static void usage() {
  fprintf( stderr, "usage: exrtotiff [options] <input.exr> <output.tiff>\n" );
  fprintf( stderr, "Supported options:\n");
  fprintf( stderr, "\t-scale scale\n" );
  fprintf( stderr, "\t-gamma gamma\n" );
  fprintf( stderr, "\t-bg bg-grey color\n");
  exit(1);
}

int main(int argc, char *argv[]) 
{
    float scale = 1.f, gamma = 2.2f;
    float bggray = -1.f;

    int argNum = 1;
    while (argNum < argc && argv[argNum][0] == '-') {
#define ARG(name, var) \
      else if (!strcmp(argv[argNum], "-" name)) { \
	if (argNum+1 == argc) \
	  usage(); \
	var = atof(argv[argNum+1]); \
	++argNum; \
      }
        if (0) ;
	ARG("scale", scale)
	ARG("gamma", gamma)
	ARG("bg", bggray)
	else 
	    usage();
      ++argNum;

    }
    if (argNum + 2 > argc) usage();

    char *inFile = argv[argNum], *outFile = argv[argNum+1]; 	  
    float *rgba;
    int xRes, yRes;
    bool hasAlpha;

    if (ReadEXR(inFile, rgba, xRes, yRes, hasAlpha)) {
	float *rgb = new float[xRes*yRes*3];
	for (int i = 0; i < xRes*yRes; ++i) {
	    for (int j = 0; j < 3; ++j) {
		rgb[3*i+j] = scale * rgba[4*i+j];
//CO		if (rgba[4*i+3] != 0.f)
//CO		    rgb[3*i+j] /= rgba[4*i+3];
		if (bggray > 0)
		    rgb[3*i+j] = rgba[4*i+3] * rgb[3*i+j] + (1.f - rgba[4*i+3]) * bggray;
	    }
	    if (bggray > 0)
		rgba[4*i+3] = 1.f;
	}

	for (int i = 0; i < xRes*yRes; ++i) {
	    for (int j = 0; j < 3; ++j) {
		rgba[4*i+j] = 255.f * powf(std::max(0.f, rgb[3*i+j]), 1.f / gamma);
                if (rgba[4*i+j] < 0.f) rgba[4*i+j] = 0.f;
                if (rgba[4*i+j] > 255.f) rgba[4*i+j] = 255.f;
//CO		if (rgba[4*i+3] != 0.f)
//CO		    rgba[4*i+j] /= rgba[4*i+3];
	    }
	    rgba[4*i+3] *= 255.f;
	}

	WriteTIFF(outFile, rgba, xRes, yRes, hasAlpha);
    }
    return 0;
}

void WriteTIFF(const char *name, float *rgba, int XRes, int YRes, bool hasAlpha) 
{
    // Open 8-bit TIFF file for writing
    TIFF *tiff = TIFFOpen(name, "w");
    if (!tiff) {
	fprintf(stderr, "Unable to open TIFF %s for writing", name);
	return;
    }

    int nChannels = hasAlpha ? 4 : 3;
    TIFFSetField(tiff, TIFFTAG_SAMPLESPERPIXEL, nChannels);
    if (hasAlpha) {
	short int extra[] = { EXTRASAMPLE_ASSOCALPHA };
	TIFFSetField(tiff, TIFFTAG_EXTRASAMPLES, (short)1, extra);
    }
    // Write image resolution information
    TIFFSetField(tiff, TIFFTAG_IMAGEWIDTH, XRes);
    TIFFSetField(tiff, TIFFTAG_IMAGELENGTH, YRes);
    TIFFSetField(tiff, TIFFTAG_BITSPERSAMPLE, 8);
    TIFFSetField(tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
    // Set Generic TIFF Fields
    TIFFSetField(tiff, TIFFTAG_ROWSPERSTRIP, 1);
    TIFFSetField(tiff, TIFFTAG_XRESOLUTION, 1.f);
    TIFFSetField(tiff, TIFFTAG_YRESOLUTION, 1.f);
    TIFFSetField(tiff, TIFFTAG_RESOLUTIONUNIT, (short)1);
    TIFFSetField(tiff, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
    TIFFSetField(tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tiff, TIFFTAG_ORIENTATION, (int)ORIENTATION_TOPLEFT);
    // Write 8-bit scanlines
    unsigned char *buf = new unsigned char[nChannels * XRes];
    for (int y = 0; y < YRes; ++y) {
	unsigned char *bufp = buf;
	for (int x = 0; x < XRes; ++x) {
	    // Pack 8-bit pixels samples into buf
	    for (int s = 0; s < nChannels; ++s)
		*bufp++ = (unsigned char)*rgba++;
	}
	TIFFWriteScanline(tiff, buf, y, 1);
    }
    // Close 8-bit TIFF file
    delete[] buf;
    TIFFClose(tiff);
}

static bool ReadEXR(const char *name, float *&rgba, int &xRes, int &yRes, bool &hasAlpha)
{
    InputFile file(name);
    Box2i dw = file.header().dataWindow();
    xRes = dw.max.x - dw.min.x + 1;
    yRes = dw.max.y - dw.min.y + 1;

    half *hrgba = new half[4 * xRes * yRes];

    // for now...
    hasAlpha = true;
    int nChannels = 4;

    half *hp = hrgba - nChannels * (dw.min.x + dw.min.y * xRes);

    FrameBuffer frameBuffer;
    frameBuffer.insert("R", Slice(HALF, (char *)hp,
				  4*sizeof(half), xRes * 4 * sizeof(half), 1, 1, 0.0));
    frameBuffer.insert("G", Slice(HALF, (char *)hp+sizeof(half),
				  4*sizeof(half), xRes * 4 * sizeof(half), 1, 1, 0.0));
    frameBuffer.insert("B", Slice(HALF, (char *)hp+2*sizeof(half),
				  4*sizeof(half), xRes * 4 * sizeof(half), 1, 1, 0.0));
    frameBuffer.insert("A", Slice(HALF, (char *)hp+3*sizeof(half),
				  4*sizeof(half), xRes * 4 * sizeof(half), 1, 1, 1.0));

    file.setFrameBuffer(frameBuffer);
    file.readPixels(dw.min.y, dw.max.y);

    rgba = new float[nChannels * xRes * yRes];
    for (int i = 0; i < nChannels * xRes * yRes; ++i)
	rgba[i] = hrgba[i];
    delete[] hrgba;

    return rgba;
}
