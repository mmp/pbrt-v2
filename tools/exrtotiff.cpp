/*
 * exrtoiff.cpp
 *
 * $Id: exrtotiff.cpp 1063 2004-05-14 01:23:50Z mmp $
 */

#include <stdio.h>
#include <stdlib.h>
#include <tiffio.h>
#include <assert.h>
#include <ImfInputFile.h>
#include <ImfChannelList.h>
#include <ImfFrameBuffer.h>
#include <half.h>
#include <algorithm>

using std::min;
using std::max;

using namespace Imf;
using namespace Imath;

static bool ReadEXR(const char *name, float *&rgba, int &xRes, int &yRes, bool &hasAlpha);
static void WriteTIFF(const char *name, float *rgba, int xRes, int yRes, bool hasAlpha);

static void usage() {
  fprintf( stderr, "usage: exrtotiff [options] <input.exr> <output.tiff>\n" );
  fprintf( stderr, "Supported options:\n");
  fprintf( stderr, "\t-scale scale\n" );
  fprintf( stderr, "\t-bloom\n" );
  fprintf( stderr, "\t-bloomscale scale [default: 0.1]\n" );
  fprintf( stderr, "\t-bloomradius radius [default: 0.01]\n" );
  fprintf( stderr, "\t-tonemap\n" );
  fprintf( stderr, "\t-repeatpix [count]\n" );
  fprintf( stderr, "\t-gamma gamma\n" );
  fprintf( stderr, "\t-bg bg-grey color\n");
  exit(1);
}

void reinhard(float *d, int xRes, int yRes) {
    const float yw[3] = { 0.212671f, 0.715160f, 0.072169f };
    // Compute world adaptation luminance, _Ywa_
    float Ywa = 0.;
    for (int i = 0; i < xRes * yRes; ++i) {
        float y = yw[0] * d[3*i] + yw[1] * d[3*i+1] + yw[2] * d[3*i+2];
//CO        if ((i % 1000) == 1) fprintf(stderr, "%f,%f,%f -> %f\n",
//CO                                     d[3*i], d[3*i+1], d[3*i+2], y);
        if (y > 1e-4f) Ywa += logf(y);
    }
    Ywa = expf(Ywa / (xRes * yRes));
//CO    Ywa /= 683.f;
//CO    fprintf(stderr, "YWA = %f\n", Ywa);
    float invY2 = 1.f / (Ywa * Ywa);
//CO    float invY2 = 1.f / (maxY * maxY);

    for (int i = 0; i < xRes * yRes; ++i) {
        float y = yw[0] * d[3*i] + yw[1] * d[3*i+1] + yw[2] * d[3*i+2];

        float s = (1.f + y * invY2) / (1.f + y);
//CO        if ((i % 1000) == 1) fprintf(stderr, "%f %f -> %f\n", y, invY2, s);
        d[3*i] *= s;
        d[3*i+1] *= s;
        d[3*i+2] *= s;
    }
}

inline float Lerp(float t, float a, float b) {
    return (1.f - t) * a + t * b;
}

void bloom(float *rgb, int xResolution, int yResolution, float bloomRadius,
           float bloomWeight) { 
    int nPix = xResolution * yResolution;
    // Possibly apply bloom effect to image
    if (bloomRadius > 0.f && bloomWeight > 0.f) {
        // Compute image-space extent of bloom effect
        int bloomSupport = int(ceilf(bloomRadius *
                                     max(xResolution, yResolution)));
        int bloomWidth = bloomSupport / 2;
        // Initialize bloom filter table
        float *bloomFilter = new float[bloomWidth * bloomWidth];
        for (int i = 0; i < bloomWidth * bloomWidth; ++i) {
            float dist = sqrtf(float(i)) / float(bloomWidth);
            bloomFilter[i] = powf(max(0.f, 1.f - dist), 4.f);
        }
        // Apply bloom filter to image pixels
        float *bloomImage = new float[3*nPix];
        for (int y = 0; y < yResolution; ++y) {
            for (int x = 0; x < xResolution; ++x) {
                // Compute bloom for pixel _(x,y)_
                // Compute extent of pixels contributing bloom
                int x0 = max(0, x - bloomWidth);
                int x1 = min(x + bloomWidth, xResolution - 1);
                int y0 = max(0, y - bloomWidth);
                int y1 = min(y + bloomWidth, yResolution - 1);
                int offset = y * xResolution + x;
                float sumWt = 0.;
                for (int by = y0; by <= y1; ++by)
                    for (int bx = x0; bx <= x1; ++bx) {
                        // Accumulate bloom from pixel $(bx,by)$
                        int dx = x - bx, dy = y - by;
                        if (dx == 0 && dy == 0) continue;
                        int dist2 = dx*dx + dy*dy;
                        if (dist2 < bloomWidth * bloomWidth) {
                            int bloomOffset = bx + by * xResolution;
                            float wt = bloomFilter[dist2];
                            sumWt += wt;
                            for (int j = 0; j < 3; ++j)
                                bloomImage[3*offset+j] += wt * rgb[3*bloomOffset+j];
                        }
                    }
                bloomImage[3*offset  ] /= sumWt;
                bloomImage[3*offset+1] /= sumWt;
                bloomImage[3*offset+2] /= sumWt;
            }
        }
        // Mix bloom effect into each pixel
        for (int i = 0; i < 3 * nPix; ++i)
            rgb[i] = Lerp(bloomWeight, rgb[i], bloomImage[i]);
        // Free memory allocated for bloom effect
        delete[] bloomFilter;
        delete[] bloomImage;
    }
}

int main(int argc, char *argv[]) 
{
    float scale = 1.f, gamma = 2.2f;
    float bggray = -1.f;
    bool tonemap = false, bloom = false;
    float bloomRadius = .01f, bloomScale = .1f;
    int repeat = 1;
    float rp = 1;

    int argNum = 1;
    while (argNum < argc && argv[argNum][0] == '-') {
#define ARG(name, var) \
      else if (!strcmp(argv[argNum], "-" name)) { \
	if (argNum+1 == argc) \
	  usage(); \
	var = atof(argv[argNum+1]); \
	++argNum; \
      }
        if (!strcmp(argv[argNum], "-tonemap")) tonemap = true;
        else if (!strcmp(argv[argNum], "-bloom")) bloom = true;
        ARG("bloomscale", bloomScale)
        ARG("bloomradius", bloomRadius)
        ARG("repeatpix", rp)
	ARG("scale", scale)
	ARG("gamma", gamma)
	ARG("bg", bggray)
	else 
	    usage();
      ++argNum;

    }
    if (argNum + 2 > argc) usage();
    repeat = int(rp);

    char *inFile = argv[argNum], *outFile = argv[argNum+1]; 	  
    float *rgba;
    int xRes, yRes;
    bool hasAlpha;

    if (ReadEXR(inFile, rgba, xRes, yRes, hasAlpha)) {
        if (repeat > 1) {
            assert(hasAlpha != 0);
            float *rscale = new float[4 * repeat * xRes * repeat * yRes];
            float *rsp = rscale;
            for (int y = 0; y < repeat * yRes; ++y) {
                int yy = y / repeat;
                for (int x = 0; x < repeat * xRes; ++x) {
                    int xx = x / repeat;
                    for (int c = 0; c < 4; ++c) 
                        *rsp++ = rgba[4 * (yy * xRes + xx) + c];
                }
            }
            xRes *= repeat;
            yRes *= repeat;
            rgba = rscale;
        }

	float *rgb = new float[xRes*yRes*3];
	for (int i = 0; i < xRes*yRes; ++i) {
	    for (int j = 0; j < 3; ++j) {
		rgb[3*i+j] = scale * rgba[4*i+j];
		if (bggray > 0)
		    rgb[3*i+j] = rgba[4*i+3] * rgb[3*i+j] + (1.f - rgba[4*i+3]) * bggray;
	    }
	    if (bggray > 0)
		rgba[4*i+3] = 1.f;
	}

        if (bloom)   ::bloom(rgb, xRes, yRes, bloomRadius, bloomScale);
        if (tonemap) reinhard(rgb, xRes, yRes);

	for (int i = 0; i < xRes*yRes; ++i) {
            float m = 0.f;
	    for (int j = 0; j < 3; ++j) {
		rgba[4*i+j] = 255.f * powf(std::max(0.f, rgb[3*i+j]), 1.f / gamma);
                m = std::max(m, rgba[4*i+j]);
            }
            if (m > 255.f) {
                for (int j = 0; j < 3; ++j)
                    rgba[4*i+j] = 255.f * (rgba[4*i+j] / m);
            }
//CO                if (rgba[4*i+j] < 0.f) rgba[4*i+j] = 0.f;
//CO                if (rgba[4*i+j] > 255.f) rgba[4*i+j] = 255.f;

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
