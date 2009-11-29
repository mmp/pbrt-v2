/*
*/

#include <stdio.h>
#include <stdlib.h>
#include <ImfInputFile.h>
#include <ImfChannelList.h>
#include <ImfFrameBuffer.h>
#include <ImfRgbaFile.h>
#include <half.h>
#include <assert.h>

using namespace Imf;
using namespace Imath;

static bool ReadEXR(const char *name, float *&rgba, int &xRes, int &yRes, bool &hasAlpha);
static void WriteEXR(const char *name, float *pixels, int xRes, int yRes);

static void usage() {
    fprintf(stderr, "usage: exrdiff [-o difffile.exr] <foo.exr> <bar.exr>\n");
    exit(1);
}

int main(int argc, char *argv[]) 
{
    const char *outfile = NULL;
    const char *imageFile1, *imageFile2;

    if (argc == 1) usage();
    if (!strcmp(argv[1], "-o")) {
        if (argc != 5)
            usage();
        outfile = argv[2];
        imageFile1 = argv[3];
        imageFile2 = argv[4];
    }
    else if (argc == 3) {
        imageFile1 = argv[1];
        imageFile2 = argv[2];
    }
    else
        usage();

    float *im1, *im2;
    int r1[2], r2[2];
    bool hasAlpha;
    if (!ReadEXR(imageFile1, im1, r1[0], r1[1], hasAlpha)) {
	printf("couldn't read image %s\n", imageFile1);
        return 1;
    }
    assert(hasAlpha);
    if (!ReadEXR(imageFile2, im2, r2[0], r2[1], hasAlpha)) {
	printf("couldn't read image %s\n", imageFile2);
	return 1;
    }
    assert(hasAlpha);
    if (r1[0] != r2[0] || r1[1] != r2[1]) {
	printf("%s/%s:\n\tresolutions don't match! (%d,%d) vs (%d,%d)\n",
               imageFile1, imageFile2, r1[0], r1[1], r2[0], r2[1]);
	return 1;
    }

    float *diffImage = NULL;
    if (outfile != NULL)
        diffImage = new float[4 * r1[0] * r1[1]];

    double sum1 = 0.f, sum2 = 0.f;
    int smallDiff = 0, bigDiff = 0;
    double mse = 0.f;
    for (int i = 0; i < 4*r1[0]*r1[1]; ++i) {
        if (diffImage) diffImage[i] = fabsf(im1[i] - im2[i]);
	if (im1[i] == 0 && im2[i] == 0) 
	    continue;

        sum1 += im1[i];
        sum2 += im2[i];
	float d = fabsf(im1[i] - im2[i]) / im1[i];
        mse += (im1[i] - im2[i]) * (im1[i] - im2[i]);
	if (d > .005) ++smallDiff;
	if (d > .05) ++bigDiff;
    }
    double avg1 = sum1 / (r1[0] * r1[1]);
    double avg2 = sum2 / (r1[0] * r1[1]);
    double avgDelta = (avg1-avg2) / std::min(avg1, avg2);
    if (bigDiff > 0 || smallDiff > 0 || fabs(avgDelta) > 1e-4) {
	printf("%s %s\n\tImages differ: %d big (%.2f%%), %d small (%.2f%%)\n"
               "\tavg 1 = %g, avg2 = %g (%f%% delta)\n"
               "\tMSE = %g\n",
               imageFile1, imageFile2,
               bigDiff, 100.f * float(bigDiff) / (4 * r1[0] * r1[1]),
               smallDiff, 100.f * float(smallDiff) / (4 * r1[0] * r1[1]),
               avg1, avg2, 100. * avgDelta,
               mse / (4. * r1[0] * r1[1]));
        if (outfile)
            WriteEXR(outfile, diffImage, r1[0], r1[1]);
        return 1;
    }

    return 0;
}

static bool ReadEXR(const char *name, float *&rgba, int &xRes, int &yRes, bool &hasAlpha)
{
    try {
    InputFile file(name);
    Box2i dw = file.header().dataWindow();
    xRes = dw.max.x - dw.min.x + 1;
    yRes = dw.max.y - dw.min.y + 1;

    half *hrgba = new half[4 * xRes * yRes];

    // for now...
    hasAlpha = true;
    int nChannels = 4;

    hrgba -= 4 * (dw.min.x + dw.min.y * xRes);
    FrameBuffer frameBuffer;
    frameBuffer.insert("R", Slice(HALF, (char *)hrgba,
				  4*sizeof(half), xRes * 4 * sizeof(half), 1, 1, 0.0));
    frameBuffer.insert("G", Slice(HALF, (char *)hrgba+sizeof(half),
				  4*sizeof(half), xRes * 4 * sizeof(half), 1, 1, 0.0));
    frameBuffer.insert("B", Slice(HALF, (char *)hrgba+2*sizeof(half),
				  4*sizeof(half), xRes * 4 * sizeof(half), 1, 1, 0.0));
    frameBuffer.insert("A", Slice(HALF, (char *)hrgba+3*sizeof(half),
				  4*sizeof(half), xRes * 4 * sizeof(half), 1, 1, 1.0));

    file.setFrameBuffer(frameBuffer);
    file.readPixels(dw.min.y, dw.max.y);

    hrgba += 4 * (dw.min.x + dw.min.y * xRes);
    rgba = new float[nChannels * xRes * yRes];
    for (int i = 0; i < nChannels * xRes * yRes; ++i)
	rgba[i] = hrgba[i];
    delete[] hrgba;
    } catch (const std::exception &e) {
        fprintf(stderr, "Unable to read image file \"%s\": %s", name, e.what());
        return NULL;
    }

    return rgba;
}

static void WriteEXR(const char *name, float *pixels, int xRes, int yRes) {
    Rgba *hrgba = new Rgba[xRes * yRes];
    for (int i = 0; i < xRes * yRes; ++i)
        hrgba[i] = Rgba(pixels[4*i], pixels[4*i+1], pixels[4*i+2], 1.);

    Box2i displayWindow(V2i(0,0), V2i(xRes-1, yRes-1));
    Box2i dataWindow = displayWindow;

    RgbaOutputFile file(name, displayWindow, dataWindow, WRITE_RGBA);
    file.setFrameBuffer(hrgba, 1, xRes);
    try {
        file.writePixels(yRes);
    }
    catch (const std::exception &e) {
        fprintf(stderr, "Unable to write image file \"%s\": %s", name,
                e.what());
    }

    delete[] hrgba;
}
