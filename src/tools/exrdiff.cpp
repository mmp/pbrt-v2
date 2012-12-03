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
    fprintf(stderr, "usage: exrdiff [-o difffile.exr] [-d diff tolerance %%] <foo.exr> <bar.exr>\n");
    exit(1);
}

int main(int argc, char *argv[]) 
{
    const char *outfile = NULL;
    const char *imageFile1 = NULL, *imageFile2 = NULL;
    float tol = 0.f;

    if (argc == 1) usage();
    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "-o")) {
            if (!argv[i+1]) usage();
            outfile = argv[i+1];
            ++i;
        }
        else if (!strcmp(argv[i], "-d")) {
            if (!argv[i+1]) usage();
            tol = atof(argv[i+1]);
            ++i;
        }
        else if (!imageFile1)
            imageFile1 = argv[i];
        else if (!imageFile2)
            imageFile2 = argv[i];
        else
            usage();
    }

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
        if ((i % 4) == 3) // alpha channel
            continue;

        sum1 += im1[i];
        sum2 += im2[i];
	float d = fabsf(im1[i] - im2[i]) / im1[i];
        mse += (im1[i] - im2[i]) * (im1[i] - im2[i]);
	if (d > .005) ++smallDiff;
	if (d > .05) ++bigDiff;
    }
    double avg1 = sum1 / (3. * r1[0] * r1[1]);
    double avg2 = sum2 / (3. * r1[0] * r1[1]);
    double avgDelta = (avg1-avg2) / std::min(avg1, avg2);
    if ((tol == 0. && (bigDiff > 0 || smallDiff > 0)) ||
        (tol > 0. && 100.f * fabs(avgDelta) > tol)) {
	printf("%s %s\n\tImages differ: %d big (%.2f%%), %d small (%.2f%%)\n"
               "\tavg 1 = %g, avg2 = %g (%f%% delta)\n"
               "\tMSE = %g, RMS = %.3f%%\n",
               imageFile1, imageFile2,
               bigDiff, 100.f * float(bigDiff) / (3 * r1[0] * r1[1]),
               smallDiff, 100.f * float(smallDiff) / (3 * r1[0] * r1[1]),
               avg1, avg2, 100. * avgDelta,
               mse / (3. * r1[0] * r1[1]),
               100. * sqrt(mse / (3. * r1[0] * r1[1])));
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
