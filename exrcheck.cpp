#include <ImfInputFile.h>
#include <ImfOutputFile.h>
#include <ImfChannelList.h>
#include <ImfFrameBuffer.h>
#include <half.h>

using namespace Imf;
using namespace Imath;

void
ReadImage(const char *name, int *width, int *height) 
{
    InputFile file(name);
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
}


int main() {
    return 0;
}
