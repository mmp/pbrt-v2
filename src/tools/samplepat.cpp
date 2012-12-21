
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


// tools/samplepat.cpp*
#include "pbrt.h"
#include "sampler.h"
#include "progressreporter.h"
#include "montecarlo.h"

// BestCandidate Sampling Constants
#define SQRT_SAMPLE_TABLE_SIZE 64
#define SAMPLE_TABLE_SIZE (SQRT_SAMPLE_TABLE_SIZE * \
                           SQRT_SAMPLE_TABLE_SIZE)

// Sample Pattern Definitions
#define BC_GRID_SIZE 40
typedef vector<int> SampleGrid[BC_GRID_SIZE][BC_GRID_SIZE];
#define GRID(v) (int((v) * BC_GRID_SIZE))

// Sample Pattern Precomputation
int line_num = 0; // make this link!
string current_file; // ditto.

// Pattern Precomputation Local Data
static float imageSamples[SAMPLE_TABLE_SIZE][2];
static float timeSamples[SAMPLE_TABLE_SIZE];
static float lensSamples[SAMPLE_TABLE_SIZE][2];

// Pattern Precomputation Utility Functions
static void addSampleToGrid(float sample[][2], int sampleNum,
                            SampleGrid *grid) {
    int u = GRID(sample[sampleNum][0]);
    int v = GRID(sample[sampleNum][1]);
    (*grid)[u][v].push_back(sampleNum);
}


inline float Wrapped1DDist(float a, float b) {
    float d = fabsf(a - b);
    if (d < 0.5f) return d;
    else return 1.f - max(a, b) + min(a, b);
}



// Pattern Precomputation Forward Declarations
void BestCandidate2D(float table[][2],
                     int count, RNG &rng, SampleGrid *grid = NULL);
static void
Redistribute2D(float samples[][2], SampleGrid &pixelGrid);
int main() {
    RNG rng;
    // Compute image sample positions
    SampleGrid pixelGrid;
    BestCandidate2D(imageSamples, SAMPLE_TABLE_SIZE, rng, &pixelGrid);

    // Compute time samples
    ProgressReporter timeProgress(SAMPLE_TABLE_SIZE, "Time samples");
    for (int i = 0; i < SAMPLE_TABLE_SIZE; ++i)
        timeSamples[i] = (i + rng.RandomFloat()) / SAMPLE_TABLE_SIZE;
    for (int currentSample = 1;
         currentSample < SAMPLE_TABLE_SIZE;
         ++currentSample) {
        // Select best time sample for current image sample
        int best = -1;

        // Find best time relative to neighbors
        float maxMinDelta = 0.;
        for (int t = currentSample; t < SAMPLE_TABLE_SIZE; ++t) {
            // Compute min delta for this time
            int gu = GRID(imageSamples[currentSample][0]);
            int gv = GRID(imageSamples[currentSample][1]);
            float minDelta = INFINITY;
            for (int du = -1; du <= 1; ++du) {
                for (int dv = -1; dv <= 1; ++dv) {
                    // Check offset from times of nearby samples

                    // Compute (u,v) grid cell to check
                    int u = gu + du, v = gv + dv;
                    if (u < 0)             u += BC_GRID_SIZE;
                    if (u >= BC_GRID_SIZE) u -= BC_GRID_SIZE;
                    if (v < 0)             v += BC_GRID_SIZE;
                    if (v >= BC_GRID_SIZE) v -= BC_GRID_SIZE;
                    for (uint32_t g = 0; g < pixelGrid[u][v].size(); ++g) {
                        int otherSample = pixelGrid[u][v][g];
                        if (otherSample < currentSample) {
                            float dt = Wrapped1DDist(timeSamples[otherSample], timeSamples[t]);
                            minDelta = min(minDelta, dt);
                        }
                    }
                }
            }

            // Update _best_ if this is best time so far
            if (minDelta > maxMinDelta) {
                maxMinDelta = minDelta;
                best = t;
            }
        }
        Assert(best != -1);
        swap(timeSamples[best], timeSamples[currentSample]);
        timeProgress.Update();
    }
    timeProgress.Done();;

    // Compute lens samples
    BestCandidate2D(lensSamples, SAMPLE_TABLE_SIZE, rng);
    Redistribute2D(lensSamples, pixelGrid);

    // Write sample table to disk
    FILE *f = fopen("sampledata.out", "w");
    if (f == NULL) {
        Error("Couldn't open sampledata.out for writing. (%s)", strerror(errno));
        return 1;
    }

    fprintf(f, "\n/* Automatically generated %dx%d sample "
            "table (%s @a %s) */\n\n",
            SQRT_SAMPLE_TABLE_SIZE, SQRT_SAMPLE_TABLE_SIZE,
            __DATE__, __TIME__);
    fprintf(f, "const float "
            "BestCandidateSampler::sampleTable[%d][5] = {\n",
            SAMPLE_TABLE_SIZE);
    for (int i = 0; i < SAMPLE_TABLE_SIZE; ++i) {
        fprintf(f, "  { ");
        fprintf(f, "%10.10ff, %10.10ff, ", imageSamples[i][0],
            imageSamples[i][1]);
        fprintf(f, "%10.10ff, ", timeSamples[i]);
        fprintf(f, "%10.10ff, %10.10ff, ", lensSamples[i][0],
            lensSamples[i][1]);
        fprintf(f, "},\n");
    }
    fprintf(f, "};\n");
    return 0;
}


void BestCandidate2D(float table[][2], int totalSamples,
                     RNG &rng, SampleGrid *grid) {
    SampleGrid localGrid;
    if (!grid) grid = &localGrid;
    ProgressReporter
        progress(totalSamples-1, "Throwing Darts");
    // Generate first 2D sample arbitrarily
    table[0][0] = rng.RandomFloat();
    table[0][1] = rng.RandomFloat();
    addSampleToGrid(table, 0, grid);
    for (int currentSample = 1; currentSample < totalSamples;
         ++currentSample) {
        // Generate next best 2D image sample
        float maxDist2 = 0.;
        int numCandidates = 500 * currentSample;
        for (int currentCandidate = 0;
             currentCandidate < numCandidates;
             ++currentCandidate) {
            // Generate a random candidate sample
            float candidate[2];
            candidate[0] = rng.RandomFloat();
            candidate[1] = rng.RandomFloat();

            // Loop over neighboring grid cells and check distances
            float sampleDist2 = INFINITY;
            int gu = GRID(candidate[0]), gv = GRID(candidate[1]);
            for (int du = -1; du <= 1; ++du) {
                for (int dv = -1; dv <= 1; ++dv) {
                    // Compute (u,v) grid cell to check
                    int u = gu + du, v = gv + dv;
                    if (u < 0)             u += BC_GRID_SIZE;
                    if (u >= BC_GRID_SIZE) u -= BC_GRID_SIZE;
                    if (v < 0)             v += BC_GRID_SIZE;
                    if (v >= BC_GRID_SIZE) v -= BC_GRID_SIZE;

                    // Update minimum squared distance from cell's samples
                    for (uint32_t g = 0; g < (*grid)[u][v].size(); ++g) {
                        int s = (*grid)[u][v][g];
                        float xdist = Wrapped1DDist(candidate[0], table[s][0]);
                        float ydist = Wrapped1DDist(candidate[1], table[s][1]);
                        float d2 = xdist*xdist + ydist*ydist;
                        sampleDist2 = min(sampleDist2, d2);
                    }
                }
            }

            // Keep this sample if it is the best one so far
            if (sampleDist2 > maxDist2) {
                maxDist2 = sampleDist2;
                table[currentSample][0] = candidate[0];
                table[currentSample][1] = candidate[1];
            }
        }
        addSampleToGrid(table, currentSample, grid);
        progress.Update();
    }
    progress.Done();
}


static void Redistribute2D(float samples[][2], SampleGrid &pixelGrid) {
    ProgressReporter progress(SAMPLE_TABLE_SIZE, "Redistribution");
    for (int currentSample = 1;
         currentSample < SAMPLE_TABLE_SIZE;
         ++currentSample) {
        // Select best lens sample for current image sample
        int best = -1;

        // Find best 2D sample relative to neighbors
        float maxMinDist2 = 0.f;
        for (int samp = currentSample; samp < SAMPLE_TABLE_SIZE; ++samp) {
            // Check distance to lens positions at nearby samples
            int gu = GRID(imageSamples[currentSample][0]);
            int gv = GRID(imageSamples[currentSample][1]);
            float minDist2 = INFINITY;
            for (int du = -1; du <= 1; ++du) {
                for (int dv = -1; dv <= 1; ++dv) {
                    // Check 2D samples in current grid cell

                    // Compute (u,v) grid cell to check
                    int u = gu + du, v = gv + dv;
                    if (u < 0)             u += BC_GRID_SIZE;
                    if (u >= BC_GRID_SIZE) u -= BC_GRID_SIZE;
                    if (v < 0)             v += BC_GRID_SIZE;
                    if (v >= BC_GRID_SIZE) v -= BC_GRID_SIZE;
                    for (uint32_t g = 0; g < pixelGrid[u][v].size(); ++g) {
                        int s2 = pixelGrid[u][v][g];
                        if (s2 < currentSample) {
                            float dx = Wrapped1DDist(samples[s2][0], samples[samp][0]);
                            float dy = Wrapped1DDist(samples[s2][1], samples[samp][1]);
                            float d2 = dx*dx + dy*dy;
                            minDist2 = min(d2, minDist2);
                        }
                    }
                }
            }

            // Update _best_ for 2D lens sample if it is best so far
            if (minDist2 > maxMinDist2) {
                maxMinDist2 = minDist2;
                best = samp;
            }
        }
        Assert(best != -1);
        swap(samples[best][0], samples[currentSample][0]);
        swap(samples[best][1], samples[currentSample][1]);
        progress.Update();
    }
    fprintf(stderr, "\n");
    progress.Done();
}


