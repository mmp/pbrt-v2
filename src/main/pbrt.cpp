
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


// main/pbrt.cpp*
#include "stdafx.h"
#include "api.h"
#include "probes.h"
#include "parser.h"
#include "parallel.h"

// Tao Du
#include "shapes/paraboloid.h"
#include <iostream>

void testParaboloid(float radius, float zmin, float zmax, float phimax) {
    // Identity transform by default.
    const Transform o2w;
    Paraboloid p(&o2w, &o2w, false, radius, zmin, zmax, phimax);
    // Get the surface area.
    const float analyticSurfaceArea = p.Area();
    std::cout << "Analytic surface area: " << analyticSurfaceArea << std::endl;

    // Numerical results.
    float numericalSurfaceArea = 0.0f;
    const int stepNum = 5000;
    std::vector<Point> points(stepNum + 1);
    for (int i = 0; i < stepNum + 1; ++i) {
        // t's range is [0, 1].
        float t = i * 1.f / stepNum; 
        // z's range is [zmin, zmax].
        float z = zmax * t + zmin * (1.f - t);
        Ray r(Point(0.f, radius * 2, z), Vector(0.f, -1.f, 0.f), 0.f);
        float tHit, rayEpsilon;
        DifferentialGeometry dg;
        p.Intersect(r, &tHit, &rayEpsilon, &dg);
        // Get the intersection point.
        points[i] = dg.p;
    }
    for (int i = 0; i < stepNum; ++i) {
        Point pCurrent = points[i];
        Point pNext = points[i + 1];
        numericalSurfaceArea += Radians(phimax) * pCurrent.y * (pNext - pCurrent).Length();
    }
    std::cout << "Numerical surface area: " << numericalSurfaceArea << std::endl;
}


// main program
int main(int argc, char *argv[]) {
    // Tao Du.
    // radius, zmin, zmax, phimax.
    testParaboloid(3.f, 1.5f, 2.5f, 250.f);
    testParaboloid(155.f, 0.f, 50.f, 172.f);
    testParaboloid(1.f, 1.f, 100.f, 360.f);
    return 0;


    Options options;
    vector<string> filenames;
    // Process command-line arguments
    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "--ncores")) options.nCores = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--outfile")) options.imageFile = argv[++i];
        else if (!strcmp(argv[i], "--quick")) options.quickRender = true;
        else if (!strcmp(argv[i], "--quiet")) options.quiet = true;
        else if (!strcmp(argv[i], "--verbose")) options.verbose = true;
        else if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h")) {
            printf("usage: pbrt [--ncores n] [--outfile filename] [--quick] [--quiet] "
                   "[--verbose] [--help] <filename.pbrt> ...\n");
            return 0;
        }
        else filenames.push_back(argv[i]);
    }

    // Print welcome banner
    if (!options.quiet) {
        printf("pbrt version %s of %s at %s [Detected %d core(s)]\n",
               PBRT_VERSION, __DATE__, __TIME__, NumSystemCores());
        printf("Copyright (c)1998-2014 Matt Pharr and Greg Humphreys.\n");
        printf("The source code to pbrt (but *not* the book contents) is covered by the BSD License.\n");
        printf("See the file LICENSE.txt for the conditions of the license.\n");
        fflush(stdout);
    }
    pbrtInit(options);
    // Process scene description
    PBRT_STARTED_PARSING();
    if (filenames.size() == 0) {
        // Parse scene from standard input
        ParseFile("-");
    } else {
        // Parse scene from input files
        for (u_int i = 0; i < filenames.size(); i++)
            if (!ParseFile(filenames[i]))
                Error("Couldn't open scene file \"%s\"", filenames[i].c_str());
    }
    pbrtCleanup();
    return 0;
}


