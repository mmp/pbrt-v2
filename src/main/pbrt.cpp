
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

#include <iostream>
#include "shapes/disk.h"


// main program
int main(int argc, char *argv[]) {
    // Tao Du.
    // Test Disk derivatives.
    const float height = 1.f;
    const float radius = 2.f;
    const float innerRadius = 1.f;
    const float phiMax = 180.f;
    const Transform object2World;
    Disk *disk = new Disk(&object2World, &object2World, false, height, radius, innerRadius, phiMax);

    // Test dpdu and dpdv.
    const Vector dir(0.f, 0.f, 1.f);
    const Point origin(1.5f, 0.f, 0.f);
    const float eps = 1e-3;
    Ray r1(origin, dir, 0.f);
    Ray r2 = r1;
    float tHit, rayEpsilon;
    DifferentialGeometry dg1, dg2;

    disk->Intersect(r1, &tHit, &rayEpsilon, &dg1);
    // Print analytical dpdu and dpdv from DifferentialGeometry.
    std::cout << "analytical dpdu: (" << dg1.dpdu.x << ", " << dg1.dpdu.y << ", " << dg1.dpdu.z << ")" << std::endl;
    std::cout << "analytical dpdv: (" << dg1.dpdv.x << ", " << dg1.dpdv.y << ", " << dg1.dpdv.z << ")" << std::endl;

    /////////////////////////////////////////
    // Test dpdu only.
    /////////////////////////////////////////
    r2 = r1;
    r2.o += Vector(0.f, eps, 0.f);
    // Intersection.
    disk->Intersect(r2, &tHit, &rayEpsilon, &dg2);
    // Print du and dv.
    float du = dg2.u - dg1.u;
    float dv = dg2.v - dg1.v;
    std::cout << "du: " << du << " dv: " << dv << std::endl;
    // Compute numerical dpdu.
    Vector ndpdu = (dg2.p - dg1.p) / du;
    std::cout << "numerical dpdu: (" << ndpdu.x << ", " << ndpdu.y << ", " << ndpdu.z << ")" << std::endl;

    /////////////////////////////////////////
    // Test dpdv only.
    /////////////////////////////////////////
    r2 = r1;
    r2.o += Vector(eps, 0.f, 0.f);
    // Intersection.
    disk->Intersect(r2, &tHit, &rayEpsilon, &dg2);
    // Print du and dv.
    du = dg2.u - dg1.u;
    dv = dg2.v - dg1.v;
    std::cout << "du: " << du << " dv: " << dv << std::endl;
    // Compute numerical dpdv.
    Vector ndpdv = (dg2.p - dg1.p) / dv;
    std::cout << "numerical dpdv: (" << ndpdv.x << ", " << ndpdv.y << ", " << ndpdv.z << ")" << std::endl;

    /////////////////////////////////////////
    // Test both dpdu and dpdv.
    /////////////////////////////////////////
    r2 = r1;
    r2.o += Vector(eps, eps, 0.f);
    // Intersection.
    disk->Intersect(r2, &tHit, &rayEpsilon, &dg2);
    // Compare p2 and p1 + dpdu * du + dpdv * dv.
    std::cout << "p2: (" << dg2.p.x << ", " << dg2.p.y << ", " << dg2.p.z << ")" << std::endl;
    du = dg2.u - dg1.u;
    dv = dg2.v - dg1.v;
    Point p2 = dg1.p + dg1.dpdu * du + dg1.dpdv * dv;
    std::cout << "p1 + dpdu * du + dpdv * dv: (" << p2.x << ", " << p2.y << ", " << p2.z << ")" << std::endl; 
    // The higher order term.
    std::cout << "relative error: " << (dg2.p - p2).Length() / (dg2.p - dg1.p).Length() * 100.f << "%" << std::endl;
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


