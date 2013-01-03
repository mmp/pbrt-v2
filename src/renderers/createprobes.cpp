
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


// renderers/createprobes.cpp*
#include "stdafx.h"
#include "renderers/createprobes.h"
#include "shapes/sphere.h"
#include "sh.h"
#include "parallel.h"
#include "scene.h"
#include "progressreporter.h"
#include "camera.h"
#include "intersection.h"
#include "integrator.h"
#include "volume.h"
#include "paramset.h"
#include "montecarlo.h"
#if defined(PBRT_IS_WINDOWS) || defined(PBRT_IS_LINUX)
#include <errno.h>
#else
#include <sys/errno.h>
#endif

// CreateRadianceProbes Local Declarations
class CreateRadProbeTask : public Task {
public:
    CreateRadProbeTask(int i, int nProbes[3], float t, const BBox &b, int lmax, bool id,
           bool ii, int nindir, ProgressReporter &p, Sample *sample,
           const vector<Point> &sp, const Scene *sc, const Renderer *ren, Spectrum *c);
    void Run();
private:
    int pointNum, nProbes[3];
    const BBox &bbox;
    int lmax, nIndirSamples;
    float time;
    ProgressReporter &prog;
    bool includeDirectInProbes, includeIndirectInProbes;
    Sample *origSample;
    const Renderer *renderer;
    const Scene *scene;
    const vector<Point> &surfacePoints;
    Spectrum *c_in;
};



// CreateRadianceProbes Method Definitions
CreateRadianceProbes::CreateRadianceProbes(SurfaceIntegrator *surf,
        VolumeIntegrator *vol, const Camera *cam, int lm, float ps, const BBox &b,
        int nindir, bool id, bool ii, float t, const string &fn) {
    lmax = lm;
    probeSpacing = ps;
    bbox = b;
    filename = fn;
    includeDirectInProbes = id;
    includeIndirectInProbes = ii;
    time = t;
    nIndirSamples = nindir;
    surfaceIntegrator = surf;
    volumeIntegrator = vol;
    camera = cam;
}


CreateRadianceProbes::~CreateRadianceProbes() {
    delete surfaceIntegrator;
    delete volumeIntegrator;
}


Spectrum CreateRadianceProbes::Li(const Scene *scene, const RayDifferential &ray,
        const Sample *sample, RNG &rng, MemoryArena &arena, Intersection *isect,
        Spectrum *T) const {
    Assert(ray.time == sample->time);
    Spectrum localT;
    if (!T) T = &localT;
    Intersection localIsect;
    if (!isect) isect = &localIsect;
    Assert(!ray.HasNaNs());
    Spectrum Lo = 0.f;
    if (scene->Intersect(ray, isect))
        Lo = surfaceIntegrator->Li(scene, this, ray, *isect, sample, rng, arena);
    else {
        for (uint32_t i = 0; i < scene->lights.size(); ++i)
           Lo += scene->lights[i]->Le(ray);
    }
    Spectrum Lv = volumeIntegrator->Li(scene, this, ray, sample, rng, T, arena);
    return *T * Lo + Lv;
}


Spectrum CreateRadianceProbes::Transmittance(const Scene *scene,
        const RayDifferential &ray, const Sample *sample, RNG &rng,
        MemoryArena &arena) const {
    return volumeIntegrator->Transmittance(scene, this, ray, sample, rng, arena);
}


void CreateRadianceProbes::Render(const Scene *scene) {
    // Compute scene bounds and initialize probe integrators
    if (bbox.pMin.x > bbox.pMax.x)
        bbox = scene->WorldBound();
    surfaceIntegrator->Preprocess(scene, camera, this);
    volumeIntegrator->Preprocess(scene, camera, this);
    Sample *origSample = new Sample(NULL, surfaceIntegrator, volumeIntegrator,
                                    scene);

    // Compute sampling rate in each dimension
    Vector delta = bbox.pMax - bbox.pMin;
    int nProbes[3];
    for (int i = 0; i < 3; ++i)
        nProbes[i] = max(1, Ceil2Int(delta[i] / probeSpacing));

    // Allocate SH coefficient vector pointers for sample points
    int count = nProbes[0] * nProbes[1] * nProbes[2];
    Spectrum **c_in = new Spectrum *[count];
    for (int i = 0; i < count; ++i)
        c_in[i] = new Spectrum[SHTerms(lmax)];

    // Compute random points on surfaces of scene

    // Create scene bounding sphere to catch rays that leave the scene
    Point sceneCenter;
    float sceneRadius;
    scene->WorldBound().BoundingSphere(&sceneCenter, &sceneRadius);
    Transform ObjectToWorld(Translate(sceneCenter - Point(0,0,0)));
    Transform WorldToObject(Inverse(ObjectToWorld));
    Reference<Shape> sph = new Sphere(&ObjectToWorld, &WorldToObject,
        true, sceneRadius, -sceneRadius, sceneRadius, 360.f);
    Reference<Material> nullMaterial = Reference<Material>(NULL);
    GeometricPrimitive sphere(sph, nullMaterial, NULL);
    vector<Point> surfacePoints;
    uint32_t nPoints = 32768, maxDepth = 32;
    surfacePoints.reserve(nPoints + maxDepth);
    Point pCamera = camera->CameraToWorld(camera->shutterOpen,
                                          Point(0, 0, 0));
    surfacePoints.push_back(pCamera);
    RNG rng;
    while (surfacePoints.size() < nPoints) {
        // Generate random path from camera and deposit surface points
        Point pray = pCamera;
        Vector dir = UniformSampleSphere(rng.RandomFloat(), rng.RandomFloat());
        float rayEpsilon = 0.f;
        for (uint32_t i = 0; i < maxDepth; ++i) {
            Ray ray(pray, dir, rayEpsilon, INFINITY, time);
        
            Intersection isect;
            if (!scene->Intersect(ray, &isect) &&
                !sphere.Intersect(ray, &isect))
                break;
        
            surfacePoints.push_back(ray(ray.maxt));
        
            DifferentialGeometry &hitGeometry = isect.dg;
            pray = isect.dg.p;
            rayEpsilon = isect.rayEpsilon;
            hitGeometry.nn = Faceforward(hitGeometry.nn, -ray.d);
        
            dir = UniformSampleSphere(rng.RandomFloat(), rng.RandomFloat());
            dir = Faceforward(dir, hitGeometry.nn);
        }
    }

    // Launch tasks to compute radiance probes at sample points
    vector<Task *> tasks;
    ProgressReporter prog(count, "Radiance Probes");
    for (int i = 0; i < count; ++i)
        tasks.push_back(new CreateRadProbeTask(i, nProbes, time,
                                   bbox, lmax, includeDirectInProbes,
                                   includeIndirectInProbes, nIndirSamples,
                                   prog, origSample, surfacePoints,
                                   scene, this, c_in[i]));
    EnqueueTasks(tasks);
    WaitForAllTasks();
    for (uint32_t i = 0; i < tasks.size(); ++i)
        delete tasks[i];
    prog.Done();

    // Write radiance probe coefficients to file
    FILE *f = fopen(filename.c_str(), "w");
    if (f) {
        if (fprintf(f, "%d %d %d\n", lmax, includeDirectInProbes?1:0, includeIndirectInProbes?1:0) < 0 ||
            fprintf(f, "%d %d %d\n", nProbes[0], nProbes[1], nProbes[2]) < 0 ||
            fprintf(f, "%f %f %f %f %f %f\n", bbox.pMin.x, bbox.pMin.y, bbox.pMin.z,
                    bbox.pMax.x, bbox.pMax.y, bbox.pMax.z) < 0) {
            Error("Error writing radiance file \"%s\" (%s)", filename.c_str(),
                  strerror(errno));
            exit(1);
        }

        for (int i = 0; i < nProbes[0] * nProbes[1] * nProbes[2]; ++i) {
            for (int j = 0; j < SHTerms(lmax); ++j) {
                fprintf(f, "  ");
                if (c_in[i][j].Write(f) == false) {
                    Error("Error writing radiance file \"%s\" (%s)", filename.c_str(),
                          strerror(errno));
                    exit(1);
                }
                fprintf(f, "\n");
            }
            fprintf(f, "\n");
        }
        fclose(f);
    }
    for (int i = 0; i < nProbes[0] * nProbes[1] * nProbes[2]; ++i)
        delete[] c_in[i];
    delete[] c_in;
    delete origSample;
}


CreateRadProbeTask::CreateRadProbeTask(int pn, int d[3], float t,
        const BBox &b, int lm, bool id,
        bool ii, int nindir, ProgressReporter &p, Sample *samp,
        const vector<Point> &sp, const Scene *sc, const Renderer *ren, Spectrum *c)
        : bbox(b), prog(p), surfacePoints(sp) {
    pointNum = pn;
    lmax = lm;
    nIndirSamples = nindir;
    time = t;
    for (int i = 0; i < 3; ++i) nProbes[i] = d[i];
    origSample = samp;
    c_in = c;
    includeDirectInProbes = id;
    includeIndirectInProbes = ii;
    scene = sc;
    renderer = ren;
}


void CreateRadProbeTask::Run() {
    // Compute region in which to compute incident radiance probes
    int sx = pointNum % nProbes[0];
    int sy = (pointNum / nProbes[0]) % nProbes[1];
    int sz = pointNum / (nProbes[0] * nProbes[1]);
    Assert(sx >= 0 && sx < nProbes[0]);
    Assert(sy >= 0 && sy < nProbes[1]);
    Assert(sz >= 0 && sz < nProbes[2]);
    float tx0 = float(sx) / nProbes[0], tx1 = float(sx+1) / nProbes[0];
    float ty0 = float(sy) / nProbes[1], ty1 = float(sy+1) / nProbes[1];
    float tz0 = float(sz) / nProbes[2], tz1 = float(sz+1) / nProbes[2];
    BBox b(bbox.Lerp(tx0, ty0, tz0), bbox.Lerp(tx1, ty1, tz1));

    // Initialize common variables for _CreateRadProbeTask::Run()_
    RNG rng(pointNum);
    Spectrum *c_probe = new Spectrum[SHTerms(lmax)];
    MemoryArena arena;
    uint32_t nFound = 0, lastVisibleOffset = 0;
    for (int i = 0; i < 256; ++i) {
        if (nFound == 32) break;
        // Try to compute radiance probe contribution at _i_th sample point

        // Compute _i_th candidate point _p_ in cell's bounding box
        float dx = RadicalInverse(i+1, 2);
        float dy = RadicalInverse(i+1, 3);
        float dz = RadicalInverse(i+1, 5);
        Point p = b.Lerp(dx, dy, dz);

        // Skip point _p_ if not indirectly visible from camera
        if (scene->IntersectP(Ray(surfacePoints[lastVisibleOffset],
                                  p - surfacePoints[lastVisibleOffset],
                                  1e-4f, 1.f, time))) {
            uint32_t j;
            // See if point is visible to any element of _surfacePoints_
            for (j = 0; j < surfacePoints.size(); ++j)
                if (!scene->IntersectP(Ray(surfacePoints[j], p - surfacePoints[j],
                                           1e-4f, 1.f, time))) {
                    lastVisibleOffset = j;
                    break;
                }
            if (j == surfacePoints.size())
                continue;
        }
        ++nFound;

        // Compute SH coefficients of incident radiance at point _p_
        if (includeDirectInProbes) {
            for (int i = 0; i < SHTerms(lmax); ++i)
                c_probe[i] = 0.f;
            SHProjectIncidentDirectRadiance(p, 0.f, time, arena, scene,
                                            true, lmax, rng, c_probe);
            for (int i = 0; i < SHTerms(lmax); ++i)
                c_in[i] += c_probe[i];
        }
        
        if (includeIndirectInProbes) {
            for (int i = 0; i < SHTerms(lmax); ++i)
                c_probe[i] = 0.f;
            SHProjectIncidentIndirectRadiance(p, 0.f, time, renderer,
                origSample, scene, lmax, rng, nIndirSamples, c_probe);
            for (int i = 0; i < SHTerms(lmax); ++i)
                c_in[i] += c_probe[i];
        }
        arena.FreeAll();
    }
    // Compute final average value for probe and cleanup
    if (nFound > 0)
        for (int i = 0; i < SHTerms(lmax); ++i)
            c_in[i] /= nFound;
    delete[] c_probe;
        prog.Update();
}


CreateRadianceProbes *CreateRadianceProbesRenderer(const Camera *camera,
        SurfaceIntegrator *surf, VolumeIntegrator *vol,
        const ParamSet &params) {
    bool includeDirect = params.FindOneBool("directlighting", true);
    bool includeIndirect = params.FindOneBool("indirectlighting", true);
    int lmax = params.FindOneInt("lmax", 4);
    int nindir = params.FindOneInt("indirectsamples", 512);
    int nbbox;
    BBox bounds;
    const float *b = params.FindFloat("bounds", &nbbox);
    if (b) {
        if (nbbox != 6) Warning("Expecting six values [x0 y0 z0 x1 y1 z1] for bounds");
        else bounds = BBox(Point(b[0], b[1], b[2]), Point(b[3], b[4], b[5]));
    }
    float probeSpacing = params.FindOneFloat("samplespacing", 1.f);
    float time = params.FindOneFloat("time", 0.f);
    string filename = params.FindOneFilename("filename", "probes.out");

    return new CreateRadianceProbes(surf, vol, camera, lmax, probeSpacing,
        bounds, nindir, includeDirect, includeIndirect, time, filename);
}


