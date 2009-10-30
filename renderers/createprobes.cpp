
/*
    pbrt source code Copyright(c) 1998-2009 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    pbrt is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.  Note that the text contents of
    the book "Physically Based Rendering" are *not* licensed under the
    GNU GPL.

    pbrt is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */


// renderers/createprobes.cpp*
#include "renderers/createprobes.h"
#include "integrators/photonmap.h"
#include "integrators/directlighting.h"
#include "integrators/emission.h"
#include "shapes/sphere.h"
#include "sh.h"
#include "parallel.h"
#include "scene.h"
#include "progressreporter.h"
#include "camera.h"
#include "intersection.h"
#include "paramset.h"
#include "montecarlo.h"

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
CreateRadianceProbes::CreateRadianceProbes(int lm, float ps, const BBox &b, int md,
        int nindir, bool id, bool ii, float t, const Point &pcam, const string &fn) {
    lmax = lm;
    probeSpacing = ps;
    bbox = b;
    filename = fn;
    includeDirectInProbes = id;
    includeIndirectInProbes = ii;
    time = t;
    nIndirSamples = nindir;
    // Create surface and volume integrators for _CreateRadianceProbes_
    if (includeIndirectInProbes) {
        int ncaus = 0, nindir = 400000;
        int nLookup = 100, maxphotondepth = 2;
        float maxdist = .5f;
        bool finalGather = false;
        int gatherSamples = 256;
        float ga = 10.f;
        if (getenv("PBRT_QUICK_RENDER")) {
            nindir = nindir / 10;
            nLookup = max(1, nLookup / 10);
            gatherSamples = max(1, gatherSamples / 4);
        }
        surfaceIntegrator = new PhotonIntegrator(ncaus, nindir, nLookup, md,
            maxphotondepth, maxdist, finalGather, gatherSamples, ga);
        volumeIntegrator = new EmissionIntegrator(1.f);
    }
    else {
        if (!includeDirectInProbes)
            Warning("Radiance probes will be including neither direct nor "
                    "indirect light?!");
        surfaceIntegrator = NULL;
        volumeIntegrator = NULL;
    }
    origSample = NULL;
    nProbes[0] = nProbes[1] = nProbes[2] = -1;
    c_in = NULL;
    pCamera = pcam;
}


CreateRadianceProbes::~CreateRadianceProbes() {
    delete surfaceIntegrator;
    delete volumeIntegrator;
    delete origSample;
    for (int i = 0; i < nProbes[0] * nProbes[1] * nProbes[2]; ++i)
        delete[] c_in[i];
    delete[] c_in;
}


Spectrum CreateRadianceProbes::Li(const Scene *scene, const RayDifferential &ray,
    const Sample *sample, MemoryArena &arena, Intersection *isect,
    Spectrum *T) const {
    Assert(ray.time == sample->time);
    Spectrum localT;
    if (!T) T = &localT;
    Intersection localIsect;
    if (!isect) isect = &localIsect;
    Assert(!ray.HasNaNs());
    Spectrum Lo = 0.f;
    if (scene->Intersect(ray, isect))
        Lo = surfaceIntegrator->Li(scene, this, ray, *isect, sample, arena);
    else {
        for (u_int i = 0; i < scene->lights.size(); ++i)
           Lo += scene->lights[i]->Le(ray);
    }
    Spectrum Lv = volumeIntegrator->Li(scene, this, ray, sample, T, arena);
    return *T * Lo + Lv;
}


Spectrum CreateRadianceProbes::Transmittance(const Scene *scene,
    const RayDifferential &ray, const Sample *sample,
    MemoryArena &arena, RNG *rng) const {
    return volumeIntegrator->Transmittance(scene, this, ray, sample, rng, arena);
}


void CreateRadianceProbes::Render(const Scene *scene) {
    // Compute scene bounds and initialize probe integrators
    if (bbox.pMin.x > bbox.pMax.x)
        bbox = scene->WorldBound();
    if (surfaceIntegrator && volumeIntegrator) {
        // This will be trouble if one of the integrators needs a Camera
        // pointer for its preprocess method...  Currently that's only the
        // irradiance cache.  Apologies if this bites you and you find this
        // comment after debugging. :-p
        surfaceIntegrator->Preprocess(scene, NULL, this);
        volumeIntegrator->Preprocess(scene, NULL, this);
        origSample = new Sample(NULL, surfaceIntegrator, volumeIntegrator,
                                scene);
    }

    // Compute sampling rate in each dimension
    Vector delta = bbox.pMax - bbox.pMin;
    for (int i = 0; i < 3; ++i)
        nProbes[i] = max(1, Ceil2Int(delta[i] / probeSpacing));

    // Allocate SH coefficient vector pointers for sample points
    int count = nProbes[0] * nProbes[1] * nProbes[2];
    c_in = new Spectrum *[count];
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
    u_int nPoints = 32768, maxDepth = 32;
    surfacePoints.reserve(nPoints + maxDepth);
    surfacePoints.push_back(pCamera);
    RNG rng;
    while (surfacePoints.size() < nPoints) {
        // Generate random path from eye and deposit surface points
        Point pray = pCamera;
        Vector dir = UniformSampleSphere(rng.RandomFloat(), rng.RandomFloat());
        float rayEpsilon = 0.f;
        for (u_int i = 0; i < maxDepth; ++i) {
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
    for (u_int i = 0; i < tasks.size(); ++i)
        delete tasks[i];
    prog.Done();

    // Write radiance probe coefficients to file
    FILE *f = fopen(filename.c_str(), "w");
    if (f) {
        if (fprintf(f, "%d %d %d\n", lmax, includeDirectInProbes?1:0, includeIndirectInProbes?1:0) < 0 ||
            fprintf(f, "%d %d %d\n", nProbes[0], nProbes[1], nProbes[2]) < 0 ||
            fprintf(f, "%f %f %f %f %f %f\n", bbox.pMin.x, bbox.pMin.y, bbox.pMin.z,
              bbox.pMax.x, bbox.pMax.y, bbox.pMax.z) < 0)
            Severe("Error writing radiance file \"%s\"", filename.c_str());
        for (int i = 0; i < nProbes[0] * nProbes[1] * nProbes[2]; ++i) {
            for (int j = 0; j < SHTerms(lmax); ++j) {
                fprintf(f, "  ");
                if (c_in[i][j].Write(f) == false)
                    Severe("Error writing radiance file \"%s\"", filename.c_str());
                fprintf(f, "\n");
            }
            fprintf(f, "\n");
        }
        fclose(f);
    }
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
    u_int nFound = 0, lastVisibleOffset = 0;
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
            u_int j;
            for (j = 0; j < surfacePoints.size(); ++j) {
                if (!scene->IntersectP(Ray(surfacePoints[j], p - surfacePoints[j],
                                           1e-4f, 1.f, time))) {
                    lastVisibleOffset = j;
                    break;
                }
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


CreateRadianceProbes *CreateRadianceProbesRenderer(const Point &pCamera,
        const ParamSet &params) {
    bool includeDirect = params.FindOneBool("directlighting", true);
    bool includeIndirect = params.FindOneBool("indirectlighting", true);
    int lmax = params.FindOneInt("lmax", 4);
    int nindir = params.FindOneInt("indirectsamples", 512);
    int maxDepth = params.FindOneInt("maxdepth", 5);
    int nbbox;
    BBox bounds;
    const float *b = params.FindFloat("bounds", &nbbox);
    if (b) {
        if (nbbox != 6) Warning("Expecting six values [x0 y0 z0 x1 y1 z1] for bounds");
        else bounds = BBox(Point(b[0], b[1], b[2]), Point(b[3], b[4], b[5]));
    }
    float probeSpacing = params.FindOneFloat("samplespacing", 1.f);
    float time = params.FindOneFloat("time", 0.f);
    string filename = params.FindOneString("filename", "probes.out");

    return new CreateRadianceProbes(lmax, probeSpacing, bounds, maxDepth,
        nindir, includeDirect, includeIndirect, time, pCamera, filename);
}


