
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


// integrators/dipolesubsurface.cpp*
#include "stdafx.h"
#include "integrators/dipolesubsurface.h"
#include "scene.h"
#include "montecarlo.h"
#include "sampler.h"
#include "progressreporter.h"
#include "intersection.h"
#include "paramset.h"
#include "reflection.h"
#include "octree.h"
#include "camera.h"
#include "floatfile.h"
struct DiffusionReflectance;

// DipoleSubsurfaceIntegrator Local Declarations
struct SubsurfaceOctreeNode {
    // SubsurfaceOctreeNode Methods
    SubsurfaceOctreeNode() {
        isLeaf = true;
        sumArea = 0.f;
        for (int i = 0; i < 8; ++i)
            ips[i] = NULL;
    }
    void Insert(const BBox &nodeBound, IrradiancePoint *ip,
                MemoryArena &arena) {
        Point pMid = .5f * nodeBound.pMin + .5f * nodeBound.pMax;
        if (isLeaf) {
            // Add _IrradiancePoint_ to leaf octree node
            for (int i = 0; i < 8; ++i) {
                if (!ips[i]) {
                    ips[i] = ip;
                    return;
                }
            }

            // Convert leaf node to interior node, redistribute points
            isLeaf = false;
            IrradiancePoint *localIps[8];
            for (int i = 0; i < 8; ++i) {
                localIps[i] = ips[i];
                children[i] = NULL;
            }
            for (int i = 0; i < 8; ++i)  {
                IrradiancePoint *ip = localIps[i];
                // Add _IrradiancePoint_ _ip_ to interior octree node
                int child = (ip->p.x > pMid.x ? 4 : 0) +
                    (ip->p.y > pMid.y ? 2 : 0) + (ip->p.z > pMid.z ? 1 : 0);
                if (!children[child])
                    children[child] = arena.Alloc<SubsurfaceOctreeNode>();
                BBox childBound = octreeChildBound(child, nodeBound, pMid);
                children[child]->Insert(childBound, ip, arena);
            }
            /* fall through to interior case to insert the new point... */
        }
        // Add _IrradiancePoint_ _ip_ to interior octree node
        int child = (ip->p.x > pMid.x ? 4 : 0) +
            (ip->p.y > pMid.y ? 2 : 0) + (ip->p.z > pMid.z ? 1 : 0);
        if (!children[child])
            children[child] = arena.Alloc<SubsurfaceOctreeNode>();
        BBox childBound = octreeChildBound(child, nodeBound, pMid);
        children[child]->Insert(childBound, ip, arena);
    }
    void InitHierarchy() {
        if (isLeaf) {
            // Init _SubsurfaceOctreeNode_ leaf from _IrradiancePoint_s
            float sumWt = 0.f;
            uint32_t i;
            for (i = 0; i < 8; ++i) {
                if (!ips[i]) break;
                float wt = ips[i]->E.y();
                E += ips[i]->E;
                p += wt * ips[i]->p;
                sumWt += wt;
                sumArea += ips[i]->area;
            }
            if (sumWt > 0.f) p /= sumWt;
            E /= i;
        }
        else {
            // Init interior _SubsurfaceOctreeNode_
            float sumWt = 0.f;
            uint32_t nChildren = 0;
            for (uint32_t i = 0; i < 8; ++i) {
                if (!children[i]) continue;
                ++nChildren;
                children[i]->InitHierarchy();
                float wt = children[i]->E.y();
                E += children[i]->E;
                p += wt * children[i]->p;
                sumWt += wt;
                sumArea += children[i]->sumArea;
            }
            if (sumWt > 0.f) p /= sumWt;
            E /= nChildren;
        }
    }
    Spectrum Mo(const BBox &nodeBound, const Point &p, const DiffusionReflectance &Rd,
                float maxError);

    // SubsurfaceOctreeNode Public Data
    Point p;
    bool isLeaf;
    Spectrum E;
    float sumArea;
    union {
        SubsurfaceOctreeNode *children[8];
        IrradiancePoint *ips[8];
    };
};


struct DiffusionReflectance {
    // DiffusionReflectance Public Methods
    DiffusionReflectance(const Spectrum &sigma_a, const Spectrum &sigmap_s,
                         float eta) {
        A = (1.f + Fdr(eta)) / (1.f - Fdr(eta));
        sigmap_t = sigma_a + sigmap_s;
        sigma_tr = Sqrt(3.f * sigma_a * sigmap_t);
        alphap = sigmap_s / sigmap_t;
        zpos = Spectrum(1.f) / sigmap_t;
        zneg = -zpos * (1.f + (4.f/3.f) * A);
    }
    Spectrum operator()(float d2) const {
        Spectrum dpos = Sqrt(Spectrum(d2) + zpos * zpos);
        Spectrum dneg = Sqrt(Spectrum(d2) + zneg * zneg);
        Spectrum Rd = (alphap / (4.f * M_PI)) *
            ((zpos * (dpos * sigma_tr + Spectrum(1.f)) *
              Exp(-sigma_tr * dpos)) / (dpos * dpos * dpos) -
             (zneg * (dneg * sigma_tr + Spectrum(1.f)) *
              Exp(-sigma_tr * dneg)) / (dneg * dneg * dneg));
        return Rd.Clamp();
    }

    // DiffusionReflectance Data
    Spectrum zpos, zneg, sigmap_t, sigma_tr, alphap;
    float A;
};



// DipoleSubsurfaceIntegrator Method Definitions
DipoleSubsurfaceIntegrator::~DipoleSubsurfaceIntegrator() {
    delete[] lightSampleOffsets;
    delete[] bsdfSampleOffsets;
}


void DipoleSubsurfaceIntegrator::RequestSamples(Sampler *sampler, Sample *sample,
        const Scene *scene) {
    // Allocate and request samples for sampling all lights
    uint32_t nLights = scene->lights.size();
    lightSampleOffsets = new LightSampleOffsets[nLights];
    bsdfSampleOffsets = new BSDFSampleOffsets[nLights];
    for (uint32_t i = 0; i < nLights; ++i) {
        const Light *light = scene->lights[i];
        int nSamples = light->nSamples;
        if (sampler) nSamples = sampler->RoundSize(nSamples);
        lightSampleOffsets[i] = LightSampleOffsets(nSamples, sample);
        bsdfSampleOffsets[i] = BSDFSampleOffsets(nSamples, sample);
    }
}


void DipoleSubsurfaceIntegrator::Preprocess(const Scene *scene,
        const Camera *camera, const Renderer *renderer) {
    if (scene->lights.size() == 0) return;
    vector<SurfacePoint> pts;
    // Get _SurfacePoint_s for translucent objects in scene
    if (filename != "") {
        // Initialize _SurfacePoint_s from file
        vector<float> fpts;
        if (ReadFloatFile(filename.c_str(), &fpts)) {
            if ((fpts.size() % 8) != 0)
                Error("Excess values (%d) in points file \"%s\"", int(fpts.size() % 8),
                      filename.c_str());
            for (u_int i = 0; i < fpts.size(); i += 8)
                pts.push_back(SurfacePoint(Point(fpts[i], fpts[i+1], fpts[i+2]),
                                           Normal(fpts[i+3], fpts[i+4], fpts[i+5]),
                                           fpts[i+6], fpts[i+7]));
        }
    }
    if (pts.size() == 0) {
        Point pCamera = camera->CameraToWorld(camera->shutterOpen,
                                              Point(0, 0, 0));
        FindPoissonPointDistribution(pCamera, camera->shutterOpen,
                                     minSampleDist, scene, &pts);
    }

    // Compute irradiance values at sample points
    RNG rng;
    MemoryArena arena;
    PBRT_SUBSURFACE_STARTED_COMPUTING_IRRADIANCE_VALUES();
    ProgressReporter progress(pts.size(), "Computing Irradiances");
    for (uint32_t i = 0; i < pts.size(); ++i) {
        SurfacePoint &sp = pts[i];
        Spectrum E(0.f);
        for (uint32_t j = 0; j < scene->lights.size(); ++j) {
            // Add irradiance from light at point
            const Light *light = scene->lights[j];
            Spectrum Elight = 0.f;
            int nSamples = RoundUpPow2(light->nSamples);
            uint32_t scramble[2] = { rng.RandomUInt(), rng.RandomUInt() };
            uint32_t compScramble = rng.RandomUInt();
            for (int s = 0; s < nSamples; ++s) {
                float lpos[2];
                Sample02(s, scramble, lpos);
                float lcomp = VanDerCorput(s, compScramble);
                LightSample ls(lpos[0], lpos[1], lcomp);
                Vector wi;
                float lightPdf;
                VisibilityTester visibility;
                Spectrum Li = light->Sample_L(sp.p, sp.rayEpsilon,
                    ls, camera->shutterOpen, &wi, &lightPdf, &visibility);
                if (Dot(wi, sp.n) <= 0.) continue;
                if (Li.IsBlack() || lightPdf == 0.f) continue;
                Li *= visibility.Transmittance(scene, renderer, NULL, rng, arena);
                if (visibility.Unoccluded(scene))
                    Elight += Li * AbsDot(wi, sp.n) / lightPdf;
            }
            E += Elight / nSamples;
        }
        irradiancePoints.push_back(IrradiancePoint(sp, E));
        PBRT_SUBSURFACE_COMPUTED_IRRADIANCE_AT_POINT(&sp, &E);
        arena.FreeAll();
        progress.Update();
    }
    progress.Done();
    PBRT_SUBSURFACE_FINISHED_COMPUTING_IRRADIANCE_VALUES();

    // Create octree of clustered irradiance samples
    octree = octreeArena.Alloc<SubsurfaceOctreeNode>();
    for (uint32_t i = 0; i < irradiancePoints.size(); ++i)
        octreeBounds = Union(octreeBounds, irradiancePoints[i].p);
    for (uint32_t i = 0; i < irradiancePoints.size(); ++i)
        octree->Insert(octreeBounds, &irradiancePoints[i], octreeArena);
    octree->InitHierarchy();
}


Spectrum DipoleSubsurfaceIntegrator::Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential &ray, const Intersection &isect,
        const Sample *sample, RNG &rng, MemoryArena &arena) const {
    Spectrum L(0.);
    Vector wo = -ray.d;
    // Compute emitted light if ray hit an area light source
    L += isect.Le(wo);

    // Evaluate BSDF at hit point
    BSDF *bsdf = isect.GetBSDF(ray, arena);
    const Point &p = bsdf->dgShading.p;
    const Normal &n = bsdf->dgShading.nn;
    // Evaluate BSSRDF and possibly compute subsurface scattering
    BSSRDF *bssrdf = isect.GetBSSRDF(ray, arena);
    if (bssrdf && octree) {
        Spectrum sigma_a  = bssrdf->sigma_a();
        Spectrum sigmap_s = bssrdf->sigma_prime_s();
        Spectrum sigmap_t = sigmap_s + sigma_a;
        if (!sigmap_t.IsBlack()) {
            // Use hierarchical integration to evaluate reflection from dipole model
            PBRT_SUBSURFACE_STARTED_OCTREE_LOOKUP(const_cast<Point *>(&p));
            DiffusionReflectance Rd(sigma_a, sigmap_s, bssrdf->eta());
            Spectrum Mo = octree->Mo(octreeBounds, p, Rd, maxError);
            FresnelDielectric fresnel(1.f, bssrdf->eta());
            Spectrum Ft = Spectrum(1.f) - fresnel.Evaluate(AbsDot(wo, n));
            float Fdt = 1.f - Fdr(bssrdf->eta());
            L += (INV_PI * Ft) * (Fdt * Mo);
            PBRT_SUBSURFACE_FINISHED_OCTREE_LOOKUP();
        }
    }
    L += UniformSampleAllLights(scene, renderer, arena, p, n,
        wo, isect.rayEpsilon, ray.time, bsdf, sample, rng, lightSampleOffsets,
        bsdfSampleOffsets);
    if (ray.depth < maxSpecularDepth) {
        // Trace rays for specular reflection and refraction
        L += SpecularReflect(ray, bsdf, rng, isect, renderer, scene, sample,
                             arena);
        L += SpecularTransmit(ray, bsdf, rng, isect, renderer, scene, sample,
                              arena);
    }
    return L;
}


Spectrum SubsurfaceOctreeNode::Mo(const BBox &nodeBound, const Point &pt,
        const DiffusionReflectance &Rd, float maxError) {
    // Compute $M_\roman{o}$ at node if error is low enough
    float dw = sumArea / DistanceSquared(pt, p);
    if (dw < maxError && !nodeBound.Inside(pt))
    {
        PBRT_SUBSURFACE_ADDED_INTERIOR_CONTRIBUTION(const_cast<SubsurfaceOctreeNode *>(this));
        return Rd(DistanceSquared(pt, p)) * E * sumArea;
    }

    // Otherwise compute $M_\roman{o}$ from points in leaf or recursively visit children
    Spectrum Mo = 0.f;
    if (isLeaf) {
        // Accumulate $M_\roman{o}$ from leaf node
        for (int i = 0; i < 8; ++i) {
            if (!ips[i]) break;
            PBRT_SUBSURFACE_ADDED_POINT_CONTRIBUTION(const_cast<IrradiancePoint *>(ips[i]));
            Mo += Rd(DistanceSquared(pt, ips[i]->p)) * ips[i]->E * ips[i]->area;
        }
    }
    else {
        // Recursively visit children nodes to compute $M_\roman{o}$
        Point pMid = .5f * nodeBound.pMin + .5f * nodeBound.pMax;
        for (int child = 0; child < 8; ++child) {
            if (!children[child]) continue;
            BBox childBound = octreeChildBound(child, nodeBound, pMid);
            Mo += children[child]->Mo(childBound, pt, Rd, maxError);
        }
    }
    return Mo;
}


DipoleSubsurfaceIntegrator *CreateDipoleSubsurfaceIntegrator(const ParamSet &params) {
    int maxDepth = params.FindOneInt("maxdepth", 5);
    float maxError = params.FindOneFloat("maxerror", .05f);
    float minDist = params.FindOneFloat("minsampledistance", .25f);
    string pointsfile = params.FindOneString("pointsfile", "");
    if (PbrtOptions.quickRender) { maxError *= 4.f; minDist *= 4.f; }
    return new DipoleSubsurfaceIntegrator(maxDepth, maxError, minDist, pointsfile);
}


