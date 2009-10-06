
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


// cameras/realistic.cpp*
#include "cameras/realistic.h"
#include "film.h"
#include "paramset.h"
#include "sampler.h"
#include "montecarlo.h"
#include "floatfile.h"

// RealisticCamera Local Declarations
struct LensElement {
    // LensElement Public Methods
    LensElement(float r, float z, float e, float ap)
        : radius(r), zpos(z), eta(e), aperture(ap) {
    }
    bool Trace(Ray *ray, float cndr) const;
    bool IsStop() const { return (radius == 0.0f && eta == 0.0f); }

    // LensElement Public Data
    float radius, zpos, eta, aperture;
};


class LensSystem {
public:
    // LensSystem Public Methods
    LensSystem(const char *file, float scale, float fstop, float focaldistance);
    void ComputeFStopLimits() {
        maxFS = 22.0f; // max aperture
        float tstop = minFS = 1.0f;
        while (!SetFStop(tstop) && tstop < maxFS)
            tstop += 0.01f;
        minFS = fstop;
    }
    bool ComputeCardinalPoints(bool disableFront = false);
    bool FocusPrecise(float fplane);
    bool ComputeExitPupil();
    void GenerateRays(Ray *fRay, Ray *bRay);
    LensElement *GetAperture() {
        for (u_int i = 0; i < lenses.size(); i++)
            if (lenses[i].IsStop())
                return &lenses[i];
        return NULL;
    }
    bool SetFStop(float stop);
    bool TraceThick(const Ray &inRay, Ray *outRay) const;
    bool TracePrecise(const Ray &inRay, Ray *outRay, bool reverse = false,
        bool disableAfterStop = false);
    void Reverse();
    float GetImagePlaneZ() const { return iPlane; }
    float GetExitPupilRadius() const { return pupilRadius; }
    float GetExitPupilZ() const { return pupilZ; }

private:
    // LensSystem Private Data
    vector<LensElement> lenses; // Vector of lens interfaces
    float f1, f2;               // Focal points for the lens
    float p1, p2;               // Locations of principle planes
    float pupilRadius;                // Exit pupil for the lens system
    float pupilZ;
    float fstop;                // Current fstop for the camera
    float maxFS;                // Maximum fstop for the camera
    float minFS;                // Minimum fstop for the camera
    float iPlane;                // Image plane
    float oPlane;                // Object focal plane
};



// RealisticCamera Method Definitions
bool LensElement::Trace(Ray *ray, float cndr) const {
    // Compute intersection of ray with spherical lens element
    Point org(ray->o.x, ray->o.y, ray->o.z - (zpos - radius));
    float A = ray->d.x*ray->d.x + ray->d.y*ray->d.y + ray->d.z*ray->d.z;
    float B = 2.f * (ray->d.x*org.x + ray->d.y*org.y + ray->d.z*org.z);
    float C = org.x*org.x + org.y*org.y + org.z*org.z - radius*radius;
    float t0, t1;
    if (!Quadratic(A, B, C, &t0, &t1))
        return false;

    // Determine which intersection should be used for lens element
    if (ray->d.z > 0.f) swap(t0, t1);
    if (radius < 0.f)   swap(t0, t1);
    float thit = t0;

    // Check to see if ray is outside element's aperture
    Point hit = (*ray)(thit);
    float halfAperture = aperture * 0.5f;
    if (hit.x * hit.x + hit.y * hit.y > halfAperture * halfAperture)
        return false;

    // Compute refracted ray at lens interface
    Point lhit = org + thit * ray->d;
    Vector N = Normalize(Vector(lhit));
    Vector I = Normalize(org - lhit);
    if (Dot(I, N) < 0.0f) N = -N;
    float p = Dot(I, N);
    float u = cndr / eta;
    float d = 1.0f - (u * u)*(1.0f - p * p);
    if (d <= 0.0f) return false;
    float g = u * p - sqrtf(d);
    ray->d = g*N - u*I;
    ray->o = hit;
    return true;
}


LensSystem::LensSystem(const char *filename, float scale, float fs,
        float focaldistance)
    : f1(0), f2(0), p1(0), p2(0), fstop(1) {
    oPlane = 100000.0f;
    iPlane = f2;

    // Read lens data into _fileValues_ from file
    vector<float> fileValues;
    if (!ReadFloatFile(filename, &fileValues))
        Error("Unable to open lens data file \"%s\"", filename);
    if ((fileValues.size() % 4) != 0)
        Error("Malformed lens data file \"%s\": %d items found "
              "not a multiple of four", filename, (int)fileValues.size());

    // Loop over the lens element specs and build the lens system
    float dist = 0.0f, zpos = 0.0f;
    for (u_int i = 0; i < fileValues.size(); i += 4) {
        float axpos = dist - zpos;
        dist -= zpos;
    
        float rad = scale * fileValues[i];
        zpos = scale * fileValues[i+1];
        float eta = fileValues[i+2];
        float aper = scale * fileValues[i+3];
        if (aper == 0. && rad == 0. && zpos == 0. && eta == 0.)
            Error("Malformed line in lens file \"%s\"", filename);
        lenses.push_back(LensElement(rad, axpos, eta, aper));
    }

    // Set initial fstop and focal distance for lens system
    Reverse();
    if (!ComputeCardinalPoints())
        Error("Failed to compute cardinal points.");
    ComputeFStopLimits();
    if (!SetFStop(fs))
        Error("Failed to adjust for fstop: %1.1f", fs);
    if (!ComputeExitPupil())
        Error("Failed to compute exit pupil.");
    if (!ComputeCardinalPoints())
        Error("Failed to compute cardinal points.");
    if (!FocusPrecise(focaldistance))
        Error("Failed to focus the camera: %1.1f.", focaldistance);
}


bool LensSystem::ComputeCardinalPoints(bool disableFront) {
    Ray ifray, ibray, ofray, obray;
    GenerateRays(&ifray, &ibray);

    // First trace through from the back side entry
    // Second trace through from the front side entry
    if (!TracePrecise(ibray, &obray, false, disableFront) ||
        !TracePrecise(ifray, &ofray, true, disableFront))
        return false;

    // Compute the focal points and principal planes
    float frp = ofray(1).y, brp = obray(1).y;

    float ftrav = (ifray.o.y - ofray.o.y) / (frp - ofray.o.y);
    float btrav = (ibray.o.y - obray.o.y) / (brp - obray.o.y);
    p2 = ofray(ftrav).z;
    p1 = obray(btrav).z;

    ftrav = (-ofray.o.y) / (frp - ofray.o.y);
    btrav = (-obray.o.y) / (brp - obray.o.y);
    f2 = ofray(ftrav).z;
    f1 = obray(btrav).z;

    return true;
}


bool LensSystem::FocusPrecise(float fplane) {
    if (fplane < 0) return false;
    if (fplane > 100000.0f) {
        iPlane = f2;
        oPlane = INFINITY;
        return true;
    }

    Ray ifray, ibray, ofray, obray;
    GenerateRays(&ifray, &ibray);

    ifray.o.z = fplane;
    ibray = ifray;
    ibray.d = Point(0,0,f1) - ifray.o;

    if (!TracePrecise(ibray, &obray, true) ||
        !TracePrecise(ifray, &ofray, true))
        return false;

    float trav = (ofray.o.y - obray.o.y) / (ofray.o.y - ofray(1).y);
    iPlane = ofray(trav).z;
    oPlane = fplane;
    return true;
}


bool LensSystem::ComputeExitPupil() {
    ComputeCardinalPoints(true);

    LensElement *stop = GetAperture();
    Point rp = Point(0, stop->aperture / 2.0f, stop->zpos - p1);

    // Image the aperture point using a thick lens
    float t = p2 - p1;
    float w = rp.z / (f2 - p2) + 1.0f;
    float zp = rp.z*(1+(t / (f2 - p2))) + t;
    pupilRadius = rp.y / w;
    pupilZ = zp / w + p1;

    return true;
}


void LensSystem::GenerateRays(Ray *fRay, Ray *bRay) {
    float maxAperture = 0.f, minAperture = INFINITY;
    for (u_int i = 0; i < lenses.size(); i++) {
        maxAperture = max(maxAperture, lenses[i].aperture);
        minAperture = min(minAperture, lenses[i].aperture);
    }

    float distZ = maxAperture / 2.0f;
    float fdist = (lenses.end()-1)->zpos + distZ;
    float distY = minAperture / 8.0f;
    float bdist = lenses.begin()->zpos - distZ;

    // First trace through from the left side entry
    *fRay = Ray(Point(0, distY, fdist), Vector(0, 0, -1), 0.f);
    *bRay = Ray(Point(0, -distY, bdist), Vector(0, 0, 1), 0.f);
}


bool LensSystem::SetFStop(float stop) {
    if (stop < minFS) return false;

    LensElement *astop = GetAperture();
    float theta = asin(1.0f / (2*stop));
    float fdiff = lenses[0].zpos - f2;
    float yval = fabs(fdiff*tan(theta));

    Vector direc = Normalize(Vector(0, yval, fdiff));
    Ray fsRay(Point(0, 0, f2), direc, 0.f);
    Ray foRay;

    bool escape = false;
    if (TracePrecise(fsRay, &foRay)) {
// trace again, still forwards, but stop at the stop...
        if (TracePrecise(fsRay, &foRay, false, true)) {
            float movx = fabs((astop->zpos - foRay.o.z) / foRay.d.z);
            float nstop = foRay(movx).y;
            astop->aperture = 2*fabs(nstop);
            escape = true;
        }
    }
    fstop = stop;
    return escape;
}


bool LensSystem::TraceThick(const Ray &inRay, Ray *outRay) const {
    Point ip = Point(inRay.o.x, inRay.o.y, inRay.o.z - p2);
    float dist1 = fabs(inRay.o.z - p1);
    float dist2 = fabs(inRay.o.z - p2);

    float pclose = (dist1 < dist2) ? p1 : p2;
    float pchit  = fabs((inRay.o.z - pclose) / inRay.d.z);
    Point pinter = inRay(pchit);

    float t = p1 - p2;
    Point pext = Point(pinter.x, pinter.y, pinter.z + fabs(t));

    // Image the aperture point using a thick lens
    float w = ip.z / (f1 - p1) + 1.0f;
    float zp = ip.z*(1+(t / (f1 - p1))) + t;
    Point outp = Point(ip.x / w, ip.y / w, zp / w + p2);
    outRay->o = pext;
    outRay->d = Normalize(outp - pext);
    return true;
}


bool LensSystem::TracePrecise(const Ray &inRay, Ray *outRay, bool reverse,
        bool disableAfterStop) {
    if (reverse)
        Reverse();

    Ray tRay = inRay;
    float cindex = 1.0f;
    bool passedStop = false;
    for (u_int i = 0; i < lenses.size(); i++) {
        if (lenses[i].IsStop()) passedStop = true;

        if (cindex == lenses[i].eta) continue;
        if (disableAfterStop && ((!reverse && passedStop) ||
                                 (reverse && !passedStop))) continue;

        if (lenses[i].IsStop()) {
            float aperZ = fabsf((lenses[i].zpos-tRay.o.z) / tRay.d.z);
            Point aperHit = tRay(aperZ);
            float haper = lenses[i].aperture / 2.0f;
            if (aperHit.x*aperHit.x + aperHit.y*aperHit.y > haper*haper) {
                if (reverse) Reverse();
                return false;
            }
        }
        else {
            // Trace lens interface
            if (lenses[i].Trace(&tRay, cindex))
                cindex = lenses[i].eta;
            else {
                if (reverse) Reverse();
                return false;
            }
        }
    }

    *outRay = tRay;
    if (reverse) Reverse();
    return true;
}


void LensSystem::Reverse() {
    reverse(lenses.begin(), lenses.end());
    vector<LensElement>::iterator Li;

    for (Li = lenses.begin(); Li != lenses.end() - 1; Li++) {
        if ((Li + 1)->IsStop()) {
            Li->eta = (Li + 2)->eta;
            Li++;
        }
        else
            Li->eta = (Li + 1)->eta;
    }
    // Lens system ends with an interface to air
    Li->eta = 1.0f;
}


RealisticCamera::~RealisticCamera() {
    delete lensSystem;
}


RealisticCamera::RealisticCamera(const AnimatedTransform &cam2world,
        float sopen, float sclose, Film *film, float focaldistance, float fstop,
        float filmdiagonal, float aspect, float scale, const string &specfile,
        const string &mode, bool constWeight)
    : Camera(cam2world, sopen, sclose, film) {
    lensSystem = new LensSystem(specfile.c_str(), scale, fstop, focaldistance);
    // Compute the width and height of the film
    float fdd = filmdiagonal * filmdiagonal;
    hfilm = sqrtf(fdd / (aspect * aspect + 1.f));
    wfilm = aspect * hfilm;
    constantWeightRays = constWeight;
    // Determine in which mode to trace rays through lens system
    if (mode == "accurate")
        accurate = true;
    else if (mode == "approximate")
        accurate = false;
    else
        Error("Cannot recognize trace mode: \"%s\"", mode.c_str());
}


float RealisticCamera::GenerateRay(const CameraSample &sample,
        Ray *ray) const {
    // Compute the camera ray point on the film plane
    float filmx = .5f - sample.ImageX / film->xResolution;
    float filmy = sample.ImageY / film->yResolution - .5f;
    float imagePlaneZ = lensSystem->GetImagePlaneZ();
    Point org = Point(filmx * wfilm, filmy * hfilm, imagePlaneZ);

    // Compute camera ray direction based on exit pupil
    float exitPupilRadius = lensSystem->GetExitPupilRadius();
    float exitPupilZ = lensSystem->GetExitPupilZ();
    float lensU, lensV;
    ConcentricSampleDisk(sample.LensU, sample.LensV, &lensU, &lensV);
    lensU *= exitPupilRadius;
    lensV *= exitPupilRadius;
    Point epoint(lensU, lensV, exitPupilZ);
    *ray = Ray(org, Normalize(epoint - org), 0.f);

    // Compute form factor _ff_ from point on film to exit pupil
    float ff = 1.f;
    if (!constantWeightRays) {
        Vector dir = Normalize(Point(0, 0, epoint.z) - ray->o);
        float costheta4 = (dir.z * dir.z) * (dir.z * dir.z);
        float Z = (epoint.z - ray->o.z);
        ff = costheta4 * M_PI * exitPupilRadius * exitPupilRadius / (Z * Z);
    }

    // Trace ray through the lens system and compute exiting ray
    Ray oray;
    if (accurate) {
        if (!lensSystem->TracePrecise(*ray, &oray)) return 0.f;
    }
    else {
        if (!lensSystem->TraceThick(*ray, &oray)) return 0.f;
    }
    oray.mint = 0.f;
    oray.maxt = INFINITY;
    oray.time = Lerp(sample.Time, ShutterOpen, ShutterClose);
    CameraToWorld(oray, ray);
    return ff;
}


RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film) {
    float shutteropen = params.FindOneFloat("shutteropen", 0.f);
    float shutterclose = params.FindOneFloat("shutterclose", 1.f);
    float focaldistance = params.FindOneFloat("focaldistance", 1e30f);
    float fstop = params.FindOneFloat("fstop", 0.f);
    float filmdiagonal = params.FindOneFloat("filmdiagonal", 20.f);
    string specfile = params.FindOneString("filename", "empty");
    string mode = params.FindOneString("mode", "accurate");
    float aspect = float(film->xResolution)/float(film->yResolution);
    bool constantWeight = params.FindOneBool("constantweightrays", "true");
    float scale = params.FindOneFloat("scale", 1.f);

    return new RealisticCamera(cam2world, shutteropen, shutterclose, film,
        focaldistance, fstop, filmdiagonal, aspect, scale, specfile, mode,
        constantWeight);
}


