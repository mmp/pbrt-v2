
/*
    pbrt source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.

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


// materials/measured.cpp*
#include "stdafx.h"
#include "materials/measured.h"
#include "paramset.h"
#include "floatfile.h"

/*
  File format descriptions:

  -- Irregularly Sampled Isotropic BRDF --

  This is the format of the BRDFs in the scenes/brdfs/ folder of the pbrt
  distribution.  This is a simple text format of numbers in a particular
  format; the hash character # is used to denote a comment that continues
  to the end of the current line.

  The first number in the file gives the number of wavelengths at which the
  reflection data was measured, numWls.  This is followed by numWls values
  that give the frequency in nm of these wavelengths.  Each BRDF
  measurement is represented by 4+numWls values.  The first two give the
  (theta,phi) angles of the incident illumination direction, the next two
  give (theta,phi) for the measured reflection direction, and the following
  numWls give the spectral coefficients for the measurement at the
  wavelength specified at the start of the file.


  -- Regular Half-Angle BRDF --
  This is the file format used in the MERL BRDF database; see http://merl.com/brdf.

  This file format is a binary format, with numbers encoded in low-endian
  form.  It represents a regular 3D tabularization of BRDF samples in RGB
  color where the dimensions indexed over are (delta phi, delta theta,
  sqrt(theta_h)).  Here, theta_h is the angle between the halfangle vector
  and the normal, and delta theta and delta phi are the offset in theta and
  phi of one of the two directions.  (Note that the offset would be the
  same for the other direction, since it's from the half-angle vector.)

  The starts with three 32-bit integers, giving the resolution of the
  overall table.  It then containes a number of samples equal to the
  product of those three integers, times 3 (for RGB).  Samples are laid out
  with delta phi the minor index, then delta theta, then sqrt(theta_h) as
  the major index.

  In the file each sample should be scaled by RGB(1500,1500,1500/1.6) of
  the original measurement.  (In order words, the sample values are scaled
  by the inverse of that as they are read in.  
*/


// MeasuredMaterial Method Definitions
static map<string, float *> loadedRegularHalfangle;
static map<string, KdTree<IrregIsotropicBRDFSample> *> loadedThetaPhi;
MeasuredMaterial::MeasuredMaterial(const string &filename,
      Reference<Texture<float> > bump) {
    bumpMap = bump;
    const char *suffix = strrchr(filename.c_str(), '.');
    regularHalfangleData = NULL;
    thetaPhiData = NULL;
    if (!suffix)
        Error("No suffix in measured BRDF filename \"%s\".  "
              "Can't determine file type (.brdf / .merl)", filename.c_str());
    if (!strcmp(suffix, ".brdf") || !strcmp(suffix, ".BRDF")) {
        // Load $(\theta, \phi)$ measured BRDF data
        if (loadedThetaPhi.find(filename) != loadedThetaPhi.end()) {
            thetaPhiData = loadedThetaPhi[filename];
            return;
        }
        
        vector<float> values;
        if (!ReadFloatFile(filename.c_str(), &values))
            Error("Unable to read BRDF data from file \"%s\"", filename.c_str());
        
        uint32_t pos = 0;
        int numWls = int(values[pos++]);
        if ((values.size() - 1 - numWls) % (4 + numWls) != 0)
            Error("Excess or insufficient data in theta, phi BRDF file \"%s\"",
                  filename.c_str());
        vector<float> wls;
        for (int i = 0; i < numWls; ++i)
            wls.push_back(values[pos++]);
        
        BBox bbox;
        vector<IrregIsotropicBRDFSample> samples;
        while (pos < values.size()) {
            float thetai = values[pos++];
            float phii = values[pos++];
            float thetao = values[pos++];
            float phio = values[pos++];
            Vector wo = SphericalDirection(sinf(thetao), cosf(thetao), phio);
            Vector wi = SphericalDirection(sinf(thetai), cosf(thetai), phii);
            Spectrum s = Spectrum::FromSampled(&wls[0], &values[pos], numWls);
            pos += numWls;
            Point p = BRDFRemap(wo, wi);
            samples.push_back(IrregIsotropicBRDFSample(p, s));
            bbox = Union(bbox, p);
        }
        loadedThetaPhi[filename] = thetaPhiData = new KdTree<IrregIsotropicBRDFSample>(samples);
    }
    else {
        // Load RegularHalfangle BRDF Data
        nThetaH = 90;
        nThetaD = 90;
        nPhiD = 180;
        
        if (loadedRegularHalfangle.find(filename) != loadedRegularHalfangle.end()) {
            regularHalfangleData = loadedRegularHalfangle[filename];
            return;
        }
        
        FILE *f = fopen(filename.c_str(), "rb");
        if (!f) {
            Error("Unable to open BRDF data file \"%s\"", filename.c_str());
            return;
        }
        int dims[3];
        if (fread(dims, sizeof(int), 3, f) != 3) {
            Error("Premature end-of-file in measured BRDF data file \"%s\"",
                  filename.c_str());
            fclose(f);
            return;
        }
        uint32_t n = dims[0] * dims[1] * dims[2];
        if (n != nThetaH * nThetaD * nPhiD)  {
            Error("Dimensions don't match\n");
            fclose(f);
            return;
        }
        
        regularHalfangleData = new float[3*n];
        const uint32_t chunkSize = 2*nPhiD;
        double *tmp = ALLOCA(double, chunkSize);
        uint32_t nChunks = n / chunkSize;
        Assert((n % chunkSize) == 0);
        float scales[3] = { 1.f/1500.f, 1.15f/1500.f, 1.66f/1500.f };
        for (int c = 0; c < 3; ++c) {
            int offset = 0;
            for (uint32_t i = 0; i < nChunks; ++i) {
                if (fread(tmp, sizeof(double), chunkSize, f) != chunkSize) {
                    Error("Premature end-of-file in measured BRDF data file \"%s\"",
                          filename.c_str());
                    delete[] regularHalfangleData;
                    regularHalfangleData = NULL;
                    fclose(f);
                    return;
                }
                for (uint32_t j = 0; j < chunkSize; ++j)
                    regularHalfangleData[3 * offset++ + c] = max(0., tmp[j] * scales[c]);
            }
        }
        
        loadedRegularHalfangle[filename] = regularHalfangleData;
        fclose(f);
    }
}


BSDF *MeasuredMaterial::GetBSDF(const DifferentialGeometry &dgGeom,
                                const DifferentialGeometry &dgShading,
                                MemoryArena &arena) const {
    // Allocate _BSDF_, possibly doing bump mapping with _bumpMap_
    DifferentialGeometry dgs;
    if (bumpMap)
        Bump(bumpMap, dgGeom, dgShading, &dgs);
    else
        dgs = dgShading;
    BSDF *bsdf = BSDF_ALLOC(arena, BSDF)(dgs, dgGeom.nn);
    if (regularHalfangleData)
        bsdf->Add(BSDF_ALLOC(arena, RegularHalfangleBRDF)
            (regularHalfangleData, nThetaH, nThetaD, nPhiD));
    else if (thetaPhiData)
        bsdf->Add(BSDF_ALLOC(arena, IrregIsotropicBRDF)(thetaPhiData));
    return bsdf;
}


MeasuredMaterial *CreateMeasuredMaterial(const Transform &xform,
        const TextureParams &mp) {
    Reference<Texture<float> > bumpMap = mp.GetFloatTexture("bumpmap", 0.f);
    return new MeasuredMaterial(mp.FindFilename("filename"), bumpMap);
}


