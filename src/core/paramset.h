
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

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CORE_PARAMSET_H
#define PBRT_CORE_PARAMSET_H

// core/paramset.h*
#include "pbrt.h"
#include "fileutil.h"
#include "geometry.h"
#include "texture.h"
#include "spectrum.h"
 #if (_MSC_VER >= 1400)
 #include <stdio.h>
 #define snprintf _snprintf
 #endif
#include <map>
using std::map;

// ParamSet Declarations
class ParamSet {
public:
    // ParamSet Public Methods
    ParamSet() { }
    void AddFloat(const string &, const float *, int nItems = 1);
    void AddInt(const string &, const int *, int nItems);
    void AddBool(const string &, const bool *, int nItems);
    void AddPoint(const string &, const Point *, int nItems);
    void AddVector(const string &, const Vector *, int nItems);
    void AddNormal(const string &, const Normal *, int nItems);
    void AddString(const string &, const string *, int nItems);
    void AddTexture(const string &, const string &);
    void AddRGBSpectrum(const string &, const float *, int nItems);
    void AddXYZSpectrum(const string &, const float *, int nItems);
    void AddBlackbodySpectrum(const string &, const float *, int nItems);
    void AddSampledSpectrumFiles(const string &, const char **, int nItems);
    void AddSampledSpectrum(const string &, const float *, int nItems);
    bool EraseInt(const string &);
    bool EraseBool(const string &);
    bool EraseFloat(const string &);
    bool ErasePoint(const string &);
    bool EraseVector(const string &);
    bool EraseNormal(const string &);
    bool EraseSpectrum(const string &);
    bool EraseString(const string &);
    bool EraseTexture(const string &);
    float FindOneFloat(const string &, float d) const;
    int FindOneInt(const string &, int d) const;
    bool FindOneBool(const string &, bool d) const;
    Point FindOnePoint(const string &, const Point &d) const;
    Vector FindOneVector(const string &, const Vector &d) const;
    Normal FindOneNormal(const string &, const Normal &d) const;
    Spectrum FindOneSpectrum(const string &,
                             const Spectrum &d) const;
    string FindOneString(const string &, const string &d) const;
    string FindOneFilename(const string &, const string &d) const;
    string FindTexture(const string &) const;
    const float *FindFloat(const string &, int *nItems) const;
    const int *FindInt(const string &, int *nItems) const;
    const bool *FindBool(const string &, int *nItems) const;
    const Point *FindPoint(const string &, int *nItems) const;
    const Vector *FindVector(const string &, int *nItems) const;
    const Normal *FindNormal(const string &, int *nItems) const;
    const Spectrum *FindSpectrum(const string &, int *nItems) const;
    const string *FindString(const string &, int *nItems) const;
    void ReportUnused() const;
    void Clear();
    string ToString() const;
    
private:
    // ParamSet Private Data
    vector<Reference<ParamSetItem<bool> > > bools;
    vector<Reference<ParamSetItem<int> > > ints;
    vector<Reference<ParamSetItem<float> > > floats;
    vector<Reference<ParamSetItem<Point> > > points;
    vector<Reference<ParamSetItem<Vector> > > vectors;
    vector<Reference<ParamSetItem<Normal> > > normals;
    vector<Reference<ParamSetItem<Spectrum> > > spectra;
    vector<Reference<ParamSetItem<string> > > strings;
    vector<Reference<ParamSetItem<string> > > textures;
    static map<string, Spectrum> cachedSpectra;
};


template <typename T> struct ParamSetItem : public ReferenceCounted {
    // ParamSetItem Public Methods
    ParamSetItem(const string &name, const T *val, int nItems = 1);
    ~ParamSetItem() {
        delete[] data;
    }

    // ParamSetItem Data
    string name;
    int nItems;
    T *data;
    mutable bool lookedUp;
};



// ParamSetItem Methods
template <typename T>
ParamSetItem<T>::ParamSetItem(const string &n, const T *v, int ni) {
    name = n;
    nItems = ni;
    data = new T[nItems];
    for (int i = 0; i < nItems; ++i) data[i] = v[i];
    lookedUp = false;
}



// TextureParams Declarations
class TextureParams {
public:
    // TextureParams Public Methods
    TextureParams(const ParamSet &geomp, const ParamSet &matp,
                  map<string, Reference<Texture<float> > > &ft,
                  map<string, Reference<Texture<Spectrum> > > &st)
        : floatTextures(ft), spectrumTextures(st),
          geomParams(geomp), materialParams(matp) {
    }
    Reference<Texture<Spectrum> > GetSpectrumTexture(const string &name,
            const Spectrum &def) const;
    Reference<Texture<float> > GetFloatTexture(const string &name,
            float def) const;
    Reference<Texture<float> > GetFloatTextureOrNull(const string &name) const;
    float FindFloat(const string &n, float d) const {
        return geomParams.FindOneFloat(n, materialParams.FindOneFloat(n, d));
    }
    string FindString(const string &n, const string &d = "") const {
           return geomParams.FindOneString(n, materialParams.FindOneString(n, d));
    }
    string FindFilename(const string &n, const string &d = "") const {
           return geomParams.FindOneFilename(n, materialParams.FindOneFilename(n, d));
    }
    int FindInt(const string &n, int d) const {
           return geomParams.FindOneInt(n, materialParams.FindOneInt(n, d));
    }
    bool FindBool(const string &n, bool d) const {
           return geomParams.FindOneBool(n, materialParams.FindOneBool(n, d));
    }
    Point FindPoint(const string &n, const Point &d) const {
           return geomParams.FindOnePoint(n, materialParams.FindOnePoint(n, d));
    }
    Vector FindVector(const string &n, const Vector &d) const {
           return geomParams.FindOneVector(n, materialParams.FindOneVector(n, d));
    }
    Normal FindNormal(const string &n, const Normal &d) const {
           return geomParams.FindOneNormal(n, materialParams.FindOneNormal(n, d));
    }
    Spectrum FindSpectrum(const string &n, const Spectrum &d) const {
           return geomParams.FindOneSpectrum(n, materialParams.FindOneSpectrum(n, d));
    }
    void ReportUnused() const {
        geomParams.ReportUnused();
        materialParams.ReportUnused();
    }
    const ParamSet &GetGeomParams() const { return geomParams; }
    const ParamSet &GetMaterialParams() const { return materialParams; }
private:
    // TextureParams Private Data
    map<string, Reference<Texture<float> > > &floatTextures;
    map<string, Reference<Texture<Spectrum> > > &spectrumTextures;
    const ParamSet &geomParams, &materialParams;
};



#endif // PBRT_CORE_PARAMSET_H
