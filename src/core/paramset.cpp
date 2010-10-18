
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


// core/paramset.cpp*
#include "stdafx.h"
#include "paramset.h"
#include "floatfile.h"
#include "textures/constant.h"

// ParamSet Macros
#define ADD_PARAM_TYPE(T, vec) \
    (vec).push_back(new ParamSetItem<T>(name, (T *)data, nItems))
#define LOOKUP_PTR(vec) \
    for (uint32_t i = 0; i < (vec).size(); ++i) \
        if ((vec)[i]->name == name) { \
            *nItems = (vec)[i]->nItems; \
            (vec)[i]->lookedUp = true; \
            return (vec)[i]->data; \
        } \
    return NULL
#define LOOKUP_ONE(vec) \
    for (uint32_t i = 0; i < (vec).size(); ++i) { \
        if ((vec)[i]->name == name && \
            (vec)[i]->nItems == 1) { \
            (vec)[i]->lookedUp = true; \
            return *((vec)[i]->data); \
}        } \
    return d

// ParamSet Methods
void ParamSet::AddFloat(const string &name, const float *data,
                        int nItems) {
    EraseFloat(name);
    floats.push_back(new ParamSetItem<float>(name, data, nItems));
}


void ParamSet::AddInt(const string &name, const int *data, int nItems) {
    EraseInt(name);
    ADD_PARAM_TYPE(int, ints);
}


void ParamSet::AddBool(const string &name, const bool *data, int nItems) {
    EraseBool(name);
    ADD_PARAM_TYPE(bool, bools);
}


void ParamSet::AddPoint(const string &name, const Point *data, int nItems) {
    ErasePoint(name);
    ADD_PARAM_TYPE(Point, points);
}


void ParamSet::AddVector(const string &name, const Vector *data, int nItems) {
    EraseVector(name);
    ADD_PARAM_TYPE(Vector, vectors);
}


void ParamSet::AddNormal(const string &name, const Normal *data, int nItems) {
    EraseNormal(name);
    ADD_PARAM_TYPE(Normal, normals);
}


void ParamSet::AddRGBSpectrum(const string &name, const float *data, int nItems) {
    EraseSpectrum(name);
    Assert(nItems % 3 == 0);
    nItems /= 3;
    Spectrum *s = new Spectrum[nItems];
    for (int i = 0; i < nItems; ++i)
        s[i] = Spectrum::FromRGB(&data[3*i]);
    spectra.push_back(new ParamSetItem<Spectrum>(name, s, nItems));
    delete[] s;
}


void ParamSet::AddXYZSpectrum(const string &name, const float *data, int nItems) {
    EraseSpectrum(name);
    Assert(nItems % 3 == 0);
    nItems /= 3;
    Spectrum *s = new Spectrum[nItems];
    for (int i = 0; i < nItems; ++i)
        s[i] = Spectrum::FromXYZ(&data[3*i]);
    spectra.push_back(new ParamSetItem<Spectrum>(name, s, nItems));
    delete[] s;
}


void ParamSet::AddBlackbodySpectrum(const string &name, const float *data,
        int nItems) {
    EraseSpectrum(name);
    Assert(nItems % 2 == 0); // temperature (K), scale, ...
    nItems /= 2;
    Spectrum *s = new Spectrum[nItems];
    float *v = new float[nCIESamples];
    for (int i = 0; i < nItems; ++i) {
        Blackbody(CIE_lambda, nCIESamples, data[2*i], v);
        s[i] = data[2*i+1] * Spectrum::FromSampled(CIE_lambda, v, nCIESamples);
    }
    spectra.push_back(new ParamSetItem<Spectrum>(name, s, nItems));
    delete[] s;
    delete[] v;
}


void ParamSet::AddSampledSpectrum(const string &name, const float *data,
        int nItems) {
    EraseSpectrum(name);
    Assert(nItems % 2 == 0);
    nItems /= 2;
    float *wl = new float[nItems], *v = new float[nItems];
    for (int i = 0; i < nItems; ++i) {
        wl[i] = data[2*i];
        v[i] = data[2*i+1];
    }
    Spectrum s = Spectrum::FromSampled(wl, v, nItems);
    spectra.push_back(new ParamSetItem<Spectrum>(name, &s, 1));
}


void ParamSet::AddSampledSpectrumFiles(const string &name, const char **names,
        int nItems) {
    EraseSpectrum(name);
    Spectrum *s = new Spectrum[nItems];
    for (int i = 0; i < nItems; ++i) {
        string fn = AbsolutePath(ResolveFilename(names[i]));
        if (cachedSpectra.find(fn) != cachedSpectra.end()) {
            s[i] = cachedSpectra[fn];
            continue;
        }

        vector<float> vals;
        if (!ReadFloatFile(fn.c_str(), &vals)) {
            Warning("Unable to read SPD file \"%s\".  Using black distribution.", 
                    fn.c_str());
            s[i] = Spectrum(0.);
        }
        else {
            if (vals.size() % 2) {
                Warning("Extra value found in spectrum file \"%s\". "
                        "Ignoring it.", fn.c_str());
            }
            vector<float> wls, v;
            for (uint32_t j = 0; j < vals.size() / 2; ++j) {
                wls.push_back(vals[2*j]);
                v.push_back(vals[2*j+1]);
            }
            s[i] = Spectrum::FromSampled(&wls[0], &v[0], wls.size());
        }
        cachedSpectra[fn] = s[i];
    }

    spectra.push_back(new ParamSetItem<Spectrum>(name, s, nItems));
    delete[] s;
}


map<string, Spectrum> ParamSet::cachedSpectra;
void ParamSet::AddString(const string &name, const string *data, int nItems) {
    EraseString(name);
    ADD_PARAM_TYPE(string, strings);
}


void ParamSet::AddTexture(const string &name, const string &value) {
    EraseTexture(name);
    textures.push_back(new ParamSetItem<string>(name, (string *)&value, 1));
}


bool ParamSet::EraseInt(const string &n) {
    for (uint32_t i = 0; i < ints.size(); ++i)
        if (ints[i]->name == n) {
            ints.erase(ints.begin() + i);
            return true;
        }
    return false;
}


bool ParamSet::EraseBool(const string &n) {
    for (uint32_t i = 0; i < bools.size(); ++i)
        if (bools[i]->name == n) {
            bools.erase(bools.begin() + i);
            return true;
        }
    return false;
}


bool ParamSet::EraseFloat(const string &n) {
    for (uint32_t i = 0; i < floats.size(); ++i)
        if (floats[i]->name == n) {
            floats.erase(floats.begin() + i);
            return true;
        }
    return false;
}


bool ParamSet::ErasePoint(const string &n) {
    for (uint32_t i = 0; i < points.size(); ++i)
        if (points[i]->name == n) {
            points.erase(points.begin() + i);
            return true;
        }
    return false;
}


bool ParamSet::EraseVector(const string &n) {
    for (uint32_t i = 0; i < vectors.size(); ++i)
        if (vectors[i]->name == n) {
            vectors.erase(vectors.begin() + i);
            return true;
        }
    return false;
}


bool ParamSet::EraseNormal(const string &n) {
    for (uint32_t i = 0; i < normals.size(); ++i)
        if (normals[i]->name == n) {
            normals.erase(normals.begin() + i);
            return true;
        }
    return false;
}


bool ParamSet::EraseSpectrum(const string &n) {
    for (uint32_t i = 0; i < spectra.size(); ++i)
        if (spectra[i]->name == n) {
            spectra.erase(spectra.begin() + i);
            return true;
        }
    return false;
}


bool ParamSet::EraseString(const string &n) {
    for (uint32_t i = 0; i < strings.size(); ++i)
        if (strings[i]->name == n) {
            strings.erase(strings.begin() + i);
            return true;
        }
    return false;
}


bool ParamSet::EraseTexture(const string &n) {
    for (uint32_t i = 0; i < textures.size(); ++i)
        if (textures[i]->name == n) {
            textures.erase(textures.begin() + i);
            return true;
        }
    return false;
}


float ParamSet::FindOneFloat(const string &name, float d) const {
    for (uint32_t i = 0; i < floats.size(); ++i)
        if (floats[i]->name == name && floats[i]->nItems == 1) {
            floats[i]->lookedUp = true;
            return *(floats[i]->data);
        }
    return d;
}


const float *ParamSet::FindFloat(const string &name, int *n) const {
    for (uint32_t i = 0; i < floats.size(); ++i)
        if (floats[i]->name == name) {
            *n = floats[i]->nItems;
            floats[i]->lookedUp = true;
            return floats[i]->data;
        }
    return NULL;
}


const int *ParamSet::FindInt(const string &name, int *nItems) const {
    LOOKUP_PTR(ints);
}


const bool *ParamSet::FindBool(const string &name, int *nItems) const {
    LOOKUP_PTR(bools);
}


int ParamSet::FindOneInt(const string &name, int d) const {
    LOOKUP_ONE(ints);
}


bool ParamSet::FindOneBool(const string &name, bool d) const {
    LOOKUP_ONE(bools);
}


const Point *ParamSet::FindPoint(const string &name, int *nItems) const {
    LOOKUP_PTR(points);
}


Point ParamSet::FindOnePoint(const string &name, const Point &d) const {
    LOOKUP_ONE(points);
}


const Vector *ParamSet::FindVector(const string &name, int *nItems) const {
    LOOKUP_PTR(vectors);
}


Vector ParamSet::FindOneVector(const string &name, const Vector &d) const {
    LOOKUP_ONE(vectors);
}


const Normal *ParamSet::FindNormal(const string &name, int *nItems) const {
    LOOKUP_PTR(normals);
}


Normal ParamSet::FindOneNormal(const string &name, const Normal &d) const {
    LOOKUP_ONE(normals);
}


const Spectrum *ParamSet::FindSpectrum(const string &name, int *nItems) const {
    LOOKUP_PTR(spectra);
}


Spectrum ParamSet::FindOneSpectrum(const string &name, const Spectrum &d) const {
    LOOKUP_ONE(spectra);
}


const string *ParamSet::FindString(const string &name, int *nItems) const {
    LOOKUP_PTR(strings);
}


string ParamSet::FindOneString(const string &name, const string &d) const {
    LOOKUP_ONE(strings);
}


string ParamSet::FindOneFilename(const string &name, const string &d) const {
    string filename = FindOneString(name, "");
    if (filename == "")
        return d;
    filename = AbsolutePath(ResolveFilename(filename));
    return filename;
}


string ParamSet::FindTexture(const string &name) const {
    string d = "";
    LOOKUP_ONE(textures);
}


void ParamSet::ReportUnused() const {
#define CHECK_UNUSED(v) \
    for (i = 0; i < (v).size(); ++i) \
        if (!(v)[i]->lookedUp) \
            Warning("Parameter \"%s\" not used", \
                (v)[i]->name.c_str())
    uint32_t i;
    CHECK_UNUSED(ints);    CHECK_UNUSED(bools);
    CHECK_UNUSED(floats);  CHECK_UNUSED(points);
    CHECK_UNUSED(vectors); CHECK_UNUSED(normals);
    CHECK_UNUSED(spectra); CHECK_UNUSED(strings);
    CHECK_UNUSED(textures);
}


void ParamSet::Clear() {
#define DEL_PARAMS(name) \
    (name).erase((name).begin(), (name).end())
    DEL_PARAMS(ints);    DEL_PARAMS(bools);
    DEL_PARAMS(floats);  DEL_PARAMS(points);
    DEL_PARAMS(vectors); DEL_PARAMS(normals);
    DEL_PARAMS(spectra); DEL_PARAMS(strings);
    DEL_PARAMS(textures);
#undef DEL_PARAMS
}


string ParamSet::ToString() const {
    string ret;
    uint32_t i;
    int j;
    string typeString;
    const int bufLen = 48*1024*1024;
    static char *buf = new char[bufLen];
    char *bufEnd = buf + bufLen;
    for (i = 0; i < ints.size(); ++i) {
        char *bufp = buf;
        *bufp = '\0';
        const Reference<ParamSetItem<int> > &item = ints[i];
        typeString = "integer ";
        // Print _ParamSetItem_ declaration, determine how many to print
        int nPrint = item->nItems;
        ret += string("\"");
        ret += typeString;
        ret += item->name;
        ret += string("\"");
        ret += string(" [");
        for (j = 0; j < nPrint; ++j)
            bufp += snprintf(bufp, bufEnd - bufp, "%d ", item->data[j]);
        ret += buf;
        ret += string("] ");
    }
    for (i = 0; i < bools.size(); ++i) {
        char *bufp = buf;
        *bufp = '\0';
        const Reference<ParamSetItem<bool> > &item = bools[i];
        typeString = "bool ";
        // Print _ParamSetItem_ declaration, determine how many to print
        int nPrint = item->nItems;
        ret += string("\"");
        ret += typeString;
        ret += item->name;
        ret += string("\"");
        ret += string(" [");
        for (j = 0; j < nPrint; ++j)
            bufp += snprintf(bufp, bufEnd - bufp, "\"%s\" ", item->data[j] ? "true" : "false");
        ret += buf;
        ret += string("] ");
    }
    for (i = 0; i < floats.size(); ++i) {
        char *bufp = buf;
        *bufp = '\0';
        const Reference<ParamSetItem<float> > &item = floats[i];
        typeString = "float ";
        // Print _ParamSetItem_ declaration, determine how many to print
        int nPrint = item->nItems;
        ret += string("\"");
        ret += typeString;
        ret += item->name;
        ret += string("\"");
        ret += string(" [");
        for (j = 0; j < nPrint; ++j)
            bufp += snprintf(bufp, bufEnd - bufp, "%.8g ", item->data[j]);
        ret += buf;
        ret += string("] ");
    }
    for (i = 0; i < points.size(); ++i) {
        char *bufp = buf;
        *bufp = '\0';
        const Reference<ParamSetItem<Point> > &item = points[i];
        typeString = "point ";
        // Print _ParamSetItem_ declaration, determine how many to print
        int nPrint = item->nItems;
        ret += string("\"");
        ret += typeString;
        ret += item->name;
        ret += string("\"");
        ret += string(" [");
        for (j = 0; j < nPrint; ++j)
            bufp += snprintf(bufp, bufEnd - bufp, "%.8g %.8g %.8g ", item->data[j].x,
                item->data[j].y, item->data[j].z);
        ret += buf;
        ret += string("] ");
    }
    for (i = 0; i < vectors.size(); ++i) {
        char *bufp = buf;
        *bufp = '\0';
        const Reference<ParamSetItem<Vector> > &item = vectors[i];
        typeString = "vector ";
        // Print _ParamSetItem_ declaration, determine how many to print
        int nPrint = item->nItems;
        ret += string("\"");
        ret += typeString;
        ret += item->name;
        ret += string("\"");
        ret += string(" [");
        for (j = 0; j < nPrint; ++j)
            bufp += snprintf(bufp, bufEnd - bufp, "%.8g %.8g %.8g ", item->data[j].x,
                item->data[j].y, item->data[j].z);
        ret += buf;
        ret += string("] ");
    }
    for (i = 0; i < normals.size(); ++i) {
        char *bufp = buf;
        *bufp = '\0';
        const Reference<ParamSetItem<Normal> > &item = normals[i];
        typeString = "normal ";
        // Print _ParamSetItem_ declaration, determine how many to print
        int nPrint = item->nItems;
        ret += string("\"");
        ret += typeString;
        ret += item->name;
        ret += string("\"");
        ret += string(" [");
        for (j = 0; j < nPrint; ++j)
            bufp += snprintf(bufp, bufEnd - bufp, "%.8g %.8g %.8g ", item->data[j].x,
                item->data[j].y, item->data[j].z);
        ret += buf;
        ret += string("] ");
    }
    for (i = 0; i < strings.size(); ++i) {
        char *bufp = buf;
        *bufp = '\0';
        const Reference<ParamSetItem<string> > &item = strings[i];
        typeString = "string ";
        // Print _ParamSetItem_ declaration, determine how many to print
        int nPrint = item->nItems;
        ret += string("\"");
        ret += typeString;
        ret += item->name;
        ret += string("\"");
        ret += string(" [");
        for (j = 0; j < nPrint; ++j)
            bufp += snprintf(bufp, bufEnd - bufp, "\"%s\" ", item->data[j].c_str());
        ret += buf;
        ret += string("] ");
    }
    for (i = 0; i < textures.size(); ++i) {
        char *bufp = buf;
        *bufp = '\0';
        const Reference<ParamSetItem<string> > &item = textures[i];
        typeString = "texture ";
        // Print _ParamSetItem_ declaration, determine how many to print
        int nPrint = item->nItems;
        ret += string("\"");
        ret += typeString;
        ret += item->name;
        ret += string("\"");
        ret += string(" [");
        for (j = 0; j < nPrint; ++j)
            bufp += snprintf(bufp, bufEnd - bufp, "\"%s\" ", item->data[j].c_str());
        ret += buf;
        ret += string("] ");
    }
    for (i = 0; i < spectra.size(); ++i) {
        char *bufp = buf;
        *bufp = '\0';
        const Reference<ParamSetItem<Spectrum> > &item = spectra[i];
        typeString = "color ";
        // Print _ParamSetItem_ declaration, determine how many to print
        int nPrint = item->nItems;
        ret += string("\"");
        ret += typeString;
        ret += item->name;
        ret += string("\"");
        ret += string(" [");
        for (j = 0; j < nPrint; ++j) {
            float rgb[3];
            item->data[j].ToRGB(rgb);
            bufp += snprintf(bufp, bufEnd - bufp, "%.8g %.8g %.8g ", rgb[0], rgb[1], rgb[2]);
        }
        ret += buf;
        ret += string("] ");
    }
    return ret;
}


// TextureParams Method Definitions
Reference<Texture<Spectrum> >
TextureParams::GetSpectrumTexture(const string &n,
                                  const Spectrum &def) const {
    string name = geomParams.FindTexture(n);
    if (name == "") name = materialParams.FindTexture(n);
    if (name != "") {
        if (spectrumTextures.find(name) != spectrumTextures.end())
            return spectrumTextures[name];
        else
            Error("Couldn't find spectrum texture named \"%s\" "
                  "for parameter \"%s\"", name.c_str(), n.c_str());
    }
    Spectrum val = geomParams.FindOneSpectrum(n,
                                materialParams.FindOneSpectrum(n, def));
    return new ConstantTexture<Spectrum>(val);
}


Reference<Texture<float> > TextureParams::GetFloatTexture(const string &n,
        float def) const {
    string name = geomParams.FindTexture(n);
    if (name == "") name = materialParams.FindTexture(n);
    if (name != "") {
        if (floatTextures.find(name) != floatTextures.end())
            return floatTextures[name];
        else
            Error("Couldn't find float texture named \"%s\" for parameter \"%s\"",
                  name.c_str(), n.c_str());
    }
    float val = geomParams.FindOneFloat(n,
        materialParams.FindOneFloat(n, def));
    return new ConstantTexture<float>(val);
}


