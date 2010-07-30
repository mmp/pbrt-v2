
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

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CORE_SPECTRUM_H
#define PBRT_CORE_SPECTRUM_H

// core/spectrum.h*
#include "pbrt.h"
#include "parallel.h"

// Spectrum Utility Declarations
static const int sampledLambdaStart = 400;
static const int sampledLambdaEnd = 700;
static const int nSpectralSamples = 30;
extern bool SpectrumSamplesSorted(const float *lambda, const float *vals, int n);
extern void SortSpectrumSamples(float *lambda, float *vals, int n);
extern float AverageSpectrumSamples(const float *lambda, const float *vals,
    int n, float lambdaStart, float lambdaEnd);
inline void XYZToRGB(const float xyz[3], float rgb[3]) {
    rgb[0] =  3.240479f*xyz[0] - 1.537150f*xyz[1] - 0.498535f*xyz[2];
    rgb[1] = -0.969256f*xyz[0] + 1.875991f*xyz[1] + 0.041556f*xyz[2];
    rgb[2] =  0.055648f*xyz[0] - 0.204043f*xyz[1] + 1.057311f*xyz[2];
}


inline void RGBToXYZ(const float rgb[3], float xyz[3]) {
    xyz[0] = 0.412453f*rgb[0] + 0.357580f*rgb[1] + 0.180423f*rgb[2];
    xyz[1] = 0.212671f*rgb[0] + 0.715160f*rgb[1] + 0.072169f*rgb[2];
    xyz[2] = 0.019334f*rgb[0] + 0.119193f*rgb[1] + 0.950227f*rgb[2];
}


enum SpectrumType { SPECTRUM_REFLECTANCE, SPECTRUM_ILLUMINANT };
extern void Blackbody(const float *wl, int n, float temp, float *vals);
extern float InterpolateSpectrumSamples(const float *lambda, const float *vals,
                                        int n, float l);

// Spectral Data Declarations
static const int nCIESamples = 471;
extern const float CIE_X[nCIESamples];
extern const float CIE_Y[nCIESamples];
extern const float CIE_Z[nCIESamples];
extern const float CIE_lambda[nCIESamples];
static const int nRGB2SpectSamples = 32;
extern const float RGB2SpectLambda[nRGB2SpectSamples];
extern const float RGBRefl2SpectWhite[nRGB2SpectSamples];
extern const float RGBRefl2SpectCyan[nRGB2SpectSamples];
extern const float RGBRefl2SpectMagenta[nRGB2SpectSamples];
extern const float RGBRefl2SpectYellow[nRGB2SpectSamples];
extern const float RGBRefl2SpectRed[nRGB2SpectSamples];
extern const float RGBRefl2SpectGreen[nRGB2SpectSamples];
extern const float RGBRefl2SpectBlue[nRGB2SpectSamples];
extern const float RGBIllum2SpectWhite[nRGB2SpectSamples];
extern const float RGBIllum2SpectCyan[nRGB2SpectSamples];
extern const float RGBIllum2SpectMagenta[nRGB2SpectSamples];
extern const float RGBIllum2SpectYellow[nRGB2SpectSamples];
extern const float RGBIllum2SpectRed[nRGB2SpectSamples];
extern const float RGBIllum2SpectGreen[nRGB2SpectSamples];
extern const float RGBIllum2SpectBlue[nRGB2SpectSamples];

// Spectrum Declarations
template <int nSamples> class CoefficientSpectrum {
public:
    // CoefficientSpectrum Public Methods
    CoefficientSpectrum(float v = 0.f) {
        for (int i = 0; i < nSamples; ++i)
            c[i] = v;
        Assert(!HasNaNs());
    }
#ifdef DEBUG
    CoefficientSpectrum(const CoefficientSpectrum &s) {
        Assert(!s.HasNaNs());
        for (int i = 0; i < nSamples; ++i)
            c[i] = s.c[i];
    }
    
    CoefficientSpectrum &operator=(const CoefficientSpectrum &s) {
        Assert(!s.HasNaNs());
        for (int i = 0; i < nSamples; ++i)
            c[i] = s.c[i];
        return *this;
    }
#endif // DEBUG
    void Print(FILE *f) const {
        fprintf(f, "[ ");
        for (int i = 0; i < nSamples; ++i) {
            fprintf(f, "%f", c[i]);
            if (i != nSamples-1) fprintf(f, ", ");
        }
        fprintf(f, "]");
    }
    CoefficientSpectrum &operator+=(const CoefficientSpectrum &s2) {
        Assert(!s2.HasNaNs());
        for (int i = 0; i < nSamples; ++i)
            c[i] += s2.c[i];
        return *this;
    }
    CoefficientSpectrum operator+(const CoefficientSpectrum &s2) const {
        Assert(!s2.HasNaNs());
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSamples; ++i)
            ret.c[i] += s2.c[i];
        return ret;
    }
    CoefficientSpectrum operator-(const CoefficientSpectrum &s2) const {
        Assert(!s2.HasNaNs());
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSamples; ++i)
            ret.c[i] -= s2.c[i];
        return ret;
    }
    CoefficientSpectrum operator/(const CoefficientSpectrum &s2) const {
        Assert(!s2.HasNaNs());
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSamples; ++i)
            ret.c[i] /= s2.c[i];
        return ret;
    }
    CoefficientSpectrum operator*(const CoefficientSpectrum &sp) const {
        Assert(!sp.HasNaNs());
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSamples; ++i)
            ret.c[i] *= sp.c[i];
        return ret;
    }
    CoefficientSpectrum &operator*=(const CoefficientSpectrum &sp) {
        Assert(!sp.HasNaNs());
        for (int i = 0; i < nSamples; ++i)
            c[i] *= sp.c[i];
        return *this;
    }
    CoefficientSpectrum operator*(float a) const {
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSamples; ++i)
            ret.c[i] *= a;
        Assert(!ret.HasNaNs());
        return ret;
    }
    CoefficientSpectrum &operator*=(float a) {
        for (int i = 0; i < nSamples; ++i)
            c[i] *= a;
        Assert(!HasNaNs());
        return *this;
    }
    friend inline
    CoefficientSpectrum operator*(float a, const CoefficientSpectrum &s) {
        Assert(!isnan(a) && !s.HasNaNs());
        return s * a;
    }
    CoefficientSpectrum operator/(float a) const {
        Assert(!isnan(a));
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSamples; ++i)
            ret.c[i] /= a;
        Assert(!ret.HasNaNs());
        return ret;
    }
    CoefficientSpectrum &operator/=(float a) {
        Assert(!isnan(a));
        for (int i = 0; i < nSamples; ++i)
            c[i] /= a;
        return *this;
    }
    bool operator==(const CoefficientSpectrum &sp) const {
        for (int i = 0; i < nSamples; ++i)
            if (c[i] != sp.c[i]) return false;
        return true;
    }
    bool operator!=(const CoefficientSpectrum &sp) const {
        return !(*this == sp);
    }
    bool IsBlack() const {
        for (int i = 0; i < nSamples; ++i)
            if (c[i] != 0.) return false;
        return true;
    }
    friend CoefficientSpectrum Sqrt(const CoefficientSpectrum &s) {
        CoefficientSpectrum ret;
        for (int i = 0; i < nSamples; ++i)
            ret.c[i] = sqrtf(s.c[i]);
        Assert(!ret.HasNaNs());
        return ret;
    }
    template <int n> friend inline CoefficientSpectrum<n> Pow(const CoefficientSpectrum<n> &s, float e);
    CoefficientSpectrum operator-() const {
        CoefficientSpectrum ret;
        for (int i = 0; i < nSamples; ++i)
            ret.c[i] = -c[i];
        return ret;
    }
    friend CoefficientSpectrum Exp(const CoefficientSpectrum &s) {
        CoefficientSpectrum ret;
        for (int i = 0; i < nSamples; ++i)
            ret.c[i] = expf(s.c[i]);
        Assert(!ret.HasNaNs());
        return ret;
    }
    CoefficientSpectrum Clamp(float low = 0, float high = INFINITY) const {
        CoefficientSpectrum ret;
        for (int i = 0; i < nSamples; ++i)
            ret.c[i] = ::Clamp(c[i], low, high);
        Assert(!ret.HasNaNs());
        return ret;
    }
    bool HasNaNs() const {
        for (int i = 0; i < nSamples; ++i)
            if (isnan(c[i])) return true;
        return false;
    }
    bool Write(FILE *f) const {
        for (int i = 0; i < nSamples; ++i)
            if (fprintf(f, "%f ", c[i]) < 0) return false;
        return true;
    }
    bool Read(FILE *f) {
        for (int i = 0; i < nSamples; ++i)
            if (fscanf(f, "%f ", &c[i]) != 1) return false;
        return true;
    }
protected:
    // CoefficientSpectrum Protected Data
    float c[nSamples];
};


class SampledSpectrum : public CoefficientSpectrum<nSpectralSamples> {
public:
    // SampledSpectrum Public Methods
    SampledSpectrum(float v = 0.f) {
        for (int i = 0; i < nSpectralSamples; ++i) c[i] = v;
    }
    SampledSpectrum(const CoefficientSpectrum<nSpectralSamples> &v)
        : CoefficientSpectrum<nSpectralSamples>(v) { }
    static SampledSpectrum FromSampled(const float *lambda,
                                       const float *v, int n) {
        // Sort samples if unordered, use sorted for returned spectrum
        if (!SpectrumSamplesSorted(lambda, v, n)) {
            vector<float> slambda(&lambda[0], &lambda[n]);
            vector<float> sv(&v[0], &v[n]);
            SortSpectrumSamples(&slambda[0], &sv[0], n);
            return FromSampled(&slambda[0], &sv[0], n);
        }
        SampledSpectrum r;
        for (int i = 0; i < nSpectralSamples; ++i) {
            // Compute average value of given SPD over $i$th sample's range
            float lambda0 = Lerp(float(i) / float(nSpectralSamples),
                                 sampledLambdaStart, sampledLambdaEnd);
            float lambda1 = Lerp(float(i+1) / float(nSpectralSamples),
                                 sampledLambdaStart, sampledLambdaEnd);
            r.c[i] = AverageSpectrumSamples(lambda, v, n, lambda0, lambda1);
        }
        return r;
    }
    static void Init() {
        // Compute XYZ matching functions for _SampledSpectrum_
        for (int i = 0; i < nSpectralSamples; ++i) {
            float wl0 = Lerp(float(i) / float(nSpectralSamples),
                             sampledLambdaStart, sampledLambdaEnd);
            float wl1 = Lerp(float(i+1) / float(nSpectralSamples),
                             sampledLambdaStart, sampledLambdaEnd);
            X.c[i] = AverageSpectrumSamples(CIE_lambda, CIE_X, nCIESamples,
                                            wl0, wl1);
            Y.c[i] = AverageSpectrumSamples(CIE_lambda, CIE_Y, nCIESamples,
                                            wl0, wl1);
            Z.c[i] = AverageSpectrumSamples(CIE_lambda, CIE_Z, nCIESamples,
                                            wl0, wl1);
            yint += Y.c[i];
        }

        // Compute RGB to spectrum functions for _SampledSpectrum_
        for (int i = 0; i < nSpectralSamples; ++i) {
            float wl0 = Lerp(float(i) / float(nSpectralSamples),
                             sampledLambdaStart, sampledLambdaEnd);
            float wl1 = Lerp(float(i+1) / float(nSpectralSamples),
                             sampledLambdaStart, sampledLambdaEnd);
            rgbRefl2SpectWhite.c[i] = AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectWhite,
                nRGB2SpectSamples, wl0, wl1);
            rgbRefl2SpectCyan.c[i] = AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectCyan,
                nRGB2SpectSamples, wl0, wl1);
            rgbRefl2SpectMagenta.c[i] = AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectMagenta,
                nRGB2SpectSamples, wl0, wl1);
            rgbRefl2SpectYellow.c[i] = AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectYellow,
                nRGB2SpectSamples, wl0, wl1);
            rgbRefl2SpectRed.c[i] = AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectRed,
                nRGB2SpectSamples, wl0, wl1);
            rgbRefl2SpectGreen.c[i] = AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectGreen,
                nRGB2SpectSamples, wl0, wl1);
            rgbRefl2SpectBlue.c[i] = AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectBlue,
                nRGB2SpectSamples, wl0, wl1);
        
            rgbIllum2SpectWhite.c[i] = AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectWhite,
                nRGB2SpectSamples, wl0, wl1);
            rgbIllum2SpectCyan.c[i] = AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectCyan,
                nRGB2SpectSamples, wl0, wl1);
            rgbIllum2SpectMagenta.c[i] = AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectMagenta,
                nRGB2SpectSamples, wl0, wl1);
            rgbIllum2SpectYellow.c[i] = AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectYellow,
                nRGB2SpectSamples, wl0, wl1);
            rgbIllum2SpectRed.c[i] = AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectRed,
                nRGB2SpectSamples, wl0, wl1);
            rgbIllum2SpectGreen.c[i] = AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectGreen,
                nRGB2SpectSamples, wl0, wl1);
            rgbIllum2SpectBlue.c[i] = AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectBlue,
                nRGB2SpectSamples, wl0, wl1);
        }
    }
    void ToXYZ(float xyz[3]) const {
        xyz[0] = xyz[1] = xyz[2] = 0.f;
        for (int i = 0; i < nSpectralSamples; ++i) {
            xyz[0] += X.c[i] * c[i];
            xyz[1] += Y.c[i] * c[i];
            xyz[2] += Z.c[i] * c[i];
        }
        xyz[0] /= yint;
        xyz[1] /= yint;
        xyz[2] /= yint;
    }
    float y() const {
        float yy = 0.f;
        for (int i = 0; i < nSpectralSamples; ++i)
            yy += Y.c[i] * c[i];
        return yy / yint;
    }
    void ToRGB(float rgb[3]) const {
        float xyz[3];
        ToXYZ(xyz);
        XYZToRGB(xyz, rgb);
    }
    RGBSpectrum ToRGBSpectrum() const;
    static SampledSpectrum FromRGB(const float rgb[3],
        SpectrumType type = SPECTRUM_REFLECTANCE);
    static SampledSpectrum FromXYZ(const float xyz[3],
            SpectrumType type = SPECTRUM_REFLECTANCE) {
        float rgb[3];
        XYZToRGB(xyz, rgb);
        return FromRGB(rgb, type);
    }
    SampledSpectrum(const RGBSpectrum &r, SpectrumType type = SPECTRUM_REFLECTANCE);
private:
    // SampledSpectrum Private Data
    static SampledSpectrum X, Y, Z;
    static float yint;
    static SampledSpectrum rgbRefl2SpectWhite, rgbRefl2SpectCyan;
    static SampledSpectrum rgbRefl2SpectMagenta, rgbRefl2SpectYellow;
    static SampledSpectrum rgbRefl2SpectRed, rgbRefl2SpectGreen;
    static SampledSpectrum rgbRefl2SpectBlue;
    static SampledSpectrum rgbIllum2SpectWhite, rgbIllum2SpectCyan;
    static SampledSpectrum rgbIllum2SpectMagenta, rgbIllum2SpectYellow;
    static SampledSpectrum rgbIllum2SpectRed, rgbIllum2SpectGreen;
    static SampledSpectrum rgbIllum2SpectBlue;
};


class RGBSpectrum : public CoefficientSpectrum<3> {
    using CoefficientSpectrum<3>::c;
public:
    // RGBSpectrum Public Methods
    RGBSpectrum(float v = 0.f) : CoefficientSpectrum<3>(v) { }
    RGBSpectrum(const CoefficientSpectrum<3> &v)
        : CoefficientSpectrum<3>(v) { }
    RGBSpectrum(const RGBSpectrum &s, SpectrumType type = SPECTRUM_REFLECTANCE) {
        *this = s;
    }
    static RGBSpectrum FromRGB(const float rgb[3],
            SpectrumType type = SPECTRUM_REFLECTANCE) {
        RGBSpectrum s;
        s.c[0] = rgb[0];
        s.c[1] = rgb[1];
        s.c[2] = rgb[2];
        Assert(!s.HasNaNs());
        return s;
    }
    void ToRGB(float *rgb) const {
        rgb[0] = c[0];
        rgb[1] = c[1];
        rgb[2] = c[2];
    }
    const RGBSpectrum &ToRGBSpectrum() const {
        return *this;
    }
    void ToXYZ(float xyz[3]) const {
        RGBToXYZ(c, xyz);
    }
    static RGBSpectrum FromXYZ(const float xyz[3],
            SpectrumType type = SPECTRUM_REFLECTANCE) {
        RGBSpectrum r;
        XYZToRGB(xyz, r.c);
        return r;
    }
    float y() const {
        const float YWeight[3] = { 0.212671f, 0.715160f, 0.072169f };
        return YWeight[0] * c[0] + YWeight[1] * c[1] + YWeight[2] * c[2];
    }
    static RGBSpectrum FromSampled(const float *lambda, const float *v,
                                   int n) {
        // Sort samples if unordered, use sorted for returned spectrum
        if (!SpectrumSamplesSorted(lambda, v, n)) {
            vector<float> slambda(&lambda[0], &lambda[n]);
            vector<float> sv(&v[0], &v[n]);
            SortSpectrumSamples(&slambda[0], &sv[0], n);
            return FromSampled(&slambda[0], &sv[0], n);
        }
        float xyz[3] = { 0, 0, 0 };
        float yint = 0.f;
        for (int i = 0; i < nCIESamples; ++i) {
            yint += CIE_Y[i];
            float val = InterpolateSpectrumSamples(lambda, v, n,
                                                   CIE_lambda[i]);
            xyz[0] += val * CIE_X[i];
            xyz[1] += val * CIE_Y[i];
            xyz[2] += val * CIE_Z[i];
        }
        xyz[0] /= yint;
        xyz[1] /= yint;
        xyz[2] /= yint;
        return FromXYZ(xyz);
    }
};



// Spectrum Inline Functions
template <int nSamples> inline CoefficientSpectrum<nSamples>
Pow(const CoefficientSpectrum<nSamples> &s, float e) {
    CoefficientSpectrum<nSamples> ret;
    for (int i = 0; i < nSamples; ++i)
        ret.c[i] = powf(s.c[i], e);
    Assert(!ret.HasNaNs());
    return ret;
}


inline Spectrum Lerp(float t, const Spectrum &s1, const Spectrum &s2) {
    return (1.f - t) * s1 + t * s2;
}



#endif // PBRT_CORE_SPECTRUM_H
