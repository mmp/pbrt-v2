
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
#define NOMINMAX 
#pragma once
#endif

#ifndef PBRT_CORE_PBRT_H
#define PBRT_CORE_PBRT_H

// core/pbrt.h*

#if defined(_WIN32) || defined(_WIN64)
#define PBRT_IS_WINDOWS
#elif defined(__linux__)
#define PBRT_IS_LINUX
#elif defined(__APPLE__)
  #define PBRT_IS_APPLE
  #if !(defined(__i386__) || defined(__amd64__))
  #define PBRT_IS_APPLE_PPC
  #else
  #define PBRT_IS_APPLE_X86
  #endif
#elif defined(__OpenBSD__)
#define PBRT_IS_OPENBSD
#endif

// Global Include Files
#include <math.h>
#include <stdlib.h>
#define _GNU_SOURCE 1
#include <stdio.h>
#include <string.h>
#include <string>
using std::string;
#include <vector>
using std::vector;
#include "error.h"
#if !defined(PBRT_IS_APPLE) && !defined(PBRT_IS_OPENBSD)
#include <malloc.h> // for _alloca, memalign
#endif
#if !defined(PBRT_IS_WINDOWS) && !defined(PBRT_IS_APPLE) && !defined(PBRT_IS_OPENBSD)
#include <alloca.h>
#endif
#include <assert.h>
#include <algorithm>
using std::min;
using std::max;
using std::swap;
using std::sort;

// Platform-specific definitions
#if defined(PBRT_IS_WINDOWS)
#include <float.h>
#define isnan _isnan
#define isinf(f) (!_finite((f)))
typedef __int8 int8_t;
typedef unsigned __int8 uint8_t;
typedef __int16 int16_t;
typedef unsigned __int16 uint16_t;
typedef __int32 int32_t;
typedef unsigned __int32 uint32_t;
typedef __int64 int64_t;
typedef unsigned __int64 uint64_t;
#pragma warning (disable : 4305) // double constant assigned to float
#pragma warning (disable : 4244) // int -> float conversion
#pragma warning (disable : 4267) // size_t -> unsigned int conversion
#endif

#if defined(PBRT_IS_LINUX) || defined(PBRT_IS_APPLE)
#include <stdint.h>
#endif // PBRT_IS_LINUX || PBRT_IS_APPLE
#if defined(PBRT_IS_WINDOWS)
#define isnan _isnan
#define isinf(f) (!_finite((f)))
#endif

// Global Macros
#define ALLOCA(TYPE, COUNT) (TYPE *)alloca((COUNT) * sizeof(TYPE))

// Global Forward Declarations
class RNG;
class Timer;
class ProgressReporter;
class MemoryArena;
template <typename T, int logBlockSize = 2> class BlockedArray;
struct Matrix4x4;
class Mutex;
class RWMutex;
class Shape;
class ParamSet;
template <typename T> struct ParamSetItem;
struct Options {
    Options() { nCores = 0;
                quickRender = quiet = openWindow = verbose = false;
                imageFile = ""; }
    int nCores;
    bool quickRender;
    bool quiet, verbose;
    bool openWindow;
    string imageFile;
};


extern Options PbrtOptions;
class TextureParams;
class Scene;
class Renderer;
class Vector;
class Point;
class Normal;
class Ray;
class RayDifferential;
class BBox;
class Transform;
struct DifferentialGeometry;
class Primitive;
struct Intersection;
class GeometricPrimitive;
template <int nSamples> class CoefficientSpectrum;
class RGBSpectrum;
class SampledSpectrum;
typedef RGBSpectrum Spectrum;
// typedef SampledSpectrum Spectrum;
class Camera;
class ProjectiveCamera;
class Sampler;
struct CameraSample;
struct Sample;
class Filter;
class Film;
class BxDF;
class BRDF;
class BTDF;
class BSDF;
class Material;
template <typename T> class Texture;
class VolumeRegion;
class BSSRDF;
class Light;
struct VisibilityTester;
class AreaLight;
struct Distribution1D;
struct Distribution2D;
struct BSDFSample;
struct BSDFSampleOffsets;
struct LightSample;
struct LightSampleOffsets;
class SurfaceIntegrator;
class Integrator;
class VolumeIntegrator;

// Global Constants
#define PBRT_VERSION "2.0.0"
#ifdef M_PI
#undef M_PI
#endif
#define M_PI       3.14159265358979323846f
#define INV_PI     0.31830988618379067154f
#define INV_TWOPI  0.15915494309189533577f
#define INV_FOURPI 0.07957747154594766788f
#ifndef INFINITY
#define INFINITY FLT_MAX
#endif
#if defined(PBRT_IS_WINDOWS)
#define alloca _alloca
#endif
#ifndef PBRT_L1_CACHE_LINE_SIZE
#define PBRT_L1_CACHE_LINE_SIZE 64
#endif
#ifndef PBRT_POINTER_SIZE
#if defined(__amd64__) || defined(_M_X64)
#define PBRT_POINTER_SIZE 8
#elif defined(__i386__) || defined(_M_IX86)
#define PBRT_POINTER_SIZE 4
#endif
#endif
#ifndef PBRT_HAS_64_BIT_ATOMICS
#if (PBRT_POINTER_SIZE == 8)
#define PBRT_HAS_64_BIT_ATOMICS
#endif
#endif // PBRT_HAS_64_BIT_ATOMICS

// Global Inline Functions
inline float Lerp(float t, float v1, float v2) {
    return (1.f - t) * v1 + t * v2;
}


inline float Clamp(float val, float low, float high) {
    if (val < low) return low;
    else if (val > high) return high;
    else return val;
}


inline int Clamp(int val, int low, int high) {
    if (val < low) return low;
    else if (val > high) return high;
    else return val;
}


inline int Mod(int a, int b) {
    int n = int(a/b);
    a -= n*b;
    if (a < 0) a += b;
    return a;
}


inline float Radians(float deg) {
    return ((float)M_PI/180.f) * deg;
}


inline float Degrees(float rad) {
    return (180.f/(float)M_PI) * rad;
}


inline float Log2(float x) {
    static float invLog2 = 1.f / logf(2.f);
    return logf(x) * invLog2;
}


inline int Floor2Int(float val);
inline int Log2Int(float v) {
    return Floor2Int(Log2(v));
}


inline bool IsPowerOf2(int v) {
    return (v & (v - 1)) == 0;
}


inline uint32_t RoundUpPow2(uint32_t v) {
    v--;
    v |= v >> 1;    v |= v >> 2;
    v |= v >> 4;    v |= v >> 8;
    v |= v >> 16;
    return v+1;
}


inline int Floor2Int(float val) {
    return (int)floorf(val);
}


inline int Round2Int(float val) {
    return Floor2Int(val + 0.5f);
}


inline int Float2Int(float val) {
    return (int)val;
}


inline int Ceil2Int(float val) {
    return (int)ceilf(val);
}


#ifdef NDEBUG
#define Assert(expr) ((void)0)
#else
#define Assert(expr) \
    ((expr) ? (void)0 : \
        Severe("Assertion \"%s\" failed in %s, line %d", \
               #expr, __FILE__, __LINE__))
#endif // NDEBUG
inline bool Quadratic(float A, float B, float C, float *t0, float *t1) {
    // Find quadratic discriminant
    float discrim = B * B - 4.f * A * C;
    if (discrim <= 0.) return false;
    float rootDiscrim = sqrtf(discrim);

    // Compute quadratic _t_ values
    float q;
    if (B < 0) q = -.5f * (B - rootDiscrim);
    else       q = -.5f * (B + rootDiscrim);
    *t0 = q / A;
    *t1 = C / q;
    if (*t0 > *t1) swap(*t0, *t1);
    return true;
}



#endif // PBRT_CORE_PBRT_H
