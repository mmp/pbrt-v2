
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

#ifndef PBRT_CORE_PBRT_H
#define PBRT_CORE_PBRT_H

// core/pbrt.h*

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
#if !defined(__APPLE__) && !defined(__OpenBSD__)
#include <malloc.h> // for _alloca, memalign
#endif
#if !defined(WIN32) && !defined(__APPLE__) && !defined(__OpenBSD__)
#include <alloca.h>
#endif
#include <assert.h>
#include <algorithm>
using std::min;
using std::max;
using std::swap;
using std::sort;

// Platform-specific definitions
#ifdef WIN32
#include <float.h>
#define isnan _isnan
#define isinf(f) (!_finite((f)))
#define int8_t __int8
#define uint8_t unsigned __int8
#define int16_t __int16
#define uint16_t unsigned __int16
#define int32_t __int32
#define uint32_t unsigned __int32
#define int64_t __int64
#define uint64_t unsigned __int64
#pragma warning (disable : 4305) // double constant assigned to float
#pragma warning (disable : 4244) // int -> float conversion
#endif // WIN32
#ifdef WIN32
#define isnan _isnan
#define isinf(f) (!_finite((f)))
#endif

// Global Macros
#define ALLOCA(TYPE, COUNT) (TYPE *)alloca((COUNT) * sizeof(TYPE))

// Global Type Declarations
typedef unsigned char u_char;
typedef unsigned short u_short;
typedef unsigned int u_int;
typedef unsigned long u_long;

// Global Forward Declarations
class Timer;
class ProgressReporter;
class MemoryArena;
template <typename T, int logBlockSize = 2> class BlockedArray;
struct Matrix4x4;
class RNG;
class Mutex;
class RWMutex;
class Shape;
class ParamSet;
template <typename T> struct ParamSetItem;
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
//typedef SampledSpectrum Spectrum;
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
#define PBRT_VERSION "2.0alpha1"
#ifdef WIN32
#define alloca _alloca
#endif
#ifndef PBRT_L1_CACHE_LINE_SIZE
#define PBRT_L1_CACHE_LINE_SIZE 64
#endif
#ifdef M_PI
#undef M_PI
#endif
#define M_PI           3.14159265358979323846f
#define INV_PI     0.31830988618379067154f
#define INV_TWOPI  0.15915494309189533577f
#ifndef INFINITY
#define INFINITY FLT_MAX
#endif

// Global Inline Functions
#ifdef NDEBUG
#define Assert(expr) ((void)0)
#else
#define Assert(expr) \
    ((expr) ? (void)0 : \
        Severe("Assertion \"%s\" failed in %s, line %d", \
               #expr, __FILE__, __LINE__))
#endif // NDEBUG
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


inline u_int RoundUpPow2(u_int v) {
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
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
