
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


// core/texture.cpp*
#include "stdafx.h"
#include "texture.h"
#include "shape.h"

// Texture Inline Functions
inline float SmoothStep(float min, float max, float value) {
    float v = Clamp((value - min) / (max - min), 0.f, 1.f);
    return v * v * (-2.f * v  + 3.f);
}



// Texture Forward Declarations
inline float Grad(int x, int y, int z, float dx, float dy, float dz);
inline float NoiseWeight(float t);

// Perlin Noise Data
#define NOISE_PERM_SIZE 256
static int NoisePerm[2 * NOISE_PERM_SIZE] = {
    151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96,
    53, 194, 233, 7, 225, 140, 36, 103, 30, 69, 142,
    // Remainder of the noise permutation table
    8, 99, 37, 240, 21, 10, 23,
       190,  6, 148, 247, 120, 234, 75, 0, 26, 197, 62, 94, 252, 219, 203, 117, 35, 11, 32, 57, 177, 33,
       88, 237, 149, 56, 87, 174, 20, 125, 136, 171, 168,  68, 175, 74, 165, 71, 134, 139, 48, 27, 166,
       77, 146, 158, 231, 83, 111, 229, 122, 60, 211, 133, 230, 220, 105, 92, 41, 55, 46, 245, 40, 244,
       102, 143, 54,  65, 25, 63, 161,  1, 216, 80, 73, 209, 76, 132, 187, 208,  89, 18, 169, 200, 196,
       135, 130, 116, 188, 159, 86, 164, 100, 109, 198, 173, 186,  3, 64, 52, 217, 226, 250, 124, 123,
       5, 202, 38, 147, 118, 126, 255, 82, 85, 212, 207, 206, 59, 227, 47, 16, 58, 17, 182, 189, 28, 42,
       223, 183, 170, 213, 119, 248, 152,  2, 44, 154, 163,  70, 221, 153, 101, 155, 167,  43, 172, 9,
       129, 22, 39, 253,  19, 98, 108, 110, 79, 113, 224, 232, 178, 185,  112, 104, 218, 246, 97, 228,
       251, 34, 242, 193, 238, 210, 144, 12, 191, 179, 162, 241,  81, 51, 145, 235, 249, 14, 239, 107,
       49, 192, 214,  31, 181, 199, 106, 157, 184,  84, 204, 176, 115, 121, 50, 45, 127,  4, 150, 254,
       138, 236, 205, 93, 222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180,
       151, 160, 137, 91, 90, 15,
       131, 13, 201, 95, 96, 53, 194, 233, 7, 225, 140, 36, 103, 30, 69, 142, 8, 99, 37, 240, 21, 10, 23,
       190,  6, 148, 247, 120, 234, 75, 0, 26, 197, 62, 94, 252, 219, 203, 117, 35, 11, 32, 57, 177, 33,
       88, 237, 149, 56, 87, 174, 20, 125, 136, 171, 168,  68, 175, 74, 165, 71, 134, 139, 48, 27, 166,
       77, 146, 158, 231, 83, 111, 229, 122, 60, 211, 133, 230, 220, 105, 92, 41, 55, 46, 245, 40, 244,
       102, 143, 54,  65, 25, 63, 161,  1, 216, 80, 73, 209, 76, 132, 187, 208,  89, 18, 169, 200, 196,
       135, 130, 116, 188, 159, 86, 164, 100, 109, 198, 173, 186,  3, 64, 52, 217, 226, 250, 124, 123,
       5, 202, 38, 147, 118, 126, 255, 82, 85, 212, 207, 206, 59, 227, 47, 16, 58, 17, 182, 189, 28, 42,
       223, 183, 170, 213, 119, 248, 152,  2, 44, 154, 163,  70, 221, 153, 101, 155, 167,  43, 172, 9,
       129, 22, 39, 253,  19, 98, 108, 110, 79, 113, 224, 232, 178, 185,  112, 104, 218, 246, 97, 228,
       251, 34, 242, 193, 238, 210, 144, 12, 191, 179, 162, 241,  81, 51, 145, 235, 249, 14, 239, 107,
       49, 192, 214,  31, 181, 199, 106, 157, 184,  84, 204, 176, 115, 121, 50, 45, 127,  4, 150, 254,
       138, 236, 205, 93, 222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180
};



// Texture Method Definitions
UVMapping2D::UVMapping2D(float ssu, float ssv, float ddu, float ddv)
    : su(ssu), sv(ssv), du(ddu), dv(ddv) { }
void UVMapping2D::Map(const DifferentialGeometry &dg,
                      float *s, float *t, float *dsdx, float *dtdx,
                      float *dsdy, float *dtdy) const {
    *s = su * dg.u + du;
    *t = sv * dg.v + dv;
    // Compute texture differentials for 2D identity mapping
    *dsdx = su * dg.dudx;
    *dtdx = sv * dg.dvdx;
    *dsdy = su * dg.dudy;
    *dtdy = sv * dg.dvdy;
}


void SphericalMapping2D::Map(const DifferentialGeometry &dg,
        float *s, float *t, float *dsdx, float *dtdx,
        float *dsdy, float *dtdy) const {
    sphere(dg.p, s, t);
    // Compute texture coordinate differentials for sphere $(u,v)$ mapping
    float sx, tx, sy, ty;
    const float delta = .1f;
    sphere(dg.p + delta * dg.dpdx, &sx, &tx);
    *dsdx = (sx - *s) / delta;
    *dtdx = (tx - *t) / delta;
    if (*dtdx > .5) *dtdx = 1.f - *dtdx;
    else if (*dtdx < -.5f) *dtdx = -(*dtdx + 1);
    sphere(dg.p + delta * dg.dpdy, &sy, &ty);
    *dsdy = (sy - *s) / delta;
    *dtdy = (ty - *t) / delta;
    if (*dtdy > .5) *dtdy = 1.f - *dtdy;
    else if (*dtdy < -.5f) *dtdy = -(*dtdy + 1);
}


void SphericalMapping2D::sphere(const Point &p, float *s, float *t) const {
    Vector vec = Normalize(WorldToTexture(p) - Point(0,0,0));
    float theta = SphericalTheta(vec);
    float phi = SphericalPhi(vec);
    *s = theta * INV_PI;
    *t = phi * INV_TWOPI;
}


void CylindricalMapping2D::Map(const DifferentialGeometry &dg,
        float *s, float *t, float *dsdx, float *dtdx,
        float *dsdy, float *dtdy) const {
    cylinder(dg.p, s, t);
    // Compute texture coordinate differentials for cylinder $(u,v)$ mapping
    float sx, tx, sy, ty;
    const float delta = .01f;
    cylinder(dg.p + delta * dg.dpdx, &sx, &tx);
    *dsdx = (sx - *s) / delta;
    *dtdx = (tx - *t) / delta;
    if (*dtdx > .5) *dtdx = 1.f - *dtdx;
    else if (*dtdx < -.5f) *dtdx = -(*dtdx + 1);
    cylinder(dg.p + delta * dg.dpdy, &sy, &ty);
    *dsdy = (sy - *s) / delta;
    *dtdy = (ty - *t) / delta;
    if (*dtdy > .5) *dtdy = 1.f - *dtdy;
    else if (*dtdy < -.5f) *dtdy = -(*dtdy + 1);
}


void PlanarMapping2D::Map(const DifferentialGeometry &dg,
        float *s, float *t, float *dsdx, float *dtdx,
        float *dsdy, float *dtdy) const {
    Vector vec = dg.p - Point(0,0,0);
    *s = ds + Dot(vec, vs);
    *t = dt + Dot(vec, vt);
    *dsdx = Dot(dg.dpdx, vs);
    *dtdx = Dot(dg.dpdx, vt);
    *dsdy = Dot(dg.dpdy, vs);
    *dtdy = Dot(dg.dpdy, vt);
}


Point IdentityMapping3D::Map(const DifferentialGeometry &dg,
                             Vector *dpdx, Vector *dpdy) const {
    *dpdx = WorldToTexture(dg.dpdx);
    *dpdy = WorldToTexture(dg.dpdy);
    return WorldToTexture(dg.p);
}


float Noise(float x, float y, float z) {
    // Compute noise cell coordinates and offsets
    int ix = Floor2Int(x), iy = Floor2Int(y), iz = Floor2Int(z);
    float dx = x - ix, dy = y - iy, dz = z - iz;

    // Compute gradient weights
    ix &= (NOISE_PERM_SIZE-1);
    iy &= (NOISE_PERM_SIZE-1);
    iz &= (NOISE_PERM_SIZE-1);
    float w000 = Grad(ix,   iy,   iz,   dx,   dy,   dz);
    float w100 = Grad(ix+1, iy,   iz,   dx-1, dy,   dz);
    float w010 = Grad(ix,   iy+1, iz,   dx,   dy-1, dz);
    float w110 = Grad(ix+1, iy+1, iz,   dx-1, dy-1, dz);
    float w001 = Grad(ix,   iy,   iz+1, dx,   dy,   dz-1);
    float w101 = Grad(ix+1, iy,   iz+1, dx-1, dy,   dz-1);
    float w011 = Grad(ix,   iy+1, iz+1, dx,   dy-1, dz-1);
    float w111 = Grad(ix+1, iy+1, iz+1, dx-1, dy-1, dz-1);

    // Compute trilinear interpolation of weights
    float wx = NoiseWeight(dx), wy = NoiseWeight(dy), wz = NoiseWeight(dz);
    float x00 = Lerp(wx, w000, w100);
    float x10 = Lerp(wx, w010, w110);
    float x01 = Lerp(wx, w001, w101);
    float x11 = Lerp(wx, w011, w111);
    float y0 = Lerp(wy, x00, x10);
    float y1 = Lerp(wy, x01, x11);
    return Lerp(wz, y0, y1);
}


float Noise(const Point &P) { return Noise(P.x, P.y, P.z); }
inline float Grad(int x, int y, int z, float dx, float dy, float dz) {
    int h = NoisePerm[NoisePerm[NoisePerm[x]+y]+z];
    h &= 15;
    float u = h<8 || h==12 || h==13 ? dx : dy;
    float v = h<4 || h==12 || h==13 ? dy : dz;
    return ((h&1) ? -u : u) + ((h&2) ? -v : v);
}


inline float NoiseWeight(float t) {
    float t3 = t*t*t;
    float t4 = t3*t;
    return 6.f*t4*t - 15.f*t4 + 10.f*t3;
}


float FBm(const Point &P, const Vector &dpdx, const Vector &dpdy,
          float omega, int maxOctaves) {
    // Compute number of octaves for antialiased FBm
    float s2 = max(dpdx.LengthSquared(), dpdy.LengthSquared());
    float foctaves = min((float)maxOctaves, 1.f - .5f * Log2(s2));
    int octaves = Floor2Int(foctaves);

    // Compute sum of octaves of noise for FBm
    float sum = 0., lambda = 1., o = 1.;
    for (int i = 0; i < octaves; ++i) {
        sum += o * Noise(lambda * P);
        lambda *= 1.99f;
        o *= omega;
    }
    float partialOctave = foctaves - octaves;
    sum += o * SmoothStep(.3f, .7f, partialOctave) * Noise(lambda * P);
    return sum;
}


float Turbulence(const Point &P, const Vector &dpdx, const Vector &dpdy,
                 float omega, int maxOctaves) {
    // Compute number of octaves for antialiased FBm
    float s2 = max(dpdx.LengthSquared(), dpdy.LengthSquared());
    float foctaves = min((float)maxOctaves, 1.f - .5f * Log2(s2));
    int octaves = Floor2Int(foctaves);

    // Compute sum of octaves of noise for turbulence
    float sum = 0., lambda = 1., o = 1.;
    for (int i = 0; i < octaves; ++i) {
        sum += o * fabsf(Noise(lambda * P));
        lambda *= 1.99f;
        o *= omega;
    }
    float partialOctave = foctaves - octaves;
    sum += o * SmoothStep(.3f, .7f, partialOctave) *
           fabsf(Noise(lambda * P));

    // finally, add in value to account for average value of fabsf(Noise())
    // (~0.2) for the remaining octaves...
    sum += (maxOctaves - foctaves) * 0.2f;

    return sum;
}



// Texture Function Definitions
float Lanczos(float x, float tau) {
    x = fabsf(x);
    if (x < 1e-5) return 1;
    if (x > 1.)    return 0;
    x *= M_PI;
    float s = sinf(x * tau) / (x * tau);
    float lanczos = sinf(x) / x;
    return s * lanczos;
}


