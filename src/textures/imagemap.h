
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

#ifndef PBRT_TEXTURES_IMAGEMAP_H
#define PBRT_TEXTURES_IMAGEMAP_H

// textures/imagemap.h*
#include "pbrt.h"
#include "texture.h"
#include "mipmap.h"
#include "paramset.h"
#include <map>

// TexInfo Declarations
struct TexInfo {
    TexInfo(const string &f, bool dt, float ma, ImageWrap wm, float sc, float ga)
        : filename(f), doTrilinear(dt), maxAniso(ma), wrapMode(wm), scale(sc), gamma(ga) { }
    string filename;
    bool doTrilinear;
    float maxAniso;
    ImageWrap wrapMode;
    float scale, gamma;
    bool operator<(const TexInfo &t2) const {
        if (filename != t2.filename) return filename < t2.filename;
        if (doTrilinear != t2.doTrilinear) return doTrilinear < t2.doTrilinear;
        if (maxAniso != t2.maxAniso) return maxAniso < t2.maxAniso;
        if (scale != t2.scale) return scale < t2.scale;
        if (gamma != t2.gamma) return gamma < t2.gamma;
        return wrapMode < t2.wrapMode;
    }
};



// ImageTexture Declarations
template <typename Tmemory, typename Treturn>
    class ImageTexture : public Texture<Treturn> {
public:
    // ImageTexture Public Methods
    ImageTexture(TextureMapping2D *m, const string &filename, bool doTri,
                 float maxAniso, ImageWrap wm, float scale, float gamma);
    Treturn Evaluate(const DifferentialGeometry &) const;
    ~ImageTexture();
    static void ClearCache() {
        typename std::map<TexInfo, MIPMap<Tmemory> *>::iterator iter;
        iter = textures.begin();
        while (iter != textures.end()) {
            delete iter->second;
            ++iter;
        }
        textures.erase(textures.begin(), textures.end());
    }
private:
    // ImageTexture Private Methods
    static MIPMap<Tmemory> *GetTexture(const string &filename,
        bool doTrilinear, float maxAniso, ImageWrap wm, float scale, float gamma);
    static void convertIn(const RGBSpectrum &from, RGBSpectrum *to,
                          float scale, float gamma) {
        *to = Pow(scale * from, gamma);
    }
    static void convertIn(const RGBSpectrum &from, float *to,
                          float scale, float gamma) {
        *to = powf(scale * from.y(), gamma);
    }
    static void convertOut(const RGBSpectrum &from, Spectrum *to) {
        float rgb[3];
        from.ToRGB(rgb);
        *to = Spectrum::FromRGB(rgb);
    }
    static void convertOut(float from, float *to) {
        *to = from;
    }

    // ImageTexture Private Data
    MIPMap<Tmemory> *mipmap;
    TextureMapping2D *mapping;
    static std::map<TexInfo, MIPMap<Tmemory> *> textures;
};


ImageTexture<float, float> *CreateImageFloatTexture(const Transform &tex2world,
        const TextureParams &tp);
ImageTexture<RGBSpectrum, Spectrum> *CreateImageSpectrumTexture(const Transform &tex2world,
        const TextureParams &tp);

#endif // PBRT_TEXTURES_IMAGEMAP_H
