
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


// textures/imagemap.cpp*
#include "stdafx.h"
#include "textures/imagemap.h"
#include "imageio.h"

// ImageTexture Method Definitions
template <typename Tmemory, typename Treturn>
ImageTexture<Tmemory, Treturn>::ImageTexture(TextureMapping2D *m,
        const string &filename, bool doTrilinear, float maxAniso,
        ImageWrap wrapMode, float scale, float gamma) {
    mapping = m;
    mipmap = GetTexture(filename, doTrilinear, maxAniso,
                        wrapMode, scale, gamma);
}


template <typename Tmemory, typename Treturn>
    ImageTexture<Tmemory, Treturn>::~ImageTexture() {
    delete mapping;
}


template <typename Tmemory, typename Treturn> MIPMap<Tmemory> *
ImageTexture<Tmemory, Treturn>::GetTexture(const string &filename,
        bool doTrilinear, float maxAniso, ImageWrap wrap,
        float scale, float gamma) {
    // Look for texture in texture cache
    TexInfo texInfo(filename, doTrilinear, maxAniso, wrap, scale, gamma);
    if (textures.find(texInfo) != textures.end())
        return textures[texInfo];
    int width, height;
    RGBSpectrum *texels = ReadImage(filename, &width, &height);
    MIPMap<Tmemory> *ret = NULL;
    if (texels) {
        // Convert texels to type _Tmemory_ and create _MIPMap_
        Tmemory *convertedTexels = new Tmemory[width*height];
        for (int i = 0; i < width*height; ++i)
            convertIn(texels[i], &convertedTexels[i], scale, gamma);
        ret = new MIPMap<Tmemory>(width, height, convertedTexels, doTrilinear,
                                  maxAniso, wrap);
        delete[] texels;
        delete[] convertedTexels;
    }
    else {
        // Create one-valued _MIPMap_
        Tmemory *oneVal = new Tmemory[1];
        oneVal[0] = powf(scale, gamma);
        ret = new MIPMap<Tmemory>(1, 1, oneVal);
        delete[] oneVal;
    }
    textures[texInfo] = ret;
    PBRT_LOADED_IMAGE_MAP(const_cast<char *>(filename.c_str()), width, height, sizeof(Tmemory), ret);
    return ret;
}


template <typename Tmemory, typename Treturn>
    std::map<TexInfo,
             MIPMap<Tmemory> *> ImageTexture<Tmemory, Treturn>::textures;
template <typename Tmemory, typename Treturn> Treturn
ImageTexture<Tmemory,
             Treturn>::Evaluate(const DifferentialGeometry &dg) const {
    float s, t, dsdx, dtdx, dsdy, dtdy;
    mapping->Map(dg, &s, &t, &dsdx, &dtdx, &dsdy, &dtdy);
    Tmemory mem = mipmap->Lookup(s, t, dsdx, dtdx, dsdy, dtdy);
    Treturn ret;
    convertOut(mem, &ret);
    return ret;
}


ImageTexture<float, float> *CreateImageFloatTexture(const Transform &tex2world,
        const TextureParams &tp) {
    // Initialize 2D texture mapping _map_ from _tp_
    TextureMapping2D *map = NULL;
    string type = tp.FindString("mapping", "uv");
    if (type == "uv") {
        float su = tp.FindFloat("uscale", 1.);
        float sv = tp.FindFloat("vscale", 1.);
        float du = tp.FindFloat("udelta", 0.);
        float dv = tp.FindFloat("vdelta", 0.);
        map = new UVMapping2D(su, sv, du, dv);
    }
    else if (type == "spherical") map = new SphericalMapping2D(Inverse(tex2world));
    else if (type == "cylindrical") map = new CylindricalMapping2D(Inverse(tex2world));
    else if (type == "planar")
        map = new PlanarMapping2D(tp.FindVector("v1", Vector(1,0,0)),
            tp.FindVector("v2", Vector(0,1,0)),
            tp.FindFloat("udelta", 0.f), tp.FindFloat("vdelta", 0.f));
    else {
        Error("2D texture mapping \"%s\" unknown", type.c_str());
        map = new UVMapping2D;
    }

    // Initialize _ImageTexture_ parameters
    float maxAniso = tp.FindFloat("maxanisotropy", 8.f);
    bool trilerp = tp.FindBool("trilinear", false);
    string wrap = tp.FindString("wrap", "repeat");
    ImageWrap wrapMode = TEXTURE_REPEAT;
    if (wrap == "black") wrapMode = TEXTURE_BLACK;
    else if (wrap == "clamp") wrapMode = TEXTURE_CLAMP;
    float scale = tp.FindFloat("scale", 1.f);
    float gamma = tp.FindFloat("gamma", 1.f);
    return new ImageTexture<float, float>(map, tp.FindFilename("filename"),
        trilerp, maxAniso, wrapMode, scale, gamma);
}



ImageTexture<RGBSpectrum, Spectrum> *CreateImageSpectrumTexture(const Transform &tex2world,
        const TextureParams &tp) {
    // Initialize 2D texture mapping _map_ from _tp_
    TextureMapping2D *map = NULL;
    string type = tp.FindString("mapping", "uv");
    if (type == "uv") {
        float su = tp.FindFloat("uscale", 1.);
        float sv = tp.FindFloat("vscale", 1.);
        float du = tp.FindFloat("udelta", 0.);
        float dv = tp.FindFloat("vdelta", 0.);
        map = new UVMapping2D(su, sv, du, dv);
    }
    else if (type == "spherical") map = new SphericalMapping2D(Inverse(tex2world));
    else if (type == "cylindrical") map = new CylindricalMapping2D(Inverse(tex2world));
    else if (type == "planar")
        map = new PlanarMapping2D(tp.FindVector("v1", Vector(1,0,0)),
            tp.FindVector("v2", Vector(0,1,0)),
            tp.FindFloat("udelta", 0.f), tp.FindFloat("vdelta", 0.f));
    else {
        Error("2D texture mapping \"%s\" unknown", type.c_str());
        map = new UVMapping2D;
    }

    // Initialize _ImageTexture_ parameters
    float maxAniso = tp.FindFloat("maxanisotropy", 8.f);
    bool trilerp = tp.FindBool("trilinear", false);
    string wrap = tp.FindString("wrap", "repeat");
    ImageWrap wrapMode = TEXTURE_REPEAT;
    if (wrap == "black") wrapMode = TEXTURE_BLACK;
    else if (wrap == "clamp") wrapMode = TEXTURE_CLAMP;
    float scale = tp.FindFloat("scale", 1.f);
    float gamma = tp.FindFloat("gamma", 1.f);
    return new ImageTexture<RGBSpectrum, Spectrum>(map, tp.FindFilename("filename"),
        trilerp, maxAniso, wrapMode, scale, gamma);
}


