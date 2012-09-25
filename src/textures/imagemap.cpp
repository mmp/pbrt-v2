
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


