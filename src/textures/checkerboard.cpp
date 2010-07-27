
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


// textures/checkerboard.cpp*
#include "stdafx.h"
#include "textures/checkerboard.h"

// CheckerboardTexture Method Definitions
Texture<float> *CreateCheckerboardFloatTexture(const Transform &tex2world,
        const TextureParams &tp) {
    int dim = tp.FindInt("dimension", 2);
    if (dim != 2 && dim != 3) {
        Error("%d dimensional checkerboard texture not supported", dim);
        return NULL;
    }
    Reference<Texture<float> > tex1 = tp.GetFloatTexture("tex1", 1.f);
    Reference<Texture<float> > tex2 = tp.GetFloatTexture("tex2", 0.f);
    if (dim == 2) {
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
        string aamode = tp.FindString("aamode", "closedform");
        return new Checkerboard2DTexture<float>(map, tex1, tex2, aamode);
    }
    else {
        // Initialize 3D texture mapping _map_ from _tp_
        TextureMapping3D *map = new IdentityMapping3D(tex2world);
        return new Checkerboard3DTexture<float>(map, tex1, tex2);
    }
}



Texture<Spectrum> *CreateCheckerboardSpectrumTexture(const Transform &tex2world,
        const TextureParams &tp) {
    int dim = tp.FindInt("dimension", 2);
    if (dim != 2 && dim != 3) {
        Error("%d dimensional checkerboard texture not supported", dim);
        return NULL;
    }
    Reference<Texture<Spectrum> > tex1 = tp.GetSpectrumTexture("tex1", 1.f);
    Reference<Texture<Spectrum> > tex2 = tp.GetSpectrumTexture("tex2", 0.f);
    if (dim == 2) {
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
        string aamode = tp.FindString("aamode", "closedform");
        return new Checkerboard2DTexture<Spectrum>(map, tex1, tex2, aamode);
    }
    else {
        // Initialize 3D texture mapping _map_ from _tp_
        TextureMapping3D *map = new IdentityMapping3D(tex2world);
        return new Checkerboard3DTexture<Spectrum>(map, tex1, tex2);
    }
}


