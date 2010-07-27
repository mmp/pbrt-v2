
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


// filters/mitchell.cpp*
#include "stdafx.h"
#include "filters/mitchell.h"
#include "paramset.h"

// Mitchell Filter Method Definitions
float MitchellFilter::Evaluate(float x, float y) const {
    return Mitchell1D(x * invXWidth) * Mitchell1D(y * invYWidth);
}


MitchellFilter *CreateMitchellFilter(const ParamSet &ps) {
    // Find common filter parameters
    float xw = ps.FindOneFloat("xwidth", 2.f);
    float yw = ps.FindOneFloat("ywidth", 2.f);
    float B = ps.FindOneFloat("B", 1.f/3.f);
    float C = ps.FindOneFloat("C", 1.f/3.f);
    return new MitchellFilter(B, C, xw, yw);
}


