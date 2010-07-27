
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

#ifndef PBRT_FILTERS_SINC_H
#define PBRT_FILTERS_SINC_H

// filters/sinc.h*
#include "filter.h"

// Sinc Filter Declarations
class LanczosSincFilter : public Filter {
public:
    // LanczosSincFilter Public Methods
    LanczosSincFilter(float xw, float yw, float t)
        : Filter(xw, yw), tau(t) { }
    float Evaluate(float x, float y) const;
    float Sinc1D(float x) const {
        x = fabsf(x);
        if (x < 1e-5) return 1.f;
        if (x > 1.)   return 0.f;
        x *= M_PI;
        float sinc = sinf(x * tau) / (x * tau);
        float lanczos = sinf(x) / x;
        return sinc * lanczos;
    }
private:
    const float tau;
};


LanczosSincFilter *CreateSincFilter(const ParamSet &ps);

#endif // PBRT_FILTERS_SINC_H
