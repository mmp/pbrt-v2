
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

#ifndef PBRT_CORE_TIMER_H
#define PBRT_CORE_TIMER_H

// core/timer.h*
#include "pbrt.h"
#if defined (PBRT_IS_WINDOWS)
#include <windows.h>
#if (_MSC_VER >= 1400)
#include <stdio.h>
#define snprintf _snprintf
#endif
#else
#include <sys/time.h>
#endif

// Timer Declarations
class Timer {
public:
    // Public Timer Methods
    Timer();
    
    void Start();
    void Stop();
    void Reset();
    
    double Time();
private:
    // Private Timer Data
    double time0, elapsed;
    bool running;
    double GetTime();
#if defined(PBRT_IS_WINDOWS)
    // Private Windows Timer Data
    LARGE_INTEGER performance_counter, performance_frequency;
    double one_over_frequency;
#else
    // Private UNIX Timer Data
    struct timeval timeofday;
#endif
};



#endif // PBRT_CORE_TIMER_H
