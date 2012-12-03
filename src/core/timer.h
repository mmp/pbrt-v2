
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
