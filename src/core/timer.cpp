
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


// core/timer.cpp*
#include "stdafx.h"
#include "timer.h"

// Timer Method Definitions
Timer::Timer()
{
#if defined( PBRT_IS_WINDOWS )
    // Windows Timer Initialization
    QueryPerformanceFrequency( &performance_frequency );
    one_over_frequency = 1.0/((double)performance_frequency.QuadPart);
#endif
    time0 = elapsed = 0;
    running = 0;
}




double Timer::GetTime()
{
#if defined( PBRT_IS_WINDOWS )
    // Windows GetTime
    QueryPerformanceCounter( &performance_counter );
    return (double) performance_counter.QuadPart * one_over_frequency;
#else
    // UNIX GetTime
    gettimeofday( &timeofday, NULL );
    return timeofday.tv_sec + timeofday.tv_usec / 1000000.0;
#endif
}



void Timer::Start()
{
    Assert( !running );
    running = 1;
    time0 = GetTime();
}



void Timer::Stop()
{
    Assert( running );
    running = 0;

    elapsed += GetTime() - time0;
}



void Timer::Reset()
{
    running = 0;
    elapsed = 0;
}



double Timer::Time()
{
    if (running) {
        Stop();
        Start();
    }
    return elapsed;
}


