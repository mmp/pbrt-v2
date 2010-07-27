
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


