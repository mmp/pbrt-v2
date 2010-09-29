
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


// core/memory.cpp*
#include "stdafx.h"
#include "memory.h"

// Memory Allocation Functions
void *AllocAligned(size_t size) {
#if defined(PBRT_IS_WINDOWS)
    return _aligned_malloc(size, PBRT_L1_CACHE_LINE_SIZE);
#elif defined (PBRT_IS_OPENBSD) || defined(PBRT_IS_APPLE)
    // Allocate excess memory to ensure an aligned pointer can be returned
    void *mem = malloc(size + (PBRT_L1_CACHE_LINE_SIZE-1) + sizeof(void*));
    char *amem = ((char*)mem) + sizeof(void*);
#if (PBRT_POINTER_SIZE == 8)
    amem += PBRT_L1_CACHE_LINE_SIZE - (reinterpret_cast<uint64_t>(amem) &
                                       (PBRT_L1_CACHE_LINE_SIZE - 1));
#else
    amem += PBRT_L1_CACHE_LINE_SIZE - (reinterpret_cast<uint32_t>(amem) &
                                       (PBRT_L1_CACHE_LINE_SIZE - 1));
#endif
    ((void**)amem)[-1] = mem;
    return amem;
#else
    return memalign(PBRT_L1_CACHE_LINE_SIZE, size);
#endif
}


void FreeAligned(void *ptr) {
    if (!ptr) return;
#if defined(PBRT_IS_WINDOWS)
    _aligned_free(ptr);
#elif defined (PBRT_IS_OPENBSD) || defined(PBRT_IS_APPLE)
    free(((void**)ptr)[-1]);
#else
    free(ptr);
#endif
}


