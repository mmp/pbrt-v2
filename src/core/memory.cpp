
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


