
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

#ifndef PBRT_CORE_MEMORY_H
#define PBRT_CORE_MEMORY_H

// core/memory.h*
#include "pbrt.h"
#include "parallel.h"

// Memory Declarations
class ReferenceCounted {
public:
    ReferenceCounted() { nReferences = 0; }
    AtomicInt32 nReferences;
private:
    ReferenceCounted(const ReferenceCounted &);
    ReferenceCounted &operator=(const ReferenceCounted &);
};


template <typename T> class Reference {
public:
    // Reference Public Methods
    Reference(T *p = NULL) {
        ptr = p;
        if (ptr) AtomicAdd(&ptr->nReferences, 1);
    }
    Reference(const Reference<T> &r) {
        ptr = r.ptr;
        if (ptr) AtomicAdd(&ptr->nReferences, 1);
    }
    Reference &operator=(const Reference<T> &r) {
        if (r.ptr) AtomicAdd(&r.ptr->nReferences, 1);
        if (ptr && AtomicAdd(&ptr->nReferences, -1) == 0) delete ptr;
        ptr = r.ptr;
        return *this;
    }
    Reference &operator=(T *p) {
        if (p) AtomicAdd(&p->nReferences, 1);
        if (ptr && AtomicAdd(&ptr->nReferences, -1) == 0) delete ptr;
        ptr = p;
        return *this;
    }
    ~Reference() {
        if (ptr && AtomicAdd(&ptr->nReferences, -1) == 0)
            delete ptr;
    }
    T *operator->() { return ptr; }
    const T *operator->() const { return ptr; }
    operator bool() const { return ptr != NULL; }
    const T *GetPtr() const { return ptr; }
private:
    T *ptr;
};


void *AllocAligned(size_t size);
template <typename T> T *AllocAligned(uint32_t count) {
    return (T *)AllocAligned(count * sizeof(T));
}


void FreeAligned(void *);
class MemoryArena {
public:
    // MemoryArena Public Methods
    MemoryArena(uint32_t bs = 32768) {
        blockSize = bs;
        curBlockPos = 0;
        currentBlock = AllocAligned<char>(blockSize);
    }
    ~MemoryArena() {
        FreeAligned(currentBlock);
        for (uint32_t i = 0; i < usedBlocks.size(); ++i)
            FreeAligned(usedBlocks[i]);
        for (uint32_t i = 0; i < availableBlocks.size(); ++i)
            FreeAligned(availableBlocks[i]);
    }
    void *Alloc(uint32_t sz) {
        // Round up _sz_ to minimum machine alignment
        sz = ((sz + 15) & (~15));
        if (curBlockPos + sz > blockSize) {
            // Get new block of memory for _MemoryArena_
            usedBlocks.push_back(currentBlock);
            if (availableBlocks.size() && sz <= blockSize) {
                currentBlock = availableBlocks.back();
                availableBlocks.pop_back();
            }
            else
                currentBlock = AllocAligned<char>(max(sz, blockSize));
            curBlockPos = 0;
        }
        void *ret = currentBlock + curBlockPos;
        curBlockPos += sz;
        return ret;
    }
    template<typename T> T *Alloc(uint32_t count = 1) {
        T *ret = (T *)Alloc(count * sizeof(T));
        for (uint32_t i = 0; i < count; ++i)
            new (&ret[i]) T();
        return ret;
    }
    void FreeAll() {
        curBlockPos = 0;
        while (usedBlocks.size()) {
    #ifndef NDEBUG
            memset(usedBlocks.back(), 0xfa, blockSize);
    #endif
            availableBlocks.push_back(usedBlocks.back());
            usedBlocks.pop_back();
        }
    }
private:
    // MemoryArena Private Data
    uint32_t curBlockPos, blockSize;
    char *currentBlock;
    vector<char *> usedBlocks, availableBlocks;
};


template <typename T, int logBlockSize> class BlockedArray {
public:
    // BlockedArray Public Methods
    BlockedArray(uint32_t nu, uint32_t nv, const T *d = NULL) {
        uRes = nu;
        vRes = nv;
        uBlocks = RoundUp(uRes) >> logBlockSize;
        uint32_t nAlloc = RoundUp(uRes) * RoundUp(vRes);
        data = AllocAligned<T>(nAlloc);
        for (uint32_t i = 0; i < nAlloc; ++i)
            new (&data[i]) T();
        if (d)
            for (uint32_t v = 0; v < vRes; ++v)
                for (uint32_t u = 0; u < uRes; ++u)
                    (*this)(u, v) = d[v * uRes + u];
    }
    uint32_t BlockSize() const { return 1 << logBlockSize; }
    uint32_t RoundUp(uint32_t x) const {
        return (x + BlockSize() - 1) & ~(BlockSize() - 1);
    }
    uint32_t uSize() const { return uRes; }
    uint32_t vSize() const { return vRes; }
    ~BlockedArray() {
        for (uint32_t i = 0; i < uRes * vRes; ++i)
            data[i].~T();
        FreeAligned(data);
    }
    uint32_t Block(uint32_t a) const { return a >> logBlockSize; }
    uint32_t Offset(uint32_t a) const { return (a & (BlockSize() - 1)); }
    T &operator()(uint32_t u, uint32_t v) {
        uint32_t bu = Block(u), bv = Block(v);
        uint32_t ou = Offset(u), ov = Offset(v);
        uint32_t offset = BlockSize() * BlockSize() * (uBlocks * bv + bu);
        offset += BlockSize() * ov + ou;
        return data[offset];
    }
    const T &operator()(uint32_t u, uint32_t v) const {
        uint32_t bu = Block(u), bv = Block(v);
        uint32_t ou = Offset(u), ov = Offset(v);
        uint32_t offset = BlockSize() * BlockSize() * (uBlocks * bv + bu);
        offset += BlockSize() * ov + ou;
        return data[offset];
    }
    void GetLinearArray(T *a) const {
        for (uint32_t v = 0; v < vRes; ++v)
            for (uint32_t u = 0; u < uRes; ++u)
                *a++ = (*this)(u, v);
    }
private:
    // BlockedArray Private Data
    T *data;
    uint32_t uRes, vRes, uBlocks;
};



#endif // PBRT_CORE_MEMORY_H
