
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


// core/probes.cpp*
#include "stdafx.h"
#include "probes.h"
#ifdef PBRT_PROBES_COUNTERS
#include "parallel.h"
#include <map>
using std::map;

// Statistics Counters Local Declarations
#ifdef PBRT_HAS_64_BIT_ATOMICS
typedef AtomicInt64 StatsCounterType;
#else
typedef AtomicInt32 StatsCounterType;
#endif
class StatsCounter {
public:
    // StatsCounter Public Methods
    StatsCounter(const string &category, const string &name);
    void operator++() {
        AtomicAdd(&num, 1);
    }
    void operator++(int) {
        AtomicAdd(&num, 1);
    }
#ifdef PBRT_HAS_64_BIT_ATOMICS
    void Max(int64_t newval) {
        int64_t oldval;
        do {
            oldval = num;
            newval = max(oldval, newval);
            if (newval == oldval) return;
        } while(AtomicCompareAndSwap(&num, newval, oldval) != oldval);
    }
    void Min(int64_t newval) {
        int64_t oldval;
        do {
            oldval = (uint32_t)num;
            newval = min(oldval, newval);
            if (newval == oldval) return;
        } while(AtomicCompareAndSwap(&num, newval, oldval) != oldval);
    }
#else
    void Max(int32_t newval) {
        int32_t oldval;
        do {
            oldval = num;
            newval = max(oldval, newval);
            if (newval == oldval) return;
        } while(AtomicCompareAndSwap(&num, newval, oldval) != oldval);
    }
    void Min(int32_t newval) {
        int32_t oldval;
        do {
            oldval = (uint32_t)num;
            newval = min(oldval, newval);
            if (newval == oldval) return;
        } while(AtomicCompareAndSwap(&num, newval, oldval) != oldval);
    }
#endif
#ifdef PBRT_HAS_64_BIT_ATOMICS
    operator int64_t() volatile { return (int64_t)num; }
#else
    operator int32_t() volatile { return (int32_t)num; }
#endif
private:
    // StatsCounter Private Data
    StatsCounterType num;
};


class StatsRatio {
public:
    // StatsRatio Public Methods
    StatsRatio(const string &category, const string &name);
    void Add(int a, int b) {
        AtomicAdd(&na, a);
        AtomicAdd(&nb, b);
    }
private:
    // StatsRatio Private Data
    StatsCounterType na, nb;
};


class StatsPercentage {
public:
    // StatsPercentage Public Methods
    void Add(int a, int b) {
        AtomicAdd(&na, a);
        AtomicAdd(&nb, b);
    }
    StatsPercentage(const string &category, const string &name);
private:
    // StatsPercentage Private Data
    StatsCounterType na, nb;
};



// Statistics Counters Definitions
static void ProbesPrintVal(FILE *f, const StatsCounterType &v);
static void ProbesPrintVal(FILE *f, const StatsCounterType &v1,
    const StatsCounterType &v2);
struct StatTracker {
    StatTracker(const string &cat, const string &n,
                StatsCounterType *pa, StatsCounterType *pb = NULL,
            bool percentage = true);
    string category, name;
    StatsCounterType *ptra, *ptrb;
    bool percentage;
};


typedef map<std::pair<string, string>, StatTracker *> TrackerMap;
static TrackerMap trackers;
static void addTracker(StatTracker *newTracker) {
    static Mutex *mutex = Mutex::Create();
    MutexLock lock(*mutex);
    std::pair<string, string> s = std::make_pair(newTracker->category, newTracker->name);
    if (trackers.find(s) != trackers.end()) {
        newTracker->ptra = trackers[s]->ptra;
        newTracker->ptrb = trackers[s]->ptrb;
        return;
    }
    trackers[s] = newTracker;
}



// Statistics Counters Function Definitions
void ProbesPrint(FILE *dest) {
    fprintf(dest, "Statistics:\n");
    TrackerMap::iterator iter = trackers.begin();
    string lastCategory;
    while (iter != trackers.end()) {
        // Print statistic
        StatTracker *tr = iter->second;
        if (tr->category != lastCategory) {
            fprintf(dest, "%s\n", tr->category.c_str());
            lastCategory = tr->category;
        }
        fprintf(dest, "    %s", tr->name.c_str());

        // Pad out to results column
        int resultsColumn = 56;
        int paddingSpaces = resultsColumn - (int) tr->name.size();
        while (paddingSpaces-- > 0)
            putc(' ', dest);
        if (tr->ptrb == NULL)
            ProbesPrintVal(dest, *tr->ptra);
        else {
            if ((int64_t)(*tr->ptrb) > 0) {
                float ratio = double((int64_t)*tr->ptra) / double((int64_t)*tr->ptrb);
                ProbesPrintVal(dest, *tr->ptra, *tr->ptrb);
                if (tr->percentage)
                    fprintf(dest, " (%3.2f%%)", 100. * ratio);
                else
                    fprintf(dest, " (%.2fx)", ratio);
            }
            else
                ProbesPrintVal(dest, *tr->ptra, *tr->ptrb);
        }
        fprintf(dest, "\n");
        ++iter;
    }
}


static void ProbesPrintVal(FILE *f, const StatsCounterType &v) {
#if defined(PBRT_IS_WINDOWS)
    LONG vv = v;
    fprintf(f, "%d", vv);
#else
#ifdef PBRT_HAS_64_BIT_ATOMICS
    int64_t vv = v;
    fprintf(f, "%lld", vv);
#else
    int32_t vv = v;
    fprintf(f, "%d", vv);
#endif
#endif
}


static void ProbesPrintVal(FILE *f, const StatsCounterType &v1,
        const StatsCounterType &v2) {
#if defined(PBRT_IS_WINDOWS)
    LONG vv1 = v1, vv2 = v2;
    fprintf(f, "%d:%d", vv1, vv2);
#else
#ifdef PBRT_HAS_64_BIT_ATOMICS
    int64_t vv1 = v1, vv2 = v2;
    fprintf(f, "%lld:%lld", vv1, vv2);
#else
    int32_t vv1 = v1, vv2 = v2;
    fprintf(f, "%d:%d", vv1, vv2);
#endif
#endif
}


void ProbesCleanup() {
    TrackerMap::iterator iter = trackers.begin();
    string lastCategory;
    while (iter != trackers.end()) {
        delete iter->second;
        ++iter;
    }
    trackers.erase(trackers.begin(), trackers.end());
}


StatTracker::StatTracker(const string &cat, const string &n,
                         StatsCounterType *pa, StatsCounterType *pb, bool p) {
    category = cat;
    name = n;
    ptra = pa;
    ptrb = pb;
    percentage = p;
}


StatsCounter::StatsCounter(const string &category, const string &name) {
    num = 0;
    addTracker(new StatTracker(category, name, &num));
}


StatsRatio::StatsRatio(const string &category, const string &name) {
    na = nb = 0;
    addTracker(new StatTracker(category, name, &na, &nb, false));
}


StatsPercentage::StatsPercentage(const string &category, const string &name) {
    addTracker(new StatTracker(category, name, &na, &nb, true));
}



// Statistics Counters Probe Declarations
static StatsCounter shapesMade("Shapes", "Total Shapes Created");
static StatsCounter trianglesMade("Shapes", "Total Triangles Created");
static StatsCounter cameraRays("Rays", "Camera Rays Traced");
static StatsCounter specularReflectionRays("Rays", "Specular Reflection Rays Traced");
static StatsCounter specularRefractionRays("Rays", "Specular Refraction Rays Traced");
static StatsCounter shadowRays("Rays", "Shadow Rays Traced");
static StatsCounter nonShadowRays("Rays", "Total Non-Shadow Rays Traced");
static StatsCounter kdTreeInteriorNodes("Kd-Tree", "Interior Nodes Created");
static StatsCounter kdTreeLeafNodes("Kd-Tree", "Interior Nodes Created");
static StatsCounter kdTreeMaxPrims("Kd-Tree", "Maximum Primitives in Leaf");
static StatsCounter kdTreeMaxDepth("Kd-Tree", "Maximum Depth of Leaf Nodes");
static StatsPercentage rayTriIntersections("Intersections", "Ray/Triangle Intersection Hits");
static StatsPercentage rayTriIntersectionPs("Intersections", "Ray/Triangle IntersectionP Hits");

// Statistics Counters Probe Definitions
void PBRT_CREATED_SHAPE(Shape *) {
    ++shapesMade;
}


void PBRT_CREATED_TRIANGLE(Triangle *) {
    ++trianglesMade;
}


void PBRT_STARTED_GENERATING_CAMERA_RAY(const CameraSample *) {
    ++cameraRays;
}



void PBRT_KDTREE_CREATED_INTERIOR_NODE(int axis, float split) {
    ++kdTreeInteriorNodes;
}



void PBRT_KDTREE_CREATED_LEAF(int nprims, int depth) {
    ++kdTreeLeafNodes;
    kdTreeMaxPrims.Max(nprims);
    kdTreeMaxDepth.Max(depth);
}



void PBRT_RAY_TRIANGLE_INTERSECTION_TEST(const Ray *, const Triangle *) {
    rayTriIntersections.Add(0, 1);
}



void PBRT_RAY_TRIANGLE_INTERSECTIONP_TEST(const Ray *, const Triangle *) {
    rayTriIntersectionPs.Add(0, 1);
}



void PBRT_RAY_TRIANGLE_INTERSECTION_HIT(const Ray *, float t) {
    rayTriIntersections.Add(1, 0);
}



void PBRT_RAY_TRIANGLE_INTERSECTIONP_HIT(const Ray *, float t) {
    rayTriIntersectionPs.Add(1, 0);
}



void PBRT_FINISHED_RAY_INTERSECTION(const Ray *, const Intersection *, int hit) {
    ++nonShadowRays;
}



void PBRT_FINISHED_RAY_INTERSECTIONP(const Ray *, int hit) {
    ++shadowRays;
}



void PBRT_STARTED_SPECULAR_REFLECTION_RAY(const RayDifferential *) {
    ++specularReflectionRays;
}



void PBRT_STARTED_SPECULAR_REFRACTION_RAY(const RayDifferential *) {
    ++specularRefractionRays;
}


#endif // PBRT_PROBES_COUNTERS
