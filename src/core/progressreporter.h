
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

#ifndef PBRT_CORE_PROGRESSREPORTER_H
#define PBRT_CORE_PROGRESSREPORTER_H

// core/progressreporter.h*
#include "pbrt.h"

// ProgressReporter Declarations
class ProgressReporter {
public:
    // ProgressReporter Public Methods
    ProgressReporter(int totalWork, const string &title,
                     int barLength = -1);
    ~ProgressReporter();
    void Update(int num = 1);
    void Done();
private:
    // ProgressReporter Private Data
    const int totalWork;
    int workDone, plussesPrinted, totalPlusses;
    Timer *timer;
    FILE *outFile;
    char *buf, *curSpace;
    Mutex *mutex;
};


extern int TerminalWidth();

#endif // PBRT_CORE_PROGRESSREPORTER_H
