
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


// core/progressreporter.cpp*
#include "stdafx.h"
#include "progressreporter.h"
#include "timer.h"
#include "parallel.h"
#if defined(PBRT_IS_WINDOWS)
#include <windows.h>
#else
#include <sys/ioctl.h>
#include <unistd.h>
#include <errno.h>
#endif // !PBRT_IS_WINDOWS

// ProgressReporter Method Definitions
ProgressReporter::ProgressReporter(int tw, const string &title, int barLength)
    : totalWork(tw) {
    if (barLength <= 0)
        barLength = TerminalWidth() - 28;
    totalPlusses = max(2, barLength - (int)title.size());
    mutex = Mutex::Create();
    plussesPrinted = 0;
    workDone = 0;
    timer = new Timer;
    timer->Start();
    outFile = stdout;
    // Initialize progress string
    const int bufLen = title.size() + totalPlusses + 64;
    buf = new char[bufLen];
    snprintf(buf, bufLen, "\r%s: [", title.c_str());
    curSpace = buf + strlen(buf);
    char *s = curSpace;
    for (int i = 0; i < totalPlusses; ++i)
        *s++ = ' ';
    *s++ = ']';
    *s++ = ' ';
    *s++ = '\0';
    if (!PbrtOptions.quiet) {
        fputs(buf, outFile);
        fflush(outFile);
    }
}


ProgressReporter::~ProgressReporter() {
    delete[] buf;
    delete timer;
    Mutex::Destroy(mutex);
}


void ProgressReporter::Update(int num) {
    if (num == 0 || PbrtOptions.quiet) return;
    MutexLock lock(*mutex);
    workDone += num;
    float percentDone = float(workDone) / float(totalWork);
    int plussesNeeded = Round2Int(totalPlusses * percentDone);
    if (plussesNeeded > totalPlusses) plussesNeeded = totalPlusses;
    while (plussesPrinted < plussesNeeded) {
        *curSpace++ = '+';
        ++plussesPrinted;
    }
    fputs(buf, outFile);
    // Update elapsed time and estimated time to completion
    float seconds = (float)timer->Time();
    float estRemaining = seconds / percentDone - seconds;
    if (percentDone == 1.f)
        fprintf(outFile, " (%.1fs)       ", seconds);
    else
        fprintf(outFile, " (%.1fs|%.1fs)  ", seconds, max(0.f, estRemaining));
    fflush(outFile);
}


void ProgressReporter::Done() {
    if (PbrtOptions.quiet) return;
    MutexLock lock(*mutex);
    while (plussesPrinted++ < totalPlusses)
        *curSpace++ = '+';
    fputs(buf, outFile);
    float seconds = (float)timer->Time();
    fprintf(outFile, " (%.1fs)       \n", seconds);
    fflush(outFile);
}


int TerminalWidth() {
#if defined(PBRT_IS_WINDOWS)
    HANDLE h = GetStdHandle(STD_OUTPUT_HANDLE);
    if (h == INVALID_HANDLE_VALUE || h == NULL) {
        fprintf(stderr, "GetStdHandle() call failed");
        return 80;
    }
    CONSOLE_SCREEN_BUFFER_INFO bufferInfo = { 0 };
    GetConsoleScreenBufferInfo(h, &bufferInfo);
    return bufferInfo.dwSize.X;
#else
    struct winsize w;
    if (ioctl(STDOUT_FILENO, TIOCGWINSZ, &w) < 0) {
        fprintf(stderr, "Error in ioctl() in TerminalWidth(): %d", errno);
        return 80;
    }
    return w.ws_col;
#endif // PBRT_IS_WINDOWS
}



