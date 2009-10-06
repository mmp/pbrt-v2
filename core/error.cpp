
/*
    pbrt source code Copyright(c) 1998-2009 Matt Pharr and Greg Humphreys.

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


// core/error.cpp*
#include "pbrt.h"

// Error Reporting Includes
#include <stdarg.h>

// Error Reporting Definitions
#define PBRT_ERROR_IGNORE 0
#define PBRT_ERROR_CONTINUE 1
#define PBRT_ERROR_ABORT 2

// Error Reporting Functions
static void processError(const char *format, va_list args,
        const char *message, int disposition) {
#ifndef WIN32
    char *errorBuf;
    if (vasprintf(&errorBuf, format, args) == -1) {
        fprintf(stderr, "vasprintf() unable to allocate memory!\n");
        abort();
    }
#else
    char errorBuf[2048];
    vsnprintf_s(errorBuf, sizeof(errorBuf), _TRUNCATE, format, args);
#endif
    // Report error
    switch (disposition) {
    case PBRT_ERROR_IGNORE:
        return;
    case PBRT_ERROR_CONTINUE:
        // Print scene file and line number, if appropriate
        extern int line_num;
        if (line_num != 0) {
            extern string current_file;
            fprintf(stderr, "%s(%d): ", current_file.c_str(), line_num);
        }
        fprintf(stderr, "%s: %s\n", message, errorBuf);
        break;
    case PBRT_ERROR_ABORT:
        // Print scene file and line number, if appropriate
        extern int line_num;
        if (line_num != 0) {
            extern string current_file;
            fprintf(stderr, "%s(%d): ", current_file.c_str(), line_num);
        }
        fprintf(stderr, "%s: %s\n", message, errorBuf);
#ifdef WIN32
        __debugbreak();
#else
        abort();
#endif
    }
#ifndef WIN32
    free(errorBuf);
#endif
}


void Info(const char *format, ...) {
    va_list args;
    va_start(args, format);
    processError(format, args, "Notice", PBRT_ERROR_CONTINUE);
    va_end(args);
}


void Warning(const char *format, ...) {
    va_list args;
    va_start(args, format);
    processError(format, args, "Warning", PBRT_ERROR_CONTINUE);
    va_end(args);
}


void Error(const char *format, ...) {
    va_list args;
    va_start(args, format);
    processError(format, args, "Error", PBRT_ERROR_CONTINUE);
    va_end(args);
}


void Severe(const char *format, ...) {
    va_list args;
    va_start(args, format);
    processError(format, args, "Fatal Error", PBRT_ERROR_ABORT);
    va_end(args);
}


