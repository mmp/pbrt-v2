
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


// core/parser.cpp*
#include "stdafx.h"
#include "parser.h"
#include "fileutil.h"

// Parsing Global Interface
bool ParseFile(const string &filename) {
    extern FILE *yyin;
    extern int yyparse(void);
    extern string current_file;
    extern int line_num;
    extern int yydebug;

    if (getenv("PBRT_YYDEBUG") != NULL)
        yydebug = 1;

    if (filename == "-")
        yyin = stdin;
    else {
        yyin = fopen(filename.c_str(), "r");
        SetSearchDirectory(DirectoryContaining(filename));
    }

    if (yyin != NULL) {
        current_file = filename;
        if (yyin == stdin) current_file = "<standard input>";
        line_num = 1;
        yyparse();
        if (yyin != stdin) fclose(yyin);
    }
    current_file = "";
    line_num = 0;
    return (yyin != NULL);
}


