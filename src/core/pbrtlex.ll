
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

/* state used for include file stuff */
%{

#define YY_MAIN 0
#define YY_NEVER_INTERACTIVE 1

#include "pbrt.h"
#include "api.h"
#include "fileutil.h"

struct ParamArray;

#if defined(PBRT_IS_WINDOWS)
#pragma warning(disable:4244)
#pragma warning(disable:4065)
#pragma warning(disable:4018)
#pragma warning(disable:4996)
#endif
#include "pbrtparse.hpp"

struct IncludeInfo {
    string filename;
    YY_BUFFER_STATE bufState;
    int lineNum;
};


vector<IncludeInfo> includeStack;

extern int line_num;
int str_pos;

void add_string_char(char c) {
    yylval.string[str_pos++] = c;
    yylval.string[str_pos] = '\0';
}


void include_push(char *filename) {
    if (includeStack.size() > 32)
        Severe("Only 32 levels of nested Include allowed in scene files.");
    IncludeInfo ii;
    extern string current_file;
    ii.filename = current_file;
    ii.bufState = YY_CURRENT_BUFFER;
    ii.lineNum = line_num;
    includeStack.push_back(ii);

    current_file = AbsolutePath(ResolveFilename(filename));
    line_num = 1;

    yyin = fopen(current_file.c_str(), "r");
    if (!yyin)
        Severe("Unable to open included scene file \"%s\"", current_file.c_str());
    yy_switch_to_buffer(yy_create_buffer(yyin, YY_BUF_SIZE));
}



void include_pop() {
    extern int line_num;
    extern string current_file;
    fclose(yyin);
    yy_delete_buffer(YY_CURRENT_BUFFER);
    yy_switch_to_buffer(includeStack.back().bufState);
    current_file = includeStack.back().filename;
    line_num = includeStack.back().lineNum;
    includeStack.pop_back();
}


%}
%option nounput
WHITESPACE [ \t\r]+
NUMBER [-+]?([0-9]+|(([0-9]+\.[0-9]*)|(\.[0-9]+)))([eE][-+]?[0-9]+)?
IDENT [a-zA-Z_][a-zA-Z_0-9]*
%x STR COMMENT INCL INCL_FILE
%%

"#" { BEGIN COMMENT; }
<COMMENT>. /* eat it up */
<COMMENT>\n { line_num++; BEGIN INITIAL; }
Accelerator             { return ACCELERATOR; }
ActiveTransform         { return ACTIVETRANSFORM; }
All                     { return ALL; }
AreaLightSource         { return AREALIGHTSOURCE; }
AttributeBegin          { return ATTRIBUTEBEGIN; }
AttributeEnd            { return ATTRIBUTEEND; }
Camera                  { return CAMERA; }
ConcatTransform         { return CONCATTRANSFORM; }
CoordinateSystem        { return COORDINATESYSTEM; }
CoordSysTransform       { return COORDSYSTRANSFORM; }
EndTime                 { return ENDTIME; }
Film                    { return FILM; }
Identity                { return IDENTITY; }
Include                 { return INCLUDE; }
LightSource             { return LIGHTSOURCE; }
LookAt                  { return LOOKAT; }
MakeNamedMaterial       { return MAKENAMEDMATERIAL; }
Material                { return MATERIAL; }
NamedMaterial           { return NAMEDMATERIAL; }
ObjectBegin             { return OBJECTBEGIN; }
ObjectEnd               { return OBJECTEND; }
ObjectInstance          { return OBJECTINSTANCE; }
PixelFilter             { return PIXELFILTER; }
Renderer                { return RENDERER; }
ReverseOrientation      { return REVERSEORIENTATION; }
Rotate                  { return ROTATE; }
Sampler                 { return SAMPLER; }
Scale                   { return SCALE; }
Shape                   { return SHAPE; }
StartTime               { return STARTTIME; }
SurfaceIntegrator       { return SURFACEINTEGRATOR; }
Texture                 { return TEXTURE; }
TransformBegin          { return TRANSFORMBEGIN; }
TransformEnd            { return TRANSFORMEND; }
TransformTimes          { return TRANSFORMTIMES; }
Transform               { return TRANSFORM; }
Translate               { return TRANSLATE; }
Volume                  { return VOLUME; }
VolumeIntegrator        { return VOLUMEINTEGRATOR; }
WorldBegin              { return WORLDBEGIN; }
WorldEnd                { return WORLDEND; }
{WHITESPACE} /* do nothing */
\n { line_num++; }
{NUMBER} {
    yylval.num = (float) atof(yytext);
    return NUM;
}


{IDENT} {
    strcpy(yylval.string, yytext);
    return ID;
}


"[" { return LBRACK; }
"]" { return RBRACK; }
\" { BEGIN STR; str_pos = 0; }
<STR>\\n {add_string_char('\n');}
<STR>\\t {add_string_char('\t');}
<STR>\\r {add_string_char('\r');}
<STR>\\b {add_string_char('\b');}
<STR>\\f {add_string_char('\f');}
<STR>\\\" {add_string_char('\"');}
<STR>\\\\ {add_string_char('\\');}
<STR>\\[0-9]{3} {
  int val = atoi(yytext+1);
  while (val > 256)
    val -= 256;
  add_string_char(val);
}


<STR>\\\n {line_num++;}
<STR>\\. { add_string_char(yytext[1]);}
<STR>\" {BEGIN INITIAL; return STRING;}
<STR>. {add_string_char(yytext[0]);}
<STR>\n {Error("Unterminated string!");}

. { Error( "Illegal character: %c (0x%x)", yytext[0], int(yytext[0])); }
%%
int yywrap() {
    if (includeStack.size() == 0) return 1;
    include_pop();
    return 0;
}



