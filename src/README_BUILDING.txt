
pbrt is designed to be easily ported to various platforms; the authors
regularly compile it on Mac OS X, Windows, and a variety of Linux variants.
pbrt users have sent patches to ensure that it compiles cleanly on FreeBSD,
OpenBSD, and other systems.  We will happily incorporate patches to make
the system build on other platforms!  (Please send patches or other notes
about issues with building pbrt to authors@pbrt.org.)

=== Building The System ===

--- Windows ---

Please see the file README_BUILDING_MSVC2008.txt or
README_BUILDING_MSVC2010.txt, as appropriate for the version of Visual
Studio you're using.

There is not currently support for building pbrt with the Cygwin gcc
compiler, but patches would be happily accepted.

--- Linux ---

On Linux and other Unix platforms, pbrt can be compiled with either the
provided Makefile or the provided SCons build files (see http://scons.org
for information about SCons).  Please see the notes below about installing
OpenEXR libraries on your system before building pbrt. 

The SCons build files build both debug and release configurations of the
system, while the Makefile only builds a release build.  See comment at the
top of the Makefile for how to modify it to do a debug build instead.

--- Mac OS X ---

Please see the notes below about installing OpenEXR libraries on your
system before building pbrt.

In addition to the Makefile and SCons files described in the "Linux"
section, there is also is also an XCode project file for Mac OS X,
pbrt.xcodeproj.


=== Build Options and Configuration ===

The remainder of this document has notes about the two main build
configuration options: using the OpenEXR and TIFF image formats or not, and
how pbrt collects and reports statistics.  It then has sections about each
of the main target platforms (Linux, Mac OS X, and Windows).

--- OpenEXR ---

If you have the OpenEXR image library installed (see http://openexr.com),
then pbrt will read and write OpenEXR format images.  A number of the
example scenes use OpenEXR image files.  If you do have OpenEXR installed,
then PBRT_HAS_OPENEXR should be #defined and the paths to the OpenEXR
headers and libraries should be set in the build rules as appropriate.

On Mac OS X and Linux, OpenEXR compiles easily from the distribution from
the OpenEXR website.  Alternatively, most package or "ports" systems
provide an OpenEXR installation.

--- TIFF ---

We provide a utility program to convert from high dynamic range EXR images
to low dynamic range TIFF images, exrtotiff.  This program includes a
rudimentary tone mapping pipeline and support for image bloom and gamma
correction.  To build this program, modify the user configuration section
appropriately in either Makefile or SCons file (depending on how you're
building the system.)  There is not currently any support for building this
from the MSVC solution file on Windows.

More comprehensive sets of programs to work with EXR images are available
from http://scanline.ca/exrtools/ and http://pfstools.sourceforge.net/.

--- Probes and Statistics ---

pbrt no longer collects runtime rendering statistics by default.  (Updating
shared statistics counters can cause a substantial performance impact with
multi-threaded execution.)  To override this, change the definition of the
PBRT_PROBES_NONE preprocessor #define to either PBRT_PROBES_DTRACE or
PBRT_PROBES_COUNTERS.

If your system supports dtrace (OSX, FreeBSD, and Solaris, currently--see
http://en.wikipedia.org/wiki/DTrace for more information), then pbrt is
instrumented to provide a large number of dtrace probes that can be
analyzed with dtrace scripts; building with dtrace support is supported if
you set the PBRT_PROBES_DTRACE preprocessor #define.  See the dtrace/
directory for a number of example dtrace scripts to use with pbrt.

Alternatively, PBRT_PROBES_COUNTERS can be set to compile the system to
gather a number of statistics with shared counters, incurring the
corresponding performance penalty.
