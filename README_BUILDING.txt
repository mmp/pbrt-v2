
On Linux/Mac/other Unix platforms, pbrt can be compiled with either the
provided Makefile or the provided scons build files.  There is also an
experimental XCode project file for Mac OS X.

Under Windows, use the provided MSVC solution file. On Windows, may also
need to install Cygwin with the flex and bison packages.

A small number of parameters may need to be set in the Makefile.  If your
system has OpenEXR installed (this is recommended), PBRT_HAS_OPENEXR should
be #defined and the paths to the OpenEXR headers and libraries should be
set in the build rules as appropriate.

For a more optimized build, you may wish to set the NDEBUG #define; this is
done for 'release' builds in the Windows MSVS, Mac XCode, and scons build
files automatically.  Doing so will disable the assertions.  However, it's
helpful to have these enabled at this point in the system's development.

Under OSX version 10.6, you may want pbrt to use Grand Central Dispatch for
managing the parallel rendering tasks.  Add a definition of
PBRT_USE_GRAND_CENTRAL_DISPATCH to do so.

Finally, this version of pbrt doesn't collect any runtime rendering
statistics by default.  (Updating shared statistics counters can cause a
substantial performance impact with multi-threaded execution.)  To override
this, change the definition of PBRT_PROBES_NONE to either
PBRT_PROBES_DTRACE or PBRT_PROBES_COUNTERS.  If your system supports dtrace
(OSX, FreeBSD, and Solaris, currently), then pbrt is instrumented to
provide a large number of dtrace probes that can be analyzed with dtrace
scripts; building with dtrace support is supported if you set
PBRT_PROBES_DTRACE.  (However, you need to manually compile the dtrace.d
file into a header file by running:

% /usr/sbin/dtrace -h -s core/dtrace.h -o core/dtrace.h

(The Makefile currently doesn't have a rule for this, but the scons and
XCode builds do.)  Finally, PBRT_PROBES_COUNTERS can be set to compile the
system to gather a number of statistics with shared counters, incurring the
corresponding performance penalty.
