
On Linux/Mac/other Unix platforms, pbrt can be compiled with the provided
Makefile.  Under Windows, use the provided MSVC solution file. You may also
need to install Cygwin with the flex and bison packages.

A small number of parameters may need to be set in the Makefile.  If your
system has OpenEXR installed (this is recommended), PBRT_HAS_OPENEXR should
be set in the DEFS line and the paths to the OpenEXR headers and libraries
in the EXRINCLUDE and EXRLIBDIR lines should be updated if needed.

The default build assumes you're compiling on a 32-bit system.  If this is
not the case, comment the 32-bit section in the Makefile and uncomment the
64-bit section.  For a more optimized build, you may wish to add -DNDEBUG
to the DEFS in the Makefile; this will disable the assertions.  However,
it's helpful to have these enabled at this point in the system's
development.

Under OSX version 10.6, you may want pbrt to use Grand Central Dispatch for
managing the parallel rendering tasks.  Add a definition of
PBRT_USE_GRAND_CENTRAL_DISPATCH in the DEFS line to do so.

Finally, this version of pbrt doesn't connect any runtime rendering
statistics.  (Updating shared statistics counters can cause a substantial
performance impact with multi-threaded execution.)  To override this,
change the definition of PBRT_STATS_NONE to either PBRT_STATS_DTRACE
or PBRT_STATS_COUNTERS.  If your system supports dtrace (OSX, FreeBSD, and
Solaris, currently), then pbrt is instrumented to provide a large number of
dtrace probes that can be analyzed with dtrace scripts; building with
dtrace support is supported if you set PBRT_STATS_DTRACE.  (However, you
need to manually compile the dtrace.d file into a header file by running:

% /usr/sbin/dtrace -h -s core/dtrace.h -o core/dtrace.h

The Makefile currently doesn't have a rule for this.)  Finally,
PBRT_STATS_COUNTERS can be set to compile the system to gather a number of
statistics with shared counters, incurring the corresponding performance
penalty.
