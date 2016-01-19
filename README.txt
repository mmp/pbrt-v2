
Welcome to the 'final' version of pbrt, version 2.  This version of the
source code corresponds to the system that is described in the second
edition of the book "Physically Based Rendering: From Theory To
Implementation", released in July 2010.

The latest version of the pbrt code, including all of the latest bug fixes,
updates, and conversion utilities is available from
http://github.com/mmp/pbrt-v2.  Please report any bugs encountered or other
issues with the system to Matt Pharr and Greg Humphreys, via the
authors@pbrt.org e-mail address.

pbrt-v2 has been ported to a number of architectures and operating systems,
including Windows (x86 and x64), Mac OS X, Linux, and OpenBSD.  It should
build on most UNIX-like systems.  Please see the file
src/README_BUILDING.txt for more information about compiling the system.

--- Organization ---

src/ : The implementation of the pbrt rendering system is in this
directory.  It includes a MSVC project to build the system as well as an
XCode project, build scripts for scons, and a Makefile as well.

scenes/ : A number of simple example scenes.

exporters/ : Scripts to export to the pbrt file format from a number of
modeling systems.  (Currently, 3ds max, Blender, StructureSynth, as well as
a Mathematica exporter)

dtrace/ : A number of scripts for gathering data about the run-time
behavior of pbrt using dtrace (http://en.wikipedia.org/wiki/DTrace).
Dtrace is only supported on Mac OSX and FreeBSD; see the
src/README_BUILDING.txt file for more information about building pbrt with
dtrace support.

--- Changes ---

The remainder of this document will summarize the major changes to the
system since the version described in the first edition of the book.  See
the file src/README_BUILDING.txt for information about how to compile the
system.

--- Incompatibilities ---

Many existing pbrt scene files will work unmodified with the second version
of the system.  The most significant user-visible change is that we have
taken this opportunity to fix the long-standing bug where LookAt
inadvertently flipped the handed-ness of the coordinate system.
You may find that scene files
that use LookAt now render flipped images or otherwise have problems with
the camera positioning.  Add "Scale -1 1 1" to the top of any scene
description files with this problem to return to the previous camera
position.

pbrt no longer has a plugin implementation based on runtime loading of DSOs
or DLLs.  All implementations of shapes, lights, etc., are just statically
compiled into the pbrt binary.  (We found that the complexity of supporting
this, both in cross-platform headaches, user issues with plugins not being
found, and so forth wasn't worth the limited benefits.)  The
PBRT_SEARCHPATH environment variable is no longer used, and the SearchPath
directive has been removed from the input file syntax.

We have removed the "shinymetal" material, replacing it with a physically
based "metal" material.  (Described further in the new features section
below.)

--- New Features ---

-- General --

pbrt is now multithreaded; performance should scale well with increasing
numbers of CPU cores.  We are particularly interested to hear any
experience running the system with >4 cores as well as any insight gained
from digging into any scalability bottlenecks.  The system attempts to
automatically determine how many CPU cores are present in the system, but
the --ncores command line argument can be used to override this.

OpenEXR is no longer required to build the system (but it is highly
recommended).  pbrt now includes code to read and write both TGA and PFM
format files; support for those file format is thus always available.  If
pbrt is compiled with PBRT_HAS_OPENEXR #defined, then OpenEXR files can be
used as well.

An accelerator based on bounding volume hierarchies has been added
(accelerators/bvh*).  This accelerator is now the default; it spends
substantially less time than the kd-tree accelerator in hierarchy
construction, while providing nearly equal performance. 

pbrt now supports full spectral rendering as a compile-time option; to
enable it, change the "typedef RGBSpectrum Spectrum" in core/pbrt.h to
"typedef SampledSpectrum Spectrum".  The number of spectral samples taken
(30 by default, leading to 10nm spacing), can be changed in the
core/spectrum.h file.

Animation is now supported via animated transformations for cameras and
shapes in the scene.  (But not light sources, however.)  See the included
example files scenes/anim-killeroos-moving.pbrt and
scenes/anim-moving-reflection.pbrt.

A rudimentary adaptive sampler is included, see samplers/adaptive.*.

-- Integrators --

The 'instant global illumination', 'extended photon map', and 'extended
infinite area light source' implementations from the author-supplied
plugins for pbrt-v1 are now part of the standard
distribution; the previous photon map and infinite area light source
implementations described in the first version of the book have been
removed.

The irradiance cache implementation has been substantially improved,
following many of the ideas from the Tabellion and Lamorlette paper from
SIGGRAPH 2004.

A subsurface scattering integrator based on the Jensen and Buhler 2002
SIGGRAPH paper has been added (integrators/dipolesubsurface.*).  Two new
materials for translucent materials have been added, "kdsubsurface" and
"subsurface".  The former implements Jensen and Buhler's intuitive approach
for setting subsurface scattering parameters.  The latter also incorporates
measured scattering parameters from a number of recent papers; many
translucent material types are available by name.  (See the source code for
details.)

An implementation of Metropolis light transport has been added (based on
the Kelemen et al 2002 Eurographics paper).  This is implemented as a new
"renderer", rather than for example a surface integrator.  See the included
example scene file scenes/metal.pbrt.

Support for a variety of forms of precomputed radiance transport has been
added.  A new renderer that computes radiance light probes in spherical
harmonics on a grid has been added (renderers/createprobes.*), and a
surface integrator that uses them is available in integrators/useprobes.*.
The diffuse precomputed radiance transfer (PRT) method of Sloan et al's 2002
SIGGRAPH paper is implemented in integrators/diffuseprt.* and a technique
for glossy PRT is in integrators/glossyprt.*.  Note that in general, one
would use pbrt to do the precomputation part of PRT and then use the
results in a real-time rendering system.  Here we have also implemented the
code that uses the results within pbrt for pedagogical purposes (and so
that we don't need to include a real-time renderer with the book!).

-- Materials --

A new 'metal' material has been added; it supports setting the spectral
index of refraction and extinction coefficients via measured data from real
metals.  See the example file scenes/metal.pbrt.  A large number of
spectral measurements are available in the scenes/spds/metals/ directory.

Support for measured BRDF data is now available via the "measured"
material.  Two BRDF representations and file formats are supported.  First
is the binary file format used by the MERL BRDF database.
(http://merl.com/brdf).  For arbitrary scattered BRDF measurements, a
general text file format is supported; see the comments in
scenes/brdfs/acryl_blue.brdf.

--- Example Scenes ---

A small number of example scenes that demonstrate some of the new features
are provided in the scenes/ directory.

anim-killeroos-moving.pbrt, anim-moving-reflection.pbrt: demonstrates
motion blur features.

bunny.pbrt: measured BRDFs

killeroo-simple.pbrt: simple scene with "Killeroo" model

metal.pbrt: Metropolis light transport, measured BRDFs

prt-teapot.pbrt: precomputed radiance transfer

ss-envmap.pbrt: subsurface scattering
