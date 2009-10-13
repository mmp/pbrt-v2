
Welcome to this early release of pbrt, version 2.  This version of the
source code corresponds to the system that will be described in the second
edition of the book "Physically Based Rendering: From Theory To
Implementation", which is expected to be available at SIGGRAPH 2010.

We are making the source code available before publication of the book in
the hopes that it will be useful to developers and users.  We in
particular hope that we will receive some feedback on the changes to the
system: comments on better implementation approaches, and performance or
portability shortcomings would be highly welcomed!  If you do have any
feedback, please either send it to the pbrt mailing list, or to
authors@pbrt.org.

We will release updated versions of the system frequently; the latest
version of the code will always be available at the github page
http://github.com/mmp/pbrt-v2.  The best way to track changes to the system
is to use the git source code management system, although tar files and zip
files are also available for download from the github page by clicking the
"download" button.

The remainder of this document will summarize the major changes to the
system.  See the file README_BUILDING.txt for information about how to
compile the system.

--- Incompatibilities ---

Many existing pbrt scene files will work unmodified with the second version
of the system.  The most significant user-visible change is that we have
taken this opportunity to fix the long-standing bug where LookAt
inadvertently flipped the handed-ness of the coordinate system
(http://pbrt.org/bugtracker/view.php?id=5).  You may find that scene files
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
the PBRT_NTHREADS environment variable can be used to override this.

OpenEXR is no longer required to build the system (but it is highly
recommended).  pbrt now includes code to read and write TGA format files;
support for that file format is thus always available.  If pbrt is compiled
with PBRT_HAS_OPENEXR #defined, then OpenEXR files can be used as well.
(With the file suffix determining which format is used.) 

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
scenes/anim-moving-reflection.pbrt.  It is likely that the method for
specifying animated transformations in input files will change in order to
be more similar to that used by lux render.

A rudimentary adaptive sampler is included, see samplers/adaptive.*.

-- Integrators --

The 'instant global illumination', 'extended photon map', and 'extended
infinite area light source' implementations from the author-supplied
plugins at http://pbrt.org/downloads.php are now standard; the previous
photon map and infinite area light source implementations described in the
first version of the book have been removed.

The irradiance cache implementation has been substantially improved,
following many of the ideas from the Tabellion and Lamorlette paper from
SIGGRAPH 2004.

A subsurface scattering integrator based on the Jensen and Buhler 2002
SIGGRAPH paper has been added (integrators/dipolesubsurface.*).  Two new
materials for translucent materials have been added, "kdsubsurface" and
"subsurface".  The former implements Jensen and Buhler's more intuitive
approach for setting subsurface scattering parameters.  The latter also
incorporates measured scattering parameters from a number of recent papers;
many translucent material types are available by name.  (See the source
code for details.)

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

metal.pbrt: Metropolis light transport, measured BRDFs

prt-teapot.pbrt: precomputed radiance transfer

ss-envmap.pbrt: subsurface scattering
