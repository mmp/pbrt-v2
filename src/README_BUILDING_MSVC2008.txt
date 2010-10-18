Building PBRT with MSVC 2008

First install the "flex" and "bison" packages from the Cygwin set of
tools (http://cygwin.com).  You need to explicitly select
"devel/bison" and "devel/flex" to be installed by the Cygwin
installer; they aren't installed by default.

The provided MSVC solution file for pbrt, src/pbrt.sln, assumes that
the Cygwin installation is c:\cygwin; if it's not installed there,
you'll need to modify their custom build rules.  Right click on "pbrt"
in the Solution Explorer and choose "Custom Build Rules".  Click on
"Bison/Flex" in the list of available rule files and choose "Modify
Rule File...".  Finally, select each of the two custom build rules in
turn and choose "Modify Build Rule...", changing the "Command Line"
property for each one to point to the actual locations of bison.exe
and flex.exe.

Next, compile the OpenEXR libraries; their source code as well as a
custom solution file to build them is provided in
src/windows/3rdparty.vs2008/3rdparty.sln.  Build whichever of the
Debug/Release, x86/x64 variants you need.  The built libraries will be
automatically stored in the directories where the main pbrt build
rules will look for them.

The system should then compile cleanly from the provided MSVC solution
file, src/pbrt.sln.  The solution file supports both 32-bit and 64-bit
builds, with both Debug and Release configurations.
