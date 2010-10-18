BUILDING PBRT:

1) Open the src/pbrt.vs2010/pbrt.sln Visual Studio solution file.
2) Select either Debug or Release, and x86 (win32) or x64.
3) Build the solution.


RUNNING PBRT:

1) Run the desired executable or program in the top-level bin/ directory
   of the pbrt distribution (<pbrt-dist>/bin).

DEBUGGING PBRT:

1) Select the Debug configuration for either the x86 or x64 platform.
2) Start the program in Visual Studio with debugging.

CLEANING UP:

1) Run cleanup.bat in the top-level directory of the pbrt distribution
   (<pbrt-dist>/cleanup.bat).

MODIFYING THE PBRT PARSING CODE:


If you need to modify the pbrt parsing code, follow these steps first:

1) Install Bison and Flex programs if they are not already installed.
   Choose either Cygwin or GnuWin32; either one will work.
   
   IMPORTANT NOTE:
   Do not install under "Program Files" or "Program Files (x86)".
   The installation path should not contain any spaces.
   Use C:\cygwin or C:\gnuwin32, for example.

   Option A: GnuWin32 <http://gnuwin32.sourceforge.net/>

      Install Bison and Flex via setup packages.

   Option B: Cygwin <http://www.cygwin.com/>

      Install Cygwin with Bison and Flex packages.
      These may not be selected by default,
      so be sure to select them.


2) Add Bison and Flex executable paths to the system PATH
   so that they can be run from any directory or location.

   a) Right-click "My Computer" or "Computer".

   b) Click Properties.

   c) Click "Advanced system settings" or "System Properties".

   d) Click the Advanced tab.

   e) Click the "Environment Variables..." button.

   f) Add the path(s) of the Flex and Bison programs to the PATH
      variable, under the system variables group box.
      Use the semi-colon (;) character to separate path items.

      Example path to add: C:\gnuwin32\bin;

   g) Press the OK button.

   h) Press the OK button.

   i) Either log off and log back on the computer, or restart it
      for the updated PATH to take effect.