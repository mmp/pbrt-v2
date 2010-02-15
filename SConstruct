# -*- mode: python -*-

import sys, platform
arch = sys.platform

Decider('MD5-timestamp')
#CacheDir('scons-cache')

#import os
#print "Pruning scons cache..."
#os.system('cd scons-cache && find . -type f -atime +1 -delete')

use_sdl = False
Export('use_sdl')

######################################################################
## Configure generic environment

def setup_nice_print(env):
    if ARGUMENTS.get('VERBOSE') != '1':
        env['YACCCOMSTR'] = "Compiling $TARGET"
        env['LEXCOMSTR'] = "Compiling $TARGET"
        env['CCCOMSTR'] = "Compiling $TARGET"
        env['CXXCOMSTR'] = "Compiling $TARGET"
        env['LINKCOMSTR'] = "Linking $TARGET"
        env['ARCOMSTR'] = "Linking $TARGET"

build_envs = { }

env = Environment(CCFLAGS = [ '-Wall', '-g' ],
                  CPPPATH = [ '#core', '#', '/usr/local/include', '/opt/local/include', '/opt/local/include/OpenEXR',
                              '/usr/local/include/OpenEXR', '/usr/include/OpenEXR', '.' ],
                  LIBPATH = ['/opt/local/lib'],
                  YACCFLAGS = [ '-d', '-v', '-t' ],
                  YACCHXXFILESUFFIX = '.hh',
#                      CC='llvm-gcc-4.2',
#                      CXX='llvm-g++-4.2',
                  ENV = { 'PATH' : [ '/usr/local/bin', '/usr/bin', '/bin', 
#                                         '/Developer/usr/llvm-gcc-4.2/bin',
                                         '/usr/sbin', '/sbin' ] })

if arch == 'darwin':
#        has_gcd = float(platform.mac_ver()[0][:4] >= 10.6)
#        if (has_gcd):
#            env.Append(CPPDEFINES = [ 'PBRT_USE_GRAND_CENTRAL_DISPATCH' ])
    if (use_sdl):
        env.Append(CPPDEFINES = [ 'PBRT_HAS_LIBSDL' ],
                   LINKFLAGS = ['-F/System/Library/Frameworks', '-framework', 'Cocoa'],
                   CPPPATH = ['/opt/local/include/SDL'])
        sdl_main = [ 'main/SDLMain.m' ]
        Export('sdl_main')
    env.Append(BUILDERS = { 'DTrace' : 
                            Builder(action = '/usr/sbin/dtrace -h -s $SOURCE -o $TARGET')})

env.Append(CPPDEFINES = [ 'PBRT_HAS_PTHREADS' ])
env.Append(CPPDEFINES = [ 'PBRT_HAS_OPENEXR' ])
parallel_libs = [ 'pthread' ]
if arch != 'darwin':
    # 32-bit build
    env.Append(CPPDEFINES = [ 'PBRT_POINTER_SIZE=4' ])
    exr_libs = [ 'Iex', 'IlmImf', 'Imath', 'Iex', 'IlmThread', 'Half' ] 
else:
    env.Append(CPPDEFINES = [ 'PBRT_POINTER_SIZE=8' ])
    exr_libs = [ 'Iex', 'IlmImf', 'Imath', 'Iex', 'IlmThread', 'Half', 'z' ] 
    env.Append(CCFLAGS = [ '-m64' ],
               LINKFLAGS = [ '-m64', '-Z' ],
               CPPDEFINES = [ 'PBRT_HAS_64_BIT_ATOMICS' ])
    
tiff_libs = [ 'tiff' ]
Export('exr_libs', 'parallel_libs', 'tiff_libs')

setup_nice_print(env)

debug_env = env.Clone()
debug_env.Append(CPPDEFINES = [ 'DEBUG' ])
if arch == 'darwin':
    debug_env.Append(CPPDEFINES = [ 'PBRT_STATS_DTRACE' ])
else:
    debug_env.Append(CPPDEFINES = [ 'PBRT_STATS_COUNTERS' ])
build_envs['debug'] = debug_env

coverage_env = debug_env.Clone()
coverage_env.Append(CCFLAGS = [ '-fprofile-arcs', '-ftest-coverage' ])
coverage_env.Append(LINKFLAGS = [ '-fprofile-arcs', '-ftest-coverage' ])
build_envs['coverage'] = coverage_env

release_env = env.Clone()
release_env.Append(CCFLAGS = [ '-O2', '-finline-functions' ])
release_env.Append(CCFLAGS = [ '-msse3', '-mfpmath=sse', '-march=nocona' ])
release_env.Append(CPPDEFINES = [ 'NDEBUG' ])
build_envs['release'] = release_env

perf_env = release_env.Clone()
perf_env.Append(CCFLAGS = [ '-finstrument-functions' ])
perf_env.Append(LIBS = [ 'Saturn' ])
build_envs['perf'] = perf_env

stats_env = release_env.Clone()
if arch == 'darwin':
    stats_env.Append(CPPDEFINES = [ 'PBRT_STATS_DTRACE' ])
else:
    stats_env.Append(CPPDEFINES = [ 'PBRT_STATS_COUNTERS' ])
build_envs['stats'] = stats_env

release_env.Append(CPPDEFINES = [ 'PBRT_STATS_NONE' ])
perf_env.Append(CPPDEFINES = [ 'PBRT_STATS_NONE' ])

###########################################################################

    # exrcheck:
    #	@echo -n Checking for EXR installation... 
    #	@$(CXX) $(CXXFLAGS) -o exrcheck exrcheck.cpp $(LIBS) || \
        #		(cat exrinstall.txt; exit 1)

for target in build_envs:
    env = build_envs[target]
    Export('env')
    output = SConscript(dirs = '.',
                        variant_dir = 'build/' + arch + '-' + target)
    env.Alias(target, output['defaults'])
    env.Default(target)
