# -*- mode: python -*-

import sys, platform
arch = sys.platform

######################################################################
# user-configurable section

has_openexr = True
exr_includes = [ '/usr/local/include', '/opt/local/include', '/opt/local/include/OpenEXR',
                 '/usr/local/include/OpenEXR', '/usr/include/OpenEXR' ]
exr_libdir = [ '/opt/local/lib' ]

has_dtrace = (arch == 'darwin')
Export('has_dtrace')

#has_gcd = float(platform.mac_ver()[0][:4] >= 10.6)
has_gcd = False

build_64bit = True

parallel_libs = [ 'pthread' ]
Export('parallel_libs')
    
tiff_libs = [ ]
# tiff_libs = [ 'tiff' ]
tiff_includes = [ ]
tiff_libdir = [ ]
Export('tiff_libs')

######################################################################
## Configure generic environment

Decider('MD5-timestamp')
#CacheDir('scons-cache')

#import os
#print "Pruning scons cache..."
#os.system('cd scons-cache && find . -type f -atime +1 -delete')

def setup_nice_print(env):
    if ARGUMENTS.get('VERBOSE') != '1':
        env['YACCCOMSTR'] = "Compiling $TARGET"
        env['LEXCOMSTR'] = "Compiling $TARGET"
        env['CCCOMSTR'] = "Compiling $TARGET"
        env['CXXCOMSTR'] = "Compiling $TARGET"
        env['LINKCOMSTR'] = "Linking $TARGET"
        env['ARCOMSTR'] = "Linking $TARGET"

env = Environment(CCFLAGS = [ '-Wall', '-g' ],
                  CPPPATH = [ '#core', '#', '.' ] + tiff_includes,
                  LIBPATH = tiff_libdir,
                  YACCFLAGS = [ '-d', '-v', '-t' ],
                  YACCHXXFILESUFFIX = '.hh',
                  ENV = { 'PATH' : [ '/usr/local/bin', '/usr/bin', '/bin', 
                                     '/usr/sbin', '/sbin' ] })
if build_64bit:
    env.Append(CCFLAGS = [ '-m64' ],
               LINKFLAGS = [ '-m64' ])

setup_nice_print(env)

if has_openexr:
    env.Append(CPPPATH = exr_includes)
    env.Append(LIBPATH = exr_libdir)
    env.Append(CPPDEFINES = [ 'PBRT_HAS_OPENEXR' ])
    if arch != 'darwin':
        exr_libs = [ 'Iex', 'IlmImf', 'Imath', 'Iex', 'IlmThread', 'Half' ] 
    else:
        exr_libs = [ 'Iex', 'IlmImf', 'Imath', 'Iex', 'IlmThread', 'Half', 'z' ] 
else:
    exr_libs = [ ]
Export('exr_libs')

if has_dtrace:
    env.Append(BUILDERS = { 'DTrace' : 
                            Builder(action = '/usr/sbin/dtrace -h -s $SOURCE -o $TARGET')})

if has_gcd:
    env.Append(CPPDEFINES = [ 'PBRT_USE_GRAND_CENTRAL_DISPATCH' ])

######################################################################

build_envs = { }

debug_env = env.Clone()
debug_env.Append(CPPDEFINES = [ 'DEBUG' ])
if has_dtrace:
    debug_env.Append(CPPDEFINES = [ 'PBRT_PROBES_DTRACE' ])
else:
    debug_env.Append(CPPDEFINES = [ 'PBRT_PROBES_NONE' ])
build_envs['debug'] = debug_env

release_env = env.Clone()
release_env.Append(CCFLAGS = [ '-O2', '-finline-functions' ])
release_env.Append(CCFLAGS = [ '-msse3', '-mfpmath=sse', '-march=nocona' ])
release_env.Append(CPPDEFINES = [ 'NDEBUG' ])
build_envs['release'] = release_env

stats_env = release_env.Clone()
if has_dtrace:
    stats_env.Append(CPPDEFINES = [ 'PBRT_PROBES_DTRACE' ])
else:
    stats_env.Append(CPPDEFINES = [ 'PBRT_PROBES_COUNTERS' ])
build_envs['stats'] = stats_env

release_env.Append(CPPDEFINES = [ 'PBRT_PROBES_NONE' ])

for target in build_envs:
    env = build_envs[target]
    Export('env')
    output = SConscript(dirs = '.',
                        variant_dir = 'build/' + arch + '-' + target)
    env.Alias(target, output['defaults'])
    env.Default(target)
