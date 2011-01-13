# -*- mode: python -*-

import sys
arch = sys.platform

Import('env', 'has_dtrace', 'parallel_libs', 'exr_libs', 'tiff_libs')

if has_dtrace == True:
    env.DTrace('core/dtrace.h', 'core/dtrace.d')
    env.Depends('core/dtrace.h', 'core/dtrace.d')

renderer_src = [ 'main/pbrt.cpp' ]

core_src = [ 'core/api.cpp',           'core/camera.cpp',         'core/diffgeom.cpp',
             'core/error.cpp',         'core/film.cpp',           'core/fileutil.cpp',
             'core/filter.cpp',
             'core/floatfile.cpp',     'core/geometry.cpp',       'core/imageio.cpp', 
             'core/integrator.cpp',    'core/intersection.cpp',   'core/light.cpp', 
             'core/material.cpp',      'core/memory.cpp',         'core/montecarlo.cpp',
             'core/paramset.cpp',      'core/parser.cpp',         'core/primitive.cpp',
             'core/parallel.cpp',      'core/probes.cpp',         'core/progressreporter.cpp', 
             'core/quaternion.cpp',    'core/reflection.cpp',     'core/renderer.cpp',
             'core/rng.cpp',           'core/sampler.cpp',        'core/scene.cpp',
             'core/sh.cpp',            'core/shrots.cpp',         'core/shape.cpp',
             'core/spectrum.cpp',      'core/texture.cpp',        'core/timer.cpp', 
             'core/transform.cpp',     'core/volume.cpp' ]

parser_flex = env.CXXFile('core/pbrtlex.ll')
parser_bison = env.CXXFile('core/pbrtparse.yy')
env.Depends('core/pbrtlex.ll', parser_bison)
parser_bison.pop()  # yuck a muck, but stop trying to link in pbrtparse.hh!
core_src = [ parser_flex, parser_bison ] + core_src


accelerators_src = [ 'accelerators/bvh.cpp', 
                     'accelerators/grid.cpp',
                     'accelerators/kdtreeaccel.cpp' ]
cameras_src = [ 'cameras/environment.cpp', 
                'cameras/orthographic.cpp', 
                'cameras/perspective.cpp' ]
film_src = [ 'film/image.cpp' ]
filters_src = [ 'filters/box.cpp',              'filters/gaussian.cpp', 
                'filters/mitchell.cpp',         'filters/sinc.cpp',
                'filters/triangle.cpp' ]
integrators_src = [ 'integrators/ambientocclusion.cpp',      'integrators/diffuseprt.cpp',
                    'integrators/dipolesubsurface.cpp',      'integrators/directlighting.cpp', 
                    'integrators/emission.cpp',              'integrators/glossyprt.cpp',
                    'integrators/igi.cpp',                   'integrators/irradiancecache.cpp', 
                    'integrators/path.cpp',                  'integrators/photonmap.cpp', 
                    'integrators/single.cpp',                'integrators/useprobes.cpp',
                    'integrators/whitted.cpp' ]
lights_src = [ 'lights/diffuse.cpp',           'lights/distant.cpp',
               'lights/goniometric.cpp',       'lights/infinite.cpp',
               'lights/point.cpp',             'lights/projection.cpp', 
               'lights/spot.cpp' ]
materials_src = [ 'materials/glass.cpp',         'materials/kdsubsurface.cpp',
                  'materials/matte.cpp',         'materials/measured.cpp',
                  'materials/metal.cpp',         'materials/mirror.cpp',
                  'materials/mixmat.cpp',        'materials/plastic.cpp', 
                  'materials/substrate.cpp',     'materials/subsurface.cpp',
                  'materials/translucent.cpp',   'materials/uber.cpp',
                  'materials/shinymetal.cpp',
                  ]
renderers_src = [ 'renderers/aggregatetest.cpp',   'renderers/createprobes.cpp',
                  'renderers/metropolis.cpp',      'renderers/samplerrenderer.cpp',
                  'renderers/surfacepoints.cpp' ]
samplers_src = [ 'samplers/adaptive.cpp',         'samplers/bestcandidate.cpp',
                 'samplers/halton.cpp',           'samplers/lowdiscrepancy.cpp', 
                 'samplers/random.cpp',           'samplers/stratified.cpp' ]
shapes_src = [ 'shapes/cone.cpp',        'shapes/cylinder.cpp',
               'shapes/disk.cpp',        'shapes/heightfield.cpp',
               'shapes/hyperboloid.cpp', 'shapes/loopsubdiv.cpp',
               'shapes/nurbs.cpp',       'shapes/paraboloid.cpp',
               'shapes/sphere.cpp',      'shapes/trianglemesh.cpp' ]
textures_src = [ 'textures/bilerp.cpp',          'textures/checkerboard.cpp',
                 'textures/constant.cpp',        'textures/dots.cpp',
                 'textures/fbm.cpp',             'textures/imagemap.cpp', 
                 'textures/marble.cpp',          'textures/mix.cpp',
                 'textures/scale.cpp',           'textures/uv.cpp',
                 'textures/windy.cpp',           'textures/wrinkled.cpp' ]
volumes_src = [ 'volumes/exponential.cpp',       'volumes/homogeneous.cpp', 
                'volumes/volumegrid.cpp' ]


lib_src = [ core_src + accelerators_src + cameras_src + film_src + filters_src +
            integrators_src + lights_src + materials_src + renderers_src +
            samplers_src + shapes_src + textures_src + volumes_src ]

output = { }
output['pbrt_lib'] = env.StaticLibrary('libpbrt', lib_src)

env_libs = [ ]
output['pbrt'] = env.Program('pbrt', renderer_src + output['pbrt_lib'],
                             LIBS = env_libs + exr_libs + parallel_libs)
if arch == 'win32':
    env.AddPostAction(output['pbrt'], 
                      'mt.exe /outputresource:"$TARGET;#1" /manifest "${TARGET}.manifest" /nologo')
output['defaults'] = [ output['pbrt' ] ]


if len(exr_libs) > 0:
    output['exrdiff'] = env.Program('exrdiff', [ 'tools/exrdiff.cpp' ], 
                                    LIBS = env_libs + exr_libs)
    output['exravg'] = env.Program('exravg', [ 'tools/exravg.cpp' ], 
                                   LIBS = env_libs + exr_libs)
    output['bsdftest'] = env.Program('bsdftest', [ 'tools/bsdftest.cpp' ] + 
                                     output['pbrt_lib'], 
                                     LIBS = env_libs + exr_libs + parallel_libs)
    output['defaults'] = output['defaults'] +  \
        [ output['exrdiff'], output['exravg'], output['bsdftest'] ]

if len(exr_libs) > 0 and len(tiff_libs) > 0:
    output['exrtotiff'] = env.Program('exrtotiff', [ 'tools/exrtotiff.cpp' ], 
                                      LIBS = env_libs + exr_libs + tiff_libs)
    output['defaults'] = output['defaults'] + output['exrtotiff']

Return('output')
