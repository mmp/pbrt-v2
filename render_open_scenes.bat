SET sceneName=anim-killeroos-moving

REM Run pbrt
bin\pbrt_x64_Release.exe --ncores 4 --outfile results\pbrt.exr scenes\%sceneName%.pbrt

REM Conver exr image to tiff image
bin\exrtotiff.exe results\pbrt.exr results\%sceneName%.tiff

REM Display tiff image
results\%sceneName%.tiff



