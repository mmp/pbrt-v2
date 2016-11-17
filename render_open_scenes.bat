SET sceneName=bunny

REM Run pbrt
bin\pbrt.exe scenes\%sceneName%.pbrt

REM Conver exr image to tiff image
bin\exrtotiff.exe pbrt.exr %sceneName%.tiff

REM Display tiff image
%sceneName%.tiff



