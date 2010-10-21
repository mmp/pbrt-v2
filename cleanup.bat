@rmdir /s /q bin
@rmdir /s /q lib
@rmdir /s /q tmp
@del /f /q /s *.user
@del /f /q /s *.suo
@del /f /q /s *.sdf
@del /f /q /s *.ncb
@cd src/pbrt.vs2010
@rmdir /s /q ipch