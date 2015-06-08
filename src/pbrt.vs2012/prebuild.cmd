@ECHO OFF
setlocal EnableDelayedExpansion
set "winbison="
set "winflex="

set "drives="
for /f "delims=" %%a in ('fsutil fsinfo drives') do @set "drives=%%a"

REM :~8 is to slice off "Drives: " returned by fsutil
for %%i in (%drives:~8%) do (
	if exist %%iNUL (
		pushd %%i
		if not defined winbison (
			for /f "delims=" %%a in (
				'dir win_bison.exe /s /b 2^>nul') do @set "winbison=%%a"
		)

		if not defined winflex (
			for /f "delims=" %%a in (
				'dir win_flex.exe /s /b 2^>nul') do @set "winflex=%%a"
		)
		popd
	)
)

nmake BISON_PROGRAM="%winbison%" FLEX_PROGRAM="%winflex%" prebuild.makefile
cd ..\3rdparty\ilmbase-1.0.2
buildLUTS.cmd

@ECHO ON