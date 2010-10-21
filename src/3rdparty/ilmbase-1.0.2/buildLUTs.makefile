tmp_dir = ..\..\..\tmp


.PHONY : toFloat.h eLut.h


toFloat.h : $(tmp_dir)\toFloat.exe
	$(tmp_dir)\toFloat.exe > toFloat.h

eLut.h : $(tmp_dir)\eLut.exe
	$(tmp_dir)\eLut.exe > eLut.h

$(tmp_dir)\toFloat.exe : toFloat.cpp
	cl /nologo /GR /EHsc toFloat.cpp /Fo$(tmp_dir)\toFloat.obj /Fe$(tmp_dir)\toFloat.exe

$(tmp_dir)\eLut.exe : eLut.cpp
	cl /nologo /GR /EHsc eLut.cpp /Fo$(tmp_dir)\eLut.obj /Fe$(tmp_dir)\eLut.exe