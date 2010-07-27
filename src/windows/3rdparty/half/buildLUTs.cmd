@cl /nologo /GR /EHsc toFloat.cpp /FetoFloat.exe
@cl /nologo /GR /EHsc eLut.cpp /FeeLut.exe

@.\toFloat.exe > toFloat.h
@.\eLut.exe > eLut.h

@del .\toFloat.obj .\toFloat.exe
@del .\eLut.exe .\eLut.obj
