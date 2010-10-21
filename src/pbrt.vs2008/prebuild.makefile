core_dir = ..\core

bison_bin = bison.exe
bison_cygwin_bin = c:\cygwin\bin\$(bison_bin)
bison_args = -d -v -t -o $(core_dir)\pbrtparse.cpp $(core_dir)\pbrtparse.yy

flex_bin = flex.exe
flex_cygwin_bin = c:\cygwin\bin\$(flex_bin)
flex_args = -o$(core_dir)\pbrtlex.cpp $(core_dir)\pbrtlex.ll


.PHONY : $(core_dir)\pbrtparse.cpp $(core_dir)\pbrtlex.cpp


$(core_dir)\pbrtparse.cpp : $(core_dir)\pbrtparse.yy
	if exist $(bison_cygwin_bin) \
	( \
		$(bison_cygwin_bin) $(bison_args) \
	) \
	else \
	( \
		$(bison_bin) $(bison_args) \
	)

$(core_dir)\pbrtlex.cpp : $(core_dir)\pbrtlex.ll
	if exist $(flex_cygwin_bin) \
	( \
		$(flex_cygwin_bin) $(flex_args) \
	) \
	else \
	( \
		$(flex_bin) $(flex_args) \
	)