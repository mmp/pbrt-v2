core_dir = ..\core

bison_args = -d -v -t -o $(core_dir)\pbrtparse.cpp $(core_dir)\pbrtparse.yy
flex_args = -o$(core_dir)\pbrtlex.cpp $(core_dir)\pbrtlex.ll
	
.PHONY : $(core_dir)\pbrtparse.cpp $(core_dir)\pbrtlex.cpp

$(core_dir)\pbrtparse.cpp $(core_dir)\pbrtparse.hpp : $(core_dir)\pbrtparse.yy
	if exist $(BISON_PROGRAM) \
	( \
		$(BISON_PROGRAM) $(bison_args) \
	) \

$(core_dir)\pbrtparse.hh : $(core_dir)\pbrtparse.hpp
	if exist $(core_dir)\pbrtparse.hh del $(core_dir)\pbrtparse.hh
	ren $(core_dir)\pbrtparse.hpp pbrtparse.hh

$(core_dir)\pbrtlex.cpp $(core_dir)\pbrtparse.hh : $(core_dir)\pbrtlex.ll
	if exist $(FLEX_PROGRAM) \
	( \
		$(FLEX_PROGRAM) $(flex_args) \
	) \