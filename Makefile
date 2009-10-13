ARCH = $(shell uname)

# user-configuration section

EXRINCLUDE=-I/usr/local/include/OpenEXR -I/usr/include/OpenEXR
EXRLIBDIR=-L/usr/local/lib

DEFS=-DPBRT_STATS_NONE -DPBRT_POINTER_SIZE=4 -DPBRT_HAS_PTHREADS -DPBRT_HAS_OPENEXR

#########################################################################

LEX=flex
YACC=bison -d -v -t
LEXLIB = -lfl

EXRLIBS=$(EXRLIBDIR) -Bstatic -lIex -lIlmImf -lIlmThread -lImath -lIex -lHalf -Bdynamic
ifeq ($(ARCH),Linux)
  EXRLIBS += -lpthread
endif
ifeq ($(ARCH),Darwin)
  EXRLIBS += -lz
endif

CC=gcc
CXX=g++
LD=$(CXX) $(OPT)
OPT=-O2 -m32 -msse2 -mfpmath=sse
INCLUDE=-I. -Icore $(EXRINCLUDE)
WARN=-Wall
CWD=$(shell pwd)
CXXFLAGS=$(OPT) $(INCLUDE) $(WARN) $(DEFS)
CCFLAGS=$(CXXFLAGS)
LIBS=$(LEXLIB) $(EXRLIBDIR) $(EXRLIBS) -lm 

LIBSRCS=$(wildcard core/*.cpp) core/pbrtlex.cpp core/pbrtparse.cpp
LIBSRCS += $(wildcard accelerators/*.cpp cameras/*.cpp film/*.cpp filters/*.cpp )
LIBSRCS += $(wildcard integrators/*.cpp lights/*.cpp materials/*.cpp renderers/*.cpp )
LIBSRCS += $(wildcard samplers/*.cpp shapes/*.cpp textures/*.cpp volumes/*.cpp)

#LIBOBJS=$(addprefix objs/, $(notdir $(LIBSRCS:.cpp=.o)))
LIBOBJS=$(addprefix objs/, $(subst /,_,$(LIBSRCS:.cpp=.o)))

HEADERS = $(wildcard */*.h)

default: bin/pbrt #tools

pbrt: bin/pbrt

dirs:
	/bin/mkdir -p bin objs

$(LIBOBJS): $(HEADERS)

.PHONY: dirs tools exrcheck

objs/libpbrt.a: $(LIBOBJS)
	@echo "Building the core rendering library (libpbrt.a)"
	@ar rcs $@ $(LIBOBJS)

objs/accelerators_%.o: accelerators/%.cpp
	@echo "Building object $@"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

objs/cameras_%.o: cameras/%.cpp
	@echo "Building object $@"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

objs/core_%.o: core/%.cpp
	@echo "Building object $@"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

objs/film_%.o: film/%.cpp
	@echo "Building object $@"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

objs/filters_%.o: filters/%.cpp
	@echo "Building object $@"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

objs/integrators_%.o: integrators/%.cpp
	@echo "Building object $@"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

objs/lights_%.o: lights/%.cpp
	@echo "Building object $@"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

objs/main_%.o: main/%.cpp
	@echo "Building object $@"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

objs/materials_%.o: materials/%.cpp
	@echo "Building object $@"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

objs/renderers_%.o: renderers/%.cpp
	@echo "Building object $@"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

objs/samplers_%.o: samplers/%.cpp
	@echo "Building object $@"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

objs/shapes_%.o: shapes/%.cpp
	@echo "Building object $@"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

objs/textures_%.o: textures/%.cpp
	@echo "Building object $@"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

objs/volumes_%.o: volumes/%.cpp
	@echo "Building object $@"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

objs/pbrt.o: main/pbrt.cpp
	@echo "Building object $@"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

bin/pbrt: dirs objs/libpbrt.a objs/pbrt.o
	@echo "Linking $@"
	@$(CXX) $(CXXFLAGS) -o $@ objs/pbrt.o objs/libpbrt.a $(LIBS)

core/pbrtlex.cpp: core/pbrtlex.ll core/pbrtparse.cpp
	@echo "Lex'ing pbrtlex.ll"
	@$(LEX) -o$@ core/pbrtlex.ll

core/pbrtparse.cpp: core/pbrtparse.yy
	@echo "YACC'ing pbrtparse.yy"
	@$(YACC) -o $@ core/pbrtparse.yy
	@if [ -e core/pbrtparse.cpp.h ]; then /bin/mv core/pbrtparse.cpp.h core/pbrtparse.hh; fi
	@if [ -e core/pbrtparse.hpp ]; then /bin/mv core/pbrtparse.hpp core/pbrtparse.hh; fi

$(RENDERER_BINARY): $(RENDERER_OBJS) $(CORE_LIB)

clean:
	rm -f objs/* bin/* core/pbrtlex.[ch]* core/pbrtparse.[ch]*
	#(cd tools && $(MAKE) clean)

objs/exrio.o: exrcheck

exrcheck:
	@echo -n Checking for EXR installation... 
	@$(CXX) $(CXXFLAGS) -o exrcheck exrcheck.cpp $(LIBS) || \
		(cat exrinstall.txt; exit 1)
