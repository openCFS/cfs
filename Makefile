# ------------- MAKEFILE FOR ADAPT_GRID ------------------------------------
# Author: LSE.
# Date:   13.09.2002.
# --------------------------------------------------------------------------
#
# ------------- Name for library
#
libname = adaptgrid
#
# ------------ Path to a place of location of the library
#
library = ../lib/$(MACHINE)/lib$(libname).so
#
# ------------ Source
#
src = Chronometer.cpp Edge.cpp GeometrySensor.cpp GridLevel.cpp Hexahedron.cpp Jacobian.cpp MeshReader.cpp MultilevelGrid.cpp Octahedron.cpp Quadrangle.cpp Tetrahedron.cpp  Vertex.cpp 
#
# ------------------------- Flags
#
CC=g++
FLAG_MKLIB= -v -shared
CFLAGS1=-O3 -m486 -c
#
# ------------------------ Rules
#
.SUFFIXES: .cpp .o
.cpp.o:
	$(CC) $(CFLAGS1) $<
#
# ------------------------ Goals
#                               ---- Make library
$(library): $(src:.cpp=.o)
	$(CC) $(FLAG_MKLIB) -o $@ $(?:.cpp=.o)
#
#                               ---- Clean
#
clean:
	-rm *.o    
#
#
