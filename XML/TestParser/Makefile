#
# Makefile for compiling the CFS - XML - Test Parser
#

TARGET = ../parse
SOURCE = TestParser
PLATFORM = LINUX
XERCESBASE = /home/data/proj/SoftwareLibraries/xerces-c
XMLSCHEMA = $(dir $(PWD))
INCPATH = $(XERCESBASE)/include
LIBPATH = $(XERCESBASE)/lib/LINUX
CXX = g++
STDOPTS = -ansi -Wall -pedantic -ggdb -I$(INCPATH) -D$(PLATFORM)
DEBUGOPTS =
CXXOPTS = $(STDOPTS)
LD = g++
LDOPTS = -L$(LIBPATH) -lxerces-c

$(TARGET):	$(SOURCE).o SPHandler.o;
		$(LD) $(LDOPTS) -o $(TARGET) $(SOURCE).o SPHandler.o

$(SOURCE).o:	$(SOURCE).cc;	$(CXX) -c $(CXXOPTS) $(SOURCE).cc

SPHandler.o:	SPHandler.cc;
		$(CXX) -c $(CXXOPTS) -DXMLSCHEMA='"$(XMLSCHEMA)"' SPHandler.cc

clean:	;	/bin/rm $(SOURCE).o SPHandler.o $(TARGET)
