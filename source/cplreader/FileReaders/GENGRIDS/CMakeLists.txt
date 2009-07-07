# ------------- MAKEFILE FOR FILEREADER FOR COUPLING -----------------------
# Author: LSE.
# Date:   05.09.2003.
# --------------------------------------------------------------------------
#

SET(CPLREADER_GENGRIDS_SRCS
  FileReader_GENGRIDS.cc
  SphericalShellGenerator.cc
  RectilinearUniformGridGenerator.cc
  ReferenceElementGenerator.cc
)

ADD_LIBRARY(cplreader_gengrids STATIC ${CPLREADER_GENGRIDS_SRCS})
