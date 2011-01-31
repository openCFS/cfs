# ------------- MAKEFILE FOR FILEREADER FOR COUPLING -----------------------
# Author: szoerner.
# Date:   13.12.2010.
# --------------------------------------------------------------------------
#

SET(CPLREADER_CFS_SRCS
  FileReader_CfsHdf5.cc
)

ADD_LIBRARY(cplreader_cfs STATIC ${CPLREADER_CFS_SRCS} )

# add sources of hdf5 reader
TARGET_LINK_LIBRARIES(cplreader_cfs
  ${CFS_HDF5_LIB}
  ${HDF5_CPP_LIBRARY}
  ${HDF5_LIBRARY})

