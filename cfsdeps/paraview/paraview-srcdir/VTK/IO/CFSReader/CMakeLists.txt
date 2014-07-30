set(Module_SRCS
  vtkCFSReader.cxx
  hdf5Reader.cc
  hdf5Common.cc
)

SET_SOURCE_FILES_PROPERTIES(
  hdf5Reader
  hdf5Common
  WRAP_EXCLUDE
)

vtk_module_library(vtkIOCFSReader ${Module_SRCS})

