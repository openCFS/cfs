# Make sure the compiler can find include files from our CFSMESHIO library.
#INCLUDE_DIRECTORIES(${CFSMESHIO_SOURCE_DIR}/include ${HDF5_INCLUDE_DIR})

# Make sure the linker can find the CFSMESHIO library once it is built.
#LINK_DIRECTORIES()

SET(COUPLED_PDE_SRCS itercoupledpde.cc coupledpdedef.cc pdecoupling.cc
  couplingmemento.cc DirectCoupledPDE.cc BasePairCoupling.cc
  PiezoCoupling.cc AcouMechCoupling.cc)

ADD_LIBRARY(coupledpde STATIC ${COUPLED_PDE_SRCS})
