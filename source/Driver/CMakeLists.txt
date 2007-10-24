# Make sure the compiler can find include files from our CFSMESHIO library.
#INCLUDE_DIRECTORIES(${CFSMESHIO_SOURCE_DIR}/include ${HDF5_INCLUDE_DIR})

# Make sure the linker can find the CFSMESHIO library once it is built.
#LINK_DIRECTORIES()

SET(DRIVER_SRCS solveStepAcoustic.cc stdSolveStep.cc baseSolveStep.cc
  solveStepElec.cc basedriver.cc transientdriver.cc staticdriver.cc
  harmonicDriver.cc multiSequenceDriver.cc assemble.cc singleDriver.cc
  iterSolveStep.cc solveStepPiezo.cc
  solveStepStokesFluid.cc eigenFrequencyDriver.cc formsContext.cc
  solveStepODE.cc solveStepMpCCI.cc solveStepMagHyst.cc)

ADD_LIBRARY(driver STATIC ${DRIVER_SRCS})

SET(TARGET_LL
  paramident
  utils-olas
  algsys-olas
  graph-olas
  utils
)

TARGET_LINK_LIBRARIES(driver ${TARGET_LL})
