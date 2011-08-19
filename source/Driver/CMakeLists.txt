SET(DRIVER_SRCS
  assemble.cc
  basedriver.cc
  baseSolveStep.cc
  eigenFrequencyDriver.cc
  formsContext.cc
  harmonicDriver.cc
  iterSolveStep.cc
  multiSequenceDriver.cc
  singleDriver.cc
  solveStepAcoustic.cc
  solveStepElec.cc
  solveStepMagEdge.cc
#  solveStepFluidMech.cc
#  solveStepMagHyst.cc
#  solveStepMech.cc
#  solveStepMpCCI.cc
#  solveStepODE.cc
#  solveStepPiezo.cc
#  solveStepMicroPiezo.cc
#  solveStepSmooth.cc
  staticdriver.cc
  stdSolveStep.cc
  transientdriver.cc
  )

ADD_LIBRARY(driver STATIC ${DRIVER_SRCS})

SET(TARGET_LL
#  paramident
  utils-olas
  algsys-olas
  graph-olas
  utils
)

TARGET_LINK_LIBRARIES(driver ${TARGET_LL})
