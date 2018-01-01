SET(DRIVER_SRCS
  Assemble.cc
  BaseDriver.cc
  EigenFrequencyDriver.cc
  FormsContexts.cc
  HarmonicDriver.cc
  MultiHarmonicDriver.cc
  MultiSequenceDriver.cc
  SingleDriver.cc
  StaticDriver.cc
  TransientDriver.cc
  InverseSourceDriver.cc
  AnalysisID.cc
  SolveSteps/BaseSolveStep.cc
  SolveSteps/IterSolveStep.cc
  SolveSteps/SolveStepElec.cc
  SolveSteps/StdSolveStep.cc
  )

ADD_LIBRARY(driver STATIC ${DRIVER_SRCS})

SET(TARGET_LL
  utils-olas
  algsys-olas
  graph-olas
  utils
  matvec
)

TARGET_LINK_LIBRARIES(driver ${TARGET_LL})
