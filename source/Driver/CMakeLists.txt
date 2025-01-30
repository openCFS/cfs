SET(DRIVER_SRCS
  Assemble.cc
  BaseDriver.cc
  EigenFrequencyDriver.cc
  EigenValueDriver.cc
  BucklingDriver.cc
  FormsContexts.cc
  HarmonicDriver.cc
  Harmonic25DDriver.cc
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
  SolveSteps/SolveStepEB.cc
  SolveSteps/StdSolveStep.cc
  SolveSteps/PrescribedSolveStep.cc
  SolveSteps/SolveStepHyst.cc
  )

ADD_LIBRARY(driver STATIC ${DRIVER_SRCS})

SET(TARGET_LL
  utils-olas
  algsys-olas
  graph-olas
  utils
  preciceadapter
  matvec
)

TARGET_LINK_LIBRARIES(driver ${TARGET_LL})
