SET(DRIVER_SRCS
  assemble.cc
  stdSolveStep.cc
  transientdriver.cc
  harmonicDriver.cc
  solveStepMicroPiezo.cc
  multiSequenceDriver.cc
  solveStepPiezo.cc
  basedriver.cc
  solveStepAcoustic.cc
  solveStepElec.cc
  formsContext.cc
  solveStepMagHyst.cc
  solveStepFluidMech.cc
  solveStepSmooth.cc
  eigenFrequencyDriver.cc
  iterSolveStep.cc
  solveStepODE.cc
  staticdriver.cc
  solveStepMech.cc
  singleDriver.cc
  baseSolveStep.cc
  )
# solveStepMpCCI.cc

ADD_LIBRARY(driver STATIC ${DRIVER_SRCS})

TARGET_LINK_LIBRARIES(driver
  algsys-olas
  cfsgeneral
  coupledpde
  datainout
  design
  domain
  elements
  forms
  logging
  materials
  matvec
  optimization
  paramh
  paramident
  pde
  scripting
  utils  
)
