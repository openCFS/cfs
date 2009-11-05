# Just redefine __builtin_isnan to ::isnan for Intel compiler for 
# assemble.cc
IF(CFS_CXX_COMPILER_NAME STREQUAL "ICC")
    GET_SOURCE_FILE_PROPERTY(OPTCONDFLAGS OptimalityCondition.cc
        COMPILE_FLAGS)
    #    MESSAGE("OPTCONDFLAGS ${OPTCONDFLAGS}")
    IF(OPTCONDFLAGS STREQUAL NOTFOUND)
      SET(OPTCONDFLAGS "")
    ENDIF(OPTCONDFLAGS STREQUAL NOTFOUND)
    
    SET_SOURCE_FILES_PROPERTIES(assemble.cc
      PROPERTIES COMPILE_FLAGS "${OPTCONDFLAGS} -D__builtin_isnan=::isnan -D__builtin_isinf=::isinf")
ENDIF(CFS_CXX_COMPILER_NAME STREQUAL "ICC")

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
  solveStepFluidMech.cc
  solveStepMagHyst.cc
  solveStepMech.cc
  solveStepMpCCI.cc
  solveStepODE.cc
  solveStepMicroPiezo.cc
  solveStepPiezo.cc
  solveStepSmooth.cc
  staticdriver.cc
  stdSolveStep.cc
  transientdriver.cc
  )

ADD_LIBRARY(driver STATIC ${DRIVER_SRCS})

SET(TARGET_LL
  paramident
  utils-olas
  algsys-olas
  graph-olas
  utils
)

TARGET_LINK_LIBRARIES(driver ${TARGET_LL})
