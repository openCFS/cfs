INCLUDE_DIRECTORIES(/home/data/programs/ansys/v100/ansys/customize/include)

SET(CMAKE_Fortran_FLAGS "-fcray-pointer -DLINUXOP64_SYS")

IF(DEBUG)
  SET(CMAKE_Fortran_FLAGS "-g ${CMAKE_Fortran_FLAGS}")
IF(DEBUG)

SET(SIMOUTANSYSRST_SRCS simOutputRST.cc ResWr.F)

ADD_LIBRARY(simoutansysrst STATIC ${SIMOUTANSYSRST_SRCS})

TARGET_LINK_LIBRARIES(simoutansysrst bin pgf90 ${GFORTRAN_LIBRARY})

SET_TARGET_PROPERTIES(simoutansysrst PROPERTIES LINKER_LANGUAGE CXX)

# Add the ANSYS ResWrDemo program for test purposes
ADD_EXECUTABLE(ResWrDemo ResWrDemo.F)

TARGET_LINK_LIBRARIES(ResWrDemo simoutansysrst)
