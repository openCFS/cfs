# this replaces the original CMakeLists.txt as we need a (static) lib only
cmake_minimum_required(VERSION 3.3)
project(dumas)

if(USE_OPENMP)
  # only MMASolver.cpp has some loops, which might be questionable
  add_compile_definitions(MMA_WITH_OPENMP)
endif()

add_library(dumas STATIC src/mma/MMASolver.cpp src/gcmma/GCMMASolver.cpp)

set(MYLIB_HEADERS src/mma/MMASolver.h src/gcmma/GCMMASolver.h)
set_target_properties(dumas PROPERTIES PUBLIC_HEADER "${MYLIB_HEADERS}")

install(TARGETS dumas PUBLIC_HEADER DESTINATION include ARCHIVE DESTINATION lib)
