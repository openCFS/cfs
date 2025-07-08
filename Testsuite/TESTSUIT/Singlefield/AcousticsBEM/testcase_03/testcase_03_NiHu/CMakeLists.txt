cmake_minimum_required (VERSION 3.15.0 FATAL_ERROR)

# set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# set the binary directory to /NiHu
message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")

project(Hello_NiHu)

## get folder's name and set names of .cc and executable accordingly
#
#  documentation: The .cc file to create NiHu's executabe has to have 
#                 the same name as the NiHu folder!
#

# Get the name of the current directory (excluding the path)
get_filename_component(folder_name ${CMAKE_CURRENT_SOURCE_DIR} NAME)

# Set the source file
set(source_file "${folder_name}.cc")

# Set the executable name
set(executable_name "${folder_name}_output")

# Create the executable
add_executable(${executable_name} ${source_file})

# add_executable(hello_world hello_world.cc)
# enable concepts for the target
target_compile_features(${executable_name} PRIVATE cxx_std_14)
target_compile_options(${executable_name} PRIVATE -fconcepts)

set_target_properties(${executable_name} PROPERTIES LINKER_LANGUAGE CXX)

## mod 24 - (1) include directories
target_include_directories(${executable_name} PUBLIC /home/candussi/Devel_RK1/FEM_BEM_20240108/BIN/build_23/cfsdeps/nihu/src/NiHu-install)
target_include_directories(${executable_name} PUBLIC /home/candussi/Devel_RK1/FEM_BEM_20240108/BIN/build_23/cfsdeps/nihu/src/NiHu-install/include/)
target_include_directories(${executable_name} PUBLIC /home/candussi/Devel_RK1/FEM_BEM_20240108/BIN/build_23/cfsdeps/nihu/src/NiHu-install/include/fmm)
target_include_directories(${executable_name} PUBLIC /home/candussi/Devel_RK1/FEM_BEM_20240108/BIN/build_23/cfsdeps/nihu/src/NiHu-install/include/library)
target_include_directories(${executable_name} PUBLIC /home/candussi/Devel_RK1/FEM_BEM_20240108/BIN/build_23/cfsdeps/nihu/src/NiHu-install/include/Eigen/eigen-3.3.9)

## additional includes - NiHu base headers
target_include_directories(${executable_name} PUBLIC /home/candussi/Devel_RK1/FEM_BEM_20231204/build_debug/cfsdeps/nihu/src/NiHu/src/core)
target_include_directories(${executable_name} PUBLIC /home/candussi/Devel_RK1/FEM_BEM_20231204/build_debug/cfsdeps/nihu/src/NiHu/src/tmp)

## mod 24 - (2) set directories' aliases
set(NIHU_EXTERNAL_LIBRARY_PATH /home/candussi/Devel_RK1/FEM_BEM_20240108/BIN/build_23/cfsdeps/nihu/src/NiHu-build/lib/nihu.a)
set(NIHU_FMM_EXTERNAL_LIBRARY_PATH /home/candussi/Devel_RK1/FEM_BEM_20240108/BIN/build_23/cfsdeps/nihu/src/NiHu-build/lib/nihu_fmm.a)
# set(NIHU_EIGEN_EXTERNAL_LIBRARY_PATH /home/candussi/Devel_RK1/FEM_BEM_20240108/BIN/build_23/cfsdeps/nihu/src/NiHu-install/include/Eigen)
# set(NIHU_EIGEN_EXTERNAL_LIBRARY_PATH /home/candussi/Devel_RK1/FEM_BEM_20240108/BIN/build_23/cfsdeps/nihu/src/NiHu-install/include/eigen-3.3.9/Eigen/Eigen)

# create a new library target for the external library
add_library(nihu STATIC IMPORTED)
add_library(nihu_fmm STATIC IMPORTED)
# add_library(weighted_residual STATIC IMPORTED)

## mod 24 - (3) create library targets for external libraries
# add_library(Eigen STATIC IMPORTED)

## mod 24 - (4) set location of external libraries for targets
# set_target_properties(Eigen PROPERTIES IMPORTED_LOCATION ${NIHU_EIGEN_EXTERNAL_LIBRARY_PATH})

## mod 24 - (5) link libraries' targets with external
# target_link_libraries(hello_world Eigen)

# set the location of the external library file for the target
set_target_properties(nihu PROPERTIES IMPORTED_LOCATION ${NIHU_EXTERNAL_LIBRARY_PATH})
set_target_properties(nihu_fmm PROPERTIES IMPORTED_LOCATION ${NIHU_FMM_EXTERNAL_LIBRARY_PATH})

# link executable or library target with the external library
target_link_libraries(${executable_name} nihu)
target_link_libraries(${executable_name} nihu_fmm)