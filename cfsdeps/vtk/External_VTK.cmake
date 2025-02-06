# VTK is an open-source software system for image processing, 3D graphics, volume rendering and visualization
# https://gitlab.kitware.com/vtk/vtk
# needed for Ensight reader
clear_depencency_variables()

# set mandatory variables for the macros in DependencyTools.cmake.
set(PACKAGE_NAME "vtk")
set(PACKAGE_VER "7.1") # must be 2 digits for include-dir to be correct
set(VTK_VERSION ${PACKAGE_VER}) # required later
set(PACKAGE_FILE "VTK-7.1.1.tar.gz") # note, one digit more and VERSION
set(PACKAGE_MD5 "daee43460f4e95547f0635240ffbc9cb")
set(DEPS_VER "-a") # set to "-a", "-b", when dependency changed with same PACKAGE_VER. Reset to "" with new PACKAGE_VER.

set(PACKAGE_MIRRORS "http://www.vtk.org/files/release/${PACKAGE_VER}/${PACKAGE_FILE}")
# add default mirrors to PACKAGE_MIRRORS or replace all with LOCAL_PACKAGE_FILE if we already have it
add_standard_mirrors_or_set_local()

# pure C
use_c_and_fortran(ON OFF)

# sets PRECOMPILED_PCKG_FILE to the full precompiled name including path
set_precompiled_pckg_file()

set(VTK_NAMES_LIBS "IOXML;IOXMLParser;expat;IOEnSight;FiltersSMP;FiltersParallel")
list(APPEND VTK_NAMES_LIBS "RenderingCore;ParallelCore;IOLegacy;IOCore;zlib;FiltersModeling")
list(APPEND VTK_NAMES_LIBS "FiltersGeometry;FiltersSources;FiltersExtraction;FiltersStatistics")
list(APPEND VTK_NAMES_LIBS "alglib;ImagingFourier;ImagingCore;FiltersGeneral;FiltersCore")
list(APPEND VTK_NAMES_LIBS "CommonExecutionModel;CommonComputationalGeometry;CommonDataModel;CommonTransforms")
list(APPEND VTK_NAMES_LIBS "CommonSystem;CommonMisc;CommonMath;CommonCore;sys")

foreach(ITEM ${VTK_NAMES_LIBS})
  list(APPEND PACKAGE_LIBRARY ${CMAKE_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}vtk${ITEM}-${VTK_VERSION}${CMAKE_STATIC_LIBRARY_SUFFIX})
endforeach()
if(UNIX)
  list(APPEND PACKAGE_LIBRARY dl)
endif()

# set hidden cache variables *_LIBRARY = PACKAGE_LIBRARY, *_INCLUDE and some defaults
# we have own VTK_INCLUDE, set the cache variable before calling set_standard_variables()
set(VTK_INCLUDE_DIR "${CMAKE_BINARY_DIR}/include/vtk-${VTK_VERSION}" CACHE FILEPATH "${PACKAGE_NAME} include dir.")
set_standard_variables()

# we use the cmake manifest
set(DEPS_INSTALL "${CMAKE_BINARY_DIR}")

# don't set compiler flags, we need to change
set_deps_args_default(OFF) 

# VTK 7.1 does not compile with msvc as C++17 -> set back to C++14
if(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
  string(REPLACE "c++17" "c++14" _CXX_FLAGS ${CFSDEPS_CXX_FLAGS})
else()
  set(_CXX_FLAGS ${CFSDEPS_CXX_FLAGS})
endif()

set(DEPS_ARGS ${DEPS_ARGS}
  -DBUILD_SHARED_LIBS:BOOL=OFF
  -DVTK_Group_Rendering:BOOL=OFF
  -DVTK_Group_StandAlone:BOOL=OFF
  -DModule_vtkFiltersParallel:BOOL=ON
  -DModule_vtkFiltersSMP:BOOL=ON
  # enable builds on systems without OpenGL
  -DVTK_RENDERING_BACKEND:STRING=None
  #if we add more file reader we should make this optional depending on USE_ENSIGHT
  -DModule_vtkIOEnSight:BOOL=ON
  -DModule_vtkIOXML:BOOL=ON
  -DModule_vtkhdf5:BOOL=OFF
  -DModule_vtkIONetCDF:BOOL=OFF
  -DModule_vtkIOParallel:BOOL=OFF
  -DvtkexodusII:BOOL=OFF
  -DCMAKE_CXX_FLAGS:STRING=${_CXX_FLAGS})

# --- it follows generic final block for cmake packages with a patch and no postinstall ---

# copy "static" license as we configure this dependency. Check if license is still valid!
file(COPY "${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/license/" DESTINATION "${CMAKE_BINARY_DIR}/license/${PACKAGE_NAME}" )

# Generate ${PACKAGE_NAME}-patch.cmake we use for our external project
generate_patches_script() # sets PATCHES_SCRIPT

# generate package creation script. We get the files from an install_manifest.txt
generate_packing_script_manifest()

# we have no postinstall, so don't call generate_postinstall_script()
assert_unset(POSTINSTALL_SCRIPT)

#dump_depencency_variables()

# do we want to use precompiled and do we already have the package?
if(${CFS_DEPS_PRECOMPILED} AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  # copy files from cache
  create_external_unpack_precompiled()
# if not, build newly and possibly pack the stuff
else()
  # standard cmake build with patch
  create_external_cmake_patched()  

  # new data just built: shall we pack and store as precompiled?
  if(${CFS_DEPS_PRECOMPILED})
    # add custom step to zip a precompiled package to the cache.
    add_external_storage_step()
  endif() 
endif()

# add project to global list of CFSDEPS
set(CFSDEPS ${CFSDEPS} ${PACKAGE_NAME})
