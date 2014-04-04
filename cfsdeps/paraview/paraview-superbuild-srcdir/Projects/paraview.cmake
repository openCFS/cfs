set(PV_EXTRA_CMAKE_ARGS ""
    CACHE STRING "Extra arguments to be passed to ParaView when configuring.")
mark_as_advanced(PV_EXTRA_CMAKE_ARGS)

# Enable CFS++ readers
set (extra_cmake_args
  -DModule_vtkIOOFFReader:BOOL=ON
  -DModule_vtkIOCFSReader:BOOL=ON
)
if (manta_ENABLED)
  list (APPEND extra_cmake_args
    -DMANTA_BUILD:PATH=${SuperBuild_BINARY_DIR}/manta/src/manta-build)
endif()
if(PV_NIGHTLY_SUFFIX)
  list (APPEND extra_cmake_args
    -DPV_NIGHTLY_SUFFIX:STRING=${PV_NIGHTLY_SUFFIX})
endif()

if (APPLE)
  # We are having issues building mpi4py with Python 2.6 on Mac OSX. Hence,
  # disable it for now.
  list (APPEND extra_cmake_args
        -DPARAVIEW_USE_SYSTEM_MPI4PY:BOOL=OFF)
endif()

set (PARAVIEW_INSTALL_DEVELOPMENT_FILES FALSE)
if (paraviewsdk_ENABLED)
  set (PARAVIEW_INSTALL_DEVELOPMENT_FILES TRUE)
endif()

#-------------------------------------------------------------------------------
# Set names of patch file and template file.
#-------------------------------------------------------------------------------
SET(PFN_TEMPL "${CMAKE_CURRENT_SOURCE_DIR}/Projects/paraview.patch.cmake.in")
SET(PFN "${SuperBuild_BINARY_DIR}/paraview/paraview-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY)

add_external_project(paraview
  DEPENDS_OPTIONAL
    # boost hdf5 qt zlib cgns
    cosmotools ffmpeg libxml3 manta matplotlib mpi numpy png python visitbridge silo
    mesa osmesa nektarreader
    ${PV_EXTERNAL_PROJECTS}

  PATCH_COMMAND 
    "${CMAKE_COMMAND}" -P "${PFN}"

  CMAKE_ARGS
    -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
    -DCMAKE_MAKE_PROGRAM:FILEPATH=${CMAKE_MAKE_PROGRAM}
    -DBUILD_SHARED_LIBS:BOOL=ON
    -DBUILD_TESTING:BOOL=OFF
    -DPARAVIEW_BUILD_PLUGIN_CoProcessingScriptGenerator:BOOL=ON
    -DPARAVIEW_BUILD_PLUGIN_EyeDomeLighting:BOOL=ON
    -DPARAVIEW_BUILD_PLUGIN_MantaView:BOOL=${manta_ENABLED}
    -DPARAVIEW_BUILD_QT_GUI:BOOL=${qt_ENABLED}
    -DQT_QMAKE_EXECUTABLE:FILEPATH=${QT_QMAKE_EXECUTABLE}
    -DPARAVIEW_ENABLE_FFMPEG:BOOL=${ffmpeg_ENABLED}
    -DPARAVIEW_ENABLE_PYTHON:BOOL=${python_ENABLED}
    -DPARAVIEW_ENABLE_COSMOTOOLS:BOOL=${cosmotools_ENABLED}
    -DPARAVIEW_USE_MPI:BOOL=${mpi_ENABLED}
    -DPARAVIEW_USE_VISITBRIDGE:BOOL=${visitbridge_ENABLED}
    -DVISIT_BUILD_READER_CGNS:BOOL=ON
    -DVISIT_BUILD_READER_Silo:BOOL=${silo_ENABLED}
    -DPARAVIEW_INSTALL_DEVELOPMENT_FILES:BOOL=${PARAVIEW_INSTALL_DEVELOPMENT_FILES}
    -DCMAKE_Fortran_COMPILER:FILEPATH=${CMAKE_Fortran_COMPILER}
    -DFFMPEG_INCLUDE_DIR:PATH=${SuperBuild_BINARY_DIR}/install/include
    # since VTK mangles all the following, I wonder if there's any point in
    # making it use system versions.
#    -DVTK_USE_SYSTEM_FREETYPE:BOOL=${ENABLE_FREETYPE}
#    -DVTK_USE_SYSTEM_LIBXML2:BOOL=${ENABLE_LIBXML2}
#    -DVTK_USE_SYSTEM_PNG:BOOL=${ENABLE_PNG}
    -DVTK_USE_SYSTEM_ZLIB:BOOL=ON
    -DVTK_USE_SYSTEM_HDF5:BOOL=ON
    -DBoost_DIR:PATH=${Boost_DIR}
    -DBoost_INCLUDE_DIR:PATH=${Boost_INCLUDE_DIR}
    -DBoost_LIBRARY_DIRS:STRING=${Boost_LIBRARY_DIRS}
    -DCGNS_INCLUDE_DIR:PATH=${CGNS_INCLUDE_DIR}
    -DCGNS_LIBRARY:FILEPATH=${CGNS_LIBRARY}
    -DHDF5_CXX_LIBRARY:FILEPATH=${HDF5_CXX_LIBRARY}
    -DHDF5_C_LIBRARY:FILEPATH=${HDF5_C_LIBRARY}
    -DHDF5_DIR:PATH=${HDF5_DIR}
    -DHDF5_HL_LIBRARY:FILEPATH=${HDF5_HL_LIBRARY}
    -DZLIB_INCLUDE_DIR:PATH=${ZLIB_INCLUDE_DIR}
    -DZLIB_LIBRARY:FILEPATH=${ZLIB_LIBRARY}

    # Web documentation
    -DPARAVIEW_BUILD_WEB_DOCUMENTATION:BOOL=${PARAVIEW_BUILD_WEB_DOCUMENTATION}

    # specify the apple app install prefix. No harm in specifying it for all
    # platforms.
    -DMACOSX_APP_INSTALL_PREFIX:PATH=<INSTALL_DIR>/Applications

  ${extra_cmake_args}

  ${PV_EXTRA_CMAKE_ARGS}
)
