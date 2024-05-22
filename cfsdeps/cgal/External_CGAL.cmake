#-------------------------------------------------------------------------------
# Computational Geometry Algorithms Library (CGAL)
# As of April 2024, CGAL is used in CFS to enable enhanced element intersection  
# computations for non-conforming interfaces and grid interpolation techniques. 
# It is also used for 3D geometry approximations on surface nodes in the 
# layerGeneration feature and the curvilinear PML formulation.
#
# CGAL is only used as a supplementary build option in CFS. It is not part of the 
# default/distributed build, due to its GPL3 licence.
#
# Project Homepage
# http://www.cgal.org/
#
# Description copied from the CGAL project homepage:
# "CGAL is an open source software project that provides easy access to efficient 
# and reliable geometric algorithms in the form of a C++ library. CGAL is used in
# various areas needing geometric computation, such as geographic information systems, 
# computer aided design, molecular biology, medical imaging, computer graphics, and robotics.
# The library offers data structures and algorithms like triangulations, Voronoi diagrams, 
# Boolean operations on polygons and polyhedra, point set processing, arrangements of curves, 
# surface and volume mesh generation, geometry processing, alpha shapes, convex hull algorithms,
# shape reconstruction, AABB and KD trees... Explore the complete list of features and 
# capabilities by visiting the CGAL Package Overview.
# The CGAL data structures and algorithms are distributed under a dual license, namely under 
# the GPL v3+ and, alternatively, under a commercial license by GeometryFactory."
#-------------------------------------------------------------------------------

# make sure not to uninetendently use another packages settings. Supports assert_set() checks. Is mandatory!
clear_depencency_variables()

# set mandatory variables for the macros in DependencyTools.cmake.
set(CFS_VERSION_MONTH "04")
set(PACKAGE_NAME "cgal")
set(PACKAGE_VER "5.6.1")
set(CGAL_VER ${PACKAGE_VER}) # for Dependencies.cc
set(PACKAGE_FILE "CGAL-${PACKAGE_VER}.tar.xz")
set(PACKAGE_MD5 "90ccd8f68894ab2f89933163afd66535")
set(DEPS_VER "") # set to "-a", "-b", when dependency changed with same PACKAGE_VER. Reset to "" with new PACKAGE_VER.

# the mirrors can point to arbitrary file names. 
set(PACKAGE_MIRRORS "https://github.com/CGAL/cgal/releases/download/v${PACKAGE_VER}/${PACKAGE_FILE}") 

# add default mirrors to PACKAGE_MIRRORS or replace all with LOCAL_PACKAGE_FILE if we already have it
add_standard_mirrors_or_set_local()

# we do not need to compile as CGAL is now header only
use_c_and_fortran(OFF OFF)

# sets PRECOMPILED_PCKG_FILE to the full precompiled name including path
set_precompiled_pckg_file()

# set hidden cache variables *_LIBRARY = PACKAGE_LIBRARY, *_INCLUDE and some defaults
set_header_only_standard_variables()

# this is the standard target for cmake projects. The files to package come from the install_manifest.txt
set(DEPS_INSTALL "${CMAKE_BINARY_DIR}")

# set DEPS_ARG with defaults for a cmake project
set_deps_args_default(ON)

# set custom build arguments
set(DEPS_ARGS 
  ${DEPS_ARGS}
  -DWITH_EXAMPLES:BOOL=OFF
  -DWITH_DEMOS:BOOL=OFF
  -DWITH_CGAL_Core:BOOL=ON
  -DWITH_CGAL_Qt5:BOOL=OFF
  -DWITH_CGAL_ImageIO:BOOL=OFF
  -DCGAL_DIR:STRING="${CMAKE_BINARY_DIR}/include/"
)

# --- it follows generic final block for cmake packages with no patch and no postinstall ---

# copy "static" license as we configure this dependency. Check if license is still valid!
file(COPY "${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/license/" DESTINATION "${CMAKE_BINARY_DIR}/license/${PACKAGE_NAME}" )

assert_unset(PATCHES_SCRIPT)

# generate package creation script. We get the files from an install_manifest.txt
generate_packing_script_manifest()

# we have no postinstall, so don't call generate_postinstall_script()
assert_unset(POSTINSTALL_SCRIPT)

# for debugging purposes, dump all dependency variables
dump_depencency_variables()

# do we want to use precompiled and do we already have the package?
if(${CFS_DEPS_PRECOMPILED} AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  # copy files from cache
  create_external_unpack_precompiled()

# if not, build newly and possibly pack the stuff
else()
  create_external_cmake()

  # new data just built: shall we pack and store as precompiled?
  if(${CFS_DEPS_PRECOMPILED})
    # add custom step to zip a precompiled package to the cache.
    add_external_storage_step()
  else()
    # without manifest (installs directly to binary dir) an without packing, we need to copy manually  
    add_install_dir_to_binary_step()  
  endif()  
endif()

# add dependencies for cgal
add_dependencies(cgal boost gmp mpfr)

# add project to global list of CFSDEPS
set(CFSDEPS ${CFSDEPS} ${PACKAGE_NAME})