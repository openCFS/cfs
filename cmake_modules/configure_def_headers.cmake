#-----------------------------------------------------------------------------
# In this CMake file all headers get configured with the current values of
# the CMake variables. The configured headers will be written to
# CFS_BUILD_DIR/include. To see what values are written to which file just
# open one of the *.hh.in header templates.
#-----------------------------------------------------------------------------

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_build_type_options.hh.in"
  "${CFS_BINARY_DIR}/include/def_build_type_options.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_config.hh.in"
  "${CFS_BINARY_DIR}/include/def_config.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_reordering.hh.in"
  "${CFS_BINARY_DIR}/include/def_reordering.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_arpack.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_arpack.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_eigen.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_eigen.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_build_quadeigensolver.hh.in"
  "${CFS_BINARY_DIR}/include/def_build_quadeigensolver.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_lis.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_lis.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_ginkgo.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_ginkgo.hh")
  
CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_petsc.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_petsc.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_hwloc.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_hwloc.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_ghost.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_ghost.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_phist_cg.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_phist_cg.hh")
  
CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_phist_ev.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_phist_ev.hh")


CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_mpi.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_mpi.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_suitesparse.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_suitesparse.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_superlu.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_superlu.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_blas.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_blas.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_pardiso.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_pardiso.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_metis.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_metis.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_ipopt.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_ipopt.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_scpip.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_scpip.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_dumas.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_dumas.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_sgp.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_sgp.hh")

configure_file("${CFS_SOURCE_DIR}/include/def_use_embedded_python.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_embedded_python.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_snopt.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_snopt.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_cgal.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_cgal.hh")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_nanoflann.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_nanoflann.hh")

IF(USE_LIBFBI)
  SET(USE_LIBFBI 1)
ELSE()
  SET(USE_LIBFBI 0)
ENDIF()
CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_libfbi.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_libfbi.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_xerces.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_xerces.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_libxml2.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_libxml2.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_feast.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_feast.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_gidpost.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_gidpost.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_ensight.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_ensight.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_cgns.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_cgns.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_build_quadeigensolver.hh.in"
  "${CFS_BINARY_DIR}/include/def_build_quadeigensolver.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_unv.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_unv.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_openmp.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_openmp.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_cuda.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_cuda.hh")


CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_sgpp.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_sgpp.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_xmlschema.hh.in"
  "${CFS_BINARY_DIR}/include/def_xmlschema.hh"
  @ONLY )

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_cfs_stats.hh.in"
  "${CFS_BINARY_DIR}/include/def_cfs_stats.hh"
  @ONLY )

