#-----------------------------------------------------------------------------
# In this CMake file all headers get configured with the current values of
# the CMake variables. The configured headers will be written to
# CFS_BUILD_DIR/include. To see what values are written to which file just
# open one of the *.hh.in header templates.
#-----------------------------------------------------------------------------

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_build_type_options.hh.in"
  "${CFS_BINARY_DIR}/include/def_build_type_options.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_arpack.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_arpack.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_ilupack.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_ilupack.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_cholmod.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_cholmod.hh")

IF(CFS_BLAS_LAPACK STREQUAL "MKL")
  SET(USE_MKL 1)
ENDIF(CFS_BLAS_LAPACK STREQUAL "MKL")

IF(CFS_BLAS_LAPACK STREQUAL "ACML")
  SET(USE_ACML 1)
ENDIF(CFS_BLAS_LAPACK STREQUAL "ACML")

IF(CFS_BLAS_LAPACK STREQUAL "GOTO")
  SET(USE_GOTO 1)
ENDIF(CFS_BLAS_LAPACK STREQUAL "GOTO")

IF(CFS_BLAS_LAPACK STREQUAL "NETLIB")
  SET(USE_NETLIB 1)
ENDIF(CFS_BLAS_LAPACK STREQUAL "NETLIB")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_blas.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_blas.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_lapack.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_lapack.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_pardiso.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_pardiso.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_metis.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_metis.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_ipopt.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_ipopt.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_scpip.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_scpip.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_snopt.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_snopt.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_knitro.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_knitro.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_interpolation.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_interpolation.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_xerces.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_xerces.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_mesh.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_mesh.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_gidpost.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_gidpost.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_hdf5.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_hdf5.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_gmsh.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_gmsh.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_gmv.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_gmv.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_unv.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_unv.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_ansysrst.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_ansysrst.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_xmlschema.hh.in"
  "${CFS_BINARY_DIR}/include/def_xmlschema.hh"
  @ONLY )

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_cfs_stats.hh.in"
  "${CFS_BINARY_DIR}/include/def_cfs_stats.hh"
  @ONLY )

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_expl_templ_inst.hh.in"
  "${CFS_BINARY_DIR}/include/def_expl_templ_inst.hh"
  @ONLY )
  
CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_cplreader.hh.in"
  "${CFS_BINARY_DIR}/include/def_cplreader.hh")
