CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_build_type_options.hh.in"
  "${CFS_BINARY_DIR}/include/def_build_type_options.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_arpack.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_arpack.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_ilupack.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_ilupack.hh")

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


CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_xerces.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_xerces.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_scripting.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_scripting.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_tcl.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_tcl.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_python.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_python.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_mesh.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_mesh.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_gidpost.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_gidpost.hh")

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_hdf5.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_hdf5.hh")

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

IF(USE_MPCCI)
  SET(MpCCI 1)
ENDIF(USE_MPCCI)

CONFIGURE_FILE("${CFS_SOURCE_DIR}/include/def_use_mpcci.hh.in"
  "${CFS_BINARY_DIR}/include/def_use_mpcci.hh")
