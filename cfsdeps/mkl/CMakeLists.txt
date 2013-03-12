PROJECT(wrap-chkstk CXX)

IF(MINGW)
  IF(CFS_BLAS_LAPACK STREQUAL MKL)
    ADD_LIBRARY(wrap-chkstk STATIC wrap_chkstk.cc)
  ENDIF()
ENDIF()
