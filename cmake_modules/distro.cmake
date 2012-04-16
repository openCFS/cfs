SET(DISTRO_SCRIPT "${CFS_SOURCE_DIR}/share/scripts/distro.sh")

EXEC_PROGRAM("${DISTRO_SCRIPT}"
  ARGS -c
  OUTPUT_VARIABLE CFS_ARCH_TEST
  RETURN_VALUE RETVAL)

LIST(GET CFS_ARCH_TEST 0 CFS_DISTRO)
LIST(GET CFS_ARCH_TEST 1 CFS_DISTRO_VER)
LIST(GET CFS_ARCH_TEST 2 CFS_ARCH)
# Determine the subarchitecture of the platform.
LIST(GET CFS_ARCH_TEST 3 CFS_SUBARCH)

# MESSAGE("CFS_ARCH_TEST: ${CFS_ARCH_TEST}")

# MESSAGE("CFS_DISTRO: ${CFS_DISTRO}")
# MESSAGE("CFS_DISTRO_VER: ${CFS_DISTRO_VER}")
# MESSAGE("CFS_ARCH: ${CFS_ARCH}")
# MESSAGE("CFS_ARCH_STR: ${CFS_ARCH_STR}")
# MESSAGE("CFS_SUBARCH: ${CFS_SUBARCH}")

# Set a few distribution specific variables (FORTRAN libs, ...)
# G2C_LIBRARY and GFORTRAN_LIBRARY are defined in
# cmake_modules/FindFortranLibs.cmake
IF(CFS_DISTRO STREQUAL "SUSE" OR
   CFS_DISTRO STREQUAL "OPENSUSE" OR
   CFS_DISTRO STREQUAL "SLE")

 # Supply a list of SUSE system packages which are needed to compile CFS
 SET(CFS_PCKG_HINT
   "blas"
   "boost(-devel)"
   "compat-g77"
   "lapack"
   "python(-devel)"
   "subversion"
   "tcl(-devel)"
   "Xerces-C(-devel)"
   "doxygen"
   "graphviz"
   )

 IF(CFS_DISTRO_VER GREATER 9.3)
   SET(CFS_PCKG_HINT ${CFS_PCKG_HINT} gfortran)
   SET(CFS_FORTRAN_LIBS ${CFS_FORTRAN_LIBS} ${GFORTRAN_LIBRARY})   
 ENDIF(CFS_DISTRO_VER GREATER 9.3)
 
ENDIF(CFS_DISTRO STREQUAL "SUSE" OR
  CFS_DISTRO STREQUAL "OPENSUSE" OR
  CFS_DISTRO STREQUAL "SLE")

IF(CFS_DISTRO STREQUAL "DEBIAN")
  SET(CFS_PCKG_HINT
   "blas-dev, refblas3-dev"
   "libg2c0-dev"
   "libboost-*-dev"
   "lapack-dev, lapack99-dev, lapack3-dev"
   "python2.4-dev"
   "tcl8.4-dev"
   "subversion"
   "libxerces25-dev, libxerces26-dev, libxerces27-dev"
   "doxygen"
   "tetex-bin"
   "graphviz"
   )

  IF(CFS_DISTRO_VER GREATER 3)
    SET(CFS_FORTRAN_LIBS ${GFORTRAN_LIBRARY})
    SET(CFS_PCKG_HINT ${CFS_PCKG_HINT} "gfortran, libgfortran1-dev")
  ELSE(CFS_DISTRO_VER GREATER 3)
    SET(CFS_FORTRAN_LIBS ${G2C_LIBRARY})
  ENDIF(CFS_DISTRO_VER GREATER 3)

ENDIF(CFS_DISTRO STREQUAL "DEBIAN")


IF(CFS_DISTRO STREQUAL "UBUNTU")

  # Supply a list of Ubuntu system packages which are needed to compile CFS
  SET(CFS_PCKG_HINT
   "refblas3-dev"
   "libboost-*-dev"
   "libg2c0-dev"
   "lapack3-dev"
   "python2.4-dev"
   "subversion"
   "tcl8.4-dev"
   "libxerces27-dev"
   "doxygen"
   "tetex-bin"
   "graphviz"
   )

  # The CFS_FORTRAN_LIBS variable has been set empirically.
  # At the moment I just have Edgy Eft (6.10)!
  IF(CFS_DISTRO_VER GREATER 6)
    SET(CFS_FORTRAN_LIBS ${GFORTRAN_LIBRARY})
    SET(CFS_PCKG_HINT ${CFS_PCKG_HINT} "libgfortran1-dev")
  ENDIF(CFS_DISTRO_VER GREATER 6)

ENDIF(CFS_DISTRO STREQUAL "UBUNTU")

IF(CFS_DISTRO STREQUAL "FEDORA" OR
   CFS_DISTRO STREQUAL "REDHAT" OR
   CFS_DISTRO STREQUAL "RHEL" OR
   CFS_DISTRO STREQUAL "CENTOS")

  SET(CFS_PCKG_HINT
   "blas, libblas.so.3"
   "boost, boost-devel"
   "compat-gcc-g77, compat-gcc-34-g77, compat-gcc-32-g77, gcc-g77"
   "gcc-gfortran"
   "lapack, lapack-devel"
   "python, python-devel"
   "tcl, tcl-devel"
   "subversion"
   "xerces-c, xerces-c-devel"
   "doxygen"
   "graphviz"
   "tetex-latex"
   )

   # Fortran support for Fedora/Redhat/RHEL never tested
   # TODO: add a sane check here
   SET(CFS_FORTRAN_LIBS ${GFORTRAN_LIBRARY})

ENDIF(CFS_DISTRO STREQUAL "FEDORA" OR
   CFS_DISTRO STREQUAL "REDHAT" OR
   CFS_DISTRO STREQUAL "RHEL" OR
   CFS_DISTRO STREQUAL "CENTOS")

IF(CFS_DISTRO STREQUAL "MANDRAKE" OR
   CFS_DISTRO STREQUAL "MANDRIVA")

  SET(CFS_PCKG_HINT
   "libblas1.1-devel libblas3-devel, libblas-devel, blas-devel"
   "libboost1, libboost, boost, libboost1-devel, libboost-devel boost-devel"
   "gcc3.3-g77"
   "libgfortran0, libgfortran, libgfortran4.0 libgfortran.so.0"
   "lapack"
   "libpython2.5-devel, python-devel, libpython-devel"
   "libtcl8.4-devel, tcl-devel, libtcl-devel"
   "subversion"
   "libxerces-c26-devel, xerces-c-devel, libxerces-c-devel"
   "doxygen"
   "tetex-latex"
   "graphviz"
   )

   # Fortran support for Mandrake/Mandriva never tested
   # TODO: add a sane check here
   SET(CFS_FORTRAN_LIBS ${G2C_LIBRARY} ${GFORTRAN_LIBRARY})

ENDIF(CFS_DISTRO STREQUAL "MANDRAKE" OR
   CFS_DISTRO STREQUAL "MANDRIVA")

IF(CFS_DISTRO STREQUAL "MACOSX")

  IF(NOT ${CFS_FORTRAN_COMPILER_NAME} STREQUAL "IFORT")
    SET(CFS_FORTRAN_LIBS ${GFORTRAN_LIBRARY})
  ENDIF(NOT ${CFS_FORTRAN_COMPILER_NAME} STREQUAL "IFORT")

  IF(CMAKE_OSX_ARCHITECTURES)
    LIST(LENGTH CMAKE_OSX_ARCHITECTURES NUM_ARCH)

    IF(NUM_ARCH GREATER 1)
      MESSAGE(FATAL_ERROR "Only binaries for one Mac OS X architecture can be built at a time!")
    ENDIF(NUM_ARCH GREATER 1)

    LIST(GET CMAKE_OSX_ARCHITECTURES 0 OSX_ARCH)
    # MESSAGE("OSX_ARCH ${OSX_ARCH}")

    IF(OSX_ARCH STREQUAL "i386")
      SET(CFS_ARCH "I386")
    ELSE(OSX_ARCH STREQUAL "i386")
      IF(OSX_ARCH STREQUAL "x86_64")
	SET(CFS_ARCH "X86_64")
      ELSE(OSX_ARCH STREQUAL "x86_64")
	MESSAGE(FATAL_ERROR "Building for architecture ${OSX_ARCH} not supported!")
      ENDIF(OSX_ARCH STREQUAL "x86_64")
    ENDIF(OSX_ARCH STREQUAL "i386")

  ELSE(CMAKE_OSX_ARCHITECTURES)
    IF(CFS_DISTRO_VER GREATER 10.5)
      SET(CMAKE_OSX_ARCHITECTURES "x86_64")
    ELSE(CFS_DISTRO_VER GREATER 10.5)
      SET(CMAKE_OSX_ARCHITECTURES "i386")
    ENDIF(CFS_DISTRO_VER GREATER 10.5)
  ENDIF(CMAKE_OSX_ARCHITECTURES)

 
ENDIF(CFS_DISTRO STREQUAL "MACOSX")

# Macro which prints a list of libraries which are needed to compile CFS++
MACRO(SUGGEST_INSTALL_PCKG)
  MESSAGE("Either try installing the following packages to compile CFS++:")
  FOREACH(PCKG ${CFS_PCKG_HINT})
    MESSAGE(${PCKG})
  ENDFOREACH(PCKG)
  MESSAGE("or try to set the corresponding advanced variables.")
ENDMACRO(SUGGEST_INSTALL_PCKG)


IF(CFS_ARCH STREQUAL "I386")
  # Set path suffix for system libs
  SET(LIB_SUFFIX "lib")
ENDIF(CFS_ARCH STREQUAL "I386")

IF(CFS_ARCH STREQUAL "X86_64")
  # Set path suffix for system libs
  SET(LIB_SUFFIX "lib64")
ENDIF(CFS_ARCH STREQUAL "X86_64")


IF(CFS_ARCH STREQUAL "IA64")
  # Set path suffix for system libs
  SET(LIB_SUFFIX "lib")
ENDIF(CFS_ARCH STREQUAL "IA64")

# Determine the system in bits
IF(CFS_ARCH STREQUAL "I386")
  SET(SYSTEM_BITS "32")
ELSE(CFS_ARCH STREQUAL "I386")
  SET(SYSTEM_BITS "64")
ENDIF(CFS_ARCH STREQUAL "I386")

SET(SYSTEM_BITS "${SYSTEM_BITS}"
  CACHE INTERNAL "Determins if this is a '32' or '64' bit system.")
SET(CFS_DISTRO "${CFS_DISTRO}"
  CACHE INTERNAL "String specifying the distribution CFS is built on.")
SET(CFS_DISTRO_VER "${CFS_DISTRO_VER}"
  CACHE INTERNAL "Version of the distribution CFS is built on.")
SET(CFS_ARCH "${CFS_ARCH}"
  CACHE INTERNAL "Architecture CFS is built on.")
SET(CFS_ARCH_STR "${CFS_DISTRO}_${CFS_DISTRO_VER}_${CFS_ARCH}"
  CACHE INTERNAL "String specifying the architecture CFS is built on.")
SET(CFS_SUBARCH "${CFS_SUBARCH}"
  CACHE INTERNAL "Sub-Architecture CFS is built on.")
