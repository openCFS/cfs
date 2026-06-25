# libxml2 is an alternative for xerces. xerces and libxml2 seem to be the only
# xml parser libs able to validate schema and apply schema defaults.
# https://gitlab.gnome.org/GNOME/libxml2
clear_depencency_variables()

# set mandatory variables for the macros in DependencyTools.cmake.
set(PACKAGE_NAME "libxml2")
set(PACKAGE_VER "2.9.4")
set(PACKAGE_FILE "${PACKAGE_NAME}-${PACKAGE_VER}.tar.gz")
set(PACKAGE_MD5 "ae249165c173b1ff386ee8ad676815f5")
set(DEPS_VER "") # set to "-a", "-b", when dependency changed with same PACKAGE_VER. Reset to "" with new PACKAGE_VER.

# the mirrors can point to arbitrary file names. 
set(PACKAGE_MIRRORS 
  "ftp://xmlsoft.org/libxml2/${PACKAGE_FILE}"
  "http://ftp.osuosl.org/pub/blfs/conglomeration/libxml2/${PACKAGE_FILE}")

# add default mirrors to PACKAGE_MIRRORS or replace all with LOCAL_PACKAGE_FILE if we already have it
add_standard_mirrors_or_set_local()

 # we'll disable fortran for ipopt as it is not needed
use_c_and_fortran(ON OFF)

# sets PRECOMPILED_PCKG_FILE to the full precompiled name including path
set_precompiled_pckg_file()

# libxml2 does not become liblibxml2.a but libxml2.a
set(PACKAGE_LIBRARY ${CMAKE_BINARY_DIR}/lib/${PACKAGE_NAME}${CMAKE_STATIC_LIBRARY_SUFFIX})

# set hidden cache variables *_LIBRARY = PACKAGE_LIBRARY, *_INCLUDE and some defaults
set_standard_variables()
# this is the standard target for configure projects (builds in source). This directory will be zipped
set(DEPS_INSTALL "${DEPS_PREFIX}/install")

set_configure_default()
set(DEPS_CONFIGURE ${DEPS_CONFIGURE} --disable-shared --without-ftp --without-html --without-http --without-icu --without-iconv --without-python --without-modules --without-lzma)    

# --- it follows generic final block for cmake packages with a patch and no postinstall ---

# copy "static" license as we configure this dependency. Check if license is still valid!
file(COPY "${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/license/" DESTINATION "${CMAKE_BINARY_DIR}/license/${PACKAGE_NAME}" )

assert_unset(PATCHES_SCRIPT)

# generate package creation script. We do not get the files from an install_manifest.txt
generate_packing_script_install_dir()

# we have no postinstall, so don't call generate_postinstall_script()
assert_unset(POSTINSTALL_SCRIPT)

#dump_depencency_variables()

# do we want to use precompiled and do we already have the package?
if(${CFS_DEPS_PRECOMPILED} AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  # copy files from cache
  create_external_unpack_precompiled()

# if not, build newly and possibly pack the stuff
else()
  # add external project step actually building an cmake package including a patch 
  # also genearate the patch script via generate_patches_script()
  create_external_configure()  

  # new data just built: shall we pack and store as precompiled?
  if(${CFS_DEPS_PRECOMPILED})
    # add custom step to zip a precompiled package to the cache.
    add_external_storage_step()
  else()
    # without manifest (installs directly to binary dir) an without packing, we need to copy manually  
    add_install_dir_to_binary_step()  
  endif()  
endif()

# add project to global list of CFSDEPS
set(CFSDEPS ${CFSDEPS} ${PACKAGE_NAME})

# ---------------------------------------------------------------------------
# preCICE-only, position-independent (-fPIC) copy of the SAME libxml2 version.
# preCICE is a shared library and embeds libxml2 statically, so it needs a -fPIC
# build (autotools: --with-pic). The default build above is left untouched; this
# variant only exists when preCICE is built from source (CFS_BUILD_PRECICE) and
# installs into its own prefix, exposed as LIBXML2_PIC_ROOT for
# cfsdeps/precice/External_PRECICE.cmake (found via FindLibXml2 / LibXml2_ROOT).
# NOTE: not validated in this environment - verify on the build machine.
# ---------------------------------------------------------------------------
if(CFS_BUILD_PRECICE AND UNIX)
  set(LIBXML2_PIC_PREFIX  "${CMAKE_BINARY_DIR}/cfsdeps/libxml2-pic")
  set(LIBXML2_PIC_SRC     "${LIBXML2_PIC_PREFIX}/src/libxml2-pic")
  set(LIBXML2_PIC_INSTALL "${LIBXML2_PIC_PREFIX}/install")
  set(LIBXML2_PIC_ROOT "${LIBXML2_PIC_INSTALL}" CACHE PATH "preCICE-only -fPIC libxml2 prefix")
  mark_as_advanced(LIBXML2_PIC_ROOT)
  set(LIBXML2_PIC_ZIP "${CFS_DEPS_CACHE_DIR}/precompiled/libxml2-pic_${PACKAGE_VER}${DEPS_VER}_${CFS_ARCH_STR}_C-${CMAKE_CXX_COMPILER_ID}-${CMAKE_CXX_COMPILER_VERSION}.tar.gz")

  if(${CFS_DEPS_PRECOMPILED} AND EXISTS "${LIBXML2_PIC_ZIP}")
    # reuse the cached -fPIC libxml2: only unpack it into the install prefix
    ExternalProject_Add(libxml2-pic
      PREFIX ${LIBXML2_PIC_PREFIX}
      DOWNLOAD_COMMAND "" CONFIGURE_COMMAND "" BUILD_COMMAND ""
      INSTALL_COMMAND ${CMAKE_COMMAND} -E make_directory ${LIBXML2_PIC_INSTALL}
      COMMAND ${CMAKE_COMMAND} -E chdir ${LIBXML2_PIC_INSTALL} ${CMAKE_COMMAND} -E tar xzf ${LIBXML2_PIC_ZIP})
  else()
    # in-source autotools build (like the default), but with --with-pic
    ExternalProject_Add(libxml2-pic
      PREFIX ${LIBXML2_PIC_PREFIX}
      SOURCE_DIR ${LIBXML2_PIC_SRC}
      BINARY_DIR ${LIBXML2_PIC_SRC}
      URL ${PACKAGE_MIRRORS}
      URL_MD5 ${PACKAGE_MD5}
      DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/${PACKAGE_NAME}
      DOWNLOAD_NAME ${PACKAGE_FILE}
      DOWNLOAD_NO_PROGRESS ON
      CONFIGURE_COMMAND ./configure --prefix=${LIBXML2_PIC_INSTALL} --disable-shared --with-pic
                        --without-ftp --without-html --without-http --without-icu --without-iconv
                        --without-python --without-modules --without-lzma
      BUILD_COMMAND make
      INSTALL_COMMAND make install
      LOG_BUILD 1)
    if(${CFS_DEPS_PRECOMPILED})
      ExternalProject_Add_Step(libxml2-pic store-precompiled
        COMMAND ${CMAKE_COMMAND} -E make_directory ${CFS_DEPS_CACHE_DIR}/precompiled
        COMMAND ${CMAKE_COMMAND} -E chdir ${LIBXML2_PIC_INSTALL} ${CMAKE_COMMAND} -E tar czf ${LIBXML2_PIC_ZIP} .
        DEPENDEES install)
    endif()
  endif()
  set(CFSDEPS ${CFSDEPS} libxml2-pic)
endif()