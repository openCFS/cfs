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

# preCICE embeds libxml2 into its shared libprecice.so, so it needs a -fPIC build.
# The lib stays static (a plain C archive embeds cleanly, nothing extra to ship);
# only the configure flags below switch. DEPS_ID gives the -fPIC variant its own
# precompiled cache name so a cached non-PIC zip is not reused.
if(CFS_BUILD_PRECICE AND UNIX)
  set(DEPS_ID "pic")
endif()

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

# --with-pic: see DEPS_ID above. --without-zlib: otherwise undefined gz* refs get
# embedded into the shared libprecice.so (which has no DT_NEEDED libz), breaking
# every consumer's link. preCICE only parses plain XML config files, and cfs
# handles .xml.gz itself via boost::iostreams, not via libxml2.
if(CFS_BUILD_PRECICE AND UNIX)
  set(DEPS_CONFIGURE ${DEPS_CONFIGURE} --with-pic --without-zlib)
endif()

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