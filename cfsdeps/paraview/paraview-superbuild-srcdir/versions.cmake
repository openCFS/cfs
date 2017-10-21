# This maintains the links for all sources used by this superbuild.
# Simply update this file to change the revision.
# One can use different revision on different platforms.
# e.g.
# if (UNIX)
#   ..
# else (APPLE)
#   ..
# endif()

include(CMakeDependentOption)

add_revision(bzip2
  URL "http://www.paraview.org/files/dependencies/bzip2-1.0.6.tar.gz"
  URL_MD5 00b516f4704d4a7cb50a1d97e6e8e15b)
# NOTE: if updating bzip2 version, fix patch in bzip2.cmake


add_revision(zlib
  URL "http://www.paraview.org/files/dependencies/zlib-1.2.7.tar.gz"
  URL_MD5 60df6a37c56e7c1366cca812414f7b85
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/zlib)
# NOTE: if updating zlib version, fix patch in zlib.cmake


add_revision(png
  URL "http://paraview.org/files/dependencies/libpng-1.4.8.tar.gz"
  URL_MD5 49c6e05be5fa88ed815945d7ca7d4aa9
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/paraview)

add_revision(freetype
  URL "http://paraview.org/files/dependencies/freetype-2.4.8.tar.gz"
  URL_MD5 "5d82aaa9a4abc0ebbd592783208d9c76"
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/paraview)

add_revision(szip
  URL "http://paraview.org/files/dependencies/szip-2.1.tar.gz"
  URL_MD5 902f831bcefb69c6b635374424acbead
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/szip)

add_revision(hdf5
  URL "http://www.paraview.org/files/dependencies/hdf5-1.8.13.tar.gz"
  URL_MD5 c03426e9e77d7766944654280b467289
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/hdf5)

add_revision(silo
  URL "http://paraview.org/files/dependencies/silo-4.9.1-bsd.tar.gz"
  URL_MD5 465d2a0a8958b088cde83fb2a5a7eeef
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/paraview)

add_revision(cgns
  URL "http://www.paraview.org/files/dependencies/cgnslib_3.1.3-4.tar.gz"
  URL_MD5 442bba32b576f3429cbd086af43fd4ae
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/cgns)

add_revision(ffmpeg
  URL "http://paraview.org/files/dependencies/ffmpeg-2.3.3.tar.bz2"
  URL_MD5 72361d3b8717b6db3ad2b9da8df7af5e
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/paraview)

add_revision(libxml2
  URL "http://paraview.org/files/dependencies/libxml2-2.7.8.tar.gz"
  URL_MD5 8127a65e8c3b08856093099b52599c86
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/paraview)

add_revision(fontconfig
  URL "http://paraview.org/files/dependencies/fontconfig-2.8.0.tar.gz"
  URL_MD5 77e15a92006ddc2adbb06f840d591c0e
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/paraview)

add_revision(qt4
  URL "http://paraview.org/files/dependencies/qt-everywhere-opensource-src-4.8.6.tar.gz"
  URL_MD5 2edbe4d6c2eff33ef91732602f3518eb
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/qt4)

if (WIN32 OR (CROSS_BUILD_STAGE STREQUAL "CROSS"))
  add_revision(qt5
    URL "http://download.qt.io/official_releases/qt/5.4/5.4.1/single/qt-everywhere-opensource-src-5.4.1.zip"
    URL_MD5 57b25c68982237abb9e18b347034e005)
else ()
  add_revision(qt5
    URL "http://download.qt.io/official_releases/qt/5.4/5.4.1/single/qt-everywhere-opensource-src-5.4.1.tar.gz"
    URL_MD5 90f3fbce38ed16e5dc2cd0909ae86ca4)
endif ()

if (WIN32 OR (CROSS_BUILD_STAGE STREQUAL "CROSS"))
  if (64bit_build)
    add_revision(python
      URL "http://www.paraview.org/files/dependencies/python+deps.tar.bz2"
      URL_MD5 "4318b8f771eda5606d9ce7f0be9f82e1")
  else ()
    add_revision(python
      URL "http://www.paraview.org/files/dependencies/python+deps-x32.tar.bz2"
      URL_MD5 "6ba441784a672e08379d23ddd61146f0")
  endif ()
elseif (CROSS_BUILD_STAGE STREQUAL "CROSS")
  add_revision(python
    URL "http://www.paraview.org/files/dependencies/Python-2.7.3.tgz"
    URL_MD5 "2cf641732ac23b18d139be077bd906cd"
    DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/python)
else()
  add_revision(python
    URL "http://paraview.org/files/dependencies/Python-2.7.2.tgz"
    URL_MD5 "0ddfe265f1b3d0a8c2459f5bf66894c7"
    DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/python)
endif()

add_revision(numpy
  URL "http://paraview.org/files/dependencies/numpy-1.8.1+cmake+static.tar.bz2"
  URL_MD5 1974dbb4bfa1509e492791a8cd225774
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/python)

add_revision(matplotlib
  URL "http://paraview.org/files/dependencies/matplotlib-1.1.1_notests.tar.gz"
  URL_MD5 30ee59119599331bf1f3b6e838fee9a8
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/python)

add_revision(boost
  URL "http://sourceforge.net/projects/boost/files/boost/1.56.0/boost_1_56_0.tar.gz/download"
  URL_MD5 8c54705c424513fa2be0042696a3a162
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/boost)

add_revision(manta
  URL "http://paraview.org/files/dependencies/manta-r2439.tar.gz"
  URL_MD5 fbf4107fe2f6d7e8a5ae3dda71805bdc
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/paraview)

if (UNIX)
  add_revision(mpi
    URL "http://paraview.org/files/dependencies/mpich2-1.4.1p1.tar.gz"
    URL_MD5 b470666749bcb4a0449a072a18e2c204
    DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/paraview)
elseif (WIN32)
  add_revision(mpi
    URL "http://www.paraview.org/files/dependencies/openmpi-1.4.4.tar.gz"
    URL_MD5 7253c2a43445fbce2bf4f1dfbac113ad
    DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/paraview)
endif()

if (CROSS_BUILD_STAGE STREQUAL "CROSS")
  add_revision(mesa
    URL "http://www.paraview.org/files/dependencies/MesaLib-7.6.1.tar.gz"
    URL_MD5 e80fabad2e3eb7990adae773d6aeacba
    DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/paraview)
else()
  add_revision(mesa
    URL "http://paraview.org/files/dependencies/MesaLib-7.11.2.tar.gz"
    URL_MD5 b9e84efee3931c0acbccd1bb5a860554
    DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/paraview)
endif()

# We stick with 7.11.2 for Mesa version for now. Newer mesa doesn't seem to
# build correctly with certain older compilers (e.g. on neser).
add_revision(osmesa
    URL "http://paraview.org/files/dependencies/MesaLib-7.11.2.tar.gz"
    URL_MD5 b9e84efee3931c0acbccd1bb5a860554  
    DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/paraview)


# ----------------------------------------------------------------------------
# You choose to download ParaView source form GIT or other URL/FILE tarball
option(ParaView_FROM_GIT "If enabled then the repository is fetched from git" ON)
cmake_dependent_option(ParaView_FROM_SOURCE_DIR "Enable to use existing ParaView source." OFF
  "NOT ParaView_FROM_GIT" OFF)

if (ParaView_FROM_GIT)
  # Download PV from GIT
  add_customizable_revision(paraview
    GIT_REPOSITORY https://gitlab.kitware.com/paraview/paraview.git
    GIT_TAG "v4.4.0"
    DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/paraview)
else()
  if (ParaView_FROM_SOURCE_DIR)
    add_customizable_revision(paraview
      SOURCE_DIR "ParaViewSource")
  else()
    # Variables to hold the URL and MD5 (optional)
    add_customizable_revision(paraview
      URL "http://www.paraview.org/files/v4.4/ParaView-v4.4.0-source.tar.gz"
      URL_MD5 "fa1569857dd680ebb4d7ff89c2227378"
      DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/paraview)
  endif()
endif()

#------------------------------------------------------------------------------
# Optional Plugins. Doesn't affect ParaView binaries at all even if missing
# or disabled.
#------------------------------------------------------------------------------

if (USE_NONFREE_COMPONENTS)
  add_revision(qhull
    GIT_REPOSITORY https://github.com/mathstuf/qhull.git
    GIT_TAG next)

  add_revision(genericio
    GIT_REPOSITORY https://kwgitlab.kitware.com/paraview/genericio.git
    GIT_TAG master)

  # Add an option to not use diy from SVN. On Debian-Etch the SVN is too old
  # to work with invalid SVN certificates.
  option(DIY_SKIP_SVN "If enabled, we simply download diy from a source tar" OFF)
  if(DIY_SKIP_SVN)
    add_revision(diy
      URL "http://paraview.org/files/dependencies/diy-src.r178.tar.gz"
      URL_MD5 4fba13aae93927d0f32dd6db0599ffcd)
  else()
    add_revision(diy
      GIT_REPOSITORY https://gitlab.kitware.com/paraview/diy.git
      GIT_TAG f5c00a034279d20cf040705d78a6b7bdb4beae43) # r178
  endif()

  add_revision(cosmotools
    GIT_REPOSITORY git://public.kitware.com/cosmotools.git
    GIT_TAG v0.13)

  add_revision(acusolve
    GIT_REPOSITORY https://kwgitlab.kitware.com/paraview/acusolvereaderplugin.git
    GIT_TAG master)

  add_revision(vistrails
    GIT_REPOSITORY https://kwgitlab.kitware.com/paraview/vistrails.git
    GIT_TAG master)
endif ()

#add_customizable_revision(vortexfinder2
#  GIT_REPOSITORY https://github.com/hguo/vortexfinder2.git
#  GIT_TAG 2bdae9bfc0f36e1013a4b41c5d25c9e6ebbf1701)
add_customizable_revision(vortexfinder2
  GIT_REPOSITORY https://github.com/tj-corona/vortexfinder2.git
  GIT_TAG paraview-superbuild)

#add_revision(mili_plugin
# URL ${CMAKE_CURRENT_SOURCE_DIR}/Externals/mili)

add_revision(portfwd
  URL "http://www.paraview.org/files/dependencies/portfwd-0.29.tar.gz"
  URL_MD5 93161c91e12b0d67ca52dc13708a2f2f
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/paraview)

add_revision(lapack
  URL "http://paraview.org/files/dependencies/lapack-3.4.2.tgz"
  URL_MD5 61bf1a8a4469d4bdb7604f5897179478
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/lapack)

add_revision(netcdf
  URL http://www.paraview.org/files/dependencies/netcdf-4.3.2.modified.tar.gz
  URL_MD5 1841196c2bfcf10246966eecf92ad0ec)

add_revision(mxml
  URL http://www.paraview.org/files/dependencies/mxml-2.9.tar.gz
  URL_MD5 e21cad0f7aacd18f942aa0568a8dee19
)

add_revision(adios
  URL http://www.paraview.org/files/dependencies/adios-1.8-439f0fb6.tar.bz2
  URL_MD5 a88701c77a7ead5daadd8d8aff70556a
)
