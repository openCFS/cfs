#=============================================================================
#
# This CMake  script sets up  a MinGW  development environment on  Windows XP,
# Vista, 7 and 8 in 32-bit and 64-bit flavors.
#
# Once the control is handed  over from the bootstrap_devel_machine.cmd script
# it is  assumed that 7-Zip and  wget.exe have been downloaded,  installed and
# put  into their  proper  places. It  is also  assumed  that the  environment
# variables for proxy servers have been properly set
#
# The  script checks  first,  if the  machine  is 32-bit  or  64-bit and  only
# installs packages  that can  run on  the machine.  Therefore, one  gets only
# 32-bit compilers and tools on a 32-bit machine.
#
# After that a  number of binary packages gets downloaded  and unpacked in the
# BASEDIR directory.  Once a basic development  environment is in place  a few
# extra  libs (iconv,  zlib, libxml2,  libxslt, TCL  and TK)  get download  as
# sources and compiled  with the 32-bit and 64-bit compilers  to obtain proper
# binaries for further development.
#
# A simple CMake  project gets generated in  example_cmake_project for testing
# the MinGW environment.
#
# In the end the environments  get written to mingw32_env.cmd, msys32_env.cmd,
# mingw64_env.cmd and msys64_env.cmd. The  environments are made available via
# the Console2 shell to which a link  will be generated on the user's desktop.
# Inside Console2 the environments can be accessed via File -> New Tab.
#
#=============================================================================

cmake_policy(SET CMP0007 NEW)

file(TO_CMAKE_PATH "$ENV{BASEDIR}" BASEDIR)
file(TO_CMAKE_PATH "$ENV{PCKGDIR}" PCKGDIR)
file(TO_CMAKE_PATH "$ENV{TMPDIR}" TMPDIR)

set(BOOTSTRAP_PCKG_MIRROR "$ENV{BOOTSTRAP_PCKG_MIRROR}")

message("BASEDIR ${BASEDIR}")
message("PCKGDIR ${PCKGDIR}")
message("TMPDIR ${TMPDIR}")
message("BOOTSTRAP_PCKG_MIRROR ${BOOTSTRAP_PCKG_MIRROR}")

message("CMAKE_COMMAND ${CMAKE_COMMAND}")
get_filename_component(CMAKE_BIN_DIR "${CMAKE_COMMAND}" DIRECTORY)
list(APPEND BIN_DIRS
  "CMake binary directory" # 32-bit
  "${CMAKE_BIN_DIR}"
  "CMake binary directory" # 64-bit
  "${CMAKE_BIN_DIR}"
  )

find_program(SEVENZIP 7z)
message("SEVENZIP ${SEVENZIP}")
get_filename_component(SEVENZIP_BIN_DIR "${SEVENZIP}" DIRECTORY)
list(APPEND BIN_DIRS
  "7-Zip binary directory" # 32-bit
  "${SEVENZIP_BIN_DIR}"
  "7-Zip binary directory" # 64-bit
  "${SEVENZIP_BIN_DIR}"
  )

#=============================================================================
# Check if machine is 32-bit or 64-bit since it makes no sense to install
# 64-bit utilities on a 32-bit machine.
#=============================================================================
set(REGKEY "HKLM\\SYSTEM\\CurrentControlSet\\Control")
set(REGKEY "${REGKEY}\\Session Manager\\Environment")
execute_process(
  COMMAND reg query "${REGKEY}" /v PROCESSOR_ARCHITECTURE
  OUTPUT_VARIABLE PROCESSOR_ARCH
  )
string(REPLACE "\n" ";" PROCESSOR_ARCH ${PROCESSOR_ARCH})
foreach(ITEM IN ITEMS ${PROCESSOR_ARCH})
  if(ITEM MATCHES "PROCESSOR_ARCHITECTURE")
    string(STRIP "${ITEM}" ITEM)
    string(REGEX REPLACE "[\t ]+" ";" ITEM "${ITEM}")
    list(GET ITEM 2 PROCESSOR_ARCH)
    break()
  endif()
endforeach()

set(X64 1)
if(PROCESSOR_ARCH STREQUAL "x86")
  set(X64 0)
endif()
message("PROCESSOR_ARCH ${PROCESSOR_ARCH} ${X64}")

#=============================================================================
# Some programs (patch.exe, install.exe etc.) can conflict with Windows UAC.
# Therefore, we need to write a manifest file for them. Here we provide a
# convenient macro for this task.
#=============================================================================
macro(WRITE_MANIFEST PROGRAM MANIFEST)
  file(WRITE "${MANIFEST}"
"<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>
<assembly xmlns=\"urn:schemas-microsoft-com:asm.v1\" manifestVersion=\"1.0\">
  <assemblyIdentity version=\"7.95.0.0\"
     processorArchitecture=\"X86\"
     name=\"${PROGRAM}\"
     type=\"win32\"/>

  <!-- Identify the application security requirements. -->
  <trustInfo xmlns=\"urn:schemas-microsoft-com:asm.v3\">
    <security>
      <requestedPrivileges>
        <requestedExecutionLevel
          level=\"asInvoker\"
          uiAccess=\"false\"/>
        </requestedPrivileges>
       </security>
  </trustInfo>
</assembly>\n")
endmacro()

#=============================================================================
# A macro for downloading a file from a number of mirrors.
#=============================================================================
macro(DOWNLOAD_FILE_FROM_MIRRORS MIRRORS FILE MD5SUM CHECKMD5)
  set(FILE ${FILE})
  set(MD5SUM ${MD5SUM})
  set(CHECKMD5 ${CHECKMD5})

  SET(PERFORM_DOWNLOAD 0)
  SET(DOWNLOAD_OKAY 0)
  if(${CHECKMD5})
    STRING(STRIP ${MD5SUM} MD5SUM)
    STRING(TOLOWER ${MD5SUM} MD5SUM)
  endif()
  SET(TIMEOUT 1200)

  if(EXISTS "${FILE}")
    if(${CHECKMD5})
      file(MD5 "${FILE}" ACTUAL_MD5SUM)
      STRING(TOLOWER ${ACTUAL_MD5SUM} ACTUAL_MD5SUM)
      
      STRING(COMPARE EQUAL ${MD5SUM} ${ACTUAL_MD5SUM} MD5_EQUAL)
      
      IF(NOT MD5_EQUAL)
	FILE(REMOVE ${FILE})
	SET(PERFORM_DOWNLOAD 1)   
      ELSE()                                                          
	SET(DOWNLOAD_OKAY 1)
      ENDIF()
    else()
      SET(DOWNLOAD_OKAY 1)
    endif()
  else()
    SET(PERFORM_DOWNLOAD 1)
  endif()

  IF(PERFORM_DOWNLOAD)
    foreach(URL IN ITEMS ${MIRRORS})
      message(
        "downloading...\n  src='${URL}'\n  dst='${FILE}'\n  timeout=${TIMEOUT}"
        )

      file(DOWNLOAD "${URL}" "${FILE}"
        INACTIVITY_TIMEOUT ${TIMEOUT}
        TIMEOUT ${TIMEOUT}
        LOG DL_LOG
        SHOW_PROGRESS
        STATUS DLSTAT
        )

      LIST(GET DLSTAT 0 DL_FAIL)
      LIST(GET DLSTAT 1 DL_MSG)

      IF(NOT DL_FAIL)
	if(${CHECKMD5})
          FILE(MD5 ${FILE} ACTUAL_MD5SUM)
          STRING(TOLOWER ${ACTUAL_MD5SUM} ACTUAL_MD5SUM)
	  
          STRING(COMPARE EQUAL ${MD5SUM} ${ACTUAL_MD5SUM} MD5_EQUAL)
	  
          IF(MD5_EQUAL)
            SET(DOWNLOAD_OKAY 1)
            BREAK()
          ELSE()
            FILE(REMOVE ${FILE})
          ENDIF()
	else()
	  break()
	endif()
      ELSE()
        MESSAGE("Download failed: ${DL_MSG}")
        FILE(REMOVE ${FILE})
      ENDIF()
    endforeach()
  ENDIF()
  if(NOT EXISTS "${FILE}")
    message(FATAL_ERROR "\nDownload of ${FILE} failed!\n")
  endif()
endmacro(DOWNLOAD_FILE_FROM_MIRRORS)

#=============================================================================
# A macro for writing out batch files for MinGW environments.
#=============================================================================
macro(write_mingw_env)
  list(LENGTH BIN_DIRS LEN_BIN_DIRS)
  math(EXPR NUM_BIN_DIRS "${LEN_BIN_DIRS}/4")

  set(MINGW_ENV_32 "@echo off

set BASEDIR=${BASEDIR_NATIVE}

set PATH=%SystemRoot%\\system32;%SystemRoot%;%SystemRoot%\\System32\\Wbem

set CMAKE_DEFAULT_GENERATOR=MinGW Makefiles

set CC=gcc
set CXX=g++
set FC=gfortran
")
  set(MINGW_ENV_64 "${MINGW_ENV_32}")

  math(EXPR LEN_BIN_DIRS "${LEN_BIN_DIRS}-1")
  foreach(IDX RANGE 0 ${LEN_BIN_DIRS} 4)
    # message("IDX ${IDX}")
    math(EXPR IDX_DESCR_32BIT "${IDX}")
    math(EXPR IDX_BINDIR_32BIT "${IDX}+1")
    math(EXPR IDX_DESCR_64BIT "${IDX}+2")
    math(EXPR IDX_BINDIR_64BIT "${IDX}+3")

    list(GET BIN_DIRS ${IDX_DESCR_32BIT} DESCR_32BIT)
    list(GET BIN_DIRS ${IDX_BINDIR_32BIT} BINDIR_32BIT)
    list(GET BIN_DIRS ${IDX_DESCR_64BIT} DESCR_64BIT)
    list(GET BIN_DIRS ${IDX_BINDIR_64BIT} BINDIR_64BIT)


    file(TO_NATIVE_PATH ${BINDIR_32BIT} BINDIR_32BIT)
    string(REPLACE "${BASEDIR_NATIVE}" "%BASEDIR%" BINDIR_32BIT "${BINDIR_32BIT}")
    file(TO_NATIVE_PATH ${BINDIR_64BIT} BINDIR_64BIT)
    string(REPLACE "${BASEDIR_NATIVE}" "%BASEDIR%" BINDIR_64BIT "${BINDIR_64BIT}")
    # message("DESCR_32BIT ${DESCR_32BIT}")
    # message("BINDIR_32BIT ${BINDIR_32BIT}")
    # message("DESCR_64BIT ${DESCR_64BIT}")
    # message("BINDIR_64BIT ${BINDIR_64BIT}")
    set(MINGW_ENV_32 "${MINGW_ENV_32}rem ${DESCR_32BIT}\n")
    set(MINGW_ENV_32 "${MINGW_ENV_32}set PATH=${BINDIR_32BIT};%PATH%\n\n")
    set(MINGW_ENV_64 "${MINGW_ENV_64}rem ${DESCR_64BIT}\n")
    set(MINGW_ENV_64 "${MINGW_ENV_64}set PATH=${BINDIR_64BIT};%PATH%\n\n")
  endforeach()

  set(MINGW_ENV_32 "${MINGW_ENV_32}set PYTHONHOME=%BASEDIR%\\mingw32\\opt\n\n")
  set(MINGW_ENV_64 "${MINGW_ENV_64}set PYTHONHOME=%BASEDIR%\\mingw64\\opt\n\n")

  file(WRITE mingw32_env.cmd "${MINGW_ENV_32}")
  file(WRITE msys32_env.cmd "@echo off
set BASEDIR=${BASEDIR_NATIVE}
call %BASEDIR%\\mingw32_env.cmd
%BASEDIR%\\Git\\git-bash.bat\n")
  if(X64)
    file(WRITE mingw64_env.cmd "${MINGW_ENV_64}")
    file(WRITE msys64_env.cmd "@echo off
set BASEDIR=${BASEDIR_NATIVE}
call %BASEDIR%\\mingw64_env.cmd
%BASEDIR%\\Git\\git-bash.bat\n")
  endif()
endmacro(write_mingw_env)

#=============================================================================
# Vim text editor.
#=============================================================================
SET(VIM_DIR "${BASEDIR}/vim")

#-----------------------------------------------------------------------------
# Download Vim
#-----------------------------------------------------------------------------
set(VIM_EXE "gvim74.exe")
if(NOT EXISTS "${PCKGDIR}/${VIM_EXE}")
  set(MIRRORS
    ftp://ftp.vim.org/pub/vim/pc/${VIM_EXE}
    ${BOOTSTRAP_PCKG_MIRROR}/${VIM_EXE}
    )

  DOWNLOAD_FILE_FROM_MIRRORS(
    "${MIRRORS}"
    "${PCKGDIR}/${VIM_EXE}"
    9f6e4aa233c1eb8a18526681308b42d6
    1
    )
endif()

#-----------------------------------------------------------------------------
# Unpack Vim
#-----------------------------------------------------------------------------
if(NOT EXISTS ${VIM_DIR})
  execute_process(COMMAND ${SEVENZIP} x -y ${PCKGDIR}/${VIM_EXE}
    WORKING_DIRECTORY ${TMPDIR}
    )

  execute_process(
    COMMAND ${CMAKE_COMMAND} -E copy_directory "\$0" "${VIM_DIR}"
    WORKING_DIRECTORY ${TMPDIR}
    )
  write_manifest(install.exe "${VIM_DIR}/install.exe.manifest")

  file(REMOVE_RECURSE ${TMPDIR})
  file(MAKE_DIRECTORY ${TMPDIR})
endif()
list(APPEND BIN_DIRS
  "Binary directory for Vim"
  "${VIM_DIR}"
  "Binary directory for Vim"
  "${VIM_DIR}"
  )

#=============================================================================
# Doxygen.
#=============================================================================
SET(DOXYGEN_DIR "${BASEDIR}/doxygen")

#-----------------------------------------------------------------------------
# Download Doxygen
#-----------------------------------------------------------------------------
set(DOXYGEN_ZIP "doxygen-1.8.6.windows.bin.zip")
if(NOT EXISTS "${PCKGDIR}/${DOXYGEN_ZIP}")
  set(MIRRORS
    ftp://ftp.stack.nl/pub/users/dimitri/${DOXYGEN_ZIP}
    ${BOOTSTRAP_PCKG_MIRROR}/${DOXYGEN_ZIP}
    )

  DOWNLOAD_FILE_FROM_MIRRORS(
    "${MIRRORS}"
    "${PCKGDIR}/${DOXYGEN_ZIP}"
    20a290d1bc17cacb9be9961c04f267b9
    1
    )
endif()

#-----------------------------------------------------------------------------
# Unpack Doxygen
#-----------------------------------------------------------------------------
if(NOT EXISTS ${DOXYGEN_DIR})
  file(MAKE_DIRECTORY "${DOXYGEN_DIR}")
  execute_process(COMMAND ${SEVENZIP} x -y ${PCKGDIR}/${DOXYGEN_ZIP}
    WORKING_DIRECTORY ${DOXYGEN_DIR}
    )
endif()
list(APPEND BIN_DIRS
  "Binary directory for Doxygen"
  "${DOXYGEN_DIR}"
  "Binary directory for Doxygen"
  "${DOXYGEN_DIR}"
  )

#=============================================================================
# Graphviz.
#=============================================================================
SET(GRAPHVIZ_DIR "${BASEDIR}/graphviz")

#-----------------------------------------------------------------------------
# Download Graphviz
#-----------------------------------------------------------------------------
set(GRAPHVIZ_ZIP "graphviz-2.36.zip")
if(NOT EXISTS "${PCKGDIR}/${GRAPHVIZ_ZIP}")
  set(MIRRORS
    http://www.graphviz.org/pub/graphviz/stable/windows/${GRAPHVIZ_ZIP}
    ${BOOTSTRAP_PCKG_MIRROR}/${GRAPHVIZ_ZIP}
    )

  DOWNLOAD_FILE_FROM_MIRRORS(
    "${MIRRORS}"
    "${PCKGDIR}/${GRAPHVIZ_ZIP}"
    105ee5de1883767a3947ce96993ffab7  
    1
    )
endif()

#-----------------------------------------------------------------------------
# Unpack Graphviz
#-----------------------------------------------------------------------------
if(NOT EXISTS ${GRAPHVIZ_DIR})
  execute_process(COMMAND ${SEVENZIP} x -y "${PCKGDIR}/${GRAPHVIZ_ZIP}"
    WORKING_DIRECTORY ${TMPDIR}
    )
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E copy_directory "release" "${GRAPHVIZ_DIR}"
    WORKING_DIRECTORY ${TMPDIR}
    )

  file(REMOVE_RECURSE ${TMPDIR})
  file(MAKE_DIRECTORY ${TMPDIR})
endif()
list(APPEND BIN_DIRS
  "Binary directory for Graphviz"
  "${GRAPHVIZ_DIR}/bin"
  "Binary directory for Graphviz"
  "${GRAPHVIZ_DIR}/bin"
  )

#=============================================================================
# MikTex.
#=============================================================================
SET(MIKTEX_DIR "${BASEDIR}/miktex-portable")

#-----------------------------------------------------------------------------
# Download MikTex
#-----------------------------------------------------------------------------
set(MIKTEX_EXE "miktex-portable-2.9.5105.exe")
if(NOT EXISTS "${PCKGDIR}/${MIKTEX_EXE}")
  set(FTP_SERVER ftp://ftp.mpi-sb.mpg.de/pub/tex/mirror/ftp.dante.de)
  set(FTP_SERVER ${FTP_SERVER}/pub/tex/systems/win32/miktex/setup)

  set(MIRRORS
    ${FTP_SERVER}/${MIKTEX_EXE}
    ${BOOTSTRAP_PCKG_MIRROR}/${MIKTEX_EXE}
    )

  DOWNLOAD_FILE_FROM_MIRRORS(
    "${MIRRORS}"
    "${PCKGDIR}/${MIKTEX_EXE}"
    8da199906dcdc0f31814c9ed1fea67e4  
    1
    )
endif()

#-----------------------------------------------------------------------------
# Unpack Miktex and install some required extra packages.
#-----------------------------------------------------------------------------
if(NOT EXISTS ${MIKTEX_DIR})
  file(MAKE_DIRECTORY "${MIKTEX_DIR}")
  execute_process(COMMAND ${SEVENZIP} x -y "${PCKGDIR}/${MIKTEX_EXE}"
    WORKING_DIRECTORY "${MIKTEX_DIR}"
    )
  file(WRITE "${TMPDIR}/miktex-extra-packages.txt"
"algorithms\nwrapfig\nkoma-script\ncolortbl\nmptopdf\nsubfig\ncaption
booktabs\nrotating\npgf\nms\nxcolor\nurl\noverpic\neepic\nattachfile
todonotes\nminted\nfancyvrb\nifplatform\nmiktex-tex4ht-bin-2.9\n
miktex-tex4ht-base-2.9\numlaute\nsubfigure")
  execute_process(COMMAND "${MIKTEX_DIR}/miktex/bin/mpm"
    --install-some=miktex-extra-packages.txt
    WORKING_DIRECTORY "${TMPDIR}"
    )
  file(REMOVE_RECURSE ${TMPDIR})
  file(MAKE_DIRECTORY ${TMPDIR})
endif()
list(APPEND BIN_DIRS
  "Binary directory for MikTex"
  "${MIKTEX_DIR}/miktex/bin"
  "Binary directory for MikTex"
  "${MIKTEX_DIR}/miktex/bin"
  )

#=============================================================================
# ImageMagick - tools for image handling.
#=============================================================================
SET(IMGMGCK_DIR "${BASEDIR}/ImageMagick-6.8.8-7")

#-----------------------------------------------------------------------------
# Download ImageMagick
#-----------------------------------------------------------------------------
set(IMGMGCK_ZIP "ImageMagick-6.8.8-7-Q16-x86-windows.zip")

if(NOT EXISTS "${PCKGDIR}/${IMGMGCK_ZIP}")
  set(MIRRORS
    http://www.imagemagick.org/download/binaries/${IMGMGCK_ZIP}
    ${BOOTSTRAP_PCKG_MIRROR}/${IMGMGCK_ZIP}
    )

  DOWNLOAD_FILE_FROM_MIRRORS(
    "${MIRRORS}"
    "${PCKGDIR}/${IMGMGCK_ZIP}"
    d13d9c654b4555a907001c4e7503514f
    1
    )
endif()

#-----------------------------------------------------------------------------
# Unpack ImageMagick
#-----------------------------------------------------------------------------
if(NOT EXISTS ${IMGMGCK_DIR})
  execute_process(COMMAND
    ${CMAKE_COMMAND} -E tar xvf ${PCKGDIR}/${IMGMGCK_ZIP}
    WORKING_DIRECTORY ${BASEDIR}
    )
endif()
list(APPEND BIN_DIRS
  "Binary directory for ImageMagick"
  "${IMGMGCK_DIR}"
  "Binary directory for ImageMagick"
  "${IMGMGCK_DIR}"
  )

#=============================================================================
# Unix utilities (UnxUtils) for Windows. Maybe these should be replaced by
# GnuWin utilities from http://gnuwin32.sourceforge.net some time.
#=============================================================================

#-----------------------------------------------------------------------------
# Set directory for UnxUtils
#-----------------------------------------------------------------------------
set(UNXUTILS_DIR "${BASEDIR}/UnxUtils")

#-----------------------------------------------------------------------------
# Download UnxUtils
#-----------------------------------------------------------------------------
set(UNXUTILS_ZIP "UnxUtils.zip")
if(NOT EXISTS "${PCKGDIR}/${UNXUTILS_ZIP}")
  # Download UnxUtils
  set(NETCOLOGNE_MIRROR http://netcologne.dl.sourceforge.net/project)
  set(NETCOLOGNE_MIRROR ${NETCOLOGNE_MIRROR}/unxutils/unxutils/current)
  set(MIRRORS
    ${NETCOLOGNE_MIRROR}/${UNXUTILS_ZIP}
    ${BOOTSTRAP_PCKG_MIRROR}/${UNXUTILS_ZIP}
    ftp://ftp.fh-hannover.de/pandora/files/linux/${UNXUTILS_ZIP}
    ftp://ftp.ill.fr/pub/cs/reg/${UNXUTILS_ZIP}
    ftp://ftp.ufanet.ru/pub/windows/unixutils/${UNXUTILS_ZIP}
    ftp://ftp.ustb.edu.cn/pub/unix_utils/${UNXUTILS_ZIP}
    http://unxutils.sourceforge.net/${UNXUTILS_ZIP}
    )

  DOWNLOAD_FILE_FROM_MIRRORS(
    "${MIRRORS}"
    "${PCKGDIR}/${UNXUTILS_ZIP}"
    59567db7fc650e1e778f91702ed109c6
    1
    )
endif()

set(UNXUPDATES_ZIP "UnxUpdates.zip")
if(NOT EXISTS "${PCKGDIR}/${UNXUPDATES_ZIP}")
  # Download UnxUpdates
  set(MIRRORS
    ftp://ftp.fh-hannover.de/pandora/files/linux/${UNXUPDATES_ZIP}
    http://unxutils.sourceforge.net/${UNXUPDATES_ZIP}
    ${BOOTSTRAP_PCKG_MIRROR}/${UNXUPDATES_ZIP}
    )

  DOWNLOAD_FILE_FROM_MIRRORS(
    "${MIRRORS}"
    "${PCKGDIR}/${UNXUPDATES_ZIP}"
    9b44ede7449f991439b32bbdd844df7b
    1
    )
endif()

set(BASE64_ZIP "base64.zip")
if(NOT EXISTS "${PCKGDIR}/${BASE64_ZIP}")
  # Download base64
  set(MIRRORS
    http://www.fourmilab.ch/webtools/base64/${BASE64_ZIP}
    ${BOOTSTRAP_PCKG_MIRROR}/${BASE64_ZIP}
    )

  DOWNLOAD_FILE_FROM_MIRRORS(
    "${MIRRORS}"
    "${PCKGDIR}/${BASE64_ZIP}"
    f4cf691250e66f6dad822e9c2ace13c6
    1
    )
endif()

set(FILE_ZIP "file-5.03-bin.zip")
if(NOT EXISTS "${PCKGDIR}/${FILE_ZIP}")
  set(MIRRORS
    http://downloads.sourceforge.net/gnuwin32/${FILE_ZIP}
    ${BOOTSTRAP_PCKG_MIRROR}/${FILE_ZIP}
    )

  DOWNLOAD_FILE_FROM_MIRRORS(
    "${MIRRORS}"
    "${PCKGDIR}/${FILE_ZIP}"
    1c003c59ee9fb04ca84a5488543a69e0
    1
    )
endif()

set(REGEX_ZIP "regex-2.7-bin.zip")
if(NOT EXISTS "${PCKGDIR}/${REGEX_ZIP}")
  set(MIRRORS
    http://kent.dl.sourceforge.net/project/gnuwin32/regex/2.7/${REGEX_ZIP}
    ${BOOTSTRAP_PCKG_MIRROR}/${REGEX_ZIP}
    )

  DOWNLOAD_FILE_FROM_MIRRORS(
    "${MIRRORS}"
    "${PCKGDIR}/${REGEX_ZIP}"
    1c097334c32f7b977093f38b105e7580
    1
    )
endif()

set(ACK_PL "ack.pl")
if(NOT EXISTS "${PCKGDIR}/${ACK_PL}")
  set(MIRRORS
    http://beyondgrep.com/ack-2.12-single-file
    ${BOOTSTRAP_PCKG_MIRROR}/${ACK_PL}
    )

  DOWNLOAD_FILE_FROM_MIRRORS(
    "${MIRRORS}"
    "${PCKGDIR}/${ACK_PL}"
    NO_MD5_CHECK
    0
    )
endif()

set(CURL_ZIP "curl-7.35.0-win32-fix1.zip")
if(NOT EXISTS "${PCKGDIR}/${CURL_ZIP}")
  set(MIRRORS
    http://www.confusedbycode.com/curl/${CURL_ZIP}
    ${BOOTSTRAP_PCKG_MIRROR}/${CURL_ZIP}
    )

  DOWNLOAD_FILE_FROM_MIRRORS(
    "${MIRRORS}"
    "${PCKGDIR}/${CURL_ZIP}"
    5f9f91ecd51ec975a65e1dc25b0b1861
    1
    )
endif()

#-----------------------------------------------------------------------------
# Download latest version of wget.exe for Windows (needed e.g. for JRE)
#-----------------------------------------------------------------------------
set(WGET_EXE "wget.exe")
set(WGET_DOWNLOAD 0)
if(EXISTS "${PCKGDIR}/${WGET_EXE}")
  execute_process(
    COMMAND "${PCKGDIR}/${WGET_EXE}" --version
    OUTPUT_VARIABLE WGET_OUTPUT
    )
  string(REPLACE "\n" ";" WGET_OUTPUT "${WGET_OUTPUT}")
  list(GET WGET_OUTPUT 0 WGET_VERSION)
  string(REPLACE "GNU Wget " "" WGET_VERSION "${WGET_VERSION}")
  if(WGET_VERSION VERSION_LESS "1.11")
    file(REMOVE "${PCKGDIR}/${WGET_EXE}")
    set(WGET_DOWNLOAD 1)
  endif()
else()
  set(WGET_DOWNLOAD 1)
endif()

if(WGET_DOWNLOAD)
  # Download wget.exe
  set(MIRRORS
    http://users.ugent.be/~bpuype/cgi-bin/fetch.pl?dl=wget/${WGET_EXE}
    ${BOOTSTRAP_PCKG_MIRROR}/${WGET_EXE}
    )

  DOWNLOAD_FILE_FROM_MIRRORS(
    "${MIRRORS}"
    "${PCKGDIR}/${WGET_EXE}"
    MD5SUM_NOT_NEEDED
    0
    )
endif()

#-----------------------------------------------------------------------------
# Unpack UnxUtils
#-----------------------------------------------------------------------------
if(NOT EXISTS ${UNXUTILS_DIR})
  file(MAKE_DIRECTORY ${UNXUTILS_DIR})
  execute_process(
    COMMAND ${SEVENZIP} x "${PCKGDIR}/${UNXUTILS_ZIP}"
    WORKING_DIRECTORY ${UNXUTILS_DIR}
    )
  execute_process(
    COMMAND ${SEVENZIP} x -y "${PCKGDIR}/${UNXUPDATES_ZIP}"
    WORKING_DIRECTORY ${UNXUTILS_DIR}/usr/local/wbin
    )
  execute_process(
    COMMAND ${SEVENZIP} x -y "${PCKGDIR}/${BASE64_ZIP}"
    WORKING_DIRECTORY ${UNXUTILS_DIR}/usr/local/wbin
    )
  execute_process(
    COMMAND ${SEVENZIP} x -y "${PCKGDIR}/${FILE_ZIP}"
    WORKING_DIRECTORY ${UNXUTILS_DIR}/usr/local
    )
  execute_process(
    COMMAND ${SEVENZIP} x -y "${PCKGDIR}/${REGEX_ZIP}"
    WORKING_DIRECTORY ${UNXUTILS_DIR}/usr/local
    )
  execute_process(
    COMMAND ${SEVENZIP} x -y "${PCKGDIR}/${CURL_ZIP}"
    WORKING_DIRECTORY ${TMPDIR}
    )
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    "${TMPDIR}/curl-7.35.0-win32/bin" "bin"
    WORKING_DIRECTORY ${UNXUTILS_DIR}/usr/local
    )
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    "${TMPDIR}/curl-7.35.0-win32/dlls" "bin"
    WORKING_DIRECTORY ${UNXUTILS_DIR}/usr/local
    )

  execute_process(
    COMMAND ${CMAKE_COMMAND} -E copy
    "${PCKGDIR}/${ACK_PL}" "${ACK_PL}"
    WORKING_DIRECTORY ${UNXUTILS_DIR}/usr/local/bin
    )

  file(REMOVE_RECURSE ${TMPDIR})
  file(MAKE_DIRECTORY ${TMPDIR})
endif()
if(WGET_DOWNLOAD)
  file(REMOVE "${UNXUTILS_DIR}/usr/local/wbin/wget.exe")
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E copy "${PCKGDIR}/${WGET_EXE}" wget.exe
    WORKING_DIRECTORY ${UNXUTILS_DIR}/usr/local/wbin
    )
endif()
set(UNXUTILS_BIN_DIR "${UNXUTILS_DIR}/usr/local/wbin")
# Write manifests for patch.exe and install.exe
write_manifest(patch.exe "${UNXUTILS_BIN_DIR}/patch.exe.manifest")
write_manifest(install.exe "${UNXUTILS_BIN_DIR}/install.exe.manifest")
list(APPEND BIN_DIRS
  "Binary directory for UnxUtils" # 32-bit
  "${UNXUTILS_BIN_DIR}"
  "Binary directory for UnxUtils" # 64-bit
  "${UNXUTILS_BIN_DIR}"
  )
list(APPEND BIN_DIRS
  "Binary directory for file.exe" # 32-bit
  "${UNXUTILS_DIR}/usr/local/bin"
  "Binary directory for file.exe" # 64-bit
  "${UNXUTILS_DIR}/usr/local/bin"
  )

#=============================================================================
# SourceGear Diffmerge tool.
#=============================================================================

#-----------------------------------------------------------------------------
# Set directory for Diffmerge
#-----------------------------------------------------------------------------
set(DIFFMERGE_DIR "${BASEDIR}/diffmerge")

#-----------------------------------------------------------------------------
# Download Diffmerge
#-----------------------------------------------------------------------------
set(DIFFMERGE_ZIP "DiffMerge_3_3_0_18513.zip")
if(NOT EXISTS "${PCKGDIR}/${DIFFMERGE_ZIP}")
  # Download Diffmerge
  set(MIRRORS
    http://download-us.sourcegear.com/DiffMerge/3.3.0/Windows/${DIFFMERGE_ZIP}
    ${BOOTSTRAP_PCKG_MIRROR}/${DIFFMERGE_ZIP}
    )

  DOWNLOAD_FILE_FROM_MIRRORS(
    "${MIRRORS}"
    "${PCKGDIR}/${DIFFMERGE_ZIP}"
    4ff5b219d3e07d70a41b8d4fdd53313e
    1
    )
endif()

#-----------------------------------------------------------------------------
# Unpack Diffmerge
#-----------------------------------------------------------------------------
if(NOT EXISTS ${DIFFMERGE_DIR})
  file(MAKE_DIRECTORY ${DIFFMERGE_DIR})
  execute_process(COMMAND ${SEVENZIP} x "${PCKGDIR}/${DIFFMERGE_ZIP}"
    WORKING_DIRECTORY ${DIFFMERGE_DIR}
    )
endif()
list(APPEND BIN_DIRS
  "Binary directory for SourceGear Diffmerge" # 32-bit
  "${DIFFMERGE_DIR}"
  "Binary directory for SourceGear Diffmerge" # 64-bit
  "${DIFFMERGE_DIR}"
  )

#=============================================================================
# ResEdit tool.
#=============================================================================

#-----------------------------------------------------------------------------
# Set directory for ResEdit
#-----------------------------------------------------------------------------
set(RESEDIT_DIR "${BASEDIR}/resedit")

#-----------------------------------------------------------------------------
# Download ResEdit
#-----------------------------------------------------------------------------
set(RESEDIT32_7Z "ResEdit-win32.7z")
if(NOT EXISTS "${PCKGDIR}/${RESEDIT32_7Z}")
  # Download ResEdit
  set(MIRRORS
    http://www.resedit.net/${RESEDIT32_7Z}
    ${BOOTSTRAP_PCKG_MIRROR}/${RESEDIT32_7Z}
    )

  DOWNLOAD_FILE_FROM_MIRRORS(
    "${MIRRORS}"
    "${PCKGDIR}/${RESEDIT32_7Z}"
    3620c72a1b4cdedf0b2df028aa3a0c0f
    0
    )
endif()
set(RESEDIT64_7Z "ResEdit-x64.7z")
if(NOT EXISTS "${PCKGDIR}/${RESEDIT64_7Z}")
  # Download ResEdit
  set(MIRRORS
    http://www.resedit.net/${RESEDIT64_7Z}
    ${BOOTSTRAP_PCKG_MIRROR}/${RESEDIT64_7Z}
    )

  DOWNLOAD_FILE_FROM_MIRRORS(
    "${MIRRORS}"
    "${PCKGDIR}/${RESEDIT64_7Z}"
    fea492fa84f5c0da8cd7ca69a6b4dc46
    0
    )
endif()

#-----------------------------------------------------------------------------
# Unpack ResEdit
#-----------------------------------------------------------------------------
if(NOT EXISTS ${RESEDIT_DIR})
  file(MAKE_DIRECTORY ${RESEDIT_DIR})
  file(MAKE_DIRECTORY ${RESEDIT_DIR}/x86)
  execute_process(COMMAND ${SEVENZIP} x "${PCKGDIR}/${RESEDIT32_7Z}"
    WORKING_DIRECTORY ${RESEDIT_DIR}/x86
    )
  if(X64)
  file(MAKE_DIRECTORY ${RESEDIT_DIR}/x64)
    execute_process(COMMAND ${SEVENZIP} x "${PCKGDIR}/${RESEDIT64_7Z}"
      WORKING_DIRECTORY ${RESEDIT_DIR}/x64
      )
  endif()
endif()
list(APPEND BIN_DIRS
  "Binary directory for 32-bit ResEdit" # 32-bit
  "${RESEDIT_DIR}/x86"
  "Binary directory for 64-bit ResEdit" # 64-bit
  "${RESEDIT_DIR}/x64"
  )

#=============================================================================
# Dependency walker tool.
#=============================================================================

#-----------------------------------------------------------------------------
# Set directory for Depends
#-----------------------------------------------------------------------------
set(DEPENDS_DIR "${BASEDIR}/depends")

#-----------------------------------------------------------------------------
# Download Depends
#-----------------------------------------------------------------------------
set(DEPENDS32_ZIP "depends22_x86.zip")
if(NOT EXISTS "${PCKGDIR}/${DEPENDS32_ZIP}")
  # Download Depends
  set(MIRRORS
    http://www.dependencywalker.com/${DEPENDS32_ZIP}
    ${BOOTSTRAP_PCKG_MIRROR}/${DEPENDS32_ZIP}
    )

  DOWNLOAD_FILE_FROM_MIRRORS(
    "${MIRRORS}"
    "${PCKGDIR}/${DEPENDS32_ZIP}"
    675ca981ddf557eb7d4550624157dbe5
    1
    )
endif()
set(DEPENDS64_ZIP "depends22_x64.zip")
if(NOT EXISTS "${PCKGDIR}/${DEPENDS64_ZIP}")
  # Download Depends
  set(MIRRORS
    http://www.dependencywalker.com/${DEPENDS64_ZIP}
    ${BOOTSTRAP_PCKG_MIRROR}/${DEPENDS64_ZIP}
    )

  DOWNLOAD_FILE_FROM_MIRRORS(
    "${MIRRORS}"
    "${PCKGDIR}/${DEPENDS64_ZIP}"
    7975054e322794cd332d5f1c00eeec5f
    1
    )
endif()

#-----------------------------------------------------------------------------
# Unpack Depends
#-----------------------------------------------------------------------------
if(NOT EXISTS ${DEPENDS_DIR})
  file(MAKE_DIRECTORY ${DEPENDS_DIR})
  file(MAKE_DIRECTORY ${DEPENDS_DIR}/x86)
  execute_process(COMMAND ${SEVENZIP} x -y "${PCKGDIR}/${DEPENDS32_ZIP}"
    WORKING_DIRECTORY ${DEPENDS_DIR}/x86
    )
  if(X64)
    file(MAKE_DIRECTORY ${DEPENDS_DIR}/x64)
    execute_process(COMMAND ${SEVENZIP} x -y "${PCKGDIR}/${DEPENDS64_ZIP}"
      WORKING_DIRECTORY ${DEPENDS_DIR}/x64
      )
  endif()
endif()
list(APPEND BIN_DIRS
  "Binary directory for 32-bit Depends" # 32-bit
  "${DEPENDS_DIR}/x86"
  "Binary directory for 64-bit Depends" # 64-bit
  "${DEPENDS_DIR}/x64"
  )

#=============================================================================
# SlikSVN Subversion client.
#=============================================================================
SET(SLIKSVN_DIR "${BASEDIR}/SlikSvn")

#-----------------------------------------------------------------------------
# Download Sliksvn
#-----------------------------------------------------------------------------
set(SLIKSVN_MSI "Slik-Subversion-1.8.5-1-win32.msi")
if(NOT EXISTS "${PCKGDIR}/${SLIKSVN_MSI}")
  # Download Sliksvn
  set(MIRRORS
    http://www.sliksvn.com/pub/${SLIKSVN_MSI}
    ${BOOTSTRAP_PCKG_MIRROR}/${SLIKSVN_MSI}
    )

  DOWNLOAD_FILE_FROM_MIRRORS(
    "${MIRRORS}"
    "${PCKGDIR}/${SLIKSVN_MSI}"
    09178b1c6e564a619c0ce30a9dd08c0e
    1
    )
endif()

#-----------------------------------------------------------------------------
# Unpack Sliksvn
#-----------------------------------------------------------------------------
if(NOT EXISTS ${SLIKSVN_DIR})
  #file(TO_NATIVE_PATH ${TMPDIR} TARGET_DIR)
  #execute_process(
  #  COMMAND msiexec /a ${PCKGDIR}/${SLIKSVN_MSI} TARGETDIR=${TARGET_DIR} /qn
  #  WORKING_DIRECTORY ${TMPDIR}
  #  )

  execute_process(
    COMMAND ${SEVENZIP} x -y "${PCKGDIR}/${SLIKSVN_MSI}"
    WORKING_DIRECTORY ${TMPDIR}
    )
  set(SLIKDLLS
    libsvn_client_1.dll   sliksvn-libsvn_client-1.dll
    libsvn_subr_1.dll     sliksvn-libsvn_subr-1.dll
    libapr_1.dll          sliksvn-libapr-1.dll
    libsvn_delta_1.dll    sliksvn-libsvn_delta-1.dll
    libsvn_wc_1.dll       sliksvn-libsvn_wc-1.dll
    libaprutil_1.dll      sliksvn-libaprutil-1.dll
    libsvn_diff_1.dll     sliksvn-libsvn_diff-1.dll
    libsvnjavahl_1.dll    libsvnjavahl-1.dll
    libsvn_fs_1.dll       sliksvn-libsvn_fs-1.dll
    libintl.dll           sliksvn-libintl.dll
    libsvn_ra_1.dll       sliksvn-libsvn_ra-1.dll
    libsasl21.dll         sliksvn-libsasl21.dll
    libsvn_repos_1.dll    sliksvn-libsvn_repos-1.dll
    ssleay32.dll          sliksvn-ssleay32.dll
    libeay32.dll          sliksvn-libeay32.dll
    DB44_20.dll           sliksvn-DB44-20-win32.dll
    F_CENTRAL_msvcp100_x86.1DEE2A86_2F57_3629_8107_A71DBB4DBED2 msvcp100.dll
    F_CENTRAL_msvcr100_x86.1DEE2A86_2F57_3629_8107_A71DBB4DBED2 msvcr100.dll
    )

  foreach(IDX RANGE 0 35 2)
    math(EXPR SRCIDX "${IDX}")
    math(EXPR DSTIDX "${IDX}+1")
    list(GET SLIKDLLS ${SRCIDX} SRC)
    list(GET SLIKDLLS ${DSTIDX} DST)
    
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy ${SRC} ${DST}
      WORKING_DIRECTORY ${TMPDIR}
      )
    file(REMOVE "${TMPDIR}/${SRC}")
  endforeach()

  execute_process(
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${TMPDIR} ${SLIKSVN_DIR}
    WORKING_DIRECTORY ${TMPDIR}
    )
  file(REMOVE_RECURSE ${TMPDIR})
  file(MAKE_DIRECTORY ${TMPDIR})
endif()
list(APPEND BIN_DIRS
  "Binary directory for Sliksvn" # 32-bit
  "${SLIKSVN_DIR}"
  "Binary directory for Sliksvn" # 64-bit
  "${SLIKSVN_DIR}"
  )

#=============================================================================
# Git client which contains a basic MSYS environment. We extend this
# environment with a MSYS compatible GNU make.exe in order to be able to
# compile automake/autoconf-based software (configure/make/make install).
#=============================================================================
SET(GIT_DIR "${BASEDIR}/Git")

#-----------------------------------------------------------------------------
# Download Git
#-----------------------------------------------------------------------------
set(GIT_7Z "PortableGit-1.9.0-preview20140217.7z")
if(NOT EXISTS "${PCKGDIR}/${GIT_7Z}")
  set(MIRRORS
    http://msysgit.googlecode.com/files/${GIT_7Z}
    ${BOOTSTRAP_PCKG_MIRROR}/${GIT_7Z}
    )

  DOWNLOAD_FILE_FROM_MIRRORS(
    "${MIRRORS}"
    "${PCKGDIR}/${GIT_7Z}"
    64880dac9f0196365a7ae45becf8af23
    1
    )
endif()
set(MAKE_LZMA "make-3.81-3-msys-1.0.13-bin.tar.lzma")
if(NOT EXISTS "${PCKGDIR}/${MAKE_LZMA}")
  set(SERVER 
  http://skylink.dl.sourceforge.net/project/mingw/MSYS/Base/make/make-3.81-3)

  set(MIRRORS
    ${SERVER}/${MAKE_LZMA}
    ${BOOTSTRAP_PCKG_MIRROR}/${MAKE_LZMA}
    )

  DOWNLOAD_FILE_FROM_MIRRORS(
    "${MIRRORS}"
    "${PCKGDIR}/${MAKE_LZMA}"
    4fe6fe4dd974c73803413a44bbd9bb1d
    1
    )
endif()
set(LIBINTL_LZMA "libintl-0.18.1.1-1-msys-1.0.17-dll-8.tar.lzma")
if(NOT EXISTS "${PCKGDIR}/${LIBINTL_LZMA}")
  set(SERVER http://netcologne.dl.sourceforge.net/project/mingw/MSYS)
  set(SERVER ${SERVER}/Base/gettext/gettext-0.18.1.1-1)

  set(MIRRORS
    ${SERVER}/${LIBINTL_LZMA}
    ${BOOTSTRAP_PCKG_MIRROR}/${LIBINTL_LZMA}
    )

  DOWNLOAD_FILE_FROM_MIRRORS(
    "${MIRRORS}"
    "${PCKGDIR}/${LIBINTL_LZMA}"
    8d3561e8f94acfe0b9139aa658bb5e34
    1
    )
endif()
set(LIBICONV_LZMA "libiconv-1.14-1-msys-1.0.17-dll-2.tar.lzma")
if(NOT EXISTS "${PCKGDIR}/${LIBICONV_LZMA}")
  set(SERVER http://skylink.dl.sourceforge.net/project/mingw/MSYS)
  set(SERVER ${SERVER}/Base/libiconv/libiconv-1.14-1)

  set(MIRRORS
    ${SERVER}/${LIBICONV_LZMA}
    ${BOOTSTRAP_PCKG_MIRROR}/${LIBICONV_LZMA}
    )

  DOWNLOAD_FILE_FROM_MIRRORS(
    "${MIRRORS}"
    "${PCKGDIR}/${LIBICONV_LZMA}"
    6ffa2b9fdb3711582d6a99bab4389de6
    1
    )
endif()

#-----------------------------------------------------------------------------
# Unpack Git
#-----------------------------------------------------------------------------
if(NOT EXISTS ${GIT_DIR})
  file(MAKE_DIRECTORY ${GIT_DIR})
  execute_process(COMMAND ${SEVENZIP} x "${PCKGDIR}/${GIT_7Z}"
    WORKING_DIRECTORY ${GIT_DIR}
    )

  set(MSYS_LZMAS ${MAKE_LZMA} ${LIBINTL_LZMA} ${LIBICONV_LZMA})
  foreach(LZMA IN ITEMS ${MSYS_LZMAS})
    execute_process(
      COMMAND ${SEVENZIP} x -y "${PCKGDIR}/${LZMA}"
      WORKING_DIRECTORY ${TMPDIR}
      )
    string(REPLACE ".lzma" "" LZMA ${LZMA})
    execute_process(
      COMMAND ${SEVENZIP} x -y "${TMPDIR}/${LZMA}"
      WORKING_DIRECTORY ${GIT_DIR}
      )
  endforeach()

  file(MAKE_DIRECTORY "${GIT_DIR}/etc/profile.d") 
  file(WRITE "${GIT_DIR}/etc/profile.d/init_env.sh"
"# normalize paths to unix paths
export BASEDIR=\"\$(cd \"\$BASEDIR\" ; pwd)\"
export PYTHONHOME=\"\$(cd \"\$PYTHONHOME\" ; pwd)\"

export PATH=\$BASEDIR/Git/bin:\$PYTHONHOME/Scripts:\$PATH

export CMAKE_DEFAULT_GENERATOR='MSYS Makefiles'\n")
  file(REMOVE_RECURSE ${TMPDIR})
  file(MAKE_DIRECTORY ${TMPDIR})
endif()

file(TO_NATIVE_PATH "${GIT_DIR}" GIT_DIR_NATIVE)
file(WRITE "${UNXUTILS_DIR}/usr/local/bin/ack.cmd"
  "@echo off
${GIT_DIR_NATIVE}\\bin\\perl.exe ${UNXUTILS_DIR}/usr/local/bin/${ACK_PL} %*    
")

list(APPEND BIN_DIRS
  "Binary directory for Git" # 32-bit
  "${GIT_DIR}/cmd"
  "Binary directory for Git" # 64-bit
  "${GIT_DIR}/cmd"
  )

#=============================================================================
# GNU compiler 4.8.x.
#=============================================================================
#-----------------------------------------------------------------------------
# Download GCC
#
# For a thorough discussion on the differences between structured exception
# handling (SEH 64-bit), setjmp - longjmp exception handling (SJLJ 32-bit and
# 64-bit) and dwarf-2 exception handling (DWARF 32-bit only)
# cf. http://qt-project.org/wiki/MinGW-64-bit
#-----------------------------------------------------------------------------
set(GCC64_7Z "x64-4.8.0-release-win32-sjlj-rev2.7z")
# We rather take the MinGW builds 32-bit compiler 4.8.1 
# release-win32-sjlj-rev1 that comes with Qt 4.8.5 to make sure everything
# fits perfectly together.
if(NOT EXISTS "${PCKGDIR}/${GCC64_7Z}")
  # Download GCC 4.8
  set(SERVER http://netcologne.dl.sourceforge.net/project/mingwbuilds)
  set(SERVER ${SERVER}/host-windows/releases/4.8.0/64-bit/threads-win32/sjlj)

  set(MIRRORS
    ${SERVER}/${GCC64_7Z}
    ${BOOTSTRAP_PCKG_MIRROR}/${GCC64_7Z}
    )

  DOWNLOAD_FILE_FROM_MIRRORS(
    "${MIRRORS}"
    "${PCKGDIR}/${GCC64_7Z}"
    5302609bd82925a8c8a41338e890a617
    1
    )
endif()

#-----------------------------------------------------------------------------
# Unpack GCC 4.8
#-----------------------------------------------------------------------------
if(X64)
  if(NOT EXISTS mingw64)
    execute_process(
      COMMAND ${SEVENZIP} x -y "${PCKGDIR}/${GCC64_7Z}"
      WORKING_DIRECTORY ${BASEDIR}
      )
  endif()
endif()

set(PYTHONHOME_32 "${BASEDIR}/mingw32/opt")
set(PYTHONHOME_64 "${BASEDIR}/mingw64/opt")

# 32-bit compiler comes with Qt4 and will be installed from this package.

# The binary directories for GCC need to be added at the end of this script in
# order to have them in front of the PATH variable.

#=============================================================================
# Qt 4.8.4.
#=============================================================================

#-----------------------------------------------------------------------------
# Download Qt4
#-----------------------------------------------------------------------------
set(QT4_MINGW64_7Z  "qt-4.8.4-x64-mingw480r2.7z")
set(QT4_MINGW32_EXE "qt-4.8.5-x86-mingw481r1.exe")
set(INNOUNP_RAR "innounp040.rar")
if(NOT EXISTS "${PCKGDIR}/${QT4_MINGW64_7Z}")
  set(SERVER http://netcologne.dl.sourceforge.net/project/qtx64/qt-x64)
  set(SERVER ${SERVER}/4.8.4/mingw-4.8.0)
  set(MIRRORS
    ${SERVER}/${QT4_MINGW64_7Z}
    ${BOOTSTRAP_PCKG_MIRROR}/${QT4_MINGW64_7Z}
    )
  DOWNLOAD_FILE_FROM_MIRRORS(
    "${MIRRORS}"        
    "${PCKGDIR}/${QT4_MINGW64_7Z}"
    9a17555a2efd3a328899db298d8da22c
    1
    )
endif()
if(NOT EXISTS "${PCKGDIR}/${QT4_MINGW32_EXE}")
  set(SERVER1 ftp://unix.hensa.ac.uk/sites/downloads.sourceforge.net)
  set(SERVER1 ${SERVER1}/q/qt/qtx64/qt-x86/4.8.5/mingw-4.8/sjlj)
  set(SERVER2 ftp://www.mirrorservice.org/sites/downloads.sourceforge.net)
  set(SERVER2 ${SERVER2}/q/qt/qtx64/qt-x86/4.8.5/mingw-4.8/sjlj)
  
  set(MIRRORS
    ${SERVER1}/${QT4_MINGW32_EXE}
    ${SERVER2}/${QT4_MINGW32_EXE}
    ${BOOTSTRAP_PCKG_MIRROR}/${QT4_MINGW32_EXE}
    )
  DOWNLOAD_FILE_FROM_MIRRORS(
    "${MIRRORS}"
    "${PCKGDIR}/${QT4_MINGW32_EXE}"
    c9519f81ebe7ad7942a83fbcb7794a65
    1
    )
endif()
if(NOT EXISTS "${PCKGDIR}/${INNOUNP_RAR}")
  # Download unpacker for InnoSetup files.
  set(SERVER http://netcologne.dl.sourceforge.net/project/innounp)
  set(SERVER ${SERVER}/innounp/innounp%200.40)
  set(MIRRORS
    ${SERVER}/${INNOUNP_RAR}
    ${BOOTSTRAP_PCKG_MIRROR}/${INNOUNP_RAR}
    )
  DOWNLOAD_FILE_FROM_MIRRORS(
    "${MIRRORS}"    
    "${PCKGDIR}/${INNOUNP_RAR}"
    b30b12f837921c3e61e6b21176804cd5
    1
    )
endif()

#-----------------------------------------------------------------------------
# Unpack Qt 4.8.4 64-bit
#-----------------------------------------------------------------------------
if(X64)
  if(NOT EXISTS qt-4.8.4-x64-mingw480r2)
    execute_process(
      COMMAND ${SEVENZIP} x -y "${PCKGDIR}/${QT4_MINGW64_7Z}"
      WORKING_DIRECTORY ${BASEDIR}
      )
  endif()
  file(WRITE "${BASEDIR}/qt-4.8.4-x64-mingw480r2/bin/qt.conf"
    "[Paths]\n  Prefix = ${BASEDIR}/qt-4.8.4-x64-mingw480r2")
endif()

#-----------------------------------------------------------------------------
# Unpack Qt 4.8.5 32-bit
#-----------------------------------------------------------------------------
if(NOT EXISTS qt-4.8.5-x86-mingw481r1)
  execute_process(
    COMMAND ${SEVENZIP} x -y "${PCKGDIR}/${INNOUNP_RAR}"
    WORKING_DIRECTORY ${TMPDIR}
    )
  execute_process(
    COMMAND innounp -x "${PCKGDIR}/${QT4_MINGW32_EXE}"
    WORKING_DIRECTORY ${TMPDIR}
    )
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    "{app}/mingw32" "${BASEDIR}/mingw32"
    WORKING_DIRECTORY ${TMPDIR}
    )
  file(REMOVE_RECURSE "${TMPDIR}/{app}/mingw32")
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    "{app}" "${BASEDIR}/qt-4.8.5-x86-mingw481r1"
    WORKING_DIRECTORY ${TMPDIR}
    )
  file(REMOVE_RECURSE ${TMPDIR})
  file(MAKE_DIRECTORY ${TMPDIR})
endif()
file(WRITE "${BASEDIR}/qt-4.8.5-x86-mingw481r1/bin/qt.conf"
  "[Paths]\n  Prefix = ${BASEDIR}/qt-4.8.5-x86-mingw481r1")

set(QT4_BIN64_DIR "${BASEDIR}/qt-4.8.4-x64-mingw480r2/bin")
set(QT4_BIN32_DIR "${BASEDIR}/qt-4.8.5-x86-mingw481r1/bin")

list(APPEND BIN_DIRS
  "Binary directory for 32-bit Qt 4.8.5"
  "${QT4_BIN32_DIR}"
  "Binary directory for 64-bit Qt 4.8.4"
  "${QT4_BIN64_DIR}"
  )

#=============================================================================
# Fix shebang lines in Python config scripts
#=============================================================================
set(PYFILES 2to3 python-config python2.7-config idle pydoc 
  python2-config smtpd.py)
foreach(PYFILE IN ITEMS ${PYFILES})
  if(X64)
    file(READ "${PYTHONHOME_64}/bin/${PYFILE}" CONTENTS)
    string(REPLACE "\n" ";" CONTENTS "${CONTENTS}")
    list(REMOVE_AT CONTENTS 0)
    string(REPLACE ":" "" PYTHONHOME_64_MSYS "${PYTHONHOME_64}")
    set(PYTHONHOME_64_SHEBANG "#!/${PYTHONHOME_64_MSYS}/bin/python2.7.exe")
    list(INSERT CONTENTS 0 "${PYTHONHOME_64_SHEBANG}")
    string(REPLACE ";" "\n" CONTENTS "${CONTENTS}")
    file(WRITE "${PYTHONHOME_64}/bin/${PYFILE}" "${CONTENTS}")
    file(TO_NATIVE_PATH "${PYTHONHOME_64}" PYTHON_NATIVE)
    file(WRITE "${PYTHONHOME_64}/bin/${PYFILE}.cmd"
    "@echo off
${PYTHON_NATIVE}\\bin\\python2.7.exe ${PYTHON_NATIVE}\\bin\\${PYFILE} %*    
    ")
  endif()

  file(READ "${PYTHONHOME_32}/bin/${PYFILE}" CONTENTS)
  string(REPLACE "\n" ";" CONTENTS "${CONTENTS}")
  list(REMOVE_AT CONTENTS 0)
  string(REPLACE ":" "" PYTHONHOME_32_MSYS "${PYTHONHOME_32}")
  set(PYTHONHOME_32_SHEBANG "#!/${PYTHONHOME_32_MSYS}/bin/python2.7.exe")
  list(INSERT CONTENTS 0 "${PYTHONHOME_32_SHEBANG}")
  string(REPLACE ";" "\n" CONTENTS "${CONTENTS}")
  file(WRITE "${PYTHONHOME_32}/bin/${PYFILE}" "${CONTENTS}")
  file(TO_NATIVE_PATH "${PYTHONHOME_32}" PYTHON_NATIVE)
  file(WRITE "${PYTHONHOME_32}/bin/${PYFILE}.cmd"
  "@echo off
${PYTHON_NATIVE}\\bin\\python2.7.exe ${PYTHON_NATIVE}\\bin\\${PYFILE} %*    
  ")
endforeach()
  
#=============================================================================
# Oracle Java JRE for Eclipse and HDFView
#=============================================================================
set(JRE_DIR "${BASEDIR}/jre1.7.0_51")
#-----------------------------------------------------------------------------
# Download JRE
#-----------------------------------------------------------------------------
set(JRE64_TGZ "jre-7u51-windows-x64.tar.gz")
set(JRE32_TGZ "jre-7u51-windows-i586.tar.gz")
if(NOT EXISTS "${PCKGDIR}/${JRE64_TGZ}")
  execute_process(
    COMMAND "${PCKGDIR}/wget.exe" --no-cookies --no-check-certificate
    --header "Cookie: gpw_e24=http%3A%2F%2Fwww.oracle.com%2F"
    "http://download.oracle.com/otn-pub/java/jdk/7u51-b13/${JRE64_TGZ}"
    WORKING_DIRECTORY ${PCKGDIR}
    RESULT_VARIABLE RETVAL
    )
  if(RETVAL)
    set(MIRRORS
      ${BOOTSTRAP_PCKG_MIRROR}/${JRE64_TGZ}
      )
    DOWNLOAD_FILE_FROM_MIRRORS(
      "${MIRRORS}"    
      "${PCKGDIR}/${JRE64_TGZ}"
      1931de2341f22408be9d6639205675c9
      1
      )    
  else()
    file(GLOB JRE64 "${PCKGDIR}/${JRE64_TGZ}*")
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy "${JRE64}" "${PCKGDIR}/${JRE64_TGZ}"
      )
    file(REMOVE "${JRE64}")
  endif()
endif()
if(NOT EXISTS "${PCKGDIR}/${JRE32_TGZ}")
  execute_process(
    COMMAND "${PCKGDIR}/wget.exe" --no-cookies --no-check-certificate
    --header "Cookie: gpw_e24=http%3A%2F%2Fwww.oracle.com%2F"
    "http://download.oracle.com/otn-pub/java/jdk/7u51-b13/${JRE32_TGZ}"
    WORKING_DIRECTORY ${PCKGDIR}
    RESULT_VARIABLE RETVAL
    )
  if(RETVAL)
    set(MIRRORS
      ${BOOTSTRAP_PCKG_MIRROR}/${JRE32_TGZ}
      )
    DOWNLOAD_FILE_FROM_MIRRORS(
      "${MIRRORS}"    
      "${PCKGDIR}/${JRE32_TGZ}"
      3921c19528d180902939b9f4c9ac92f1
      1
      )    
  else()
    file(GLOB JRE32 "${PCKGDIR}/${JRE32_TGZ}*")
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy "${JRE32}" "${PCKGDIR}/${JRE32_TGZ}"
      )
    file(REMOVE "${JRE32}")
  endif()
endif()

#-----------------------------------------------------------------------------
# Unpack JRE
#-----------------------------------------------------------------------------
if(NOT EXISTS "${JRE_DIR}")
  file(MAKE_DIRECTORY ${JRE_DIR})
  if(X64)
    execute_process(
      COMMAND ${SEVENZIP} x "${PCKGDIR}/${JRE64_TGZ}"
      WORKING_DIRECTORY ${TMPDIR}
      )
    execute_process(
      COMMAND ${SEVENZIP} x jre-7u51-windows-x64.tar
      WORKING_DIRECTORY ${TMPDIR}
      )
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy_directory
      ${TMPDIR}/jre1.7.0_51 ${JRE_DIR}/x64
      WORKING_DIRECTORY ${TMPDIR}
      )
    file(REMOVE_RECURSE ${TMPDIR})
    file(MAKE_DIRECTORY ${TMPDIR})
  endif()
  execute_process(
    COMMAND ${SEVENZIP} x "${PCKGDIR}/${JRE32_TGZ}"
    WORKING_DIRECTORY ${TMPDIR}
    )
  execute_process(
    COMMAND ${SEVENZIP} x jre-7u51-windows-i586.tar
    WORKING_DIRECTORY ${TMPDIR}
    )
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${TMPDIR}/jre1.7.0_51 ${JRE_DIR}/x86
    WORKING_DIRECTORY ${TMPDIR}
    )
  file(REMOVE_RECURSE ${TMPDIR})
  file(MAKE_DIRECTORY ${TMPDIR})
endif()

list(APPEND BIN_DIRS
  "Binary directory for 32-bit JRE"
  "${JRE_DIR}/x86/bin"
  "Binary directory for 64-bit JRE"
  "${JRE_DIR}/x64/bin"
  )

#=============================================================================
# Eclipse Kepler SR1 with CDT
#=============================================================================
set(ECLIPSE_DIR "${BASEDIR}/eclipse")
#-----------------------------------------------------------------------------
# Download Eclipse
#-----------------------------------------------------------------------------
set(ECLIPSE64_ZIP "eclipse-cpp-kepler-SR1-win32-x86_64.zip")
set(ECLIPSE32_ZIP "eclipse-cpp-kepler-SR1-win32.zip")
if(NOT EXISTS "${PCKGDIR}/${ECLIPSE64_ZIP}")
  set(SERVER1 http://ftp.fau.de/eclipse/technology)
  set(SERVER1 ${SERVER1}/epp/downloads/release/kepler/SR1)
  set(SERVER2 http://ftp.halifax.rwth-aachen.de/eclipse/technology)
  set(SERVER2 ${SERVER2}/epp/downloads/release/kepler/SR1)
  
  set(MIRRORS
    ${SERVER1}/${ECLIPSE64_ZIP}
    ${SERVER2}/${ECLIPSE64_ZIP}
    ${BOOTSTRAP_PCKG_MIRROR}/${ECLIPSE64_ZIP}
    )
  DOWNLOAD_FILE_FROM_MIRRORS(
    "${MIRRORS}"
    "${PCKGDIR}/${ECLIPSE64_ZIP}"
    d35cf915b31f7e1068b264cff25119c1
    1
    )
endif()
if(NOT EXISTS "${PCKGDIR}/${ECLIPSE32_ZIP}")
  set(SERVER1 http://ftp.wh2.tu-dresden.de/pub/mirrors/eclipse/technology)
  set(SERVER1 ${SERVER1}/epp/downloads/release/kepler/SR1)
  set(SERVER2 http://ftp.fau.de/eclipse/technology)
  set(SERVER2 ${SERVER2}/epp/downloads/release/kepler/SR1)
  
  set(MIRRORS
    ${SERVER1}/${ECLIPSE32_ZIP}
    ${SERVER2}/${ECLIPSE32_ZIP}
    ${BOOTSTRAP_PCKG_MIRROR}/${ECLIPSE32_ZIP}
    )
  DOWNLOAD_FILE_FROM_MIRRORS(
    "${MIRRORS}"
    "${PCKGDIR}/${ECLIPSE32_ZIP}"
    9cca1051a86681491b2921bab91aea6f
    1
    )
endif()

#-----------------------------------------------------------------------------
# Unpack Eclipse
#-----------------------------------------------------------------------------
if(NOT EXISTS "${ECLIPSE_DIR}")
  file(MAKE_DIRECTORY "${ECLIPSE_DIR}")
  if(X64)
    execute_process(
      COMMAND ${SEVENZIP} x -y "${PCKGDIR}/${ECLIPSE64_ZIP}"
      WORKING_DIRECTORY ${ECLIPSE_DIR}
      )
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy_directory eclipse x64
      WORKING_DIRECTORY ${ECLIPSE_DIR}
      )
    file(REMOVE_RECURSE ${ECLIPSE_DIR}/eclipse)
  endif()
  execute_process(
    COMMAND ${SEVENZIP} x -y "${PCKGDIR}/${ECLIPSE32_ZIP}"
    WORKING_DIRECTORY ${ECLIPSE_DIR}
    )
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E copy_directory eclipse x86
    WORKING_DIRECTORY ${ECLIPSE_DIR}
    )
  file(REMOVE_RECURSE ${ECLIPSE_DIR}/eclipse)
endif()

list(APPEND BIN_DIRS
  "Binary directory for 32-bit Eclipse"
  "${ECLIPSE_DIR}/x86"
  "Binary directory for 64-bit Eclipse"
  "${ECLIPSE_DIR}/x64"
  )

#=============================================================================
# HDFView HDF5 viewer
#=============================================================================
set(HDFVIEW_DIR "${BASEDIR}/hdfview")
#-----------------------------------------------------------------------------
# Download hdfview
#-----------------------------------------------------------------------------
set(HDFVIEW64_ZIP "hdf-java-2.10-win64.zip")
set(HDFVIEW32_ZIP "hdf-java-2.10-win32.zip")
if(NOT EXISTS "${PCKGDIR}/${HDFVIEW64_ZIP}")
  set(MIRRORS
    http://hdfgroup.org/ftp/HDF5/hdf-java/current/bin/${HDFVIEW64_ZIP}
    ${BOOTSTRAP_PCKG_MIRROR}/${HDFVIEW64_ZIP}
    )
  DOWNLOAD_FILE_FROM_MIRRORS(
    "${MIRRORS}"
    "${PCKGDIR}/${HDFVIEW64_ZIP}"
    b5c3154ee63485770b3ff9cd24d03152
    1
    )
endif()
if(NOT EXISTS "${PCKGDIR}/${HDFVIEW32_ZIP}")
  set(MIRRORS
    http://hdfgroup.org/ftp/HDF5/hdf-java/current/bin/${HDFVIEW32_ZIP}
    ${BOOTSTRAP_PCKG_MIRROR}/${HDFVIEW32_ZIP}
    )
  DOWNLOAD_FILE_FROM_MIRRORS(
    "${MIRRORS}"
    "${PCKGDIR}/${HDFVIEW32_ZIP}"
    b1e3649753f54db89dffc788cbd93b45
    1
    )
endif()

#-----------------------------------------------------------------------------
# Unpack HDFView
#-----------------------------------------------------------------------------
if(NOT EXISTS "${HDFVIEW_DIR}")
  file(MAKE_DIRECTORY "${HDFVIEW_DIR}")
  if(X64)
    execute_process(
      COMMAND ${SEVENZIP} x -y "${PCKGDIR}/${HDFVIEW64_ZIP}"
      WORKING_DIRECTORY ${TMPDIR}
      )
    execute_process(
      COMMAND ${SEVENZIP} x HDF-JAVA-2.10.0-win64.exe
      WORKING_DIRECTORY ${TMPDIR}
      )
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy_directory
      ${TMPDIR}/\$_OUTDIR ${HDFVIEW_DIR}/x64
      WORKING_DIRECTORY ${TMPDIR}
      )
    file(READ "${HDFVIEW_DIR}/x64/bin/hdfview.bat" HDFVIEW_BAT)
    string(REPLACE "C:/Program Files/HDF_Group/HDF-JAVA/2.10.0"
      "${HDFVIEW_DIR}/x64" HDFVIEW_BAT "${HDFVIEW_BAT}")
    string(REPLACE "%INSTALLDIR%\\bin\\JRE\\bin"
      "${JRE_DIR}/x64/bin" HDFVIEW_BAT "${HDFVIEW_BAT}")
    file(WRITE "${HDFVIEW_DIR}/x64/bin/hdfview.bat" "new_${HDFVIEW_BAT}")
    file(REMOVE_RECURSE ${TMPDIR})
    file(MAKE_DIRECTORY ${TMPDIR})
  endif()
  execute_process(
    COMMAND ${SEVENZIP} x -y "${PCKGDIR}/${HDFVIEW32_ZIP}"
    WORKING_DIRECTORY ${TMPDIR}
    )
  execute_process(
    COMMAND ${SEVENZIP} x HDF-JAVA-2.10.0-win32.exe
    WORKING_DIRECTORY ${TMPDIR}
    )
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${TMPDIR}/\$_OUTDIR ${HDFVIEW_DIR}/x86
    WORKING_DIRECTORY ${TMPDIR}
    )
  file(READ "${HDFVIEW_DIR}/x86/bin/hdfview.bat" HDFVIEW_BAT)
  string(REPLACE "C:/Program Files/HDF_Group/HDF-JAVA/2.10.0"
    "${HDFVIEW_DIR}/x86" HDFVIEW_BAT "${HDFVIEW_BAT}")
  string(REPLACE "%INSTALLDIR%\\bin\\JRE\\bin"
    "${JRE_DIR}/x86/bin" HDFVIEW_BAT "${HDFVIEW_BAT}")
  file(WRITE "${HDFVIEW_DIR}/x86/bin/hdfview.bat" "new_${HDFVIEW_BAT}")
  file(REMOVE_RECURSE ${TMPDIR})
  file(MAKE_DIRECTORY ${TMPDIR})
endif()

list(APPEND BIN_DIRS
  "Binary directory for 32-bit HDFView"
  "${HDFVIEW_DIR}/x86/bin"
  "Binary directory for 64-bit HDFView"
  "${HDFVIEW_DIR}/x64/bin"
  )

#=============================================================================
# CMake with configurable default generator. On Windows CMake searches for
# an installed MS Visual Studio by default. This can be very annoying if one
# forgets to tell CMake the proper generator on the command line. Therefore
# we patch the CMake source, so that the default generator can be specified
# in an environment variable.
#=============================================================================
set(CMAKE_DIR "${BASEDIR}/cmake-with-default-generator")
#-----------------------------------------------------------------------------
# Download CMake sources
#-----------------------------------------------------------------------------
set(CMAKE_TGZ "cmake-2.8.12.2.tar.gz")
if(NOT EXISTS "${PCKGDIR}/${CMAKE_TGZ}")
  set(MIRRORS
    http://www.cmake.org/files/v2.8/${CMAKE_TGZ}
    ${BOOTSTRAP_PCKG_MIRROR}/${CMAKE_TGZ}
    )
  DOWNLOAD_FILE_FROM_MIRRORS(
    "${MIRRORS}"
    "${PCKGDIR}/${CMAKE_TGZ}"
    17c6513483d23590cbce6957ec6d1e66
    1
    )
endif()

#-----------------------------------------------------------------------------
# Patch for making default generator configurable by the environment
# variable CMAKE_DEFAULT_GENERATOR.
#-----------------------------------------------------------------------------
set(CMAKE_DEFAULT_GEN_PATCH 
"--- Source/cmake.cxx.orig	Tue Nov  5 19:07:23 2013
+++ Source/cmake.cxx	Thu Feb 27 18:31:30 2014
@@ -2234,47 +2234,53 @@
 #if defined(__BORLANDC__) && defined(_WIN32)
       this->SetGlobalGenerator(new cmGlobalBorlandMakefileGenerator);
 #elif defined(_WIN32) && !defined(__CYGWIN__) && !defined(CMAKE_BOOT_MINGW)
-      std::string installedCompiler;
-      // Try to find the newest VS installed on the computer and
-      // use that as a default if -G is not specified
-      const std::string vsregBase =
-        \"[HKEY_LOCAL_MACHINE\\\\SOFTWARE\\\\Microsoft\\\\\";
-      std::vector<std::string> vsVerions;
-      vsVerions.push_back(\"VisualStudio\\\\\");
-      vsVerions.push_back(\"VCExpress\\\\\");
-      vsVerions.push_back(\"WDExpress\\\\\");
-      struct VSRegistryEntryName
-      {
-        const char* MSVersion;
-        const char* GeneratorName;
-      };
-      VSRegistryEntryName version[] = {
-        {\"6.0\", \"Visual Studio 6\"},
-        {\"7.0\", \"Visual Studio 7\"},
-        {\"7.1\", \"Visual Studio 7 .NET 2003\"},
-        {\"8.0\", \"Visual Studio 8 2005\"},
-        {\"9.0\", \"Visual Studio 9 2008\"},
-        {\"10.0\", \"Visual Studio 10\"},
-        {\"11.0\", \"Visual Studio 11\"},
-        {\"12.0\", \"Visual Studio 12\"},
-        {0, 0}};
-      for(int i=0; version[i].MSVersion != 0; i++)
-        {
-        for(size_t b=0; b < vsVerions.size(); b++)
-          {
-          std::string reg = vsregBase + vsVerions[b] + version[i].MSVersion;
-          reg += \";InstallDir]\";
-          cmSystemTools::ExpandRegistryValues(reg,
-                                              cmSystemTools::KeyWOW64_32);
-          if (!(reg == \"/registry\"))
-            {
-            installedCompiler = version[i].GeneratorName;
-            break;
-            }
-          }
-        }
+      const char* genName = getenv(\"CMAKE_DEFAULT_GENERATOR\");
+      if(!genName) {      
+    	  std::string installedCompiler;
+    	  // Try to find the newest VS installed on the computer and
+    	  // use that as a default if -G is not specified
+    	  const std::string vsregBase =
+    			  \"[HKEY_LOCAL_MACHINE\\\\SOFTWARE\\\\Microsoft\\\\\";
+    	  std::vector<std::string> vsVerions;
+    	  vsVerions.push_back(\"VisualStudio\\\\\");
+    	  vsVerions.push_back(\"VCExpress\\\\\");
+    	  vsVerions.push_back(\"WDExpress\\\\\");
+    	  struct VSRegistryEntryName
+    	  {
+    		  const char* MSVersion;
+    		  const char* GeneratorName;
+    	  };
+    	  VSRegistryEntryName version[] = {
+    			  {\"6.0\", \"Visual Studio 6\"},
+    			  {\"7.0\", \"Visual Studio 7\"},
+    			  {\"7.1\", \"Visual Studio 7 .NET 2003\"},
+    			  {\"8.0\", \"Visual Studio 8 2005\"},
+    			  {\"9.0\", \"Visual Studio 9 2008\"},
+    			  {\"10.0\", \"Visual Studio 10\"},
+    			  {\"11.0\", \"Visual Studio 11\"},
+    			  {\"12.0\", \"Visual Studio 12\"},
+    			  {0, 0}};
+    	  for(int i=0; version[i].MSVersion != 0; i++)
+    	  {
+    		  for(size_t b=0; b < vsVerions.size(); b++)
+    		  {
+    			  std::string reg = vsregBase + vsVerions[b] + version[i].MSVersion;
+    			  reg += \";InstallDir]\";
+    			  cmSystemTools::ExpandRegistryValues(reg,
+                                              	  cmSystemTools::KeyWOW64_32);
+    			  if (!(reg == \"/registry\"))
+    			  {
+    				  installedCompiler = version[i].GeneratorName;
+    				  break;
+    			  }
+    		  }
+    	  }
+    	  
+    	  genName = installedCompiler.c_str();
+      }
+      
       cmGlobalGenerator* gen
-        = this->CreateGlobalGenerator(installedCompiler.c_str());
+        = this->CreateGlobalGenerator(genName);
       if(!gen)
         {
         gen = new cmGlobalNMakeMakefileGenerator;"
)
file(WRITE "${PCKGDIR}/cmake-default-gen.patch" "${CMAKE_DEFAULT_GEN_PATCH}")
#-----------------------------------------------------------------------------
# Make sure, that we have unix line endings in patch.
#-----------------------------------------------------------------------------
configure_file(
  "${PCKGDIR}/cmake-default-gen.patch"
  "${TMPDIR}/cmake-default-gen.patch"
  @ONLY NEWLINE_STYLE UNIX)
configure_file(
  "${TMPDIR}/cmake-default-gen.patch"
  "${PCKGDIR}/cmake-default-gen.patch"
  @ONLY NEWLINE_STYLE UNIX)

#-----------------------------------------------------------------------------
# Unpack, patch, build and install CMake
#-----------------------------------------------------------------------------
if(NOT EXISTS "${CMAKE_DIR}")
  file(MAKE_DIRECTORY "${CMAKE_DIR}")
  set(PATH_SAVE $ENV{PATH}) # save old PATH
  set(ENV{PATH} "${BASEDIR}\\Git\\bin;$ENV{PATH}")     # for patch.exe
  set(ENV{PATH} "${BASEDIR}\\mingw32\\bin;$ENV{PATH}") # for GNU compilers
  set(ENV{PATH} "${QT4_BIN32_DIR};$ENV{PATH}")         # for qmake.exe
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E tar xvf "${PCKGDIR}/${CMAKE_TGZ}"
    WORKING_DIRECTORY ${TMPDIR}
    )
  execute_process(
    COMMAND patch --binary -p0
    -i "${TMPDIR}/cmake-default-gen.patch"
    WORKING_DIRECTORY ${TMPDIR}/cmake-2.8.12.2
    )
  file(MAKE_DIRECTORY "${TMPDIR}/build")
  execute_process(
    COMMAND ${CMAKE_COMMAND}
    -G "MinGW Makefiles"
    -DCMAKE_C_COMPILER:FILEPATH=${BASEDIR}/mingw32/bin/gcc.exe
    -DCMAKE_CXX_COMPILER:FILEPATH=${BASEDIR}/mingw32/bin/g++.exe
    -DCMAKE_RC_COMPILER:FILEPATH=${BASEDIR}/mingw32/bin/windres.exe
    -DCMAKE_MAKE_PROGRAM:FILEPATH=${BASEDIR}/mingw32/bin/mingw32-make.exe
    -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_DIR}
    -DBUILD_QtDialog:BOOL=ON
    ../cmake-2.8.12.2
    WORKING_DIRECTORY ${TMPDIR}/build
    )
  execute_process(
    COMMAND mingw32-make -j install
    WORKING_DIRECTORY ${TMPDIR}/build
    )
  file(REMOVE_RECURSE ${TMPDIR})
  file(MAKE_DIRECTORY ${TMPDIR})

  set(QTDLLS libgcc_s_sjlj-1.dll libstdc++-6.dll qtcore4.dll qtgui4.dll)
  foreach(DLL IN ITEMS ${QTDLLS})
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy ${DLL} ${CMAKE_DIR}/bin/${DLL}
      WORKING_DIRECTORY ${QT4_BIN32_DIR}
      )    
  endforeach()

  set(ENV{PATH} ${PATH_SAVE}) # reset old PATH
  file(REMOVE_RECURSE ${TMPDIR})
  file(MAKE_DIRECTORY ${TMPDIR})
endif()
file(REMOVE "${TMPDIR}/cmake-default-gen.patch")

list(APPEND BIN_DIRS
  "Binary directory for CMake with configurable default generator"
  "${CMAKE_DIR}/bin"
  "Binary directory for CMake with configurable default generator"
  "${CMAKE_DIR}/bin"
  )

file(TO_NATIVE_PATH ${BASEDIR} BASEDIR_NATIVE)

# Add binary directories for GCC. This needs to be done here, in order to have
# them in front of the PATH variable.
list(APPEND BIN_DIRS
  "Binary directory for 32-bit GCC 4.8"
  "${BASEDIR}/mingw32/bin"
  "Binary directory for 64-bit GCC 4.8"
  "${BASEDIR}/mingw64/bin"
  "Binary directory for 32-bit Python 2.7"
  "${PYTHONHOME_32}/bin"
  "Binary directory for 64-bit Python 2.7"
  "${PYTHONHOME_64}/bin"
  )

#=============================================================================
# Download and build a few very basic libs (zlib, iconv, libxml2, libxslt,
# TCL) from source, so that they match the GNU compiler.
#=============================================================================
SET(EXTRALIBS_DIR "${BASEDIR}/extralibs")

#-----------------------------------------------------------------------------
# Download sources for extra libs
#-----------------------------------------------------------------------------
set(ZLIB_TGZ "zlib-1.2.8.tar.gz")
set(ICONV_TGZ "libiconv-1.14.tar.gz")
set(LIBXML2_TGZ "libxml2-2.9.1.tar.gz")
set(LIBXSLT_TGZ "libxslt-1.1.28.tar.gz")
set(PYGMENTS_TGZ "Pygments-1.6.tar.gz")
set(TCL_ZIP "tcl861-src.zip")
set(TK_ZIP "tk861-src.zip")

if(NOT EXISTS "${PCKGDIR}/${ZLIB_TGZ}")
  set(MIRRORS
    ftp://ftp.pl.pgpi.org/vol/rzm1/GraphicsMagick/delegates/${ZLIB_TGZ}
    ftp://ftp.uwsg.indiana.edu/linux/gentoo/distfiles/${ZLIB_TGZ}
    ${BOOTSTRAP_PCKG_MIRROR}/${CMAKE_TGZ}
    )
  
  DOWNLOAD_FILE_FROM_MIRRORS(
    "${MIRRORS}"
    "${PCKGDIR}/${ZLIB_TGZ}"
    44d667c142d7cda120332623eab69f40
    1
    )
endif()
if(NOT EXISTS "${PCKGDIR}/${ICONV_TGZ}")
  set(MIRRORS
    http://ftp.gnu.org/gnu/libiconv/${ICONV_TGZ}
    ${BOOTSTRAP_PCKG_MIRROR}/${ICONV_TGZ}
    )
  
  DOWNLOAD_FILE_FROM_MIRRORS(
    "${MIRRORS}"
    "${PCKGDIR}/${ICONV_TGZ}"
    e34509b1623cec449dfeb73d7ce9c6c6
    1
    )
endif()
if(NOT EXISTS "${PCKGDIR}/${LIBXML2_TGZ}")
  set(MIRRORS
    ftp://xmlsoft.org/libxslt/${LIBXML2_TGZ}
    ${BOOTSTRAP_PCKG_MIRROR}/${LIBXML2_TGZ}
    )
  
  DOWNLOAD_FILE_FROM_MIRRORS(
    "${MIRRORS}"
    "${PCKGDIR}/${LIBXML2_TGZ}"
    9c0cfef285d5c4a5c80d00904ddab380
    1
    )
endif()
if(NOT EXISTS "${PCKGDIR}/${LIBXSLT_TGZ}")
  set(MIRRORS
    ftp://xmlsoft.org/libxslt/${LIBXSLT_TGZ}
    ${BOOTSTRAP_PCKG_MIRROR}/${LIBXSLT_TGZ}
    )
  
  DOWNLOAD_FILE_FROM_MIRRORS(
    "${MIRRORS}"
    "${PCKGDIR}/${LIBXSLT_TGZ}"
    9667bf6f9310b957254fdcf6596600b7
    1
    )
endif()
if(NOT EXISTS "${PCKGDIR}/${PYGMENTS_TGZ}")
  set(MD5 a18feedf6ffd0b0cc8c8b0fbdb2027b1)
  set(MIRRORS
    https://pypi.python.org/packages/source/P/Pygments/${PYGMENTS_TGZ}\#md5=${MD5}
    ${BOOTSTRAP_PCKG_MIRROR}/${PYGMENTS_TGZ}
    )
  
  DOWNLOAD_FILE_FROM_MIRRORS(
    "${MIRRORS}"
    "${PCKGDIR}/${PYGMENTS_TGZ}"
    a18feedf6ffd0b0cc8c8b0fbdb2027b1
    1
    )
endif()
if(NOT EXISTS "${PCKGDIR}/${TCL_ZIP}")
  set(MIRRORS
    http://garr.dl.sourceforge.net/project/tcl/Tcl/8.6.1/${TCL_ZIP}
    ${BOOTSTRAP_PCKG_MIRROR}/${TCL_ZIP}
    )
  
  DOWNLOAD_FILE_FROM_MIRRORS(
    "${MIRRORS}"
    "${PCKGDIR}/${TCL_ZIP}"
    3573f2e35e7bc2dd2fe4d00be16ad767
    1
    )
endif()
if(NOT EXISTS "${PCKGDIR}/${TK_ZIP}")
  set(MIRRORS
    http://garr.dl.sourceforge.net/project/tcl/Tcl/8.6.1/${TK_ZIP}
    ${BOOTSTRAP_PCKG_MIRROR}/${TK_ZIP}
    )
  
  DOWNLOAD_FILE_FROM_MIRRORS(
    "${MIRRORS}"
    "${PCKGDIR}/${TK_ZIP}"
    8899c37f9f3c8bd6b712fc8ab89d1184
    1
    )
endif()
  
#-----------------------------------------------------------------------------
# Unpack and patch sources and build extra libs.
#-----------------------------------------------------------------------------
macro(build_extra_libs BITS ARCH)
  set(BITS "${BITS}")
  set(ARCH "${ARCH}")

  execute_process(
    COMMAND ${CMAKE_COMMAND} -E tar xvf "${PCKGDIR}/${ZLIB_TGZ}"
    WORKING_DIRECTORY ${TMPDIR}
    )
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E tar xvf "${PCKGDIR}/${ICONV_TGZ}"
    WORKING_DIRECTORY ${TMPDIR}
    )
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E tar xvf "${PCKGDIR}/${LIBXML2_TGZ}"
    WORKING_DIRECTORY ${TMPDIR}
    )
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E tar xvf "${PCKGDIR}/${LIBXSLT_TGZ}"
    WORKING_DIRECTORY ${TMPDIR}
    )
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E tar xvf "${PCKGDIR}/${PYGMENTS_TGZ}"
    WORKING_DIRECTORY ${TMPDIR}
    )
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E tar xvf "${PCKGDIR}/${TCL_ZIP}"
    WORKING_DIRECTORY ${TMPDIR}
    )
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E tar xvf "${PCKGDIR}/${TK_ZIP}"
    WORKING_DIRECTORY ${TMPDIR}
    )

  set(XML_PATCHES 
"--- libxml2-2.9.1/win32/Makefile.mingw.orig	Wed Feb 27 06:08:51 2013
+++ libxml2-2.9.1/win32/Makefile.mingw	Sun Mar 02 02:25:47 2014
@@ -40,7 +40,7 @@
 
 # The compiler and its options.
 CC = gcc.exe
-CFLAGS += -DWIN32 -D_WINDOWS -D_MBCS -DNOLIBTOOL 
+CFLAGS += -DWIN32 -D_WIN32 -D_WINDOWS -D_MBCS -DNOLIBTOOL 
 CFLAGS += -I\$(XML_SRCDIR) -I\$(XML_SRCDIR)/include -I\$(INCPREFIX) \$(INCLUDE)
 ifneq (\$(WITH_THREADS),no)
 CFLAGS += -D_REENTRANT
--- libxslt-1.1.28/win32/Makefile.mingw.orig	Tue Sep 04 15:26:23 2012
+++ libxslt-1.1.28/win32/Makefile.mingw	Sun Mar 02 20:00:49 2014
@@ -49,7 +49,7 @@
 
 # The compiler and its options.
 CC = gcc.exe
-CFLAGS += -DWIN32 -D_WINDOWS -D_MBCS
+CFLAGS += -DWIN32 -D_WIN32 -D_WINDOWS -D_MBCS
 CFLAGS += -I\$(BASEDIR) -I\$(XSLT_SRCDIR) -I\$(INCPREFIX)
 
 # The linker and its options.
@@ -72,6 +72,7 @@
 
 # Libxslt object files.
 XSLT_OBJS = \$(XSLT_INTDIR)/attributes.o\\
+        \$(XSLT_INTDIR)/attrvt.o\\
 	\$(XSLT_INTDIR)/documents.o\\
 	\$(XSLT_INTDIR)/extensions.o\\
 	\$(XSLT_INTDIR)/extra.o\\
@@ -93,6 +94,7 @@
 
 # Static libxslt object files.
 XSLT_OBJS_A = \$(XSLT_INTDIR_A)/attributes.o\\
+        \$(XSLT_INTDIR_A)/attrvt.o\\
 	\$(XSLT_INTDIR_A)/documents.o\\
 	\$(XSLT_INTDIR_A)/extensions.o\\
 	\$(XSLT_INTDIR_A)/extra.o\\
@@ -173,17 +175,17 @@
 	cmd.exe /C if not exist \$(INCPREFIX)\\\$(EXSLT_BASENAME) mkdir \$(INCPREFIX)\\\$(EXSLT_BASENAME)
 	cmd.exe /C if not exist \$(BINPREFIX) mkdir \$(BINPREFIX)
 	cmd.exe /C if not exist \$(LIBPREFIX) mkdir \$(LIBPREFIX)
-	cmd.exe /C copy \$(XSLT_SRCDIR)\\*.h \$(INCPREFIX)\\\$(XSLT_BASENAME)
-	cmd.exe /C copy \$(EXSLT_SRCDIR)\\*.h \$(INCPREFIX)\\\$(EXSLT_BASENAME)
-	cmd.exe /C copy \$(BINDIR)\\\$(XSLT_SO) \$(SOPREFIX)
-	cmd.exe /C copy \$(BINDIR)\\\$(XSLT_A) \$(LIBPREFIX)
-	cmd.exe /C copy \$(BINDIR)\\\$(XSLT_IMP) \$(LIBPREFIX)
-	cmd.exe /C copy \$(BINDIR)\\\$(EXSLT_SO) \$(SOPREFIX)
-	cmd.exe /C copy \$(BINDIR)\\\$(EXSLT_A) \$(LIBPREFIX)
-	cmd.exe /C copy \$(BINDIR)\\\$(EXSLT_IMP) \$(LIBPREFIX)
+	cp \$(XSLT_SRCDIR)\\\\*.h \$(INCPREFIX)\\\$(XSLT_BASENAME)
+	cp \$(EXSLT_SRCDIR)\\\\*.h \$(INCPREFIX)\\\$(EXSLT_BASENAME)
+	cp \$(BINDIR)\\\$(XSLT_SO) \$(BINPREFIX)
+	cp \$(BINDIR)\\\$(XSLT_A) \$(LIBPREFIX)
+	cp \$(BINDIR)\\\$(XSLT_IMP) \$(LIBPREFIX)
+	cp \$(BINDIR)\\\$(EXSLT_SO) \$(BINPREFIX)
+	cp \$(BINDIR)\\\$(EXSLT_A) \$(LIBPREFIX)
+	cp \$(BINDIR)\\\$(EXSLT_IMP) \$(LIBPREFIX)
 
 install : install-libs
-	cmd.exe /C copy \$(BINDIR)\\*.exe \$(BINPREFIX)
+	cmd.exe /C copy \$(BINDIR)\\\\*.exe \$(BINPREFIX)
 
 install-dist : install
 
@@ -280,8 +282,8 @@
 	cmd.exe /C if not exist \$(UTILS_INTDIR) mkdir \$(UTILS_INTDIR)
 
 # An implicit rule for xsltproc and friends.
-APPLIBS = \$(LIBS)
-APPLIBS += -llibxml2 -l\$(XSLT_BASENAME) -l\$(EXSLT_BASENAME)
+APPLIBS = -l\$(XSLT_BASENAME) -l\$(EXSLT_BASENAME) -llibxml2
+APPLIBS += \$(LIBS)
 APP_LDFLAGS = \$(LDFLAGS)
 APP_LDFLAGS += -Wl,--major-image-version,\$(LIBXSLT_MAJOR_VERSION)
 APP_LDFLAGS += -Wl,--minor-image-version,\$(LIBXSLT_MINOR_VERSION)
@@ -290,7 +292,7 @@
 APP_LDFLAGS += -Bstatic
 \$(BINDIR)/%.exe : \$(UTILS_SRCDIR)/%.c
 	\$(CC) \$(CFLAGS) -o \$(subst .c,.o,\$(UTILS_INTDIR)/\$(<F)) -c \$< 
-	\$(LD) \$(APP_LDFLAGS) -o \$@ \$(APPLIBS) \$(subst .c,.o,\$(UTILS_INTDIR)/\$(<F))
+	\$(LD) \$(APP_LDFLAGS) -o \$@ \$(subst .c,.o,\$(UTILS_INTDIR)/\$(<F)) \$(APPLIBS)
 else
 \$(BINDIR)/%.exe : \$(UTILS_SRCDIR)/%.c 
 	\$(CC) \$(CFLAGS) -o \$(subst .c,.o,\$(UTILS_INTDIR)/\$(<F)) -c \$< \n")
  file(WRITE "${PCKGDIR}/libxml2-libxslt.patch" "${XML_PATCHES}")
  #---------------------------------------------------------------------------
  # Make sure, that we have unix line endings in patch.
  #---------------------------------------------------------------------------
  configure_file(
    "${PCKGDIR}/libxml2-libxslt.patch"
    "${TMPDIR}/libxml2-libxslt.patch"
    @ONLY NEWLINE_STYLE UNIX)
  configure_file(
    "${TMPDIR}/libxml2-libxslt.patch"
    "${PCKGDIR}/libxml2-libxslt.patch"
    @ONLY NEWLINE_STYLE UNIX)

  execute_process(
    COMMAND ${BASEDIR}/Git/bin/patch.exe -p0 --binary
    -i ${PCKGDIR}/libxml2-libxslt.patch
    WORKING_DIRECTORY ${TMPDIR}
    )  
  file(REMOVE "${TMPDIR}/tcl8.6.1/compat/zlib/win32/zdll.lib")
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E copy
    zlib1.dll "${UNXUTILS_BIN_DIR}/../bin/zlib1.dll"
    WORKING_DIRECTORY ${TMPDIR}/tcl8.6.1/compat/zlib/win32
    )  
  file(REMOVE "${TMPDIR}/tcl8.6.1/compat/zlib/win64/zdll.lib")

  set(BUILD_CMD "@echo off

set EXTRADIR=extralibs
set BASEDIR_NATIVE=${BASEDIR_NATIVE}
set BASEDIR_CMAKE=${BASEDIR}

call \"%BASEDIR_NATIVE%\\mingw@BITS@_env.cmd\"

echo Building zlib...
cd /d %BASEDIR_NATIVE%\\tmp\\zlib-1.2.8
mkdir build@BITS@ && cd build@BITS@
cmake -DCMAKE_INSTALL_PREFIX=%BASEDIR_CMAKE%/%EXTRADIR%/@ARCH@ ..
make install
cd /d %BASEDIR_NATIVE%\\%EXTRADIR%\\x86\\bin
copy libzlib.dll \"%BASEDIR_NATIVE%\\UnxUtils\\usr\\local\\bin\\zlib1.dll\"
cd /d %BASEDIR_NATIVE%\\tmp\\tcl8.6.1\\compat\\zlib\\win32
copy %BASEDIR_NATIVE%\\%EXTRADIR%\\@ARCH@\\lib\\libzlib.dll.a zdll.lib
cd /d %BASEDIR_NATIVE%\\tmp\\tcl8.6.1\\compat\\zlib\\win64
copy %BASEDIR_NATIVE%\\%EXTRADIR%\\@ARCH@\\lib\\libzlib.dll.a zdll.lib

echo Building libiconv...
cd /d %BASEDIR_NATIVE%\\tmp\\libiconv-1.14
>build.sh ECHO ./configure --prefix=\"\$BASEDIR/%EXTRADIR%/@ARCH@\"
>>build.sh ECHO make install
%BASEDIR_NATIVE%\\Git\\bin\\bash.exe --login -i build.sh

echo Building libxml2...
cd /d %BASEDIR_NATIVE%\\tmp\\libxml2-2.9.1\\win32
cscript configure.js prefix=%BASEDIR_NATIVE%\\%EXTRADIR%\\@ARCH@ zlib=yes python=yes static=yes compiler=mingw
make -f Makefile.mingw install
cd /d %BASEDIR_NATIVE%\\%EXTRADIR%\\@ARCH@\\include
move libxml2\\libxml .
rmdir /s/q libxml2 

echo Building libxslt...
cd /d %BASEDIR_NATIVE%\\tmp\\libxslt-1.1.28\\win32
cscript configure.js prefix=%BASEDIR_NATIVE%\\%EXTRADIR%\\@ARCH@ include=%BASEDIR_NATIVE%\\%EXTRADIR%\\@ARCH@\\include lib=%BASEDIR_NATIVE%\\%EXTRADIR%\\@ARCH@\\lib zlib=yes  static=yes compiler=mingw
make -f Makefile.mingw install

echo Building Pygments...
cd /d %BASEDIR_NATIVE%\\tmp\\Pygments-1.6
python setup.py build
python setup.py install
set PYG_CMD=%PYTHONHOME%\\bin\\pygmentize.cmd
>%PYG_CMD% ECHO ^@echo off
>>%PYG_CMD% ECHO \"%PYTHONHOME%\\bin\\python.exe\" \"%PYTHONHOME%\\Scripts\\pygmentize\" %%*

set PATH=%BASEDIR_NATIVE%\\bin;%PATH%

cd /d %BASEDIR_NATIVE%\\tmp\\tcl8.6.1\\win
>build.sh ECHO ./configure --prefix=\"\$BASEDIR/%EXTRADIR%/@ARCH@\"
>>build.sh ECHO make install
%BASEDIR_NATIVE%\\Git\\bin\\bash.exe --login -i build.sh

cd /d %BASEDIR_NATIVE%\\tmp\\tk8.6.1\\win
>build.sh ECHO ./configure --prefix=\"\$BASEDIR/%EXTRADIR%/@ARCH@\"
>>build.sh ECHO make install
copy *.h %BASEDIR_NATIVE%\\%EXTRADIR%\\@ARCH@\\include
copy ..\\generic\\*.h %BASEDIR_NATIVE%\\%EXTRADIR%\\@ARCH@\\include
%BASEDIR_NATIVE%\\Git\\bin\\bash.exe --login -i build.sh
\n")
  string(CONFIGURE "${BUILD_CMD}" BUILD_CMD32 @ONLY)
  file(WRITE "${TMPDIR}/build.cmd" "${BUILD_CMD32}")

  execute_process(
    COMMAND cmd /C build.cmd
    WORKING_DIRECTORY ${TMPDIR}
    )

  execute_process(
    COMMAND ${CMAKE_COMMAND} -E copy tclsh86.exe tclsh.exe
    WORKING_DIRECTORY "${EXTRALIBS_DIR}/${ARCH}/bin"
    )  
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E copy wish86.exe wish.exe
    WORKING_DIRECTORY "${EXTRALIBS_DIR}/${ARCH}/bin"
    )  

  file(REMOVE_RECURSE ${TMPDIR})
  file(MAKE_DIRECTORY ${TMPDIR})
endmacro(build_extra_libs)

list(APPEND BIN_DIRS
  "32-bit Binary directory for extra libs (iconv, zlib, libxml2, TCL)"
  "${EXTRALIBS_DIR}/x86/bin"
  "64-bit Binary directory for extra libs (iconv, zlib, libxml2, TCL)"
  "${EXTRALIBS_DIR}/x64/bin"
  )

# Write intermediate version of MinGW environments.
write_mingw_env()

if(NOT EXISTS ${EXTRALIBS_DIR})
  build_extra_libs(32 x86)
  if(X64)
    build_extra_libs(64 x64)
  endif()
endif()

endif()

#=============================================================================
# Console2 - an improved windows command line shell.
#=============================================================================
SET(CONSOLE2_DIR "${BASEDIR}/Console2")

#-----------------------------------------------------------------------------
# Download Console2
#-----------------------------------------------------------------------------
set(CONSOLE2_ZIP "Console-2.00b148-Beta_32bit.zip")
if(NOT EXISTS "${PCKGDIR}/${CONSOLE2_ZIP}")
  set(MIRRORS
    http://netcologne.dl.sourceforge.net/project/console/console-devel/2.00/Console-2.00b148-Beta_32bit.zip
    ftp://ftp.pl.pgpi.org/vol/rzm1/GraphicsMagick/delegates/${ZLIB_TGZ}
    ftp://ftp.uwsg.indiana.edu/linux/gentoo/distfiles/${ZLIB_TGZ}
    ${BOOTSTRAP_PCKG_MIRROR}/${CMAKE_TGZ}
    )
  
  DOWNLOAD_FILE_FROM_MIRRORS(
    "${MIRRORS}"
    "${PCKGDIR}/${CONSOLE2_ZIP}"
    5167a73fc962eac4b0ef9f7958e8c847
    1
    )
endif()

set(MSYS_ICO "msys.ico")
if(NOT EXISTS "${PCKGDIR}/${MSYS_ICO}")
  execute_process(
    COMMAND "${PCKGDIR}/wget.exe" --no-check-certificate
    https://raw.github.com/OpenModelica/OMDev/master/tools/msys/${MSYS_ICO}
    WORKING_DIRECTORY ${PCKGDIR}
    RESULT_VARIABLE RETVAL
    )
  if(RETVAL)
    set(MIRRORS
      ${BOOTSTRAP_PCKG_MIRROR}/${MSYS_ICO}
      )
    DOWNLOAD_FILE_FROM_MIRRORS(
      "${MIRRORS}"
      "${PCKGDIR}/${MSYS_ICO}"
      6999c504c2fddc295fb9a96a55baa40b
      1
      )
  endif()
endif()

#-----------------------------------------------------------------------------
# Unpack Console2
#-----------------------------------------------------------------------------
if(NOT EXISTS ${CONSOLE2_DIR})
  execute_process(COMMAND ${SEVENZIP} x -y "${PCKGDIR}/${CONSOLE2_ZIP}"
    WORKING_DIRECTORY ${TMPDIR}
    )
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E copy_directory "Console2" "${CONSOLE2_DIR}"
    WORKING_DIRECTORY ${TMPDIR}
    )
  file(MAKE_DIRECTORY "$ENV{APPDATA}/Console")
  execute_process(
    COMMAND ${SEVENZIP} x -y $ENV{SystemRoot}\\system32\\cmd.exe
    WORKING_DIRECTORY ${TMPDIR}
    )
  execute_process(
    COMMAND "${IMGMGCK_DIR}/convert.exe"
    3.ico 4.ico 5.ico
    -alpha on -colors 256
    "${CONSOLE2_DIR}/cmd.ico"
    WORKING_DIRECTORY ${TMPDIR}/.rsrc/ICON
    )
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E copy
    "${PCKGDIR}/msys.ico" "${CONSOLE2_DIR}/msys.ico"
    WORKING_DIRECTORY ${TMPDIR}
    )
  file(TO_NATIVE_PATH "${CONSOLE2_DIR}" CONSOLE2_NATIVE_DIR)
  set(CONSOLE_XML
    "
<?xml version=\"1.0\"?>
<settings>
  <console change_refresh=\"10\" refresh=\"100\" rows=\"25\" columns=\"80\"
           buffer_rows=\"500\" buffer_columns=\"0\" shell=\"\"
           init_dir=\"\" start_hidden=\"0\" save_size=\"0\">
    <colors>
      <color id=\"0\" r=\"7\" g=\"54\" b=\"66\"/>
      <color id=\"1\" r=\"38\" g=\"139\" b=\"210\"/>
      <color id=\"2\" r=\"133\" g=\"153\" b=\"0\"/>
      <color id=\"3\" r=\"42\" g=\"161\" b=\"152\"/>
      <color id=\"4\" r=\"220\" g=\"50\" b=\"47\"/>
      <color id=\"5\" r=\"211\" g=\"54\" b=\"130\"/>
      <color id=\"6\" r=\"181\" g=\"137\" b=\"0\"/>
      <color id=\"7\" r=\"238\" g=\"232\" b=\"213\"/>
      <color id=\"8\" r=\"42\" g=\"161\" b=\"152\"/>
      <color id=\"9\" r=\"131\" g=\"148\" b=\"150\"/>
      <color id=\"10\" r=\"88\" g=\"110\" b=\"117\"/>
      <color id=\"11\" r=\"147\" g=\"161\" b=\"161\"/>
      <color id=\"12\" r=\"203\" g=\"75\" b=\"22\"/>
      <color id=\"13\" r=\"108\" g=\"113\" b=\"196\"/>
      <color id=\"14\" r=\"101\" g=\"123\" b=\"131\"/>
      <color id=\"15\" r=\"253\" g=\"246\" b=\"227\"/>
    </colors>
  </console>
  <appearance>
    <font name=\"Courier New\" size=\"10\" bold=\"0\" italic=\"0\"
          smoothing=\"0\">
      <color use=\"0\" r=\"0\" g=\"0\" b=\"0\"/>
    </font>
    <window title=\"Console\" icon=\"\" use_tab_icon=\"1\"
            use_console_title=\"0\" show_cmd=\"1\" show_cmd_tabs=\"1\"
            use_tab_title=\"1\" trim_tab_titles=\"20\"
            trim_tab_titles_right=\"0\"/>
    <controls show_menu=\"1\" show_toolbar=\"1\" show_statusbar=\"1\"
              show_tabs=\"1\" hide_single_tab=\"1\" show_scrollbars=\"1\"
              flat_scrollbars=\"0\" tabs_on_bottom=\"0\"/>
    <styles caption=\"1\" resizable=\"1\" taskbar_button=\"1\" border=\"1\"
            inside_border=\"2\" tray_icon=\"1\">
      <selection_color r=\"255\" g=\"255\" b=\"255\"/>
    </styles>
    <position x=\"-1\" y=\"-1\" dock=\"-1\" snap=\"0\" z_order=\"0\"
              save_position=\"0\"/>
    <transparency type=\"0\" active_alpha=\"255\" inactive_alpha=\"255\"
                  r=\"0\" g=\"0\" b=\"0\"/>
  </appearance>
  <behavior>
    <copy_paste copy_on_select=\"0\" clear_on_copy=\"1\" no_wrap=\"1\"
                trim_spaces=\"1\" copy_newline_char=\"0\"
                sensitive_copy=\"1\"/>
    <scroll page_scroll_rows=\"0\"/>
    <tab_highlight flashes=\"3\" stay_highligted=\"1\"/>
  </behavior>
  <hotkeys use_scroll_lock=\"1\">
    <hotkey ctrl=\"1\" shift=\"0\" alt=\"0\" extended=\"0\" code=\"83\" command=\"settings\"/>
    <hotkey ctrl=\"0\" shift=\"0\" alt=\"0\" extended=\"0\" code=\"112\" command=\"help\"/>
    <hotkey ctrl=\"0\" shift=\"0\" alt=\"1\" extended=\"0\" code=\"115\" command=\"exit\"/>
    <hotkey ctrl=\"1\" shift=\"0\" alt=\"0\" extended=\"0\" code=\"112\" command=\"newtab1\"/>
    <hotkey ctrl=\"1\" shift=\"0\" alt=\"0\" extended=\"0\" code=\"113\" command=\"newtab2\"/>
    <hotkey ctrl=\"1\" shift=\"0\" alt=\"0\" extended=\"0\" code=\"114\" command=\"newtab3\"/>
    <hotkey ctrl=\"1\" shift=\"0\" alt=\"0\" extended=\"0\" code=\"115\" command=\"newtab4\"/>
    <hotkey ctrl=\"1\" shift=\"0\" alt=\"0\" extended=\"0\" code=\"116\" command=\"newtab5\"/>
    <hotkey ctrl=\"1\" shift=\"0\" alt=\"0\" extended=\"0\" code=\"117\" command=\"newtab6\"/>
    <hotkey ctrl=\"1\" shift=\"0\" alt=\"0\" extended=\"0\" code=\"118\" command=\"newtab7\"/>
    <hotkey ctrl=\"1\" shift=\"0\" alt=\"0\" extended=\"0\" code=\"119\" command=\"newtab8\"/>
    <hotkey ctrl=\"1\" shift=\"0\" alt=\"0\" extended=\"0\" code=\"120\" command=\"newtab9\"/>
    <hotkey ctrl=\"1\" shift=\"0\" alt=\"0\" extended=\"0\" code=\"121\" command=\"newtab10\"/>
    <hotkey ctrl=\"1\" shift=\"0\" alt=\"0\" extended=\"0\" code=\"49\" command=\"switchtab1\"/>
    <hotkey ctrl=\"1\" shift=\"0\" alt=\"0\" extended=\"0\" code=\"50\" command=\"switchtab2\"/>
    <hotkey ctrl=\"1\" shift=\"0\" alt=\"0\" extended=\"0\" code=\"51\" command=\"switchtab3\"/>
    <hotkey ctrl=\"1\" shift=\"0\" alt=\"0\" extended=\"0\" code=\"52\" command=\"switchtab4\"/>
    <hotkey ctrl=\"1\" shift=\"0\" alt=\"0\" extended=\"0\" code=\"53\" command=\"switchtab5\"/>
    <hotkey ctrl=\"1\" shift=\"0\" alt=\"0\" extended=\"0\" code=\"54\" command=\"switchtab6\"/>
    <hotkey ctrl=\"1\" shift=\"0\" alt=\"0\" extended=\"0\" code=\"55\" command=\"switchtab7\"/>
    <hotkey ctrl=\"1\" shift=\"0\" alt=\"0\" extended=\"0\" code=\"56\" command=\"switchtab8\"/>
    <hotkey ctrl=\"1\" shift=\"0\" alt=\"0\" extended=\"0\" code=\"57\" command=\"switchtab9\"/>
    <hotkey ctrl=\"1\" shift=\"0\" alt=\"0\" extended=\"0\" code=\"48\" command=\"switchtab10\"/>
    <hotkey ctrl=\"1\" shift=\"0\" alt=\"0\" extended=\"0\" code=\"9\" command=\"nexttab\"/>
    <hotkey ctrl=\"1\" shift=\"1\" alt=\"0\" extended=\"0\" code=\"9\" command=\"prevtab\"/>
    <hotkey ctrl=\"1\" shift=\"0\" alt=\"0\" extended=\"0\" code=\"87\" command=\"closetab\"/>
    <hotkey ctrl=\"1\" shift=\"0\" alt=\"0\" extended=\"0\" code=\"82\" command=\"renametab\"/>
    <hotkey ctrl=\"1\" shift=\"0\" alt=\"0\" extended=\"1\" code=\"45\" command=\"copy\"/>
    <hotkey ctrl=\"1\" shift=\"0\" alt=\"0\" extended=\"1\" code=\"46\" command=\"clear_selection\"/>
    <hotkey ctrl=\"0\" shift=\"1\" alt=\"0\" extended=\"1\" code=\"45\" command=\"paste\"/>
    <hotkey ctrl=\"0\" shift=\"0\" alt=\"0\" extended=\"0\" code=\"0\" command=\"stopscroll\"/>
    <hotkey ctrl=\"0\" shift=\"0\" alt=\"0\" extended=\"0\" code=\"0\" command=\"scrollrowup\"/>
    <hotkey ctrl=\"0\" shift=\"0\" alt=\"0\" extended=\"0\" code=\"0\" command=\"scrollrowdown\"/>
    <hotkey ctrl=\"0\" shift=\"0\" alt=\"0\" extended=\"0\" code=\"0\" command=\"scrollpageup\"/>
    <hotkey ctrl=\"0\" shift=\"0\" alt=\"0\" extended=\"0\" code=\"0\" command=\"scrollpagedown\"/>
    <hotkey ctrl=\"0\" shift=\"0\" alt=\"0\" extended=\"0\" code=\"0\" command=\"scrollcolleft\"/>
    <hotkey ctrl=\"0\" shift=\"0\" alt=\"0\" extended=\"0\" code=\"0\" command=\"scrollcolright\"/>
    <hotkey ctrl=\"0\" shift=\"0\" alt=\"0\" extended=\"0\" code=\"0\" command=\"scrollpageleft\"/>
    <hotkey ctrl=\"0\" shift=\"0\" alt=\"0\" extended=\"0\" code=\"0\" command=\"scrollpageright\"/>
    <hotkey ctrl=\"1\" shift=\"1\" alt=\"0\" extended=\"0\" code=\"112\" command=\"dumpbuffer\"/>
    <hotkey ctrl=\"0\" shift=\"0\" alt=\"0\" extended=\"0\" code=\"0\" command=\"activate\"/>
  </hotkeys>
  <mouse>
    <actions>
      <action ctrl=\"0\" shift=\"0\" alt=\"0\" button=\"1\" name=\"copy\"/>
      <action ctrl=\"0\" shift=\"1\" alt=\"0\" button=\"1\" name=\"select\"/>
      <action ctrl=\"0\" shift=\"0\" alt=\"0\" button=\"3\" name=\"paste\"/>
      <action ctrl=\"1\" shift=\"0\" alt=\"0\" button=\"1\" name=\"drag\"/>
      <action ctrl=\"0\" shift=\"0\" alt=\"0\" button=\"2\" name=\"menu\"/>
    </actions>
  </mouse>
  <tabs>\n")
  if(X64)
    set(CONSOLE_XML "${CONSOLE_XML}\n
    <tab title=\"CMD MinGW 64-bit Environment\"
         icon=\"${CONSOLE2_NATIVE_DIR}\\cmd.ico\"
         use_default_icon=\"0\">
      <console shell=\"cmd /K ${BASEDIR_NATIVE}\\mingw64_env.cmd\"
               init_dir=\"${BASEDIR_NATIVE}\" run_as_user=\"0\" user=\"\"/>
      <cursor style=\"0\" r=\"255\" g=\"255\" b=\"255\"/>
      <background type=\"0\" r=\"0\" g=\"0\" b=\"0\">
	<image file=\"\" relative=\"0\" extend=\"0\" position=\"0\">
	  <tint opacity=\"0\" r=\"0\" g=\"0\" b=\"0\"/>
	</image>
      </background>
    </tab>
    <tab title=\"MSYS MinGW 64-bit Environment\"
         icon=\"${CONSOLE2_NATIVE_DIR}\\msys.ico\"
         use_default_icon=\"0\">
      <console shell=\"cmd /K ${BASEDIR_NATIVE}\\msys64_env.cmd\"
               init_dir=\"${BASEDIR_NATIVE}\" run_as_user=\"0\" user=\"\"/>
      <cursor style=\"0\" r=\"255\" g=\"255\" b=\"255\"/>
      <background type=\"0\" r=\"0\" g=\"0\" b=\"0\">
	<image file=\"\" relative=\"0\" extend=\"0\" position=\"0\">
	  <tint opacity=\"0\" r=\"0\" g=\"0\" b=\"0\"/>
	</image>
      </background>
    </tab>\n")
  endif()  
  set(CONSOLE_XML "${CONSOLE_XML}\n
    <tab title=\"CMD MinGW 32-bit Environment\"
         icon=\"${CONSOLE2_NATIVE_DIR}\\cmd.ico\"
         use_default_icon=\"0\">
      <console shell=\"cmd /K ${BASEDIR_NATIVE}\\mingw32_env.cmd\"
               init_dir=\"${BASEDIR_NATIVE}\" run_as_user=\"0\" user=\"\"/>
      <cursor style=\"0\" r=\"255\" g=\"255\" b=\"255\"/>
      <background type=\"0\" r=\"0\" g=\"0\" b=\"0\">
	<image file=\"\" relative=\"0\" extend=\"0\" position=\"0\">
	  <tint opacity=\"0\" r=\"0\" g=\"0\" b=\"0\"/>
	</image>
      </background>
    </tab>
    <tab title=\"MSYS MinGW 32-bit Environment\"
         icon=\"${CONSOLE2_NATIVE_DIR}\\msys.ico\"
         use_default_icon=\"0\">
      <console shell=\"cmd /K ${BASEDIR_NATIVE}\\msys32_env.cmd\"
               init_dir=\"${BASEDIR_NATIVE}\" run_as_user=\"0\" user=\"\"/>
      <cursor style=\"0\" r=\"255\" g=\"255\" b=\"255\"/>
      <background type=\"0\" r=\"0\" g=\"0\" b=\"0\">
	<image file=\"\" relative=\"0\" extend=\"0\" position=\"0\">
	  <tint opacity=\"0\" r=\"0\" g=\"0\" b=\"0\"/>
	</image>
      </background>
    </tab>
    <tab title=\"Console2\" use_default_icon=\"1\">
      <console shell=\"\" init_dir=\"\" run_as_user=\"0\" user=\"\"/>
      <cursor style=\"0\" r=\"255\" g=\"255\" b=\"255\"/>
      <background type=\"0\" r=\"0\" g=\"0\" b=\"0\">
	<image file=\"\" relative=\"0\" extend=\"0\" position=\"0\">
	  <tint opacity=\"0\" r=\"0\" g=\"0\" b=\"0\"/>
	</image>
      </background>
    </tab>
  </tabs>
</settings>")
  set(CONSOLE_XML_FILE "$ENV{APPDATA}/Console/console.xml")
  if(EXISTS "${CONSOLE_XML_FILE}")
    file(MD5 "${CONSOLE_XML_FILE}" CONSOLE_MD5) 
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy
      "${CONSOLE_XML_FILE}" "${CONSOLE_XML_FILE}_${CONSOLE_MD5}"
      WORKING_DIRECTORY ${TMPDIR}
      )
  endif()

  file(WRITE "${CONSOLE_XML_FILE}" "${CONSOLE_XML}")

  set(CONSOLE_LNK
    "
L_Welcome_MsgBox_Message_Text = _
    \"A shortcut to Console2 will be created on your desktop.\" & vbcrlf & _
    vbcrlf & \"The MinGW environments are available inside Console2 \" & _
    vbcrlf & \"via File -> New Tab -> MinGW [32|64]-bit Environment.\"
L_Welcome_MsgBox_Title_Text = _
    \"Creation of shortcut to Console2\"
Call Welcome()
Dim WSHShell
Set WSHShell = _
   WScript.CreateObject(\"WScript.Shell\")
Dim MyShortcut, MyDesktop, DesktopPath
' Read desktop path using WshSpecialFolders object
DesktopPath = _
   WSHShell.SpecialFolders(\"Desktop\")
' Create a shortcut object on the desktop
Set MyShortcut = _
   WSHShell.CreateShortcut( _
   DesktopPath & \"\\Console2.lnk\")
' Set shortcut object properties and save it
MyShortcut.TargetPath = _
   WSHShell.ExpandEnvironmentStrings( _
   \"${BASEDIR_NATIVE}\\Console2\\Console.exe\")
MyShortcut.WorkingDirectory = _
   WSHShell.ExpandEnvironmentStrings( _
   \"${BASEDIR_NATIVE}\")
MyShortcut.WindowStyle = 4
MyShortcut.IconLocation = _
   WSHShell.ExpandEnvironmentStrings( _
   \"${BASEDIR_NATIVE}\\Console2\\Console.exe, 0\")
MyShortcut.Save
WScript.Echo _
   \"A shortcut to Console2 now exists on your Desktop.\"
Sub Welcome()
   Dim intDoIt
   intDoIt = MsgBox(L_Welcome_MsgBox_Message_Text, _
      vbOKOnly + vbInformation, _
      L_Welcome_MsgBox_Title_Text )
   If intDoIt = vbCancel Then
      WScript.Quit
   End If
End Sub")
  file(WRITE "${TMPDIR}/link.vbs" "${CONSOLE_LNK}")
  execute_process(
    COMMAND cscript /Nologo link.vbs
    WORKING_DIRECTORY ${TMPDIR}
    )

  file(REMOVE_RECURSE ${TMPDIR})
  file(MAKE_DIRECTORY ${TMPDIR})
endif()

#=============================================================================
# Create a simple CMake project for testing the environment.
#=============================================================================
SET(EXAMPLE_CMAKE_PROJ_DIR "${BASEDIR}/example_cmake_project")
if(NOT EXISTS "${EXAMPLE_CMAKE_PROJ_DIR}")
  file(MAKE_DIRECTORY "${EXAMPLE_CMAKE_PROJ_DIR}")

  file(WRITE "${EXAMPLE_CMAKE_PROJ_DIR}/CMakeLists.txt"
    "PROJECT(example_project C CXX Fortran)

CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

#=============================================================================
# C, C++ and Fortran hello world programs.
#=============================================================================
ADD_EXECUTABLE(example_c example.c)

ADD_EXECUTABLE(example_cpp example.cpp)

ADD_EXECUTABLE(example_fortran example.f90)

#=============================================================================
# Check for version control systems.
#=============================================================================
FIND_PACKAGE(Git)
FIND_PACKAGE(Subversion)

#=============================================================================
# A simple Qt4 program.
#=============================================================================
# https://qt-project.org/quarterly/view/using_cmake_to_build_qt_projects
FIND_PACKAGE(Qt4 REQUIRED)

include(\${QT_USE_FILE})
add_definitions(\${QT_DEFINITIONS})
include_directories(\${QT_INCLUDE_DIR})
set(LIBRARIES \${LIBRARIES} \${QT_LIBRARIES})

# Local files
set(HEADERS MainWindow.hpp)
set(SOURCES MainWindow.cpp Main.cpp)

# Create Qt MOC files
qt4_wrap_cpp(HEADERS_MOC \${HEADERS})

# Create executable
add_executable(qtapp \${HEADERS_MOC} \${SOURCES})
target_link_libraries(qtapp \${LIBRARIES})
set_target_properties(qtapp PROPERTIES OUTPUT_NAME qtapp)

#=============================================================================
# An extended Python interpreter.
#=============================================================================

# Since FindPython.cmake does not respect PYTHONHOME on Windows we have to
# specify the path to the headers and libs by hand.
file(TO_CMAKE_PATH \"\$ENV{PYTHONHOME}\" PYTHONHOME)

set(PYTHON_INCLUDE_DIR \"\${PYTHONHOME}/include/python2.7\" CACHE PATH
  \"Path to where Python.h is found\")
set(PYTHON_LIBRARY \"\${PYTHONHOME}/lib/python2.7/config/libpython2.7.a\" CACHE PATH
  \"Path to Python library\")

FIND_PACKAGE(PythonLibs)
INCLUDE_DIRECTORIES(\${PYTHON_INCLUDE_PATH})

ADD_EXECUTABLE(mypython pymodule.c)
target_link_libraries (mypython \${PYTHON_LIBRARIES})
")


  file(WRITE "${EXAMPLE_CMAKE_PROJ_DIR}/Main.cpp"
"#include \"MainWindow.hpp\"
#include <QApplication>

int main(int argc, char** argv)
{
    QApplication app(argc, argv);
    MainWindow window;
    window.show();
    return app.exec();
}\n")

  file(WRITE "${EXAMPLE_CMAKE_PROJ_DIR}/MainWindow.hpp"
"#pragma once
#include <QMainWindow>

class MainWindow : public QMainWindow
{
    Q_OBJECT
public:
    MainWindow(void);
};\n")

  file(WRITE "${EXAMPLE_CMAKE_PROJ_DIR}/MainWindow.cpp"
"#include \"MainWindow.hpp\"
#include <QLabel>

MainWindow::MainWindow(void)
{
    QLabel* label = new QLabel(\"Success!\");
    setCentralWidget(label);
}\n")

  file(WRITE "${EXAMPLE_CMAKE_PROJ_DIR}/pymodule.c"
"#include <Python.h>                         

static PyObject *
spam_simon(PyObject *self, PyObject *args)
{
    const char *command;
    int sts = 42;

    if (!PyArg_ParseTuple(args, \"s\", &command))
        return NULL;

    printf(\"Jipijajei Schweinebacke!\\nArgument: %s\\n\", command);

    return Py_BuildValue(\"i\", sts);
}

static PyMethodDef SpamMethods[] = {
    {\"simon\",  spam_simon, METH_VARARGS,
     \"Call Simon's method.\"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC
initspam(void)
{
    (void) Py_InitModule(\"spam\", SpamMethods);
}

int
main(int argc, char *argv[])
{
char *args[] = {\"embed\", \"hello\", 0};

    /* Pass argv[0] to the Python interpreter */
    Py_SetProgramName(argv[0]);

    /* Initialize the Python interpreter.  Required. */
    Py_Initialize();

    /* Add a static module */
    initspam();



PySys_SetArgvEx(2, args, 0);

printf(\"Hello, brave new world\\n\\n\");

PyRun_SimpleString(\"import sys\\n\");
PyRun_SimpleString(\"print(sys.builtin_module_names)\\n\");
PyRun_SimpleString(\"print(sys.modules.keys())\\n\");
PyRun_SimpleString(\"print(sys.executable)\\n\");
PyRun_SimpleString(\"print(sys.argv)\\n\");

PyRun_SimpleString(\"import spam\\n\");
PyRun_SimpleString(\"print(spam.simon(\\\"FRANZ-XAVER\\\"))\\n\");

printf(\"\\nGoodbye, cruel world\\n\");

Py_Exit(0);
}\n")

  file(WRITE "${EXAMPLE_CMAKE_PROJ_DIR}/example.c"
    "#include <stdio.h>

int main(int argc, char** argv) {
  printf(\"Hello world!\\n\");
  return 0;
}
")

  file(WRITE "${EXAMPLE_CMAKE_PROJ_DIR}/example.cpp"
    "#include <iostream>

int main(int argc, char** argv) {
  std::cout << \"Hello world!\" << std::endl;
  return 0;
}
")

  file(WRITE "${EXAMPLE_CMAKE_PROJ_DIR}/example.f90"
    "program hello
   print *, \"Hello World!\"
end program hello
")
endif()


write_mingw_env()

