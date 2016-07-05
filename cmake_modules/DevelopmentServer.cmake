#=============================================================================
# Global defines for CFS++ development server.
#=============================================================================

#-----------------------------------------------------------------------------
# CFS_(D)evelopment_(S)erver_HOSTNAME
#-----------------------------------------------------------------------------
SET(CFS_DS_HOSTNAME "cfs.mdmt.tuwien.ac.at")

#-----------------------------------------------------------------------------
# Base directory for public FTP server.
#-----------------------------------------------------------------------------
SET(CFS_DS_FTP "ftp://${CFS_DS_HOSTNAME}")

#-----------------------------------------------------------------------------
# Base directory for secure FTP server.
#-----------------------------------------------------------------------------
SET(CFS_DS_SFTP "sftp://${CFS_DS_HOSTNAME}")

#-----------------------------------------------------------------------------
# Base directory for restricted HTTPS server.
#-----------------------------------------------------------------------------
SET(CFS_DS_HTTPS "https://${CFS_DS_HOSTNAME}")

#-----------------------------------------------------------------------------
# Base directory of WEBDAV server.
#-----------------------------------------------------------------------------
SET(CFS_DS_WEBDAV "${CFS_DS_HTTPS}/files")

# ----------------------------------------------------------------------------
# Alternative mirror from the optimization group - plain http!
# ----------------------------------------------------------------------------
SET(CFS_FAU_MIRROR "http://eamc080.eam.uni-erlangen.de/cfsdeps")

#-----------------------------------------------------------------------------
# Base directory of Subversion server.
#-----------------------------------------------------------------------------
SET(CFS_DS_SVN "${CFS_DS_HTTPS}/svn")

#-----------------------------------------------------------------------------
# Base directory CFS_DEPS_CACHE_DIR mirror.
#-----------------------------------------------------------------------------
SET(CFS_DS_CFSDEPS "${CFS_DS_FTP}")

#-----------------------------------------------------------------------------
# Drop site for nightly build/test results.
# The port is within TU-Wien 44088 but from outside (FAU) 1022 and the port needs
# to be opened for the sending IP.
# Hence FAU needs to configure with CFS_DS_PORT
#-----------------------------------------------------------------------------
SET(CFS_DS_PORT_DEFAULT "44088")
IF(NOT CFS_DS_PORT)
  SET(CFS_DS_PORT ${CFS_DS_PORT_DEFAULT})
ENDIF()
SET(CFS_DS_CDASH_DROP_SITE "${CFS_DS_HOSTNAME}:${CFS_DS_PORT}")
