# ===========================================================================
#  Switch to mirror since official server is no more available for hardy.
# ===========================================================================
FILE(READ /etc/apt/sources.list SOURCES_LIST)
STRING(REPLACE
  "us.archive.ubuntu.com"
  "ftp.uni-muenster.de/pub/mirrors/ftp.ubuntu.com"
  SOURCES_LIST
  ${SOURCES_LIST})
STRING(REPLACE
  "security.ubuntu.com"
  "ftp.uni-muenster.de/pub/mirrors/ftp.ubuntu.com"
  SOURCES_LIST
  ${SOURCES_LIST})
FILE(WRITE sources.list ${SOURCES_LIST})

EXECUTE_PROCESS(
  COMMAND sudo mv sources.list /etc/apt/sources.list
  )

EXECUTE_PROCESS(
  COMMAND sudo apt-get update
  )

# ===========================================================================
#  Install latest tk-dev package needed for cgnsview
# ===========================================================================
EXECUTE_PROCESS(
  COMMAND sudo apt-get install -y tk-dev
  )
