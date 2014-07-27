# Auto    generate   bootstrap_devel_machine.txt   from    distro.sh   and
# bootstrap_devel_machine.sh

FILE(READ "distro.sh" DISTRO_SH)
STRING(REPLACE ";" "__XXXXXXXXX__" DISTRO_SH "${DISTRO_SH}")
STRING(REPLACE "\n" ";" DISTRO_SH ${DISTRO_SH})

SET(printout on)
SET(OUTPUT "")

foreach(line IN LISTS DISTRO_SH)
  if(line MATCHES "^Help")
    SET(printout off)
  endif()
  if(line MATCHES "^# Set default sub-architecture")
    SET(printout on)
  endif()

  if(line MATCHES "while :")
    SET(printout off)
  endif()
  
  if(printout)
    foreach(part IN LISTS line)
      LIST(APPEND OUTPUT "${part}")
    endforeach(part)
    
#    MESSAGE("${line}")
  endif(printout)
  
endforeach(line)


FILE(READ "bootstrap_devel_machine.sh" BOOTSTRAP_SH)
STRING(REPLACE ";" "__XXXXXXXXX__" BOOTSTRAP_SH "${BOOTSTRAP_SH}")
STRING(REPLACE "\n" ";" BOOTSTRAP_SH ${BOOTSTRAP_SH})

foreach(line IN LISTS BOOTSTRAP_SH)
  if(line MATCHES "^# This script")
    SET(printout on)
  endif()

  if(printout)
    foreach(part IN LISTS line)
      LIST(APPEND OUTPUT "${part}")
    endforeach(part)
    
#    MESSAGE("${line}")
  endif(printout)
  
endforeach(line)

STRING(REPLACE ";" "\n" OUTPUT "${OUTPUT}")
STRING(REPLACE "__XXXXXXXXX__" ";" OUTPUT ${OUTPUT})

# FILE(WRITE "out.txt" "${OUTPUT}")
MESSAGE("${OUTPUT}")

