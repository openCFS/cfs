if(USE_EMBEDDED_PYTHON)
  # see FindPrograms.cmake on how the PYTHON_ stuff is set
  set(TARGET_LL ${TARGET_LL} ${PYTHON_LIBRARY})
  # already in FindCFSDEPS.cmake: include_directories(${PYTHON_INCLUDE_DIR})
  #include_directories(${PYTHON_SITE_PACKAGES_DIR})
  
  set(SIMINPUTPYTHON_SRCS SimInputPython.cc)

  add_library(siminputpython STATIC ${SIMINPUTPYTHON_SRCS})
  # required for numpy
  include_directories(${PYTHON_SITE_PACKAGES_DIR})
  
  target_link_libraries(siminputpython ${TARGET_LL})
endif()
