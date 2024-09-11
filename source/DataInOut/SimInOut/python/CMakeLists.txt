if(USE_EMBEDDED_PYTHON)
  # see FindPrograms.cmake on how the Python_* stuff is set
  list(APPEND TARGET_LL ${Python_LIBRARIES})
  
  set(SIMINPUTPYTHON_SRCS SimInputPython.cc)

  add_library(siminputpython STATIC ${SIMINPUTPYTHON_SRCS})
  
  target_link_libraries(siminputpython ${TARGET_LL})
endif()
