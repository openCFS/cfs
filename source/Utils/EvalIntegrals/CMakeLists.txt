SET(EVALINTEGRAL_SRC
  BiotSavart.cc
  )

ADD_LIBRARY(evalintegral STATIC ${EVALINTEGRAL_SRC})

SET(TARGET_LL
  mathparser
  datainout
  )


TARGET_LINK_LIBRARIES(evalintegral
   ${TARGET_LL}
  )


