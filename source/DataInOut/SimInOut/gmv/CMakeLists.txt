SET(SIMINOUTGMV_SRCS SimInputGMV.cc gmvread.cc SimOutGMV.cc)
ADD_LIBRARY(siminoutgmv STATIC ${SIMINOUTGMV_SRCS})

ADD_DEPENDENCIES(siminoutgmv boost)
