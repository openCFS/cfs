SET(ODEDESCR_SRCS
    Gilmore.cc
    Gilmoredimlos.cc
    KellerMiksis.cc
    LinearKellerMiksis.cc
    RPNNP.cc
   )

ADD_LIBRARY(odedescr STATIC ${ODEDESCR_SRCS})
