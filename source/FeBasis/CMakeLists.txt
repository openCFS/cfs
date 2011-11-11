SET(FEBASIS_SRC 
  BaseFE.cc
  FeSpace.cc
  FeFunctions.cc
  H1/FeSpaceH1.cc
  H1/FeSpaceH1Nodal.cc
  H1/FeSpaceH1Hi.cc
  HCurl/FeSpaceHCurlHi.cc
  H1/H1Elems.cc
  H1/H1ElemsLagExpl.cc
  H1/H1ElemsLagVar.cc
  H1/H1ElemsHi.cc
  HCurl/HCurlElemsHi.cc
  HCurl/HCurlElems.cc
)

ADD_LIBRARY(febasis STATIC ${FEBASIS_SRC})

TARGET_LINK_LIBRARIES(febasis
  paramh
  logging
  matvec
  cfsgeneral
  elemmapping
  domain
  results
  algsys-olas
  pde
  utils
  integration
)
