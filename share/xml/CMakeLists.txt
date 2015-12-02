SET(SRC_DIR ${CFS_SOURCE_DIR}/share/xml/CFS-Simulation)
SET(DST_DIR ${CFS_BINARY_DIR}/share/xml/CFS-Simulation)

ADD_CUSTOM_TARGET(share_xml_cfssim ALL
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/CFS.xsd ${DST_DIR}/CFS.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_Analysis.xsd ${DST_DIR}/Schemas/CFS_Analysis.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_Coupling.xsd ${DST_DIR}/Schemas/CFS_Coupling.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_DirectCoupling.xsd ${DST_DIR}/Schemas/CFS_DirectCoupling.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_Domain.xsd ${DST_DIR}/Schemas/CFS_Domain.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_Documentation.xsd ${DST_DIR}/Schemas/CFS_Documentation.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_FePolynomial.xsd ${DST_DIR}/Schemas/CFS_FePolynomial.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_FileFormats.xsd ${DST_DIR}/Schemas/CFS_FileFormats.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_IntegrationScheme.xsd ${DST_DIR}/Schemas/CFS_IntegrationScheme.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_LinearSystem.xsd ${DST_DIR}/Schemas/CFS_LinearSystem.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_Misc.xsd ${DST_DIR}/Schemas/CFS_Misc.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDEacoustic.xsd ${DST_DIR}/Schemas/CFS_PDEacoustic.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDEacousticMixed.xsd ${DST_DIR}/Schemas/CFS_PDEacousticMixed.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDEbasic.xsd ${DST_DIR}/Schemas/CFS_PDEbasic.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDEelectrostatic.xsd ${DST_DIR}/Schemas/CFS_PDEelectrostatic.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDEfluidMech.xsd ${DST_DIR}/Schemas/CFS_PDEfluidMech.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDEheatConduction.xsd ${DST_DIR}/Schemas/CFS_PDEheatConduction.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDEmagedge.xsd ${DST_DIR}/Schemas/CFS_PDEmagedge.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDEmagnetic.xsd ${DST_DIR}/Schemas/CFS_PDEmagnetic.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDEmechanic.xsd ${DST_DIR}/Schemas/CFS_PDEmechanic.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDEelectricConduction.xsd ${DST_DIR}/Schemas/CFS_PDEelectricConduction.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDEwaterWave.xsd ${DST_DIR}/Schemas/CFS_PDEwaterWave.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDEtest.xsd ${DST_DIR}/Schemas/CFS_PDEtest.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDEsmooth.xsd ${DST_DIR}/Schemas/CFS_PDEsmooth.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PostProc.xsd ${DST_DIR}/Schemas/CFS_PostProc.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_Preconds.xsd ${DST_DIR}/Schemas/CFS_Preconds.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_Solvers.xsd ${DST_DIR}/Schemas/CFS_Solvers.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_Optimization.xsd ${DST_DIR}/Schemas/CFS_Optimization.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_ErsatzMaterial.xsd ${DST_DIR}/Schemas/CFS_ErsatzMaterial.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDELatticeBoltzmann.xsd ${DST_DIR}/Schemas/CFS_PDELatticeBoltzmann.xsd
  COMMENT "Updating XML schemata for CFS++ simulation..."
)
  
SET(SRC_DIR ${CFS_SOURCE_DIR}/share/xml/CFS-Material)
SET(DST_DIR ${CFS_BINARY_DIR}/share/xml/CFS-Material)

ADD_CUSTOM_TARGET(share_xml_cfsmat ALL
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/CFS_Material.xsd ${DST_DIR}/CFS_Material.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_MATacoustic.xsd ${DST_DIR}/Schemas/CFS_MATacoustic.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_MATelectric.xsd ${DST_DIR}/Schemas/CFS_MATelectric.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_MATflow.xsd ${DST_DIR}/Schemas/CFS_MATflow.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_MATelecConduction.xsd ${DST_DIR}/Schemas/CFS_MATelecConduction.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_MATtest.xsd ${DST_DIR}/Schemas/CFS_MATtest.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_MATheatConduction.xsd ${DST_DIR}/Schemas/CFS_MATheatConduction.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_MATmagnetic.xsd ${DST_DIR}/Schemas/CFS_MATmagnetic.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_MATmagStrict.xsd ${DST_DIR}/Schemas/CFS_MATmagStrict.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_MATmechanical.xsd ${DST_DIR}/Schemas/CFS_MATmechanical.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_MATpiezo.xsd ${DST_DIR}/Schemas/CFS_MATpiezo.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_MATpyroelectric.xsd ${DST_DIR}/Schemas/CFS_MATpyroelectric.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_MATthermoelastic.xsd ${DST_DIR}/Schemas/CFS_MATthermoelastic.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_Misc.xsd ${DST_DIR}/Schemas/CFS_Misc.xsd
  COMMENT "Updating XML schemata for CFS++ materials..."
)

SET(SRC_DIR ${CFS_SOURCE_DIR}/share/xml/config)
SET(DST_DIR ${CFS_BINARY_DIR}/share/xml/config)
ADD_CUSTOM_TARGET(share_xml_config ALL
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/cfstool.xml ${DST_DIR}/cfstool.xml
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/cfstool.json ${DST_DIR}/cfstool.json
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/log.xml ${DST_DIR}/log.xml
  COMMENT "Updating XML & JSON config files for cfstool and logging..."
)

SET(SRC_DIR ${CFS_SOURCE_DIR}/share/xml/CFS-Dat)
SET(DST_DIR ${CFS_BINARY_DIR}/share/xml/CFS-Dat)
ADD_CUSTOM_TARGET(share_xml_dat ALL
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/CFS_Dat.xsd ${DST_DIR}/CFS_Dat.xsd
  COMMENT "Updating XML schemata for CFS-Dat..."
)
