SET(SRC_DIR ${CFS_SOURCE_DIR}/share/xml/CFS-Simulation)
SET(DST_DIR ${CFS_BINARY_DIR}/share/xml/CFS-Simulation)

ADD_CUSTOM_TARGET(share_xml_cfssim ALL
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/CFS.xsd ${DST_DIR}/CFS.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_Analysis.xsd ${DST_DIR}/Schemas/CFS_Analysis.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_Coupling.xsd ${DST_DIR}/Schemas/CFS_Coupling.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_DirectCoupling.xsd ${DST_DIR}/Schemas/CFS_DirectCoupling.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_Domain.xsd ${DST_DIR}/Schemas/CFS_Domain.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_Documentation.xsd ${DST_DIR}/Schemas/CFS_Documentation.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_ErsatzMaterial.xsd ${DST_DIR}/Schemas/CFS_ErsatzMaterial.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_FileFormats.xsd ${DST_DIR}/Schemas/CFS_FileFormats.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_IntegRules.xsd ${DST_DIR}/Schemas/CFS_IntegRules.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_LinearSystem.xsd ${DST_DIR}/Schemas/CFS_LinearSystem.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_MeshDeformations.xsd ${DST_DIR}/Schemas/CFS_MeshDeformations.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_Misc.xsd ${DST_DIR}/Schemas/CFS_Misc.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_MpCCI.xsd ${DST_DIR}/Schemas/CFS_MpCCI.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_Optimization.xsd ${DST_DIR}/Schemas/CFS_Optimization.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDEacoustic.xsd ${DST_DIR}/Schemas/CFS_PDEacoustic.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDEacousticMixed.xsd ${DST_DIR}/Schemas/CFS_PDEacousticMixed.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDEbasic.xsd ${DST_DIR}/Schemas/CFS_PDEbasic.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDEelectrostatic.xsd ${DST_DIR}/Schemas/CFS_PDEelectrostatic.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDEfluidMech.xsd ${DST_DIR}/Schemas/CFS_PDEfluidMech.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDEheatConduction.xsd ${DST_DIR}/Schemas/CFS_PDEheatConduction.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDEmagedge.xsd ${DST_DIR}/Schemas/CFS_PDEmagedge.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDEmagnetic.xsd ${DST_DIR}/Schemas/CFS_PDEmagnetic.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDEmagneticScalar.xsd ${DST_DIR}/Schemas/CFS_PDEmagneticScalar.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDEmechanic.xsd ${DST_DIR}/Schemas/CFS_PDEmechanic.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDEmpcci.xsd ${DST_DIR}/Schemas/CFS_PDEmpcci.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDEsmooth.xsd ${DST_DIR}/Schemas/CFS_PDEsmooth.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDEstokesFluid.xsd ${DST_DIR}/Schemas/CFS_PDEstokesFluid.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_ParamIdent.xsd ${DST_DIR}/Schemas/CFS_ParamIdent.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PostProc.xsd ${DST_DIR}/Schemas/CFS_PostProc.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_Preconds.xsd ${DST_DIR}/Schemas/CFS_Preconds.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_Solvers.xsd ${DST_DIR}/Schemas/CFS_Solvers.xsd
  COMMENT "Updating XML schemata for CFS++ simulation..."
)

SET(SRC_DIR ${CFS_SOURCE_DIR}/share/xml/CFS-Material)
SET(DST_DIR ${CFS_BINARY_DIR}/share/xml/CFS-Material)

ADD_CUSTOM_TARGET(share_xml_cfsmat ALL
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/CFS_Material.xsd ${DST_DIR}/CFS_Material.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_MATacoustic.xsd ${DST_DIR}/Schemas/CFS_MATacoustic.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_MATelectric.xsd ${DST_DIR}/Schemas/CFS_MATelectric.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_MATflow.xsd ${DST_DIR}/Schemas/CFS_MATflow.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_MATheatConduction.xsd ${DST_DIR}/Schemas/CFS_MATheatConduction.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_MATmagnetic.xsd ${DST_DIR}/Schemas/CFS_MATmagnetic.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_MATmechanical.xsd ${DST_DIR}/Schemas/CFS_MATmechanical.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_MATpiezo.xsd ${DST_DIR}/Schemas/CFS_MATpiezo.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_MATpyroelectric.xsd ${DST_DIR}/Schemas/CFS_MATpyroelectric.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_MATthermoelastic.xsd ${DST_DIR}/Schemas/CFS_MATthermoelastic.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_Misc.xsd ${DST_DIR}/Schemas/CFS_Misc.xsd
  COMMENT "Updating XML schemata for CFS++ materials..."
)
