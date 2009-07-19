ADD_CUSTOM_TARGET(share_xml ALL
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CFS_SOURCE_DIR}/share/xml/CFS-Simulation/CFS.xsd ${CFS_BINARY_DIR}/share/xml/CFS-Simulation/CFS.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CFS_SOURCE_DIR}/share/xml/CFS-Simulation/Schemas/CFS_Analysis.xsd ${CFS_BINARY_DIR}/share/xml/CFS-Simulation/Schemas/CFS_Analysis.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CFS_SOURCE_DIR}/share/xml/CFS-Simulation/Schemas/CFS_Coupling.xsd ${CFS_BINARY_DIR}/share/xml/CFS-Simulation/Schemas/CFS_Coupling.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CFS_SOURCE_DIR}/share/xml/CFS-Simulation/Schemas/CFS_DirectCoupling.xsd ${CFS_BINARY_DIR}/share/xml/CFS-Simulation/Schemas/CFS_DirectCoupling.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CFS_SOURCE_DIR}/share/xml/CFS-Simulation/Schemas/CFS_Domain.xsd ${CFS_BINARY_DIR}/share/xml/CFS-Simulation/Schemas/CFS_Domain.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CFS_SOURCE_DIR}/share/xml/CFS-Simulation/Schemas/CFS_ErsatzMaterial.xsd ${CFS_BINARY_DIR}/share/xml/CFS-Simulation/Schemas/CFS_ErsatzMaterial.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CFS_SOURCE_DIR}/share/xml/CFS-Simulation/Schemas/CFS_FileFormats.xsd ${CFS_BINARY_DIR}/share/xml/CFS-Simulation/Schemas/CFS_FileFormats.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CFS_SOURCE_DIR}/share/xml/CFS-Simulation/Schemas/CFS_IntegRules.xsd ${CFS_BINARY_DIR}/share/xml/CFS-Simulation/Schemas/CFS_IntegRules.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CFS_SOURCE_DIR}/share/xml/CFS-Simulation/Schemas/CFS_LinearSystem.xsd ${CFS_BINARY_DIR}/share/xml/CFS-Simulation/Schemas/CFS_LinearSystem.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CFS_SOURCE_DIR}/share/xml/CFS-Simulation/Schemas/CFS_MeshDeformations.xsd ${CFS_BINARY_DIR}/share/xml/CFS-Simulation/Schemas/CFS_MeshDeformations.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CFS_SOURCE_DIR}/share/xml/CFS-Simulation/Schemas/CFS_Misc.xsd ${CFS_BINARY_DIR}/share/xml/CFS-Simulation/Schemas/CFS_Misc.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CFS_SOURCE_DIR}/share/xml/CFS-Simulation/Schemas/CFS_MpCCI.xsd ${CFS_BINARY_DIR}/share/xml/CFS-Simulation/Schemas/CFS_MpCCI.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CFS_SOURCE_DIR}/share/xml/CFS-Simulation/Schemas/CFS_Optimization.xsd ${CFS_BINARY_DIR}/share/xml/CFS-Simulation/Schemas/CFS_Optimization.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CFS_SOURCE_DIR}/share/xml/CFS-Simulation/Schemas/CFS_PDEacoustic.xsd ${CFS_BINARY_DIR}/share/xml/CFS-Simulation/Schemas/CFS_PDEacoustic.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CFS_SOURCE_DIR}/share/xml/CFS-Simulation/Schemas/CFS_PDEacousticMixed.xsd ${CFS_BINARY_DIR}/share/xml/CFS-Simulation/Schemas/CFS_PDEacousticMixed.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CFS_SOURCE_DIR}/share/xml/CFS-Simulation/Schemas/CFS_PDEbasic.xsd ${CFS_BINARY_DIR}/share/xml/CFS-Simulation/Schemas/CFS_PDEbasic.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CFS_SOURCE_DIR}/share/xml/CFS-Simulation/Schemas/CFS_PDEelectrostatic.xsd ${CFS_BINARY_DIR}/share/xml/CFS-Simulation/Schemas/CFS_PDEelectrostatic.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CFS_SOURCE_DIR}/share/xml/CFS-Simulation/Schemas/CFS_PDEfluidMech.xsd ${CFS_BINARY_DIR}/share/xml/CFS-Simulation/Schemas/CFS_PDEfluidMech.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CFS_SOURCE_DIR}/share/xml/CFS-Simulation/Schemas/CFS_PDEheatConduction.xsd ${CFS_BINARY_DIR}/share/xml/CFS-Simulation/Schemas/CFS_PDEheatConduction.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CFS_SOURCE_DIR}/share/xml/CFS-Simulation/Schemas/CFS_PDEmagedge.xsd ${CFS_BINARY_DIR}/share/xml/CFS-Simulation/Schemas/CFS_PDEmagedge.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CFS_SOURCE_DIR}/share/xml/CFS-Simulation/Schemas/CFS_PDEmagnetic.xsd ${CFS_BINARY_DIR}/share/xml/CFS-Simulation/Schemas/CFS_PDEmagnetic.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CFS_SOURCE_DIR}/share/xml/CFS-Simulation/Schemas/CFS_PDEmechanic.xsd ${CFS_BINARY_DIR}/share/xml/CFS-Simulation/Schemas/CFS_PDEmechanic.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CFS_SOURCE_DIR}/share/xml/CFS-Simulation/Schemas/CFS_PDEmpcci.xsd ${CFS_BINARY_DIR}/share/xml/CFS-Simulation/Schemas/CFS_PDEmpcci.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CFS_SOURCE_DIR}/share/xml/CFS-Simulation/Schemas/CFS_PDEsmooth.xsd ${CFS_BINARY_DIR}/share/xml/CFS-Simulation/Schemas/CFS_PDEsmooth.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CFS_SOURCE_DIR}/share/xml/CFS-Simulation/Schemas/CFS_PDEstokesFluid.xsd ${CFS_BINARY_DIR}/share/xml/CFS-Simulation/Schemas/CFS_PDEstokesFluid.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CFS_SOURCE_DIR}/share/xml/CFS-Simulation/Schemas/CFS_ParamIdent.xsd ${CFS_BINARY_DIR}/share/xml/CFS-Simulation/Schemas/CFS_ParamIdent.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CFS_SOURCE_DIR}/share/xml/CFS-Simulation/Schemas/CFS_PostProc.xsd ${CFS_BINARY_DIR}/share/xml/CFS-Simulation/Schemas/CFS_PostProc.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CFS_SOURCE_DIR}/share/xml/CFS-Simulation/Schemas/CFS_Preconds.xsd ${CFS_BINARY_DIR}/share/xml/CFS-Simulation/Schemas/CFS_Preconds.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CFS_SOURCE_DIR}/share/xml/CFS-Simulation/Schemas/CFS_Solvers.xsd ${CFS_BINARY_DIR}/share/xml/CFS-Simulation/Schemas/CFS_Solvers.xsd

  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CFS_SOURCE_DIR}/share/xml/CFS-Material/CFS_Material.xsd ${CFS_BINARY_DIR}/share/xml/CFS-Material/CFS_Material.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CFS_SOURCE_DIR}/share/xml/CFS-Material/Schemas/CFS_MATacoustic.xsd ${CFS_BINARY_DIR}/share/xml/CFS-Material/Schemas/CFS_MATacoustic.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CFS_SOURCE_DIR}/share/xml/CFS-Material/Schemas/CFS_MATelectric.xsd ${CFS_BINARY_DIR}/share/xml/CFS-Material/Schemas/CFS_MATelectric.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CFS_SOURCE_DIR}/share/xml/CFS-Material/Schemas/CFS_MATflow.xsd ${CFS_BINARY_DIR}/share/xml/CFS-Material/Schemas/CFS_MATflow.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CFS_SOURCE_DIR}/share/xml/CFS-Material/Schemas/CFS_MATheatConduction.xsd ${CFS_BINARY_DIR}/share/xml/CFS-Material/Schemas/CFS_MATheatConduction.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CFS_SOURCE_DIR}/share/xml/CFS-Material/Schemas/CFS_MATmagnetic.xsd ${CFS_BINARY_DIR}/share/xml/CFS-Material/Schemas/CFS_MATmagnetic.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CFS_SOURCE_DIR}/share/xml/CFS-Material/Schemas/CFS_MATmechanical.xsd ${CFS_BINARY_DIR}/share/xml/CFS-Material/Schemas/CFS_MATmechanical.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CFS_SOURCE_DIR}/share/xml/CFS-Material/Schemas/CFS_MATpiezo.xsd ${CFS_BINARY_DIR}/share/xml/CFS-Material/Schemas/CFS_MATpiezo.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CFS_SOURCE_DIR}/share/xml/CFS-Material/Schemas/CFS_MATpyroelectric.xsd ${CFS_BINARY_DIR}/share/xml/CFS-Material/Schemas/CFS_MATpyroelectric.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CFS_SOURCE_DIR}/share/xml/CFS-Material/Schemas/CFS_MATthermoelastic.xsd ${CFS_BINARY_DIR}/share/xml/CFS-Material/Schemas/CFS_MATthermoelastic.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CFS_SOURCE_DIR}/share/xml/CFS-Material/Schemas/CFS_Misc.xsd ${CFS_BINARY_DIR}/share/xml/CFS-Material/Schemas/CFS_Misc.xsd

  COMMENT "Updating contents of ${CFS_BINARY_DIR}/share/xml..."
)
