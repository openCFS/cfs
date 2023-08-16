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
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_HystereticSystem.xsd ${DST_DIR}/Schemas/CFS_HystereticSystem.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_LinearSystem.xsd ${DST_DIR}/Schemas/CFS_LinearSystem.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_Misc.xsd ${DST_DIR}/Schemas/CFS_Misc.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDEacoustic.xsd ${DST_DIR}/Schemas/CFS_PDEacoustic.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDEacousticMixed.xsd ${DST_DIR}/Schemas/CFS_PDEacousticMixed.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDEbasic.xsd ${DST_DIR}/Schemas/CFS_PDEbasic.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDEelectrostatic.xsd ${DST_DIR}/Schemas/CFS_PDEelectrostatic.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDEfluidMech.xsd ${DST_DIR}/Schemas/CFS_PDEfluidMech.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDEfluidMechLin.xsd ${DST_DIR}/Schemas/CFS_PDEfluidMechLin.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDEheatConduction.xsd ${DST_DIR}/Schemas/CFS_PDEheatConduction.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDEmagedge.xsd ${DST_DIR}/Schemas/CFS_PDEmagedge.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDEsplit.xsd ${DST_DIR}/Schemas/CFS_PDEsplit.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDEmagnetic.xsd ${DST_DIR}/Schemas/CFS_PDEmagnetic.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDEmechanic.xsd ${DST_DIR}/Schemas/CFS_PDEmechanic.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDEelectricConduction.xsd ${DST_DIR}/Schemas/CFS_PDEelectricConduction.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDEelecQuasistatic.xsd ${DST_DIR}/Schemas/CFS_PDEelecQuasistatic.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDEwaterWave.xsd ${DST_DIR}/Schemas/CFS_PDEwaterWave.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDEtest.xsd ${DST_DIR}/Schemas/CFS_PDEtest.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDEsmooth.xsd ${DST_DIR}/Schemas/CFS_PDEsmooth.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PostProc.xsd ${DST_DIR}/Schemas/CFS_PostProc.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_Preconds.xsd ${DST_DIR}/Schemas/CFS_Preconds.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_Solvers.xsd ${DST_DIR}/Schemas/CFS_Solvers.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_ErsatzMaterial.xsd ${DST_DIR}/Schemas/CFS_ErsatzMaterial.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PDELatticeBoltzmann.xsd ${DST_DIR}/Schemas/CFS_PDELatticeBoltzmann.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_Optimization.xsd ${DST_DIR}/Schemas/CFS_Optimization.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_PythonKernel.xsd ${DST_DIR}/Schemas/CFS_PythonKernel.xsd
  COMMENT "Updating XML schemata for openCFS simulation..."
)

SET(SRC_DIR ${CFS_SOURCE_DIR}/share/xml/CFS-Material)
SET(DST_DIR ${CFS_BINARY_DIR}/share/xml/CFS-Material)

ADD_CUSTOM_TARGET(share_xml_cfsmat ALL
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/CFS_Material.xsd ${DST_DIR}/CFS_Material.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_MATacoustic.xsd ${DST_DIR}/Schemas/CFS_MATacoustic.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_MATelectric.xsd ${DST_DIR}/Schemas/CFS_MATelectric.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_MATflow.xsd ${DST_DIR}/Schemas/CFS_MATflow.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_MATelecConduction.xsd ${DST_DIR}/Schemas/CFS_MATelecConduction.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_MATelecQuasistatic.xsd ${DST_DIR}/Schemas/CFS_MATelecQuasistatic.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_MATtest.xsd ${DST_DIR}/Schemas/CFS_MATtest.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_MATheatConduction.xsd ${DST_DIR}/Schemas/CFS_MATheatConduction.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_MATmagnetic.xsd ${DST_DIR}/Schemas/CFS_MATmagnetic.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_MATmagStrict.xsd ${DST_DIR}/Schemas/CFS_MATmagStrict.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_MATmechanical.xsd ${DST_DIR}/Schemas/CFS_MATmechanical.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_MATsmooth.xsd ${DST_DIR}/Schemas/CFS_MATsmooth.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_MAThysteresis.xsd ${DST_DIR}/Schemas/CFS_MAThysteresis.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_MATpiezo.xsd ${DST_DIR}/Schemas/CFS_MATpiezo.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_MATpyroelectric.xsd ${DST_DIR}/Schemas/CFS_MATpyroelectric.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_MATthermoelastic.xsd ${DST_DIR}/Schemas/CFS_MATthermoelastic.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFS_Misc.xsd ${DST_DIR}/Schemas/CFS_Misc.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Examples/MaterialDataBase.xml ${DST_DIR}/Examples/MaterialDataBase.xml
  COMMENT "Updating XML schemata for openCFS materials..."
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
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFSDat_BaseFilter.xsd ${DST_DIR}/Schemas/CFSDat_BaseFilter.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFSDat_BinOpFilt.xsd ${DST_DIR}/Schemas/CFSDat_BinOpFilt.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFSDat_Tensor.xsd ${DST_DIR}/Schemas/CFSDat_Tensor.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFSDat_DatMeshInput.xsd ${DST_DIR}/Schemas/CFSDat_DatMeshInput.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFSDat_DatMeshOutput.xsd ${DST_DIR}/Schemas/CFSDat_DatMeshOutput.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFSDat_DefOn.xsd ${DST_DIR}/Schemas/CFSDat_DefOn.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFSDat_GaussDerivative.xsd ${DST_DIR}/Schemas/CFSDat_GaussDerivative.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFSDat_VolumeMultiplication.xsd ${DST_DIR}/Schemas/CFSDat_VolumeMultiplication.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFSDat_Differentiation.xsd ${DST_DIR}/Schemas/CFSDat_Differentiation.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFSDat_DifferentiationType.xsd ${DST_DIR}/Schemas/CFSDat_DifferentiationType.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFSDat_FilterList.xsd ${DST_DIR}/Schemas/CFSDat_FilterList.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFSDat_GridType.xsd ${DST_DIR}/Schemas/CFSDat_GridType.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFSDat_Interpolation.xsd ${DST_DIR}/Schemas/CFSDat_Interpolation.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFSDat_InterpolationRBF.xsd ${DST_DIR}/Schemas/CFSDat_InterpolationRBF.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFSDat_InterpolationType.xsd ${DST_DIR}/Schemas/CFSDat_InterpolationType.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFSDat_NearestNeighbour.xsd ${DST_DIR}/Schemas/CFSDat_NearestNeighbour.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFSDat_FEBased.xsd ${DST_DIR}/Schemas/CFSDat_FEBased.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFSDat_ResultDef.xsd ${DST_DIR}/Schemas/CFSDat_ResultDef.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFSDat_RotDomain.xsd ${DST_DIR}/Schemas/CFSDat_RotDomain.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFSDat_RotSubstantialTimeDeriv.xsd ${DST_DIR}/Schemas/CFSDat_RotSubstantialTimeDeriv.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFSDat_SimplePoint.xsd ${DST_DIR}/Schemas/CFSDat_SimplePoint.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFSDat_FftFilter.xsd ${DST_DIR}/Schemas/CFSDat_FftFilter.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFSDat_TDeriv1.xsd ${DST_DIR}/Schemas/CFSDat_TDeriv1.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFSDat_FIRFilter.xsd ${DST_DIR}/Schemas/CFSDat_FIRFilter.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFSDat_TMean.xsd ${DST_DIR}/Schemas/CFSDat_TMean.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFSDat_Aeroacoustic.xsd ${DST_DIR}/Schemas/CFSDat_Aeroacoustic.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFSDat_SNGR.xsd ${DST_DIR}/Schemas/CFSDat_SNGR.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFSDat_TimeType.xsd ${DST_DIR}/Schemas/CFSDat_TimeType.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFSDat_AeroacousticType.xsd ${DST_DIR}/Schemas/CFSDat_AeroacousticType.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFSDat_PostLighthill.xsd ${DST_DIR}/Schemas/CFSDat_PostLighthill.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFSDat_TemporalBlending.xsd ${DST_DIR}/Schemas/CFSDat_TemporalBlending.xsd
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/Schemas/CFSDat_DiffScheme.xsd ${DST_DIR}/Schemas/CFSDat_DiffScheme.xsd
  COMMENT "Updating XML schemata for CFS-Dat..."
)

install(DIRECTORY "${CFS_BINARY_DIR}/share/xml" DESTINATION "share/" USE_SOURCE_PERMISSIONS PATTERN "CMakeFiles" EXCLUDE PATTERN "Makefile" EXCLUDE PATTERN "*.dir" EXCLUDE PATTERN "*.cmake" EXCLUDE)
