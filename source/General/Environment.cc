#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>

#include "Environment.hh"
#include "Utils/tools.hh"
#include "Domain/Domain.hh"
#include "def_use_blas.hh"
#include "def_use_openmp.hh"

#ifdef USE_OPENMP
  #include <omp.h>
#endif


using std::to_string;

namespace CoupledField {

  // Define global objects 
  Domain* domain = NULL;
  // Initialisation of some global pointers
  //WriteInfo *Info = NULL;

  // Definition of finite element type mappings
  static EnumTuple complexPartTuples[] = {
      EnumTuple(Global::INTEGER,  "INTEGER (does not belong into ComplexPart!)"),
      EnumTuple(Global::REAL,     "RealPart"),
      EnumTuple(Global::IMAG,     "ImaginaryPart"),
      EnumTuple(Global::COMPLEX,  "Complex")
  };

  Enum<Global::ComplexPart> Global::complexPart = \
  Enum<Global::ComplexPart>("Parts of a complex number",
                              sizeof(complexPartTuples) / sizeof(EnumTuple),
                              complexPartTuples);




  /** instantiate Enum<> objects, Fill the content in SetEnvironmentEnums() below */
  Enum<SolutionType> SolutionTypeEnum;
  Enum<MaterialType> MaterialTypeEnum;
  Enum<MaterialClass> MaterialClassEnum;
  Enum<MaterialTensorNotation> tensorNotation;
  Enum<ApproxCurveType> ApproxCurveTypeEnum;
  Enum<NonLinMethodType> NonLinMethodTypeEnum;
  Enum<FEMatrixType> feMatrixType;

  UInt MAX_NUM_FE_MATRICES;
  //give default value of 1 for the OMP threads
  UInt CFS_NUM_THREADS = 1;

  void SetEnvironmentEnums()
  {
    // SolutionType is the largest type - other types follow below

    SolutionTypeEnum.SetName("SolutionTypeEnum");
    //mechanics
    SolutionTypeEnum.Add(MECH_DISPLACEMENT, "mechDisplacement");
    SolutionTypeEnum.Add(MECH_ELEM_POROSITY, "mechElemPorosity");
    SolutionTypeEnum.Add(MECH_NORMAL_DISPLACEMENT, "mechNormalDisplacement");

    SolutionTypeEnum.Add(MECH_ACCELERATION, "mechAcceleration");
    SolutionTypeEnum.Add(MECH_VELOCITY, "mechVelocity");
    SolutionTypeEnum.Add(MECH_RHS_LOAD, "mechRhsLoad");

    SolutionTypeEnum.Add(MECH_STRESS, "mechStress");
    SolutionTypeEnum.Add(MECH_THERMAL_STRESS, "mechThermalStress");
    SolutionTypeEnum.Add(MECH_STRAIN, "mechStrain");
    SolutionTypeEnum.Add(MECH_IRR_STRAIN, "mechIrrStrain");
    SolutionTypeEnum.Add(MECH_IRR_STRESS, "mechIrrStress");
    SolutionTypeEnum.Add(MECH_THERMAL_STRAIN, "mechThermalStrain");
    SolutionTypeEnum.Add(MECH_STRUCT_INTENSTIY, "mechStructIntensity");
    SolutionTypeEnum.Add(MECH_NORMAL_STRUCT_INTENSITY, "mechNormalStructIntensity");
    SolutionTypeEnum.Add(MECH_PRINCIPAL_STRESS, "mechPrincipalStress");
    SolutionTypeEnum.Add(MECH_PRINCIPAL_STRESS_MIN, "mechPrincipalStressMin");
    SolutionTypeEnum.Add(MECH_PRINCIPAL_STRESS_MAX, "mechPrincipalStressMax");
    SolutionTypeEnum.Add(MECH_PRINCIPAL_STRESS_MED, "mechPrincipalStressMed");
    SolutionTypeEnum.Add(MECH_PRINCIPAL_STRESS_MIN_SCAL, "mechPrincipalStressMinScalar");
    SolutionTypeEnum.Add(MECH_PRINCIPAL_STRESS_MAX_SCAL, "mechPrincipalStressMaxScalar");
    SolutionTypeEnum.Add(MECH_PRINCIPAL_STRESS_MED_SCAL, "mechPrincipalStressMedScalar");
    SolutionTypeEnum.Add(MECH_PRINCIPAL_STRAIN, "mechPrincipalStrain");
    SolutionTypeEnum.Add(MECH_PRINCIPAL_STRAIN_MIN, "mechPrincipalStrainMin");
    SolutionTypeEnum.Add(MECH_PRINCIPAL_STRAIN_MAX, "mechPrincipalStrainMax");
    SolutionTypeEnum.Add(MECH_PRINCIPAL_STRAIN_MED, "mechPrincipalStrainMed");
    SolutionTypeEnum.Add(MECH_PRINCIPAL_STRAIN_MIN_SCAL, "mechPrincipalStrainMinScalar");
    SolutionTypeEnum.Add(MECH_PRINCIPAL_STRAIN_MAX_SCAL, "mechPrincipalStrainMaxScalar");
    SolutionTypeEnum.Add(MECH_PRINCIPAL_STRAIN_MED_SCAL, "mechPrincipalStrainMedScalar");
    SolutionTypeEnum.Add(VON_MISES_STRESS, "vonMisesStress");
    SolutionTypeEnum.Add(VON_MISES_STRAIN, "vonMisesStrain");
    SolutionTypeEnum.Add(MECH_KIN_ENERGY_DENS, "mechKinEnergyDensity");
    SolutionTypeEnum.Add(MECH_DEFORM_ENERGY_DENS, "mechDeformEnergyDensity");
    SolutionTypeEnum.Add(MECH_TOTAL_ENERGY_DENS, "mechTotalEnergyDensity");
    SolutionTypeEnum.Add(MECH_KIN_ENERGY, "mechKinEnergy");
    SolutionTypeEnum.Add(MECH_DEFORM_ENERGY, "mechDeformEnergy");
    SolutionTypeEnum.Add(MECH_TOTAL_ENERGY, "mechTotalEnergy");
    SolutionTypeEnum.Add(MECH_POWER, "mechPower");
    SolutionTypeEnum.Add(MECH_DEF_SURF_VOLUME, "mechDisplacedSurfVolume");
    SolutionTypeEnum.Add(MECH_FORCE, "mechForce");
    SolutionTypeEnum.Add(MECH_NORMAL_STRESS, "mechNormalStress");
    SolutionTypeEnum.Add(MECH_DYADIC_STRAIN, "mechDyadicStrain");
    SolutionTypeEnum.Add(MECH_DYADIC_STRAIN_SUM, "mechDyadicStrainSum");
    SolutionTypeEnum.Add(MECH_QUAD_DISP, "mechQuadDisplacement");
    SolutionTypeEnum.Add(MECH_QUAD_DISP_SUM, "mechQuadDisplacementSum");


    SolutionTypeEnum.Add(MECH_PSEUDO_DENSITY, "mechPseudoDensity"); // shall be replaced by PSEUDO_DENSITY in the future
    SolutionTypeEnum.Add(PSEUDO_DENSITY, "pseudoDensity");
    SolutionTypeEnum.Add(PHYSICAL_PSEUDO_DENSITY, "physicalPseudoDensity");
    SolutionTypeEnum.Add(MECH_SHAPE, "mechShape");
    SolutionTypeEnum.Add(MECH_TENSOR_TRACE, "mechTensorTrace");
    SolutionTypeEnum.Add(MECH_TENSOR, "mechTensor");
    SolutionTypeEnum.Add(MECH_TENSOR_HILL_MANDEL, "mechTensorHillMandel");
    SolutionTypeEnum.Add(MECH_ELEM_VOL, "mechElemVol");
    
    SolutionTypeEnum.Add(MECH_VELOCITY_ELEM, "mechVelocityElem");


    //electrostatics / elctric current conduction
    SolutionTypeEnum.Add(ELEC_POTENTIAL, "elecPotential");
    SolutionTypeEnum.Add(ELEC_FIELD_INTENSITY, "elecFieldIntensity");
    SolutionTypeEnum.Add(ELEC_FIELD_INTENSITY_SURF, "elecFieldIntensitySurf");
    SolutionTypeEnum.Add(ELEC_FIELD_INTENSITY_TRANSVERSAL, "elecFieldIntensityTransversal");
    SolutionTypeEnum.Add(ELEC_FIELD_INTENSITY_LONGITUDINAL, "elecFieldIntensityLongitudinal");
    SolutionTypeEnum.Add(ELEC_CURRENT_FIELD_INTENSITY, "elecCurrentDensitySurf");
    SolutionTypeEnum.Add(ELEC_CURRENT_SURF, "elecCurrentSurf");
    SolutionTypeEnum.Add(DISPLACEMENT_CURRENT_FIELD_INTENSITY, "displacementCurrentDensity");
    SolutionTypeEnum.Add(DISPLACEMENT_CURRENT_SURF, "displacementCurrent");
    SolutionTypeEnum.Add(ELECTRIC_AND_DISPLACEMENT_CURRENT_DENSITY, "electricAndDisplacemetCurrentDensity");
    SolutionTypeEnum.Add(ELEC_POLARIZATION, "elecPolarization");
    SolutionTypeEnum.Add(ELEC_CURRENT_DENSITY, "elecCurrentDensity");
    SolutionTypeEnum.Add(ELEC_GRAD_V_INT, "elecGradVInt");
    SolutionTypeEnum.Add(ELEC_NORMAL_CURRENT_DENSITY, "elecNormalCurrentDensity");
    SolutionTypeEnum.Add(DISPLACEMENT_NORMAL_CURRENT_DENSITY, "displacementNormalCurrentDensity");
    SolutionTypeEnum.Add(ELECTRIC_AND_DISPLACEMENT_NORMAL_CURRENT_DENSITY,"elecAndDisplacementNormalCurrentDensity");
    SolutionTypeEnum.Add(ELEC_CURRENT, "elecCurrent");
    SolutionTypeEnum.Add(ELEC_AND_DISPLACEMENT_CURRENT, "elecAndDisplacementCurrent");
    SolutionTypeEnum.Add(ELEC_POWER_DENSITY, "elecPowerDensity");
    SolutionTypeEnum.Add(ELEC_POWER, "elecPower");
    SolutionTypeEnum.Add(ELEC_FORCE, "elecForce");
    SolutionTypeEnum.Add(ELEC_FORCE_VWP, "elecForceVWP");
    SolutionTypeEnum.Add(ELEC_CHARGE, "elecCharge");
    SolutionTypeEnum.Add(ELEC_SURFACE_CHARGE_DENSITY, "elecSurfaceChargeDensity");
    SolutionTypeEnum.Add(ELEC_FLUX_DENSITY, "elecFluxDensity");
    SolutionTypeEnum.Add(ELEC_ENERGY, "elecEnergy");
    SolutionTypeEnum.Add(ELEC_ENERGY_DENSITY, "elecEnergyDensity");
    SolutionTypeEnum.Add(ELEC_FORCE_DENSITY, "elecForceDensity");
    SolutionTypeEnum.Add(ELEC_RHS_LOAD, "elecRhsLoad");
    SolutionTypeEnum.Add(ELEC_ELEM_PERMITTIVITY, "elecElemPermittivity");
    SolutionTypeEnum.Add(ELEC_COND_TENSOR, "ElectCond");

    SolutionTypeEnum.Add(ELEC_PSEUDO_POLARIZATION, "elecPseudoPolarization");
    SolutionTypeEnum.Add(ELEC_PHYSICAL_PSEUDO_DENSITY, "elecPhysicalPseudoDensity");
    SolutionTypeEnum.Add(ELEC_TENSOR, "elecTensor");
    SolutionTypeEnum.Add(ELEC_TENSOR_TRACE, "elecTensorTrace");

    //smoothing PDE
    SolutionTypeEnum.Add(SMOOTH_DISPLACEMENT, "smoothDisplacement");
    SolutionTypeEnum.Add(SMOOTH_VELOCITY, "smoothVelocity");
    // SolutionTypeEnum.Add(SMOOTH_ACCELERATION, "smoothAcceleration");
    SolutionTypeEnum.Add(SMOOTH_ZERO_PRESSURE, "smoothZeroStress");
    SolutionTypeEnum.Add(SMOOTH_STRAIN, "smoothStrain");
    SolutionTypeEnum.Add(SMOOTH_CONTACT_FORCE_DENSITY, "smoothContactForceDensity");
    SolutionTypeEnum.Add(SMOOTH_CONTACT_FORCE, "smoothContactForce");
    SolutionTypeEnum.Add(SMOOTH_DEFORM_ENERGY_DENS, "smoothDeformEnergyDensity");
    SolutionTypeEnum.Add(SMOOTH_DEFORM_ENERGY, "smoothDeformEnergy");
    
    //acoustics
    SolutionTypeEnum.Add(ACOU_PRESSURE, "acouPressure");
    SolutionTypeEnum.Add(ACOU_ACCELERATION, "acouAcceleration");
    SolutionTypeEnum.Add(ACOU_POTENTIAL, "acouPotential");
    SolutionTypeEnum.Add(ACOU_VELOCITY, "acouVelocity");
    SolutionTypeEnum.Add(ACOU_POSITION, "acouPosition");
    SolutionTypeEnum.Add(ACOU_NORMAL_VELOCITY, "acouNormalVelocity");
    SolutionTypeEnum.Add(ACOU_PRESSURE_DERIV_1, "acouPressureD1");
    SolutionTypeEnum.Add(ACOU_PRESSURE_DERIV_2, "acouPressureD2");
    SolutionTypeEnum.Add(ACOU_FORCE, "acouForce");
    SolutionTypeEnum.Add(ACOU_POTENTIAL_DERIV_1, "acouPotentialD1");
    SolutionTypeEnum.Add(ACOU_POTENTIAL_DERIV_2, "acouPotentialD2");
    SolutionTypeEnum.Add(ACOU_RHS_LOAD, "acouRhsLoad");
    SolutionTypeEnum.Add(ACOU_RHS_LOAD_DENSITY, "acouRhsLoadDensity");
    SolutionTypeEnum.Add(ACOU_RHS_LOAD_DENSITY_VECTOR, "acouRhsLoadDensityVector");
    SolutionTypeEnum.Add(ACOU_RHS_LOADP, "acouRhsLoadP");
    SolutionTypeEnum.Add(ACOU_APE_RHS_LOAD, "apeRhsLoad");
    SolutionTypeEnum.Add(ACOU_VORTEX_RHS_LOAD, "vortexRhsLoad");
    SolutionTypeEnum.Add(ACOU_DIV_LH_TENSOR, "acouDivLighthillTensor");
    SolutionTypeEnum.Add(ACOU_RHSVAL, "acouRHSval");
    SolutionTypeEnum.Add(ACOUSURF_RHSVAL, "acouSurfRHSval");
    SolutionTypeEnum.Add(ACOU_ELEM_SPEED_OF_SOUND,"acouSpeedOfSound");
    SolutionTypeEnum.Add(ACOU_POWERDENSITY, "acouPowerDensity");
    SolutionTypeEnum.Add(ACOU_POWER, "acouPower");
    SolutionTypeEnum.Add(ACOU_POWER_PLANEWAVE, "acouPowerPlaneWave");
    SolutionTypeEnum.Add(ACOU_INTENSITY, "acouIntensity");
    SolutionTypeEnum.Add(ACOU_NORMAL_INTENSITY, "acouNormalIntensity");
    SolutionTypeEnum.Add(ACOU_NORMAL_INTENSITY_PLANEWAVE, "acouNormalIntensityPlaneWave");
    SolutionTypeEnum.Add(ACOU_SURFINTENSITY, "acouSurfIntensity");
    SolutionTypeEnum.Add(ACOU_SURFPRESSURE, "acouSurfPressure");
    SolutionTypeEnum.Add(ACOU_POT_ENERGY, "acouPotEnergy");
    SolutionTypeEnum.Add(ACOU_KIN_ENERGY, "acouKinEnergy");
    SolutionTypeEnum.Add(ACOU_PMLAUXVEC,"acouPmlAuxVec");
    SolutionTypeEnum.Add(ACOU_PMLAUXSCALAR, "acouPmlAuxScalar");

    SolutionTypeEnum.Add(ACOU_MIXED_MASS_LOAD, "acouMixedMassLoad");
    SolutionTypeEnum.Add(ACOU_MIXED_MOMENTUM_LOAD, "acouMixedMomentumLoad");
    SolutionTypeEnum.Add(ACOU_LAMB_RHS, "acouLambRhs");

    // Auxiliary variables for TDEF formulation (currently limited to 15 poles)
    for (unsigned int i = 0; i < 15; i++)
      SolutionTypeEnum.Add((SolutionType)(ACOU_TDEF_PHI_C_1 + i), "acouTDEFPhiC_" + std::to_string(i));

    for (unsigned int i = 0; i < 15; i++)
      SolutionTypeEnum.Add((SolutionType)(ACOU_TDEF_PSI_C_1 + i), "acouTDEFPsiC_" + std::to_string(i));

    for (unsigned int i = 0; i < 15; i++)
      SolutionTypeEnum.Add((SolutionType)(ACOU_TDEF_PHI_V_1 + i), "acouTDEFPhiV_" + std::to_string(i));

    for (unsigned int i = 0; i < 15; i++)
      SolutionTypeEnum.Add((SolutionType)(ACOU_TDEF_PSI_V_1 + i), "acouTDEFPsiV_" + std::to_string(i));

    //Splitting
    SolutionTypeEnum.Add(SPLIT_SCALAR, "splitScalar");
    SolutionTypeEnum.Add(SPLIT_VECTOR, "splitVector");
    SolutionTypeEnum.Add(SPLIT_RHS_LOAD, "splitRhsLoad");
    SolutionTypeEnum.Add(SPLIT_SCALAR_VELOCITY, "splitScalVel");
    SolutionTypeEnum.Add(SPLIT_VECTOR_VELOCITY, "splitVectVel");
    SolutionTypeEnum.Add(SPLIT_POT_ENERGY, "splitPotEnergy");
    SolutionTypeEnum.Add(SPLIT_LAMB, "splitLamb");
    SolutionTypeEnum.Add(SPLIT_DIVLAMB, "splitDivLamb");

    //water waves
    SolutionTypeEnum.Add(WATER_PRESSURE, "waterPressure");
    SolutionTypeEnum.Add(WATER_POSITION, "waterPosition");
    SolutionTypeEnum.Add(WATER_PMLAUXVEC,"waterPmlAuxVec");
    SolutionTypeEnum.Add(WATER_PMLAUXSCALAR, "waterPmlAuxScalar");
    SolutionTypeEnum.Add(WATER_RHS_LOAD, "waterRhsLoad");
    SolutionTypeEnum.Add(WATER_PRES_TENS, "waterPressureTensor");
    SolutionTypeEnum.Add(WATER_SURFACE_TRACTION, "waterSurfaceTraction");
    SolutionTypeEnum.Add(WATER_SURFACE_FORCE, "waterSurfaceForce");
    SolutionTypeEnum.Add(WATER_SURFACE_TORQUE, "waterSurfaceTorque");
    SolutionTypeEnum.Add(WATER_TDT, "waterTorqueDensityTensor");
    SolutionTypeEnum.Add(WATER_SURFACE_TORQUE_DENSITY, "waterSurfaceTorqueDensity");

    //magnetics
    SolutionTypeEnum.Add(MAG_POTENTIAL, "magPotential");
    SolutionTypeEnum.Add(MAG_POTENTIAL_ADJ, "magPotentialAdj");
    SolutionTypeEnum.Add(MAG_POTENTIAL_DERIV1, "magPotentialD1");
    SolutionTypeEnum.Add(MAG_POTENTIAL_ADJ_DERIV1, "magPotentialAdjD1");
    SolutionTypeEnum.Add(MAG_TOTAL_POTENTIAL, "magTotalPotential");
    SolutionTypeEnum.Add(MAG_REDUCED_POTENTIAL, "magReducedPotential");
    SolutionTypeEnum.Add(MAG_RHS_LOAD, "magRhsLoad");
    SolutionTypeEnum.Add(FLUX_INDUCED_STRAIN, "fluxIndStrain");

    SolutionTypeEnum.Add(MAG_FLUX_DENSITY, "magFluxDensity");
    SolutionTypeEnum.Add(MAG_FLUX_DENSITY_SURF, "magFluxDensitySurf");
    SolutionTypeEnum.Add(MAG_FLUX, "magFlux");
    SolutionTypeEnum.Add(MAG_NORMAL_FLUX_DENSITY, "magNormalFluxDensity");
    SolutionTypeEnum.Add(MAG_AVERAGED_FLUX_DENSITY, "magAveragedFluxDensity"); 
    SolutionTypeEnum.Add(MAG_FIELD_INTENSITY, "magFieldIntensity");
    SolutionTypeEnum.Add(MAG_AVERAGED_FIELD_INTENSITY, "magAveragedFieldIntensity");       
    SolutionTypeEnum.Add(MAG_EDDY_CURRENT_DENSITY, "magEddyCurrentDensity");
    SolutionTypeEnum.Add(MAG_COIL_CURRENT_DENSITY, "magCoilCurrentDensity");
    SolutionTypeEnum.Add(MAG_TOTAL_CURRENT_DENSITY, "magTotalCurrentDensity");
    SolutionTypeEnum.Add(MAG_JOULE_LOSS_POWER_DENSITY, "magJouleLossPowerDensity");
    SolutionTypeEnum.Add(MAG_JOULE_LOSS_POWER_DENSITY_ON_NODES, "magJouleLossPowerDensityOnNodes");
    SolutionTypeEnum.Add(MAG_JOULE_LOSS_POWER, "magJouleLossPower");
    SolutionTypeEnum.Add(MAG_POTENTIAL_DIV, "magPotentialDiv");
    SolutionTypeEnum.Add(MAG_POTENTIAL_GRAD, "magPotentialGrad");
    SolutionTypeEnum.Add(MAG_CURL_ADJ, "magCurlAdj");
    SolutionTypeEnum.Add(MAG_FORCE_LORENTZ_DENSITY, "magForceLorentzDensity");
    SolutionTypeEnum.Add(MAG_FORCE_LORENTZ_DENSITY_STATIC, "magForceLorentzDensityStatic");
    SolutionTypeEnum.Add(MAG_FORCE_LORENTZ_DENSITY_HARMONIC, "magForceLorentzDensityHarmonic");    
    SolutionTypeEnum.Add(MAG_FORCE_MAXWELL_DENSITY, "magForceMaxwellDensity");
    SolutionTypeEnum.Add(MAG_NORMALFORCE_MAXWELL_DENSITY, "magNormalForceMaxwellDensity");
    SolutionTypeEnum.Add(MAG_TANGENTIALFORCE_MAXWELL_DENSITY, "magTangentialForceMaxwellDensity");
    SolutionTypeEnum.Add(MAG_EDDY_POWER_DENSITY, "magEddyPowerDensity");
    SolutionTypeEnum.Add(MAG_ENERGY_DENSITY, "magEnergyDensity");
    SolutionTypeEnum.Add(MAG_CORE_LOSS_DENSITY, "magCoreLossDensity");
    SolutionTypeEnum.Add(MAG_CORE_LOSS, "magCoreLoss");
    SolutionTypeEnum.Add(MAG_GRAD_ADJ_PARAM, "magGradAdjParam");
    SolutionTypeEnum.Add(MAG_GRAD_ADJ_PARAM1, "magGradAdjParam1");
    SolutionTypeEnum.Add(MAG_GRAD_ADJ_PARAM2, "magGradAdjParam2");
    SolutionTypeEnum.Add(MAG_GRAD_ADJ_PARAM3, "magGradAdjParam3");
    SolutionTypeEnum.Add(MAG_GRAD_ADJ_PARAM4, "magGradAdjParam4");

    SolutionTypeEnum.Add(MAG_FORCE_VWP, "magForceVWP");
    SolutionTypeEnum.Add(MAG_FORCE_LORENTZ, "magForceLorentz");
    SolutionTypeEnum.Add(MAG_FORCE_MAXWELL, "magForceMaxwell");
    SolutionTypeEnum.Add(MAG_ENERGY, "magEnergy");
    SolutionTypeEnum.Add(MAG_EDDY_POWER, "magEddyPower");
    SolutionTypeEnum.Add(MAG_EDDY_CURRENT, "magEddyCurrent");
    SolutionTypeEnum.Add(MAG_EDDY_CURRENT1, "magEddyCurrent1");
    SolutionTypeEnum.Add(MAG_EDDY_CURRENT2, "magEddyCurrent2");
    SolutionTypeEnum.Add(MAG_TOTAL_CURRENT, "magTotalCurrent");
    SolutionTypeEnum.Add(MAG_ELEM_PERMEABILITY, "magElemPermeability");
    SolutionTypeEnum.Add(MAG_ELEM_RELUCTIVITY, "magElemReluctivity");
    SolutionTypeEnum.Add(MAG_MAGNETIZATION, "magMagnetization");
    SolutionTypeEnum.Add(MAG_POLARIZATION, "magPolarization");

    // for magnetic coil optimization
    SolutionTypeEnum.Add(RHS_PSEUDO_DENSITY, "rhsPseudoDensity");
    SolutionTypeEnum.Add(PHYSICAL_RHS_PSEUDO_DENSITY, "physicalRhsPseudoDensity");
    SolutionTypeEnum.Add(ELEC_POTENTIAL_DERIV_1, "elecPotentialD1");

    // magnetic - coil results
    SolutionTypeEnum.Add(COIL_CURRENT, "coilCurrent");
    SolutionTypeEnum.Add(COIL_CURRENT_DERIV1, "coilCurrentD1");
    SolutionTypeEnum.Add(COIL_VOLTAGE, "coilVoltage");
    SolutionTypeEnum.Add(COIL_VOLTAGE_INTEGRAL, "coilVoltageIntegral");

    SolutionTypeEnum.Add(COIL_INDUCED_VOLTAGE, "coilInducedVoltage");
    SolutionTypeEnum.Add(COIL_INDUCTANCE, "coilInductance");
    SolutionTypeEnum.Add(COIL_LINKED_FLUX, "coilLinkedFlux");

    //heat conduction
    SolutionTypeEnum.Add(HEAT_TEMPERATURE, "heatTemperature");
    SolutionTypeEnum.Add(HEAT_MEAN_TEMPERATURE, "heatMeanTemperature");
    SolutionTypeEnum.Add(HEAT_TEMPERATURE_D1, "heatTemperatureD1");
    SolutionTypeEnum.Add(HEAT_FLUX_DENSITY, "heatFluxDensity");
    SolutionTypeEnum.Add(HEAT_FLUX_INTENSITY, "heatFluxIntensity");
    SolutionTypeEnum.Add(HEAT_FLUX, "heatFlux");
    SolutionTypeEnum.Add(HEAT_RHS_LOAD, "heatRhsLoad");
    SolutionTypeEnum.Add(HEAT_SOURCE_DENSITY, "heatSourceDensity");
    SolutionTypeEnum.Add(HEAT_CONDUCTIVITY_TENSOR_HOM,"heatConductivityTensorHomogenized");

    //fluidMech
    SolutionTypeEnum.Add(FLUID_FORCE, "fluidForce");
    SolutionTypeEnum.Add(MEAN_FLUIDMECH_VELOCITY, "meanFluidMechVelocity");
    SolutionTypeEnum.Add(MEAN_FLUIDMECH_VELOCITY_NORMAL, "meanFluidMechVelocityNormal");
    SolutionTypeEnum.Add(DIV_MEAN_FLUIDMECH_VELOCITY, "divMeanFluidMechVelocity");
    SolutionTypeEnum.Add(FLUIDMECH_VORTICITY, "fluidMechVorticity");
    SolutionTypeEnum.Add(FLUIDMECH_VELOCITY, "fluidMechVelocity");
    SolutionTypeEnum.Add(FLUIDMECH_MESH_VELOCITY, "fluidMechMeshVelocity");
    SolutionTypeEnum.Add(FLUIDMECH_MESH_VELOCITY_NODE, "fluidMechMeshVelocityNode");
    SolutionTypeEnum.Add(FLUIDMECH_MESH_VELOCITY_ELEM, "fluidMechMeshVelocityElem");
    SolutionTypeEnum.Add(FLUIDMECH_TOTAL_VELOCITY_ELEM, "fluidMechTotalVelocityElem");
    SolutionTypeEnum.Add(FLUIDMECH_NORMAL_VELOCITY, "fluidMechNormalVelocity");
    SolutionTypeEnum.Add(FLUIDMECH_PRESSURE, "fluidMechPressure");
    SolutionTypeEnum.Add(FLUIDMECH_ZERO_PRESSURE, "fluidMechZeroPressure");
    SolutionTypeEnum.Add(FLUIDMECH_VELOCITY_DERIV_1, "fluidMechVelocity_deriv1");
    SolutionTypeEnum.Add(FLUIDMECH_PRESSURE_DERIV_1, "fluidMechPressure_deriv1");
    SolutionTypeEnum.Add(FLUIDMECH_VELOCITY_DERIV_2, "fluidMechVelocity_deriv2");
    SolutionTypeEnum.Add(FLUIDMECH_PRESSURE_DERIV_2, "fluidMechPressure_deriv2");

    SolutionTypeEnum.Add(FLUIDMECH_PRESSURE_TIME_DERIV_1, "fluidMechPressure_timeDeriv1");
    SolutionTypeEnum.Add(FLUIDMECH_PRESSURE_TIME_DERIV_2, "fluidMechPressure_timeDeriv2");

    SolutionTypeEnum.Add(FLUIDMECH_FORCE, "fluidMechForce");
    SolutionTypeEnum.Add(FLUIDMECH_DENSITY, "fluidMechDensity");
    SolutionTypeEnum.Add(FLUIDMECH_TEMP, "fluidMechTemp");
    SolutionTypeEnum.Add(FLUIDMECH_TKE, "fluidMechTKE");
    SolutionTypeEnum.Add(FLUIDMECH_TDR, "fluidMechTDR");
    SolutionTypeEnum.Add(FLUIDMECH_TEF, "fluidMechTEF");
    SolutionTypeEnum.Add(FLUIDMECH_STRESS, "fluidMechStress");
    SolutionTypeEnum.Add(FLUIDMECH_COMP_STRESS, "fluidMechCompressibleStress");
    SolutionTypeEnum.Add(FLUIDMECH_VISC_STRESS, "fluidMechViscousStress");
    SolutionTypeEnum.Add(FLUIDMECH_PRES_TENS, "fluidMechPressureTensor");
    SolutionTypeEnum.Add(FLUIDMECH_TOTAL_STRESS, "fluidMechTotalStress");
    SolutionTypeEnum.Add(FLUIDMECH_SURFACE_TRACTION, "fluidMechSurfaceTraction");
    SolutionTypeEnum.Add(FLUIDMECH_STRAINRATE, "fluidMechStrainRate");
    SolutionTypeEnum.Add(FLUIDMECH_WVT, "fluidMechWVT");
    SolutionTypeEnum.Add(FLUIDMECH_WVT_DENSITY, "fluidMechWVTDensity");
    SolutionTypeEnum.Add(FLUIDMECH_WVT_PHI, "fluidMechWVTPhi");
    SolutionTypeEnum.Add(FLUIDMECH_WVT_DENSITY_PHI, "fluidMechWVTDensityPhi");
    SolutionTypeEnum.Add(FLUIDMECH_WVT_U1, "fluidMechWVT_u1");
    SolutionTypeEnum.Add(FLUIDMECH_WVT_U2, "fluidMechWVT_u2");
    SolutionTypeEnum.Add(FLUIDMECH_WVT_F, "fluidMechWVT_f");

    SolutionTypeEnum.Add(FLUIDMECH_VISCOUS_DISS_POWER_DENS, "fluidMechViscousDissipationDensity");
    SolutionTypeEnum.Add(FLUIDMECH_VISCOUS_DISS_POWER_DENS_DIV, "fluidMechViscousDissipationDensityDivergencePart");
    SolutionTypeEnum.Add(FLUIDMECH_VISCOUS_DISS_POWER_DENS_STRAIN, "fluidMechViscousDissipationDensityTotalStrainPart");
    SolutionTypeEnum.Add(FLUIDMECH_VISCOUS_DISS_POWER, "fluidMechViscousDissipation");

    SolutionTypeEnum.Add(FLUIDMECH_INTENSITY, "fluidMechIntensity");
    SolutionTypeEnum.Add(FLUIDMECH_INTENSITY_PRESSURE_ONLY, "fluidMechIntensityPressureOnly");
    SolutionTypeEnum.Add(FLUIDMECH_SURFINTENSITY, "fluidMechSurfaceIntensity");
    SolutionTypeEnum.Add(FLUIDMECH_SURFINTENSITY_PRESSURE_ONLY, "fluidMechSurfaceIntensityPressureOnly");
    SolutionTypeEnum.Add(FLUIDMECH_NORMAL_INTENSITY, "fluidMechNormalIntensity");
    SolutionTypeEnum.Add(FLUIDMECH_NORMAL_INTENSITY_PRESSURE_ONLY, "fluidMechNormalIntensityPressureOnly");
    SolutionTypeEnum.Add(FLUIDMECH_POWER, "fluidMechPower");
    SolutionTypeEnum.Add(FLUIDMECH_POWER_PRESSURE_ONLY, "fluidMechPowerPressureOnly");
    SolutionTypeEnum.Add(FLUIDMECH_SURFIMPEDANCE, "fluidMechSurfaceImpedance");
    SolutionTypeEnum.Add(FLUIDMECH_IMPEDANCE, "fluidMechImpedance");

    SolutionTypeEnum.Add(LAMBDA_K, "lambda_k");
    SolutionTypeEnum.Add(VOLUME, "volume");

    // TEST PDE
    SolutionTypeEnum.Add(TEST_DOF, "testDof");
    SolutionTypeEnum.Add(TEST_FIELD, "testField");
    SolutionTypeEnum.Add(TEST_RHS_LOAD, "testRhsLoad");

    // optimization
    SolutionTypeEnum.Add(HOMOGENIZED_TENSOR, "homogenizedTensor");
    // the actual result type is given in result descriptions
    // in the xml file in the optimization element.
    for(unsigned int i = 0; i < 66; i++)
      SolutionTypeEnum.Add( (SolutionType) (OPT_RESULT_1 + i), "optResult_" + std::to_string(i+1));
      // SolutionTypeEnum.Add(OPT_RESULT_1, "optResult_1");

    // independent
    SolutionTypeEnum.Add(LAGRANGE_MULT, "lagrangeMultiplier");
    SolutionTypeEnum.Add(LAGRANGE_MULT_DERIV_1, "lagrangeMultiplierD1");
    SolutionTypeEnum.Add(LAGRANGE_MULT_DERIV_2, "lagrangeMultiplierD2");
    SolutionTypeEnum.Add(LAGRANGE_MULT_1, "lagrangeMultiplier1");
    // evaluates the spacial gradient of the solution at the nodes.
    // common for all PDEs, no unit
    SolutionTypeEnum.Add(GRAD_ACOU_SOLUTION, "gradAcousticSolution"); // independent on acoustic formulation
    SolutionTypeEnum.Add(GRAD_ELEC_POTENTIAL, "gradElecPotential");
    SolutionTypeEnum.Add(GRAD_ELEC_POTENTIAL_DERIV1, "gradElecPotentialD1");
    SolutionTypeEnum.Add(ELEM_DENSITY, "density");
    // PML parameters for all formulations
    SolutionTypeEnum.Add(PML_DAMP_FACTOR, "pmlDampFactor"); // eigenvalues of the PML tensor, stacked in a vector 
                                                            // (=diagonal entries for classic and shifted PML, and 
                                                            // diagonal entries of the 'inner' tensor for curvilinear PML)
    SolutionTypeEnum.Add(PML_DISTANCE, "pmlDistance"); // distance between point and the PML interface 
    // PML parameters for curvilinear PML only
    SolutionTypeEnum.Add(PML_TENSOR, "pmlTensor"); // whole PML operation tensor of the curvilinear PML
    SolutionTypeEnum.Add(PML_DETERMINANT, "pmlDeterminant"); // the determinant of the PML tensor

    // General (grid related) results
    SolutionTypeEnum.Add(ELEM_LOC_DIR, "localDirection");
    SolutionTypeEnum.Add(JACOBIAN, "jacobian");
    SolutionTypeEnum.Add(ASPECT_RATIO, "aspectRatio");

    //LBM velocity
    SolutionTypeEnum.Add(LBM_PROBABILITY_DISTRIBUTION, "LBMProbabilityDistribution");
    SolutionTypeEnum.Add(LBM_VELOCITY, "LBMVelocity");
    SolutionTypeEnum.Add(LBM_DENSITY, "LBMDensity");
    SolutionTypeEnum.Add(LBM_PRESSURE, "LBMPressure");
    SolutionTypeEnum.Add(LBM_PHYSICAL_PSEUDO_DENSITY, "LBMPhysicalPseudoDensity");

    //Geometry
    SolutionTypeEnum.Add(NODE_NORMAL, "NodeNormal");
    SolutionTypeEnum.Add(COOSY_X,"CooSy-default-x");
    SolutionTypeEnum.Add(COOSY_Y,"CooSy-default-y");
    SolutionTypeEnum.Add(SURFACE_NORMAL,"SurfaceNormal");
    SolutionTypeEnum.Add(AREA,"area");
    SolutionTypeEnum.Add(ETA, "eta");
    SolutionTypeEnum.Add(XI, "xi");

    // ==== Initialization of Material Constants ====
    MaterialTypeEnum.Add( NO_MATERIAL, "noMaterial" );

    // -- ACOUSTIC --
    MaterialTypeEnum.Add( ACOU_BULK_MODULUS, "Acoustic_Bulk_Modulus" );
    MaterialTypeEnum.Add( ACOU_SOUND_SPEED, "Acoustic_Sound_Speed" );
    MaterialTypeEnum.Add( ACOU_BOVERA, "BoverA" );
    MaterialTypeEnum.Add( ACOU_IMPEDANCE_REAL_VAL, "acouImpedanceRealValue");
    MaterialTypeEnum.Add( ACOU_IMPEDANCE_IMAG_VAL, "acouImpedanceImagValue");
    MaterialTypeEnum.Add( ACOU_TDEF_INVDENS_CONST, "acouTDEF_invDens_Const");
    MaterialTypeEnum.Add( ACOU_TDEF_INVDENS_A, "acouTDEF_invDens_CoefA");
    MaterialTypeEnum.Add( ACOU_TDEF_INVDENS_ALPHA, "acouTDEF_invDens_CoefAlpha");
    MaterialTypeEnum.Add( ACOU_TDEF_INVDENS_B, "acouTDEF_invDens_CoefB");
    MaterialTypeEnum.Add( ACOU_TDEF_INVDENS_BETA, "acouTDEF_invDens_CoefBeta");
    MaterialTypeEnum.Add( ACOU_TDEF_INVDENS_C, "acouTDEF_invDens_CoefC");
    MaterialTypeEnum.Add( ACOU_TDEF_INVDENS_GAMMA, "acouTDEF_invDens_CoefGamma");
    MaterialTypeEnum.Add( ACOU_TDEF_INVBLK_CONST, "acouTDEF_invBlk_Const");
    MaterialTypeEnum.Add( ACOU_TDEF_INVBLK_A, "acouTDEF_invBlk_CoefA");
    MaterialTypeEnum.Add( ACOU_TDEF_INVBLK_ALPHA, "acouTDEF_invBlk_CoefAlpha");
    MaterialTypeEnum.Add( ACOU_TDEF_INVBLK_B, "acouTDEF_invBlk_CoefB");
    MaterialTypeEnum.Add( ACOU_TDEF_INVBLK_BETA, "acouTDEF_invBlk_CoefBeta");
    MaterialTypeEnum.Add( ACOU_TDEF_INVBLK_C, "acouTDEF_invBlk_CoefC");
    MaterialTypeEnum.Add( ACOU_TDEF_INVBLK_GAMMA, "acouTDEF_invBlk_CoefGamma");
//    MaterialTypeEnum.Add( IMP_HOLE_DIAM, "holeDiam");
//    MaterialTypeEnum.Add( IMP_PLATE_THICKNESS, "plateThick");
//    MaterialTypeEnum.Add( IMP_END_CORRECTION, "impEndCorrection");
//    MaterialTypeEnum.Add( POROSITY, "sigma");
//    MaterialTypeEnum.Add( BETA, "beta");
//    MaterialTypeEnum.Add( MPP_VOLUME_DEPTH, "mppVolDepth");
//    MaterialTypeEnum.Add( FLOW_MACH_NUMBER, "flowMachNr");
//    MaterialTypeEnum.Add( ACOU_ALPHA, "AcousticAlpha" );
//    MaterialTypeEnum.Add( FRACTIONAL_ALG, "FractionalAlg" );
//    MaterialTypeEnum.Add( FRACTIONAL_MEMORY, "FractionalMemory" );
//    MaterialTypeEnum.Add( FRACTIONAL_INTERPOL, "FractionalInterpol" );
//    MaterialTypeEnum.Add( FRACTIONAL_EXPONENT, "FractionalExponent" );

    // -- Electrostatic --
    MaterialTypeEnum.Add( ELEC_PERMITTIVITY_TENSOR, "Electric_Permittivity_Tensdor" );
    MaterialTypeEnum.Add( ELEC_PERMITTIVITY_SCALAR, "Electric_Permittivity_Scalar" );
    MaterialTypeEnum.Add( ELEC_PERMITTIVITY_1, "Electric_Permittivity_1" );
    MaterialTypeEnum.Add( ELEC_PERMITTIVITY_2, "Electric_Permittivity_2" );
    MaterialTypeEnum.Add( ELEC_PERMITTIVITY_3, "Electric_Permittivity_3" );

    // -- Electrostatic Jiles Hysterese
    MaterialTypeEnum.Add( ELEC_PS_JILES, "elec_Ps_Jiles" );
    MaterialTypeEnum.Add( ELEC_ALPHA_JILES, "elec_alpha_Jiles" );
    MaterialTypeEnum.Add( ELEC_A_JILES, "elec_a_Jiles" );
    MaterialTypeEnum.Add( ELEC_K_JILES, "elec_k_Jiles" );
    MaterialTypeEnum.Add( ELEC_C_JILES, "elec_c_Jiles" );

    // -- Electric Conduction --
    MaterialTypeEnum.Add( ELEC_CONDUCTIVITY_TENSOR, "Electric_Conductivity_Tensor" );
    MaterialTypeEnum.Add( ELEC_CONDUCTIVITY_SCALAR, "Electric_Conductivity_Scalar" );
    MaterialTypeEnum.Add( ELEC_CONDUCTIVITY_1, "Electric_Conductivity_1" );
    MaterialTypeEnum.Add( ELEC_CONDUCTIVITY_2, "Electric_Conductivity_2" );
    MaterialTypeEnum.Add( ELEC_CONDUCTIVITY_3, "Electric_Conductivity_3" );

    // -- Flow --
    MaterialTypeEnum.Add( FLUID_ADIABATIC_EXPONENT, "Flow_Adiabatic_Exponent");
    MaterialTypeEnum.Add( FLUID_DYNAMIC_VISCOSITY, "Flow_Dynamic_Viscosity" );
    MaterialTypeEnum.Add( FLUID_KINEMATIC_VISCOSITY, "Flow_Kinematic_Viscosity" );
    MaterialTypeEnum.Add( FLUID_BULK_VISCOSITY, "Flow_Bulk_Viscosity" );
    MaterialTypeEnum.Add( FLUID_BULK_MODULUS, "Flow_Bulk_Modulus" );

    // -- Heat Conduction --
    MaterialTypeEnum.Add( HEAT_CONDUCTIVITY_TENSOR, "Heat_Conductivity_Tensor" );
    MaterialTypeEnum.Add( HEAT_CONDUCTIVITY_SCALAR, "Heat_Conductivity_Scalar" );
    MaterialTypeEnum.Add( HEAT_CONDUCTIVITY_1, "Heat_Conductivity_1" );
    MaterialTypeEnum.Add( HEAT_CONDUCTIVITY_2, "Heat_Conductivity_2" );
    MaterialTypeEnum.Add( HEAT_CONDUCTIVITY_3, "Heat_Conductivity_3" );
    MaterialTypeEnum.Add( HEAT_CAPACITY, "Heat_Capacity" );
    MaterialTypeEnum.Add( HEAT_REF_TEMPERATURE, "Heat_Reference_Temperature" );

    // -- Magnetic --
    MaterialTypeEnum.Add( MAG_PERMEABILITY_TENSOR, "Magnetic_Permeability_Tensor" );
    MaterialTypeEnum.Add( MAG_PERMEABILITY_SCALAR, "Magnetic_Permeability_Scalar" );
    MaterialTypeEnum.Add( MAG_PERMEABILITY_1, "Magnetic_Permeability_1" );
    MaterialTypeEnum.Add( MAG_PERMEABILITY_2, "Magnetic_Permeability_2" );
    MaterialTypeEnum.Add( MAG_PERMEABILITY_3, "Magnetic_Permeability_3" );
    MaterialTypeEnum.Add( MAG_RELUCTIVITY_TENSOR, "Magnetic_Reluctivity_Tensor" );
    MaterialTypeEnum.Add( MAG_RELUCTIVITY_SCALAR, "Magnetic_Reluctivity_Scalar" );
    MaterialTypeEnum.Add( MAG_RELUCTIVITY_DERIV, "Magnetic_Reluctivity_Derivative" );
    MaterialTypeEnum.Add( MAG_RELUCTIVITY_DERIV_P1, "Magnetic_Reluctivity_DerivativeP1" );
    MaterialTypeEnum.Add( MAG_RELUCTIVITY_DERIV_P2, "Magnetic_Reluctivity_DerivativeP2" );
    MaterialTypeEnum.Add( MAG_RELUCTIVITY_DERIV_P3, "Magnetic_Reluctivity_DerivativeP3" );
    MaterialTypeEnum.Add( MAG_RELUCTIVITY_DERIV_P4, "Magnetic_Reluctivity_DerivativeP4" );    
    MaterialTypeEnum.Add( MAG_CONDUCTIVITY_TENSOR, "Magnetic_Conductivity_Tensor" );
    MaterialTypeEnum.Add( MAG_CONDUCTIVITY_SCALAR, "Magnetic_Conductivity_Scalar" );
    MaterialTypeEnum.Add( MAG_PERMITTIVITY_SCALAR, "Magnetic_Permittivity_Scalar" );
    MaterialTypeEnum.Add( MAG_PERMITTIVITY_TENSOR, "Magnetic_Permittivity_Tensor" );
    MaterialTypeEnum.Add( MAG_PERMITTIVITY_1, "Magnetic_Permittivity_1" );
    MaterialTypeEnum.Add( MAG_PERMITTIVITY_2, "Magnetic_Permittivity_2" );
    MaterialTypeEnum.Add( MAG_PERMITTIVITY_3, "Magnetic_Permittivity_3" );
    MaterialTypeEnum.Add( MAG_CONDUCTIVITY_1, "Magnetic_Conductivity_1" );
    MaterialTypeEnum.Add( MAG_CONDUCTIVITY_2, "Magnetic_Conductivity_2" );
    MaterialTypeEnum.Add( MAG_CONDUCTIVITY_3, "Magnetic_Conductivity_3" );
//    MaterialTypeEnum.Add( MAGNETIZATION, "Magnetization");
    MaterialTypeEnum.Add( MAG_CORE_LOSS_PER_MASS, "Magnetic_Core_Loss_Per_Mass" );
    MaterialTypeEnum.Add( MAG_BH_VALUES, "Magnetic_BH_Curve");
    MaterialTypeEnum.Add( MAG_BH_VALUES_1, "Magnetic_BH_Curve_1");
    MaterialTypeEnum.Add( MAG_BH_VALUES_2, "Magnetic_BH_Curve_2");
    MaterialTypeEnum.Add( MAG_BH_VALUES_3, "Magnetic_BH_Curve_3");
    MaterialTypeEnum.Add( MAG_BH_DATA_ACCURACY, "Magnetic_BH_Data_Accuracy");
    MaterialTypeEnum.Add( MAG_BH_MAX_APPROX_VAL, "Magnetic_BH_Max_Approx_Value");
    // Magnetic EB Hysteresis Parameters
    MaterialTypeEnum.Add(MAG_PS_EB, "mag_Ps_EB");
    MaterialTypeEnum.Add(MAG_A_EB, "mag_A_EB");
    MaterialTypeEnum.Add(MAG_MU0_EB, "mag_mu0_EB");
    MaterialTypeEnum.Add(MAG_NUMS_EB, "mag_numS_EB");
    MaterialTypeEnum.Add(MAG_CHI_FACTOR_EB, "mag_chi_factor_EB");

    // -- Mechanical --
    MaterialTypeEnum.Add( MECH_STIFFNESS_TENSOR, "Mechanic_Stiffness_Tensor" );
    MaterialTypeEnum.Add( MECH_KMODULUS, "Mechanic_Bulk_Modulus" );
    MaterialTypeEnum.Add( MECH_LAME_MU, "Mechanic_Lame_Mu" );
    MaterialTypeEnum.Add( MECH_LAME_LAMBDA, "Mechanic_Lame_Lambda" );
    MaterialTypeEnum.Add( MECH_EMODULUS, "Mechanic_Elasticity_Modulus" );
    MaterialTypeEnum.Add( MECH_EMODULUS_1, "Mechanic_Elasticity_Modulus_1" );
    MaterialTypeEnum.Add( MECH_EMODULUS_2, "Mechanic_Elasticity_Modulus_2" );
    MaterialTypeEnum.Add( MECH_EMODULUS_3, "Mechanic_Elasticity_Modulus_3" );
    MaterialTypeEnum.Add( MECH_POISSON, "Mechanic_Poisson_Ratio" );
    MaterialTypeEnum.Add( MECH_POISSON_3, "Mechanic_Poisson_Ratio_3" );
    MaterialTypeEnum.Add( MECH_POISSON_12, "Mechanic_PoissonRation_12" );
    MaterialTypeEnum.Add( MECH_POISSON_23, "Mechanic_PoissonRation_23" );
    MaterialTypeEnum.Add( MECH_POISSON_13, "Mechanic_PoissonRation_13" );
    MaterialTypeEnum.Add( MECH_GMODULUS, "Mechanic_Shear_Modulus" );
    MaterialTypeEnum.Add( MECH_GMODULUS_3, "Mechanic_Shear_Modulus_3" );
    MaterialTypeEnum.Add( MECH_GMODULUS_23, "Mechanic_Shear_Modulus_23" );
    MaterialTypeEnum.Add( MECH_GMODULUS_13, "Mechanic_Shear_Modulus_13" );
    MaterialTypeEnum.Add( MECH_GMODULUS_12, "Mechanic_Shear_Modulus_12" );
    // Viscoelasticity
    MaterialTypeEnum.Add( MECH_BULK_RELAX_TIMES, "Mechanic_Bulk_Relaxation_Times"),
    MaterialTypeEnum.Add( MECH_BULK_RELAX_MODULI, "Mechanic_Bulk_Relaxation_Moduli"),
    MaterialTypeEnum.Add( MECH_VISCO_BULK_INITIAL, "Mechanic_Visco_Initial_Bulk_Modulus"),
    MaterialTypeEnum.Add( MECH_SHEAR_RELAX_TIMES, "Mechanic_Shear_Relaxation_Times"),
    MaterialTypeEnum.Add( MECH_SHEAR_RELAX_MODULI, "Mechanic_Shear_Relaxation_Moduli"),
    MaterialTypeEnum.Add( MECH_VISCO_SHEAR_INITIAL, "Mechanic_Visco_Initial_Shear_Modulus"),
    MaterialTypeEnum.Add( MECH_VISCO_STIFFNESS_TENSOR, "Mechanic_Visco_Stiffness_Tensor"),
    MaterialTypeEnum.Add( MECH_VISCO_LONGTERM_TENSOR, "Mechanic_Visco_Longterm_Tensor"),
    // Themal Expansion
    MaterialTypeEnum.Add( MECH_THERMAL_EXPANSION_TENSOR, "Mechanic_Thermal_Expansion_Tensor" );
    MaterialTypeEnum.Add( MECH_THERMAL_EXPANSION_SCALAR, "Mechanic_Thermal_Expansion_Scalar" );
    MaterialTypeEnum.Add( MECH_THERMAL_EXPANSION_1, "Mechanic_Thermal_Expansion_1" );
    MaterialTypeEnum.Add( MECH_THERMAL_EXPANSION_2, "Mechanic_Thermal_Expansion_2" );
    MaterialTypeEnum.Add( MECH_THERMAL_EXPANSION_3, "Mechanic_Thermal_Expansion_3" );
    MaterialTypeEnum.Add( MECH_TE_REFTEMPERATURE, "Mechanic_Thermal_Expansion_Reference_Temperature");
//    MaterialTypeEnum.Add( PYROCOEFFICIENT_TENSOR, "Pyrocoefficient_Tensor" );

    // -- Test PDE --
    MaterialTypeEnum.Add( TEST_ALPHA, "Test_Alpha" );
    MaterialTypeEnum.Add( TEST_BETA, "Test_BETA" );

    // -- Piezo-electric --
    MaterialTypeEnum.Add( PIEZO_TENSOR, "PiezoTensor" );
    // Piezo Micro Model
    MaterialTypeEnum.Add( MEAN_TEMPERATURE, "meanTemperature" );
    MaterialTypeEnum.Add( SPON_POLARIZATION, "sponPolarization" );
    MaterialTypeEnum.Add( SPON_STRAIN, "sponStrain" );
    MaterialTypeEnum.Add( DRIVING_FORCE_90, "Driving_Force_90");
    MaterialTypeEnum.Add( DRIVING_FORCE_180, "Driving_Force_180");
    MaterialTypeEnum.Add( RATE_CONSTANT, "rateConstant" );
    MaterialTypeEnum.Add( VISCO_PLASTIC_INDEX, "viscoPlasticIndex" );
    MaterialTypeEnum.Add( SATURATION_INDEX, "saturationIndex" );
    MaterialTypeEnum.Add( VOLUME_FRAC_INIT, "volumeFracInit" );
    MaterialTypeEnum.Add( EFIELD0, "Efield0" );
    MaterialTypeEnum.Add( STRESS0, "Stress0" );
    MaterialTypeEnum.Add( DCOUPLE0, "dCouple0" );
    MaterialTypeEnum.Add( SCALE_FORCE_ELEC, "scaleForceElec" );
    MaterialTypeEnum.Add( SCALE_FORCE_MECH, "scaleForceMech" );
    MaterialTypeEnum.Add( SCALE_FORCE_COUPLE, "scaleForceCouple" );

    // -- Magnetostriction --
    MaterialTypeEnum.Add( MAGSTRICT_RELUCTIVITY,"Magstrict_reluctivity");
    MaterialTypeEnum.Add( MAGNETOSTRICTION_TENSOR_h, "Magnetostriction_Tensor_h" );
    MaterialTypeEnum.Add( MAGNETOSTRICTION_TENSOR_h_mech, "Magnetostriction_Tensor_h_mech" );
    MaterialTypeEnum.Add( MAGNETOSTRICTION_TENSOR_h_mag, "Magnetostriction_Tensor_h_mag" );
    MaterialTypeEnum.Add( A_JILES, "aJiles" );
    MaterialTypeEnum.Add( ALPHA_JILES, "alphaJiles" );
    MaterialTypeEnum.Add( K_JILES, "kJiles" );
    MaterialTypeEnum.Add( C_JILES, "cJiles" );
    MaterialTypeEnum.Add( JILES_TEST, "dummyParamForJiles" );
    MaterialTypeEnum.Add( Y_REMANENCE, "Yremanence" );

    // -- General --
    MaterialTypeEnum.Add( DENSITY, "Density" );

    // Rayleigh Damping
    MaterialTypeEnum.Add( RAYLEIGH_ALPHA, "Rayleigh_Alpha" );
    MaterialTypeEnum.Add( RAYLEIGH_BETA, "Rayleigh_Beta" );
    MaterialTypeEnum.Add( LOSS_TANGENS_DELTA, "Loss_TangensDelta" );

    // General Material Nonlinearity
    MaterialTypeEnum.Add( NONLIN_DEPENDENCY, "nonLinDependency" );

    // -- Hysteresis --
    MaterialTypeEnum.Add( X_SATURATION, "Xsaturation" );
    MaterialTypeEnum.Add( Y_SATURATION, "Ysaturation" );
    MaterialTypeEnum.Add( S_SATURATION, "Ssaturation" );
    MaterialTypeEnum.Add( PREISACH_WEIGHTS, "preisachWeights" );
    MaterialTypeEnum.Add( PREISACH_WEIGHTS_DIM, "preisachWeights_dim" );
    MaterialTypeEnum.Add( PREISACH_WEIGHTS_CONSTVALUE, "preisachWeights_constValue" );
    MaterialTypeEnum.Add( PREISACH_WEIGHTS_TYPE, "preisachWeights_type" );
    MaterialTypeEnum.Add( PREISACH_WEIGHTS_MUDAT_A, "preisachWeights_mudatA" );
    MaterialTypeEnum.Add( PREISACH_WEIGHTS_MUDAT_H, "preisachWeights_mudath" );
    MaterialTypeEnum.Add( PREISACH_WEIGHTS_MUDAT_H2, "preisachWeights_mudath2" );
    MaterialTypeEnum.Add( PREISACH_WEIGHTS_MUDAT_SIGMA, "preisachWeights_mudatsigma" );
    MaterialTypeEnum.Add( PREISACH_WEIGHTS_MUDAT_SIGMA2, "preisachWeights_mudatsigma2" );
    MaterialTypeEnum.Add( PREISACH_WEIGHTS_MUDAT_ETA, "preisachWeights_mudateta" );
    MaterialTypeEnum.Add( PREISACH_WEIGHTS_ANHYSTCOUNTINGTOOUTPUTSAT, "preisachWeights_anhystCountsToSat" );
    MaterialTypeEnum.Add( PREISACH_WEIGHTS_MUDAT_PARAMSFORHALFRANGE, "preisachWeights_mudat_paramsforhalfrange" );
    MaterialTypeEnum.Add( PREISACH_WEIGHTS_FOR_MAYERGOYZ_VECTOR, "preisachWeights_forMayergoyzVector" );
    MaterialTypeEnum.Add( PREISACH_WEIGHTS_ANHYST_ONLY, "preisachWeights_anhyst_ONLY" );
    MaterialTypeEnum.Add( PREISACH_WEIGHTS_ANHYST_A, "preisachWeights_anhystA" );
    MaterialTypeEnum.Add( PREISACH_WEIGHTS_ANHYST_B, "preisachWeights_anhystB" );
    MaterialTypeEnum.Add( PREISACH_WEIGHTS_ANHYST_C, "preisachWeights_anhystC" );
    MaterialTypeEnum.Add( PREISACH_WEIGHTS_ANHYST_D, "preisachWeights_anhystD" );
    MaterialTypeEnum.Add( PREISACH_WEIGHTS_ANHYST_PARAMSFORHALFRANGE, "preisachWeights_anhyst_paramsforhalfrange" );
    MaterialTypeEnum.Add( PREISACH_WEIGHTS_TENSOR, "preisachWeights_givenTensor" );
    MaterialTypeEnum.Add( PREISACH_MAYERGOYZ_NUM_DIR, "preisach_mayergoyz_num_dir" );
    MaterialTypeEnum.Add( PREISACH_MAYERGOYZ_ISOTROPIC, "preisach_mayergoyz_isotropic" );
    MaterialTypeEnum.Add( PREISACH_MAYERGOYZ_CLIPOUTPUT, "preisach_mayergoyz_clip_output" );
    MaterialTypeEnum.Add( MAYERGOYZ_STARTAXIS_X, "preisach_mayergoyz_start_x" );
    MaterialTypeEnum.Add( MAYERGOYZ_STARTAXIS_Y, "preisach_mayergoyz_start_y" );
    MaterialTypeEnum.Add( MAYERGOYZ_STARTAXIS_Z, "preisach_mayergoyz_start_z" );
    MaterialTypeEnum.Add( MAYERGOYZ_LOSSPARAM_A, "preisach_mayergoyz_lossa" );
    MaterialTypeEnum.Add( MAYERGOYZ_LOSSPARAM_B, "preisach_mayergoyz_lossb" );
    MaterialTypeEnum.Add( MAYERGOYZ_USEABSDPHI, "preisach_mayergoyz_useAbsDPhi" );
    MaterialTypeEnum.Add( MAYERGOYZ_NORMALIZEXINEXP, "preisach_mayergoyz_normalizeXinExp" );
    MaterialTypeEnum.Add( MAYERGOYZ_RESTRICTIONOFPSI, "preisach_mayergoyz_restrictionOfPsi" );
    MaterialTypeEnum.Add( MAYERGOYZ_SCALINGOFXINEXP, "preisach_mayergoyz_scalingOfXinExp" );
    MaterialTypeEnum.Add( HYST_TYPE_IS_PREISACH, "isPreisachType" );
    MaterialTypeEnum.Add( PREISACH_PRESCRIBEOUTPUT, "preisach_prescribe_output" );
    MaterialTypeEnum.Add( PREISACH_SCALEINITIALSTATE, "preisach_scale_initial_state" );
    MaterialTypeEnum.Add( SCALETOSAT, "preisach_scaletosat" );
    MaterialTypeEnum.Add( SCALETOSAT_STRAIN, "preisach_scaletosat_strain" );
    MaterialTypeEnum.Add( P_DIRECTION, "Pdirection" );
    MaterialTypeEnum.Add( P_DIRECTION_X, "PdirectionX" );
    MaterialTypeEnum.Add( P_DIRECTION_Y, "PdirectionY" );
    MaterialTypeEnum.Add( P_DIRECTION_Z, "PdirectionZ" );
    MaterialTypeEnum.Add( HYST_MODEL, "hystModel" );
    MaterialTypeEnum.Add( HYSTERESIS_DIM, "PreisachDim" );
    MaterialTypeEnum.Add( HYST_STRAIN_FORM, "strainForm" );

    MaterialTypeEnum.Add( HYST_IRRSTRAINS, "irrStrainForm" );
    MaterialTypeEnum.Add( HYST_IRRSTRAIN_C1, "irrStrainMuDatC1" );
    MaterialTypeEnum.Add( HYST_IRRSTRAIN_C2, "irrStrainMuDatC2" );
    MaterialTypeEnum.Add( HYST_IRRSTRAIN_C3, "irrStrainMuDatC3" );
    MaterialTypeEnum.Add( HYST_IRRSTRAIN_CI, "irrStrainMuDatCI" );
    MaterialTypeEnum.Add( HYST_IRRSTRAIN_CI_SIZE, "irrStrainMuDatCI_size" );
    MaterialTypeEnum.Add( HYST_IRRSTRAIN_D0, "irrStrainMuDatD0" );
    MaterialTypeEnum.Add( HYST_IRRSTRAIN_D1, "irrStrainMuDatD1" );
    MaterialTypeEnum.Add( HYST_IRRSTRAIN_SCALETOSAT, "irrStrainMuDat_scaleToSat" );
    MaterialTypeEnum.Add( HYST_IRRSTRAIN_PARAMSFORHALFRANGE, "irrStrain_paramsforhalfrange" );
    MaterialTypeEnum.Add( ROT_RESISTANCE, "RotResistance" );
    MaterialTypeEnum.Add( PRINT_PREISACH, "printOut" );
    MaterialTypeEnum.Add( PRINT_PREISACH_RESOLUTION, "bmpResolution" );
    MaterialTypeEnum.Add( IS_TESTING, "isTesting" );
    MaterialTypeEnum.Add( EVAL_VERSION, "evalVersion" );
    MaterialTypeEnum.Add( ANG_DISTANCE, "angularDistance" );
    MaterialTypeEnum.Add( ANG_CLIPPING, "angularClipping" );
    MaterialTypeEnum.Add( ANG_RESOLUTION, "angularResolution" );
    MaterialTypeEnum.Add( AMP_RESOLUTION, "amplitudeResolution" );

    MaterialTypeEnum.Add( MAX_NUM_IT_HYST_INV, "hystInvOuterNumIt" );
    MaterialTypeEnum.Add( VEC_HYST_INV_METHOD, "hystInvMethod" );
    MaterialTypeEnum.Add( RES_TOL_H_HYST_INV, "residualTolH" );
    MaterialTypeEnum.Add( RES_TOL_B_HYST_INV, "residualTolB" );
    MaterialTypeEnum.Add( RES_TOL_H_HYST_INV_ISREL, "residualTolH_isrel" );
    MaterialTypeEnum.Add( RES_TOL_B_HYST_INV_ISREL, "residualTolB_isrel" );

    MaterialTypeEnum.Add( ALPHA_REG_HYST_INV, "alphaRegStart" );
    MaterialTypeEnum.Add( ALPHA_REG_MIN_HYST_INV, "alphaRegMin" );
    MaterialTypeEnum.Add( ALPHA_REG_MAX_HYST_INV, "alphaRegMax" );
    MaterialTypeEnum.Add( MAX_NUM_REG_IT_HYST_INV, "hystInvRegNumIt" );
    MaterialTypeEnum.Add( ALPHA_LS_MIN_HYST_INV, "alphaLSMin" );
    MaterialTypeEnum.Add( ALPHA_LS_MAX_HYST_INV, "alphaLSMax" );
    MaterialTypeEnum.Add( MAX_NUM_LS_IT_HYST_INV, "hystInvLSNumIt" );
    MaterialTypeEnum.Add( STOP_INV_LS_AT_LOCAL_MIN, "stopInvLineSearchAtLocalMin" );
    MaterialTypeEnum.Add( JAC_RESOLUTION_HYST_INV, "jacobiResolution" );
    MaterialTypeEnum.Add( JAC_IMPLEMENTATION_HYST_INV, "jacobiImplementation" );
    MaterialTypeEnum.Add( TRUST_LOW_HYST_INV, "trustRegionLow" );
    MaterialTypeEnum.Add( TRUST_MID_HYST_INV, "trustRegionMid" );
    MaterialTypeEnum.Add( TRUST_HIGH_HYST_INV, "trustRegionHigh" );

    MaterialTypeEnum.Add( HYST_INV_PROJLM_MU, "hystInv_MU" );
    MaterialTypeEnum.Add( HYST_INV_PROJLM_RHO, "hystInv_RHO" );
    MaterialTypeEnum.Add( HYST_INV_PROJLM_BETA, "hystInv_BETA" );
    MaterialTypeEnum.Add( HYST_INV_PROJLM_SIGMA, "hystInv_SIGMA" );
    MaterialTypeEnum.Add( HYST_INV_PROJLM_GAMMA, "hystInv_GAMMA" );
    MaterialTypeEnum.Add( HYST_INV_PROJLM_TAU, "hystInv_TAU" );
    MaterialTypeEnum.Add( HYST_INV_PROJLM_C, "hystInv_C" );
    MaterialTypeEnum.Add( HYST_INV_PROJLM_P, "hystInv_P" );

    MaterialTypeEnum.Add( HYST_INV_FP_SAFETYFACTOR, "hystInv_FP_safetyfactor" );
    MaterialTypeEnum.Add( HYST_LOCAL_INVERSION_PRINT_WARNINGS, "hystInv_printWarnings" );

    MaterialTypeEnum.Add( INITIAL_STATE, "initialStates" );
    MaterialTypeEnum.Add( INITIAL_STATE_X, "initialStatesX" );
    MaterialTypeEnum.Add( INITIAL_STATE_Y, "initialStatesY" );
    MaterialTypeEnum.Add( INITIAL_STATE_Z, "initialStatesZ" );

    MaterialTypeEnum.Add( PRESCRIBED_MAGNETIZATION, "prescribedMagnetization" );
    MaterialTypeEnum.Add( PRESCRIBED_MAGNETIZATION_X, "prescribedMagnetization_x" );
    MaterialTypeEnum.Add( PRESCRIBED_MAGNETIZATION_Y, "prescribedMagnetization_y" );
    MaterialTypeEnum.Add( PRESCRIBED_MAGNETIZATION_Z, "prescribedMagnetization_z" );

    MaterialTypeEnum.Add( X_SATURATION_STRAIN, "XsaturationStrain" );
    MaterialTypeEnum.Add( Y_SATURATION_STRAIN, "YsaturationStrain" );
    MaterialTypeEnum.Add( PREISACH_WEIGHTS_STRAIN, "preisachWeights_strain" );
    MaterialTypeEnum.Add( PREISACH_WEIGHTS_DIM_STRAIN, "preisachWeights_dim_strain" );
    MaterialTypeEnum.Add( PREISACH_WEIGHTS_CONSTVALUE_STRAIN, "preisachWeights_constValue_strain" );
    MaterialTypeEnum.Add( PREISACH_WEIGHTS_TYPE_STRAIN, "preisachWeights_type_strain" );
    MaterialTypeEnum.Add( PREISACH_WEIGHTS_MUDAT_A_STRAIN, "preisachWeights_mudatA_STRAIN" );
    MaterialTypeEnum.Add( PREISACH_WEIGHTS_MUDAT_H_STRAIN, "preisachWeights_mudath_STRAIN" );
    MaterialTypeEnum.Add( PREISACH_WEIGHTS_MUDAT_H2_STRAIN, "preisachWeights_mudath2_STRAIN" );
    MaterialTypeEnum.Add( PREISACH_WEIGHTS_MUDAT_SIGMA_STRAIN, "preisachWeights_mudatsigma_STRAIN" );
    MaterialTypeEnum.Add( PREISACH_WEIGHTS_MUDAT_SIGMA2_STRAIN, "preisachWeights_mudatsigma2_STRAIN" );
    MaterialTypeEnum.Add( PREISACH_WEIGHTS_MUDAT_ETA_STRAIN, "preisachWeights_mudateta_STRAIN" );
    MaterialTypeEnum.Add( PREISACH_WEIGHTS_ANHYSTCOUNTINGTOOUTPUTSAT_STRAIN, "preisachWeights_anhystCountsToSat_STRAIN" );
    MaterialTypeEnum.Add( PREISACH_WEIGHTS_MUDAT_PARAMSFORHALFRANGE_STRAIN, "preisachWeights_mudat_paramsforhalfrange_STRAIN" );

    MaterialTypeEnum.Add( PREISACH_WEIGHTS_FOR_MAYERGOYZ_VECTOR_STRAIN, "preisachWeights_forMayergoyzVector_strain" );
    MaterialTypeEnum.Add( PREISACH_WEIGHTS_TENSOR_STRAIN, "preisachWeights_givenTensor_strain" );

    MaterialTypeEnum.Add( S_DIRECTION, "Sdirection" );
    MaterialTypeEnum.Add( S_DIRECTION_X, "SdirectionX" );
    MaterialTypeEnum.Add( S_DIRECTION_Y, "SdirectionY" );
    MaterialTypeEnum.Add( S_DIRECTION_Z, "SdirectionZ" );
    MaterialTypeEnum.Add( HYST_MODEL_STRAIN, "hystModel_strain" );
    MaterialTypeEnum.Add( HYSTERESIS_DIM_STRAIN, "PreisachDim_strain" );
    MaterialTypeEnum.Add( ROT_RESISTANCE_STRAIN, "RotResistance_strain" );
    MaterialTypeEnum.Add( EVAL_VERSION_STRAIN, "evalVersion_strain" );
    MaterialTypeEnum.Add( ANG_DISTANCE_STRAIN, "angularDistance_strain" );
    MaterialTypeEnum.Add( ANG_CLIPPING_STRAIN, "angularClipping_strain" );
    MaterialTypeEnum.Add( ANG_RESOLUTION_STRAIN, "angularResolution_strain" );
    MaterialTypeEnum.Add( AMP_RESOLUTION_STRAIN, "amplitudeResolution_strain" );
    MaterialTypeEnum.Add( IRRSTRAIN_REUSE_P, "reusePolarizationForStrain" );

    MaterialTypeEnum.Add( PREISACH_WEIGHTS_ANHYST_ONLY_STRAIN, "preisachWeights_anhyst_ONLY_strain" );
    MaterialTypeEnum.Add( PREISACH_WEIGHTS_ANHYST_A_STRAIN, "preisachWeights_anhystA_strain" );
    MaterialTypeEnum.Add( PREISACH_WEIGHTS_ANHYST_B_STRAIN, "preisachWeights_anhystB_strain" );
    MaterialTypeEnum.Add( PREISACH_WEIGHTS_ANHYST_C_STRAIN, "preisachWeights_anhystC_strain" );
    MaterialTypeEnum.Add( PREISACH_WEIGHTS_ANHYST_D_STRAIN, "preisachWeights_anhystD_strain" );
    MaterialTypeEnum.Add( PREISACH_WEIGHTS_ANHYST_PARAMSFORHALFRANGE_STRAIN, "preisachWeights_anhyst_paramsforhalfrange_strain" );

    MaterialTypeEnum.Add( PREISACH_MAYERGOYZ_NUM_DIR_STRAIN, "preisach_mayergoyz_num_dir_strain" );
    MaterialTypeEnum.Add( PREISACH_MAYERGOYZ_ISOTROPIC_STRAIN, "preisach_mayergoyz_isotropic_strain" );
    MaterialTypeEnum.Add( PREISACH_MAYERGOYZ_CLIPOUTPUT_STRAIN, "preisach_mayergoyz_clip_output_strain" );
    MaterialTypeEnum.Add( MAYERGOYZ_STARTAXIS_X_STRAIN, "preisach_mayergoyz_start_x_strain" );
    MaterialTypeEnum.Add( MAYERGOYZ_STARTAXIS_Y_STRAIN, "preisach_mayergoyz_start_y_strain" );
    MaterialTypeEnum.Add( MAYERGOYZ_STARTAXIS_Z_STRAIN, "preisach_mayergoyz_start_z_strain" );
    MaterialTypeEnum.Add( MAYERGOYZ_LOSSPARAM_A_STRAIN, "preisach_mayergoyz_lossa_strain" );
    MaterialTypeEnum.Add( MAYERGOYZ_LOSSPARAM_B_STRAIN, "preisach_mayergoyz_lossb_strain" );
    MaterialTypeEnum.Add( MAYERGOYZ_USEABSDPHI_STRAIN, "preisach_mayergoyz_useAbsDPhi_strain" );
    MaterialTypeEnum.Add( MAYERGOYZ_NORMALIZEXINEXP_STRAIN, "preisach_mayergoyz_normalizeXinExp_strain" );
    MaterialTypeEnum.Add( MAYERGOYZ_RESTRICTIONOFPSI_STRAIN, "preisach_mayergoyz_restrictionOfPsi_strain" );
    MaterialTypeEnum.Add( MAYERGOYZ_SCALINGOFXINEXP_STRAIN, "preisach_mayergoyz_scalingOfXinExp_strain" );
    MaterialTypeEnum.Add( HYST_TYPE_IS_PREISACH_STRAIN, "isPreisachType_strain" );

    MaterialTypeEnum.Add( HYST_COUPLING_DEFINED, "hystCouplingDefined" );

    MaterialTypeEnum.Add( TRACE_FORCE_CENTRALDIFF, "traceForceCentralDifferences" );
    MaterialTypeEnum.Add( TRACE_FORCE_RETRACING, "traceForceRetracing" );
    MaterialTypeEnum.Add( TRACE_JAC_RESOLUTION, "traceJacResolution" );

    // -- Smooth --
    MaterialTypeEnum.Add( SMOOTH_STIFFNESS_TENSOR, "Smooth_Stiffness_Tensor" );
    MaterialTypeEnum.Add( SMOOTH_KMODULUS, "Smooth_Bulk_Modulus" );
    MaterialTypeEnum.Add( SMOOTH_LAME_MU, "Smooth_Lame_Mu" );
    MaterialTypeEnum.Add( SMOOTH_LAME_LAMBDA, "Smooth_Lame_Lambda" );
    MaterialTypeEnum.Add( SMOOTH_EMODULUS, "Smooth_Elasticity_Modulus" );
    MaterialTypeEnum.Add( SMOOTH_EMODULUS_1, "Smooth_Elasticity_Modulus_1" );
    MaterialTypeEnum.Add( SMOOTH_EMODULUS_2, "Smooth_Elasticity_Modulus_2" );
    MaterialTypeEnum.Add( SMOOTH_EMODULUS_3, "Smooth_Elasticity_Modulus_3" );
    MaterialTypeEnum.Add( SMOOTH_POISSON, "Smooth_Poisson_Ratio" );
    MaterialTypeEnum.Add( SMOOTH_POISSON_3, "Smooth_Poisson_Ratio_3" );
    MaterialTypeEnum.Add( SMOOTH_POISSON_12, "Smooth_PoissonRation_12" );
    MaterialTypeEnum.Add( SMOOTH_POISSON_23, "Smooth_PoissonRation_23" );
    MaterialTypeEnum.Add( SMOOTH_POISSON_13, "Smooth_PoissonRation_13" );
    MaterialTypeEnum.Add( SMOOTH_GMODULUS, "Smooth_Shear_Modulus" );

    // ==== Initialization of Material Classes ====
    MaterialClassEnum.Add(NO_CLASS, "no material class");
    MaterialClassEnum.Add(TESTMAT, "test");
    MaterialClassEnum.Add(ELECTROMAGNETIC, "electromagnetic");
    MaterialClassEnum.Add(ELECTROSTATIC, "electrostatic");
    MaterialClassEnum.Add(ACOUSTIC, "acoustic");
    MaterialClassEnum.Add(FLOW, "fluid");
    MaterialClassEnum.Add(MECHANIC, "mechanical");
    MaterialClassEnum.Add(PIEZO, "piezo");
    MaterialClassEnum.Add(THERMIC, "thermal");
    MaterialClassEnum.Add(PYROELECTRIC, "pyroelectric");
    MaterialClassEnum.Add(THERMOELASTIC, "thermoelastic");
    MaterialClassEnum.Add(MAGNETOSTRICTIVE, "magnetostrictive");
    MaterialClassEnum.Add(ELECTRICCONDUCTION, "conductive");
    MaterialClassEnum.Add(SMOOTH, "smooth");

    // ==== Initialization of Material tensor notations ====
    tensorNotation.SetName("MaterialTensorNotation");
    tensorNotation.Add(NO_NOTATION, "no-notation");
    tensorNotation.Add(HILL_MANDEL, "hill_mandel");
    tensorNotation.Add(VOIGT, "voigt");

    // ==== Initialization of Matrix Types ====
    feMatrixType.Add( NOTYPE, "none" );
    feMatrixType.Add( SYSTEM, "system" );
    feMatrixType.Add( STIFFNESS, "stiffness");
    feMatrixType.Add( DAMPING, "damping" );
    feMatrixType.Add( DAMPING_AUX, "damping_aux" );
    feMatrixType.Add( CONVECTION, "convection");
    feMatrixType.Add( MASS, "mass" );
    feMatrixType.Add( AUXILIARY, "auxiliary" );
    feMatrixType.Add( MASS_UPDATE, "mass_update" );
    feMatrixType.Add( STIFFNESS_UPDATE, "stiffness_update");
    feMatrixType.Add( DAMPING_UPDATE, "damping_update" );
    feMatrixType.Add( SYSTEM_HYSTFREE, "system_hystfree" );
    feMatrixType.Add( SYSTEM_FIXPOINT, "system_fixpoint" );
    feMatrixType.Add( SYSTEM_FD_JACOBIAN, "system_fd_Jacobian" );
    feMatrixType.Add( SYSTEM_DELTAMAT_JACOBIAN, "system_deltaMat_Jacobian" );
    feMatrixType.Add( GEOMETRIC_STIFFNESS, "geometric_stiffness" );
    feMatrixType.Add( BACKUP, "backup" );

    // ==== Initialization of ApproxCurveTypes ====
    ApproxCurveTypeEnum.Add( NO_APPROX_TYPE , "No approximation" );
    ApproxCurveTypeEnum.Add( LIN_INTERPOLATE, "LinInterpolate" );
    ApproxCurveTypeEnum.Add( BILIN_INTERPOLATE, "BiLinInterpolate" );
    ApproxCurveTypeEnum.Add( TRILIN_INTERPOLATE, "TriLinInterpolate" );
    ApproxCurveTypeEnum.Add( CUBIC_SPLINES, "CubicSplines" );
    ApproxCurveTypeEnum.Add( SMOOTH_SPLINES, "smoothSplines" );
    ApproxCurveTypeEnum.Add( ANALYTIC, "analytic" );


    MAX_NUM_FE_MATRICES = feMatrixType.map.size() - 1;

    // ==== Initialization of NonLinMethodEnum ====
    NonLinMethodTypeEnum.Add( FIXEDPOINT, "fixPoint" );
    NonLinMethodTypeEnum.Add( NEWTON, "newton" );
    NonLinMethodTypeEnum.Add( HYST_FIXPOINT_IT, "HYST_fixPoint_IT" );
    NonLinMethodTypeEnum.Add( HYST_FIXPOINT_TS, "HYST_fixPoint_TS" );
    NonLinMethodTypeEnum.Add( HYST_DELTAMAT_IT, "HYST_deltaMat_IT" );
    NonLinMethodTypeEnum.Add( HYST_DELTAMAT_TS, "HYST_deltaMat_TS" );
    NonLinMethodTypeEnum.Add( HYST_DEBUG, "HYST_debug" );
  }

  void AddGenericSolution(std::string name, Domain* domain) {
    // use the given name and assign a generic result based on the internal counter
    // afterwards, return the new solutionType
    SolutionTypeEnum.Add( (SolutionType) (GENERIC_RESULT_0 + domain->GetGenericResultIndex()), name);
    // increment the counter
    domain->IncrementGenericResultIndex();
  }

  std::string GetSolAsString(std::string name) {
    SolutionType solType = SolutionTypeEnum.Parse(name);
    return SolutionTypeEnum.ToString(solType);
  }

  // SolutionType
  std::string MapSolTypeToUnit(SolutionType solType)
  {
    switch(solType)
    {

      case ACOU_FORCE:
      case FLUIDMECH_WVT_F:
      case FLUIDMECH_FORCE:
      case MAG_FORCE_VWP:
      case MECH_RHS_LOAD:
      case ELEC_FORCE:
      case SMOOTH_CONTACT_FORCE:
      case WATER_SURFACE_FORCE:
        return "N";
        break;

      case ACOU_POWER:
      case ACOU_POWER_PLANEWAVE:
      case ELEC_POWER:
      case HEAT_FLUX:
      case MAG_EDDY_POWER:
      case MAG_CORE_LOSS:
      case FLUIDMECH_VISCOUS_DISS_POWER:
      case FLUIDMECH_POWER:
      case FLUIDMECH_POWER_PRESSURE_ONLY:
        return "W";
         break;

      case ACOU_ENERGY:
      case ACOU_POT_ENERGY:
      case ACOU_KIN_ENERGY:
      case ELEC_ENERGY:
      case MAG_ENERGY:
      case MECH_KIN_ENERGY:
      case MECH_DEFORM_ENERGY:
      case MECH_TOTAL_ENERGY:
      case SMOOTH_DEFORM_ENERGY:
        return "Ws";
        break;

      case ACOU_NORMAL_INTENSITY_PLANEWAVE:
      case ACOU_INTENSITY:
      case ACOU_NORMAL_INTENSITY:
      case ACOU_SURFINTENSITY:
      case FLUIDMECH_INTENSITY:
      case FLUIDMECH_INTENSITY_PRESSURE_ONLY:
      case FLUIDMECH_SURFINTENSITY:
      case FLUIDMECH_SURFINTENSITY_PRESSURE_ONLY:
      case FLUIDMECH_NORMAL_INTENSITY:
      case FLUIDMECH_NORMAL_INTENSITY_PRESSURE_ONLY:
      case HEAT_FLUX_DENSITY:
      case HEAT_FLUX_INTENSITY:
        return "W/m^2";
        break;

      case HEAT_CONDUCTIVITY_TENSOR_HOM:
        return "W/(mK)";
        break;

      case SPLIT_POT_ENERGY:
      case ELEC_ENERGY_DENSITY:
      case SMOOTH_DEFORM_ENERGY_DENS:
        return "Ws/m^3";
        break;

      case ACOU_POSITION:
      case WATER_POSITION:
      case MECH_DISPLACEMENT:
      case MECH_NORMAL_DISPLACEMENT:
      case SMOOTH_DISPLACEMENT:
        return "m";
        break;

      case ACOU_NORMAL_VELOCITY:
      case ACOU_VELOCITY:
      case ACOU_ELEM_SPEED_OF_SOUND:
      case FLUIDMECH_VELOCITY:
      case MEAN_FLUIDMECH_VELOCITY:
      case MEAN_FLUIDMECH_VELOCITY_NORMAL:
      case FLUIDMECH_MESH_VELOCITY:
      case FLUIDMECH_MESH_VELOCITY_NODE:
      case FLUIDMECH_MESH_VELOCITY_ELEM:
      case FLUIDMECH_TOTAL_VELOCITY_ELEM:
      case MECH_VELOCITY:
      case MECH_VELOCITY_ELEM:
      case SMOOTH_VELOCITY:
      case SPLIT_VECTOR_VELOCITY:
      case SPLIT_SCALAR_VELOCITY:
        return "m/s";
        break;

      case MECH_ACCELERATION:
      // case SMOOTH_ACCELERATION:
        return "m/s^2";
        break;

      case ACOU_ACCELERATION:
      case ACOU_POTENTIAL:
      case SPLIT_SCALAR:
      case SPLIT_VECTOR:
        return "m^2/s";
        break;

      case ACOU_POTENTIAL_DERIV_1:
        return "m^2/s^2";
        break;

      case ACOU_POTENTIAL_DERIV_2:
      case FLUIDMECH_TEF:
        return "m^2/s^3";
        break;

      case MECH_DEF_SURF_VOLUME:
        return "m^3";
        break;
      case SPLIT_RHS_LOAD:
        return "m^3/s";
        break;

      case ACOU_PRESSURE:
      case ACOU_SURFPRESSURE:
      case FLUIDMECH_STRESS:
      case FLUIDMECH_PRESSURE:
      case FLUIDMECH_COMP_STRESS:
      case FLUIDMECH_VISC_STRESS:
      case FLUIDMECH_PRES_TENS:
      case FLUIDMECH_TOTAL_STRESS:
      case FLUIDMECH_SURFACE_TRACTION:
      case FLUIDMECH_ZERO_PRESSURE:
      case WATER_PRESSURE:
      case WATER_PRES_TENS:
      case SMOOTH_ZERO_PRESSURE:
      case WATER_SURFACE_TRACTION:
        return "Pa";
        break;

      case FLUIDMECH_PRESSURE_TIME_DERIV_1:
      case ACOU_PRESSURE_DERIV_1:
        return "Pa/s";
        break;

      case ACOU_PRESSURE_DERIV_2:
      case FLUIDMECH_PRESSURE_TIME_DERIV_2:
        return "Pa/s^2";
        break;

      case FLUIDMECH_PRESSURE_DERIV_1:
        return "Pa/m";
        break;

      case FLUIDMECH_PRESSURE_DERIV_2:
        return "Pa/m^2";
        break;

      case ACOU_DIV_LH_TENSOR:
        return "kg m^-2 s^-2";
        break;

      case ACOU_RHS_LOAD:
      case ACOU_RHS_LOADP:
      case ACOU_APE_RHS_LOAD:
      case ACOU_VORTEX_RHS_LOAD:
        return "kg m^-3 s^-2";
        break;

      case ACOU_RHS_LOAD_DENSITY:
        return "kg m^-6 s^-2";
        break;
        
      case ACOU_RHS_LOAD_DENSITY_VECTOR:
        return "kg m^-5 s^-2";
        break;

      case SPLIT_LAMB:
        return "kg/(ms)^2";
        break;

      case SPLIT_DIVLAMB:
        return "kg/(m^3s^2)";
        break;

      case ELEC_CHARGE:
      case ELEC_RHS_LOAD:
        return "C";
        break;

      case ELEC_SURFACE_CHARGE_DENSITY:
        return "C/m^2";
        break;

      case ELEC_POLARIZATION:
      case ELEC_FLUX_DENSITY:
        return "C/m^2";
        break;

      case ELEC_POWER_DENSITY:
      case FLUIDMECH_VISCOUS_DISS_POWER_DENS:
        return "W/m^3";
        break;

      case ELEC_FORCE_DENSITY:
      case MAG_FORCE_MAXWELL_DENSITY:
      case MAG_NORMALFORCE_MAXWELL_DENSITY:
      case MAG_TANGENTIALFORCE_MAXWELL_DENSITY:
      case SMOOTH_CONTACT_FORCE_DENSITY:
        return "N/m^2";

      case DIV_MEAN_FLUIDMECH_VELOCITY:
      case FLUIDMECH_VELOCITY_DERIV_1:
      case FLUIDMECH_VORTICITY:
      case FLUIDMECH_TDR:
      case FLUIDMECH_STRAINRATE:
        return "1/s";
        break;

      case FLUIDMECH_VELOCITY_DERIV_2:
        return "1/s^2";
        break;

      case FLUIDMECH_TEMP:
      case HEAT_TEMPERATURE:
      case HEAT_MEAN_TEMPERATURE:
        return "K";
        break;

      case FLUIDMECH_TKE:
        return "J";
        break;

      case WATER_SURFACE_TORQUE:
        return "Nm";
        break;

      case WATER_TDT:
      case WATER_SURFACE_TORQUE_DENSITY:
        return "N/m";
        break;

      case FLUIDMECH_WVT:
        return "kg m^-2 s^-2";
        break;

      case FLUIDMECH_WVT_DENSITY:
        return "kg m^-2 s^-2 m s^-1";
        break;

      case FLUIDMECH_DENSITY:
      case ELEM_DENSITY:
        return "kg/m^3";
        break;

      case FLUIDMECH_WVT_PHI:
        return "kg m^-3 s^-1";
        break;

      case FLUIDMECH_WVT_DENSITY_PHI:
        return "kg m^-3 s^-1 m s^-1";
        break;

      case FLUIDMECH_WVT_U1:
      case FLUIDMECH_WVT_U2:
        return "m s^-1";
        break;

      case HEAT_TEMPERATURE_D1:
        return "K/s";
        break;

      case HEAT_SOURCE_DENSITY:
        return "W/m^3";
        break;

      case COIL_CURRENT:
      case ELEC_GRAD_V_INT:
      case ELEC_CURRENT:
      case MAG_TOTAL_POTENTIAL:
      case MAG_REDUCED_POTENTIAL:
        return "A";
        break;
        
      case MAG_FIELD_INTENSITY:
      case MAG_AVERAGED_FIELD_INTENSITY:
      case MAG_MAGNETIZATION:
      case MAG_POTENTIAL_GRAD:
        return "A/m";
        break;

      case ELEC_NORMAL_CURRENT_DENSITY:
      case ELEC_CURRENT_DENSITY:
        return "A/m^2";
        break;

      case COIL_CURRENT_DERIV1:
        return "A/s";
        break;

      case MAG_RHS_LOAD:
        return "Am";
        break;

      case MAG_EDDY_CURRENT_DENSITY:
      case MAG_TOTAL_CURRENT_DENSITY:
        return "A/m^2";
        break;

      case MAG_ELEM_RELUCTIVITY:
        return "Am/Vs";

      case ELEC_ELEM_PERMITTIVITY:
        return "(As)/(Vm)";

      case ELEC_POTENTIAL:
      case COIL_VOLTAGE:
        return "V";
        break;

      case ELEC_FIELD_INTENSITY:
      case ELEC_FIELD_INTENSITY_SURF:
      case MAG_POTENTIAL_DERIV1:
      case MAG_POTENTIAL_ADJ_DERIV1:
        return "V/m";
        break;

      case MAG_FLUX:
      case COIL_LINKED_FLUX:
      case COIL_VOLTAGE_INTEGRAL:
        return "Vs";
        break;

      case MAG_POTENTIAL:
        return "Vs/m";
        break;

      case MAG_POTENTIAL_ADJ:
        return "Vs/m";
        break;

      case MAG_FLUX_DENSITY:
      case MAG_AVERAGED_FLUX_DENSITY:
      case MAG_FLUX_DENSITY_SURF:
      case MAG_CURL_ADJ:
      case MAG_NORMAL_FLUX_DENSITY:
        return "Vs/m^2";
        break;

      case ELEC_POTENTIAL_DERIV_1:
        return "V/s";
        break;


      case MAG_ELEM_PERMEABILITY:
        return "Vs/Am";

      case MAG_POLARIZATION:
        return "Vs/m^2";
        

      case MAG_CORE_LOSS_DENSITY:
        return "W/kg";
        break;


      case MECH_STRESS:
      case MECH_IRR_STRESS:
      case MECH_NORMAL_STRESS:
      case MECH_THERMAL_STRESS:
      case MECH_PRINCIPAL_STRESS:
      case MECH_PRINCIPAL_STRESS_MIN:
      case MECH_PRINCIPAL_STRESS_MAX:
      case MECH_PRINCIPAL_STRESS_MED:
      case MECH_PRINCIPAL_STRESS_MIN_SCAL:
      case MECH_PRINCIPAL_STRESS_MED_SCAL:
      case MECH_PRINCIPAL_STRESS_MAX_SCAL:
      return "N/m^2";
        break;

      case ACOU_MIXED_MASS_LOAD:
      case ACOU_MIXED_MOMENTUM_LOAD:
      case ACOU_LAMB_RHS:
        return "-";
        break;

      case HEAT_RHS_LOAD:
        return "?";
        break;

      case ELEC_PSEUDO_POLARIZATION:
      case ELEC_PHYSICAL_PSEUDO_DENSITY:
      case FLUX_INDUCED_STRAIN:
      case LBM_PHYSICAL_PSEUDO_DENSITY:
      case MECH_STRAIN:
      case MECH_IRR_STRAIN:
      case MECH_PRINCIPAL_STRAIN:
      case MECH_PRINCIPAL_STRAIN_MIN:
      case MECH_PRINCIPAL_STRAIN_MAX:
      case MECH_PRINCIPAL_STRAIN_MED:
      case MECH_PRINCIPAL_STRAIN_MIN_SCAL:
      case MECH_PRINCIPAL_STRAIN_MAX_SCAL:
      case MECH_PRINCIPAL_STRAIN_MED_SCAL:
      case MECH_THERMAL_STRAIN:
      case MECH_ELEM_VOL:
      case MECH_ELEM_POROSITY:
      case MECH_PSEUDO_DENSITY:
      case PSEUDO_DENSITY:
      case PHYSICAL_PSEUDO_DENSITY:
        return "";
        break;
      case VOLUME:
        return "m^3";
        break;

      default:
        return "unknown";
        break;
    }
  }

  
  /** the String2Enum/Enum2String stuff is depreciated and shall be realized as Enum<> */
  template <class TYPE>
  void String2Enum(const std::string &in, TYPE &out)
  {
    EXCEPTION("Not implemented");
  }

  template <class TYPE>
  void Enum2String(const TYPE &in, std::string &out)
  {
    EXCEPTION("Not implemented");
  }

  // ComplexFormat
  template<>
  void String2Enum<ComplexFormat>(const std::string &in, ComplexFormat &out) {
   if ( in == "amplPhase" ) {
      out = AMPLITUDE_PHASE;
    } else if (in == "realImag" ) {
      out = REAL_IMAG;
    } else {
     EXCEPTION( "ComplexFormat '" << in << "' not known!" );
   }
  }

  template<>
  void Enum2String<ComplexFormat>(const ComplexFormat &in, std::string &out) {

    switch( in ) {

    case AMPLITUDE_PHASE:
      out = "amplPhase";
      break;

    case REAL_IMAG:
      out = "realImag";
      break;

    default:
      EXCEPTION( "ComplexFormat " << in << " not known!") ;
    }
  }

  // FreqSamplingType
  template<>
  void String2Enum<FreqSamplingType>( const std::string &in,
                                      FreqSamplingType &out ) {

    if ( in == "noSamplingType" ) {
      out = NO_SAMPLING_TYPE;
    }
    else if ( in == "linear" ) {
      out = LINEAR_SAMPLING;
    }
    else if ( in == "log" ) {
      out = LOG_SAMPLING;
    }
    else if ( in == "reverseLog" ) {
      out = REVERSE_LOG_SAMPLING;
    }
    else {
      EXCEPTION( "'" << in << "' cannot be converted into an '"
                 << "FreqSamplingType' item!" );
    }
  }

  template<>
  void Enum2String<FreqSamplingType>( const FreqSamplingType &in,
                                      std::string &out ) {

    switch( in ) {

    case NO_SAMPLING_TYPE:
      out = "noSamplingType";
      break;

    case LINEAR_SAMPLING:
      out = "linear";
      break;

    case LOG_SAMPLING:
      out = "log";
      break;

    case REVERSE_LOG_SAMPLING:
      out = "reverseLog";
      break;

    default:
      EXCEPTION( "Found no conversion for supplied FreqSamplingType value!" );
    }
  }

  template<>
  void Enum2String<SubTensorType>(const SubTensorType &in, std::string &out) {
    switch(in) {
      case NO_TENSOR:
        out = "noTensor";
        break;
    case PLANE_STRAIN:
      out = "planeStrain";
      break;
    case PLANE_STRESS:
      out = "planeStress";
      break;
    case PLANE:
      out = "plane";
      break;
    case AXI:
      out = "axi";
      break;
    case FULL:
      out = "3d";
      break;
    default:
      EXCEPTION("No conversion found for your 'SubTensorType'");
    }
  }

  template<>
  void String2Enum<SubTensorType>( const std::string &in, SubTensorType &out ) {
    if ( in == "noTensor" ) {
       out = NO_TENSOR;
    }
    if ( in == "planeStrain" ) {
      out = PLANE_STRAIN;
    }
    else if ( in == "planeStress" ) {
      out = PLANE_STRESS;
    }
    else if ( in == "plane" ) {
      out = PLANE;
    }
    else if ( in == "axi" ) {
      out = AXI;
    }
    else if ( in == "3d" || in == "2.5d") {
      out = FULL;
    }
    else {
      EXCEPTION("No conversion from string to 'SubTensorType' found" );
    }
  }



  template<>
  void Enum2String<MaterialClass>(const MaterialClass &in,
                                  std::string &out) {
    switch(in) {
      case NO_CLASS:
        out = "No MaterialClass";
        break;
      case TESTMAT:
        out = "testmat";
        break;
      case ELECTROMAGNETIC:
        out = "magnetic";
        break;
      case ELECTROMAGNETIC_DARWIN:
        out = "magnetic";
        break;
      case ELECTROSTATIC:
        out = "electric";
        break;
      case ELECTRICCONDUCTION:
        out = "elecConduction";
        break;
      case ELECQUASISTATIC:
        out = "elecQuasistatic";
        break;
      case ACOUSTIC:
        out = "acoustic";
        break;
      case FLOW:
        out = "flow";
        break;
      case MECHANIC:
        out = "mechanical";
        break;
      case PIEZO:
        out = "piezo";
        break;
      case THERMIC:
        out = "heatConduction";
        break;
      case PYROELECTRIC:
        out = "pyroelectric";
        break;
      case THERMOELASTIC:
        out = "thermoelastic";
        break;
      case MAGNETOSTRICTIVE:
        out = "magnetoStrictive";
        break;
      case SMOOTH:
        out = "smooth";
        break;
      default:
        EXCEPTION("No conversion found for your 'MaterialClass'" );
    }
  }


  template<>
  void String2Enum<MaterialClass>( const std::string &in, MaterialClass &out ) {

    if ( in == "No MaterialClass" ) {
      out = NO_CLASS;
    }
    else if ( in == "testmat" ) {
       out = TESTMAT;
     }
    else if ( in == "electromagnetic" ) {
      out = ELECTROMAGNETIC;
    }
    else if ( in == "electromagnetic_darwin" ) {
	  out = ELECTROMAGNETIC_DARWIN;
	}
    else if ( in == "elecConduction" ) {
      out = ELECTRICCONDUCTION;
    }
    else if ( in == "elecQuasistatic" ) {
          out = ELECQUASISTATIC;
    }
    else if ( in == "electromagnetic" ) {
      out = ELECTROMAGNETIC;
    }
    else if ( in == "fluid" ) {
      out = ACOUSTIC;
    }
    else if ( in == "flow" ) {
      out = FLOW;
    }
    else if ( in == "mechanic" ) {
      out = MECHANIC;
    }
    else if ( in == "piezo" ) {
      out = PIEZO;
    }
    else if ( in == "thermic" ) {
      out = THERMIC;
    }
    else if ( in == "pyroelectric" ) {
      out = PYROELECTRIC;
    }
    else if ( in == "thermoelastic" ) {
      out = THERMOELASTIC;
    }
    else if ( in == "magnetoStrictive" ) {
      out = MAGNETOSTRICTIVE;
    }
    else if ( in == "smooth" ) {
      out = SMOOTH;
    }
    else {
      EXCEPTION( "'" << in << "' cannot be converted into an '"
                 << "MaterialClass' item!" );
    }
  }


  template<>
  void String2Enum<Directions>( const std::string &in, Directions &out ) {

    if ( in == "X" ) {
      out = X;
    }
    else if ( in == "Y" ) {
      out = Y;
    }
    else if ( in == "Z" ) {
      out = Z;
    }
    else if ( in == "radXY" ) {
      out = radXY;
    }
    else if ( in == "radXZ" ) {
      out = radXZ;
    }
    else if ( in == "radYZ" ) {
      out = radYZ;
    }
    else {
      EXCEPTION( "'" << in << "' cannot be converted into an '"
                 << "DataType' item!" );
    }
  }

  template<>
  void Enum2String<Directions>(const Directions &in,
                               std::string &out) {
    switch(in) {
      case X:
        out = "X";
        break;
      case Y:
        out = "Y";
        break;
      case Z:
        out = "Z";
        break;
      case radXY:
        out = "radXY";
        break;
      case radXZ:
        out = "radXZ";
        break;
      case radYZ:
        out = "radYZ";
        break;

      default:
        EXCEPTION("No conversion found for your 'MaterialClass'" );
    }
  }

  // localInversionFlag
  template<>
  void String2Enum<localInversionFlag>( const std::string &in,
                                      localInversionFlag &out ) {

    if ( in == "LevenbergMarquardt" ) {
      out = LOCAL_LEVENBERGMARQUARDT;
    }
    else if ( in == "Newton" ) {
      out = LOCAL_NEWTON;
    }
    else if ( in == "JacobiFreeNewton" ) {
      out = LOCAL_JACOBIFREENEWTON;
    }
    else if ( in == "ProjectedLevenbergMarquardt" ) {
      out = LOCAL_PROJECTEDLM;
    }
    else if ( in == "EverettBased-scalarOnly" ) {
      out = LOCAL_EVERETTBASED;
    }
    else if ( in == "Fixpoint" ) {
      out = LOCAL_FIXPOINT;
    }
    else if ( in == "NoLocalInversion" ){
      out = LOCAL_NOINVERSION;
    }    
    else {
      EXCEPTION( "'" << in << "' cannot be converted into an '"
                 << "localInversionFlag' item!" );
    }
  }

  template<>
  void Enum2String<localInversionFlag>( const localInversionFlag &in,
                                      std::string &out ) {

    switch( in ) {

    case LOCAL_LEVENBERGMARQUARDT:
      out = "LevenbergMarquardt";
      break;

    case LOCAL_NEWTON:
      out = "Newton";
      break;

    case LOCAL_JACOBIFREENEWTON:
      out = "JacobiFreeNewton";
      break;

    case LOCAL_PROJECTEDLM:
      out = "ProjectedLevenbergMarquardt";
      break;
      
    case LOCAL_EVERETTBASED:
      out = "EverettBased-scalarOnly";
      break;

    case LOCAL_FIXPOINT:
      out = "Fixpoint";
      break;
    
    case LOCAL_NOINVERSION:
      out = "NoLocalInversion";
      break;
      
    default:
      EXCEPTION( "Found no conversion for supplied localInversionFlag value!" );
    }
  }
  
  template<>
  void String2Enum<NonLinType>( const std::string &in, NonLinType &out )
  {
    if ( in == "noNonLinearity" || in == "" ) {
      out = NO_NONLINEARITY;
    } else if( in == "westervelt" ) {
      out = WESTERVELT;;
    } else if( in == "kuznetsov") {
      out = KUZNETSOV;
    } else if( in == "variableSOS_CN1") {
      out = VARIABLE_SOS_CN1;
    } else if( in == "variableSOS_CN2") {
      out = VARIABLE_SOS_CN2;
    } else if( in == "variableSOS_CN2Mean") {
      out = VARIABLE_SOS_CN2Mean;
    } else if( in == "material") {
      out = MATERIAL;
    } else if( in == "geometric") {
      out = GEOMETRIC;
    } else if( in == "hysteresis") {
      out = HYSTERESIS;
    } else if( in == "piezoMicroHF") {
      out = PIEZO_MICRO_HF;
    } else if( in == "permeability") {
      out = PERMEABILITY;
     } else if( in == "reluctivity") {
      out = RELUCTIVITY; 
    } else if( in == "reluctivity_magstrict"){
      out = RELUCTIVITY_MAGSTRICT;
    } else if( in == "heatConductivity") {
      out = NLHEAT_CONDUCTIVITY;
    } else if( in == "heatCapacity") {
      out = NLHEAT_CAPACITY;
    } else if( in == "density") {
      out = NLHEAT_DENSITY;
    } else if( in == "heatSource") {
      out = NLHEAT_SOURCE;
    } else if( in == "thermalRadiation") {
      out = THERMAL_RADIATION;
    } else if( (in == "elecConductivity") || (in == "electricConductivity") ) {
       out = NLELEC_CONDUCTIVITY;
    } else if( in == "elecBiPole") {
       out = NLELEC_BIPOLE;
    } else if( in == "elecBiPoleTempDep") {
       out = NLELEC_BIPOLE_TEMP_DEP;
    } else if( in == "elecTriPole") {
       out = NLELEC_TRIPOLE;
    } else if( in == "elecTriPoleTempDep") {
       out = NLELEC_TRIPOLE_TEMP_DEP;
    } else if( in == "elecPermittivity") {
       out = NLELEC_PERMITTIVITY;
    } else {
      EXCEPTION( "'" << in << "' cannot be converted into an "
                 << "'NonLinType' item!" );
    }
  }


  template<>
  void Enum2String<NonLinType>( const NonLinType &in, std::string& out ) {

    switch(in) {
      case NO_NONLINEARITY:
        out = "noNonLinearity";
        break;
      case WESTERVELT:
        out = "westervelt";
        break;
      case KUZNETSOV:
        out = "kuznetsov";
        break;
      case VARIABLE_SOS_CN1:
        out = "variableSOS_CN1";
        break;
      case VARIABLE_SOS_CN2:
        out = "variableSOS_CN2";
        break;
      case VARIABLE_SOS_CN2Mean:
        out = "variableSOS_CN2Mean";
        break;
      case MATERIAL:
        out = "material";
        break;
      case GEOMETRIC:
        out = "geometric";
        break;
      case HYSTERESIS:
        out = "hysteresis";
        break;
      case NLELEC_CONDUCTIVITY:
        out = "elecConductivity";
        break;
      case NLELEC_BIPOLE:
        out = "elecBiPole";
        break;
      case NLELEC_BIPOLE_TEMP_DEP:
        out = "elecBiPoleTempDep";
        break;
      case NLELEC_TRIPOLE:
        out = "elecTriPole";
        break;
      case NLELEC_TRIPOLE_TEMP_DEP:
        out = "elecTriPoleTempDep";
        break;
      case PIEZO_MICRO_HF:
        out = "piezoMicroHF";
        break;
      case PERMEABILITY:
        out = "permeability";
        break;
      case RELUCTIVITY:
        out = "reluctivity";
        break;  
      case RELUCTIVITY_MAGSTRICT:
      out = "reluctivity_magstrict";
      break;
      case NLHEAT_CONDUCTIVITY:
        out = "heatConductivity";
        break;
      case NLHEAT_CAPACITY:
        out = "density";
        break;
      case NLHEAT_DENSITY:
        out = "heat";
        break;
      case NLHEAT_SOURCE:
        out = "heatSource";
        break;
      case THERMAL_RADIATION:
        out = "thermalRadiation";
        break;
      case NLELEC_PERMITTIVITY:
        out = "elecPermittivity";
        break;
      default:
        EXCEPTION( "No conversion found for 'NonLinType' " << in );
    }
  }

  template <>
  void String2Enum<DampingType>(const std::string &in, DampingType &out)
  {
    if (in == "none") {
      out = NONE;
    }
    else if (in == "rayleigh") {
      out = RAYLEIGH;
    }
    else if (in == "abc") {
      out = ABCDAMP;
    }
    else if (in == "thermoViscous") {
      out = THERMOVISCOUS;
    }
    else if (in == "fractional") {
      out = FRACTIONAL;
    }
    else if (in == "fractiona_gl") {
      out = FRACTIONAL_GL;
    }
    else if (in == "fractional_blank") {
      out = FRACTIONAL_BLANK;
    }
    else if (in == "fractional_gl_int") {
      out = FRACTIONAL_GL_INT;
    }
    else if (in == "fractional_blank_int") {
      out = FRACTIONAL_BLANK_INT;
    }
    else if (in == "pml") {
      out = PML;
    }
    else if (in == "dampLayer") {
      out = DAMPLAYER;
    }
    else if (in == "mapping") {
      out = MAPPING;
    }
    else if (in == "adaptedLossTangensDelta") {
      out = ADAPTED_LOSS_TANGENS_DELTA;
    }
    else if (in == "globalRayleigh") {
      out = GLOBAL_RAYLEIGH;
    }
    else {
      EXCEPTION("'" << in << "' cannot be converted into a 'DampingType' item!");
    }
  }

  template <>
  void Enum2String<DampingType>(const DampingType &in, std::string &out)
  {
    switch (in) {
    case NONE:
      out = "none";
      break;
    case RAYLEIGH:
      out = "rayleigh";
      break;
    case ABCDAMP:
      out = "abc";
      break;
    case THERMOVISCOUS:
      out = "thermoViscous";
      break;
    case FRACTIONAL:
      out = "fractional";
      break;
    case FRACTIONAL_GL:
      out = "fractiona_gl";
      break;
    case FRACTIONAL_BLANK:
      out = "fractional_blank";
      break;
    case FRACTIONAL_GL_INT:
      out = "fractional_gl_int";
      break;
    case FRACTIONAL_BLANK_INT:
      out = "fractional_blank_int";
      break;
    case PML:
      out = "pml";
      break;
    case DAMPLAYER:
      out = "dampLayer";
      break;
    case MAPPING:
      out = "mapping";
      break;
    case ADAPTED_LOSS_TANGENS_DELTA:
      out = "adaptedLossTangensDelta";
      break;
    case GLOBAL_RAYLEIGH:
      out = "globalRayleigh";
      break;
    default:
      EXCEPTION("No conversion found for 'DapmingType' " << in);
    }
  }

  // ****************************************************************
  // ****************************************************************
  //              THE OLAS ENVIRONMENT STARTS HERE
  // ****************************************************************
  // ****************************************************************

  // Specialisation for StopCritType
  template<>
  void Enum2String<StopCritType>(const StopCritType &in,
                                 std::string &out) {
    switch( in ) {
      case NOSTOPCRITTYPE:
        out = "no stopping criterion";
        break;
      case ABSNORM:
        out = "absNorm";
        break;
      case RELNORM_RHS:
        out = "relNormRHS";
        break;
      case RELNORM_RES0:
        out = "relNormRes0";
        break;

      default:
        EXCEPTION( "No string value found for the specified value of the "
            << "enumeration datatype StopCritType.\n"
            << "Seems to indicate a missing case implementation!" );
    }
  }

  // Specialisation for AMG interpolation type
  template<>
  void Enum2String<AMGInterpolationType>(const AMGInterpolationType &in,
                                         std::string &out) {
    switch( in ) {
      case AMG_INTERPOLATION_CONSTANT:
        out = "constant";
        break;
      case AMG_INTERPOLATION_SIMPLE_WEIGHTED:
        out = "simpleWeighted";
        break;
      case AMG_INTERPOLATION_SMOOTHED_SCALING:
        out = "smoothedScaling";
        break;
      case AMG_INTERPOLATION_DEVELOP:
        out = "develop";
        break;
      default:
        EXCEPTION( "No string value found for the specified value of the "
            << "enumeration datatype AMGInterpolationType.\n"
            << "Seems to indicate a missing case implementation!" );
    }
  }

  // specialisation for AMG smoother type
  template<>
  void Enum2String<AMGSmootherType>(const AMGSmootherType &in,
                                    std::string &out) {
    switch( in ) {
      case AMG_SMOOTHER_GAUSSSEIDEL:
        out = "GaussSeidel";
        break;
      case AMG_SMOOTHER_DAMPED_JACOBI:
        out = "Jacobi";
        break;
      default:
        EXCEPTION( "No string value found for the specified value of the "
            << "enumeration datatype AMGSmootherType.\n"
            << "Seems to indicate a missing case implementation!"; );
    }
  }

  // Specialisation for FEMatrixType
  template<>
  void Enum2String<FEMatrixType>(const FEMatrixType &in,
                                 std::string &out) {
    switch( in ) {
      case NOTYPE:
        out = "no_fe_matrix";
        break;
      case SYSTEM:
        out = "system";
        break;
      case STIFFNESS:
        out = "stiffness";
        break;
      case DAMPING:
        out = "damping";
        break;
      case DAMPING_AUX:
        out = "damping_aux";
        break;
      case CONVECTION:
        out = "convection";
        break;
      case MASS:
        out = "mass";
        break;
      case AUXILIARY:
        out = "auxiliary";
        break;
      case STIFFNESS_UPDATE:
        out = "stiffness_update";
        break;
      case DAMPING_UPDATE:
        out = "damping_update";
        break;
      case MASS_UPDATE:
        out = "MASS_update";
        break;
      case SYSTEM_HYSTFREE:
        out = "system_hystfree";
        break;
      case SYSTEM_FIXPOINT:
        out = "system_fixpoint";
        break;
      case SYSTEM_FD_JACOBIAN:
        out = "system_fd_Jacobian";
        break;
      case SYSTEM_DELTAMAT_JACOBIAN:
        out = "system_deltaMat_Jacobian";
        break;
      case BACKUP:
        out = "backup";
        break;
      case GEOMETRIC_STIFFNESS:
        out = "geometric_stiffness";
        break;

      default:
        EXCEPTION( "No string value found for the specified value of the "
            << "enumeration datatype FEMatrixTypeType.\n"
            << "Seems to indicate a missing case implementation!" );
    }
  }

  // Specialisation for IDBCType
  template<>
  void Enum2String<IDBCType>(const IDBCType &in,
                             std::string &out) {
    switch( in ) {

      case IDBC_NOTYPE:
        out = "notype";
        break;
      case IDBC_ELIMINATION:
        out = "elimination";
        break;
      case IDBC_PENALTY:
        out = "penalty";
        break;
      default:
        EXCEPTION( "No string value found for the specified value of the "
            << "enumeration datatype IDBCType.\n"
            << "Seems to indicate a missing case implementation!" );
    }
  }

  // -----------------------------------------------
  //  Implementation of string conversion routines
  // -----------------------------------------------

  // Specialisation for StopCritType
  template<>
  void String2Enum<StopCritType>( const std::string &in, StopCritType &out ) {

    if ( in == "noStopCritType" ) {
      out = NOSTOPCRITTYPE;
    }
    else if ( in == "absNorm" ) {
      out = ABSNORM;
    }
    else if ( in == "relNormRHS" ) {
      out = RELNORM_RHS;
    }
    else if ( in == "relNormRes0" ) {
      out = RELNORM_RES0;
    }
    else {
      EXCEPTION( "No enumeration value found in StopCritType for '"
          << in << "'\n A missing case implementation?" );
    }
  }


  // map specialisation for interpolation types
  template<>
  void String2Enum<AMGInterpolationType>( const std::string &in,
                                          AMGInterpolationType &out ) {
    if ( in == "constant" ) {
      out = AMG_INTERPOLATION_CONSTANT;
    } else if( in == "simpleWeighted" ) {
      out = AMG_INTERPOLATION_SIMPLE_WEIGHTED;
    } else if( in == "develop" ) {
      out = AMG_INTERPOLATION_DEVELOP;
    } else {
      EXCEPTION( "No enumeration value found in AMGInterpolationType "
          << "for '" << in << "'\n A missing case implementation?" );
    }
  }

  // map specialisation for interpolation types
  template <>
  void String2Enum<AMGSmootherType>( const std::string &in,
                                     AMGSmootherType   &out ) {
    if ( in == "GaussSeidel" ) {
      out = AMG_SMOOTHER_GAUSSSEIDEL;
    } else if( in == "Jacobi" ) {
      out = AMG_SMOOTHER_DAMPED_JACOBI;
    } else {
      EXCEPTION( "No enumeration value found in AMGSmoothertype "
          << "for '" << in << "'\n A missing case implementation?" );
    }
  }

  // Specialisation for FEMatrixType
  template<>
  void String2Enum<FEMatrixType>( const std::string &in,
                                  FEMatrixType &out ) {

    if ( in == "nomatrixtype" )
      out = NOTYPE;
    else if ( in == "system" )
      out = SYSTEM;
    else if ( in == "stiffness" )
      out = STIFFNESS;
    else if ( in == "damping" )
      out = DAMPING;
    else if ( in == "convection" )
      out = CONVECTION;
    else if ( in == "mass" )
      out = MASS;
    else if ( in == "mass_update" )
      out = MASS_UPDATE;
    else if ( in == "stiffness_update" )
      out = STIFFNESS_UPDATE;
    else if ( in == "damping_update" )
      out = DAMPING_UPDATE;
    else if ( in == "system_hystfree")
      out = SYSTEM_HYSTFREE;
    else if ( in == "system_fixpoint")
      out = SYSTEM_FIXPOINT;
    else if ( in == "system_fd_Jacobian")
      out = SYSTEM_FD_JACOBIAN;
    else if ( in == "system_deltaMat_Jacobian")
      out = SYSTEM_DELTAMAT_JACOBIAN;
    else if ( in == "backup")
      out = BACKUP;
    else if ( in == "geometric_stiffness")
          out = GEOMETRIC_STIFFNESS;
    else {
      EXCEPTION( "String '" << in << "' cannot be converted to item of "
                 << "'FEMatrixType'!" );
    }
  }

  // Specialisation for IDBCType
  template<>
  void String2Enum<IDBCType>( const std::string &in, IDBCType &out ) {

    if ( in == "notype" )
      out = IDBC_NOTYPE;
    else if ( in == "elimination" )
      out = IDBC_ELIMINATION;
    else if ( in == "penalty" )
      out = IDBC_PENALTY;
    else {
      EXCEPTION( "String '" << in << "' cannot be converted to item of "
                 << "'IDBCType'!" );
    }
  }


  const char* GetBlasThreadsEnvVariable(bool full)
  {
// see configure_def_headers.cmake and def_use_blas_hh.in
#if defined(USE_MKL)
  return full ? "MKL_NUM_THREADS" : "MKL";
#elif defined(USE_ACCELERATE)
  return full ? "VECLIB_MAXIMUM_THREADS" : "VECLIB";
#else
  return full ? "DUMMY_NUM_THREADS" : "DUMMY";
#endif
  }

  bool HasBlasThreadsEnvVariable()
  {
#if defined(USE_MKL) || defined(USE_ACCELERATE)
  return true;
#else
  return false;
#endif
  }

  // on windows we have no setenv() but putenv() which demands a non-local string space
  #if defined(USE_OPENMP) && defined(WIN32)
    static char omp_threads[1024];
    static char other_threads[1024]; // for Windows on MKL_NUM_THREADS
  #endif

  void SetNumberOfThreads(int cfs, bool homogenize, bool quiet)
  {
    using std::cout;

    CFS_NUM_THREADS = cfs; // this is a global int variable

    // on almost all Linux and Windows systems this is MKL_NUM_THREADS.
    const char* otherenv = GetBlasThreadsEnvVariable(true); // DUMMY_NUM_THREADS for openblas and netlib
    const std::string otherstr = std::string(otherenv);

#ifdef USE_OPENMP
    assert(cfs >= 1); // programOptions has this as fallback
    // our variables are cfs (command line or CFS_NUM_THREADS) and omp and mkl from environment
    int omp = getenv("OMP_NUM_THREADS") != NULL ? atoi(getenv("OMP_NUM_THREADS")) : -1;

    // here mostly for MKL_NUM_THREADS or VECLIB_MAXIMUM_THREADS or not used (openblas and netlib) and only DUMMY_NUM_THREADS is used
    // we set for openblas the dummy as we already have OMP_NUM_THREADS and don't want to double
    int other = getenv(otherenv) != NULL ? atoi(getenv(otherenv)) : -1;

    // set the threads via environment variables, this is how e.g. the external libs (lis, cholmod, mkl) read it.
    // note that this change is only for this process and child processes, it does not change the system settings.
    if(omp <= 0) {
      omp = cfs; // for later comparison
      omp_set_num_threads(omp); // e.g. for ginkgo we cannot rely on the environment variable
      #ifndef WIN32
        setenv("OMP_NUM_THREADS",to_string(omp).c_str(),1); // libs read it there
      #else
        // there is no setenv() on Windows and ancient POSIX putenv needs a static string source
        strncpy(omp_threads, string("OMP_NUM_THREADS=" + to_string(omp)).c_str(), 1022);
        putenv(omp_threads);
      #endif
    }
    if(other <= 0) {
      other = cfs;
      #ifndef WIN32
        setenv(otherenv,to_string(other).c_str(),1); // 1 = overwrite but it shall not have been there before
      #else
        strncpy(other_threads, string(otherstr + "=" + to_string(other)).c_str(), 1022);
        putenv(other_threads);
     #endif
    }

    // now some hints if OMP and MKL/ACCELERATE were set to something different than CFS
    std::string msg = "Hint: openCFS threads (CFS_NUM_THREADS=" + to_string(cfs) + ") differs from ";
    unsigned int org = msg.size();
    if(cfs != omp && cfs != other)
      msg+= "OMP_NUM_THREADS=" + to_string(omp) + " and " + otherstr + "=" + to_string(other);
    if(cfs != omp && cfs == other)
      msg+= "OMP_NUM_THREADS=" + to_string(omp);
    if(cfs == omp && cfs != other)
      msg+= otherstr + "=" + to_string(other);
    if(msg.size() != org && !quiet)
      cout << ">> " << msg << std::endl;
#endif
  }

} // end of namespace




