#ifndef ENVIRONMENT_TPYES_HH
#define ENVIRONMENT_TPYES_HH

/** To help avoiding cyclic inclusion, we strictly keep here:
 *  - typedef (enums)
 *  - defines which cannot be avoided
 *  - strictly no include, also not from std stuff
 *  - type which need includes (e.g. shared_ptr) need to go to Environment.hh
 *  - class definitions are ok
 * Many enums have a corresponding Enum<TYPE> object, this is to be defined in EnvironmentEnum.hh
 */

#define REFACTOR WARN("Commented out due to refactoring");

namespace CoupledField {

  // forward class declaration
  class Domain;
  class CFSMessenger;
  class LogConfigurator;
  
  // Type definitions for regions
  typedef int RegionIdType;
  static const RegionIdType NO_REGION_ID = -1;
  static const RegionIdType ALL_REGIONS = -2;

  static const double NORM_EPS = 1e-6;  // needed e.g. for lower bounds of norms in iteration loops
  static const double EPS = 1e-12;     // value for absolute precision (needed e.g. for lower bounds of norms in iteration loops)

  // Abbreviation for global handle
  static const int MathParser_GLOB_HANDLER = 0;


  //! number of threads for parallel CFS loops
  extern unsigned int CFS_NUM_THREADS;

  /** options for StdVector/Vector/Matrix::ToString().
   * TS_MATLAB and TS_PYTHON allow copy and paste to matlab/python
   * See the ::ToString() for behavior for the format */
  typedef enum {TS_PLAIN, TS_MATLAB, TS_PYTHON, TS_INFO, TS_NONZEROS} ToStringFormat;

  //! specifications of Lapack routines for different types of system matrices in
  //! matrix.solveWithLapack
  //! Z - Complex valued matrix
  //! SY - symmetric matrix
  //! GE - general matrix
  //! HE - hermitian matrix
  typedef enum {ZSYSVMTX, ZGESVMTX, ZHESVMTX} lapackSysMatType;


  //! logging configurator
  extern LogConfigurator* logConf;
  
  //! Global pointer to domain object
  extern Domain* domain;

  //! Damping type
  enum DampingType{NONE=0, RAYLEIGH=1, ABCDAMP=2, THERMOVISCOUS=3,
    FRACTIONAL=4, FRACTIONAL_GL=5, FRACTIONAL_BLANK=6,
    FRACTIONAL_GL_INT=7, FRACTIONAL_BLANK_INT=8,
    PML = 9, DAMPLAYER = 10, MAPPING = 11, 
    ADAPTED_LOSS_TANGENS_DELTA = 12, GLOBAL_RAYLEIGH = 13};

  //! Type of nonlinearity for certain pdes
  typedef enum { NO_NONLINEARITY, WESTERVELT, KUZNETSOV, VARIABLE_SOS_CN1,
    VARIABLE_SOS_CN2, VARIABLE_SOS_CN2Mean,
    MATERIAL, GEOMETRIC, HYSTERESIS, PIEZO_MICRO_HF, PERMEABILITY, RELUCTIVITY_MAGSTRICT,
    NLHEAT_CONDUCTIVITY, NLHEAT_CAPACITY, NLHEAT_DENSITY, NLHEAT_SOURCE, THERMAL_RADIATION, NLELEC_CONDUCTIVITY, /* temperature dep */
    NLELEC_BIPOLE, NLELEC_BIPOLE_TEMP_DEP, NLELEC_TRIPOLE, NLELEC_TRIPOLE_TEMP_DEP, NLELEC_PERMITTIVITY} NonLinType;

  //! Terminal Connector enum
  typedef enum { CATHODE, ANODE, DRAIN, SOURCE, GATE, BODY, COLLECTOR,
    EMITTER, BASE} TerminalConnector;

  //! Describes all possible solution types in a CFS simulation
  typedef enum{
    NO_SOLUTION_TYPE, INVALID_SOLUTION_TYPE,
    // ============
    //  MECHANICAL
    // ============
    // -- primary results --
    MECH_DISPLACEMENT, MECH_ACCELERATION, MECH_VELOCITY,MECH_RHS_LOAD,

    // --- flux / derived quantities --
    MECH_STRESS, MECH_PRINCIPAL_STRESS,
    MECH_PRINCIPAL_STRESS_MIN, MECH_PRINCIPAL_STRESS_MAX, MECH_PRINCIPAL_STRESS_MED,
    MECH_PRINCIPAL_STRESS_MIN_SCAL, MECH_PRINCIPAL_STRESS_MAX_SCAL, MECH_PRINCIPAL_STRESS_MED_SCAL,
    MECH_STRAIN, MECH_PRINCIPAL_STRAIN,
    MECH_PRINCIPAL_STRAIN_MIN, MECH_PRINCIPAL_STRAIN_MAX, MECH_PRINCIPAL_STRAIN_MED,
    MECH_PRINCIPAL_STRAIN_MIN_SCAL, MECH_PRINCIPAL_STRAIN_MAX_SCAL, MECH_PRINCIPAL_STRAIN_MED_SCAL,
    MECH_THERMAL_STRAIN, MECH_STRUCT_INTENSTIY,
    MECH_NORMAL_STRUCT_INTENSITY, VON_MISES_STRESS,
    VON_MISES_STRAIN, MECH_KIN_ENERGY_DENS, MECH_DEFORM_ENERGY_DENS,
    MECH_TOTAL_ENERGY_DENS, MECH_NORMAL_STRESS, MECH_THERMAL_STRESS,
    MECH_DYADIC_STRAIN, MECH_QUAD_DISP,
    MECH_NORMAL_DISPLACEMENT, MECH_NORMAL_VELOCITY,
    // -- integrated quantities --
    MECH_KIN_ENERGY, MECH_DEFORM_ENERGY, MECH_TOTAL_ENERGY,
    MECH_POWER, MECH_DEF_SURF_VOLUME, MECH_WEIGHT, MECH_DYADIC_STRAIN_SUM, MECH_QUAD_DISP_SUM,
    // -- coupling quantities --
    MECH_FORCE, MECH_VELOCITY_ELEM,

    // -- optimization properties
    MECH_PSEUDO_DENSITY, // this is historical. We shall switch to PSEUDO_DENSITY as in mag
    PSEUDO_DENSITY, // this shall be used for all physical - is the design variable in topology optimization
    PHYSICAL_PSEUDO_DENSITY, // the filtered and penalized PSEUDO_DENSITY as used for the simulation
    MECH_SHAPE, MECH_TENSOR_TRACE, MECH_TENSOR, MECH_TENSOR_HILL_MANDEL,MECH_ELEM_VOL, MECH_ELEM_POROSITY,

    // ===============
    //  ELECTROSTATIC / ELECTRIC CURRENT
    // ===============
    ELEC_POTENTIAL, ELEC_POTENTIAL_DERIV1, ELEC_FIELD_INTENSITY, ELEC_FIELD_INTENSITY_SURF,
    ELEC_FORCE , ELEC_FORCE_VWP, ELEC_FORCE_DENSITY,
    ELEC_CHARGE, ELEC_SURFACE_CHARGE_DENSITY, ELEC_FLUX_DENSITY, ELEC_RHS_LOAD,
    ELEC_ENERGY, ELEC_ENERGY_DENSITY, ELEC_POLARIZATION, ELEC_ELEM_PERMITTIVITY,
    ELEC_CURRENT_DENSITY, ELEC_FIELD_INTENSITY_TRANSVERSAL, ELEC_FIELD_INTENSITY_LONGITUDINAL,
	  ELEC_GRAD_V_INT, ELEC_NORMAL_CURRENT_DENSITY, ELEC_CURRENT, ELEC_CURRENT_SURF, 
    ELEC_CURRENT_FIELD_INTENSITY, ELEC_POWER, ELEC_POWER_DENSITY,
	  DISPLACEMENT_CURRENT_FIELD_INTENSITY, DISPLACEMENT_CURRENT_SURF, ELEC_COND_TENSOR,
    DISPLACEMENT_NORMAL_CURRENT_DENSITY, ELECTRIC_AND_DISPLACEMENT_CURRENT_DENSITY,
    ELECTRIC_AND_DISPLACEMENT_NORMAL_CURRENT_DENSITY, DISPLACEMENT_CURRENT, 
    ELEC_AND_DISPLACEMENT_CURRENT, 

    // -- optimization properties
    ELEC_PSEUDO_POLARIZATION, ELEC_PHYSICAL_PSEUDO_DENSITY, ELEC_TENSOR_TRACE, ELEC_TENSOR,

    // =================
    //  SMOOTH
    // =================
    SMOOTH_DISPLACEMENT, SMOOTH_VELOCITY, SMOOTH_ZERO_PRESSURE, SMOOTH_STRAIN,
    SMOOTH_CONTACT_FORCE_DENSITY, SMOOTH_CONTACT_FORCE, SMOOTH_DEFORM_ENERGY_DENS, SMOOTH_DEFORM_ENERGY,

    // ==========
    //  ACOUSTIC
    // ==========
    ACOU_POTENTIAL, ACOU_PRESSURE, ACOU_SURFPRESSURE, ACOU_PRESSURE_DERIV_1,
    ACOU_PRESSURE_DERIV_2,ACOU_VELOCITY,ACOU_POSITION, ACOU_NORMAL_VELOCITY,ACOU_ACCELERATION, ACOU_FORCE,
    ACOU_POWERDENSITY, ACOU_POWER, ACOU_INTENSITY, ACOU_NORMAL_INTENSITY,
    ACOU_POTENTIAL_DERIV_1, ACOU_POTENTIAL_DERIV_2, ACOU_RHS_LOAD, ACOU_RHS_LOAD_DENSITY, ACOU_RHS_LOAD_DENSITY_VECTOR,ACOU_RHS_LOAD_DENSITY_VECTOR_IN_NORMAL,
    ACOU_RHSVAL, ACOUSURF_RHSVAL, ACOU_ENERGY, ACOU_POWER_PLANEWAVE,ACOU_NORMAL_INTENSITY_PLANEWAVE,
    ACOU_POT_ENERGY, ACOU_KIN_ENERGY,
    ACOU_SURFINTENSITY, ACOU_DIV_LH_TENSOR,
    ACOU_PMLAUXVEC, ACOU_PMLAUXSCALAR,
    ACOU_ELEM_SPEED_OF_SOUND, ACOU_IMPEDANCE,

    // ==========
    //  SPLITTING
    // ==========
    SPLIT_SCALAR, SPLIT_VECTOR, SPLIT_RHS_LOAD, SPLIT_SCALAR_VELOCITY, SPLIT_VECTOR_VELOCITY,
    SPLIT_POT_ENERGY, SPLIT_LAMB, SPLIT_DIVLAMB,

    // ==========
    //  WATER WAVES
    // ==========
    WATER_PRESSURE, WATER_PMLAUXVEC, WATER_PMLAUXSCALAR, WATER_RHS_LOAD,WATER_POSITION, WATER_PRES_TENS,
    WATER_SURFACE_TRACTION, WATER_SURFACE_FORCE, WATER_SURFACE_TORQUE, WATER_TDT, WATER_SURFACE_TORQUE_DENSITY,

    // =========
    // AEROACOUSTIC SOURCE TERMS
    // =========
    ACOU_RHS_LOADP, ACOU_APE_RHS_LOAD, ACOU_VORTEX_RHS_LOAD,

    // ==========
    //  ACOUSTIC MIXED
    // ==========
    ACOU_MIXED_MASS_LOAD, ACOU_MIXED_MOMENTUM_LOAD, ACOU_LAMB_RHS,

    // ==========
    //  MAGNETIC
    // ==========
    // -- primary results --
    MAG_POTENTIAL, MAG_POTENTIAL_DERIV1, MAG_TOTAL_POTENTIAL, MAG_REDUCED_POTENTIAL, MAG_RHS_LOAD,

    // --- flux / derived quantities --
    MAG_FLUX_DENSITY, MAG_FLUX, MAG_NORMAL_FLUX_DENSITY, MAG_FLUX_DENSITY_SURF, MAG_FIELD_INTENSITY, MAG_EDDY_CURRENT_DENSITY,
    MAG_COIL_CURRENT_DENSITY, MAG_TOTAL_CURRENT_DENSITY, MAG_POTENTIAL_DIV, MAG_FORCE_LORENTZ_DENSITY,
    MAG_FORCE_LORENTZ_DENSITY_STATIC, MAG_FORCE_LORENTZ_DENSITY_HARMONIC, 

    MAG_FORCE_MAXWELL_DENSITY, MAG_NORMALFORCE_MAXWELL_DENSITY, MAG_TANGENTIALFORCE_MAXWELL_DENSITY,
    MAG_EDDY_POWER_DENSITY, MAG_ENERGY_DENSITY, MAG_CORE_LOSS_DENSITY,
    MAG_JOULE_LOSS_POWER_DENSITY, MAG_JOULE_LOSS_POWER_DENSITY_ON_NODES, MAG_JOULE_LOSS_POWER,ELEC_POTENTIAL_DERIV_1,

    // -- integrated quantities --
    MAG_FORCE_VWP, MAG_FORCE_LORENTZ, MAG_FORCE_MAXWELL,
    MAG_ENERGY, MAG_EDDY_POWER, MAG_EDDY_CURRENT, MAG_EDDY_CURRENT1, MAG_EDDY_CURRENT2, MAG_TOTAL_CURRENT, MAG_CORE_LOSS,


    // -- coil quantities --
    COIL_INDUCTANCE, COIL_INDUCED_VOLTAGE, COIL_CURRENT, COIL_VOLTAGE, COIL_LINKED_FLUX, COIL_CURRENT_DERIV1,

    // -- material related results --
    MAG_ELEM_PERMEABILITY, MAG_MAGNETIZATION,FLUX_INDUCED_STRAIN,MAG_POLARIZATION,MAG_ELEM_RELUCTIVITY,

    // -- magnetic topology optimization
    // for pseudo density see mech
    RHS_PSEUDO_DENSITY, PHYSICAL_RHS_PSEUDO_DENSITY,


    // =================
    //  HEAT CONDUCTION
    // =================
    HEAT_TEMPERATURE, HEAT_TEMPERATURE_D1, HEAT_RHS_LOAD, HEAT_SOURCE_DENSITY, HEAT_FLUX_DENSITY, HEAT_FLUX_INTENSITY, HEAT_FLUX,
    HEAT_MEAN_TEMPERATURE, HEAT_CONDUCTIVITY_TENSOR_HOM,
    // ===============
    //  FLUIDMECHANIC
    // ===============
    FLUID_FORCE,
    LAMBDA_K, MEAN_FLUIDMECH_VELOCITY, FLUIDMECH_VELOCITY, FLUIDMECH_MESH_VELOCITY,
    FLUIDMECH_MESH_VELOCITY_NODE, FLUIDMECH_MESH_VELOCITY_ELEM, FLUIDMECH_TOTAL_VELOCITY_ELEM,
    FLUIDMECH_NORMAL_VELOCITY,
    FLUIDMECH_PRESSURE, FLUIDMECH_ZERO_PRESSURE,
    DIV_MEAN_FLUIDMECH_VELOCITY,MEAN_FLUIDMECH_VELOCITY_NORMAL,
    FLUIDMECH_VELOCITY_DERIV_1, FLUIDMECH_PRESSURE_DERIV_1, FLUIDMECH_PRESSURE_TIME_DERIV_1,
    FLUIDMECH_VELOCITY_DERIV_2, FLUIDMECH_PRESSURE_DERIV_2, FLUIDMECH_PRESSURE_TIME_DERIV_2,
    FLUIDMECH_DENSITY, FLUIDMECH_FORCE, FLUIDMECH_TEMP,
    FLUIDMECH_TKE, FLUIDMECH_TDR, FLUIDMECH_TEF,
    FLUIDMECH_STRESS, FLUIDMECH_COMP_STRESS, FLUIDMECH_VISC_STRESS, FLUIDMECH_PRES_TENS, FLUIDMECH_TOTAL_STRESS,
    FLUIDMECH_SURFACE_TRACTION,
    FLUIDMECH_STRAINRATE,
    FLUIDMECH_ENERGY, FLUIDMECH_STABILPARAM,
    FLUIDMECH_WVT, FLUIDMECH_WVT_DENSITY, FLUIDMECH_VORTICITY,
    FLUIDMECH_WVT_PHI, FLUIDMECH_WVT_DENSITY_PHI,
    FLUIDMECH_WVT_U1, FLUIDMECH_WVT_U2, FLUIDMECH_WVT_F,
    FLUIDMECH_VISCOUS_DISS_POWER_DENS, FLUIDMECH_VISCOUS_DISS_POWER_DENS_DIV, FLUIDMECH_VISCOUS_DISS_POWER_DENS_STRAIN,
    FLUIDMECH_VISCOUS_DISS_POWER,
    FLUIDMECH_INTENSITY, FLUIDMECH_SURFINTENSITY, FLUIDMECH_NORMAL_INTENSITY, FLUIDMECH_POWER,
    FLUIDMECH_INTENSITY_PRESSURE_ONLY, FLUIDMECH_SURFINTENSITY_PRESSURE_ONLY, FLUIDMECH_NORMAL_INTENSITY_PRESSURE_ONLY, 
    FLUIDMECH_POWER_PRESSURE_ONLY, 
    FLUIDMECH_SURFIMPEDANCE, FLUIDMECH_IMPEDANCE,


    // ===============
    //  TESTPDE
    // ===============
    TEST_DOF, TEST_FIELD, TEST_RHS_LOAD,

    // ======
    //  MISC
    // ======
    HOMOGENIZED_TENSOR,
    OPT_RESULT_1, OPT_RESULT_2, OPT_RESULT_3, OPT_RESULT_4, OPT_RESULT_5, OPT_RESULT_6,
    OPT_RESULT_7, OPT_RESULT_8, OPT_RESULT_9, OPT_RESULT_10, OPT_RESULT_11, OPT_RESULT_12,
    OPT_RESULT_13, OPT_RESULT_14, OPT_RESULT_15, OPT_RESULT_16, OPT_RESULT_17, OPT_RESULT_18,
    OPT_RESULT_19, OPT_RESULT_20,OPT_RESULT_21,OPT_RESULT_22,OPT_RESULT_23,OPT_RESULT_24,OPT_RESULT_25,
    OPT_RESULT_26,OPT_RESULT_27,OPT_RESULT_28,OPT_RESULT_29,OPT_RESULT_30,OPT_RESULT_31, OPT_RESULT_32, OPT_RESULT_33,
    OPT_RESULT_34, OPT_RESULT_35, OPT_RESULT_36, OPT_RESULT_37, OPT_RESULT_38, OPT_RESULT_39,
    OPT_RESULT_40, OPT_RESULT_41, OPT_RESULT_42, OPT_RESULT_43, OPT_RESULT_44, OPT_RESULT_45,
    OPT_RESULT_46, OPT_RESULT_47, OPT_RESULT_48, OPT_RESULT_49, OPT_RESULT_50, OPT_RESULT_51,
    OPT_RESULT_52, OPT_RESULT_53, OPT_RESULT_54, OPT_RESULT_55, OPT_RESULT_56, OPT_RESULT_57,
    OPT_RESULT_58, OPT_RESULT_59, OPT_RESULT_60, OPT_RESULT_61, OPT_RESULT_62, OPT_RESULT_63,
    OPT_RESULT_64, OPT_RESULT_65, OPT_RESULT_66,

    LAGRANGE_MULT, LAGRANGE_MULT_DERIV_1, LAGRANGE_MULT_DERIV_2, LAGRANGE_MULT_1,
    THERMOMECH_FORCE, THERMOELEC_FORCE,
    GRAD_ACOU_SOLUTION, GRAD_ELEC_POTENTIAL, GRAD_ELEC_POTENTIAL_DERIV1,
    ELEM_DENSITY,

    // ======
    // PML
    // ======
    PML_DAMP_FACTOR, PML_TENSOR, PML_DETERMINANT, PML_DISTANCE,

    // ==========
    //  GEOMETRY
    // ==========
    ELEM_LOC_DIR, JACOBIAN, ASPECT_RATIO, VOLUME, NODE_NORMAL,COOSY_X,COOSY_Y,SURFACE_NORMAL,AREA, ETA,XI,

    // === NonFEM LBM results ===
    LBM_PROBABILITY_DISTRIBUTION, LBM_VELOCITY, LBM_DENSITY,
    LBM_PRESSURE, LBM_PHYSICAL_PSEUDO_DENSITY,

    // for hysteresis
    MECH_IRR_STRESS, MECH_IRR_STRAIN,

    // ==========
    //  GENERIC
    // ==========
    // these types can be used by generic postproc functions
    // don't change the naming convention (_X has to be the last part of the name, where X is any integer number)
    // if you need something else, adapt the code so that "resNr" is calculated correctly!!
    GENERIC_RESULT_0, GENERIC_RESULT_1, GENERIC_RESULT_2, GENERIC_RESULT_3, GENERIC_RESULT_4,
    GENERIC_RESULT_5, GENERIC_RESULT_6, GENERIC_RESULT_7, GENERIC_RESULT_8, GENERIC_RESULT_9,

  } SolutionType;

  //! describes the possible material types
  /* Note: For general tensor properties (3x3 tensor), we have the following
   * naming convention:
   *
   *   - PROPERTY_TENSOR:  Denotes the full 3x3 tensor which will be always
   *                       available, independent of the definition type
   *                       (isotropic, orthotropic, transversly isotropic, tensor)
   *   - PROPERTY_SCALAR:  The scalar property in case of isotropic or
   *                       transversly isotropic (in the 1-2-plane) definition.
   *   - PROPERTY_1/_2/_3: Denotes the orthotropic scalar values in the 3 main
   *                       directions. In the transversly isotropic case, only
   *                       the _3 component is used.
   */
  typedef enum{
    NO_MATERIAL,

    // ==========
    //  ACOUSTIC
    // ==========
    ACOU_BULK_MODULUS,
    ACOU_SOUND_SPEED, ACOU_IMPEDANCE_REAL_VAL, ACOU_IMPEDANCE_IMAG_VAL, ACOU_BOVERA,

    // These seem to be unused
    /*ACOU_ALPHA, IMP_HOLE_DIAM, IMP_PLATE_THICKNESS, FLOW_MACH_NUMBER,
    IMP_END_CORRECTION, POROSITY, BETA, MPP_VOLUME_DEPTH, FRACTIONAL_ALG,
    FRACTIONAL_MEMORY, FRACTIONAL_INTERPOL, FRACTIONAL_EXPONENT,*/

    // ===============
    //  ELECTROSTATIC
    // ===============
    ELEC_PERMITTIVITY_TENSOR, ELEC_PERMITTIVITY_SCALAR,
    ELEC_PERMITTIVITY_1, ELEC_PERMITTIVITY_2, ELEC_PERMITTIVITY_3,

    //Parameter for Hysteresis Model
    ELEC_PS_JILES, ELEC_ALPHA_JILES, ELEC_A_JILES, ELEC_K_JILES, ELEC_C_JILES,

    // =====================
    //  ELECTRIC CONDUCTION
    // =====================
    ELEC_CONDUCTIVITY_TENSOR, ELEC_CONDUCTIVITY_SCALAR,
    ELEC_CONDUCTIVITY_1, ELEC_CONDUCTIVITY_2, ELEC_CONDUCTIVITY_3,

    // =======
    //  FLUID
    // =======
    FLUID_ADIABATIC_EXPONENT, FLUID_DYNAMIC_VISCOSITY, FLUID_KINEMATIC_VISCOSITY,
    FLUID_BULK_VISCOSITY, FLUID_BULK_MODULUS,

    // =================
    //  HEAT CONDUCTION
    // =================
    HEAT_CONDUCTIVITY_TENSOR, HEAT_CONDUCTIVITY_SCALAR,
    HEAT_CONDUCTIVITY_1, HEAT_CONDUCTIVITY_2, HEAT_CONDUCTIVITY_3,
    HEAT_CAPACITY, HEAT_REF_TEMPERATURE,

    // ==========
    //  MAGNETIC
    // ==========
    MAG_PERMEABILITY_TENSOR, MAG_PERMEABILITY_SCALAR,
    MAG_PERMEABILITY_1, MAG_PERMEABILITY_2, MAG_PERMEABILITY_3,
    MAG_RELUCTIVITY_TENSOR, MAG_RELUCTIVITY_SCALAR, MAG_RELUCTIVITY_DERIV,
    MAG_CONDUCTIVITY_TENSOR, MAG_CONDUCTIVITY_SCALAR, MAG_PERMITTIVITY_SCALAR, MAG_PERMITTIVITY_TENSOR,
	  MAG_PERMITTIVITY_1, MAG_PERMITTIVITY_2, MAG_PERMITTIVITY_3,
    MAG_CONDUCTIVITY_1, MAG_CONDUCTIVITY_2, MAG_CONDUCTIVITY_3,
    //MAGNETIZATION,
    MAG_CORE_LOSS_PER_MASS,
    MAG_BH_VALUES, MAG_BH_VALUES_1, MAG_BH_VALUES_2, MAG_BH_VALUES_3,
    MAG_BH_DATA_ACCURACY, MAG_BH_MAX_APPROX_VAL,
    // -- Magnetic EB Hysteresis Parameters
    MAG_PS_EB, MAG_A_EB, MAG_MU0_EB, MAG_NUMS_EB, MAG_CHI_FACTOR_EB,
    
    // ============
    //  MECHANICAL
    // ============
    MECH_STIFFNESS_TENSOR, MECH_KMODULUS, MECH_LAME_MU, MECH_LAME_LAMBDA,
    MECH_EMODULUS, MECH_EMODULUS_1, MECH_EMODULUS_2, MECH_EMODULUS_3,
    MECH_POISSON, MECH_POISSON_3, MECH_POISSON_12, MECH_POISSON_23, MECH_POISSON_13,
    MECH_GMODULUS, MECH_GMODULUS_3, MECH_GMODULUS_23, MECH_GMODULUS_13, MECH_GMODULUS_12,

    // -- Viscoelasticity --
    MECH_BULK_RELAX_TIMES, MECH_BULK_RELAX_MODULI, MECH_VISCO_BULK_INITIAL,
    MECH_SHEAR_RELAX_TIMES, MECH_SHEAR_RELAX_MODULI, MECH_VISCO_SHEAR_INITIAL,
    MECH_VISCO_STIFFNESS_TENSOR, MECH_VISCO_LONGTERM_TENSOR,

    // -- Thermal expansion --
    MECH_THERMAL_EXPANSION_TENSOR, MECH_THERMAL_EXPANSION_SCALAR,
    MECH_THERMAL_EXPANSION_1, MECH_THERMAL_EXPANSION_2, MECH_THERMAL_EXPANSION_3,
    MECH_TE_REFTEMPERATURE,

    // Seems to be unused
    //PYROCOEFFICIENT_TENSOR,

    // ==========
    //  TEST PDE
    // ==========
    TEST_ALPHA, TEST_BETA,

    // ===============
    //  PIEZOELECTRIC
    // ===============
    PIEZO_TENSOR,

    MEAN_TEMPERATURE, SPON_POLARIZATION, SPON_STRAIN, DRIVING_FORCE_90,
    DRIVING_FORCE_180, RATE_CONSTANT, VISCO_PLASTIC_INDEX, SATURATION_INDEX,
    VOLUME_FRAC_INIT, EFIELD0, STRESS0, DCOUPLE0,
    SCALE_FORCE_ELEC, SCALE_FORCE_MECH, SCALE_FORCE_COUPLE, CORE_LOSS,

    // prescribed magnetization for permanent magnets
    PRESCRIBED_MAGNETIZATION,PRESCRIBED_MAGNETIZATION_X,
    PRESCRIBED_MAGNETIZATION_Y,PRESCRIBED_MAGNETIZATION_Z,

    // ==================
    //  MAGNETOSTRICTION
    // ==================
    MAGSTRICT_RELUCTIVITY, MAGNETOSTRICTION_TENSOR_h,
    MAGNETOSTRICTION_TENSOR_h_mech, MAGNETOSTRICTION_TENSOR_h_mag,
    A_JILES, ALPHA_JILES, K_JILES, C_JILES, JILES_TEST, Y_REMANENCE,

    // =========
    //  GENERAL
    // =========
    DENSITY,

    // =========
    // Rayleigh damping
    // =========
    RAYLEIGH_ALPHA, RAYLEIGH_BETA, LOSS_TANGENS_DELTA,

    // -- General Material Nonlinearity --
    NONLIN_DEPENDENCY,
    // These seem to be unused
    //NONLIN_COEFFICIENT, NONLIN_APPROXIMATION_TYPE,
    //NONLIN_DATA_NAME, DATA_ACCURACY, MAX_APPROX_VAL,

    /*
     * ENUM for Hysteresis with Preisach like models
     * New: in order to allow for a separate modeling of irreversible, hysteretic
     * strain via a second hyst operator, most enoms have to be defined twice
     * (i.e. one for hyst operator 1 and one for the additional hyst operator).
     * To avoid double code (due to different enum names during set/get functions)
     * divide the followin enums in two blocks that are separated by a fixed offset
     * (e.g. 100). Doing so, it should be equivalent to write
     * material->GetScalar(paramSet.muDat_A_, PREISACH_WEIGHTS_MUDAT_A + 100, Global::REAL);
     * instead of
     * material->GetScalar(paramSet.muDat_A_, PREISACH_WEIGHTS_MUDAT_A_STRAIN, Global::REAL);
     */
    /*
     * polarization; start at 300
     */
    X_SATURATION=300, Y_SATURATION,
    PREISACH_WEIGHTS_TYPE, PREISACH_WEIGHTS, PREISACH_WEIGHTS_DIM,
    PREISACH_WEIGHTS_CONSTVALUE,
    PREISACH_WEIGHTS_MUDAT_A, PREISACH_WEIGHTS_MUDAT_H, PREISACH_WEIGHTS_MUDAT_H2,
    PREISACH_WEIGHTS_MUDAT_SIGMA, PREISACH_WEIGHTS_MUDAT_SIGMA2, PREISACH_WEIGHTS_MUDAT_ETA,
    PREISACH_WEIGHTS_MUDAT_PARAMSFORHALFRANGE,
    PREISACH_WEIGHTS_TENSOR,
    PREISACH_WEIGHTS_ANHYST_A, PREISACH_WEIGHTS_ANHYST_B, PREISACH_WEIGHTS_ANHYST_C,
    PREISACH_WEIGHTS_ANHYST_ONLY,
    PREISACH_WEIGHTS_ANHYST_D,PREISACH_WEIGHTS_ANHYST_PARAMSFORHALFRANGE,
    PREISACH_WEIGHTS_ANHYSTCOUNTINGTOOUTPUTSAT,
    PREISACH_WEIGHTS_FOR_MAYERGOYZ_VECTOR,
    P_DIRECTION, P_DIRECTION_X, P_DIRECTION_Y, P_DIRECTION_Z,
    HYST_MODEL, HYSTERESIS_DIM,
    ROT_RESISTANCE, EVAL_VERSION, ANG_DISTANCE, ANG_CLIPPING, ANG_RESOLUTION, AMP_RESOLUTION,
    SCALETOSAT,
    PREISACH_MAYERGOYZ_NUM_DIR,PREISACH_MAYERGOYZ_ISOTROPIC,PREISACH_MAYERGOYZ_CLIPOUTPUT,
    MAYERGOYZ_STARTAXIS_X,MAYERGOYZ_STARTAXIS_Y,MAYERGOYZ_STARTAXIS_Z,MAYERGOYZ_LOSSPARAM_A,MAYERGOYZ_LOSSPARAM_B,
    MAYERGOYZ_USEABSDPHI,MAYERGOYZ_NORMALIZEXINEXP,MAYERGOYZ_RESTRICTIONOFPSI,MAYERGOYZ_SCALINGOFXINEXP,
    HYST_TYPE_IS_PREISACH,
    /*
     * strain; offset 100
     */
    X_SATURATION_STRAIN=400, Y_SATURATION_STRAIN,
    PREISACH_WEIGHTS_TYPE_STRAIN, PREISACH_WEIGHTS_STRAIN, PREISACH_WEIGHTS_DIM_STRAIN,
    PREISACH_WEIGHTS_CONSTVALUE_STRAIN,
    PREISACH_WEIGHTS_MUDAT_A_STRAIN, PREISACH_WEIGHTS_MUDAT_H_STRAIN, PREISACH_WEIGHTS_MUDAT_H2_STRAIN,
    PREISACH_WEIGHTS_MUDAT_SIGMA_STRAIN, PREISACH_WEIGHTS_MUDAT_SIGMA2_STRAIN, PREISACH_WEIGHTS_MUDAT_ETA_STRAIN,
    PREISACH_WEIGHTS_MUDAT_PARAMSFORHALFRANGE_STRAIN,
    PREISACH_WEIGHTS_TENSOR_STRAIN,
    PREISACH_WEIGHTS_ANHYST_A_STRAIN, PREISACH_WEIGHTS_ANHYST_B_STRAIN, PREISACH_WEIGHTS_ANHYST_C_STRAIN,
    PREISACH_WEIGHTS_ANHYST_ONLY_STRAIN,
    PREISACH_WEIGHTS_ANHYST_D_STRAIN,PREISACH_WEIGHTS_ANHYST_PARAMSFORHALFRANGE_STRAIN,
    PREISACH_WEIGHTS_ANHYSTCOUNTINGTOOUTPUTSAT_STRAIN,
    PREISACH_WEIGHTS_FOR_MAYERGOYZ_VECTOR_STRAIN,
    S_DIRECTION, S_DIRECTION_X, S_DIRECTION_Y, S_DIRECTION_Z,
    HYST_MODEL_STRAIN, HYSTERESIS_DIM_STRAIN,
    ROT_RESISTANCE_STRAIN, EVAL_VERSION_STRAIN, ANG_DISTANCE_STRAIN, ANG_CLIPPING_STRAIN, ANG_RESOLUTION_STRAIN, AMP_RESOLUTION_STRAIN,
    SCALETOSAT_STRAIN,
    PREISACH_MAYERGOYZ_NUM_DIR_STRAIN,PREISACH_MAYERGOYZ_ISOTROPIC_STRAIN,PREISACH_MAYERGOYZ_CLIPOUTPUT_STRAIN,
    MAYERGOYZ_STARTAXIS_X_STRAIN,MAYERGOYZ_STARTAXIS_Y_STRAIN,MAYERGOYZ_STARTAXIS_Z_STRAIN,
    MAYERGOYZ_LOSSPARAM_A_STRAIN,MAYERGOYZ_LOSSPARAM_B_STRAIN,
    MAYERGOYZ_USEABSDPHI_STRAIN,MAYERGOYZ_NORMALIZEXINEXP_STRAIN,MAYERGOYZ_RESTRICTIONOFPSI_STRAIN,MAYERGOYZ_SCALINGOFXINEXP_STRAIN,
    HYST_TYPE_IS_PREISACH_STRAIN,
    /*
     * additional parameter that are not double
     */
    // strain related
    S_SATURATION, IRRSTRAIN_REUSE_P, HYST_COUPLING_DEFINED,
    HYST_STRAIN_FORM,
    HYST_IRRSTRAINS, HYST_IRRSTRAIN_C1, HYST_IRRSTRAIN_C2, HYST_IRRSTRAIN_C3,
    HYST_IRRSTRAIN_CI, HYST_IRRSTRAIN_CI_SIZE, HYST_IRRSTRAIN_D0, HYST_IRRSTRAIN_D1,
    HYST_IRRSTRAIN_SCALETOSAT, HYST_IRRSTRAIN_PARAMSFORHALFRANGE,
    // initial state
    PREISACH_SCALEINITIALSTATE, PREISACH_PRESCRIBEOUTPUT,
    INITIAL_STATE, INITIAL_STATE_X, INITIAL_STATE_Y, INITIAL_STATE_Z,
    // model related
    PRINT_PREISACH, PRINT_PREISACH_RESOLUTION, IS_TESTING,
    // Parameter for local inversion of hyst operator
    MAX_NUM_IT_HYST_INV, RES_TOL_H_HYST_INV, RES_TOL_B_HYST_INV, RES_TOL_H_HYST_INV_ISREL, RES_TOL_B_HYST_INV_ISREL,VEC_HYST_INV_METHOD,
    MAX_NUM_REG_IT_HYST_INV, ALPHA_REG_HYST_INV, ALPHA_REG_MIN_HYST_INV, ALPHA_REG_MAX_HYST_INV,
    TRUST_LOW_HYST_INV, TRUST_MID_HYST_INV, TRUST_HIGH_HYST_INV,
    JAC_RESOLUTION_HYST_INV, JAC_IMPLEMENTATION_HYST_INV, MAX_NUM_LS_IT_HYST_INV, ALPHA_LS_MIN_HYST_INV, ALPHA_LS_MAX_HYST_INV,
    STOP_INV_LS_AT_LOCAL_MIN, HYST_INV_PROJLM_MU, HYST_INV_PROJLM_RHO, HYST_INV_PROJLM_BETA, HYST_INV_PROJLM_SIGMA, HYST_INV_PROJLM_GAMMA,
    HYST_INV_PROJLM_TAU, HYST_INV_PROJLM_C, HYST_INV_PROJLM_P, HYST_INV_FP_SAFETYFACTOR, HYST_LOCAL_INVERSION_PRINT_WARNINGS,
    // For tracing of hyst operator (debugging and finetuning parameters)
    TRACE_FORCE_CENTRALDIFF,TRACE_FORCE_RETRACING,TRACE_JAC_RESOLUTION,

    // ============
    //  SMOOTH
    // ============
    SMOOTH_STIFFNESS_TENSOR, SMOOTH_KMODULUS, 
    SMOOTH_LAME_MU, SMOOTH_LAME_LAMBDA,
    SMOOTH_EMODULUS, 
    SMOOTH_EMODULUS_1, SMOOTH_EMODULUS_2, SMOOTH_EMODULUS_3,
    SMOOTH_POISSON,
    SMOOTH_POISSON_3, SMOOTH_POISSON_12, SMOOTH_POISSON_23, SMOOTH_POISSON_13,
    SMOOTH_GMODULUS,
    SMOOTH_GMODULUS_3, SMOOTH_GMODULUS_23, SMOOTH_GMODULUS_13, SMOOTH_GMODULUS_12
    
  } MaterialType;
    
  typedef enum{NO_TENSOR = -1, FULL, PLANE_STRAIN, PLANE_STRESS, PLANE, AXI} SubTensorType;

  //! Denote the available material classes
  typedef enum{
    NO_CLASS, TESTMAT, ELECTROMAGNETIC, ELECTROMAGNETIC_DARWIN, ELECTROSTATIC, ACOUSTIC, FLOW, MECHANIC,
    PIEZO, THERMIC, PYROELECTRIC, THERMOELASTIC, MAGNETOSTRICTIVE, 
    ELECTRICCONDUCTION, ELECQUASISTATIC, SMOOTH
  } MaterialClass;

  //! material tensor notation
  typedef enum { NO_NOTATION, VOIGT, HILL_MANDEL } MaterialTensorNotation;

  // type of approximation / interpolation
  typedef enum{ NO_APPROX_TYPE, LIN_INTERPOLATE, BILIN_INTERPOLATE, TRILIN_INTERPOLATE, CUBIC_SPLINES, SMOOTH_SPLINES, ANALYTIC } ApproxCurveType;

  //! material parameter to be approximated / interpolated
  typedef enum{ GENERIC, MAGNETIC_MAT_BH, HEAT_MAT_CONDUCTIVITY, HEAT_MAT_CAPACITY, ELEC_MAT_CONDUCTIVITY }  ApproxMaterialCurves;

  //! Enumeration for directions
  //! direction of various fields
  //! "Rad" means readial, following two letters indicate the stress-plane
  enum  Directions {X=0, Y=1, Z=2, radXY=3, radXZ=4, radYZ=5};

  //! nonlinear method definition
  typedef enum {FIXEDPOINT=1, NEWTON=2, HYST_DEBUG=3, HYST_FIXPOINT_IT=4, HYST_FIXPOINT_TS=5, HYST_DELTAMAT_IT=6, HYST_DELTAMAT_TS=7 } NonLinMethodType ;
  
  //! output format for complex numbers
  typedef enum {REAL_IMAG, AMPLITUDE_PHASE} ComplexFormat;

  //! Data type for specification of frequency sampling approach

  //! This enumeration data type is used for distinguishing the different
  //! approaches for sampling the frequency domain when performing a
  //! harmonic analysis with openCFS. The type currently allows the following
  //! values
  //! - NO_SAMPLING_TYPE
  //! - LINEAR_SAMPLING
  //! - LOG_SAMPLING
  //! - REVERSE_LOG_SAMPLING
  typedef enum { NO_SAMPLING_TYPE, LINEAR_SAMPLING, LOG_SAMPLING,
                 REVERSE_LOG_SAMPLING } FreqSamplingType;

  // ****************************************************************
  // ****************************************************************
  //              THE OLAS ENVIRONMENT STARTS HERE
  // ****************************************************************
  // ****************************************************************

  //! Standard type for identifiying a PDE in the algebraic system
  //! \note This type is guaranteed to be of an Integer or
  //! enumeration type, so it can be used as array index
  typedef int FeFctIdType;
  
  //! Standard default identifier
  static const FeFctIdType NO_FCT_ID = -1;

  //! Standard default identifier
  static const FeFctIdType PSEUDO_FCT_ID = -999;

  //! Standard type for identifiying a sparsity pattern of a matrix object
  //! in combination with the PatternPool class.
  typedef unsigned int PatternIdType;
  //! Value to allow consistent storage of unavailable pattern identifier
  static const PatternIdType NO_PATTERN_ID = 0;


  // ---- Enumerations for general algebraic system information ----
  //! Type of finite element matrix

  //! This enumeration data type describes the type of FE-matrices which
  //! are used to assemble an effective system matrix.
  //! The enumeration contains the following values:
  //! - NOTYPE
  //! - SYSTEM
  //! - STIFFNESS
  //! - DAMPING
  //! - CONVECTION
  //! - MASS
  //! - AUXILIARY - e.g. for radiation optimization used in SurfaceNormalInt
  //! - STIFFNESS_UPDATE - utilized if only some parts needs reassemble
  //! - DAMPING_UPDATE - utilized if only some parts needs reassemble
  //! - MASS_UPDATE - utilized if only some parts needs reassemble
  //! For hysteresis, several versions of the system matrices are utlized:
  //! - SYSTEM_HYSTFREE - effectuve system matrix without consider hysteresis; may contain other non-linearities, though
  //!                     this matrix is also utilized for FIXPOINT scheme
  //! - SYSTEM_FD_JACOBIAN - effective system matrix representing an approximation to the Jacobian of the residual
  //!                    (i.e. K_FD_Jac = K_eff + df_hyst/du), using a finite difference approximation of df_hyst/du
  //! - SYSTEM_DELTAMAT_JACOBIAN - as FD_JACOBIAN, but df_hyst/du is approximated using the deltaMat approach
  //!                    (i.e. no full FD_JACOBIAN but difference quotient between current and previous values)
  //! - GEOMETRIC_STIFFNESS - used in Buckling Analysis (also called stress stiffness matrix)
  //! - BACKUP - currently unused? Will be now used to store a copy of the system matrix when we reduce the number of non-zero elements

  typedef enum {NOTYPE, SYSTEM = 1, STIFFNESS, DAMPING, DAMPING_AUX, CONVECTION, MASS, AUXILIARY, STIFFNESS_UPDATE, DAMPING_UPDATE , MASS_UPDATE,
    SYSTEM_HYSTFREE, SYSTEM_FIXPOINT, SYSTEM_FD_JACOBIAN, SYSTEM_DELTAMAT_JACOBIAN, BACKUP, GEOMETRIC_STIFFNESS}
  FEMatrixType;

  //! Maximal number of different FE matrix types

  //! This macro specifies the maximal number
  //! of different types of Finite-Element matrices that can occur in a
  //! openCFS run, i.e. the number of possible choices from FEMatrixType
  //! besides the NOTYPE value.
  extern unsigned int MAX_NUM_FE_MATRICES;

  //! Type of algebraic system

  //! This enumeration data type describes the type of the algebraic system.
  //! To be more precise is describes the type of linear system objects that
  //! are assembled and managed by the linear system. These can either be
  //! StdMatrix / StdVector objects or SBM_Matrix / SBM_Vector objects. Thus,
  //! the current value are:
  //! - NOALGSYSTYPE
  //! - SBM_SYSTEM
  //! - STANDARD_SYSTEM
  typedef enum { NOALGSYSTYPE, SBM_SYSTEM, STANDARD_SYSTEM } AlgSysType;

  //! Type of Eigenvalue Problem
  typedef enum {NO_TYPE, REAL_SYMMETRIC, REAL_GENERAL, COMPLEX_SYMMETRIC, COMPLEX_HERMITIAN, COMPLEX_GENERAL} EigenValueProblemType;

  //! Type of cycle used for algebraic multigrid preconditioner

  //! This enumeration data type describes the type of solutions
  //! cycle which is used by the AMG(algebraic multigrid) preconditioner.
  //! It can take one of the following values:
  //! - NOCYCLE
  //! - VCYCLE
  //! - WCYCLE
  typedef enum {NOCYCLE, VCYCLE, WCYCLE} CycleType;

  //! Switch, needed for the specialized AMG-methods

  //! This enumeration data type describes the version of the special AMG-precond/solver
  //! - SCALAR: one dof per node
  //! - VECTORIAL: 2d or 3d dof's per node
  //! - EDGE: dof's defined on edges (HCurl elements)
  typedef enum { SCALAR = 0, VECTORIAL = 1, EDGE = 2} AMGType;


  //! Type of stopping criterion used by iterative solvers

  //! This enumeration data type describes the different stopping criteria
  //! supported by the iterative solvers. It can take one of the following
  //! values
  //! - NOSTOPCRIT
  //! - ABSNORM
  //! - RELNORM_RHS
  //! - RELNORM_RES0
  typedef enum { NOSTOPCRITTYPE, ABSNORM, RELNORM_RHS,
                 RELNORM_RES0 } StopCritType;

  //! Parallel type of the matrix.

  //! This enumeration data type describes the parallelisation type of the
  //! matrix.
  //! <center>
  //! <table border="1" cellpadding="10">
  //! <tr>
  //!   <td align="center"><b>Enumeration value</b></td>
  //!   <td align="center"><b>Description</b></td>
  //! </tr>
  //! <tr><td>NOPARMATRIXTYPE</td><td>No type was specified.</td></tr>
  //! <tr><td>SEQMATRIX</td><td>This is a sequential matrix object.</td></tr>
  //! <tr><td>PARMATRIX</td><td>This is a ParMatrix object. The type of the
  //! submatrices used for storing the internal and coupling parts can be
  //! obtained from MatrixStorageType.</td></tr>
  //! <tr><td>PARCRSMATRIX</td><td>This is a ParCRS_Matrix object.</td></tr>
  //! </table>
  //! </center>
  typedef enum { NOPARMATRIXTYPE, SEQMATRIX, PARMATRIX, PARCRSMATRIX } ParMatrixType;

  //! \name enumeration types for (parallel) AMG
  //@{
  //! interpolation type
  typedef enum { AMG_INTERPOLATION_CONSTANT,
                 AMG_INTERPOLATION_SIMPLE_WEIGHTED,
                 AMG_INTERPOLATION_SMOOTHED_SCALING,
                 AMG_INTERPOLATION_DEVELOP
               } AMGInterpolationType;
  //! smoother type
  typedef enum { AMG_SMOOTHER_GAUSSSEIDEL,
                 AMG_SMOOTHER_DAMPED_JACOBI
               } AMGSmootherType;
  //@}

  //! Type of treatment of inhomogeneous Dirichlet boundary conditions

  //! Enumeration data type describing the approach used for the treatment
  //! of inhomogeneous Dirichlet boundary conditions. Can currently take one
  //! of the following values
  //! - IDBC_NOTYPE
  //! - IDBC_ELIMINATION
  //! - IDBC_PENALTY
  typedef enum { IDBC_NOTYPE, IDBC_ELIMINATION, IDBC_PENALTY } IDBCType;

//  //! Enums for Hysteresis (solveStepHyst, coefFunctionHyst, Hysteresis)
  typedef enum { SOLVE_NOTSET = -1, SOLVE_GLOBAL_FIXPOINT_B = 0, SOLVE_LOCAL_FIXPOINT_B = 1, SOLVE_JACOBIAN = 2, SOLVE_BGM = 3, SOLVE_DELTAMAT_IT = 4,
    SOLVE_DELTAMAT_TS = 5, SOLVE_DELTAMAT_TS_TOWARDS_IT = 6, SOLVE_BGM_ENHANCED = 7, SOLVE_GLOBAL_FIXPOINT_H = 8, SOLVE_LOCAL_FIXPOINT_B_v2 = 9, SOLVE_CHORD = 10,
    SOLVE_LOCAL_FIXPOINT_H = 11} solveFlag;
  typedef enum { NOTSET = -1, HYSTFREE = 0, FDJACOBIAN = 1, DELTAMATJACOBIAN_LAST_TS = 2, DELTAMATJACOBIAN_LAST_IT = 3, REUSED = 4,
    HYSTFREE_GLOBALLY_SCALED_B_VERSION = 5, HYSTFREE_LOCALLY_SCALED_B_VERSION = 6, HYSTFREE_SCALED_H_VERSION = 7,HYSTFREE_LOCALLY_SCALED_H_VERSION = 8} matrixFlag;
  typedef enum { NOTSET_VEC = -2, PREVIOUS_TS_RES = -1, CURRENT_RES = 0, PREVIOUS_IT_RES = 1, RES_REM = 2 } vectorFlag;
  typedef enum { INVALID_LINESEARCH = -2, NO_LINESEARCH = -1, LS_TRIAL_AND_ERROR = 0, LS_GOLDENSECTION = 1, LS_QUADINTERPOL = 2, LS_VECBASEDPOLY = 3,
    LS_THREEPOINTPOLY = 4, LS_BACKTRACKING_ARMIJO = 5, LS_BACKTRACKING_ARMIJO_NONMONOTONIC = 6, LS_BACKTRACKING_GRADFREE = 7 } linesearchFlag;
  typedef enum { NOFP = 0, FP_GLOBAL_B = 1, FP_LOCAL_B_DIAG = 2, FP_LOCAL_B_FULLJAC = 21, FP_GLOBAL_H = 3, FP_LOCAL_H_DIAG = 4, FP_LOCAL_H_FULLJAC = 41 } fixpointFlag;

    //! local inversion methods for hysteresis operator
  typedef enum {LOCAL_LEVENBERGMARQUARDT=0, LOCAL_NEWTON=1, LOCAL_JACOBIFREENEWTON=2, LOCAL_PROJECTEDLM=3, LOCAL_EVERETTBASED=4, LOCAL_FIXPOINT=5, LOCAL_NOINVERSION=6, LOCAL_NOTIMPLEMENTED=10} localInversionFlag;
  
} // end of namespace

#endif // ENVIRONMENT_TPYES_HH
