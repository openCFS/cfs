// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>
#include <string>

#include "environment.hh"
#include "Utils/tools.hh"
#include "Domain/domain.hh"
#include "DataInOut/Logging/cfslog.hh"


// Since OLAS uses a separate namespace for
// writing out data, two different declarations
// have to be made
namespace OutInfo{

  std::ostream *debug    = NULL;
  std::ostream *cla      = NULL;
  std::ostream *memtrace = NULL;

  // Generate string stream for generation of error messages
//  std::stringstream *error = new std::stringstream();

  // Generate string stream for generation of warning messages
  std::stringstream *warning = new std::stringstream();

#ifdef MEMTRACE
  double sumdmem = 0;
  double sumimem = 0;
#endif
}

namespace CoupledField {


  LogConfigurator * logConf = new LogConfigurator();

#ifdef USE_SCRIPTING
  CFSMessenger * messenger = NULL;
#endif

  Domain * domain = NULL;


  BaseFE * ptQ1     = NULL;
  BaseFE * ptQ2     = NULL;
  BaseFE * ptQ9     = NULL;
  BaseFE * ptTet1   = NULL;
  BaseFE * ptTet2   = NULL;
  BaseFE * ptL1     = NULL;
  BaseFE * ptL2     = NULL;
  BaseFE * ptTr1    = NULL;
  BaseFE * ptTr2    = NULL;
  BaseFE * ptHexa1  = NULL;
  BaseFE * ptHexa2  = NULL;
  BaseFE * ptHexa27 = NULL;
  BaseFE * ptPyra1  = NULL;
  BaseFE * ptPyra2  = NULL;
  BaseFE * ptWedge1 = NULL;
  BaseFE * ptWedge2 = NULL;

  // Initialisation of some global pointers
  WriteInfo *Info = NULL;

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


  // --------------------------------------------
  //  Implementation of enum conversion routines
  // --------------------------------------------
  // General template function
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

  // CouplingInputType
  template<>
  void String2Enum<CouplingInputType>( const std::string &in,
                                       CouplingInputType &out ) {

    if (in == "Coordinate-Displacement")
      out = COORD;
    else if (in == "RHS")
      out = RHS;
    else if (in == "DirichletInhom")
      out = ID_BC;
    else if (in == "materialParam")
      out = MAT;
    else if (in == "gridVel")
      out = GRID_VEL;
    else {
      EXCEPTION( "'" << in << "' cannot be converted into item of "
                 << "'CouplingInputType'!" );
    }
  }

  template<>
  void Enum2String<CouplingInputType>(const CouplingInputType &in,
                                      std::string &out) {

    switch(in) {

    case COORD:
      out = "Coordinate-Displacement";
      break;
    case RHS:
      out = "RHS";
      break;
    case ID_BC:
      out = "DirichletInhom";
      break;
    case MAT:
      out = "materialParam";
      break;
    case GRID_VEL:
      out = "gridVel";
      break;
    default:  
      EXCEPTION("No conversion found for your 'CouplingInputType'" );
    }
  }

  // IntegrationMethod
  template<>
  void String2Enum<IntegrationMethod>( const std::string &in, IntegrationMethod &out ) {
    if ( in == "Gauss_standard" )
      out = ECONOMICAL;

    else if ( in == "Classic" )
      out = CLASSICAL;

    else if ( in == "Lobatto" )
      out = LOBATTO;

    else if ( in == "Chebyshev" )
      out = CHEBYSHEV;

    else if ( in == "Experimental" )
      out = EXPERIMENTAL;

    else if ( in == "Cartesian" )
      out = CARTESIAN;

    else if ( in == "Special" )
      out = SPECIAL;

    else if ( in == "Undefined" )
      out = UNDEFINED;


    else {
      EXCEPTION( "'" << in << "' cannot be converted into item of "
                 << "'IntegrationMethod'!" );
    }
  }

  template<>
  void Enum2String<IntegrationMethod>( const IntegrationMethod &in, std::string &out ) {
    switch(in) {
    case ECONOMICAL:
      out = "Gauss_standard";
      break;

    case CLASSICAL:
      out = "Classic";
      break;

    case LOBATTO:
      out = "Lobatto";
      break;

    case CHEBYSHEV:
      out = "Chebyshev";
      break;

    case EXPERIMENTAL:
      out = "Experimental";
      break;

    case CARTESIAN:
      out = "Cartesian";
      break;

    case SPECIAL:
      out = "Special";
      break;

    case UNDEFINED:
       out = "Undefined";
       break;

    default:
      EXCEPTION( "No conversion found for your 'IntegrationMethod'" );
    }
  }



  // BubbleDynType
  template<>
  void String2Enum<BubbleDynType>( const std::string &in,
                                   BubbleDynType &out ) {

    if ( in == "KellerMiksis" )
      out = KELLERMIKSIS;
    else if ( in == "Gilmore" )
      out = GILMORE;
    else {
      EXCEPTION( "'" << in << "' cannot be converted into item of "
                 << "'BubbleDynType'!" );
    }
  }

  template<>
  void Enum2String<BubbleDynType>( const BubbleDynType &in,
                                   std::string &out ) {
    switch(in) {

    case KELLERMIKSIS:
      out = "KellerMiksis";
      break;
    case GILMORE:
      out = "Gilmore";
      break;
    default:
      EXCEPTION( "No conversion found for your 'BubbleDynType'" );
    }
  }

  // CouplingOutputType
  template<>
  void String2Enum<CouplingOutputType>(const std::string &in,
                                       CouplingOutputType &out) {

    if (in == "nodal")
      out = NODE;
    else if (in == "elem")
      out = ELEM;
    else {
      EXCEPTION( "'" << in << "' cannot be converted into item of "
                 <<"'CouplingOutputType'!" );
    }
  }

  template<>
  void Enum2String<CouplingOutputType>(const CouplingOutputType &in,
                                       std::string &out) {
    switch(in) {

    case NODE:
      out = "nodal";
      break;
    case ELEM:
      out = "elem";
      break;
    default:
      EXCEPTION("No conversion found for your 'CouplingOutputType'" );
    }
  }

  // CouplingRegionType
  template<>
  void String2Enum<CouplingRegionType>(const std::string &in,
                                       CouplingRegionType &out) {
    if (in == "region")
      out = REGION;
    else if (in == "node")
      out = NODES;
    else if (in == "interface")
      out = SURFACE;
    else {
      EXCEPTION( "'" << in << "' cannot be converted into item of "
                 << "'CouplingRegionType'!" );
    }
  }

  template<>
  void Enum2String<CouplingRegionType>(const CouplingRegionType &in,
                                       std::string &out) {

    switch(in) {

    case REGION:
      out = "region";
      break;
    case NODES:
      out = "node";
      break;
    case SURFACE:
      out = "interface";
      break;
    default:
      EXCEPTION( "No conversion found for your 'CouplingRegionType'" );
    }
  }

  // NormType
  template<>
  void String2Enum<NormType>(const std::string &in, NormType &out) {

    if (in == "no")
      out = NO_NORM;
    else if (in == "rel")
      out = L2REL;
    else if (in == "abs")
      out = L2ABS;
    else {
      EXCEPTION( "'" << in << "' cannot be converted into item of "
                 << "'NormType'!" );
    }
  }

  template<>
  void Enum2String<NormType>(const NormType &in, std::string &out) {

    switch(in) {
    case NO_NORM:
      out = "no";
      break;
    case L2REL:
      out = "rel";
      break;
    case L2ABS:
      out = "abs";
      break;
    default:
      EXCEPTION("No conversion found for your 'NormType'");
    }
  }

  // SolutionType
  std::string MapSolTypeToUnit(SolutionType solType)
  {
    switch(solType)
    {

    case ACOU_FORCE:
      return "N";
      break;

    case ACOU_INTENSITY:
      return "W/m^2";
      break;

    case ACOU_POTENTIAL:
      return "m^2/s";
      break;

    case ACOU_POTENTIAL_DERIV_1:
      return "m^2/s^2";
      break;

    case ACOU_POTENTIAL_DERIV_2:
      return "m^2/s^3";
      break;

    case ACOU_POWER:
      return "W";
      break;

    case ACOU_PRESSURE:
      return "Pa";
      break;

    case ACOU_PRESSURE_DERIV_1:
      return "Pa/s";
      break;

    case ACOU_PRESSURE_DERIV_2:
      return "Pa/s^2";
      break;

    case ACOU_PRESSUREXYZ:
      return "Pa";
      break;

    case ACOU_RHS_LOAD:
      return "kg m^-3 s^-2";
      break;

    case ACOU_SURFINTENSITY:
      return "W/m^2";
      break;

    case ACOU_VELOCITY:
      return "m/s";
      break;


    case ELEC_CHARGE:
      return "C";
      break;

    case ELEC_ENERGY:
      return "Ws";
      break;

    case ELEC_POLARIZATION:
      return "C/m^2";
      break;

    case ELEC_FLUX_DENSITY:
      return "C/m^2";
      break;

    case ELEC_PSEUDO_POLARIZATION:
      return "";
      break;

    case ELEC_FIELD_INTENSITY:
      return "V/m";
      break;

    case ELEC_POTENTIAL:
      return "V";
      break;

    case ELEC_RHS_LOAD:
      return "C";
      break;

    case FLUIDMECH_VELOCITY:
      return "m/s";
      break;

    case FLUIDMECH_PRESSURE:
      return "Pa";
      break;

    case FLUIDMECH_DENSITY:
      return "kg/m^3";
      break;

    case FLUIDMECH_TKE:
      return "J";
      break;

    case HEAT_TEMPERATURE:
      return "K";
      break;

    case HEAT_RHS_LOAD:
      return "?";
      break;

    case MAG_FLUX_DENSITY:
      return "Vs/m^2";
      break;

    case MAG_HFIELD:
      return "A/m";
      break;

    case MAG_EDDY_CURRENT:
      return "A/m^2";
      break;

    case MAG_POTENTIAL:
      return "Vs/m";
      break;

    case MAG_RHS_LOAD:
      return "Am";
      break;

    case MAG_EDDY_POWER:
      return "W";
      break;

    case MAG_ENERGY:
      return "Ws";
      break;

    case MAG_FORCE_VWP:
      return "N";
      break;

    case MECH_DEF_VOLUME:
      return "m^3";
      break;

    case MECH_DISPLACEMENT:
      return "m";
      break;

    case MECH_VELOCITY:
      return "m/s";
      break;

    case MECH_ACCELERATION:
      return "m/s^2";
      break;

    case MECH_RHS_LOAD:
      return "N";
      break;

    case MECH_ENERGY:
      return "Ws";
      break;

    case MECH_PSEUDO_DENSITY:
      return "";
      break;

    case MECH_STRESS:
      return "N/m^2";
      break;

    case MECH_STRAIN:
      return "";
      break;

    case MECH_STRAIN_IRR:
      return "N/m^2";
      break;

    case SMOOTH_DISPLACEMENT:
      return "m";
      break;

    case SMOOTH_VELOCITY:
      return "m/s";
      break;

    default:
      return "unknown";
      break;
    }
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

    std::cout << "is " << in << std::endl;
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
  void Enum2String<MaterialType>(const MaterialType &in, std::string &out) {
    switch(in) {
    case NO_MATERIAL:
      out = "noMaterial";
      break;
    case MAG_PERMEABILITY:
      out = "Magnetic_permability";
      break;
    case MAG_PERMEABILITY_1:
      out = "Magnetic_permability_1";
      break;
    case MAG_PERMEABILITY_2:
      out = "Magnetic_permability_2";
      break;
    case MAG_PERMEABILITY_3:
      out = "Magnetic_permability_3";
      break;
    case MAG_RELUCTIVITY:
      out = "Magnetic_reluctivity";
      break;
    case MAG_CONDUCTIVITY:
      out = "Magnetic_Conductiuvity";
      break;
    case ELEC_PERMITTIVITY:
      out = "Electric_Permittivity";
      break;
    case MECH_STIFFNESS_TENSOR:
      out = "MechanicStiffnessTensor";
      break;
    case COEFF_STRAIN_IRREVERSIBLE:
      out = "Coeff_Strain_Irreversible";
      break;
    case MECH_EMODULUS:
      out = "Mechanic_Emodulus";
      break;
    case MECH_EMODULUS_X:
      out = "Mechanic_Emodulus_X";
      break;
    case MECH_EMODULUS_Y:
      out = "Mechanic_Emodulus_Y";
      break;
    case MECH_EMODULUS_Z:
      out = "Mechanic_Emodulus_Z";
      break;
    case MECH_POISSON:
      out = "Mechanic_PoissonRation";
      break;
    case MECH_POISSON_XY:
      out = "Mechanic_PoissonRation_XY";
      break;
    case MECH_POISSON_YZ:
      out = "Mechanic_PoissonRation_YZ";
      break;
    case MECH_POISSON_XZ:
      out = "Mechanic_PoissonRation_XZ";
      break;
    case MECH_KMODULUS:
      out = "Mechanic_Kmodulus";
      break;
    case MECH_GMODULUS:
      out = "Mechanic_Gmodulus";
      break;
    case MECH_GMODULUS_YZ:
      out = "Mechanic_Gmodulus_YZ";
      break;
    case MECH_GMODULUS_ZX:
      out = "Mechanic_Gmodulus_ZX";
      break;
    case MECH_GMODULUS_XY:
      out = "Mechanic_Gmodulus_XY";
      break;
    case MECH_LAME_MU:
      out = "Mechanic_LameMu";
      break;
    case MECH_LAME_LAMBDA:
      out = "Mechanic_LameLambda";
      break;
    case RAYLEIGH_ALPHA:
      out = "Rayleigh_Alpha";
      break;
    case RAYLEIGH_BETA:
      out = "Rayleigh_Beta";
      break;
    case RAYLEIGH_FREQUENCY:
      out = "Rayleigh_Frequency";
      break;
    case RAYLEIGH_DELTA_FREQ:
      out = "Rayleigh_DeltaFreq";
      break;
    case LOSS_TANGENS_DELTA:
      out = "Loss_TangensDelta";
      break;
     case DENSITY:
      out = "Density";
      break;
    case ACOU_BULK_MODULUS:
      out = "AcousticBulkModulus";
      break;
    case ACOU_SOUND_SPEED:
      out = "Acoustic_SoundSpeed";
      break;
    case BOVERA:
      out = "BoverA";
      break;
    case ACOU_ALPHA:
      out = "AcousticAlpha";
      break;
    case FRACTIONAL_ALG:
      out = "FractionalAlg";
      break;
    case FRACTIONAL_MEMORY:
      out = "FractionalMemory";
      break;
    case FRACTIONAL_INTERPOL:
      out = "FractionalInterpol";
      break;
    case FRACTIONAL_EXPONENT:
      out = "FractionalExponent";
      break;
    case HEAT_CONDUCTIVITY:
      out = "HeatConductivity";
      break;
    case HEAT_CONDUCTIVITY_TENSOR:
      out = "HeatConductivity_Tensor";
      break;
    case HEAT_CAPACITY:
      out = "HeatCapacity";
      break;
    case THERMAL_EXPANSION_TENSOR:
      out = "thermalExpansion";
      break;
    case DYNAMIC_VISCOSITY:
      out = "dynamicViscosity";
      break;
    case KINEMATIC_VISCOSITY:
      out = "kinematicViscosity";
      break;
    case PIEZO_TENSOR:
      out = "PiezoTensor";
      break;
    case SPON_POLARIZATION:
      out = "sponPolarization";
      break;
    case SPON_STRAIN:
      out = "sponStrain";
      break;
    case EFIELD0:
      out = "Efield0";
      break;
    case STRESS0:
      out = "Stress0";
      break;
    case DCOUPLE0:
      out = "dCouple0";
      break;
    case RATE_CONSTANT:
      out = "rateConstant";
      break;
    case VISCO_PLASTIC_INDEX:
      out = "viscoPlasticIndex";
      break;
    case SATURATION_INDEX:
      out = "saturationIndex";
      break;
    case  SCALE_FORCE_ELEC:
      out = "scaleForceElec";
      break;
    case  SCALE_FORCE_MECH:
      out = "scaleForceMech";
      break;
    case  SCALE_FORCE_COUPLE:
      out = "scaleForceCouple";
      break;
    case VOLUME_FRAC_INIT:
      out = "volumeFracInit";
      break;
    case MEAN_TEMPERATURE:
      out = "meanTemperature";
      break;

    case X_SATURATION:
      out = "Xsaturation";
      break;
    case Y_SATURATION:
      out = "Ysaturation";
      break;
    case Y_REMANENCE:
      out = "Yremanence";
      break;
    case PREISACH_WEIGHTS:
      out = "preisachWeights";
      break;
    case A_JILES:
      out = "aJiles";
      break;
    case ALPHA_JILES:
      out = "alphaJiles";
      break;
    case K_JILES:
      out = "kJiles";
      break;
    case C_JILES:
      out = "cJiles";
      break;
    case P_DIRECTION:
      out = "Pdirection";
      break;
    case HYST_MODEL:
      out = "hystModel";
      break;
    case NONLIN_COEFFICIENT:
      out = "nonLinCoefficient";
      break;
    case NONLIN_DEPENDENCY:
      out = "nonLinDependency";
      break;
    case NONLIN_APPROXIMATION_TYPE:
      out = "nonLinApproximationType";
      break;
    case NONLIN_DATA_NAME:
      out = "nonLinDataName";
      break;
    case DATA_ACCURACY:
      out = "dataAccuracy";
      break;
    case MAX_APPROX_VAL:
      out = "maxApproxVal";
      break;
    case PYROCOEFFICIENT_TENSOR:
      out = "Pyrocoefficient_Tensor";
      break;
    default:
      EXCEPTION("No conversion found for your 'DataType'");
    }
  }

  template<>
  void String2Enum<MaterialType>( const std::string &in, MaterialType &out ) {
    if ( in == "noMaterial" ) {
      out = NO_MATERIAL;
    }
    else if ( in == "Magnetic_permability" ) {
      out = MAG_PERMEABILITY;
    }
    else if ( in == "Magnetic_permability_1" ) {
      out = MAG_PERMEABILITY_1;
    }
    else if ( in == "Magnetic_permability_2" ) {
      out = MAG_PERMEABILITY_2;
    }
    else if ( in == "Magnetic_permability_3" ) {
      out = MAG_PERMEABILITY_3;
    }
    else if ( in == "Magnetic_reluctivity" ) {
      out = MAG_RELUCTIVITY;
    }
    else if ( in == "Magnetic_Conductiuvity" ) {
      out = MAG_CONDUCTIVITY;
    }
    else if ( in == "Electric_Permittivity" ) {
      out = ELEC_PERMITTIVITY;
    }
    else if ( in == "MechanicStiffnessTensor" ) {
      out = MECH_STIFFNESS_TENSOR;
    }
    else if ( in == "Coeff_Strain_Irreversibel" ) {
      out = COEFF_STRAIN_IRREVERSIBLE;
    }
    else if ( in == "Mechanic_Emodulus" ) {
      out = MECH_EMODULUS;
    }
    else if ( in == "Mechanic_Emodulus_X" ) {
      out = MECH_EMODULUS_X;
    }
    else if ( in == "Mechanic_Emodulus_Y" ) {
      out = MECH_EMODULUS_Y;
    }
    else if ( in == "Mechanic_Emodulus_Z" ) {
      out = MECH_EMODULUS_Z;
    }
    else if ( in == "Mechanic_PoissonRation" ) {
      out = MECH_POISSON;
    }
    else if ( in == "Mechanic_PoissonRation_XY" ) {
      out = MECH_POISSON_XY;
    }
    else if ( in == "Mechanic_PoissonRation_YZ" ) {
      out = MECH_POISSON_YZ;
    }
    else if ( in == "Mechanic_PoissonRation_XZ" ) {
      out = MECH_POISSON_XZ;
    }
    else if ( in == "Mechanic_Kmodulus" ) {
      out = MECH_KMODULUS;
    }
    else if ( in == "Mechanic_Gmodulus" ) {
      out = MECH_GMODULUS;
    }
    else if ( in == "Mechanic_Gmodulus_YZ" ) {
      out = MECH_GMODULUS_YZ;
    }
    else if ( in == "Mechanic_Gmodulus_ZX" ) {
      out = MECH_GMODULUS_ZX;
    }
    else if ( in == "Mechanic_Gmodulus_XY" ) {
      out = MECH_GMODULUS_XY;
    }
    else if ( in == "Mechanic_LameMu" ) {
      out = MECH_LAME_MU;
    }
    else if ( in == "Mechanic_LameLambda" ) {
      out = MECH_LAME_LAMBDA;
    }
    else if ( in == "Rayleigh_Alpha" ) {
      out = RAYLEIGH_ALPHA;
    }
    else if ( in == "Rayleigh_Beta" ) {
      out = RAYLEIGH_BETA;
    }
    else if ( in == "Rayleigh_Frequency" ) {
      out = RAYLEIGH_FREQUENCY;
    }
    else if ( in == "Rayleigh_DeltaFreq" ) {
      out = RAYLEIGH_DELTA_FREQ;
    }
    else if ( in == "Loss_TangensDelta" ) {
      out = LOSS_TANGENS_DELTA;
    }
    else if ( in == "Density" ) {
      out = DENSITY;
    }
    else if ( in == "AcousticBulkModulus" ) {
      out = ACOU_BULK_MODULUS;
    }
    else if ( in == "Acoustic_SoundSpeed" ) {
      out = ACOU_SOUND_SPEED;
    }
    else if ( in == "BoverA" ) {
      out = BOVERA;
    }
    else if ( in == "AcousticAlpha" ) {
      out = ACOU_ALPHA;
    }
    else if ( in == "FractionalAlg" ) {
      out = FRACTIONAL_ALG;
    }
    else if ( in == "FractionalMemory" ) {
      out = FRACTIONAL_MEMORY;
    }
    else if ( in == "FractionalInterpol" ) {
      out = FRACTIONAL_INTERPOL;
    }
    else if ( in == "FractionalExponent" ) {
      out = FRACTIONAL_EXPONENT;
    }
    else if ( in == "HeatConductivity" ) {
      out = HEAT_CONDUCTIVITY;
    }
    else if ( in == "HeatConductivity_Tensor" ) {
      out = HEAT_CONDUCTIVITY_TENSOR;
    }
    else if ( in == "HeatCapacity" ) {
      out = HEAT_CAPACITY;
    }
    else if ( in == "thermalExpansion_Tensor" ) {
      out = THERMAL_EXPANSION_TENSOR;
    }
    else if ( in == "dynamicViscosity" ) {
      out = DYNAMIC_VISCOSITY;
    }
    else if ( in == "kinematicViscosity" ) {
      out = KINEMATIC_VISCOSITY;
    }
    else if ( in == "PiezoTensor" ) {
      out = PIEZO_TENSOR;
    }
    else if ( in == "Xsaturation" ) {
      out = X_SATURATION;
    }
    else if ( in == "Ysaturation" ) {
      out = Y_SATURATION;
    }
    else if ( in == "Yremanence" ) {
      out = Y_REMANENCE;
    }
    else if ( in == "sponPolarization" ) {
      out = SPON_POLARIZATION;
    }
    else if ( in == "sponStrain" ) {
      out = SPON_STRAIN;
    }
    else if ( in == "Efield0" ) {
      out = EFIELD0;
    }
    else if ( in == "Stress0" ) {
      out = STRESS0;
    }
    else if ( in == "dCouple0" ) {
      out = DCOUPLE0;
    }
    else if ( in == "rateConstant" ) {
      out = RATE_CONSTANT;
    }
    else if ( in == "viscoPlasticIndex" ) {
      out = VISCO_PLASTIC_INDEX;
    }
    else if ( in == "saturationIndex" ) {
      out = SATURATION_INDEX;
    }
    else if ( in == "scaleForceElec" ) {
      out =  SCALE_FORCE_ELEC;
    }
    else if ( in == "scaleForceMech" ) {
      out =  SCALE_FORCE_MECH;
    }
    else if ( in == "scaleForceCouple" ) {
      out =  SCALE_FORCE_COUPLE;
    }
    else if ( in == "volumeFracInit" ) {
      out = VOLUME_FRAC_INIT;
    }
    else if ( in == "meanTemperature" ) {
      out =  MEAN_TEMPERATURE;
    }

    else if ( in == "preisachWeights" ) {
      out = PREISACH_WEIGHTS;
    }
    else if ( in == "aJiles" ) {
      out = A_JILES;
    }
    else if ( in == "alphaJiles" ) {
      out = ALPHA_JILES;
    }
    else if ( in == "kJiles" ) {
      out = K_JILES;
    }
    else if ( in == "cJiles" ) {
      out = C_JILES;
    }
    else if ( in == "Pdirection" ) {
      out = P_DIRECTION;
    }
    else if ( in == "hystModel" ) {
      out = HYST_MODEL;
    }
    else if ( in == "nonLinCoefficient" ) {
      out = NONLIN_COEFFICIENT;
    }
    else if ( in == "nonLinDependency" ) {
      out = NONLIN_DEPENDENCY;
    }
    else if ( in == "nonLinApproximationType" ) {
      out = NONLIN_APPROXIMATION_TYPE;
    }
    else if ( in == "nonLinDataName" ) {
      out = NONLIN_DATA_NAME;
    }
    else if ( in == "dataAccuracy" ) {
      out = DATA_ACCURACY;
    }
    else if ( in == "maxApproxVal" ) {
      out = MAX_APPROX_VAL;
    }
    else if ( in == "Pyrocoefficient_Tensor" ) {
      out = PYROCOEFFICIENT_TENSOR;
    }
    else {
      EXCEPTION("No conversion from string to 'MaterialType' found");
    }
  }


  template<>
  void Enum2String<SubTensorType>(const SubTensorType &in, std::string &out) {
    switch(in) {
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
    else if ( in == "3d" ) {
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
    case ELECTROMAGNETIC:
      out = "magnetic";
      break;
    case ELECTROSTATIC:
      out = "electric";
      break;
    case FLUID:
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

    default:
      EXCEPTION("No conversion found for your 'MaterialClass'" );
    }
  }

  template<>
  void String2Enum<MaterialClass>( const std::string &in, MaterialClass &out ) {

    if ( in == "No MaterialClass" ) {
      out = NO_CLASS;
    }
    else if ( in == "electromagnetic" ) {
      out = ELECTROMAGNETIC;
    }
    else if ( in == "electrostatic" ) {
      out = ELECTROSTATIC;
    }
    else if ( in == "electromagnetic" ) {
      out = ELECTROMAGNETIC;
    }
    else if ( in == "fluid" ) {
      out = FLUID;
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

  template<>
  void String2Enum<NonLinType>( const std::string &in, NonLinType &out ) {

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
    case PIEZO_MICRO_HF:
      out = "piezoMicroHF";
      break;
    case PERMEABILITY:
      out = "permeability";
      break;
    default:
      EXCEPTION( "No conversion found for 'NonLinType' " << in );
    }
  }

  template<>
  void String2Enum<DampingType>( const std::string &in, DampingType &out ) {

    if( in == "none" ) {
      out = NONE;
    } else if( in == "rayleigh" ) {
      out = RAYLEIGH;
    } else if( in == "abc") {
      out = ABCDAMP;
    } else if( in == "thermoViscous") {
      out = THERMOVISCOUS;
    } else if( in == "fractional") {
      out = FRACTIONAL;
    } else if( in == "fractiona_gl") {
      out = FRACTIONAL_GL;
    } else if( in == "fractional_blank") {
      out = FRACTIONAL_BLANK;
    } else if( in == "fractional_gl_int") {
      out = FRACTIONAL_GL_INT;
    } else if( in == "fractional_blank_int") {
      out = FRACTIONAL_BLANK_INT;
    } else if( in == "pml" ) {
      out = PML;
    } else if( in == "dampLayer" ) {
      out = DAMPLAYER;
    } else {
      EXCEPTION( "'" << in << "' cannot be converted into an "
                 << "'DampingType' item!" );
    }
  }

  template<>
  void Enum2String<DampingType>( const DampingType &in, std::string& out ) {
    switch(in) {
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
  default:
    EXCEPTION( "No conversion found for 'DapmingType' " << in );
    }
  }

  // ****************************************************************
  // ****************************************************************
  //              THE OLAS ENVIRONMENT STARTS HERE
  // ****************************************************************
  // ****************************************************************


  // Specialisation for SolverType
  template<>
  void Enum2String<SolverType>(const SolverType &in,
                                      std::string &out) {
    switch( in ) {
    case NOSOLVER:
      out = "no solver";
      break;
    case RICHARDSON:
      out = "Richardson";
      break;
    case DIAGSOLVER:
      out = "diagsolver";
      break;
    case CG:
      out = "cg";
      break;
    case GMRES:
      out = "gmres";
      break;
    case MINRES:
      out = "minres";
      break;
    case SYMMLQ:
      out = "symmlq";
      break;
    case LAPACK_LU:
      out = "lapackLU";
      break;
    case LAPACK_LL:
      out = "lapackLL";
      break;
    case LU_SOLVER:
      out = "directLU";
      break;
    case LDL_SOLVER:
      out = "directLDL";
      break;
    case LDL_SOLVER2:
      out = "directLDL2";
      break;
    case PARDISO:
      out = "pardiso";
      break;
    case ILUPACK_SOLVER:
      out = "ilupack";
      break;
    case CHOLMOD:
      out = "cholmod";
      break;



    default:
      EXCEPTION( "No string value found for the specified value of the "
           << "enumeration datatype SolverType.\n"
           << "Seems to indicate a missing case implementation!" );
    }
  }

  // Specialisation for EigenSolverType
  template<>
  void Enum2String<EigenSolverType>(const EigenSolverType &in,
                                      std::string &out) {
    switch( in ) {
    case NOEIGENSOLVER:
      out = "no eigensolver";
      break;
    case ARPACK:
      out = "arpack";
      break;
    case SUBSPACE:
      out = "subspace";
      break;
    default:
      EXCEPTION( "No string value found for the specified value of the "
               << "enumeration datatype EigenSolverType.\n"
               << "Seems to indicate a missing case implementation!" );
    }
  }


  // Specialisation for PrecondType
  template<>
  void Enum2String<PrecondType>(const PrecondType &in,
                                      std::string &out) {
    switch( in ) {
    case NOPRECOND:
      out = "no precond";
      break;
    case ID:
      out = "Id";
      break;
    case MG:
      out = "MG";
      break;
    case JACOBI:
      break;
    case SSOR:
      out = "SSOR";
      break;
    case ILU0:
      out = "ILU0";
      break;
    case ILUTP:
      out = "ILUTP";
      break;
    case ILUK:
      out = "ILUK";
      break;
    case ILDL0:
      out = "ILDL0";
      break;
    case ILDLK:
      out = "ILDLK";
      break;
    case ILDLTP:
      out = "ILDLTP";
      break;
    case ILDLCN:
      out = "ILDLCN";
      break;
    case IC0:
      out = "IC0";
      break;

    default:
      EXCEPTION( "No string value found for the specified value of the "
           << "enumeration datatype PrecondType.\n"
           << "Seems to indicate a missing case implementation!" );
    }
  }

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

  // Specialisation for ReorderingType
  template<>
  void Enum2String<ReorderingType>(const ReorderingType &in,
                                        std::string &out) {
    switch( in ) {
    case NOREORDERING:
      out = "noReordering";
      break;
    case SLOAN:
      out = "Sloan";
      break;
    case METIS:
      out = "Metis";
      break;
    case MINIMUM_DEGREE:
      out = "minimumDegree";
      break;
    case NESTED_DISSECTION:
      out = "nestedDissection";
      break;
    default:
      EXCEPTION( "No string value found for the specified value of the "
           << "enumeration datatype ReorderingType.\n"
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
    case CONVECTION:
      out = "convection";
      break;
    case MASS:
      out = "mass";
      break;
    case AUXILIARY:
      out = "auxiliary";
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

  // Specialisation for SolverType
  template<>
  void String2Enum<SolverType>( const std::string &in, SolverType &out ) {

    if ( in == "no solver" ) {
      out = NOSOLVER;
    }
    else if ( in == "Richardson" ) {
      out = RICHARDSON;
    }
    else if ( in == "diagsolver" ) {
      out = DIAGSOLVER;
    }
    else if ( in == "cg" ) {
      out = CG;
    }
    else if ( in == "gmres" ) {
      out = GMRES;
    }
    else if ( in == "minres" ) {
      out = MINRES;
    }
    else if ( in == "symmlq" ) {
      out = SYMMLQ;
    }
    else if ( in == "lapackLL" ) {
      out = LAPACK_LL;
    }
    else if ( in == "directLU" ) {
      out = LU_SOLVER;
    }
    else if ( in == "directLDL" ) {
      out = LDL_SOLVER;
    }
    else if ( in == "directLDL2" ) {
      out = LDL_SOLVER2;
    }
    else if ( in == "pardiso" ) {
      out = PARDISO;
    }
    else if ( in == "ilupack" ) {
      out = ILUPACK_SOLVER;
    }
    else if ( in == "cholmod") {
      out = CHOLMOD;
    }
    else {
      EXCEPTION( "No enumeration value found in SolverType for '"
           << in << "'\n A missing case implementation?" );
    }
  }

 // Specialisation for EigenSolverType
  template<>
  void String2Enum<EigenSolverType>( const std::string &in, EigenSolverType &out ) {

    if ( in == "no eigensolver" ) {
      out = NOEIGENSOLVER;
    }
    else if ( in == "arpack" ) {
      out = ARPACK;
    }
    else if ( in == "subspace" ) {
      out = SUBSPACE;
    }
    else {
      EXCEPTION( "No enumeration value found in EigenSolverType for '"
               << in << "'\n A missing case implementation?" );
    }
  }
  // Specialisation for PrecondType
  template<>
  void String2Enum<PrecondType>( const std::string &in, PrecondType &out ) {

    if ( in == "noPrecond" ) {
      out = NOPRECOND;
    }
    else if ( in == "Id" ) {
      out = ID;
    }
    else if ( in == "MG" ) {
      out = MG;
    }
    else if ( in == "Jacobi" ) {
      out = JACOBI;
    }
    else if ( in == "SSOR" ) {
      out = SSOR;
    }
    else if ( in == "ILU0" ) {
      out = ILU0;
    }
    else if ( in == "ILUTP" ) {
      out = ILUTP;
    }
    else if ( in == "ILUK" ) {
      out = ILUK;
    }
    else if ( in == "ILDL0" ) {
      out = ILDL0;
    }
    else if ( in == "ILDLK" ) {
      out = ILDLK;
    }
    else if ( in == "ILDLTP" ) {
      out = ILDLTP;
    }
    else if ( in == "ILDLCN" ) {
      out = ILDLCN;
    }
    else if ( in == "IC0" ) {
      out = IC0;
    }
    else {
      EXCEPTION( "No enumeration value found in PrecondType for '"
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
                                     AMGSmootherType   &out )
  {
    if ( in == "GaussSeidel" ) {
      out = AMG_SMOOTHER_GAUSSSEIDEL;
    } else if( in == "Jacobi" ) {
      out = AMG_SMOOTHER_DAMPED_JACOBI;
    } else {
      EXCEPTION( "No enumeration value found in AMGSmoothertype "
       << "for '" << in << "'\n A missing case implementation?" );
    }
  }

  // Specialisation for ReorderingType
  template<>
  void String2Enum<ReorderingType>( const std::string &in,
                    ReorderingType &out ) {

    if ( in == "noReordering" )
      out = NOREORDERING;
    else if ( in == "Sloan" )
      out = SLOAN;
    else if ( in == "Metis" )
      out = METIS;
    else if ( in == "minimumDegree" )
      out = MINIMUM_DEGREE;
    else if ( in == "nestedDissection" )
      out = NESTED_DISSECTION;
    else {
      EXCEPTION( "String '" << in << "' cannot be converted to item of "
           << "'ReorderingType'!"; );
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
  
  void SetEnvironmentEnums(){
    // SolutionType
    
    SolutionTypeEnum.SetName("SolutionTypeEnum");
    //mechanics
    SolutionTypeEnum.Add(MECH_DISPLACEMENT, "mechDisplacement");
    SolutionTypeEnum.Add(MECH_ACCELERATION, "mechAcceleration");
    SolutionTypeEnum.Add(MECH_VELOCITY, "mechVelocity");
    SolutionTypeEnum.Add(MECH_FORCE, "mechForce");
    SolutionTypeEnum.Add(MECH_STRESS, "mechStress");
    SolutionTypeEnum.Add(MECH_STRAIN, "mechStrain");
    SolutionTypeEnum.Add(MECH_STRAIN_IRR, "mechStrainIrr");
    SolutionTypeEnum.Add(MECH_ENERGY, "mechEnergy");
    SolutionTypeEnum.Add(MECH_DEF_VOLUME, "volumeAboveDefSurf");
    SolutionTypeEnum.Add(MECH_RHS_LOAD, "mechRhsLoad");
    SolutionTypeEnum.Add(MECH_PSEUDO_DENSITY, "mechPseudoDensity");
    SolutionTypeEnum.Add(MECH_SHAPE, "mechShape");
    //electrostatics
    SolutionTypeEnum.Add(ELEC_POTENTIAL, "elecPotential");
    SolutionTypeEnum.Add(ELEC_FIELD_INTENSITY, "elecFieldIntensity");
    SolutionTypeEnum.Add(ELEC_POLARIZATION, "elecPolarization");
    SolutionTypeEnum.Add(ELEC_PSEUDO_POLARIZATION, "elecPseudoPolarization");
    SolutionTypeEnum.Add(ELEC_FORCE_VWP, "elecForceVWP");
    SolutionTypeEnum.Add(ELEC_INTERFACE_FORCE, "elecInterfaceForce");
    SolutionTypeEnum.Add(ELEC_CHARGE, "elecCharge");
    SolutionTypeEnum.Add(ELEC_FLUX_DENSITY, "elecFluxDensity");
    SolutionTypeEnum.Add(ELEC_ENERGY, "elecEnergy");
    SolutionTypeEnum.Add(ELEC_RHS_LOAD, "elecRhsLoad");
    //smoothing PDE
    SolutionTypeEnum.Add(SMOOTH_DISPLACEMENT, "smoothDisplacement");
    SolutionTypeEnum.Add(SMOOTH_VELOCITY, "smoothVelocity");
    SolutionTypeEnum.Add(GRID_VELOCITY, "gridVelocity");
    SolutionTypeEnum.Add(SMOOTH_STRAIN, "smoothStrain");
    //acoustics
    SolutionTypeEnum.Add(ACOU_PRESSURE, "acouPressure");
    SolutionTypeEnum.Add(ACOU_POTENTIAL, "acouPotential");
    SolutionTypeEnum.Add(ACOU_VELOCITY, "acouVelocity");
    SolutionTypeEnum.Add(ACOU_PRESSURE_DERIV_1, "acouPressureD1");
    SolutionTypeEnum.Add(ACOU_PRESSURE_DERIV_2, "acouPressureD2");
    SolutionTypeEnum.Add(ACOU_FORCE, "acouForce");
    SolutionTypeEnum.Add(ACOU_POTENTIAL_DERIV_1, "acouPotentialD1");
    SolutionTypeEnum.Add(ACOU_POTENTIAL_DERIV_2, "acouPotentialD2");
    SolutionTypeEnum.Add(ACOU_RHS_LOAD, "acouRhsLoad");
    SolutionTypeEnum.Add(ACOU_RHSVAL, "acouRHSval");
    SolutionTypeEnum.Add(ACOUSURF_RHSVAL, "acouSurfRHSval");
    SolutionTypeEnum.Add(ACOU_SOUND_SPEEED, "acouSoundSpeed");
    SolutionTypeEnum.Add(ACOU_BUBBLE_RHS_VAL, "acouBubbleRhsVal");
    SolutionTypeEnum.Add(ACOU_POT_NRBC, "acouPotNRBC");
    SolutionTypeEnum.Add(NRBC_PHI, "nrbcPhi");
    SolutionTypeEnum.Add(ACOU_PRESSUREXYZ, "acouPressureXYZ");
    SolutionTypeEnum.Add(ACOU_POWERDENSITY, "acouPowerDensity");
    SolutionTypeEnum.Add(ACOU_POWER, "acouPower");
    SolutionTypeEnum.Add(ACOU_INTENSITY, "acouIntensity");
    SolutionTypeEnum.Add(ACOU_SURFINTENSITY, "acouSurfIntensity");
    //magnetics
    SolutionTypeEnum.Add(MAG_POTENTIAL, "magPotential");
    SolutionTypeEnum.Add(MAG_FLUX_DENSITY, "magFluxDensity");
    SolutionTypeEnum.Add(MAG_POTENTIAL_DIV, "magPotentialDiv");
    SolutionTypeEnum.Add(MAG_HFIELD, "magHfield");
    SolutionTypeEnum.Add(MAG_EDDY_CURRENT, "magEddyCurrent");
    SolutionTypeEnum.Add(MAG_FORCE_VWP, "magForceVWP");
    SolutionTypeEnum.Add(MAG_FORCE_LORENTZ, "magForceLorentz");
    SolutionTypeEnum.Add(MAG_ENERGY, "magEnergy");
    SolutionTypeEnum.Add(MAG_EDDY_POWER, "magEddyPower");
    SolutionTypeEnum.Add(MAG_RHS_LOAD, "magRhsLoad");
    //heat conduction
    SolutionTypeEnum.Add(HEAT_TEMPERATURE, "heatTemperature");
    SolutionTypeEnum.Add(HEAT_RHS_LOAD, "heatRhsLoad");
    //mpcci
    SolutionTypeEnum.Add(FLUID_FORCE, "fluidForce");
    //fluidMech
    SolutionTypeEnum.Add(FLUIDMECH_VELOCITY, "fluidMechVelocity");
    SolutionTypeEnum.Add(FLUIDMECH_PRESSURE, "fluidMechPressure");
    SolutionTypeEnum.Add(FLUIDMECH_FORCE, "fluidMechForce");
    SolutionTypeEnum.Add(FLUIDMECH_DENSITY, "fluidMechDensity");
    SolutionTypeEnum.Add(FLUIDMECH_TKE, "fluidMechTKE");
    SolutionTypeEnum.Add(LAMBDA_K, "lambda_k");
    // bubble
    SolutionTypeEnum.Add(BUBBLE_RADIUS, "bubbleRadius");
    SolutionTypeEnum.Add(BUBBLE_RADIUS_DERIV_1, "bubbleRadiusD1");
    // optimization
    SolutionTypeEnum.Add(HOMOGENIZED_TENSOR, "homogenizedTensor");
    // the actual result type is given in result descriptions
    // in the xml file in the optimization element.
    SolutionTypeEnum.Add(OPT_RESULT_1, "optResult_1");
    SolutionTypeEnum.Add(OPT_RESULT_2, "optResult_2");
    SolutionTypeEnum.Add(OPT_RESULT_3, "optResult_3");
    SolutionTypeEnum.Add(OPT_RESULT_4, "optResult_4");
    SolutionTypeEnum.Add(OPT_RESULT_5, "optResult_5");
    SolutionTypeEnum.Add(OPT_RESULT_6, "optResult_6");
    SolutionTypeEnum.Add(OPT_RESULT_7, "optResult_7");
    SolutionTypeEnum.Add(OPT_RESULT_8, "optResult_8");
    SolutionTypeEnum.Add(OPT_RESULT_9, "optResult_9");
    // independent
    SolutionTypeEnum.Add(LAGRANGE_MULT, "LagrangeMultiplier");
  }
  
  Enum<SolutionType> SolutionTypeEnum;

}

