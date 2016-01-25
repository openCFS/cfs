#include <iostream>
#include <fstream>
#include <string>

#include "Environment.hh"
#include "Utils/tools.hh"
#include "Domain/Domain.hh"

// Since OLAS uses a separate namespace for
// writing out data, two different declarations
// have to be made
namespace OutInfo{

  std::ostream *cla      = NULL;
}

namespace CoupledField {

  // Define global objects 
  Domain * domain = NULL;
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
        
      case ACOU_NORMAL_INTENSITY:
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

      case WATER_PRESSURE:
        return "Pa";
        break;

      case ACOU_ACCELERATION:
        return "m/s^2";
        break;

      case ACOU_PRESSURE_DERIV_1:
        return "Pa/s";
        break;

      case ACOU_PRESSURE_DERIV_2:
        return "Pa/s^2";
        break;

      case ACOU_RHS_LOAD:
        return "kg m^-3 s^-2";
        break;

      case ACOU_RHS_LOADP:
        return "kg m^-3 s^-2";
        break;

      case ACOU_APE_RHS_LOAD:
        return "kg m^-3 s^-2";
        break;

      case ACOU_VORTEX_RHS_LOAD:
        return "kg m^-3 s^-2";
        break;

      case ACOU_RHS_LOAD_DENSITY:
        return "kg m^-6 s^-2";
        break;

      case ACOU_DIV_LH_TENSOR:
        return "kg m^-2 s^-2";
        break;

      case ACOU_VELOCITY:
        return "m/s";
        break;

      case ACOU_SURFINTENSITY:
        return "W/m^2";
        break;

      case ACOU_ENERGY:
        return "Ws";
        break;

      case ACOU_POT_ENERGY:
        return "Ws";
        break;

      case ACOU_KIN_ENERGY:
        return "Ws";
        break;

      case ACOU_ELEM_SPEED_OF_SOUND:
        return "m/s";
        break;

      case ACOU_MIXED_MASS_LOAD:
        return "-";
        break;

      case ACOU_MIXED_MOMENTUM_LOAD:
        return "-";
        break;

      case ACOU_LAMB_RHS:
        return "-";
        break;

      case ELEC_CHARGE:
        return "C";
        break;

      case ELEC_CHARGE_DENSITY:
        return "C/m^3";
        break;

      case ELEC_ENERGY:
        return "Ws";
        break;
        
      case ELEC_ENERGY_DENSITY:
        return "Ws/m^3";
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

      case ELEC_CURRENT_DENSITY:
        return "A/m^2";
        break;

      case ELEC_NORMAL_CURRENT_DENSITY:
        return "A/m^2";
        break;

      case ELEC_CURRENT:
        return "A";
        break;

      case ELEC_POWER_DENSITY:
        return "W/m^3";
        break;

      case ELEC_POWER:
        return "W";
        break;

      case FLUIDMECH_VELOCITY:
        return "m/s";
        break;

      case MEAN_FLUIDMECH_VELOCITY:
        return "m/s";
        break;

      case FLUIDMECH_VELOCITY_DERIV_1:
        return "1/s";
        break;

      case FLUIDMECH_VELOCITY_DERIV_2:
        return "1/s^2";
        break;

      case MEAN_FLUIDMECH_VELOCITY_NORMAL:
        return "m/s";
        break;

      case DIV_MEAN_FLUIDMECH_VELOCITY:
        return "1/s";
        break;

      case FLUIDMECH_PRESSURE:
        return "Pa";
        break;

      case FLUIDMECH_PRESSURE_DERIV_1:
        return "Pa/m";
        break;

      case FLUIDMECH_PRESSURE_DERIV_2:
        return "Pa/m^2";
        break;

      case FLUIDMECH_PRESSURE_TIME_DERIV_1:
        return "Pa/s";
        break;

      case FLUIDMECH_PRESSURE_TIME_DERIV_2:
        return "Pa/s^2";
        break;

      case FLUIDMECH_DENSITY:
        return "kg/m^3";
        break;



      case FLUIDMECH_TKE:
        return "J";
        break;

      case FLUIDMECH_STRESS:
        return "Pa";
        break;

      case FLUIDMECH_STRAINRATE:
        return "1/s";
        break;

      case FLUIDMECH_WVT:
        return "kg m^-2 s^-2";
        break;

      case FLUIDMECH_WVT_DENSITY:
        return "kg m^-2 s^-2 m s^-1";
        break;

      case FLUIDMECH_WVT_PHI:
        return "kg m^-3 s^-1";
        break;

      case FLUIDMECH_WVT_DENSITY_PHI:
        return "kg m^-3 s^-1 m s^-1";
        break;

      case FLUIDMECH_WVT_U1:
        return "m s^-1";
        break;

      case FLUIDMECH_WVT_U2:
        return "m s^-1";
        break;

      case FLUIDMECH_WVT_F:
        return "N";
        break;

      case HEAT_TEMPERATURE:
        return "K";
        break;
        
      case HEAT_MEAN_TEMPERATURE:
        return "K";
        break;

      case HEAT_TEMPERATURE_D1:
        return "K/s";
        break;

      case HEAT_FLUX_DENSITY:
        return "W/m^2";
        break;

      case HEAT_RHS_LOAD:
        return "?";
        break;

      case HEAT_SOURCE_DENSITY:
        return "W/m^3";
        break;

      case MAG_FLUX_DENSITY:
      case MAG_NORMAL_FLUX_DENSITY:
      case COIL_LINKED_FLUX:
        return "Vs/m^2";
        break;

      case MAG_FLUX:
        return "Vs";
        break;

      case MAG_FIELD_INTENSITY:
        return "A/m";
        break;

      case MAG_EDDY_CURRENT_DENSITY:
      case MAG_TOTAL_CURRENT_DENSITY:
        return "A/m^2";
        break;
        
      case MAG_POTENTIAL:
        return "Vs/m";
        break;
        
      case MAG_POTENTIAL_DERIV1:
        return "V/m";
        break;

      case MAG_SCALAR_POTENTIAL:
      case COIL_CURRENT:
        return "A";
        break;

      case COIL_CURRENT_DERIV1:
        return "A/s";
        break;

      case MAG_RHS_LOAD:
        return "Am";
        break;

      case MAG_ELEM_PERMEABILITY:
        return "Vs/Am";

      case MAG_EDDY_POWER:
      case MAG_CORE_LOSS:
        return "W";
        break;

      case MAG_CORE_LOSS_DENSITY:
        return "W/kg";
        break;

      case MAG_ENERGY:
        return "Ws";
        break;

      case MAG_FORCE_VWP:
        return "N";
        break;

      case MECH_DEF_SURF_VOLUME:
        return "m^3";
        break;

      case MECH_DISPLACEMENT:
        return "m";
        break;

      case MECH_NORMAL_DISPLACEMENT:
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

      case MECH_KIN_ENERGY:
      case MECH_DEFORM_ENERGY: 
      case MECH_TOTAL_ENERGY:
        return "Ws";
        break;

      case MECH_STRESS:
        return "N/m^2";
        break;

      case MECH_THERMAL_STRESS:
        return "N/m^2";
        break;

      case MECH_STRAIN:
        return "";
        break;

      case MECH_THERMAL_STRAIN:
         return "";
         break;

      case MECH_NORMAL_STRESS:
        return "N/m^2";
        break;

      case SMOOTH_DISPLACEMENT:
        return "m";
        break;

      case SMOOTH_VELOCITY:
        return "m/s";
        break;

      case MECH_PSEUDO_DENSITY:
      case PHYSICAL_PSEUDO_DENSITY:
      case ELEC_PHYSICAL_PSEUDO_DENSITY:
      case LBM_PHYSICAL_PSEUDO_DENSITY:
        return "";
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
      case TESTMAT:
        out = "testmat";
        break;
      case ELECTROMAGNETIC:
        out = "magnetic";
        break;
      case ELECTROSTATIC:
        out = "electric";
        break;
      case ELECTRICCONDUCTION:
        out = "elecConduction";
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
      case MAGNETOSTRICTIVE:
        out = "magnetoStrictive";
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
    else if ( in == "elecConduction" ) {
      out = ELECTRICCONDUCTION;
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
    else if ( in == "magnetoStrictive" ) {
      out = MAGNETOSTRICTIVE;
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
    } else if( in == "heatConductivity") {
      out = NLHEAT_CONDUCTIVITY;
    } else if( in == "heatCapacity") {
      out = NLHEAT_CAPACITY;
    } else if( in == "elecConductivity") {
       out = NLELEC_CONDUCTIVITY;
    } else if( in == "elecBiPole") {
       out = NLELEC_BIPOLE;
    } else if( in == "elecBiPoleTempDep") {
       out = NLELEC_BIPOLE_TEMP_DEP;
    } else if( in == "elecTriPole") {
       out = NLELEC_TRIPOLE;
    } else if( in == "elecTriPoleTempDep") {
       out = NLELEC_TRIPOLE_TEMP_DEP;
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
      case NLHEAT_CONDUCTIVITY:
        out = "heatConductivity";
        break;
      case NLHEAT_CAPACITY:
        out = "heatCapacity";
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
    SolutionTypeEnum.Add(MECH_NORMAL_DISPLACEMENT, "mechNormalDisplacement");
    
    SolutionTypeEnum.Add(MECH_ACCELERATION, "mechAcceleration");
    SolutionTypeEnum.Add(MECH_VELOCITY, "mechVelocity");
    SolutionTypeEnum.Add(MECH_RHS_LOAD, "mechRhsLoad");    
    
    SolutionTypeEnum.Add(MECH_STRESS, "mechStress");
    SolutionTypeEnum.Add(MECH_THERMAL_STRESS, "mechThermalStress");
    SolutionTypeEnum.Add(MECH_STRAIN, "mechStrain");
    SolutionTypeEnum.Add(MECH_THERMAL_STRAIN, "mechThermalStrain");
    SolutionTypeEnum.Add(MECH_STRUCT_INTENSTIY, "mechStructIntensity");
    SolutionTypeEnum.Add(MECH_NORMAL_STRUCT_INTENSITY, "mechNormalStructIntensity");
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


    SolutionTypeEnum.Add(MECH_PSEUDO_DENSITY, "mechPseudoDensity");
    SolutionTypeEnum.Add(PHYSICAL_PSEUDO_DENSITY, "physicalPseudoDensity");
    SolutionTypeEnum.Add(MECH_SHAPE, "mechShape");
    SolutionTypeEnum.Add(MECH_TENSOR_TRACE, "mechTensorTrace");
    SolutionTypeEnum.Add(MECH_TENSOR, "mechTensor");
    SolutionTypeEnum.Add(MECH_TENSOR_HILL_MANDEL, "mechTensorHillMandel");

    //electrostatics / elctric current conduction
    SolutionTypeEnum.Add(ELEC_POTENTIAL, "elecPotential");
    SolutionTypeEnum.Add(ELEC_FIELD_INTENSITY, "elecFieldIntensity");
    SolutionTypeEnum.Add(ELEC_CURRENT_DENSITY, "elecCurrentDensity");
    SolutionTypeEnum.Add(ELEC_NORMAL_CURRENT_DENSITY, "elecNormalCurrentDensity");
    SolutionTypeEnum.Add(ELEC_CURRENT, "elecCurrent");
    SolutionTypeEnum.Add(ELEC_POWER_DENSITY, "elecPowerDensity");
    SolutionTypeEnum.Add(ELEC_POWER, "elecPower");
    SolutionTypeEnum.Add(ELEC_FORCE_VWP, "elecForceVWP");
    SolutionTypeEnum.Add(ELEC_CHARGE, "elecCharge");
    SolutionTypeEnum.Add(ELEC_CHARGE_DENSITY, "elecChargeDensity");
    SolutionTypeEnum.Add(ELEC_FLUX_DENSITY, "elecFluxDensity");
    SolutionTypeEnum.Add(ELEC_ENERGY, "elecEnergy");
    SolutionTypeEnum.Add(ELEC_ENERGY_DENSITY, "elecEnergyDensity");
    SolutionTypeEnum.Add(ELEC_RHS_LOAD, "elecRhsLoad");

    SolutionTypeEnum.Add(ELEC_PSEUDO_POLARIZATION, "elecPseudoPolarization");
    SolutionTypeEnum.Add(ELEC_PHYSICAL_PSEUDO_DENSITY, "elecPhysicalPseudoDensity");
    SolutionTypeEnum.Add(ELEC_TENSOR, "elecTensor");
    SolutionTypeEnum.Add(ELEC_TENSOR_TRACE, "elecTensorTrace");

    //smoothing PDE
    SolutionTypeEnum.Add(SMOOTH_DISPLACEMENT, "smoothDisplacement");
    SolutionTypeEnum.Add(SMOOTH_VELOCITY, "smoothVelocity");
    SolutionTypeEnum.Add(GRID_VELOCITY, "gridVelocity");
    SolutionTypeEnum.Add(SMOOTH_STRAIN, "smoothStrain");
    //acoustics
    SolutionTypeEnum.Add(ACOU_PRESSURE, "acouPressure");
    SolutionTypeEnum.Add(ACOU_ACCELERATION, "acouAcceleration");
    SolutionTypeEnum.Add(ACOU_POTENTIAL, "acouPotential");
    SolutionTypeEnum.Add(ACOU_VELOCITY, "acouVelocity");
    SolutionTypeEnum.Add(ACOU_NORMAL_VELOCITY, "acouNormalVelocity");
    SolutionTypeEnum.Add(ACOU_PRESSURE_DERIV_1, "acouPressureD1");
    SolutionTypeEnum.Add(ACOU_PRESSURE_DERIV_2, "acouPressureD2");
    SolutionTypeEnum.Add(ACOU_FORCE, "acouForce");
    SolutionTypeEnum.Add(ACOU_POTENTIAL_DERIV_1, "acouPotentialD1");
    SolutionTypeEnum.Add(ACOU_POTENTIAL_DERIV_2, "acouPotentialD2");
    SolutionTypeEnum.Add(ACOU_RHS_LOAD, "acouRhsLoad");
    SolutionTypeEnum.Add(ACOU_RHS_LOAD_DENSITY, "acouRhsLoadDensity");
    SolutionTypeEnum.Add(ACOU_RHS_LOADP, "acouRhsLoadP");
    SolutionTypeEnum.Add(ACOU_APE_RHS_LOAD, "apeRhsLoad");
    SolutionTypeEnum.Add(ACOU_VORTEX_RHS_LOAD, "vortexRhsLoad");
    SolutionTypeEnum.Add(ACOU_DIV_LH_TENSOR, "acouDivLighthillTensor");
    SolutionTypeEnum.Add(ACOU_RHSVAL, "acouRHSval");
    SolutionTypeEnum.Add(ACOUSURF_RHSVAL, "acouSurfRHSval");
    SolutionTypeEnum.Add(ACOU_ELEM_SPEED_OF_SOUND,"acouSpeedOfSound");
    SolutionTypeEnum.Add(ACOU_POWERDENSITY, "acouPowerDensity");
    SolutionTypeEnum.Add(ACOU_POWER, "acouPower");
    SolutionTypeEnum.Add(ACOU_INTENSITY, "acouIntensity");
    SolutionTypeEnum.Add(ACOU_NORMAL_INTENSITY, "acouNormalIntensity");
    SolutionTypeEnum.Add(ACOU_SURFINTENSITY, "acouSurfIntensity");
    SolutionTypeEnum.Add(ACOU_POT_ENERGY, "acouPotEnergy");
    SolutionTypeEnum.Add(ACOU_KIN_ENERGY, "acouKinEnergy");
    SolutionTypeEnum.Add(ACOU_PMLAUXVEC,"acouPmlAuxVec");
    SolutionTypeEnum.Add(ACOU_PMLAUXSCALAR, "acouPmlAuxScalar");
    SolutionTypeEnum.Add(ACOU_PSEUDO_DENSITY, "acouPseudoDensity");

    SolutionTypeEnum.Add(ACOU_MIXED_MASS_LOAD, "acouMixedMassLoad");
    SolutionTypeEnum.Add(ACOU_MIXED_MOMENTUM_LOAD, "acouMixedMomentumLoad");
    SolutionTypeEnum.Add(ACOU_LAMB_RHS, "acouLambRhs");

    //water waves
    SolutionTypeEnum.Add(WATER_PRESSURE, "waterPressure");
    SolutionTypeEnum.Add(WATER_PMLAUXVEC,"waterPmlAuxVec");
    SolutionTypeEnum.Add(WATER_PMLAUXSCALAR, "waterPmlAuxScalar");
    SolutionTypeEnum.Add(WATER_RHS_LOAD, "waterRhsLoad");

    //magnetics
    SolutionTypeEnum.Add(MAG_POTENTIAL, "magPotential");
    SolutionTypeEnum.Add(MAG_POTENTIAL_DERIV1, "magPotentialD1");
    SolutionTypeEnum.Add(MAG_SCALAR_POTENTIAL, "magScalarPotential");
    SolutionTypeEnum.Add(MAG_RHS_LOAD, "magRhsLoad");
    
    SolutionTypeEnum.Add(MAG_FLUX_DENSITY, "magFluxDensity");
    SolutionTypeEnum.Add(MAG_FLUX, "magFlux");
    SolutionTypeEnum.Add(MAG_NORMAL_FLUX_DENSITY, "magNormalFluxDensity");
    SolutionTypeEnum.Add(MAG_FIELD_INTENSITY, "magFieldIntensity");
    SolutionTypeEnum.Add(MAG_EDDY_CURRENT_DENSITY, "magEddyCurrentDensity");
    SolutionTypeEnum.Add(MAG_COIL_CURRENT_DENSITY, "magCoilCurrentDensity");
    SolutionTypeEnum.Add(MAG_TOTAL_CURRENT_DENSITY, "magTotalCurrentDensity");
    SolutionTypeEnum.Add(MAG_POTENTIAL_DIV, "magPotentialDiv");
    SolutionTypeEnum.Add(MAG_FORCE_LORENTZ_DENSITY, "magForceLorentzDensity");
    SolutionTypeEnum.Add(MAG_EDDY_POWER_DENSITY, "magEddyPowerDensity");
    SolutionTypeEnum.Add(MAG_ENERGY_DENSITY, "magEnergyDensity");
    SolutionTypeEnum.Add(MAG_CORE_LOSS_DENSITY, "magCoreLossDensity");
    SolutionTypeEnum.Add(MAG_CORE_LOSS, "magCoreLoss");
    
    SolutionTypeEnum.Add(MAG_FORCE_VWP, "magForceVWP");
    SolutionTypeEnum.Add(MAG_FORCE_LORENTZ, "magForceLorentz");
    SolutionTypeEnum.Add(MAG_ENERGY, "magEnergy");
    SolutionTypeEnum.Add(MAG_EDDY_POWER, "magEddyPower");
    SolutionTypeEnum.Add(MAG_EDDY_CURRENT, "magEddyCurrent");
    SolutionTypeEnum.Add(MAG_ELEM_PERMEABILITY, "magElemPermeability");

    // magnetic - coil results
    SolutionTypeEnum.Add(COIL_CURRENT, "coilCurrent");
    SolutionTypeEnum.Add(COIL_CURRENT_DERIV1, "coilCurrentD1");
    
    SolutionTypeEnum.Add(COIL_INDUCED_VOLTAGE, "coilInducedVoltage");
    SolutionTypeEnum.Add(COIL_INDUCTANCE, "coilInductance");
    SolutionTypeEnum.Add(COIL_LINKED_FLUX, "coilLinkedFlux");
    
    //heat conduction
    SolutionTypeEnum.Add(HEAT_TEMPERATURE, "heatTemperature");
    SolutionTypeEnum.Add(HEAT_MEAN_TEMPERATURE, "heatMeanTemperature");
    SolutionTypeEnum.Add(HEAT_TEMPERATURE_D1, "heatTemperatureD1");
    SolutionTypeEnum.Add(HEAT_FLUX_DENSITY, "heatFluxDensity");
    SolutionTypeEnum.Add(HEAT_RHS_LOAD, "heatRhsLoad");
    SolutionTypeEnum.Add(HEAT_SOURCE_DENSITY, "heatSourceDensity");
    //fluidMech
    SolutionTypeEnum.Add(FLUID_FORCE, "fluidForce");
    SolutionTypeEnum.Add(MEAN_FLUIDMECH_VELOCITY, "meanFluidMechVelocity");
    SolutionTypeEnum.Add(MEAN_FLUIDMECH_VELOCITY_NORMAL, "meanFluidMechVelocityNormal");
    SolutionTypeEnum.Add(DIV_MEAN_FLUIDMECH_VELOCITY, "divMeanFluidMechVelocity");
    SolutionTypeEnum.Add(FLUIDMECH_VELOCITY, "fluidMechVelocity");
    SolutionTypeEnum.Add(FLUIDMECH_PRESSURE, "fluidMechPressure");
    SolutionTypeEnum.Add(FLUIDMECH_VELOCITY_DERIV_1, "fluidMechVelocity_deriv1");
    SolutionTypeEnum.Add(FLUIDMECH_PRESSURE_DERIV_1, "fluidMechPressure_deriv1");
    SolutionTypeEnum.Add(FLUIDMECH_VELOCITY_DERIV_2, "fluidMechVelocity_deriv2");
    SolutionTypeEnum.Add(FLUIDMECH_PRESSURE_DERIV_2, "fluidMechPressure_deriv2");

    SolutionTypeEnum.Add(FLUIDMECH_PRESSURE_TIME_DERIV_1, "fluidMechPressure_timeDeriv1");
    SolutionTypeEnum.Add(FLUIDMECH_PRESSURE_TIME_DERIV_2, "fluidMechPressure_timeDeriv2");

    SolutionTypeEnum.Add(FLUIDMECH_FORCE, "fluidMechForce");
    SolutionTypeEnum.Add(FLUIDMECH_DENSITY, "fluidMechDensity");
    SolutionTypeEnum.Add(FLUIDMECH_TKE, "fluidMechTKE");
    SolutionTypeEnum.Add(FLUIDMECH_STRESS, "fluidMechStress");
    SolutionTypeEnum.Add(FLUIDMECH_STRAINRATE, "fluidMechStrainRate");
    SolutionTypeEnum.Add(FLUIDMECH_WVT, "fluidMechWVT");
    SolutionTypeEnum.Add(FLUIDMECH_WVT_DENSITY, "fluidMechWVTDensity");
    SolutionTypeEnum.Add(FLUIDMECH_WVT_PHI, "fluidMechWVTPhi");
    SolutionTypeEnum.Add(FLUIDMECH_WVT_DENSITY_PHI, "fluidMechWVTDensityPhi");
    SolutionTypeEnum.Add(FLUIDMECH_WVT_U1, "fluidMechWVT_u1");
    SolutionTypeEnum.Add(FLUIDMECH_WVT_U2, "fluidMechWVT_u2");
    SolutionTypeEnum.Add(FLUIDMECH_WVT_F, "fluidMechWVT_f");

    SolutionTypeEnum.Add(LAMBDA_K, "lambda_k");

    // TEST PDE
    SolutionTypeEnum.Add(TEST_DOF, "testDof");
    SolutionTypeEnum.Add(TEST_FIELD, "testField");
    SolutionTypeEnum.Add(TEST_RHS_LOAD, "testRhsLoad");

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
    SolutionTypeEnum.Add(OPT_RESULT_10, "optResult_10");
    SolutionTypeEnum.Add(OPT_RESULT_11, "optResult_11");
    SolutionTypeEnum.Add(OPT_RESULT_12, "optResult_12");
    SolutionTypeEnum.Add(OPT_RESULT_13, "optResult_13");
    SolutionTypeEnum.Add(OPT_RESULT_14, "optResult_14");
    SolutionTypeEnum.Add(OPT_RESULT_15, "optResult_15");
    SolutionTypeEnum.Add(OPT_RESULT_16, "optResult_16");
    SolutionTypeEnum.Add(OPT_RESULT_17, "optResult_17");
    SolutionTypeEnum.Add(OPT_RESULT_18, "optResult_18");
    SolutionTypeEnum.Add(OPT_RESULT_19, "optResult_19");
    SolutionTypeEnum.Add(OPT_RESULT_20, "optResult_20");
    // independent
    SolutionTypeEnum.Add(LAGRANGE_MULT, "lagrangeMultiplier");
    SolutionTypeEnum.Add(LAGRANGE_MULT_DERIV_1, "lagrangeMultiplierD1");
    SolutionTypeEnum.Add(LAGRANGE_MULT_DERIV_2, "lagrangeMultiplierD2");
    // evaluates the spacial gradient of the solution at the nodes.
    // common for all PDEs, no unit
    SolutionTypeEnum.Add(GRAD_ACOU_SOLUTION, "gradAcousticSolution"); // independent on acoustic formulation
    SolutionTypeEnum.Add(GRAD_ELEC_POTENTIAL, "gradElecPotential");
    SolutionTypeEnum.Add(ELEM_DENSITY, "density");
    SolutionTypeEnum.Add(PML_DAMP_FACTOR, "pmlDampFactor");

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

    // ==== Initialization of Material Constants ====
    MaterialTypeEnum.Add( NO_MATERIAL, "noMaterial" );
    MaterialTypeEnum.Add( MAG_PERMEABILITY, "Magnetic_permeability" );
    MaterialTypeEnum.Add( MAG_PERMEABILITY_1, "Magnetic_permeability_1" ); 
    MaterialTypeEnum.Add( MAG_PERMEABILITY_2, "Magnetic_permeability_2" ); 
    MaterialTypeEnum.Add( MAG_PERMEABILITY_3, "Magnetic_permeability_3" );
    MaterialTypeEnum.Add( MAG_PERMEABILITYCURVES, "Magnetic_permeability_curves" );
    MaterialTypeEnum.Add( MAG_RELUCTIVITY, "Magnetic_reluctivity" );
    MaterialTypeEnum.Add( MAG_RELUCTIVITY_DERIV, "Magnetic_reluctivity_derivative" );
    MaterialTypeEnum.Add( MAG_CONDUCTIVITY, "Magnetic_Conductivity" ); 
    MaterialTypeEnum.Add( ELEC_PERMITTIVITY, "Electric_Permittivity" );
    MaterialTypeEnum.Add( ELEC_CONDUCTIVITY, "Electric_Conductivity" );
    MaterialTypeEnum.Add( ELEC_CONDUCTIVITY_TENSOR, "Electric_conductivity_tensor" );
    MaterialTypeEnum.Add( MECH_STIFFNESS_TENSOR, "MechanicStiffnessTensor" );
    MaterialTypeEnum.Add( COEFF_STRAIN_IRREVERSIBLE, "Coeff_Strain_Irreversible" ); 
    MaterialTypeEnum.Add( MECH_EMODULUS, "Mechanic_Emodulus" );
    MaterialTypeEnum.Add( MECH_EMODULUS_X, "Mechanic_Emodulus_X" ); 
    MaterialTypeEnum.Add( MECH_EMODULUS_Y, "Mechanic_Emodulus_Y" );
    MaterialTypeEnum.Add( MECH_EMODULUS_Z, "Mechanic_Emodulus_Z" ); 
    MaterialTypeEnum.Add( MECH_POISSON, "Mechanic_PoissonRation" ); 
    MaterialTypeEnum.Add( MECH_POISSON_XY, "Mechanic_PoissonRation_XY" ); 
    MaterialTypeEnum.Add( MECH_POISSON_YZ, "Mechanic_PoissonRation_YZ" ); 
    MaterialTypeEnum.Add( MECH_POISSON_XZ, "Mechanic_PoissonRation_XZ" ); 
    MaterialTypeEnum.Add( MECH_KMODULUS, "Mechanic_Kmodulus" ); 
    MaterialTypeEnum.Add( MECH_GMODULUS, "Mechanic_Gmodulus" ); 
    MaterialTypeEnum.Add( MECH_GMODULUS_YZ, "Mechanic_Gmodulus_YZ" ); 
    MaterialTypeEnum.Add( MECH_GMODULUS_ZX, "Mechanic_Gmodulus_ZX" ); 
    MaterialTypeEnum.Add( MECH_GMODULUS_XY, "Mechanic_Gmodulus_XY" );
    MaterialTypeEnum.Add( MECH_LAME_MU, "Mechanic_LameMu" );
    MaterialTypeEnum.Add( MECH_LAME_LAMBDA, "Mechanic_LameLambda" ); 
    MaterialTypeEnum.Add( MECH_TEC, "Mechanic_ThermalExpCoef" );
    MaterialTypeEnum.Add( MECH_TEC1, "Mechanic_ThermalExpCoef1" );
    MaterialTypeEnum.Add( MECH_TEC2, "Mechanic_ThermalExpCoef2" );
    MaterialTypeEnum.Add( MECH_TEC3, "Mechanic_ThermalExpCoef3" );
    MaterialTypeEnum.Add( MECH_TEC_VECTOR, "Mechanic_TEC_Vector" );
    MaterialTypeEnum.Add( MECH_TEC_VECTORPLANE, "Mechanic_TEC_VectorPlane" );
    MaterialTypeEnum.Add( MECH_TEC_VECTORAXI, "Mechanic_TEC_VectorAxi" );
    MaterialTypeEnum.Add( MECH_VISCOALPHA_VECTOR, "Mechanic_ViscoAlpha_Vec" );
    MaterialTypeEnum.Add( MECH_VISCOG_VECTOR, "Mechanic_ViscoG_Vec" );
    MaterialTypeEnum.Add( MECH_VISCOK_VECTOR, "Mechanic_ViscoK_Vec" );


    MaterialTypeEnum.Add( MECH_STIFFTENSOR_TEC_VECTOR, "Mechanic_StiffTEC_Vector" );
    MaterialTypeEnum.Add( MECH_STIFFTENSOR_TEC_VECTORPLANE, "Mechanic_StiffTEC_VectorPlane" );
    MaterialTypeEnum.Add( MECH_STIFFTENSOR_TEC_VECTORAXI, "Mechanic_StiffTEC_VectorAxi" );

    MaterialTypeEnum.Add( MECH_TEC_REFTEMPERATURE, "refTemperature");
    MaterialTypeEnum.Add( RAYLEIGH_ALPHA, "Rayleigh_Alpha" ); 
    MaterialTypeEnum.Add( RAYLEIGH_BETA, "Rayleigh_Beta" ); 
    MaterialTypeEnum.Add( RAYLEIGH_FREQUENCY, "Rayleigh_Frequency" ); 
    MaterialTypeEnum.Add( LOSS_TANGENS_DELTA, "Loss_TangensDelta" ); 
    MaterialTypeEnum.Add( DENSITY, "Density" );
    MaterialTypeEnum.Add( ACOU_BULK_MODULUS, "AcousticBulkModulus" ); 
    MaterialTypeEnum.Add( ACOU_SOUND_SPEED, "Acoustic_SoundSpeed" );
    MaterialTypeEnum.Add( IMP_HOLE_DIAM, "holeDiam");
    MaterialTypeEnum.Add( IMP_PLATE_THICKNESS, "plateThick");
    MaterialTypeEnum.Add( IMP_END_CORRECTION, "impEndCorrection");
    MaterialTypeEnum.Add( POROSITY, "sigma");
    MaterialTypeEnum.Add( BETA, "beta");
    MaterialTypeEnum.Add( MPP_VOLUME_DEPTH, "mppVolDepth");
    MaterialTypeEnum.Add( FLOW_MACH_NUMBER, "flowMachNr");
    MaterialTypeEnum.Add( ACOU_IMPEDANCE_REAL_VAL, "acouImpedanceRealValue");
    MaterialTypeEnum.Add( ACOU_IMPEDANCE_IMAG_VAL, "acouImpedanceImagValue");
    MaterialTypeEnum.Add( BOVERA, "BoverA" ); 
    MaterialTypeEnum.Add( ACOU_ALPHA, "AcousticAlpha" ); 
    MaterialTypeEnum.Add( FRACTIONAL_ALG, "FractionalAlg" ); 
    MaterialTypeEnum.Add( FRACTIONAL_MEMORY, "FractionalMemory" ); 
    MaterialTypeEnum.Add( FRACTIONAL_INTERPOL, "FractionalInterpol" );
    MaterialTypeEnum.Add( FRACTIONAL_EXPONENT, "FractionalExponent" ); 
    MaterialTypeEnum.Add( HEAT_CONDUCTIVITY, "HeatConductivity" ); 
    MaterialTypeEnum.Add( HEAT_CONDUCTIVITY_TENSOR, "HeatConductivity_Tensor" );
    MaterialTypeEnum.Add( MAGNETOSTRICTION_TENSOR_h, "Magnetostriction_Tensor_h" ); 
    MaterialTypeEnum.Add( HEAT_CAPACITY, "HeatCapacity" ); 
    MaterialTypeEnum.Add( THERMAL_EXPANSION_TENSOR, "thermalExpansion" ); 
    MaterialTypeEnum.Add( DYNAMIC_VISCOSITY, "dynamicViscosity" ); 
    MaterialTypeEnum.Add( KINEMATIC_VISCOSITY, "kinematicViscosity" );
    MaterialTypeEnum.Add( PIEZO_TENSOR, "PiezoTensor" ); 
    MaterialTypeEnum.Add( SPON_POLARIZATION, "sponPolarization" ); 
    MaterialTypeEnum.Add( SPON_STRAIN, "sponStrain" ); 
    MaterialTypeEnum.Add( EFIELD0, "Efield0" ); 
    MaterialTypeEnum.Add( STRESS0, "Stress0" ); 
    MaterialTypeEnum.Add( DCOUPLE0, "dCouple0" );
    MaterialTypeEnum.Add( RATE_CONSTANT, "rateConstant" );
    MaterialTypeEnum.Add( VISCO_PLASTIC_INDEX, "viscoPlasticIndex" ); 
    MaterialTypeEnum.Add( SATURATION_INDEX, "saturationIndex" ); 
    MaterialTypeEnum.Add( SCALE_FORCE_ELEC, "scaleForceElec" ); 
    MaterialTypeEnum.Add( SCALE_FORCE_MECH, "scaleForceMech" );
    MaterialTypeEnum.Add( SCALE_FORCE_COUPLE, "scaleForceCouple" ); 
    MaterialTypeEnum.Add( VOLUME_FRAC_INIT,"volumeFracInit" ); 
    MaterialTypeEnum.Add( MEAN_TEMPERATURE, "meanTemperature" ); 
    MaterialTypeEnum.Add( X_SATURATION, "Xsaturation" ); 
    MaterialTypeEnum.Add( Y_SATURATION, "Ysaturation" ); 
    MaterialTypeEnum.Add( Y_REMANENCE, "Yremanence" );
    MaterialTypeEnum.Add( PREISACH_WEIGHTS, "preisachWeights" ); 
    MaterialTypeEnum.Add( A_JILES, "aJiles" ); 
    MaterialTypeEnum.Add( ALPHA_JILES, "alphaJiles" ); 
    MaterialTypeEnum.Add( K_JILES, "kJiles" );
    MaterialTypeEnum.Add( C_JILES, "cJiles" ); 
    MaterialTypeEnum.Add( P_DIRECTION, "Pdirection" ); 
    MaterialTypeEnum.Add( HYST_MODEL, "hystModel" ); 
    MaterialTypeEnum.Add( NONLIN_COEFFICIENT, "nonLinCoefficient" ); 
    MaterialTypeEnum.Add( NONLIN_DEPENDENCY, "nonLinDependency" );
    MaterialTypeEnum.Add( NONLIN_APPROXIMATION_TYPE, "nonLinApproximationType" );
    MaterialTypeEnum.Add( NONLIN_DATA_NAME, "nonLinDataName" ); 
    MaterialTypeEnum.Add( DATA_ACCURACY, "dataAccuracy" ); 
    MaterialTypeEnum.Add( MAX_APPROX_VAL, "maxApproxVal" ); 
    MaterialTypeEnum.Add( PYROCOEFFICIENT_TENSOR, "Pyrocoefficient_Tensor" ); 
    MaterialTypeEnum.Add( TEST_ALPHA, "Test_Alpha" );
    MaterialTypeEnum.Add( TEST_BETA, "Test_BETA" );
    MaterialTypeEnum.Add( CORE_LOSS, "coreLossPerMass" );

    // ==== Initialization of Matrix Types ====
    feMatrixType.Add( NOTYPE, "none" );
    feMatrixType.Add( SYSTEM, "system" );
    feMatrixType.Add( STIFFNESS, "stiffness");
    feMatrixType.Add( DAMPING, "damping" );
    feMatrixType.Add( CONVECTION, "convection");
    feMatrixType.Add( MASS, "mass" );
    feMatrixType.Add( AUXILIARY, "auxiliary" );

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
  }


  void SetNumberOfThreads(UInt numThreads){
    NUM_CFS_THREADS = numThreads;
  }

  Enum<SolutionType> SolutionTypeEnum;
  Enum<MaterialType> MaterialTypeEnum;
  Enum<ApproxCurveType> ApproxCurveTypeEnum;
  Enum<NonLinMethodType> NonLinMethodTypeEnum;
  Enum<FEMatrixType> feMatrixType;
  UInt MAX_NUM_FE_MATRICES;
  //give default value of 1 for the OMP threads
  UInt NUM_CFS_THREADS = 1;
}



