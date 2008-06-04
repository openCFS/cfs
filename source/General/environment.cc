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

#ifndef INTEGLIB
  std::ostream *debug    = NULL;
  std::ostream *cla      = NULL;
  std::ostream *memtrace = NULL;
#endif

  // Generate string stream for generation of error messages
  std::stringstream *error = new std::stringstream();

  // Generate string stream for generation of warning messages
  std::stringstream *warning = new std::stringstream();
}


namespace CoupledField {


  LogConfigurator * logConf = new LogConfigurator();

#ifdef PROFILING
  Profiler * profiler = NULL;
#endif  

#ifdef USE_SCRIPTING
  CFSMessenger * messenger = NULL;
#endif

  Domain * domain = NULL;


  BaseFE * ptQ1     = NULL;
  BaseFE * ptQ2     = NULL;
  BaseFE * ptTet1   = NULL;
  BaseFE * ptTet2   = NULL;
  BaseFE * ptL1     = NULL;
  BaseFE * ptL2     = NULL;
  BaseFE * ptTr1    = NULL;
  BaseFE * ptTr2    = NULL;
  BaseFE * ptHexa1  = NULL;
  BaseFE * ptHexa2  = NULL;
  BaseFE * ptPyra1  = NULL;
  BaseFE * ptPyra2  = NULL;
  BaseFE * ptWedge1 = NULL;
  BaseFE * ptWedge2 = NULL;

  // Initialisation of some global pointers
  WriteInfo *Info = NULL;
  BaseCommandLineHandler *commandLine = NULL;


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

  // ScalarType
  template<>
  void Enum2String<EntryType::ScalarType>( const EntryType::ScalarType &in,
                                           std::string &out ) {
    
    switch(in) {
    case EntryType::NOENTRYTYPE:
      out = "No Entry Type";
      break;
    case EntryType::DOUBLE:
      out = "Double";
      break;
    case EntryType::COMPLEX:
      out = "Complex";
      break;
    case EntryType::INTEGER:
      out = "Integer";
      break;
    case EntryType::UINT:
      out = "Unsigned Integer";
      break;
    default:
      EXCEPTION("No conversion found for your 'EntryType::ScalarType'" );
    }
  }

  template<>
  void String2Enum<EntryType::ScalarType>( const std::string &in,  
                                           EntryType::ScalarType &out ) {
    
    if ( in == "No Entry Type" ) {
      out = EntryType::NOENTRYTYPE;
    }
    else if ( in == "Double" ) {
      out = EntryType::DOUBLE;
    }
    else if ( in == "Complex" ) {
      out = EntryType::COMPLEX;
    }
    else if ( in == "Integer" ) {
      out = EntryType::INTEGER;
    }
    else if ( in == "Unsigned Integer" ) {
      out = EntryType::UINT;
    }
    else {
      EXCEPTION( "'" << in << "' cannot be converted into an '"
                 << "EntryType::ScalarType' item!" );
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
  template<>
  void String2Enum<SolutionType>(const std::string &in, SolutionType &out) {

    //mechanics
    if (in == "mechDisplacement")
      out = MECH_DISPLACEMENT;
    else if (in == "mechAcceleration")
      out = MECH_ACCELERATION;
    else if (in == "mechVelocity")
      out = MECH_VELOCITY;
    else if (in == "mechForce")
      out = MECH_FORCE;
    else if (in == "mechStress")
      out = MECH_STRESS;
    else if (in == "mechStrain")
      out = MECH_STRAIN;
    else if (in == "mechStrainIrr")
      out = MECH_STRAIN_IRR;
    else if (in == "mechEnergy")
      out = MECH_ENERGY;
    else if (in == "volumeAboveDefSurf" ) 
      out = MECH_DEF_VOLUME;
    else if (in == "mechRhsLoad" ) 
      out = MECH_RHS_LOAD;
    else if (in == "mechPseudoDensity")
      out = MECH_PSEUDO_DENSITY;
    //electrostatics
    else if (in == "elecPotential")
      out = ELEC_POTENTIAL;
    else if (in == "elecFieldIntensity")
      out = ELEC_FIELD_INTENSITY;
    else if (in == "elecPolarization")
      out = ELEC_POLARIZATION;
    else if (in == "elecPseudoPolarization")
      out = ELEC_PSEUDO_POLARIZATION;
    else if (in == "elecForceVWP")
      out = ELEC_FORCE_VWP;
    else if (in == "elecInterfaceForce")
      out = ELEC_INTERFACE_FORCE;
    else if (in == "elecCharge")
      out = ELEC_CHARGE;
    else if (in == "elecFluxDensity")
      out = ELEC_FLUX_DENSITY;
    else if (in == "elecEnergy")
      out = ELEC_ENERGY;
    else if (in == "elecRhsLoad")
      out = ELEC_RHS_LOAD;

    //smoothing PDE
    else if (in == "smoothDisplacement")
      out = SMOOTH_DISPLACEMENT;

    //acoustics
    else if (in == "acouPressure")
      out = ACOU_PRESSURE;
    else if (in == "acouPotential")
      out = ACOU_POTENTIAL;
    else if (in == "acouVelocity")
      out = ACOU_VELOCITY;
    else if (in == "acouPressureD1")
      out = ACOU_PRESSURE_DERIV_1;
    else if (in == "acouPressureD2")
      out = ACOU_PRESSURE_DERIV_2;
    else if (in == "acouForce")
      out = ACOU_FORCE;
    else if (in == "acouPotentialD1")
      out = ACOU_POTENTIAL_DERIV_1;
    else if (in == "acouPotentialD2")
      out = ACOU_POTENTIAL_DERIV_2;
    else if (in == "acouRhsLoad")
      out = ACOU_RHS_LOAD;
    else if (in == "acouSoundSpeed")
      out = ACOU_SOUND_SPEEED;
    else if (in == "acouBubbleRhsVal")
      out = ACOU_BUBBLE_RHS_VAL;
    else if (in == "acouPotNRBC")
      out = ACOU_POT_NRBC;
    else if (in == "nrbcPhi")
      out = NRBC_PHI;
    else if (in == "acouPressureXYZ")
      out = ACOU_PRESSUREXYZ;
    else if (in == "acouPowerDensity")
      out = ACOU_POWERDENSITY;
    else if (in == "acouPower")
      out = ACOU_POWER;
    else if (in == "acouIntensity")
      out = ACOU_INTENSITY;
    else if (in == "acouSurfIntensity")
      out = ACOU_SURFINTENSITY;
      
    //magnetics
    else if (in == "magPotential")
      out = MAG_POTENTIAL;
    else if (in == "magFluxDensity")
      out = MAG_FLUX_DENSITY;
    else if (in == "magHfield")
      out = MAG_HFIELD;
    else if (in == "magEddyCurrent")
      out = MAG_EDDY_CURRENT;
    else if (in == "magForceVWP")
      out = MAG_FORCE_VWP;
    else if (in == "magForceLorentz")
      out = MAG_FORCE_LORENTZ;
    else if (in == "magEnergy")
      out = MAG_ENERGY;
    else if (in == "magEddyPower")
      out = MAG_EDDY_POWER;
    else if (in == "magRhsLoad")
      out = MAG_RHS_LOAD;

    //heat conduction
    else if (in == "heatTemperature")
      out = HEAT_TEMPERATURE;
    else if (in == "heatRhsLoad")
      out = HEAT_RHS_LOAD;

    //mpcci
    else if (in == "fluidForce")
      out = FLUID_FORCE;

    //stokesFluid
    else if (in == "stokesFluidVelPresVort")
      out = STOKESFLUID_VEL_PRES_VORT;
    else if (in == "stokesFluidVelocity")
      out = STOKESFLUID_VELOCITY;
    else if (in == "stokesFluidPressure")
      out = STOKESFLUID_PRESSURE;
    else if (in == "stokesFluidVorticity")
      out = STOKESFLUID_VORTICITY;
    else if (in == "stokesFluidForce")
      out = STOKESFLUID_FORCE;

    //fluidMech
    else if (in == "fluidMechVelocity")
      out = FLUIDMECH_VELOCITY;
    else if (in == "fluidMechPressure")
      out = FLUIDMECH_PRESSURE;
    else if (in == "fluidMechForce")
      out = FLUIDMECH_FORCE;
    else if (in == "fluidMechDensity")
      out = FLUIDMECH_DENSITY;
    else if (in == "fluidMechTKE")
      out = FLUIDMECH_TKE;
    else if (in == "lambda_k")
      out = LAMBDA_K;
    
    // bubble
    else if (in == "bubbleRadius")
      out = BUBBLE_RADIUS;
    else if (in == "bubbleRadiusD1")
      out = BUBBLE_RADIUS_DERIV_1;
    else if (in == "bubbleValues")
      out = MAG_FLUX_DENSITY;

    // the actual result type is given in result descriptions 
    // in the xml file in the optimization element.
    else if (in == "optResult_1")
      out = OPT_RESULT_1;
    else if (in == "optResult_2")
      out = OPT_RESULT_2;
    else if (in == "optResult_3")
      out = OPT_RESULT_3;

    // independent
    else if (in == "LagrangeMultiplier")
      out = LAGRANGE_MULT;
      
    else {
      EXCEPTION( "'" << in << "' cannot be converted into item of "
               << "'SolutionType'!" );
    }
  }

  template<>
  void Enum2String<SolutionType>(const SolutionType &in, std::string &out)
  { 

    switch (in)
      {
        //mechanics
      case MECH_DISPLACEMENT:
        out = "mechDisplacement";
        break;
      case MECH_ACCELERATION:
        out = "mechAcceleration";
        break;
      case MECH_VELOCITY:
        out = "mechVelocity";
        break;
      case MECH_FORCE:
        out = "mechForce";
        break;
      case MECH_STRESS:
        out = "mechStress";
        break;
      case MECH_STRAIN:
        out = "mechStrain";
        break;
      case MECH_STRAIN_IRR:
        out = "mechStrainIrr";
        break;
      case MECH_ENERGY:
        out = "mechEnergy";
        break;
      case MECH_PSEUDO_DENSITY:
        out = "mechPseudoDensity";
        break;

      case MECH_DEF_VOLUME:
        out = "mechDeformedVolume";
        break;

      case MECH_RHS_LOAD:
        out = "mechRhsLoad";
        break;
        
        //electrostatics
      case ELEC_POTENTIAL:
        out = "elecPotential";
        break;
      case ELEC_FIELD_INTENSITY:
        out = "elecFieldIntensity";
        break;
      case ELEC_POLARIZATION:
        out = "elecPolarization";
        break;
      case ELEC_PSEUDO_POLARIZATION:
        out = "elecPseudoPolarization";
        break;
      case ELEC_FORCE_VWP: 
        out = "elecForceVWP";
        break;
      case ELEC_INTERFACE_FORCE:
        out = "elecInterfaceForce";
        break; 
      case ELEC_CHARGE:
        out = "elecCharge";
        break;
      case ELEC_FLUX_DENSITY:
        out = "elecFluxDensity";
        break; 
      case ELEC_ENERGY:
        out = "elecEnergy";
        break;
      case ELEC_RHS_LOAD:
        out = "elecRhsLoad";
        break;
      
        //smoothing PDE  
      case SMOOTH_DISPLACEMENT:
        out = "smoothDisplacement";
        break;
      
        //acoustics
      case ACOU_POTENTIAL:
        out = "acouPotential";
        break;
      case ACOU_PRESSURE:
        out = "acouPressure";
        break;
      case ACOU_VELOCITY:
        out = "acouVelocity";
        break;
      case ACOU_PRESSURE_DERIV_1:
        out = "acouPressureD1";
        break;
      case ACOU_PRESSURE_DERIV_2:
        out = "acouPressureD2";
        break;
      case ACOU_FORCE:
        out = "acouForce";
        break;
      case ACOU_POTENTIAL_DERIV_1:
        out = "acouPotentialD1";
        break;
      case ACOU_POTENTIAL_DERIV_2:
        out = "acouPotentialD2";
        break;
      case ACOU_RHS_LOAD:
        out = "acouRhsLoad";
        break;
      case ACOU_SOUND_SPEEED:
        out = "acouSoundSpeed";
        break;
      case ACOU_BUBBLE_RHS_VAL:
        out = "acouBubbleRhsVal";
        break;
      case ACOU_POT_NRBC:
        out = "acouPotNRBC";
        break;
      case NRBC_PHI:
        out = "nrbcPhi";
        break;
      case ACOU_PRESSUREXYZ:
        out = "acouPressureXYZ";
        break;
      case ACOU_POWERDENSITY:
        out = "acouPowerDensity";
        break;
      case ACOU_POWER:
        out = "acouPower";
        break;
      case ACOU_INTENSITY:
        out = "acouIntensity";
        break;
      case ACOU_SURFINTENSITY:
        out = "acouSurfIntensity";
        break;
 
        //magnetics  
      case MAG_POTENTIAL:
        out = "magPotential";
        break;
      case MAG_FLUX_DENSITY:
        out = "magFluxDensity";
        break;
      case MAG_HFIELD:
        out = "magHfield";
        break;
      case MAG_EDDY_CURRENT:
        out = "magEddyCurrent";
        break;
      case MAG_FORCE_VWP:
        out = "magForceVWP";
        break;
      case MAG_FORCE_LORENTZ:
        out = "magForceLorentz";
        break;
      case MAG_ENERGY:
        out = "magEnergy";
        break;
      case MAG_EDDY_POWER:
        out = "magEddyPower";
        break;
      case MAG_RHS_LOAD:
        out = "magRhsLoad";
        break;
        
        //heat conduction
      case HEAT_TEMPERATURE:
        out = "heatTemperature";
        break;

      case HEAT_RHS_LOAD:
        out = "heatRhsLoad";
        break;

        //mpcci PDE  
      case FLUID_FORCE:
        out = "fluidForce";
        break;

        //stokesFluid
      case STOKESFLUID_VEL_PRES_VORT:
        out = "stokesFluidVelPresVort";
        break;
      case STOKESFLUID_VELOCITY:
        out = "stokesFluidVelocity";
        break;
      case STOKESFLUID_PRESSURE:
        out = "stokesFluidPressure";
        break;
      case STOKESFLUID_VORTICITY:
        out = "stokesFluidVorticity";
        break;
      case STOKESFLUID_FORCE:
        out = "stokesFluidForce";
        break;

        //fluidMech
      case FLUIDMECH_VELOCITY:
        out = "fluidMechVelocity";
        break;
      case FLUIDMECH_PRESSURE:
        out = "fluidMechPressure";
        break;
      case FLUIDMECH_FORCE:
        out = "fluidMechForce";
        break;
      case FLUIDMECH_DENSITY:
        out = "fluidMechDensity";
        break;
      case FLUIDMECH_TKE:
        out = "fluidMechTKE";
        break;
      case LAMBDA_K:
        out = "lambda_k";
        break;
        
        // bubble
      case BUBBLE_RADIUS:
        out = "bubbleRadius";
        break;
      case BUBBLE_RADIUS_DERIV_1:
        out = "bubbleRadiusD1";
        break;
      case BUBBLE_VOLUME_FRAC:
        out = "bubbleVolumeFrac";
        break;
      
      // write design element data from optimization
      case OPT_RESULT_1:
        out = "optResult_1";
        break;
      case OPT_RESULT_2:
        out = "optResult_2";
        break;
      case OPT_RESULT_3:
        out = "optResult_3";
        break;
      
        // independent
      case LAGRANGE_MULT:
        out = "LagrangeMultiplier";
	break;

      default:
        EXCEPTION( "Wrong type of solution or 'SolutionType2String' not "
                   << "implemented for this type of solution" );
      }
  } 

  // FEType
  template<>
  void Enum2String<FEType>( const FEType &in, std::string &out ) {
      out = ELEM_TYPE_NAMES[in];
  }

  template<>
  void String2Enum<FEType>( const std::string &in, FEType &out ) {
    EXCEPTION( "String2Enum not implemented for FEType" );
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
  void Enum2String<DataType>(const DataType &in, std::string &out) {
    switch(in) {
    case COMPLEX:
      out = "Complex";
      break;
    case REAL:
      out = "Real";
      break;
    case IMAG:
      out = "ImaginaryPart";
      break;
    default:  
      EXCEPTION("No conversion found for your 'DataType'");
    }
  }


  template<>
  void String2Enum<DataType>( const std::string &in, DataType &out ) {

    if ( in == "Complex" ) {
      out = COMPLEX;
    }
    else if ( in == "Real" ) {
      out = REAL;
    }
    else if ( in == "ImaginaryPart" ) {
      out = IMAG;
    }
    else {
      EXCEPTION( "'" << in << "' cannot be converted into an '"
                 << "DataType' item!" );
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
    case PIEZO_TENSOR:
      out = "PiezoTensor";
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
  
  }
  
