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
  std::ostream *trace    = NULL;
  std::ostream *debug    = NULL;
  std::ostream *cla      = NULL;
  std::ostream *memtrace = NULL;
  std::ostream *data     = NULL;
#endif

  // Generate string stream for generation of error messages
  std::stringstream *error = new std::stringstream();

  // Generate string stream for generation of warning messages
  std::stringstream *warning = new std::stringstream();
}


namespace CoupledField {


  bool PrintGridOnly = false;

  LogConfigurator * logConf = new LogConfigurator();

#ifdef PROFILING
  Profiler * profiler = NULL;
#endif  

#ifdef USE_SCRIPTING
  CFSMessenger * messenger = NULL;
#endif

  Domain * domain = NULL;

  Flags * flags=NULL;

  

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
  BaseParamHandler *params = NULL;
  BaseCommandLineHandler *commandLine = NULL;


  // --------------------------------------------
  //  Implementation of enum conversion routines 
  // --------------------------------------------
  // General template function
  template <class TYPE> 
  void String2Enum(const std::string &in, TYPE &out)
  {
    Error("Not implemented", __FILE__, __LINE__);
  }

  template <class TYPE>
  void Enum2String(const TYPE &in, std::string &out)
  {
    Error("Not implemented", __FILE__, __LINE__);
  }
  
  // AnalysisType
  template<>
  void String2Enum<AnalysisType>(const std::string &in, AnalysisType &out) {

    if (in == "static")
      out = STATIC;
    else if (in == "transient")
      out = TRANSIENT;
    else if (in == "harmonic")
      out = HARMONIC;
    else if (in == "eigenFrequency")
      out = EIGENFREQUENCY;
    else if(in == "multiSequence")
      out = MULTI_SEQUENCE;
    else if (in == "paramIdent")
      // since the parameter identification process lives in frequency domain
      out = HARMONIC;
    else if (in=="transientHarmonic")
      out = TRANSIENTHARMONIC;
    else {
      (*error) << "'" << in << "' cannot be converted into item of "
               << "'AnalysisType'!";
      Error( __FILE__, __LINE__ );
    } 
  }

  template<> 
  void Enum2String<AnalysisType>(const AnalysisType &in, std::string &out) {
    switch(in) {
    case STATIC:
      out = "static";
      break; 
    case TRANSIENT:
      out = "transient";
      break;
    case HARMONIC:
      out = "harmonic";
      break;
    case EIGENFREQUENCY:
      out = "eigenFrequency";
      break;
    case MULTI_SEQUENCE:
      out = "multiSequence";
      break;
    case TRANSIENTHARMONIC:
      out = "transientHarmonic";
      break;
    default:  
      Error("No conversion found for your 'AnalysisType'", __FILE__, __LINE__);
    }
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
      (*error) << "'" << in << "' cannot be converted into item of "
               << "'CouplingInputType'!";
      Error( __FILE__, __LINE__);
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
      Error("No conversion found for your 'CouplingInputType'",
            __FILE__, __LINE__);
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
      (*error) << "'" << in << "' cannot be converted into item of "
               << "'IntegrationMethod'!";
      Error( __FILE__, __LINE__);
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
      Error( "No conversion found for your 'IntegrationMethod'",
             __FILE__, __LINE__ );
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
      (*error) << "'" << in << "' cannot be converted into item of "
               << "'BubbleDynType'!";
      Error( __FILE__, __LINE__);
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
      Error( "No conversion found for your 'BubbleDynType'",
             __FILE__, __LINE__ );
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
      (*error) << "'" << in << "' cannot be converted into item of "
               <<"'CouplingOutputType'!";
      Error( __FILE__, __LINE__ );
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
      Error("No conversion found for your 'CouplingOutputType'",
            __FILE__, __LINE__ );
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
      (*error) << "'" << in << "' cannot be converted into item of "
               << "'CouplingRegionType'!";
      Error( __FILE__, __LINE__ );
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
      Error( "No conversion found for your 'CouplingRegionType'",
             __FILE__, __LINE__ );
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
      (*error) << "'" << in << "' cannot be converted into item of "
               << "'NormType'!";
      Error( __FILE__, __LINE__ );
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
      Error("No conversion found for your 'EntryType::ScalarType'"
            , __FILE__, __LINE__);
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
      (*error) << "'" << in << "' cannot be converted into an '"
               << "EntryType::ScalarType' item!";
      Error( __FILE__, __LINE__);
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
      Error("No conversion found for your 'NormType'", __FILE__, __LINE__);
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
    else if (in == "mechEnergy")
      out = MECH_ENERGY;

    //electrostatics
    else if (in == "elecPotential")
      out = ELEC_POTENTIAL;
    else if (in == "elecFieldIntensity")
      out = ELEC_FIELD_INTENSITY;
    else if (in == "elecPolarization")
      out = ELEC_POLARIZATION;
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

    //smoothing PDE
    else if (in == "smoothDisplacement")
      out = SMOOTH_DISPLACEMENT;

    //acoustics
    else if (in == "acouPressure")
      out = ACOU_PRESSURE;
    else if (in == "acouPotential")
      out = ACOU_POTENTIAL;
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
    else if (in == "acouRHSval")
      out = ACOU_RHSVAL;
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
      
    //magnetics
    else if (in == "magPotential")
      out = MAG_POTENTIAL;
    else if (in == "magFluxDensity")
      out = MAG_FLUX_DENSITY;
    else if (in == "magEddyCurrent")
      out = MAG_EDDY_CURRENT;
    else if (in == "magForceVWP")
      out = MAG_FORCE_VWP;
    else if (in == "magForceLorentz")
      out = MAG_FORCE_LORENTZ;
    else if (in == "magEnergy")
      out = MAG_ENERGY;

    //heat conduction
    else if (in == "heatTemperature")
      out = HEAT_TEMPERATURE;

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

    // bubble
    else if (in == "bubbleRadius")
      out = BUBBLE_RADIUS;
    else if (in == "bubbleRadiusD1")
      out = BUBBLE_RADIUS_DERIV_1;
    else if (in == "bubbleValues")
      out = MAG_FLUX_DENSITY;
	           
    else {
      (*error) << "'" << in << "' cannot be converted into item of "
               << "'SolutionType'!";
      Error( __FILE__, __LINE__ );
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

      case MECH_ENERGY:
        out = "mechEnergy";
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
      case ACOU_RHSVAL:
        out = "acouRHSval";
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
 
        //magnetics  
      case MAG_POTENTIAL:
        out = "magPotential";
        break;
      case MAG_FLUX_DENSITY:
        out = "magFluxDensity";
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

        //heat conduction
      case HEAT_TEMPERATURE:
        out = "heatTemperature";
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
      
      default:
        (*error) << "Wrong type of solution or 'SolutionType2String' not "
                 << "implemented for this type of solution";
        Error( __FILE__, __LINE__);
      }
  } 

  // FEType
  template<>
  void Enum2String<FEType>( const FEType &in, std::string &out ) {

    switch(in) {
    case NOFETYPE:
      out = "noFEType";
      break;
    case  LINE:
      out = "line";
      break;
    case  TRIA:
      out = "triangle";
      break;
    case  QUAD:
      out = "quad";
      break;
    case  TET:
      out = "tetraeder";
      break;
    case  HEX:
      out = "hexaeder";
      break;
    case  PYR:
      out = "pyramid";
      break;
    case  WED:
      out = "wedge";
      break;
    default:
      (*error) << "No conversion found for your 'FEType'";
      Error( __FILE__, __LINE__ );
    }
  }

  template<>
  void String2Enum<FEType>( const std::string &in, FEType &out ) {
    (*error) << "String2Enum not implemented for FEType";
    Error( __FILE__, __LINE__ );
  }

  // ComplexFormat
  template<>
  void String2Enum<ComplexFormat>(const std::string &in, ComplexFormat &out) {
    Error( "Not implemented", __FILE__, __LINE__ );
  }

  template<>
  void Enum2String<ComplexFormat>(const ComplexFormat &in, std::string &out) {
    Error( "Not implemented", __FILE__, __LINE__ );
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
      (*error) << "'" << in << "' cannot be converted into an '"
               << "FreqSamplingType' item!";
      Error( __FILE__, __LINE__);
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
      (*error) << "Found no conversion for supplied FreqSamplingType value!";
      Error( __FILE__, __LINE__ );
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
      Error("No conversion found for your 'DataType'", __FILE__, __LINE__);
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
      (*error) << "'" << in << "' cannot be converted into an '"
               << "DataType' item!";
      Error( __FILE__, __LINE__);
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
    case HEAT_CAPACITY:
      out = "HeatCapacity";
      break;
    case DYNAMIC_VISCOSITY:
      out = "dynamicViscosity";
      break;
    case PIEZO_TENSOR:
      out = "PiezoTensor";
      break;
    case E_SATURATION:
      out = "Esaturation";
      break;
    case P_SATURATION:
      out = "Psaturation";
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
    default:  
      Error("No conversion found for your 'DataType'", __FILE__, __LINE__);
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
    else if ( in == "Magnetic_permability" ) {
      out = MAG_PERMEABILITY;
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
    else if ( in == "HeatCapacity" ) {
      out = HEAT_CAPACITY;
    }
    else if ( in == "dynamicViscosity" ) {
      out = DYNAMIC_VISCOSITY;
    }
    else if ( in == "PiezoTensor" ) {
      out = PIEZO_TENSOR;
    }
    else if ( in == "Esaturation" ) {
      out = E_SATURATION;
    }
    else if ( in == "Psaturation" ) {
      out = P_SATURATION;
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
    else {
      Error("No conversion from string to 'MaterialType' found", __FILE__, __LINE__);
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
      Error("No conversion found for your 'SubTensorType'", __FILE__, __LINE__);
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
      Error("No conversion from string to 'SubTensorType' found", 
            __FILE__, __LINE__);
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
      
    default:  
      Error("No conversion found for your 'MaterialClass'", 
            __FILE__, __LINE__);
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
    else {
      (*error) << "'" << in << "' cannot be converted into an '"
               << "MaterialClass' item!";
      Error( __FILE__, __LINE__);
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
      (*error) << "'" << in << "' cannot be converted into an '"
               << "DataType' item!";
      Error( __FILE__, __LINE__);
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
      Error("No conversion found for your 'MaterialClass'", 
            __FILE__, __LINE__);
    }
  }
  
}
