#include <iostream>
#include <fstream>
#include <string>

#include "environment.hh"
#include "Utils/tools.hh"


// Since OLAS uses a separate namespace for 
// writing out data, two different declarations
// have to be made
namespace OutInfo{
  std::ostream * trace    = NULL;
  std::ostream * debug    = NULL;
  std::ostream * cla      = NULL;
  std::ostream * memtrace = NULL;
  std::ostream * data     = NULL;
}


namespace CoupledField {

  // Generate string stream for generation of error messages
  std::stringstream *error = new std::stringstream();

  // Generate string stream for generation of warning messages
  std::stringstream *warning = new std::stringstream();

  Boolean PrintGridOnly = FALSE;

  Flags * flags=NULL;

  BaseFE * ptQ1     = NULL;
  BaseFE * ptQ2     = NULL;
  BaseFE * ptTet1   = NULL;
  BaseFE * ptL1     = NULL;
  BaseFE * ptL2     = NULL;
  BaseFE * ptTr1    = NULL;
  BaseFE * ptTr2    = NULL;
  BaseFE * ptHexa1  = NULL;
  BaseFE * ptHexa2  = NULL;
  BaseFE * ptPyra1  = NULL;
  BaseFE * ptWedge1 = NULL;
  BaseFE * ptWedge2 = NULL;

  WriteInfo * Info = NULL;
  //  ConfFile * conf           = NULL;
  BaseParamHandler * params = NULL;


  // --------------------------------------------------------------
  //  Bubbledynamic stuff: used to choose between different typs of
  //  bubbledynamic models  
  // --------------------------------------------------------------
  BubbleDynType bubbleDyn;

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
    else if (in == "eigenfrequency")
      out = EIGENFREQUENCY;
    else if(in == "multiSequence")
      out = MULTI_SEQUENCE;
    else if (in == "paramIdent" || in == "multiHarmonic")
      // since the parameter identification process lives in freqeuncy domain
      out = HARMONIC;
    else if(in == "bubbleDynamic")
      out = BUBBLEDYNAMIC;
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
      out = "eigenfrequency";
      break;
    case MULTI_SEQUENCE:
      out = "multiSequence";
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
    else if (in == "acouForce")
      out = ACOU_FORCE;
    else if (in == "acouPotentialD1")
      out = ACOU_POTENTIAL_DERIV_1;
    else if (in == "acouPotentialD2")
      out = ACOU_POTENTIAL_DERIV_2;

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
    //energy

    //mpcci
    else if (in == "fluidForce")
      out = FLUID_FORCE;

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
    case ACOU_FORCE:
      out = "acouForce";
      break;
    case ACOU_POTENTIAL_DERIV_1:
      out = "acouPotentialD1";
      break;
    case ACOU_POTENTIAL_DERIV_2:
      out = "acouPotentialD2";
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
    //mpcci PDE  
    case FLUID_FORCE:
      out = "fluidForce";
      break;
      
    default:
      (*error) << "Wrong type of solution or 'SolutionType2String' not "
	       << "implemented for this type of solution";
      Error( __FILE__, __LINE__);
    }
  } 

  // EQNType
  template<>
  void String2Enum<EQNType>(const std::string &in, EQNType &out) {

    if (in == "scalar")
      out = NODE_SCALAR;
    else if (in == "block")
      out = NODE_BLOCK;
    else if (in == "sclarBlock")
      out = NODE_SCALAR_BLOCK;
    else if (in == "superBlock")
      out = NODE_SUPER_BLOCK;
    else {
      (*error) << "'" << in << "' cannot be converted into item of 'EQNType'!";
      Error( __FILE__, __LINE__);
    }
  }

  template<>
  void Enum2String<EQNType>(const EQNType &in, std::string &out) {

    switch(in) {
    case NODE_SCALAR:
      out = "scalar";
      break; 
    case NODE_BLOCK:
      out = "block";
      break;
    case NODE_SCALAR_BLOCK:
      out = "scalarBlock";
      break;
    case NODE_SUPER_BLOCK:
      out = "superBlock";
      break;
    default:  
      Error("No conversion found for your 'EQNType'", __FILE__, __LINE__);
    }
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

}
