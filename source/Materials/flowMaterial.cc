// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <map>
#include <set>
#include <string>
#include <utility>

#include "General/exception.hh"
#include "Materials/baseMaterial.hh"
#include "flowMaterial.hh"
#include "limits.h"
#include "math.h"
#include "stdlib.h"



namespace CoupledField
{

  // ***********************
  //   Default Constructor
  // ***********************
  FlowMaterial::FlowMaterial() : BaseMaterial() {


    materialDatabaseName_ = "FluidFlow";
    //set the allowed material parameters
    isAllowed_.insert( DENSITY );
    isAllowed_.insert( DYNAMIC_VISCOSITY );
    isAllowed_.insert( KINEMATIC_VISCOSITY );
  }

  FlowMaterial::~FlowMaterial() {


  }

  void FlowMaterial::Finalize() {
  // Trigger calculation of kinematic or dynamic viscosity
   ComputeAllViscosities();
  }
 
  void FlowMaterial::SetScalar( Double param, MaterialType matType, 
				Global::ComplexPart dataType ) {


    //check, if allowed
    if (  isAllowed_.find( matType ) == isAllowed_.end() ) {
      std::string dim = "scalar";
      matTypeNotAllowed( matType, dim );
    }
    else {
      isSet_.insert( matType );

      Complex val;
      if ( dataType == Global::REAL ) {
	val = Complex ( param, 0.0 );
      }
      else if (dataType == Global::IMAG ) {
	val = Complex ( 0.0, param );
	isComplex_.insert( matType );
      }
      else {
	std::string msg = "SetScalar-Double";
	dataTypeNotAllowed4SetGet( dataType, msg );
      }
      
      scalarParams_[matType] = val;
    }
  }


  void FlowMaterial::GetScalar( Double& param, MaterialType matType, 
				Global::ComplexPart dataType ) const {


    scalarMap::const_iterator pos;
    pos = scalarParams_.find( matType );

    if ( pos == scalarParams_.end() ) {
      std::string dim = "scalar";
      matTypeNotInDataBase( matType, dim );
    }
    else {
      Complex val = pos->second;
      if ( dataType == Global::REAL ) {
	param = val.real();
      }
      else if ( dataType == Global::IMAG ) {
	param = val.imag();
      }
      else {
	std::string msg = "GetScalar-Double";
	dataTypeNotAllowed4SetGet( dataType, msg );
      }
    }
  }
  
  void FlowMaterial::ComputeAllViscosities(){

    Double density, dynamicViscosity, kinematicViscosity;
    if (IsSet(DENSITY))
      GetScalar(density,DENSITY,Global::REAL);
    else
      EXCEPTION("No fluid density is specified in the material file!");
    

    if (IsSet(DYNAMIC_VISCOSITY)) {
      GetScalar(dynamicViscosity,DYNAMIC_VISCOSITY,Global::REAL);
      kinematicViscosity=dynamicViscosity/density;
      SetScalar( kinematicViscosity, KINEMATIC_VISCOSITY, Global::REAL );

    }
    else if (IsSet(KINEMATIC_VISCOSITY)){
      GetScalar(kinematicViscosity,KINEMATIC_VISCOSITY,Global::REAL);
      dynamicViscosity=kinematicViscosity*density;
      SetScalar( dynamicViscosity, DYNAMIC_VISCOSITY, Global::REAL );
    }
    else
      EXCEPTION("No fluid viscosity is specified in the material file!");
  }


}
