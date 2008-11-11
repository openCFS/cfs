// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <limits.h>
#include <string>

#include "flowMaterial.hh"



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
				DataType dataType ) {


    //check, if allowed
    if (  isAllowed_.find( matType ) == isAllowed_.end() ) {
      std::string dim = "scalar";
      matTypeNotAllowed( matType, dim );
    }
    else {
      isSet_.insert( matType );

      Complex val;
      if ( dataType == REAL ) {
	val = Complex ( param, 0.0 );
      }
      else if (dataType == IMAG ) {
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
				DataType dataType ) const {


    scalarMap::const_iterator pos;
    pos = scalarParams_.find( matType );

    if ( pos == scalarParams_.end() ) {
      std::string dim = "scalar";
      matTypeNotInDataBase( matType, dim );
    }
    else {
      Complex val = pos->second;
      if ( dataType == REAL ) {
	param = val.real();
      }
      else if ( dataType == IMAG ) {
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
      GetScalar(density,DENSITY,REAL);
    else
      Error("No fluid density is specified in the material file!",__FILE__,__LINE__);
    

    if (IsSet(DYNAMIC_VISCOSITY)) {
      GetScalar(dynamicViscosity,DYNAMIC_VISCOSITY,REAL);
      kinematicViscosity=dynamicViscosity/density;
      SetScalar( kinematicViscosity, KINEMATIC_VISCOSITY, REAL );

    }
    else if (IsSet(KINEMATIC_VISCOSITY)){
      GetScalar(kinematicViscosity,KINEMATIC_VISCOSITY,REAL);
      dynamicViscosity=kinematicViscosity*density;
      SetScalar( dynamicViscosity, DYNAMIC_VISCOSITY, REAL );
    }
    else
      Error("No fluid viscosity is specified in the material file!",__FILE__,__LINE__);
  }


}
