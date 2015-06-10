// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <map>
#include <set>
#include <string>
#include <utility>

#include "Materials/BaseMaterial.hh"
#include "AcousticMaterial.hh"
#include "limits.h"
#include "math.h"
#include "stdlib.h"

namespace CoupledField
{

  // ***********************
  //   Default Constructor
  // ***********************
  AcousticMaterial::AcousticMaterial(MathParser* mp,
                                     CoordSystem * defaultCoosy) 
  : BaseMaterial(mp, defaultCoosy) {

    materialDatabaseName_ = "Acoustic";

    //set the allowed material parameters
    isAllowed_.insert( DENSITY );
    isAllowed_.insert( ACOU_BULK_MODULUS );
    isAllowed_.insert( ACOU_SOUND_SPEED );
    isAllowed_.insert( KINEMATIC_VISCOSITY );
    isAllowed_.insert( ACOU_ALPHA );
    isAllowed_.insert( FRACTIONAL_EXPONENT );
    isAllowed_.insert( RAYLEIGH_ALPHA );
    isAllowed_.insert( RAYLEIGH_BETA );
    isAllowed_.insert( RAYLEIGH_FREQUENCY);
    isAllowed_.insert( LOSS_TANGENS_DELTA);
    isAllowed_.insert( BOVERA );
    isAllowed_.insert( IMP_HOLE_DIAM );
    isAllowed_.insert( IMP_PLATE_THICKNESS );
    isAllowed_.insert( IMP_END_CORRECTION );
    isAllowed_.insert( FLOW_MACH_NUMBER );
    isAllowed_.insert( ACOU_IMPEDANCE_REAL_VAL );
    isAllowed_.insert( ACOU_IMPEDANCE_IMAG_VAL );
    isAllowed_.insert( POROSITY );
    isAllowed_.insert( BETA );
  }

  AcousticMaterial::~AcousticMaterial() {


  }



  void AcousticMaterial::SetScalar( Double param, MaterialType matType, 
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
        dataTypeNotAllowed4SetGet ( dataType, msg );
      }

      scalarParams_[matType] = val;
    }
  }

  void AcousticMaterial::SetScalar( const std::string& param, MaterialType matType,
      Global::ComplexPart dataType ) {
    if (  isAllowed_.find( matType ) == isAllowed_.end() ) {
      std::string dim = "string";
      matTypeNotAllowed( matType, dim );
    } else {
      stringParams_[matType] = param;
    }
  }



  void AcousticMaterial::GetScalar( Double& param, MaterialType matType, 
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
  void AcousticMaterial::GetScalar( std::string& param, MaterialType matType,
              Global::ComplexPart dataType ) const {

    stringMap::const_iterator pos;
    pos = stringParams_.find( matType );
    std::string value;

    if ( pos == stringParams_.end() ) {
      std::string dim = "scalar";
      matTypeNotInDataBase( matType, dim );
    }
    else {
      param=pos->second;
    }
  }

}
