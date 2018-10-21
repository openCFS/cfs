// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <map>
#include <set>
#include <string>
#include <utility>

#include "Domain/CoefFunction/CoefFunction.hh"
#include "Domain/CoefFunction/CoefXpr.hh"
#include "General/Exception.hh"
#include "Materials/BaseMaterial.hh"
#include "FlowMaterial.hh"
#include "limits.h"
#include "math.h"
#include "stdlib.h"



namespace CoupledField
{

  // ***********************
  //   Default Constructor
  // ***********************
  FlowMaterial::FlowMaterial(MathParser* mp,
                             CoordSystem * defaultCoosy) 
: BaseMaterial(mp, defaultCoosy) {


    materialDatabaseName_ = "FluidFlow";
    //set the allowed material parameters
    isAllowed_.insert( DENSITY );
    isAllowed_.insert( DYNAMIC_VISCOSITY );
    isAllowed_.insert( KINEMATIC_VISCOSITY );
    isAllowed_.insert( BULK_VISCOSITY );
    isAllowed_.insert( REF_PRESSURE );
    isAllowed_.insert( REF_TEMPERATURE );
    isAllowed_.insert( ADIABATIC_EXPONENT );
  }

  FlowMaterial::~FlowMaterial() {


  }

  void FlowMaterial::Finalize() {
  // Trigger calculation of kinematic or dynamic viscosity
   ComputeAllViscosities();
  }
 
 
  
  void FlowMaterial::ComputeAllViscosities(){

    PtrCoefFct density, dynamicViscosity, kinematicViscosity;
    if (IsSet(DENSITY))
      density = GetScalCoefFnc(DENSITY,Global::REAL);
    else
      EXCEPTION("No fluid density is specified in the material file!");
    

    if (IsSet(DYNAMIC_VISCOSITY) && IsSet(KINEMATIC_VISCOSITY)) {
      EXCEPTION("Dynamic and kinematic fluid viscosities found in mat file!");
    }
    
    if (IsSet(DYNAMIC_VISCOSITY)) {
      dynamicViscosity = GetScalCoefFnc(DYNAMIC_VISCOSITY,Global::REAL);
      kinematicViscosity = 
          CoefFunction::Generate(mp_, Global::REAL,CoefXprBinOp(mp_, dynamicViscosity, density,
                                              CoefXpr::OP_DIV ) );

      SetCoefFct( KINEMATIC_VISCOSITY, kinematicViscosity );
    }
    else if (IsSet(KINEMATIC_VISCOSITY)){
      kinematicViscosity = GetScalCoefFnc(KINEMATIC_VISCOSITY,Global::REAL);
      dynamicViscosity = 
          CoefFunction::Generate(mp_, Global::REAL,CoefXprBinOp(mp_, kinematicViscosity, density,
                                              CoefXpr::OP_MULT ) );
                                              
      SetCoefFct( DYNAMIC_VISCOSITY, dynamicViscosity );
    }
    else {
      EXCEPTION("No fluid viscosity is specified in the material file!");
    }

    // In the end, generate the tensors from the scalar values
    GenerateTensor( DYNAMIC_VISCOSITY, dynamicViscosity );
    GenerateTensor( KINEMATIC_VISCOSITY, kinematicViscosity );
  }

  void FlowMaterial::GenerateTensor( MaterialType matType, PtrCoefFct scalar ) {
   
    PtrCoefFct tens;
    
    PtrCoefFct z = CoefFunction::Generate(mp_, Global::REAL,"0");
    PtrCoefFct twice = 
        CoefFunction::Generate(mp_, Global::REAL, 
                               CoefXprBinOp(mp_, "2.0", scalar, CoefXpr::OP_MULT ));
    
    StdVector<PtrCoefFct> vals(36);
    vals.Init(z);
    vals[0+0*6] = twice;
    vals[1+1*6] = twice;
    vals[2+2*6] = twice;
    vals[3+3*6] = scalar;
    vals[4+4*6] = scalar;
    vals[5+5*6] = scalar;
    tens = CoefFunction::Generate(mp_, Global::REAL, 6, 6, vals );
    SetCoefFct(matType, tens);
    
    PtrCoefFct check = tensorCoef_[matType];
   }

}
