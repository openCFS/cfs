// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "FlowMaterial.hh"

#include "Domain/CoefFunction/CoefFunction.hh"
#include "Domain/CoefFunction/CoefXpr.hh"
#include "General/Exception.hh"
#include "Materials/BaseMaterial.hh"



namespace CoupledField
{

  // ***********************
  //   Default Constructor
  // ***********************
  FlowMaterial::FlowMaterial(MathParser* mp,
                             CoordSystem * defaultCoosy) 
  : BaseMaterial(FLOW, mp, defaultCoosy)
  {
    //set the allowed material parameters
    isAllowed_.insert( DENSITY );
    isAllowed_.insert( FLUID_DYNAMIC_VISCOSITY );
    isAllowed_.insert( FLUID_KINEMATIC_VISCOSITY );
    isAllowed_.insert( FLUID_BULK_VISCOSITY );
    isAllowed_.insert( FLUID_REF_PRESSURE );
    isAllowed_.insert( FLUID_ADIABATIC_EXPONENT );
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
    

    if (IsSet(FLUID_DYNAMIC_VISCOSITY) && IsSet(FLUID_KINEMATIC_VISCOSITY)) {
      EXCEPTION("Dynamic and kinematic fluid viscosities found in mat file!");
    }
    
    if (IsSet(FLUID_DYNAMIC_VISCOSITY)) {
      dynamicViscosity = GetScalCoefFnc(FLUID_DYNAMIC_VISCOSITY,Global::REAL);
      kinematicViscosity = 
          CoefFunction::Generate(mp_, Global::REAL,CoefXprBinOp(mp_, dynamicViscosity, density,
                                              CoefXpr::OP_DIV ) );

      SetCoefFct( FLUID_KINEMATIC_VISCOSITY, kinematicViscosity );
    }
    else if (IsSet(FLUID_KINEMATIC_VISCOSITY)){
      kinematicViscosity = GetScalCoefFnc(FLUID_KINEMATIC_VISCOSITY,Global::REAL);
      dynamicViscosity = 
          CoefFunction::Generate(mp_, Global::REAL,CoefXprBinOp(mp_, kinematicViscosity, density,
                                              CoefXpr::OP_MULT ) );
                                              
      SetCoefFct( FLUID_DYNAMIC_VISCOSITY, dynamicViscosity );
    }
    else {
      EXCEPTION("No fluid viscosity is specified in the material file!");
    }

    // In the end, generate the tensors from the scalar values
    GenerateTensor( FLUID_DYNAMIC_VISCOSITY, dynamicViscosity );
    GenerateTensor( FLUID_KINEMATIC_VISCOSITY, kinematicViscosity );
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
