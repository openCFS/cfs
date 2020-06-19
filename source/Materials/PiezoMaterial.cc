// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "PiezoMaterial.hh"

namespace CoupledField
{
  
  // ***********************
  //   Default Constructor
  // ***********************
  PiezoMaterial::PiezoMaterial(MathParser* mp,
          CoordSystem * defaultCoosy) 
  : BaseMaterial(PIEZO, mp, defaultCoosy)
  {
    //set the allowed material parameters
    isAllowed_.insert( PIEZO_TENSOR );
    isAllowed_.insert( X_SATURATION );
    isAllowed_.insert( Y_SATURATION );
    isAllowed_.insert( S_SATURATION );
    isAllowed_.insert( A_JILES );
    isAllowed_.insert( ALPHA_JILES );
    isAllowed_.insert( K_JILES );
    isAllowed_.insert( C_JILES );
    isAllowed_.insert( JILES_TEST );
    isAllowed_.insert( P_DIRECTION );
    isAllowed_.insert( P_DIRECTION_X );
    isAllowed_.insert( P_DIRECTION_Y );
    isAllowed_.insert( P_DIRECTION_Z );
    isAllowed_.insert( HYST_MODEL );
    isAllowed_.insert( HYST_TYPE_IS_PREISACH_STRAIN );
    isAllowed_.insert( HYST_TYPE_IS_PREISACH );
    isAllowed_.insert( HYSTERESIS_DIM );
    isAllowed_.insert( ROT_RESISTANCE );
    isAllowed_.insert( EVAL_VERSION );
    isAllowed_.insert( PRINT_PREISACH );
    isAllowed_.insert( PRINT_PREISACH_RESOLUTION );
    isAllowed_.insert( IS_TESTING );
    isAllowed_.insert( ANG_DISTANCE );
    isAllowed_.insert( ANG_CLIPPING );
    isAllowed_.insert( ANG_RESOLUTION );
    isAllowed_.insert( AMP_RESOLUTION );
    isAllowed_.insert( INITIAL_STATE );
    isAllowed_.insert( INITIAL_STATE_X );
    isAllowed_.insert( INITIAL_STATE_Y );
    isAllowed_.insert( INITIAL_STATE_Z );
    isAllowed_.insert( NONLIN_DEPENDENCY );
    
    // micro-piezoelectric model
    isAllowed_.insert( MEAN_TEMPERATURE );
    isAllowed_.insert( SPON_POLARIZATION );
    isAllowed_.insert( SPON_STRAIN );
    isAllowed_.insert( EFIELD0 );
    isAllowed_.insert( STRESS0 );
    isAllowed_.insert( DCOUPLE0 );
    isAllowed_.insert( DRIVING_FORCE_90 );
    isAllowed_.insert( DRIVING_FORCE_180 );
    isAllowed_.insert( RATE_CONSTANT );
    isAllowed_.insert( VISCO_PLASTIC_INDEX );
    isAllowed_.insert( SATURATION_INDEX );
    isAllowed_.insert( SCALE_FORCE_ELEC );
    isAllowed_.insert( SCALE_FORCE_MECH );
    isAllowed_.insert( SCALE_FORCE_COUPLE );
    isAllowed_.insert( VOLUME_FRAC_INIT );
  }
  
  PiezoMaterial::~PiezoMaterial() {
  }
  
}
