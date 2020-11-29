#include "ElectroStaticMaterial.hh"

#include "Domain/CoefFunction/CoefFunction.hh"
#include "Materials/Models/Preisach.hh"

namespace CoupledField
{

  // ***********************
  //   Default Constructor
  // ***********************
  ElectroStaticMaterial::ElectroStaticMaterial(MathParser* mp,
                                               CoordSystem * defaultCoosy) 
  : BaseMaterial(ELECTROSTATIC, mp, defaultCoosy)
  {
    //set the allowed material parameters
    isAllowed_.insert( ELEC_PERMITTIVITY_SCALAR );
    isAllowed_.insert( ELEC_PERMITTIVITY_TENSOR );
    isAllowed_.insert( ELEC_PERMITTIVITY_1 );
    isAllowed_.insert( ELEC_PERMITTIVITY_2 );
    isAllowed_.insert( ELEC_PERMITTIVITY_3 );
    isAllowed_.insert( NONLIN_DEPENDENCY );

    isAllowed_.insert( X_SATURATION );
    isAllowed_.insert( Y_SATURATION );
    isAllowed_.insert( Y_REMANENCE );

    isAllowed_.insert( PREISACH_WEIGHTS );
    isAllowed_.insert( PREISACH_WEIGHTS_DIM );
    isAllowed_.insert( PREISACH_WEIGHTS_CONSTVALUE );
    isAllowed_.insert( PREISACH_WEIGHTS_TYPE );
    isAllowed_.insert( PREISACH_WEIGHTS_MUDAT_A );
    isAllowed_.insert( PREISACH_WEIGHTS_MUDAT_H );
    isAllowed_.insert( PREISACH_WEIGHTS_MUDAT_H2 );
    isAllowed_.insert( PREISACH_WEIGHTS_MUDAT_SIGMA );
    isAllowed_.insert( PREISACH_WEIGHTS_MUDAT_SIGMA2 );
    isAllowed_.insert( PREISACH_WEIGHTS_MUDAT_ETA );
    isAllowed_.insert( PREISACH_WEIGHTS_ANHYSTCOUNTINGTOOUTPUTSAT );
    isAllowed_.insert( PREISACH_WEIGHTS_MUDAT_PARAMSFORHALFRANGE );

    isAllowed_.insert( PREISACH_WEIGHTS_TENSOR );
    isAllowed_.insert( PREISACH_WEIGHTS_FOR_MAYERGOYZ_VECTOR );
    isAllowed_.insert( PREISACH_WEIGHTS_ANHYST_ONLY );
    isAllowed_.insert( PREISACH_WEIGHTS_ANHYST_PARAMSFORHALFRANGE );
    isAllowed_.insert( PREISACH_WEIGHTS_ANHYST_D );
    isAllowed_.insert( PREISACH_WEIGHTS_ANHYST_A );
    isAllowed_.insert( PREISACH_WEIGHTS_ANHYST_B );
    isAllowed_.insert( PREISACH_WEIGHTS_ANHYST_C );
    isAllowed_.insert( PREISACH_MAYERGOYZ_NUM_DIR );
    isAllowed_.insert( PREISACH_MAYERGOYZ_ISOTROPIC );
    isAllowed_.insert( PREISACH_MAYERGOYZ_CLIPOUTPUT );

    isAllowed_.insert( MAYERGOYZ_STARTAXIS_X );
    isAllowed_.insert( MAYERGOYZ_STARTAXIS_Y );
    isAllowed_.insert( MAYERGOYZ_STARTAXIS_Z );
    isAllowed_.insert( MAYERGOYZ_STARTAXIS_X_STRAIN );
    isAllowed_.insert( MAYERGOYZ_STARTAXIS_Y_STRAIN );
    isAllowed_.insert( MAYERGOYZ_STARTAXIS_Z_STRAIN );
    isAllowed_.insert( MAYERGOYZ_LOSSPARAM_A );
    isAllowed_.insert( MAYERGOYZ_LOSSPARAM_B );
    isAllowed_.insert( MAYERGOYZ_LOSSPARAM_A_STRAIN );
    isAllowed_.insert( MAYERGOYZ_LOSSPARAM_B_STRAIN );
    
    isAllowed_.insert( MAYERGOYZ_USEABSDPHI );
    isAllowed_.insert( MAYERGOYZ_USEABSDPHI_STRAIN );
    isAllowed_.insert( MAYERGOYZ_NORMALIZEXINEXP );
    isAllowed_.insert( MAYERGOYZ_NORMALIZEXINEXP_STRAIN );
    isAllowed_.insert( MAYERGOYZ_RESTRICTIONOFPSI );
    isAllowed_.insert( MAYERGOYZ_RESTRICTIONOFPSI_STRAIN );
    isAllowed_.insert( MAYERGOYZ_SCALINGOFXINEXP );
    isAllowed_.insert( MAYERGOYZ_SCALINGOFXINEXP_STRAIN );
    
    isAllowed_.insert( PREISACH_PRESCRIBEOUTPUT );
    isAllowed_.insert( PREISACH_SCALEINITIALSTATE );
    isAllowed_.insert( SCALETOSAT );
    isAllowed_.insert( SCALETOSAT_STRAIN );
    isAllowed_.insert( X_SATURATION );
    isAllowed_.insert( Y_SATURATION );
    isAllowed_.insert( Y_REMANENCE );
    isAllowed_.insert( A_JILES );
    isAllowed_.insert( ALPHA_JILES );
    isAllowed_.insert( K_JILES );
    isAllowed_.insert( C_JILES );
    isAllowed_.insert( JILES_TEST );
    isAllowed_.insert( P_DIRECTION );
    isAllowed_.insert( P_DIRECTION_X );
    isAllowed_.insert( P_DIRECTION_Y );
    isAllowed_.insert( P_DIRECTION_Z );
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
    isAllowed_.insert( HYSTERESIS_DIM );
    isAllowed_.insert( HYST_STRAIN_FORM );
    isAllowed_.insert( HYST_IRRSTRAINS );
    isAllowed_.insert( HYST_IRRSTRAIN_C1 );
    isAllowed_.insert( HYST_IRRSTRAIN_C2 );
    isAllowed_.insert( HYST_IRRSTRAIN_C3 );
    isAllowed_.insert( HYST_IRRSTRAIN_CI );
    isAllowed_.insert( HYST_IRRSTRAIN_CI_SIZE );
    isAllowed_.insert( HYST_IRRSTRAIN_D0 );
    isAllowed_.insert( HYST_IRRSTRAIN_D1 );
    isAllowed_.insert( HYST_IRRSTRAIN_SCALETOSAT );
    isAllowed_.insert( HYST_IRRSTRAIN_PARAMSFORHALFRANGE );
    isAllowed_.insert( ROT_RESISTANCE );
    isAllowed_.insert( HYST_MODEL );
    isAllowed_.insert( HYST_COUPLING_DEFINED );
    isAllowed_.insert( X_SATURATION_STRAIN );
    isAllowed_.insert( Y_SATURATION_STRAIN );
    isAllowed_.insert( S_SATURATION );
    isAllowed_.insert( PREISACH_WEIGHTS_FOR_MAYERGOYZ_VECTOR_STRAIN );
    isAllowed_.insert( PREISACH_WEIGHTS_TENSOR_STRAIN );
    isAllowed_.insert( PREISACH_WEIGHTS_MUDAT_A_STRAIN );
    isAllowed_.insert( PREISACH_WEIGHTS_MUDAT_H_STRAIN );
    isAllowed_.insert( PREISACH_WEIGHTS_MUDAT_H2_STRAIN );
    isAllowed_.insert( PREISACH_WEIGHTS_MUDAT_SIGMA_STRAIN );
    isAllowed_.insert( PREISACH_WEIGHTS_MUDAT_SIGMA2_STRAIN );
    isAllowed_.insert( PREISACH_WEIGHTS_MUDAT_ETA_STRAIN );
    isAllowed_.insert( PREISACH_WEIGHTS_ANHYSTCOUNTINGTOOUTPUTSAT_STRAIN );
    isAllowed_.insert( PREISACH_WEIGHTS_MUDAT_PARAMSFORHALFRANGE_STRAIN );
    isAllowed_.insert( PREISACH_WEIGHTS_TYPE_STRAIN );
    isAllowed_.insert( PREISACH_WEIGHTS_CONSTVALUE_STRAIN );
    isAllowed_.insert( PREISACH_WEIGHTS_DIM_STRAIN );
    isAllowed_.insert( PREISACH_WEIGHTS_STRAIN );
    isAllowed_.insert( HYST_TYPE_IS_PREISACH_STRAIN );
    isAllowed_.insert( HYST_TYPE_IS_PREISACH );
    isAllowed_.insert( ANG_DISTANCE_STRAIN );
    isAllowed_.insert( EVAL_VERSION_STRAIN );
    isAllowed_.insert( ROT_RESISTANCE_STRAIN );
    isAllowed_.insert( HYSTERESIS_DIM_STRAIN );
    isAllowed_.insert( HYST_MODEL_STRAIN );
    isAllowed_.insert( S_DIRECTION_Z );
    isAllowed_.insert( S_DIRECTION_Y );
    isAllowed_.insert( S_DIRECTION_X );
    isAllowed_.insert( S_DIRECTION );
    isAllowed_.insert( PREISACH_WEIGHTS_ANHYST_PARAMSFORHALFRANGE_STRAIN );
    isAllowed_.insert( PREISACH_WEIGHTS_ANHYST_D_STRAIN );
    isAllowed_.insert( PREISACH_WEIGHTS_ANHYST_ONLY_STRAIN );
    isAllowed_.insert( PREISACH_WEIGHTS_ANHYST_C_STRAIN );
    isAllowed_.insert( PREISACH_WEIGHTS_ANHYST_B_STRAIN );
    isAllowed_.insert( PREISACH_WEIGHTS_ANHYST_A_STRAIN );
    isAllowed_.insert( ANG_CLIPPING_STRAIN );
    isAllowed_.insert( ANG_RESOLUTION_STRAIN );
    isAllowed_.insert( AMP_RESOLUTION_STRAIN );
    isAllowed_.insert( PREISACH_MAYERGOYZ_NUM_DIR_STRAIN );
    isAllowed_.insert( PREISACH_MAYERGOYZ_ISOTROPIC_STRAIN );
    isAllowed_.insert( PREISACH_MAYERGOYZ_CLIPOUTPUT_STRAIN );
    isAllowed_.insert( IRRSTRAIN_REUSE_P );
    isAllowed_.insert( TRACE_FORCE_CENTRALDIFF );
    isAllowed_.insert( TRACE_FORCE_RETRACING );
    isAllowed_.insert( TRACE_JAC_RESOLUTION );
  }

  ElectroStaticMaterial::~ElectroStaticMaterial() {
  }
  
  void ElectroStaticMaterial::Finalize() {
    MaterialType orthoProps[3] = {
        ELEC_PERMITTIVITY_1,
        ELEC_PERMITTIVITY_2,
        ELEC_PERMITTIVITY_3
    };
    CalcFull3x3Tensor(ELEC_PERMITTIVITY_SCALAR, orthoProps, ELEC_PERMITTIVITY_TENSOR);
  }


  //============================== Hysteresis =====================================
//
//  Double ElectroStaticMaterial::ComputeScalarDiffVal( UInt nrElem, Double Xval ) {
//
//    Double matDiff, eps;
//
//    UInt idx = globalElem2Local_[nrElem];
//    Double Ycurrent = hyst_->computeValueAndUpdate(Xval, idx);
//
//    //    std::cout << "epsDiff: " << " Xval=" << Xval << "  Yval=" << Ycurrent << std::endl;
//
//    //compute differential material parameter
//    Double dX = Xval - Xprevious_[idx];
//    Double dY = Ycurrent -Yprevious_[idx];
//
//    if ( (abs(dY) < 1e-12) || (abs(dX) < 1e-10) ) {
//      GetScalar(eps,ELEC_PERMITTIVITY,Global::REAL);
//      matDiff = eps;
//    }
//    else {
//      matDiff = dY / dX;
//    }
//
//    return matDiff;
//  }
//
//
//  void ElectroStaticMaterial::SetPreviousHystVal( UInt nrElem, Double Xval ) {
//
//    UInt idx = globalElem2Local_[nrElem];
//
//    Xprevious_[idx] = Xval;
//    Yprevious_[idx] = hyst_->computeValueAndUpdate( Xval, idx );
//  }
//
//
//  Double ElectroStaticMaterial::ComputeScalarHystVal( UInt nrElem, Double Xval ) {
//
//    UInt idx    = globalElem2Local_[nrElem];
//    Double Yval = hyst_->computeValueAndUpdate( Xval, idx );
//
//    return Yval;
//  }
}
