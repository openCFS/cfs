#include "ElectroMagneticMaterial.hh"

#include "Domain/ElemMapping/Elem.hh"
#include "Domain/ElemMapping/EntityLists.hh"
#include "Utils/mathParser/mathParser.hh"
#include "Utils/SmoothSpline.hh"
#include "Materials/Models/Preisach.hh"
#include "Domain/CoefFunction/CoefFunctionApprox.hh"
#include "Domain/CoefFunction/CoefFunctionCompound.hh"
#include "Domain/CoefFunction/CoefXpr.hh"
#include "Domain/CoefFunction/CoefFunctionConst.hh"

namespace CoupledField
{

  // ***********************
  //   Default Constructor
  // ***********************
  ElectroMagneticMaterial::ElectroMagneticMaterial(MathParser* mp,
                                                   CoordSystem * defaultCoosy,
                                                   bool isDarwin)
  : BaseMaterial(ELECTROMAGNETIC, mp, defaultCoosy)
  {
    this->isDarwin_ = isDarwin;

    //set the allowed material parameters
    isAllowed_.insert( MAG_PERMEABILITY_TENSOR );
    isAllowed_.insert( MAG_PERMEABILITY_SCALAR );
    isAllowed_.insert( MAG_PERMEABILITY_1 );
    isAllowed_.insert( MAG_PERMEABILITY_2 );
    isAllowed_.insert( MAG_PERMEABILITY_3 );
    isAllowed_.insert( MAG_PERMITTIVITY_TENSOR );
    isAllowed_.insert( MAG_PERMITTIVITY_SCALAR );
    isAllowed_.insert( MAG_PERMITTIVITY_1 );
    isAllowed_.insert( MAG_PERMITTIVITY_2 );
    isAllowed_.insert( MAG_PERMITTIVITY_3 );
    isAllowed_.insert( MAG_RELUCTIVITY_TENSOR );
    isAllowed_.insert( MAG_RELUCTIVITY_SCALAR );
    isAllowed_.insert( MAG_RELUCTIVITY_DERIV );
    isAllowed_.insert( MAG_RELUCTIVITY_DERIV_P1 );
    isAllowed_.insert( MAG_ANHYST_DERIV_P1 );
    isAllowed_.insert( MAG_ANHYST_DERIV_P2 );
    isAllowed_.insert( MAG_ANHYST_DERIV_P3 );
    isAllowed_.insert( MAG_ANHYST_DERIV_P4 );        
    isAllowed_.insert( MAG_RELUCTIVITY_DERIV_P2 );
    isAllowed_.insert( MAG_RELUCTIVITY_DERIV_P3 );
    isAllowed_.insert( MAG_RELUCTIVITY_DERIV_P4 );            
    isAllowed_.insert( MAG_CONDUCTIVITY_TENSOR );
    isAllowed_.insert( MAG_CONDUCTIVITY_SCALAR );
    isAllowed_.insert( MAG_CONDUCTIVITY_1 );
    isAllowed_.insert( MAG_CONDUCTIVITY_2 );
    isAllowed_.insert( MAG_CONDUCTIVITY_3 );
    isAllowed_.insert( ELEC_PERMITTIVITY_TENSOR );
    isAllowed_.insert( MAG_CORE_LOSS_PER_MASS );
    isAllowed_.insert( MAG_BH_VALUES );
    isAllowed_.insert( MAG_BH_VALUES_1 );
    isAllowed_.insert( MAG_BH_VALUES_2 );
    isAllowed_.insert( MAG_BH_VALUES_3 );
    isAllowed_.insert( MAG_BH_DATA_ACCURACY );
    isAllowed_.insert( MAG_BH_MAX_APPROX_VAL );
    isAllowed_.insert( DENSITY );
    isAllowed_.insert( NONLIN_DEPENDENCY );

    // -- Energy based vectorhysteresis
    isAllowed_.insert( MAG_PS_EB );
    isAllowed_.insert( MAG_A_EB );
    isAllowed_.insert( MAG_MSAT_PACEJKA_EB );
    isAllowed_.insert( MAG_A_PACEJKA_EB );
    isAllowed_.insert( MAG_B_PACEJKA_EB );
    isAllowed_.insert( MAG_C_PACEJKA_EB );
    isAllowed_.insert( MAG_MU0_EB );
    isAllowed_.insert( MAG_NUMS_EB );//these two are the old implementation but they are still needed for the invEBHyst
    isAllowed_.insert( MAG_CHI_FACTOR_EB ); //these two are the old implementation but they are still needed for the invEBHyst
    isAllowed_.insert( MAG_MSM_AS );
    isAllowed_.insert( MAG_APPROX_TYPE );
    isAllowed_.insert( MAG_MSM_K1 );
    isAllowed_.insert( MAG_MSM_K2 );
    isAllowed_.insert( MAG_MSM_LAMBDA100 );
    isAllowed_.insert( MAG_MSM_LAMBDA111 );
    isAllowed_.insert( MAG_MSM_PS );
    isAllowed_.insert( MAG_JACOBIAN_METHOD_EB );
    isAllowed_.insert( MAG_ANHYST_TYPE_EB );
    isAllowed_.insert( MAG_ANHYST_FORMULA_EB );
    isAllowed_.insert( MAG_PINNING_FORCES_WEIGHTS_EB );
    isAllowed_.insert( MAG_WEIGHTS_FILE_PATH_EB );
    

    // -- inverse Energy based vectorhysteresis
    isAllowed_.insert( MAG_ANHYST_TYPE_INVEB );
    isAllowed_.insert( MAG_ANHYST_FORMULA_INVEB ); 
    isAllowed_.insert( MAG_JS_INVEB );
    isAllowed_.insert( MAG_MS_INVEB );
    isAllowed_.insert( MAG_PA_INVEB );
    isAllowed_.insert( MAG_PB_INVEB );
    isAllowed_.insert( MAG_PC_INVEB );
    isAllowed_.insert( MAG_A_INVEB );
    isAllowed_.insert( MAG_P0_INVEB );
    isAllowed_.insert( MAG_P1_INVEB );
    isAllowed_.insert( MAG_P2_INVEB );
    isAllowed_.insert( MAG_PINNING_FORCES_WEIGHTS_INVEB );
    isAllowed_.insert( MAG_WEIGHTS_FILE_PATH_INVEB );
    isAllowed_.insert( MAG_LOOKUP_TABLE_FILE_INVEB );
    isAllowed_.insert( MAG_JACOBIAN_METHOD_INVEB );



    isAllowed_.insert( PRESCRIBED_MAGNETIZATION );
    isAllowed_.insert( PRESCRIBED_MAGNETIZATION_X );
    isAllowed_.insert( PRESCRIBED_MAGNETIZATION_Y );
    isAllowed_.insert( PRESCRIBED_MAGNETIZATION_Z );
    
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

    // inversion for vec hysteresis
    isAllowed_.insert( MAX_NUM_IT_HYST_INV );
    isAllowed_.insert( VEC_HYST_INV_METHOD );
    isAllowed_.insert( RES_TOL_H_HYST_INV );
    isAllowed_.insert( RES_TOL_B_HYST_INV );
    isAllowed_.insert( RES_TOL_H_HYST_INV_ISREL );
    isAllowed_.insert( RES_TOL_B_HYST_INV_ISREL );
    isAllowed_.insert( ALPHA_REG_HYST_INV );
    isAllowed_.insert( ALPHA_REG_MIN_HYST_INV );
    isAllowed_.insert( ALPHA_REG_MAX_HYST_INV );
    isAllowed_.insert( MAX_NUM_REG_IT_HYST_INV );
    isAllowed_.insert( ALPHA_LS_MIN_HYST_INV );
    isAllowed_.insert( ALPHA_LS_MAX_HYST_INV );
    isAllowed_.insert( MAX_NUM_LS_IT_HYST_INV );
    isAllowed_.insert( STOP_INV_LS_AT_LOCAL_MIN );
    isAllowed_.insert( JAC_RESOLUTION_HYST_INV );
    isAllowed_.insert( JAC_IMPLEMENTATION_HYST_INV );
    isAllowed_.insert( TRUST_LOW_HYST_INV );
    isAllowed_.insert( TRUST_MID_HYST_INV );
    isAllowed_.insert( TRUST_HIGH_HYST_INV );

    isAllowed_.insert( HYST_INV_PROJLM_MU );
    isAllowed_.insert( HYST_INV_PROJLM_RHO );
    isAllowed_.insert( HYST_INV_PROJLM_BETA );
    isAllowed_.insert( HYST_INV_PROJLM_SIGMA );
    isAllowed_.insert( HYST_INV_PROJLM_GAMMA );
    isAllowed_.insert( HYST_INV_PROJLM_TAU );
    isAllowed_.insert( HYST_INV_PROJLM_C );
    isAllowed_.insert( HYST_INV_PROJLM_P );

    isAllowed_.insert( HYST_INV_FP_SAFETYFACTOR );
    isAllowed_.insert( HYST_LOCAL_INVERSION_PRINT_WARNINGS );

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
    isAllowed_.insert( HYST_TYPE_IS_PREISACH_STRAIN );
    isAllowed_.insert( HYST_TYPE_IS_PREISACH );
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

    isAllowed_.insert( MAYERGOYZ_STARTAXIS_X_STRAIN );
    isAllowed_.insert( MAYERGOYZ_STARTAXIS_Y_STRAIN );
    isAllowed_.insert( MAYERGOYZ_STARTAXIS_Z_STRAIN );

    isAllowed_.insert( IRRSTRAIN_REUSE_P );

    isAllowed_.insert( MAGNETOSTRICTION_TENSOR_h_mag );
    isAllowed_.insert( TRACE_FORCE_CENTRALDIFF );
    isAllowed_.insert( TRACE_FORCE_RETRACING );
    isAllowed_.insert( TRACE_JAC_RESOLUTION );
  }


  ElectroMagneticMaterial::~ElectroMagneticMaterial() {
  }
  
  void ElectroMagneticMaterial::Finalize() {
    ComputeFullMuTensor();

    MaterialType orthoProps[3] = {
        MAG_CONDUCTIVITY_1,
        MAG_CONDUCTIVITY_2,
        MAG_CONDUCTIVITY_3
    };
    CalcFull3x3Tensor(MAG_CONDUCTIVITY_SCALAR, orthoProps, MAG_CONDUCTIVITY_TENSOR);

    if(isDarwin_){
      MaterialType orthoProps2[3] = {
          MAG_PERMITTIVITY_1,
          MAG_PERMITTIVITY_2,
          MAG_PERMITTIVITY_3
      };
    CalcFull3x3Tensor(MAG_PERMITTIVITY_SCALAR, orthoProps2, MAG_PERMITTIVITY_TENSOR);
    }
  }

  // Calculate full permeability and reluctivity tensors from scalar values
  void ElectroMagneticMaterial::ComputeFullMuTensor() {

    PtrCoefFct muTensor, relucTensor;
    StdVector<PtrCoefFct> tensorComp(9), relucComp(9);

    // depending on symmetry, calculate full 3x3 permeability tensor
    SymmetryType symType = GetSymmetryType(MAG_PERMEABILITY_TENSOR);

    if (symType == GENERAL) {

      // in this case we have already the full permeability tensor
      muTensor = GetTensorCoefFnc( MAG_PERMEABILITY_TENSOR, FULL, Global::COMPLEX);

      // Now we have the full mu-tensor, so we can invert the matrix
      // and store the reluctivity tensor.
      // Attention: This currently only works with constant expressions!
      if( muTensor->GetDependency() != CoefFunction::CONSTANT ) {
        EXCEPTION( "The magnetic permeability must be constant!");
      } else {
        shared_ptr< CoefFunctionConst<Complex> > reluc(new CoefFunctionConst<Complex>());
        Matrix<Complex> permMat, invMat;
        shared_ptr< CoefFunctionConst<Complex> > permConst =
            dynamic_pointer_cast< CoefFunctionConst<Complex> >(muTensor);
        LocPointMapped lpm;
        permConst->GetTensor( permMat, lpm);
        permMat.Invert( invMat );
        reluc->SetTensor( invMat );
        SetCoefFct( MAG_RELUCTIVITY_TENSOR, reluc );
      }
    }
    else if (symType == ISOTROPIC) {

      PtrCoefFct isoMu = GetScalCoefFnc( MAG_PERMEABILITY_SCALAR, Global::COMPLEX );
      // set diagonal entries
      tensorComp[0] = isoMu;
      tensorComp[4] = isoMu;
      tensorComp[8] = isoMu;
      muTensor = CoefFunction::Generate( mp_, Global::COMPLEX, 3, 3, tensorComp );
      SetCoefFct( MAG_PERMEABILITY_TENSOR, muTensor );

      // Compute reluctivity (scalar variable)
      PtrCoefFct isoRel =
          CoefFunction::Generate(mp_, Global::COMPLEX,
                                 CoefXprUnaryOp(mp_, isoMu, CoefXpr::OP_INV) );

      // Compute reluctivity
      SetCoefFct( MAG_RELUCTIVITY_SCALAR, isoRel );

      // Compute reluctivity (full tensor)
      relucComp[0] = isoRel;
      relucComp[4] = isoRel;
      relucComp[8] = isoRel;
      relucTensor = CoefFunction::Generate( mp_, Global::COMPLEX, 3, 3, relucComp );
      SetCoefFct( MAG_RELUCTIVITY_TENSOR, relucTensor );
    }
    else if (symType == TRANS_ISOTROPIC) {

      PtrCoefFct mu = GetScalCoefFnc( MAG_PERMEABILITY_SCALAR, Global::COMPLEX );
      PtrCoefFct mu3 = GetScalCoefFnc( MAG_PERMEABILITY_3, Global::COMPLEX );
      tensorComp[0] = mu;
      tensorComp[4] = mu;
      tensorComp[8] = mu3;
      muTensor = CoefFunction::Generate( mp_, Global::COMPLEX, 3, 3, tensorComp );
      SetCoefFct( MAG_PERMEABILITY_TENSOR, muTensor );

      // Compute reluctivity
      PtrCoefFct rel, rel3;
      rel = CoefFunction::Generate(mp_, Global::COMPLEX,
                                   CoefXprUnaryOp(mp_, mu, CoefXpr::OP_INV) );
      rel3 = CoefFunction::Generate(mp_, Global::COMPLEX,
                                         CoefXprUnaryOp(mp_, mu3, CoefXpr::OP_INV) );
      relucComp[0] = rel;
      relucComp[4] = rel;
      relucComp[8] = rel3;
      relucTensor = CoefFunction::Generate( mp_, Global::COMPLEX, 3, 3, relucComp );
      SetCoefFct( MAG_RELUCTIVITY_TENSOR, relucTensor );
    }
    else if (symType == ORTHOTROPIC) {

      PtrCoefFct mu1 = GetScalCoefFnc( MAG_PERMEABILITY_1, Global::COMPLEX );
      PtrCoefFct mu2 = GetScalCoefFnc( MAG_PERMEABILITY_2, Global::COMPLEX );
      PtrCoefFct mu3 = GetScalCoefFnc( MAG_PERMEABILITY_3, Global::COMPLEX );
      tensorComp[0] = mu1;
      tensorComp[4] = mu2;
      tensorComp[8] = mu3;
      muTensor = CoefFunction::Generate( mp_, Global::COMPLEX, 3, 3, tensorComp );
      SetCoefFct( MAG_PERMEABILITY_TENSOR, muTensor );

      // Compute reluctivity
      PtrCoefFct rel1, rel2, rel3;
      rel1 = CoefFunction::Generate(mp_, Global::COMPLEX,
                                   CoefXprUnaryOp(mp_, mu1, CoefXpr::OP_INV) );
      rel2 = CoefFunction::Generate(mp_, Global::COMPLEX,
                                         CoefXprUnaryOp(mp_, mu2, CoefXpr::OP_INV) );
      rel3 = CoefFunction::Generate(mp_, Global::COMPLEX,
                                         CoefXprUnaryOp(mp_, mu3, CoefXpr::OP_INV) );
      relucComp[0] = rel1;
      relucComp[4] = rel2;
      relucComp[8] = rel3;
      relucTensor = CoefFunction::Generate( mp_, Global::COMPLEX, 3, 3, relucComp );
      SetCoefFct( MAG_RELUCTIVITY_TENSOR, relucTensor );
    }
    else {
      EXCEPTION( "Calculation of full permeability matrix for symmetryType '"
          << SymmetryTypeEnum.ToString(symType) << "' not implemented!" );
    }

    // Set symmetry type also for reluctivity
    SetSymmetryType( MAG_RELUCTIVITY_TENSOR, symType );
  }

//  void ElectroMagneticMaterial::InitApproxCurves() {
//
//    // check, if we need to approx BH curve
//    if (  needApproxMatCurves_.find( MAG_PERMEABILITY ) != needApproxMatCurves_.end() ) {
//      std::string nlfnc = GetNonlinFileName(MAG_PERMEABILITY);
//      nlinFncBH_ = new SmoothSpline(nlfnc, MAG_PERMEABILITY);
//
//      //get accuracy of approximation
//      Double dataAccuracy;
//      GetScalar( dataAccuracy, DATA_ACCURACY, Global::REAL );
//
//      //get maximal approximation value
//      Double maxApproxVal;
//      GetScalar( maxApproxVal, MAX_APPROX_VAL, Global::REAL );
//
//      nlinFncBH_->SetAccuracy( dataAccuracy );
//      nlinFncBH_->SetMaxY( maxApproxVal );   //maximal value of magnetic induction B
//      nlinFncBH_->CalcBestParameter();
//      nlinFncBH_->CalcApproximation();
//      nlinFncBH_->Print(); 
//    }
//  }

  /*
   * done in base class
  void ElectroMagneticMaterial::InitHyst( UInt numElemSD, shared_ptr<ElemList> actSDList,
                                          bool isInverse, bool computeHystInverse ) {

    isHystInverse_      = isInverse;
    computeHystInverse_ = computeHystInverse;

    std::string val = stringParams_[HYST_MODEL];
    if ( val != "preisach" ) {
      EXCEPTION( "Currently we just support Preisach Hysteresis Model" );
    }
    else {
      isHysteresis_ = true;
      dim_ = 2;

      // we use vector hysteresisi model
      hyst_ = NULL;

      std::cout << "computeHystInverse: " << computeHystInverse_  
          << " isHystInverse: " << isHystInverse_ << std::endl;  

      GetScalar(Xsat_, X_SATURATION, Global::REAL);
      GetScalar(Ysat_, Y_SATURATION, Global::REAL);
      Matrix<Double> weights;
      GetTensor(weights,  PREISACH_WEIGHTS, Global::REAL);
      bool isVirgin = true;

      hyst_ = new Preisach(numElemSD, Xsat_, Ysat_, weights, isVirgin);

      vecHyst_ = new Hysteresis* [dim_];
      vecHyst_[0] = new Preisach(numElemSD, Xsat_, Ysat_, weights, isVirgin);
      vecHyst_[1] = new Preisach(numElemSD, Xsat_, Ysat_, weights, isVirgin);

      // set map: global to local element number
      EntityIterator it = actSDList->GetIterator();
      UInt iel = 0;
      UInt globalElNr;
      for ( it.Begin(); !it.IsEnd(); it++, iel++) {
        globalElNr = it.GetElem()->elemNum;
        globalElem2Local_[globalElNr] = iel;
      }
    }

    //allocate memory for previous results, needed for the
    //effective material parameter formulation
    vecXprevious_.Resize(dim_,numElemSD);
    vecYprevious_.Resize(dim_,numElemSD);
    vecXprevious_.Init();
    vecYprevious_.Init();

    //
    Double startMatdiff;
    GetScalar( startMatdiff, MAG_RELUCTIVITY, Global::REAL );
    matDiffprevious_.Resize( numElemSD );
    matDiffprevious_.Init( startMatdiff );

    if (  computeHystInverse_ ) {
      //for inverse hysteresis
      vecXact_.Resize(dim_,numElemSD);
      vecYact_.Resize(dim_,numElemSD);
      vecXact_.Init();
      vecYact_.Init();
    }
  }
*/
//
//  void ElectroMagneticMaterial::SetPreviousHystVal( UInt nrElem, Vector<Double>& valVec) {
//
//    UInt idx = globalElem2Local_[nrElem];
//
//    if ( isHystInverse_ ) {
//      vecYprevious_[0][idx] = valVec[0];
//      vecYprevious_[1][idx] = valVec[1];
//      vecXprevious_[0][idx] = vecHyst_[0]->computeValueAndUpdate( valVec[0], idx );
//      vecXprevious_[1][idx] = vecHyst_[1]->computeValueAndUpdate( valVec[1], idx );
//    }
//    else if ( computeHystInverse_ ) {
//      Vector<Double> newX(2), newY(2);
//      ComputeInverseScalar( idx, 0, valVec[0], newX[0] );
//      ComputeInverseScalar( idx, 1, valVec[1], newX[1]  );
//
//      //! now perform also an updating
//      newY[0] = vecHyst_[0]->computeValueAndUpdate( newX[0], idx);
//      newY[1] = vecHyst_[1]->computeValueAndUpdate( newX[1], idx);
//
//      Vector<Double> dX(2), dY(2);
//      dX[0] = newX[0] - vecXprevious_[0][idx];
//      dX[1] = newX[1] - vecXprevious_[1][idx];
//
//      dY[0] = newY[0] - vecYprevious_[0][idx];
//      dY[1] = newY[1] - vecYprevious_[1][idx];
//
//      //Double dB = dY[0]*dY[0] + dY[1]*dY[1];
//      Double dB = dY[1]*dY[1];
//      if ( dB > 1e-5) {
//        Double newMatDiff = ComputeMatDiff( dX, dY, idx );
//        matDiffprevious_[idx] = newMatDiff;
//      }
//
//      vecXprevious_[0][idx] = newX[0];
//      vecXprevious_[1][idx] = newX[1];
//      vecYprevious_[0][idx] = newY[0];
//      vecYprevious_[1][idx] = newY[1];
//
//      vecXact_[0][idx] = vecXprevious_[0][idx];
//      vecXact_[1][idx] = vecXprevious_[1][idx];
//      vecYact_[0][idx] = vecYprevious_[0][idx];
//      vecYact_[1][idx] = vecYprevious_[1][idx];
//    }
//    else {
//      vecXprevious_[0][idx] = valVec[0];
//      vecXprevious_[1][idx] = valVec[1];
//      vecYprevious_[0][idx] = vecHyst_[0]->computeValueAndUpdate( valVec[0], idx );
//      vecYprevious_[1][idx] = vecHyst_[1]->computeValueAndUpdate( valVec[1], idx );
//    }
//  }
//
//
//  void ElectroMagneticMaterial::ComputeScalarDiffValues( UInt nrElem,
//                                                         Vector<Double>& valVec,
//                                                         Vector<Double>& scalarValues ) {
//
//    Vector<Double> Ycurrent(dim_);
//    Vector<Double> Xcurrent(dim_);
//
//    UInt idx = globalElem2Local_[nrElem];
//
//    //   std::cout << "elNr=" << nrElem << "  idx=" << idx << std::endl;
//
//    if ( isHystInverse_ ) {
//      Ycurrent[0] = valVec[0];
//      Ycurrent[1] = valVec[1];
//      Xcurrent[0] = vecHyst_[0]->computeValueAndUpdate(valVec[0], idx);
//      Xcurrent[1] = vecHyst_[1]->computeValueAndUpdate(valVec[1], idx);
//     // std::cout << "Ycurrent:\n" << Ycurrent << "\n Xcurrent: \n" << Ycurrent << std::endl;
//    }
//    else if ( computeHystInverse_ ) {
//      ComputeInverseScalar( idx, 0, valVec[0], Xcurrent[0] );
//      ComputeInverseScalar( idx, 1, valVec[1], Xcurrent[1] );
//      Ycurrent[0] = vecHyst_[0]->computeValueAndUpdate(Xcurrent[0], idx);
//      Ycurrent[1] = vecHyst_[1]->computeValueAndUpdate(Xcurrent[1], idx);
//      //       Ycurrent[0] = valVec[0];
//      //       Ycurrent[1] = valVec[1];
//    }
//    else {
//      Xcurrent[0] = valVec[0];
//      Xcurrent[1] = valVec[1];
//      Ycurrent[0] = vecHyst_[0]->computeValueAndUpdate(valVec[0], idx);
//      Ycurrent[1] = vecHyst_[1]->computeValueAndUpdate(valVec[1], idx);
//    }
//
//    //compute differential material parameter
//    Vector<Double> dX(dim_), dY(dim_);
//    dX[0] = Xcurrent[0] - vecXprevious_[0][idx];
//    dX[1] = Xcurrent[1] - vecXprevious_[1][idx];
//
//    dY[0] = Ycurrent[0] - vecYprevious_[0][idx];
//    dY[1] = Ycurrent[1] - vecYprevious_[1][idx];
//
//  //  std::cout << "dX \n" << dX << "  DdY:\n" << dY << std::endl;
//    Double dB = dY[0]*dY[0] + dY[1]*dY[1];
//
//    scalarValues.Init();
//    if ( dB > 1e-10 ) {
//      scalarValues[0] = ( dX[0] * dY[0] + dX[1] * dY[1] ) / dB ;
//      scalarValues[1] = ( dX[0] * dY[1] + dX[1] * dY[0] ) / dB ;
//    }
//  }
//
//
//  Double ElectroMagneticMaterial::ComputeScalarDiffVal( UInt nrElem, Vector<Double>& valVec) {
//
//    // COMPWARNING: unused variable Double matDiff, eps;
//    Vector<Double> Ycurrent(dim_);
//    Vector<Double> Xcurrent(dim_);
//
//    UInt idx = globalElem2Local_[nrElem];
//
//    if ( isHystInverse_ ) {
//      Ycurrent[0] = valVec[0];
//      Ycurrent[1] = valVec[1];
//      Xcurrent[0] = vecHyst_[0]->computeValueAndUpdate(valVec[0], idx);
//      Xcurrent[1] = vecHyst_[1]->computeValueAndUpdate(valVec[1], idx);
//    }
//    else if ( computeHystInverse_ ) {
//      ComputeInverseScalar( idx, 0, valVec[0], Xcurrent[0] );
//      ComputeInverseScalar( idx, 1, valVec[1], Xcurrent[1] );
//      Ycurrent[0] = valVec[0];
//      Ycurrent[1] = valVec[1];
//    }
//    else {
//      Xcurrent[0] = valVec[0];
//      Xcurrent[1] = valVec[1];
//      Ycurrent[0] = vecHyst_[0]->computeValueAndUpdate(valVec[0], idx);
//      Ycurrent[1] = vecHyst_[1]->computeValueAndUpdate(valVec[1], idx);
//    }
//
//    //    std::cout << "Hx=" << Xcurrent[0] << "  Hy= " << Xcurrent[1] << std::endl;
//    //    std::cout << "Bx=" << Ycurrent[0] << "  By= " << Ycurrent[1] << std::endl;
//    //compute differential material parameter
//    Vector<Double> dX(dim_), dY(dim_);
//    dX[0] = Xcurrent[0] - vecXprevious_[0][idx];
//    dX[1] = Xcurrent[1] - vecXprevious_[1][idx];
//
//    dY[0] = Ycurrent[0] - vecYprevious_[0][idx];
//    dY[1] = Ycurrent[1] - vecYprevious_[1][idx];
//    Double newMatDiff = ComputeMatDiff( dX, dY, idx );
//
//    return newMatDiff;
//  }
//
//  Double ElectroMagneticMaterial::ComputeMatDiff( Vector<Double>& dX, Vector<Double>& dY, UInt idx ) {
//
//    Double matDiff, eps;
//    Double dB = dY[0]*dY[0] + dY[1]*dY[1];
//    //Double dB = dY[1]*dY[1];
//
//    if ( dB < 1e-12 ) {
//      matDiff = matDiffprevious_[idx];
//      //matDiff = eps;
//      //std::cout << "Startnu: " << matDiff << std::endl;
//    }
//    else {
//      matDiff = ( dX[0] * dY[0] + dX[1] * dY[1] ) / dB;
//      //matDiff = ( dX[1] * dY[1] ) / dB;
//    }
//
//    //    std::cout << "dB=" << dB  <<  "  dnu=" << matDiff << std::endl << std::endl;
//    //    std::cout << "  dnu=" << matDiff << std::endl << std::endl;
//
//    if ( matDiff < 50.0 ) {
//      GetScalar(eps,MAG_RELUCTIVITY,Global::REAL);
//      matDiff = eps;
//    }
//
//    return matDiff;
//  }
//
//  void ElectroMagneticMaterial::ComputeVectorHystVal( UInt nrElem, Vector<Double>& in,
//                                                      Vector<Double>& out ) {
//
//    UInt idx = globalElem2Local_[nrElem];
//
//    //    std::cout << "elNr=" << nrElem << "  idx=" << idx << std::endl;
//
//    if ( isHystInverse_ ) {
//      out[0] = vecHyst_[0]->computeValueAndUpdate( in[0], idx );
//      out[1] = vecHyst_[1]->computeValueAndUpdate( in[1], idx );
//    }
//    else if ( computeHystInverse_ ) {
//      ComputeInverseScalar( idx, 0, in[0], out[0] );
//      ComputeInverseScalar( idx, 1, in[1], out[1] );
//    }
//    else {
//      out[0] = vecHyst_[0]->computeValueAndUpdate( in[0], idx );
//      out[1] = vecHyst_[1]->computeValueAndUpdate( in[1], idx );
//    }
//  }
//
//  void ElectroMagneticMaterial::GetVectorHystVal( UInt nrElem, Vector<Double>& Val ) {
//
//    UInt idx = globalElem2Local_[nrElem];
//
//    if ( isHystInverse_ ) {
//      Val[0] = vecHyst_[0]->getValue( idx );
//      Val[1] = vecHyst_[1]->getValue( idx );
//    }
//    else if ( computeHystInverse_ ) {
//      Val[0] = vecXprevious_[0][idx];
//      Val[1] = vecXprevious_[1][idx];
//    }
//    else {
//      Val[0] = vecHyst_[0]->getValue( idx );
//      Val[1] = vecHyst_[1]->getValue( idx );
//    }
//  }
//
//  Double ElectroMagneticMaterial::GetScalarHystVal( UInt nrElem ) {
//
//    EXCEPTION("ElectroMagneticMaterial::GetScalarHystVal makes no sense");
//
//    // COMPWARNING: unused variable UInt idx = globalElem2Local_[nrElem];
//    Double Yval = 0.0; // = hyst_->getValue( idx );
//
//    return Yval;
//  }
//
//
//  //====================================== INVERSE HYST =======================================
//
//
//  void ElectroMagneticMaterial::ComputeInverseScalar( UInt idxEl, UInt comp, Double Yin,
//                                                      Double& Xout ) {
//
//
//    Double eps = 1e-3;
//    Double  dH = vecHyst_[0]->GetIncX();
//
//  //  std::cout << "Yin= " << Yin << std::endl;
//
//    if ( ( abs(Yin) + 0.05*Ysat_ ) > Ysat_ ) {
//      Xout = Xsat_;
//    }
//    else {
//      Double Hs, Ho, Hu, Hact, Bs, Bact, dB;
//      bool found = false;
//
//      //compute starting values
//      Hs = vecXact_[comp][idxEl];
//      Bs = vecHyst_[comp]->computeValueAndUpdate( Hs, idxEl, false);
//
//    //  std::cout << "Start Bs: " << Bs << "  Hs=" << Hs <<  std::endl;
//      if  ( abs(Bs - Yin) < eps ) {
//        found = true;
//        Xout  = Hs;
//     //   std::cout << "Direct found " << std::endl;
//      }
//      else if ( (Bs - Yin) > eps ) {
//        Ho = Hs;
//        Hu = Hs;
//     //   std::cout << "Fix Ho= " << Ho << std::endl;
//        do {
//          Hu  -= dH;
//          Bact = vecHyst_[comp]->computeValueAndUpdate( Hu, idxEl, false);
//          dB   = Bact - Yin;
//        } while ( dB > 0 );
//      }
//      else {
//        Hu = Hs;
//        Ho = Hs;
//      //  std::cout << "Fix Hu=  " << Hu << std::endl;
//        do {
//          Ho  += dH;
//          Bact = vecHyst_[comp]->computeValueAndUpdate( Ho, idxEl, false);
//          //  std::cout << "Compute Ho: " << Ho << "  Bact=" << Bact << std::endl;
//          dB   = Bact - Yin;
//        } while ( dB < 0 );
//      }
//
//      if ( found == false ) {
//     //   std::cout << "Do iter: Bin=" << Yin << "  Bs=" << std::endl;
//        do {
//          Hact = ( Ho + Hu ) * 0.5;
//          Bact = vecHyst_[comp]->computeValueAndUpdate( Hact, idxEl, false);
//          dB   = Bact - Yin;
//
//          if ( dB < 0 )
//            Hu = Hact;
//          else
//            Ho = Hact;
//
//       //   std::cout << "newB =" << Bact << "  Hact=" << Hact << "  Ho=" << Ho << "   Hu=" << Hu << std::endl;
//
//        } while ( abs(dB) > eps && abs(Ho-Hu) > abs(Ho)*1e-4 );
//
//        Xout = Hact;
//      }
//    }
//
//    vecXact_[comp][idxEl] = Xout;
//    vecYact_[comp][idxEl] = Yin;
//
//    // update
//    //vecHyst_[comp]->updateMinMaxList( Xout, idxEl );
//
//    //     if ( found )
//    //       std::cout << " Hval = " << Xout << "  Bval=" << Yin << "   Bs=" << Bs << std::endl;
//    //     else
//    //       std::cout << " Hval = " << Xout << "  Bval=" << Yin << "  Bact=" << Bact << std::endl;
//    //   }
//
//  }
//


  PtrCoefFct ElectroMagneticMaterial::GetScalCoefFncNonLin(MaterialType matType,
                                                           Global::ComplexPart matDataType,
                                                           PtrCoefFct fluxCoef,
                                                           PtrCoefFct tempCoef ) {
    //This method allocates the objects handling the nonlinear BH curve; thereby, we allow
    //approximation with smooth splines and analytically defined functions
    //
    //Please note: in the nonlinear bilinear form, we need the reluctivity (=1/permeability)
    //             therefore, we switch between permability and reluctivity quite often
    //             The analytic defined functions in the material file are
    //             reluctivity(magFluxDensity) = nu(B)
    //
    //The core loss factor is also handled here because it needs to be approximated.
	  //Also the temperature-dependent conductivity is processed here via a call to the BaseMaterial

    // Ensure that only MAG_RELUCTIVITY or CORE_LOSS are queried
    if( matType != MAG_RELUCTIVITY_SCALAR &&
        matType != MAG_RELUCTIVITY_DERIV_P1 &&
        matType != MAG_RELUCTIVITY_DERIV_P2 &&
        matType != MAG_RELUCTIVITY_DERIV_P3 &&
        matType != MAG_RELUCTIVITY_DERIV_P4 &&
        matType != MAG_PERMEABILITY_SCALAR &&
        matType != MAG_CORE_LOSS_PER_MASS &&
        matType != MAG_CONDUCTIVITY_SCALAR) {
      EXCEPTION("Scalar nonlinearity for magnetic materials only allowed for MAG_RELUCTIVITY, MAG_PERMEABILITY_SCALAR and CORE_LOSS!"
          << "MAG_RELUCTIVITY_DERIV must be queried using GetTensorCoefFncNonLin.");
    }

    // Ensure that only real-valued parameters are used
    if( matDataType != Global::REAL ) {
      EXCEPTION( "Only real-valued nonlinear parameters are supported");
    }

    PtrCoefFct ret;

    // TODO: There is so much C&P code, please refactor that!

    if( matType == MAG_RELUCTIVITY_SCALAR || matType == MAG_RELUCTIVITY_DERIV_P1 || matType == MAG_RELUCTIVITY_DERIV_P2 || 
        matType == MAG_RELUCTIVITY_DERIV_P3 || matType == MAG_RELUCTIVITY_DERIV_P4){
      // -----------
      // RELUCTIVITY
      // -----------
      // check if material is isotropic or anisotropic
      if( nonlinIsoParams_.find(MAG_PERMEABILITY_SCALAR) != nonlinIsoParams_.end() ) {
        // ---------------------------
        // ISOTROPIC VERSION
        // ---------------------------
        // check, if nonlinear curve was already calculated
        MatDescriptorNl & matNl = nonlinIsoParams_[MAG_PERMEABILITY_SCALAR];

        //Here we really approximate H(B); see book Kaltenbacher, 2nd, 125ff
        if( matNl.approxType == SMOOTH_SPLINES ) {
          // Check, if smooth spline approximation was already created
          // and initialized
          if( !matNl.approxData ) {
            SmoothSpline * sp = new SmoothSpline( matNl.fileName, MAG_PERMEABILITY_SCALAR );
            sp->SetAccuracy( matNl.measAccuracy );
            sp->SetMaxY( matNl.maxVal );
            sp->CalcBestParameter();
            sp->CalcApproximation();
            sp->Print();
            matNl.approxData = sp;
          }

          ApproxData * sp = matNl.approxData;
          // get linear starting value
          Double startVal = 0.0;
          this->GetScalar( startVal, matType, Global::REAL );
          shared_ptr<CoefFunctionApprox> coef( new CoefFunctionApprox());
          coef->Init( startVal, sp, fluxCoef);
          ret = coef;
          
        }
        else if( matNl.approxType == ANALYTIC ) {
          // this is for describing the reluctivity directly in the xml as analytic formula
          // idea: the string from the xml describes a function with the same notation as
          // described in CoefFunctionCompound.hh
          // basically, all occurences of B_R are replaced with the CoefFunction fluxDensAbs
          // note: a good starting value for B->0 works miracles!

          // get Euclidean norm of B
          CoefXprUnaryOp fluxDensAbsOp = CoefXprUnaryOp( mp_, fluxCoef, CoefXpr::OP_NORM );
          PtrCoefFct fluxDensAbs = CoefFunction::Generate( mp_, Global::REAL, fluxDensAbsOp );

          // get function of B
          std::string nuStr = matNl.analyticExpr;
          shared_ptr<CoefFunctionCompound<Double> > nuFnc(new CoefFunctionCompound<Double>(mp_));
          std::map<std::string,PtrCoefFct> symbolsNu;
          symbolsNu["B"] = fluxDensAbs;
          nuStr.insert(0,"( ");
          nuStr.append(" )");
          nuFnc->SetScalar(nuStr,symbolsNu);
          ret = nuFnc;
        }
      }
      else if( nonlinAnisoParams_.find(MAG_PERMEABILITY_SCALAR) != nonlinAnisoParams_.end() ) {
        // ---------------------------
        // ANISOTROPIC VERSION: here we allow for different BH-curves as a function of the angle!
        // ---------------------------
        StdVector<MatDescriptorNl> & matNl = nonlinAnisoParams_[MAG_PERMEABILITY_SCALAR];
        UInt numCurves = matNl.GetSize();
        StdVector<Double> angles(numCurves);
        StdVector<Double> zScalings(numCurves);
        StdVector<shared_ptr<CoefFunction> > approx(numCurves);
        Double startValAveraged = 0.0;

        // Loop over all entries
        for( UInt i = 0; i < matNl.GetSize(); ++i ) {
          MatDescriptorNl & actNl = matNl[i];
          angles[i] = actNl.angle;
          zScalings[i] = actNl.zScaling;
          
          // check for scaling factor == 0
          if (zScalings[i] == 0.0) {
            WARN("z-scaling factor of nonlinear curve for angle " << angles[i] << " is set to zero!\n"
              << " consideration of nonlinear curves will be disabled in z-direction which leads to 2D anisotropic behavior!")
          }              

          //Here we really approximate H(B); see book Kaltenbacher, 2nd, 125ff
          if( actNl.approxType == SMOOTH_SPLINES ) {
            // Check, if smooth spline approximation was already created
            // and initialized
            if( !actNl.approxData ) {
              SmoothSpline * sp = new SmoothSpline( actNl.fileName, MAG_PERMEABILITY_SCALAR );
              sp->SetAccuracy( actNl.measAccuracy );
              sp->SetMaxY( actNl.maxVal );
              sp->CalcBestParameter();
              sp->CalcApproximation();
              sp->Print();
              actNl.approxData = sp;
            }

            ApproxData * sp = actNl.approxData;
            // get linear starting value

            Double startVal;
            this->GetScalar( startVal, matType, Global::REAL );
            shared_ptr<CoefFunctionApprox> coef( new CoefFunctionApprox());
            coef->Init( startVal, sp, fluxCoef);

            //compute an averaged starting value
            startValAveraged += startVal / (Double)numCurves;

            //store in array
            approx[i] = coef;
          }
          else if( actNl.approxType == ANALYTIC ) {
            // this is for describing the reluctivity directly in the xml as analytic formula
            // idea: the string from the xml describes a function with the same notation as
            // described in CoefFunctionCompound.hh
            // basically, all occurences of B_R are replaced with the CoefFunction fluxDensAbs
            // note: a good starting value for B->0 works miracles!

            // get Euclidean norm of B
            CoefXprUnaryOp fluxDensAbsOp = CoefXprUnaryOp( mp_, fluxCoef, CoefXpr::OP_NORM );
            PtrCoefFct fluxDensAbs = CoefFunction::Generate( mp_, Global::REAL, fluxDensAbsOp );

            // get function of B
            std::string nuStr = actNl.analyticExpr;
            shared_ptr<CoefFunctionCompound<Double> > nuFnc(new CoefFunctionCompound<Double>(mp_));
            std::map<std::string,PtrCoefFct> symbolsNu;
            symbolsNu["B"] = fluxDensAbs;

            nuStr.insert(0,"( ");
            nuStr.append(" )");
            nuFnc->SetScalar(nuStr,symbolsNu);

            //compute an averaged starting value directly from the string
            Double B_init = 0.0;
            unsigned int handle = mp_->GetNewHandle();
            mp_->RegisterExternalVar(handle,"B_R",&B_init);
            mp_->SetExpr(handle,nuStr);
            Double nuInit = mp_->Eval(handle);
            startValAveraged += nuInit / (Double)numCurves;

            //store in array
            approx[i] = nuFnc;
          }
        }
        // -------------------------
        // Insertion sort algorithm: we sort the BH-curves starting at smallest
        //                           specified angle
        // ------------------------
        Double compAngle;
        Double compZScaling;
        shared_ptr<CoefFunction> compApprox;
        UInt j;
        for( UInt i = 1; i < numCurves; i++ ) {
          compAngle = angles[i];
          compZScaling = zScalings[i];
          compApprox = approx[i];
          j = i;
          while( ( j > 0 ) && ( angles[j - 1] > compAngle ) ) {
            angles[j] = angles[j - 1];
            zScalings[j] = zScalings[j - 1];
            approx[j] = approx[j - 1];
            j = j - 1;
          }
          angles[j] = compAngle;
          zScalings[j] = compZScaling;
          approx[j] = compApprox;
        }

        // allocate the coef-Function for handling the ansiotropy
        shared_ptr<CoefFunctionApproxAniso> coef( new CoefFunctionApproxAniso());
        coef->Init( startValAveraged, approx, angles, zScalings, fluxCoef );
        baseCoefAniso_ = coef;
        ret = coef;
      }
      else if( nonlinIsoTempDependBHParams_.find(MAG_PERMEABILITY_SCALAR) != nonlinIsoTempDependBHParams_.end() ) {
        // ---------------------------
        // ISOTROPIC TEMP-DEPEND BH CURVES VERSION: here we allow for different BH-curves as a function of temperature!
        // ---------------------------
        StdVector<MatDescriptorNl> & matNl = nonlinIsoTempDependBHParams_[MAG_PERMEABILITY_SCALAR];
        UInt numCurves = matNl.GetSize();
        StdVector<Double> temperatures(numCurves);
        StdVector<shared_ptr<CoefFunction> > approx(numCurves);
        Double startValAveraged = 0.0;

        // Loop over all entries
        for( UInt i = 0; i < matNl.GetSize(); ++i ) {
          MatDescriptorNl & actNl = matNl[i];
          temperatures[i] = actNl.temperature;
          
          //Here we really approximate H(B); see book Kaltenbacher, 2nd, 125ff
          if( actNl.approxType == SMOOTH_SPLINES ) {
            // Check, if smooth spline approximation was already created
            // and initialized
            if( !actNl.approxData ) {
              SmoothSpline * sp = new SmoothSpline( actNl.fileName, MAG_PERMEABILITY_SCALAR );
              sp->SetAccuracy( actNl.measAccuracy );
              sp->SetMaxY( actNl.maxVal );
              sp->CalcBestParameter();
              sp->CalcApproximation();
              sp->Print();
              actNl.approxData = sp;
            }

            ApproxData * sp = actNl.approxData;
            // get linear starting value

            Double startVal;
            this->GetScalar( startVal, matType, Global::REAL );
            shared_ptr<CoefFunctionApprox> coef( new CoefFunctionApprox());
            coef->Init( startVal, sp, fluxCoef);

            //compute an averaged starting value
            startValAveraged += startVal / (Double)numCurves;

            //store in array
            approx[i] = coef;
          }
          else if( actNl.approxType == ANALYTIC ) {
            // this is for describing the reluctivity directly in the xml as analytic formula
            // idea: the string from the xml describes a function with the same notation as
            // described in CoefFunctionCompound.hh
            // basically, all occurences of B_R are replaced with the CoefFunction fluxDensAbs
            // note: a good starting value for B->0 works miracles!

            // get Euclidean norm of B
            CoefXprUnaryOp fluxDensAbsOp = CoefXprUnaryOp( mp_, fluxCoef, CoefXpr::OP_NORM );
            PtrCoefFct fluxDensAbs = CoefFunction::Generate( mp_, Global::REAL, fluxDensAbsOp );

            // get function of B
            std::string nuStr = actNl.analyticExpr;
            shared_ptr<CoefFunctionCompound<Double> > nuFnc(new CoefFunctionCompound<Double>(mp_));
            std::map<std::string,PtrCoefFct> symbolsNu;
            symbolsNu["B"] = fluxDensAbs;

            nuStr.insert(0,"( ");
            nuStr.append(" )");
            nuFnc->SetScalar(nuStr,symbolsNu);

            //compute an averaged starting value directly from the string
            Double B_init = 0.0;
            unsigned int handle = mp_->GetNewHandle();
            mp_->RegisterExternalVar(handle,"B_R",&B_init);
            mp_->SetExpr(handle,nuStr);
            Double nuInit = mp_->Eval(handle);
            startValAveraged += nuInit / (Double)numCurves;

            //store in array
            approx[i] = nuFnc;
          }
        }
        // -------------------------
        // Insertion sort algorithm: we sort the BH-curves starting with the smallest temperature
        // ------------------------
        Double compTemperature;
        shared_ptr<CoefFunction> compApprox;
        UInt j;
        for( UInt i = 1; i < numCurves; i++ ) {
          compTemperature = temperatures[i];
          compApprox = approx[i];
          j = i;
          while( ( j > 0 ) && ( temperatures[j - 1] > compTemperature ) ) {
            temperatures[j] = temperatures[j - 1];
            approx[j] = approx[j - 1];
            j = j - 1;
          }
          temperatures[j] = compTemperature;
          approx[j] = compApprox;
        }

        shared_ptr<CoefFunctionApproxIsotropicTemperatureDependent> coef( new CoefFunctionApproxIsotropicTemperatureDependent());
        // initialize the coef function not only with the flux density but also with the temperature
        coef->Init( startValAveraged, approx, temperatures, fluxCoef, tempCoef ); 
        baseCoefIsoTempDependBH_ = coef;
        ret = coef;
      }

      else {
        EXCEPTION( "No nonlinear definition found for material type '"
           << MaterialTypeEnum.ToString(matType) << "'");
      }

    } else if( matType == MAG_CORE_LOSS_PER_MASS ){
      //-----------
      // CORE_LOSS
      //-----------
      if ( nonlinIsoParams_.find(MAG_CORE_LOSS_PER_MASS) != nonlinIsoParams_.end() ) {
        MatDescriptorNl & matNl = nonlinIsoParams_[MAG_CORE_LOSS_PER_MASS];
        if( matNl.approxType == SMOOTH_SPLINES ) {
          // Check, if smooth spline approximation was already created
          // and initialized
          if( !matNl.approxData ) {
            SmoothSpline * sp = new SmoothSpline( matNl.fileName, MAG_CORE_LOSS_PER_MASS );
            sp->SetAccuracy( matNl.measAccuracy );
            sp->SetMaxY( matNl.maxVal );
            sp->CalcBestParameter();
            sp->CalcApproximation();
            sp->Print();
            matNl.approxData = sp;
          }
          ApproxData * sp = matNl.approxData;
          // get linear starting value
          Double startVal = 0.0;
          this->GetScalar( startVal, matType, Global::REAL );
          shared_ptr<CoefFunctionApprox> coef( new CoefFunctionApprox() );
          coef->Init( startVal, sp, fluxCoef);
          ret = coef;
        } else if( matNl.approxType == LIN_INTERPOLATE ) {
          if ( !matNl.approxData ) {
            LinInterpolate * li = new LinInterpolate( matNl.fileName, MAG_CORE_LOSS_PER_MASS );
            matNl.approxData = li;
          }
          ApproxData * li = matNl.approxData;
          shared_ptr<CoefFunctionApprox> coef( new CoefFunctionApprox() );
          coef->Init( 0.0, li, fluxCoef );
          ret = coef;
        }
      }
      else {
        // since the core loss is an optional parameter (as well as density)
        // we have to guarantee here to return something, otherwise
        // the ResultFunctorIntegrate throws
        // checking for IsSet in the PDE is not enough, which seems odd
        ret = CoefFunction::Generate( mp_, Global::REAL, "0.0" );
      }
    } else if( matType == MAG_CONDUCTIVITY_SCALAR){
    	ret = BaseMaterial::GetScalCoefFncNonLin(matType, matDataType, fluxCoef);
    } else if( matType == MAG_PERMEABILITY_SCALAR){
      // -----------
      // PERMEABILITY
      // -----------
      // check if material is isotropic or anisotropic
      if( nonlinIsoParams_.find(MAG_PERMEABILITY_SCALAR) != nonlinIsoParams_.end() ) {
        // ---------------------------
        // ISOTROPIC VERSION
        // ---------------------------
        // check, if nonlinear curve was already calculated
        MatDescriptorNl & matNl = nonlinIsoParams_[MAG_PERMEABILITY_SCALAR];

        //Here we really approximate H(B); see book Kaltenbacher, 2nd, 125ff
        if( matNl.approxType == SMOOTH_SPLINES ) {
          EXCEPTION("mu(H) is only implemented for the analytic definition of mu(H)");
        }else if( matNl.approxType == ANALYTIC ) {
          // this is for describing the permeability directly in the xml as analytic formula
          // idea: the string from the xml describes a function with the same notation as
          // described in CoefFunctionCompound.hh
          // basically, all occurences of H_R are replaced with the CoefFunction fluxDensAbs
          // ALTHOUGH fluxDens does not mean H, in this mu(H) evaluation it means field intensity!!!
          // note: a good starting value for H->0 works miracles!

          // get Euclidean norm of H
          CoefXprUnaryOp fluxDensAbsOp = CoefXprUnaryOp( mp_, fluxCoef, CoefXpr::OP_NORM );
          PtrCoefFct fluxDensAbs = CoefFunction::Generate( mp_, Global::REAL, fluxDensAbsOp );

          // get function of H
          std::string muStr = matNl.analyticExpr;
          shared_ptr<CoefFunctionCompound<Double> > muFnc(new CoefFunctionCompound<Double>(mp_));
          std::map<std::string,PtrCoefFct> symbolsMu;
          symbolsMu["H"] = fluxDensAbs;
          muStr.insert(0,"( ");
          muStr.append(" )");
          muFnc->SetScalar(muStr,symbolsMu);
          ret = muFnc;
        }
      }else{
        EXCEPTION("mu(H) in the MagneticScalarPotentialPDE is only implemented for the isotropic case");
      }

    }

    return ret;
  }

PtrCoefFct ElectroMagneticMaterial::
GetScalCoefFncNonLin_MagStrict(MaterialType matType,
                               Global::ComplexPart matDataType,
                               PtrCoefFct mechStrain )
{
     //This method allocates the objects handling a nonlinear nu(S) curve; thereby, we allow
     //approximation with smooth splines
     //
     //Please note: in comparison to the nonlinear treatment of the BH curve
     // 		  both the analyitc function as well as the points to be interpolated must
     //		  must specify the reluctivity in dependency on the mechanical strain
     //		  as only a scalar value is returned, we assume the isotropic case 
     //		  for the dependency on S we currently consider the signed maximum value of S (+/- max(abs(S_i)) where +/- depends on sign of max component)

     // Ensure that only MAG_RELUCTIVITY or MAG_RELUCTIVITY_DERIV are queried
    // std::cout << "ElectroMagneticMaterial: GetNonLinFnc" << std::endl;
     if( matType != MAG_RELUCTIVITY_SCALAR  ) {
       EXCEPTION("Scalar Nonlinearity for magnetic materials only allowed for MAG_RELUCTIVITY!");
     }
     
     // Ensure that only real-valued parameters are used
     if( matDataType != Global::REAL ) {
       EXCEPTION( "Only real-valued nonlinear parameters are supported");
     }
     PtrCoefFct ret;
     
     // check if material is isotropic or anisotropic
     if( nonlinIsoParams_.find(MAG_PERMEABILITY_SCALAR) != nonlinIsoParams_.end() ) {
       
       // ---------------------------
       // ISOTROPIC VERSION
       // ---------------------------
       // check, if nonlinear curve was already calculated
       MatDescriptorNl & matNl = nonlinIsoParams_[MAG_PERMEABILITY_SCALAR];

       //Here we approximate nu(S) from data points
       if( matNl.approxType == SMOOTH_SPLINES ) {	 
		EXCEPTION("Currently nu(S) is currently not implemented with splines (did not find approx)");
       }
       else if( matNl.approxType == ANALYTIC ) {

	//	std::cout << "Use analytic nonlin type for reluctivity" << std::endl;
         // this is for describing the reluctivity directly in the xml as analytic formula
         // idea: the string from the xml describes a function with the same notation as
         // described in CoefFunctionCompound.hh
         // basically, all occurences of B_R are replaced with the CoefFunction fluxDensAbs
         // note: a good starting value for B->0 works miracles!

         // get Euclidean norm of S
         CoefXprUnaryOp strainAbsOp = CoefXprUnaryOp( mp_, mechStrain, CoefXpr::OP_NORM );
         PtrCoefFct mechStrainAbs = CoefFunction::Generate( mp_, Global::REAL, strainAbsOp );

         // get function of S
         std::string nuStr = matNl.analyticExpr;
         shared_ptr<CoefFunctionCompound<Double> > nuFnc(new CoefFunctionCompound<Double>(mp_));
         std::map<std::string,PtrCoefFct> symbolsNu;
         symbolsNu["s"] = mechStrain;
         symbolsNu["S"] = mechStrainAbs;

         nuStr.insert(0,"( ");
         nuStr.append(" )");
         nuFnc->SetScalar(nuStr,symbolsNu);
         return(nuFnc);
       }

     } else if( nonlinAnisoParams_.find(MAG_PERMEABILITY_SCALAR) != nonlinAnisoParams_.end() ) {
       
       EXCEPTION("Currently nu(S) is only implemented for the isotropic case");  
     }

     else {
       EXCEPTION( "No nonlinear definition found for material type '"
           << MaterialTypeEnum.ToString(matType) << "'");
     }

     return ret;
   }

  
  PtrCoefFct ElectroMagneticMaterial:: 
  GetScalCoefFncNonLinDerivParam(MaterialType matType,
                                 Global::ComplexPart matDataType,
                                 PtrCoefFct fluxCoef) {
    //This method allocates the objects handling the drivatives of the
    //reluctivity  or magnetic anhysteresis curve w.r.t parameters
    //needed in adjoint magnetic PDE

    // do check
    if( matType != MAG_RELUCTIVITY_DERIV_P1 &&
        matType != MAG_RELUCTIVITY_DERIV_P2 &&
        matType != MAG_RELUCTIVITY_DERIV_P3 &&
        matType != MAG_RELUCTIVITY_DERIV_P4 &&
        matType != MAG_ANHYST_DERIV_P1 &&
        matType != MAG_ANHYST_DERIV_P2 &&
        matType != MAG_ANHYST_DERIV_P3 &&
        matType != MAG_ANHYST_DERIV_P4 
      ) {
      EXCEPTION("Just derivative of Scalar reluctivity w.r.t. paremeters P1 - P4 or of anhysteresis curve w.r.t. P1 - P4 allowed!");
    }

    // Ensure that only real-valued parameters are used
    if ( matDataType != Global::REAL ) {
      EXCEPTION( "Only real-valued nonlinear parameters are supported");
    }

    PtrCoefFct ret;

    // check, if nonlinear curve was already calculated
    MatDescriptorNl & matNl = nonlinIsoParams_[MAG_PERMEABILITY_SCALAR];
    if ( matNl.approxType == ANALYTIC ) {
      // this is for describing the reluctivity directly in the xml as analytic formula
      // basically, all occurences of B_R or H_R are replaced with the CoefFunction fluxDensAbs

      // get Euclidean norm of B or H
      CoefXprUnaryOp fluxDensAbsOp = CoefXprUnaryOp( mp_, fluxCoef, CoefXpr::OP_NORM );
      PtrCoefFct fluxDensAbs = CoefFunction::Generate( mp_, Global::REAL, fluxDensAbsOp );

      shared_ptr<CoefFunctionCompound<Double> > parFnc(new CoefFunctionCompound<Double>(mp_));
      std::map<std::string,PtrCoefFct> symbols;

      // get function of B or H
      std::string fncStr;
      if ( matType == MAG_RELUCTIVITY_DERIV_P1 ) {
        fncStr = matNl.analyticExprDerivP1;
        symbols["B"] = fluxDensAbs;
      } 
      else if ( matType == MAG_RELUCTIVITY_DERIV_P2 ) {
        fncStr = matNl.analyticExprDerivP2;
        symbols["B"] = fluxDensAbs;
      }
      else if ( matType == MAG_RELUCTIVITY_DERIV_P3 ) {
        fncStr = matNl.analyticExprDerivP3;
        symbols["B"] = fluxDensAbs;
      } 
      else if ( matType == MAG_RELUCTIVITY_DERIV_P4 ) {
        fncStr = matNl.analyticExprDerivP4;
        symbols["B"] = fluxDensAbs;
      } 
      else if ( matType == MAG_ANHYST_DERIV_P1) {
        fncStr = matNl.analyticExprDerivP1;
        symbols["H"] = fluxDensAbs;
      }   
      else if ( matType == MAG_ANHYST_DERIV_P2) {
        fncStr = matNl.analyticExprDerivP2;
        symbols["H"] = fluxDensAbs;
      }    
      else if ( matType == MAG_ANHYST_DERIV_P3) {
        fncStr = matNl.analyticExprDerivP3;
        symbols["H"] = fluxDensAbs;
      }
      else if ( matType == MAG_ANHYST_DERIV_P4) {
        fncStr = matNl.analyticExprDerivP4;
        symbols["H"] = fluxDensAbs;
      }         
      else {
        EXCEPTION("matType has to be MAG_RELUCTIVITY_DERIV_P1 - P4 or MAG_ANHYST_DERIV_P1 - P2");
      }
      
      fncStr.insert(0,"( ");
      fncStr.append(" )");
      parFnc->SetScalar(fncStr,symbols);
      ret = parFnc;
    } else {
      EXCEPTION("Just ANALYTIC functions are allowed")
    }
    return ret;
  }


  PtrCoefFct ElectroMagneticMaterial::
  GetTensorCoefFncNonLin( MaterialType matType,
                          SubTensorType type,
                          Global::ComplexPart matDataType,
                          PtrCoefFct dependency )
  {
    //
    //This method allocates the objects handling the derivative of the reluctivity w.r.t.
    //the magnetic flux density ( nu'(B) ); therefore it is called to bulid up the nonlinear
    //bilinear form for the geometric stiffness matrix
    //
    //Please note: in the nonlinear bilinear form, we need the derivative of the
    //             reluctivity (=1/permeability); therefore, we switch between
    //             permability and reluctivity quite often
    //             The analytic defined functions in the material file diretcly define
    //             the derivative of the reluctivity as a function of mag. flux density!
    //

       // Ensure that only MAG_RELUCTIVITY or MAG_RELUCTIVITY_DERIV are queried
       if( matType != MAG_RELUCTIVITY_TENSOR && matType != MAG_RELUCTIVITY_DERIV ) {
         EXCEPTION("Nonlinearity for magnetic materials only allowed for "
             << "MAG_RELUCTIVITY_SCALAR or MAG_RELUCTIVITY_DERIV" );
       }
       
       // Ensure that only real-valued parameters are used
       if( matDataType != Global::REAL ) {
         EXCEPTION( "Only real-valued nonlinear parameters are supported");
       }
       PtrCoefFct ret;
       
       UInt dimDMat = (type == FULL) ? 3 : 2;
       
       // check if material is isotropic or anisotropic
       if( nonlinIsoParams_.find(MAG_PERMEABILITY_SCALAR) != nonlinIsoParams_.end() ) {
         
         // Check, if MAG_RELUCTIVITY is queried
         if( matType == MAG_RELUCTIVITY_SCALAR ) {
           EXCEPTION("An isotropic nonlinear MAG_RELUCTIVITY_SCALAR must be queried using "\
                     "GetScalCoefFncNonLin");
         }
         
         // ---------------------------
         // ISOTROPIC VERSION
         // ---------------------------
         // check, if nonlinear curve was already calculated
         MatDescriptorNl & matNl = nonlinIsoParams_[MAG_PERMEABILITY_SCALAR];

         if( matNl.approxType == SMOOTH_SPLINES ) {

           //Here we really approximate H(B); see book Kaltenbacher, 2nd, 125ff
           if( !matNl.approxData ) {
             SmoothSpline * sp = new SmoothSpline( matNl.fileName, MAG_PERMEABILITY_SCALAR );
             sp->SetAccuracy( matNl.measAccuracy );
             sp->SetMaxY( matNl.maxVal );
             sp->CalcBestParameter();
             sp->CalcApproximation();
             sp->Print();
             matNl.approxData = sp;
           }

           //now we need the object "CoefFunctionApproxDeriv", which returns by
           //calling "B^T [ e_B^T * nu' * |B| * e_B] B", so it is a tensor!!
           ApproxData * sp = matNl.approxData;
           shared_ptr<CoefFunctionApproxDeriv> coef( new CoefFunctionApproxDeriv());
           coef->Init( sp, dimDMat, dependency );
           ret = coef;

         } else if( matNl.approxType == ANALYTIC ) {
           //Here, we obtain " nu' " be evaluating the analytical defined function in
           //the material file, and then we have to build the tensor due to
           // "B^T [ e_B^T * nu' * |B| * e_B] B"
           // dperchto: should be [ nu' * |B| * e_B * e_B^T ] otherwise it would be scalar
           //                                   3x1   1x3
           // implemented as [ nu' / |B| * B * B^T ] which needs one mathematical operation less

           // get Euclidean norm of B
           CoefXprUnaryOp fluxDensAbsOp = CoefXprUnaryOp( mp_, dependency, CoefXpr::OP_NORM );
           PtrCoefFct fluxDensAbs = CoefFunction::Generate( mp_, Global::REAL, fluxDensAbsOp );

           // get function of B
           std::string dnudBStr = matNl.analyticExprDeriv;
           shared_ptr<CoefFunctionCompound<Double> > scalFnc(new CoefFunctionCompound<Double>(mp_));
           std::map<std::string,PtrCoefFct> symbolsNu;
           symbolsNu["B"] = fluxDensAbs;

           // avoid divisions by zero, for very small B a fixed-step iteration is performed
           dnudBStr.insert(0,"( ( B_R lt 1e-6 ) ? ( 0.0 ) : ( ( ");
           dnudBStr.append(" ) / B_R ) )");
           scalFnc->SetScalar(dnudBStr,symbolsNu);

           shared_ptr<CoefFunctionCompound<Double> > dnudBTens(new CoefFunctionCompound<Double>(mp_));
           std::map<std::string,PtrCoefFct> symbolsTens;
           symbolsTens["a"] = dependency;
           StdVector<std::string> bbStr;
           if ( dimDMat == 3 ) {
             bbStr = "( a_0_R * a_0_R )" , "( a_0_R * a_1_R )" , "( a_0_R * a_2_R )" ,
                     "( a_0_R * a_1_R )" , "( a_1_R * a_1_R )" , "( a_1_R * a_2_R )" ,
                     "( a_0_R * a_2_R )" , "( a_1_R * a_2_R )" , "( a_2_R * a_2_R )" ;
           } else {
             bbStr = "( a_0_R * a_0_R )" , "( a_0_R * a_1_R )" ,
                     "( a_0_R * a_1_R )" , "( a_1_R * a_1_R )" ;
           }
           dnudBTens->SetTensor(bbStr,dimDMat,dimDMat,symbolsTens);

           CoefXprTensScalOp dnudBOp = CoefXprTensScalOp( mp_, dnudBTens, scalFnc, CoefXpr::OP_MULT );
           PtrCoefFct dnudBFnc = CoefFunction::Generate( mp_, Global::REAL, dnudBOp );
           return(dnudBFnc);
         }

       } else if( nonlinAnisoParams_.find(MAG_PERMEABILITY_SCALAR) != nonlinAnisoParams_.end() ) {
         // ---------------------------
         // ANISOTROPIC VERSION
         // ---------------------------
         StdVector<MatDescriptorNl> & matNl = nonlinAnisoParams_[MAG_PERMEABILITY_SCALAR];

         UInt numCurves = matNl.GetSize();
         StdVector<Double> angles(numCurves);
         StdVector<Double> zScalings(numCurves);
         StdVector<shared_ptr<CoefFunction> > approx(numCurves);
         // Loop over all entries
         for( UInt i = 0; i < matNl.GetSize(); ++i ) {
           MatDescriptorNl & actNl = matNl[i];
           angles[i] = actNl.angle;
           zScalings[i] = actNl.zScaling;
           
           // check for scaling factor == 0
           if (zScalings[i] == 0.0) {
             WARN("z-scaling factor of nonlinear curve for angle " << angles[i] << " is set to zero!\n"
               << " consideration of nonlinear curves will be disabled in z-direction which leads to 2D anisotropic behavior!")
           }

           //Here we really approximate H(B); see book Kaltenbacher, 2nd, 125ff
           if( actNl.approxType == SMOOTH_SPLINES ) {
             // Check, if smooth spline approximation was already created
             // and initialized
             if( actNl.approxData == NULL) {
               SmoothSpline * sp = new SmoothSpline( actNl.fileName, MAG_PERMEABILITY_SCALAR );
               sp->SetAccuracy( actNl.measAccuracy );
               sp->SetMaxY( actNl.maxVal );
               sp->CalcBestParameter();
               sp->CalcApproximation();
               sp->Print();
               actNl.approxData = sp;
             }
             //now we need the object "CoefFunctionApproxDeriv", which returns by
             //calling coef->getScalar( nuPrime, lpm) the derivative of the reluctivity;
             //Please note: In this case, we do not need the tensor (as in the isotropic case),
             //since the method "CoefFunctionApproxDerivAniso" (see below) will do the job;
             //That's why the object "CoefFunctionApproxDeriv" has the method "GetScalar"!
             //
             ApproxData * sp = actNl.approxData;
             shared_ptr<CoefFunctionApproxDeriv> coef( new CoefFunctionApproxDeriv());
             coef->Init( sp, dimDMat, dependency );

             approx[i] = coef;
           }
           else if( actNl.approxType == ANALYTIC ) {
             //Get the analytic expression for nu'(B)
             CoefXprUnaryOp fluxDensAbsOp = CoefXprUnaryOp( mp_, dependency, CoefXpr::OP_NORM );
             PtrCoefFct fluxDensAbs = CoefFunction::Generate( mp_, Global::REAL, fluxDensAbsOp );

             // get function of B
             std::string nuStr = actNl.analyticExprDeriv;
             shared_ptr<CoefFunctionCompound<Double> > nuFncDeriv(new CoefFunctionCompound<Double>(mp_));
             std::map<std::string,PtrCoefFct> symbolsNu;
             symbolsNu["B"] = fluxDensAbs;

             nuStr.insert(0,"( ");
             nuStr.append(" )");
             nuFncDeriv->SetScalar(nuStr,symbolsNu);

             //store in array
             approx[i] = nuFncDeriv;
           }
         }

         // -------------------------
         // Insertion sort algorithm: we sort the BH-curves starting at smallest
         //                           specified angle
         // ------------------------
         Double compAngle;
         Double compZScaling;
         shared_ptr<CoefFunction> compApprox;
         UInt j;
         for( UInt i = 1; i < numCurves; i++ ) {
           compAngle = angles[i];
           compZScaling = zScalings[i];
           compApprox = approx[i];
           j = i;
           while( ( j > 0 ) && ( angles[j - 1] > compAngle ) ) {
             angles[j] = angles[j - 1];
             zScalings[j] = zScalings[j - 1];
             approx[j] = approx[j - 1];
             j = j - 1;
           }
           angles[j] = compAngle;
           zScalings[j] = compZScaling;
           approx[j] = compApprox;
         }
         
         if( matType == MAG_RELUCTIVITY_TENSOR ) {
           // get linear starting value
           Double startVal = 0.0;
           this->GetScalar( startVal, matType, Global::REAL );
           shared_ptr<CoefFunctionApproxAniso> coef( new CoefFunctionApproxAniso());
           coef->Init( startVal, approx, angles, zScalings, dependency );
           ret = coef;
         }
         else if (matType == MAG_RELUCTIVITY_DERIV ) {
           //used for the bilinear form of the Newton method
           shared_ptr<CoefFunctionApproxDerivAniso> coef( new CoefFunctionApproxDerivAniso());
           coef->Init( approx, angles, zScalings, dimDMat, dependency, baseCoefAniso_ );
           ret = coef;
         }

       } else {
         EXCEPTION( "No nonlinear definition found for material type '"
             << MaterialTypeEnum.ToString(matType) << "'");
       }

       return ret;
     }

}
