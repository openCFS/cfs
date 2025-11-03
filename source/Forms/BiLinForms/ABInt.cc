// =====================================================================================
// 
// 
//    Description:  Implementation file for ABIntegrators
//                  TAKE care:
//                  This file is just for inclusion in the header file!
// 
//        Version:  1.0
//        Created:  11/03/2011 08:04:43 PM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Andreas Hueppe (AHU), andreas.hueppe@uni-klu.ac.at
//        Company:  Universitaet Klagenfurt
// 
// =====================================================================================
#include "ABInt.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Domain/ElemMapping/SurfElem.hh"

namespace CoupledField{

template< class COEF_DATA_TYPE, class B_DATA_TYPE>
  ABInt<COEF_DATA_TYPE, B_DATA_TYPE>::
    ABInt(BaseBOperator * aOp, BaseBOperator * bOp, 
          PtrCoefFct scalCoef, MAT_DATA_TYPE factor, bool coordUpdate ): 
      BBInt<COEF_DATA_TYPE, B_DATA_TYPE>(bOp, scalCoef, factor, coordUpdate ) {
  this->type_ = BiLinearForm::AB_INT;
  this->name_ = "ABInt"; // no idea, why this doesn't compile: type.ToString(type_); do
  this->aOperator_ = aOp;
  this->solDependent_ = false;

  // Note: In general the AB-Integrator is not symmetric, as is should 
  // get only used, if the A- and B differential operators are distinct.
  this->isSymmetric_ = false;
}


template< class COEF_DATA_TYPE, class B_DATA_TYPE>
  void ABInt<COEF_DATA_TYPE, B_DATA_TYPE>::
    CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                       EntityIterator& ent1,
                       EntityIterator& ent2 ) {
  // Extract physical element
  const Elem* ptElem = ent1.GetElem();

  MAT_DATA_TYPE fac(0.0);

  // Obtain FE element from feSpace and integration scheme
  IntegOrder order1, order2;
  IntScheme::IntegMethod method1, method2;
  BaseFE* ptFeA = this->ptFeSpace1_->GetFe( ent1, method1, order1 );
  BaseFE* ptFeB = this->ptFeSpace2_->GetFe( ent2, method2, order2 );

  const UInt nrFncsA = ptFeA->GetNumFncs();
  const UInt nrFncsB = ptFeB->GetNumFncs();

  // Get shape map from grid
  shared_ptr<ElemShapeMap> esm = ent1.GetGrid()->GetElemShapeMap( ptElem, this->coordUpdate_ );

  // Get integration points
  StdVector<LocPoint> intPoints;
  StdVector<Double> weights;
  this->intScheme_->GetIntPoints( Elem::GetShapeType(ptElem->type), 
      method1, order1, method2, order2, intPoints, weights );

  elemMat.Resize( nrFncsA * aOperator_->GetDimDof(), nrFncsB * this->bOperator_->GetDimDof() );
  elemMat.Init();

  // Loop over all integration points
  LocPointMapped lpm;
  const UInt numIntPts = intPoints.GetSize();
  for( UInt iIntPts = 0; iIntPts < numIntPts; ++iIntPts  ) {

    // Calculate for each integration point the LocPointMapped
    lpm.Set( intPoints[iIntPts], esm, weights[iIntPts] );

    // Calculate A-matrix (first differential operator)
    this->aOperator_->CalcOpMat( aMat_, lpm, ptFeA );

    // Calculate B-matrix (second differential operator)
    this->bOperator_->CalcOpMat( this->bMat_, lpm, ptFeB );

    // Calculate scalar factor
    this->coefScalar_->GetScalar(fac, lpm);
    fac *= MAT_DATA_TYPE(lpm.jacDet * weights[iIntPts]); 

#ifdef NDEBUG
    aMat_.Mult_Blas(this->bMat_,elemMat,true,false,this->factor_*fac,1.0);
#else
    elemMat += Transpose(aMat_) * this->bMat_ * this->factor_*fac;
#endif
    //std::cout << "A mat " << this->aMat_.ToString() << std::endl;
    //std::cout << "B mat " << this->bMat_.ToString() << std::endl;
    //std::cout << "Elem mat " << elemMat.ToString() << std::endl;
  }
}


template< class COEF_DATA_TYPE, class B_DATA_TYPE>
  ABIntLem<COEF_DATA_TYPE, B_DATA_TYPE>::
    ABIntLem(BaseBOperator * aOp, BaseBOperator * bOp, 
          PtrCoefFct scalCoef, MAT_DATA_TYPE factor,
          bool coordUpdate ): 
          ABInt<COEF_DATA_TYPE,B_DATA_TYPE>(aOp, bOp, scalCoef, factor, coordUpdate) {
  this->type_ = BiLinearForm::AB_INT;
  this->name_ = "ABIntLem"; 
  this->aOperator_ = aOp;
  this->solDependent_ = false;

  // Note: In general the AB-Integrator is not symmetric, as is should 
  // get only used, if the A- and B differential operators are distinct.
  this->isSymmetric_ = false;
}



template< class COEF_DATA_TYPE, class B_DATA_TYPE>
  void ABIntLem<COEF_DATA_TYPE, B_DATA_TYPE>::
    CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                       EntityIterator& ent1,
                       EntityIterator& ent2 ) {
  // Extract physical element
  const Elem* ptElem1 = ent1.GetElem();
  const Elem* ptElem2 = ent2.GetElem();


  const SurfElem *sElem1 = nullptr;
  const SurfElem *sElem2 = nullptr;
  const Elem *ptVolElem1 = nullptr;
  const Elem *ptVolElem2 = nullptr;


  BaseFE* ptFeA;
  BaseFE* ptFeB;

  MAT_DATA_TYPE fac(0.0);

  // Obtain FE element from feSpace and integration scheme
  IntegOrder order1, order2;
  IntScheme::IntegMethod method1, method2;
  
  // Check if we have to evaluate the volume or surface element in the case of special conditions for the surface
  if ( this->GetUseVolEqnA() )
  {
    // get surface element
    sElem1 = ent1.GetSurfElem();
    // get volume element
    ptVolElem1 = sElem1->ptVolElems[0];

    ptFeA = this->ptFeSpace1_->GetFe(ptVolElem1->elemNum);
    
    this->ptFeSpace1_->GetIntegration(ptFeA, ptVolElem1->regionId, method1, order1);
  }
  else
  {
    ptFeA = this->ptFeSpace1_->GetFe(ent1, method1, order1);
  }

  if ( this->GetUseVolEqnB() )
  {
    // get surface element
    sElem2 = ent2.GetSurfElem();
    // get volume element
    ptVolElem2 = sElem2->ptVolElems[0];

    ptFeB = this->ptFeSpace2_->GetFe(ptVolElem2->elemNum);

    this->ptFeSpace2_->GetIntegration(ptFeB, ptVolElem2->regionId, method2, order2);
  }
  else
  {
    ptFeB = this->ptFeSpace2_->GetFe(ent2, method2, order2);
  }



  const UInt nrFncsA = ptFeA->GetNumFncs();
  const UInt nrFncsB = ptFeB->GetNumFncs();
  

  StdVector<LocPoint> intPoints;
  StdVector<Double> weights;
  shared_ptr<ElemShapeMap> esm;

  // for LEM coupling we have a special operator that passes us the element matrix from
  // an integrator from the FEM side
  // we search for that one and use it for the evaluation of said integrator
  bool aOpIsAllocOp = false;
  if (this->aOperator_->GetDimElem()==0){
    // if aOp is the alloc operator, then we use the integration points from the B side
    aOpIsAllocOp = true;
    // Get shape map from grid
    esm = ent2.GetGrid()->GetElemShapeMap( ptElem2, this->coordUpdate_ );
    // Get integration points
    this->intScheme_->GetIntPoints( Elem::GetShapeType(ptElem2->type), 
        method2, order2, intPoints, weights );
  } else if (this->bOperator_->GetDimElem()==0){
    // if bOp is the alloc operator, then we use the integration points from the A side
    aOpIsAllocOp = false;
    // Get shape map from grid
    esm = ent1.GetGrid()->GetElemShapeMap( ptElem1, this->coordUpdate_ );
    // Get integration points
    this->intScheme_->GetIntPoints( Elem::GetShapeType(ptElem1->type), 
        method1, order1, intPoints, weights );
  } else {
    EXCEPTION("ABIntLem shall only be used with a FEM allocation operator");
  }
  
  
  elemMat.Resize( nrFncsA * this->aOperator_->GetDimDof(), nrFncsB * this->bOperator_->GetDimDof() );
  elemMat.Init();

  //Matrix<MAT_DATA_TYPE> elemMat2;
  //elemMat2.Resize( nrFncsA * this->aOperator_->GetDimDof(), nrFncsB * this->bOperator_->GetDimDof() );
  //elemMat2.Init();

  // Loop over all integration points
  LocPointMapped lpm;
  const UInt numIntPts = intPoints.GetSize();
  for( UInt iIntPts = 0; iIntPts < numIntPts; ++iIntPts  ) {

    // Calculate for each integration point the LocPointMapped
    lpm.Set( intPoints[iIntPts], esm, weights[iIntPts] );
    //lpm.SetWithSurface( intPoints[iIntPts], esm, this->volRegions_, weights[iIntPts] );

    // Decide which FE function to pass (we always need the one of the FEM part for the correct number of equations)
    if( aOpIsAllocOp ){
      // if we are dealing with network coupling, we might evaluate a surface within an
      // volume element integrator
      // therefore, we have the ability to manually override the integrator and tell it
      // to be a surface integrator
      
      // set bOp to be a surface operator
      //this->bOperator_->OverrideIsSurfOperator(this->overrideSurfaceInt_);

      // Calculate A-matrix (first differential operator)
      this->aOperator_->CalcOpMat( this->aMat_, lpm, ptFeA );

      // Calculate B-matrix (second differential operator)
      this->bOperator_->CalcOpMat( this->bMat_, lpm, ptFeB );

      // reset it
      //this->bOperator_->OverrideIsSurfOperator(false);

    } else {
      // if we are dealing with network coupling, we might evaluate a surface within an
      // volume element integrator
      // therefore, we have the ability to manually override the integrator and tell it
      // to be a surface integrator
      
      //this->aOperator_->OverrideIsSurfOperator(this->overrideSurfaceInt_);

      // Calculate A-matrix (first differential operator)
      this->aOperator_->CalcOpMat( this->aMat_, lpm, ptFeA );

      // Calculate B-matrix (second differential operator)
      this->bOperator_->CalcOpMat( this->bMat_, lpm, ptFeB );

      // reset it
      //this->aOperator_->OverrideIsSurfOperator(false);
    }
    
    //std::cout << "A mat " << this->aMat_.ToString() << std::endl;
    //std::cout << "B mat " << this->bMat_.ToString() << std::endl;
    //std::cout << "weight " << weights[iIntPts] << std::endl;
    // Calculate scalar factor
    this->coefScalar_->GetScalar(fac, lpm);
    fac *= MAT_DATA_TYPE(lpm.jacDet * weights[iIntPts]); 
    //std::cout << "fac " << fac << std::endl;
    //std::cout << "factor_ " << this->factor_ << std::endl;

#ifdef NDEBUG
    aMat_.Mult_Blas(this->bMat_,elemMat,true,false,this->factor_*fac,1.0);
#else
    elemMat += Transpose(this->aMat_) * this->bMat_ * this->factor_*fac;
#endif
    //std::cout << "Elem mat " << elemMat.ToString() << std::endl;

    //elemMat2 = Transpose(this->aMat_) * this->bMat_ * this->factor_*fac;
    //std::cout << "Elem mat2 " << elemMat2.ToString() << std::endl;
  }
}

template< class COEF_DATA_TYPE, class B_DATA_TYPE>
  SurfaceABInt<COEF_DATA_TYPE, B_DATA_TYPE>::
    SurfaceABInt( BaseBOperator * aOp, BaseBOperator * bOp,
                PtrCoefFct scalCoef, MAT_DATA_TYPE factor,
                const std::set<RegionIdType>& volRegions,
                bool coordUpdate ): 
      ABInt<COEF_DATA_TYPE,B_DATA_TYPE>(aOp, bOp, scalCoef, factor, coordUpdate) {
  this->name_ = "SurfaceABInt";
  volRegions_ = volRegions;
  this->isSymmetric_ = false;
}

template< class COEF_DATA_TYPE, class B_DATA_TYPE>
  SurfaceABInt<COEF_DATA_TYPE, B_DATA_TYPE>::
    SurfaceABInt( BaseBOperator * aOp, BaseBOperator * bOp,
                  const std::map< RegionIdType, PtrCoefFct >& regionCoefs, MAT_DATA_TYPE factor,
                  const std::set<RegionIdType>& volRegions, bool coordUpdate): 
      ABInt<COEF_DATA_TYPE,B_DATA_TYPE>(aOp, bOp,regionCoefs.begin()->second, factor, coordUpdate )
{
  this->name_ = "SurfaceABInt";
  volRegions_ = volRegions;
  this->isSymmetric_ = false;
  regionCoefs_ = regionCoefs;
}


template< class COEF_DATA_TYPE, class B_DATA_TYPE>
  void SurfaceABInt<COEF_DATA_TYPE, B_DATA_TYPE>::
    CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                       EntityIterator& ent1,
                       EntityIterator& ent2 ) {
  // Extract physical element
  const Elem* ptElem1 = ent1.GetElem();
  const Elem* ptElem2 = ent2.GetElem();

  MAT_DATA_TYPE fac(0.0);

  // Obtain FE element from feSpace and integration scheme

  const SurfElem *sElem1 = nullptr;
  const SurfElem *sElem2 = nullptr;
  const Elem *ptVolElem1 = nullptr;
  const Elem *ptVolElem2 = nullptr;

  BaseFE *ptFeA = nullptr;
  BaseFE *ptFeB = nullptr;
  IntegOrder order1, order2;
  IntScheme::IntegMethod method1, method2;

  // Check if we have to evaluate the volume or surface element in the case of special conditions for the surface
  if ( this->GetUseVolEqnA() )
  {
    // get surface element
    sElem1 = ent1.GetSurfElem();
    // get volume element
    ptVolElem1 = sElem1->ptVolElems[0];

    ptFeA = this->ptFeSpace1_->GetFe(ptVolElem1->elemNum);
    
    this->ptFeSpace1_->GetIntegration(ptFeA, ptVolElem1->regionId, method1, order1);
  }
  else
  {
    ptFeA = this->ptFeSpace1_->GetFe(ent1, method1, order1);
  }

  if ( this->GetUseVolEqnB() )
  {
    // get surface element
    sElem2 = ent2.GetSurfElem();
    // get volume element
    ptVolElem2 = sElem2->ptVolElems[0];

    ptFeB = this->ptFeSpace2_->GetFe(ptVolElem2->elemNum);

    this->ptFeSpace2_->GetIntegration(ptFeB, ptVolElem2->regionId, method2, order2);
  }
  else
  {
    ptFeB = this->ptFeSpace2_->GetFe(ent2, method2, order2);
  }


  const UInt nrFncsA = ptFeA->GetNumFncs();
  const UInt nrFncsB = ptFeB->GetNumFncs();

  // Get shape map from grid
  shared_ptr<ElemShapeMap> esm1 = ent1.GetGrid()->GetElemShapeMap( ptElem1, this->coordUpdate_ );
  shared_ptr<ElemShapeMap> esm2 = ent1.GetGrid()->GetElemShapeMap( ptElem2, this->coordUpdate_ );

  // Get integration points
  StdVector<LocPoint> intPoints;
  StdVector<Double> weights;
  this->intScheme_->GetIntPoints( Elem::GetShapeType(ptElem1->type), 
      method1, order1, method2, order2, intPoints, weights );

  elemMat.Resize( nrFncsA * this->aOperator_->GetDimDof(), nrFncsB * this->bOperator_->GetDimDof() );
  elemMat.Init();

  // Loop over all integration points
  LocPointMapped lpm1,lpm2;
  const UInt numIntPts = intPoints.GetSize();
  for( UInt iIntPts = 0; iIntPts < numIntPts; ++iIntPts ) {

    // Calculate for each integration point the LocPointMapped
    lpm1.SetWithSurface( intPoints[iIntPts], esm1, this->volRegions_, weights[iIntPts] );
    lpm2.SetWithSurface( intPoints[iIntPts], esm2, this->volRegions_, weights[iIntPts] );

    // Calculate A-matrix (first differential operator)
    this->aOperator_->CalcOpMat( this->aMat_, lpm1, ptFeA );

    // Calculate B-matrix (second differential operator)
    this->bOperator_->CalcOpMat( this->bMat_, lpm2, ptFeB );

    // Calculate scalar factor
    if(!regionCoefs_.empty()) {
      regionCoefs_[lpm1.lpmVol->ptEl->regionId]->GetScalar(fac, lpm1);
    } else {
      this->coefScalar_->GetScalar(fac, lpm1);
    }

    fac *= MAT_DATA_TYPE(lpm1.jacDet * weights[iIntPts]);

#ifdef NDEBUG
    this->aMat_.Mult_Blas(this->bMat_, elemMat, true, false,
        this->factor_*fac, 1.0);
#else
    elemMat += Transpose(this->aMat_) * this->bMat_ * this->factor_*fac;
#endif

  }
}


DEFINE_LOG(nitscheInt, "Forms.NitscheInt")
template< class COEF_DATA_TYPE, class B_DATA_TYPE>
  SurfaceNitscheABInt<COEF_DATA_TYPE, B_DATA_TYPE>::
    SurfaceNitscheABInt(BaseBOperator * aOp, BaseBOperator * bOp,
                        PtrCoefFct scalCoef, MAT_DATA_TYPE factor,
                        BiLinearForm::CouplingDirection cplDir,
                        bool coordUpdate, bool isSym, bool isPenalty ): 
      ABInt<COEF_DATA_TYPE,B_DATA_TYPE>(aOp, bOp, scalCoef, factor, coordUpdate) {
  
  this->name_ = "SurfaceNitscheABInt";
  this->cplDirection_ = cplDir;
  this->isSymmetric_ = isSym;
  this->isPenalty_ = isPenalty;
}


template< class COEF_DATA_TYPE, class B_DATA_TYPE>
  void SurfaceNitscheABInt<COEF_DATA_TYPE, B_DATA_TYPE>::
    CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                       EntityIterator& ent1,
                       EntityIterator& ent2 ) {
  /*do now the following>
    1. Perform some security checks
    2. fill a volregion set according to the coupling direction
       in order to determine a consistent normal direction
    3. obtain pointers to associated volume elements
    4. loop over int points of surface element and do a standard integration
    5. FE spaces are set according to which pdes are coupled. E.g., in the SinglePDE, 
       ptFeSpace1_ and ptFeSpace2_ are set to points to the same type of feSpace (of the unknown).
       For inter-PDE coupling, it is a bit different. 
       Hamideh: ptFeSpace1_ and ptFeSpace2_ points to Primary and Secondary unknowns.
                The true feSpace for test function and  unknow will be taken from coupling direction
                For example in the case of LinFlowMechNitsche coupling:
                Primary --> Fespace1_(mechanics displacement)
                Secondary --> FeSpace2_(LinFlow Velocity)*/

  // Extract MortarNcSurfElem element
  const NcSurfElem* ptNcElem = ent1.GetNcSurfElem();
  const MortarNcSurfElem * ptMortarElem = dynamic_cast<const MortarNcSurfElem*>(ptNcElem);
  
  // Primary-side and secondary-side surface elements
  const SurfElem* ptSurfPrimary = ptMortarElem->ptPrimary;
  const SurfElem* ptSurfSecondary = ptMortarElem->ptSecondary;
  
  // Obtain FE element from feSpace and integration scheme
  BaseFE* ptFePrimary = this->ptFeSpace1_->GetFe(ptSurfPrimary->elemNum);
  BaseFE* ptFeSecondary =  this->ptFeSpace2_->GetFe(ptSurfSecondary->elemNum);

  // get integration context of primary and secondary element (method and order)
  // Integration order must be the same for all Nitsche integrators, regardless of coupling direction.
  IntegOrder orderPrimary, orderSecondary;
  IntScheme::IntegMethod methodPrimary, methodSecondary;
  this->ptFeSpace1_->GetIntegration(ptFePrimary, ptSurfPrimary->ptVolElems[0]->regionId, methodPrimary, orderPrimary);
  this->ptFeSpace2_->GetIntegration(ptFeSecondary, ptSurfSecondary->ptVolElems[0]->regionId, methodSecondary, orderSecondary);

  // assign vol elems and fe functions according to which term is constructed (PrimPrim, PrimSec, SecPrim, SecSec)
  // this is more or less a switch case statement which terms to consider and where to acquire the equation numbers
  bool usePrimary1, usePrimary2;
  Elem* ptVolElem1 = nullptr;
  Elem* ptVolElem2 = nullptr;
  BaseFE* ptFe1 = nullptr;
  BaseFE* ptFe2 = nullptr;
  this->SetCouplingDirection(ptSurfPrimary, ptSurfSecondary, ptVolElem1, ptVolElem2, ptFe1, ptFe2, usePrimary1, usePrimary2);

  // when assembling the penalty terms, we need the penaltyFactor=+-beta/h. Otherwise, it is normally +-1.
  // the +-beta or +-1 needs to be passed correctly via the constructor and is stored in the factor_ member.
  // this is done, e.g., in the SinglePDE class, where the Nitsche integrators are created.
  MAT_DATA_TYPE penaltyFactor = this->factor_; // -1.0
  if (this->isPenalty_)
    ComputePenalty(ptSurfPrimary, ptSurfSecondary, ent1, penaltyFactor); // beta/h

  // initialize element matrix
  const UInt nrFncs1 = ptFe1->GetNumFncs();
  const UInt nrFncs2 = ptFe2->GetNumFncs();
  elemMat.Resize( nrFncs1*this->aOperator_->GetDimDof(), nrFncs2*this->bOperator_->GetDimDof());
  elemMat.Init();

  // Get shape map of the mortar element from grid
  shared_ptr<ElemShapeMap> esmNc = ent1.GetGrid()->GetElemShapeMap(ptMortarElem, this->coordUpdate_);

  // Get local integration points and weights
  StdVector<LocPoint> intPoints;
  StdVector<Double> weights;
  this->intScheme_->GetIntPoints(Elem::GetShapeType(ptMortarElem->type), methodPrimary, orderPrimary, methodSecondary, orderSecondary, intPoints, weights);

#ifndef NDEBUG
  assert(ent1.GetSurfElem() == ent2.GetSurfElem());
  assert(ptMortarElem != nullptr);
  assert(ptSurfPrimary != nullptr);
  assert(ptSurfSecondary != nullptr);
  assert(ptFePrimary != nullptr);
  assert(ptFeSecondary != nullptr);
  LOG_DBG2(nitscheInt) << "Region primary: " << ptMortarElem->ptPrimary->ptVolElems[0]->regionId;
  LOG_DBG2(nitscheInt) << "Region secondary: " << ptMortarElem->ptSecondary->ptVolElems[0]->regionId;
  LOG_DBG2(nitscheInt) << "elemNum primary: " << ptMortarElem->ptPrimary->ptVolElems[0]->elemNum;
  LOG_DBG2(nitscheInt) << "elemNum secondary: " << ptMortarElem->ptSecondary->ptVolElems[0]->elemNum;
  LOG_DBG2(nitscheInt) << "NcSurfElem #" << ptMortarElem->elemNum
                       << ", type " << Elem::feType.ToString(ptMortarElem->type)
                       << "\n\tprimary parent #" << ptSurfPrimary->elemNum
                       << "\n\tsecondary parent #" << ptSurfSecondary->elemNum;
  LOG_DBG2(nitscheInt) << "\tconnectivity: " << ptMortarElem->connect.ToString();
#endif

  MAT_DATA_TYPE fac(0.0);
  Vector<Double> globIntPoint;
  LocPointMapped lpm1, lpm2;
  // Loop over all integration points
  for( UInt iIntPts = 0; iIntPts < intPoints.GetSize(); ++iIntPts ) {
    // Calculate for each integration point the LocPointMapped
    lpm1.SetWithNitscheSurface( intPoints[iIntPts], esmNc, usePrimary1, weights[iIntPts] );
    lpm2.SetWithNitscheSurface( intPoints[iIntPts], esmNc, usePrimary2, weights[iIntPts] );
    // Calculate A-matrix (first differential operator)
    this->aOperator_->CalcOpMat( this->aMat_, lpm1, ptFe1);
    // Calculate B-matrix (second differential operator)
    this->bOperator_->CalcOpMat( this->bMat_, lpm2, ptFe2);

    // Calculate scalar factor. Should not matter if we use lpm1 or lpm2 here.
    this->coefScalar_->GetScalar(fac, lpm1);

    fac *= MAT_DATA_TYPE(lpm1.jacDet * weights[iIntPts]);
#ifdef NDEBUG
    this->aMat_.Mult_Blas(this->bMat_, elemMat, true, false, penaltyFactor*fac, 1.0);
#else
    elemMat += Transpose(this->aMat_) * this->bMat_ * penaltyFactor *fac;
    LOG_DBG3(nitscheInt) << "Element matrix: \n" << elemMat;
    LOG_DBG3(nitscheInt) << "penaltyFactor: "<< penaltyFactor;
    LOG_DBG3(nitscheInt) << "fac: " << fac << std::endl << std::endl;
#endif
  }
}

template< class COEF_DATA_TYPE, class B_DATA_TYPE>
  void SurfaceNitscheABInt<COEF_DATA_TYPE, B_DATA_TYPE>::
    ComputePenalty(const SurfElem* ptSurfPrimary, const SurfElem* ptSurfSecondary, 
                   EntityIterator& ent1, MAT_DATA_TYPE& penaltyFactor) { 
  // get shape maps of surface elements
  shared_ptr<ElemShapeMap> esm1T = ent1.GetGrid()->GetElemShapeMap(ptSurfPrimary, this->coordUpdate_);
  shared_ptr<ElemShapeMap> esm2T = ent1.GetGrid()->GetElemShapeMap(ptSurfSecondary, this->coordUpdate_);
  //obtain pointer to basis functions of the surface elements
  BaseFE* SFe1 = this->ptFeSpace1_->GetFe(ptSurfPrimary->elemNum);
  BaseFE* SFe2 = this->ptFeSpace2_->GetFe(ptSurfSecondary->elemNum);

  // get the minimum length of the primary element
  Double min,max;
  esm1T->GetMaxMinEdgeLength(min,max);
  // scale by iso order...
  UInt isoOrder1 = SFe1->GetIsoOrder();
  if(isoOrder1 == 0){
    StdVector<UInt> aIsoOrder;
    SFe1->GetAnisoOrder(aIsoOrder);
    // set maximum occurring order
    for(UInt iOrder = 0; iOrder < aIsoOrder.GetSize(); ++iOrder)
      isoOrder1 = (aIsoOrder[iOrder]>isoOrder1)? aIsoOrder[iOrder] : isoOrder1;
    // that's a bit dirty because zero order Nedelec elements are neither zero nor first order...
    // but here we set it to first, otherwise we divide by zero
    if(isoOrder1 == 0 && this->ptFeSpace1_->GetSpaceType() == this->ptFeSpace1_->HCURL)
      isoOrder1 = 1;
  }
  // scale minumum edge length by iso order
  MAT_DATA_TYPE surface1(min/(Double)isoOrder1);

  // same for secondary element...
  esm2T->GetMaxMinEdgeLength(min,max);
  UInt isoOrder2 = SFe2->GetIsoOrder();
  if(isoOrder2 == 0){
    StdVector<UInt> aIsoOrder;
    SFe2->GetAnisoOrder(aIsoOrder);
    // set maximum occurring order
    for(UInt iOrder = 0; iOrder < aIsoOrder.GetSize(); ++iOrder)
      isoOrder2 = (aIsoOrder[iOrder]>isoOrder2)? aIsoOrder[iOrder] : isoOrder2;
    // that's a bit dirty because zero order Nedelec elements are neither zero nor first order...
    // but here we set it to first, otherwise we divide by zero
    if(isoOrder2 == 0 && this->ptFeSpace1_->GetSpaceType() == this->ptFeSpace1_->HCURL)
      isoOrder2 = 1;
  }
  // scale minumum edge length by iso order
  MAT_DATA_TYPE surface2(min/(Double)isoOrder2);
  // scale Nitsche factor by mean of shortest edge lengths
  MAT_DATA_TYPE tmp(2.0);
  // compute 1/h
  penaltyFactor = this->factor_ * (tmp/(surface1+surface2));
}

template< class COEF_DATA_TYPE, class B_DATA_TYPE>
  void SurfaceNitscheABInt<COEF_DATA_TYPE, B_DATA_TYPE>::
    SetCouplingDirection(const SurfElem* ptSurfPrimary, const SurfElem* ptSurfSecondary, Elem*& ptVolElem1, Elem*& ptVolElem2, 
                         BaseFE*& ptFe1, BaseFE*& ptFe2, bool& usePrimary1, bool& usePrimary2) {
  switch(this->cplDirection_) {
    case BiLinearForm::PRIM_PRIM:
      usePrimary1 = true;
      usePrimary2 = true;
      break;
    case BiLinearForm::SEC_SEC:
      usePrimary1 = false;
      usePrimary2 = false;
      break;
    case BiLinearForm::PRIM_SEC:
      usePrimary1 = true;
      usePrimary2 = false;
      break;
    case BiLinearForm::SEC_PRIM:
      usePrimary1 = false;
      usePrimary2 = true;
      break;
    default:
      EXCEPTION("Undefined coupling direction");
      break;
  }
  // set the volume elements according to the current coupling term
  ptVolElem1 = (usePrimary1) ? ptSurfPrimary->ptVolElems[0] : ptSurfSecondary->ptVolElems[0]; // take pointer volume element of test function
  ptVolElem2 = (usePrimary2) ? ptSurfPrimary->ptVolElems[0] : ptSurfSecondary->ptVolElems[0]; // take pointer volume element of unknown
  // Obtain FE element from FeSpace for operators according to the current coupling term (different FeSpaces are only used for inter-PDE coupling)
  ptFe1 = (usePrimary1) ? this->ptFeSpace1_->GetFe(ptVolElem1->elemNum) : this->ptFeSpace2_->GetFe(ptVolElem1->elemNum);
  ptFe2 = (usePrimary2) ? this->ptFeSpace1_->GetFe(ptVolElem2->elemNum) : this->ptFeSpace2_->GetFe(ptVolElem2->elemNum);
}


DEFINE_LOG(mortarInt, "Forms.MortarInt")
template< class COEF_DATA_TYPE, class B_DATA_TYPE>
  SurfaceMortarABInt<COEF_DATA_TYPE, B_DATA_TYPE>::
    SurfaceMortarABInt(BaseBOperator * aOp, BaseBOperator * bOp,
                       PtrCoefFct scalCoef, MAT_DATA_TYPE factor,
                       RegionIdType primaryVolRegion, RegionIdType secondaryVolRegion,
                       bool coplanar, bool coordUpdate,
                       BiLinearForm::CouplingDirection cplDirection): 
      ABInt<COEF_DATA_TYPE, B_DATA_TYPE>( aOp, bOp, scalCoef, factor, coordUpdate),
      isCoplanar_(coplanar), cplDirection_(cplDirection) {
  this->name_ = "SurfaceMortarABInt";
  this->isSymmetric_ = true;
  
  this->primaryVolRegion_ = primaryVolRegion;
  this->secondaryVolRegion_ = secondaryVolRegion;
  
  // create set of volume regions
  volRegions_.insert(primaryVolRegion_);
  volRegions_.insert(secondaryVolRegion_);
}


template< class COEF_DATA_TYPE, class B_DATA_TYPE>
  void SurfaceMortarABInt<COEF_DATA_TYPE, B_DATA_TYPE>::
    CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                       EntityIterator& ent1,
                       EntityIterator& ent2 ) {

  // Extract MortarNcSurfElem element
  const NcSurfElem* ptNcElem = ent1.GetNcSurfElem();
  const MortarNcSurfElem* ptMortarElem = dynamic_cast<const MortarNcSurfElem*>(ptNcElem);
  
  // Primary-side and secondary-side surface elements
  const SurfElem* ptSurfPrimary = ptMortarElem->ptPrimary;
  const SurfElem* ptSurfSecondary = ptMortarElem->ptSecondary;
  
  // Obtain FE element from feSpace and integration scheme
  BaseFE* ptFePrimary = this->ptFeSpaceField_->GetFe(ptSurfPrimary->elemNum);
  BaseFE* ptFeSecondary = this->ptFeSpaceLM_->GetFe(ptSurfSecondary->elemNum);
  
  // get integration context of primary and secondary element (method and order)
  IntegOrder orderPrimary, orderSecondary;
  IntScheme::IntegMethod methodPrimary, methodSecondary;
  this->ptFeSpaceField_->GetIntegration(ptFePrimary, ptSurfPrimary->ptVolElems[0]->regionId, methodPrimary, orderPrimary);
  this->ptFeSpaceField_->GetIntegration(ptFeSecondary, ptSurfSecondary->ptVolElems[0]->regionId, methodSecondary, orderSecondary);

  // initialize element matrix
  const UInt nrFncsPrimary = ptFePrimary->GetNumFncs();
  const UInt nrFncsSecondary = ptFeSecondary->GetNumFncs();
  elemMat.Resize( nrFncsPrimary * this->aOperator_->GetDimDof(), nrFncsSecondary * this->bOperator_->GetDimDof() );
  elemMat.Init();

  // Get shape maps from grid
  shared_ptr<ElemShapeMap> esmNc = ent1.GetGrid()->GetElemShapeMap( ptMortarElem, this->coordUpdate_ );
  shared_ptr<ElemShapeMap> esmPrimary = ent1.GetGrid()->GetElemShapeMap( ptSurfPrimary, this->coordUpdate_ );
  shared_ptr<ElemShapeMap> esmSecondary = ent1.GetGrid()->GetElemShapeMap( ptSurfSecondary, this->coordUpdate_ );

  // Get local integration points and weights
  StdVector<LocPoint> intPoints;
  StdVector<Double> weights;
  this->intScheme_->GetIntPoints(Elem::GetShapeType(ptMortarElem->type), methodPrimary, orderPrimary, methodSecondary, orderSecondary, intPoints, weights);

#ifndef NDEBUG
  assert(ent1.GetSurfElem() == ent2.GetSurfElem());
  assert(ptMortarElem != nullptr);
  assert(ptSurfPrimary != nullptr);
  assert(ptSurfSecondary != nullptr);
  assert(ptFePrimary != nullptr);
  assert(ptFeSecondary != nullptr);
  LOG_DBG2(mortarInt) << "Region primary: " << ptMortarElem->ptPrimary->ptVolElems[0]->regionId;
  LOG_DBG2(mortarInt) << "Region secondary: " << ptMortarElem->ptSecondary->ptVolElems[0]->regionId;
  LOG_DBG2(mortarInt) << "elemNum primary: " << ptMortarElem->ptPrimary->ptVolElems[0]->elemNum;
  LOG_DBG2(mortarInt) << "elemNum secondary: " << ptMortarElem->ptSecondary->ptVolElems[0]->elemNum;
  LOG_DBG2(mortarInt) << "NcSurfElem #" << ptMortarElem->elemNum
                       << ", type " << Elem::feType.ToString(ptMortarElem->type)
                       << "\n\tprimary parent #" << ptSurfPrimary->elemNum
                       << "\n\tsecondary parent #" << ptSurfSecondary->elemNum;
  LOG_DBG2(mortarInt) << "\tconnectivity: " << ptMortarElem->connect.ToString();
#endif

  MAT_DATA_TYPE fac(0.0);
  Vector<Double> globIntPoint;
  LocPoint ipPrimary, ipSecondary;
  LocPointMapped lpmNc, lpmPrimary, lpmSecondary;
  // Loop over all integration points
  for( UInt iIntPts = 0; iIntPts < intPoints.GetSize(); ++iIntPts ) {
    // perform back projection of integration points...
    // Calculate global coordinates of integration point
    esmNc->Local2Global(globIntPoint, intPoints[iIntPts]);

    // Calculate local coordinates of integration point in secondary element
    esmSecondary->Global2Local(ipSecondary.coord, globIntPoint);

    // Projection for curved interfaces
    if ( !isCoplanar_ ) {
      /* Instead of projecting the integration point back into the primary
       * element, we use a better approach: Orthogonal projection is a linear
       * operation, so any point in the projected primary element has the same
       * local coordinates as its orthogonal projection in the original primary
       * element. So we perform the Global2Local mapping in the projected
       * element, which has been computed already by the intersection
       * algorithm. Sometimes the intersection is the projected primary element
       * itself, in which case we need no coordinate mapping at all.  
       */ 
      if ( ptMortarElem->projectedPrimary != nullptr) {
        shared_ptr<ElemShapeMap> esmProj = ent1.GetGrid()->GetElemShapeMap(ptMortarElem->projectedPrimary.get(), this->coordUpdate_);
        esmProj->Global2Local(ipPrimary.coord, globIntPoint);
      } else { // projectedPrimary == nullptr means it is the MortarElem itself
        ipPrimary.coord = intPoints[iIntPts].coord;
      }
    } else {
      // Calculate local coordinates of integration point in primary element
      esmPrimary->Global2Local(ipPrimary.coord, globIntPoint); 
    }

#ifndef NDEBUG
    LOG_DBG3(mortarInt) << "Integration point #" << iIntPts+1
      << "\n\tglobal coordinates: " << globIntPoint.ToString()
      << "\n\tlocal coordinates: " << intPoints[iIntPts].coord.ToString()
      << "\n\tlocal coordinates in primary: " << ipPrimary.coord.ToString()
      << "\n\tlocal coordinares in secondary: " << ipSecondary.coord.ToString();
#endif

    // Calculate for each integration point the LocPointMapped
    lpmNc.Set( intPoints[iIntPts], esmNc, weights[iIntPts] );
    lpmPrimary.SetWithSurface( ipPrimary, esmPrimary, this->volRegions_, weights[iIntPts] );
    lpmSecondary.SetWithSurface( ipSecondary, esmSecondary, this->volRegions_, weights[iIntPts] );

    // Calculate A-matrix (first differential operator)
    // Calculate B-matrix (second differential operator)
    if (cplDirection_ == BiLinearForm::PRIM_SEC) {
      this->aOperator_->CalcOpMat( this->aMat_, lpmPrimary, ptFePrimary );
      this->bOperator_->CalcOpMat( this->bMat_, lpmSecondary, ptFeSecondary );
    }
    else if (cplDirection_ == BiLinearForm::SEC_PRIM) {
      this->aOperator_->CalcOpMat( this->aMat_, lpmSecondary, ptFeSecondary );
      this->bOperator_->CalcOpMat( this->bMat_, lpmPrimary, ptFePrimary );
    } else {
      EXCEPTION("Cannot calculate operator matrices for the current coupling direction: " + boost::lexical_cast<std::string>(cplDirection_));
    }

    // Calculate scalar factor
    this->coefScalar_->GetScalar(fac, lpmSecondary);
    fac *= MAT_DATA_TYPE(lpmNc.jacDet * weights[iIntPts]);

#ifdef NDEBUG
    this->aMat_.Mult_Blas(this->bMat_, elemMat, true, false, this->factor_*fac, 1.0);
#else
    assert(esmPrimary->CoordIsInsideElem(ipPrimary.coord, EPS));
    assert(esmSecondary->CoordIsInsideElem(ipSecondary.coord, EPS));
    elemMat += Transpose(this->aMat_) * this->bMat_ * this->factor_*fac;
    LOG_DBG3(mortarInt) << "Element matrix: \n" << elemMat;
#endif
  }
}


template< class COEF_DATA_TYPE, class B_DATA_TYPE>
  void SurfaceMortarABInt<COEF_DATA_TYPE, B_DATA_TYPE>::
    SetFeSpace( shared_ptr<FeSpace> feSpace1, shared_ptr<FeSpace> feSpace2 ) {
  
  shared_ptr<BaseFeFunction> feFunc1 = feSpace1->GetFeFunction().lock(),
                             feFunc2 = feSpace2->GetFeFunction().lock();
  // IsLagrSurfSpace check is included to be able to work with other kinds of Lagrange multipliers
  if (feFunc1->GetResultInfo()->resultType == LAGRANGE_MULT)
  {
    ptFeSpaceLM_ = feFunc1->GetFeSpace();
    ptFeSpaceField_ = feFunc2->GetFeSpace();
  }
  else if (feFunc2->GetResultInfo()->resultType == LAGRANGE_MULT)
  {
    ptFeSpaceLM_ = feFunc2->GetFeSpace();
    ptFeSpaceField_ = feFunc1->GetFeSpace();
  }
  else
  {
    EXCEPTION("FeSpace for Lagrange multiplier is not defined");
  }
  
  this->intScheme_ = ptFeSpaceLM_->GetIntScheme();
  assert(ptFeSpaceLM_ != nullptr);
  assert(ptFeSpaceField_ != nullptr);
}

DEFINE_LOG(mortarIntMA, "Forms.MortarIntMA")
template<class COEF_DATA_TYPE, class B_DATA_TYPE>
  SurfaceMortarABIntMA<COEF_DATA_TYPE, B_DATA_TYPE>::
    SurfaceMortarABIntMA(BaseBOperator * aOp, BaseBOperator * bOp,
                         PtrCoefFct scalCoef, MAT_DATA_TYPE factor,
                         bool coplanar, bool coordUpdate): 
      SurfaceABInt<COEF_DATA_TYPE, B_DATA_TYPE>(aOp, bOp, scalCoef, factor,
                                                std::set<RegionIdType>(), coordUpdate) {
  //  isCoplanar_(coplanar), cplDirection_(cplDirection)
  this->isSymmetric_ = false;
  this->isCoplanar_ = coplanar;
}

//mechanical-acoustic coupling on non-conforming grids!
template<class COEF_DATA_TYPE, class B_DATA_TYPE>
  SurfaceMortarABIntMA<COEF_DATA_TYPE, B_DATA_TYPE>::
    SurfaceMortarABIntMA(BaseBOperator * aOp, BaseBOperator * bOp,
                         const std::map< RegionIdType, PtrCoefFct >& regionCoefs,
                         MAT_DATA_TYPE factor, bool coplanar, bool coordUpdate): 
      SurfaceABInt<COEF_DATA_TYPE, B_DATA_TYPE>(aOp, bOp, regionCoefs, factor,
                                                std::set<RegionIdType>(), coordUpdate) {
  this->name_ = "SurfaceMortarABIntMA";
  this->regionCoefs_ = regionCoefs;
  this->isSymmetric_ = false;
  this->isCoplanar_  = coplanar;
}


template< class COEF_DATA_TYPE, class B_DATA_TYPE>
  void SurfaceMortarABIntMA<COEF_DATA_TYPE, B_DATA_TYPE>::
    CalcElementMatrix(Matrix<MAT_DATA_TYPE>& elemMat,
                      EntityIterator& ent1,
                      EntityIterator& ent2) {

  // Extract MortarNcSurfElem element
  const NcSurfElem* ptNcElem = ent1.GetNcSurfElem();
  const MortarNcSurfElem* ptMortarElem = dynamic_cast<const MortarNcSurfElem*>(ptNcElem);

  // Primary-side and secondary-side surface elements
  const SurfElem* ptSurfPrimary = ptMortarElem->ptPrimary;
  const SurfElem* ptSurfSecondary = ptMortarElem->ptSecondary;

  // Obtain FE element from feSpace and integration scheme
  BaseFE* ptFePrimary = this->ptFeSpacePrimary_->GetFe(ptSurfPrimary->elemNum);
  BaseFE* ptFeSecondary = this->ptFeSpaceSecondary_->GetFe(ptSurfSecondary->elemNum);

  // get integration context of primary and secondary element (method and order)
  IntegOrder orderPrimary, orderSecondary;
  IntScheme::IntegMethod methodPrimary, methodSecondary;
  this->ptFeSpacePrimary_->GetIntegration(ptFePrimary, ptSurfPrimary->ptVolElems[0]->regionId, methodPrimary, orderPrimary);
  this->ptFeSpaceSecondary_->GetIntegration(ptFeSecondary, ptSurfSecondary->ptVolElems[0]->regionId, methodSecondary, orderSecondary);

  // initialize element matrix
  const UInt nrFncsPrimary = ptFePrimary->GetNumFncs();
  const UInt nrFncsSecondary = ptFeSecondary->GetNumFncs();
  Matrix<MAT_DATA_TYPE> result; // store separately so we can transpose later in case of doTranspose_ == true
  result.Resize( nrFncsPrimary * this->ptPrimaryOp_->GetDimDof(), nrFncsSecondary * this->ptSecondaryOp_->GetDimDof() );
  result.Init();

  // Get shape maps from grid
  shared_ptr<ElemShapeMap> esmNc = ent1.GetGrid()->GetElemShapeMap( ptMortarElem, this->coordUpdate_ );
  shared_ptr<ElemShapeMap> esmPrimary = ent1.GetGrid()->GetElemShapeMap( ptSurfPrimary, this->coordUpdate_ );
  shared_ptr<ElemShapeMap> esmSecondary = ent1.GetGrid()->GetElemShapeMap( ptSurfSecondary, this->coordUpdate_ );

  // Get local integration points and weights
  StdVector<LocPoint> intPoints;
  StdVector<Double> weights;
  this->intScheme_->GetIntPoints(Elem::GetShapeType(ptMortarElem->type), methodPrimary, orderPrimary, methodSecondary, orderSecondary, intPoints, weights);

  // assign or update volRegions_ so that they can get assigned to the Set functions later on
  // it should actually not be necessary to set the volRegions_ every time the matrix gets reassembled but that's how it is
  this->volRegions_.insert(ptSurfPrimary->ptVolElems[0]->regionId);
  this->volRegions_.insert(ptSurfSecondary->ptVolElems[0]->regionId);

#ifndef NDEBUG
  assert(ent1.GetSurfElem() == ent2.GetSurfElem());
  assert(ptMortarElem != nullptr);
  assert(ptSurfPrimary != nullptr);
  assert(ptSurfSecondary != nullptr);
  assert(ptFePrimary != nullptr);
  assert(ptFeSecondary != nullptr);
  LOG_DBG2(mortarIntMA) << "Region primary: " << ptMortarElem->ptPrimary->ptVolElems[0]->regionId;
  LOG_DBG2(mortarIntMA) << "Region secondary: " << ptMortarElem->ptSecondary->ptVolElems[0]->regionId;
  LOG_DBG2(mortarIntMA) << "elemNum primary: " << ptMortarElem->ptPrimary->ptVolElems[0]->elemNum;
  LOG_DBG2(mortarIntMA) << "elemNum secondary: " << ptMortarElem->ptSecondary->ptVolElems[0]->elemNum;
  LOG_DBG2(mortarIntMA) << "NcSurfElem #" << ptMortarElem->elemNum
                       << ", type " << Elem::feType.ToString(ptMortarElem->type)
                       << "\n\tprimary parent #" << ptSurfPrimary->elemNum
                       << "\n\tsecondary parent #" << ptSurfSecondary->elemNum;
  LOG_DBG2(mortarIntMA) << "\tconnectivity: " << ptMortarElem->connect.ToString();
#endif

  MAT_DATA_TYPE fac(0.0);
  Vector<Double> globIntPoint;
  LocPoint ipPrimary, ipSecondary;
  LocPointMapped lpmNc, lpmPrimary, lpmSecondary;
  const UInt numIntPts = intPoints.GetSize();
  // Loop over all integration points
  for( UInt iIntPts = 0; iIntPts < numIntPts; ++iIntPts ) {
    // Calculate global coordinates of integration point
    esmNc->Local2Global(globIntPoint, intPoints[iIntPts]);

    // Calculate local coordinates of integration point in secondary element
    esmSecondary->Global2Local(ipSecondary.coord, globIntPoint);

    // Projection for curved interfaces
    if ( !isCoplanar_ ) {
      /* Instead of projecting the integration point back into the primary
       * element, we use a better approach: Orthogonal projection is a linear
       * operation, so any point in the projected primary element has the same
       * local coordinates as its orthogonal projection in the original primary
       * element. So we perform the Global2Local mapping in the projected
       * element, which has been computed already by the intersection
       * algorithm. Sometimes the intersection is the projected primary element
       * itself, in which case we need no coordinate mapping at all.
       */
      if ( ptMortarElem->projectedPrimary != nullptr) {
        shared_ptr<ElemShapeMap> esmProj = ent1.GetGrid()->GetElemShapeMap(ptMortarElem->projectedPrimary.get(), this->coordUpdate_);
        esmProj->Global2Local(ipPrimary.coord, globIntPoint);
      } else { // projectedPrimary == nullptr means it is the MortarElem itself
        ipPrimary.coord = intPoints[iIntPts].coord;
      }
    } else {
      // Calculate local coordinates of integration point in primary element
      esmPrimary->Global2Local(ipPrimary.coord, globIntPoint);
    }

    // Calculate for each integration point the LocPointMapped
    lpmNc.Set( intPoints[iIntPts], esmNc, weights[iIntPts] );
    lpmPrimary.SetWithSurface( ipPrimary, esmPrimary, this->volRegions_, weights[iIntPts] );
    lpmSecondary.SetWithSurface( ipSecondary, esmSecondary, this->volRegions_, weights[iIntPts] );

    // compute operator matrices directly as the lpms were set correctly beforehand
    // Calculate A-matrix (first differential operator)
    this->ptPrimaryOp_->CalcOpMat( this->aMat_, lpmPrimary, ptFePrimary );
    // Calculate B-matrix (second differential operator)
    this->ptSecondaryOp_->CalcOpMat( this->bMat_, lpmSecondary, ptFeSecondary );
    // Calculate scalar factor...
    // if the object has been constructed passing the regionCoefs_ try to set the factors accordingly
    if (!this->regionCoefs_.empty()) {
      if (this->regionCoefs_.find(ptSurfSecondary->ptVolElems[0]->regionId) != this->regionCoefs_.end()) {
        this->regionCoefs_[ptSurfSecondary->ptVolElems[0]->regionId]->GetScalar(fac, lpmSecondary);
      } else if (this->regionCoefs_.find(ALL_REGIONS) != this->regionCoefs_.end()) {
        this->regionCoefs_[ALL_REGIONS]->GetScalar(fac, lpmSecondary);
      } else {
        EXCEPTION("Could not find coefficient function for region " << ent1.GetGrid()->GetRegion().ToString(ptSurfSecondary->ptVolElems[0]->regionId));
      }
    } else {
      this->coefScalar_->GetScalar(fac, lpmSecondary);
    }
    fac *= MAT_DATA_TYPE(lpmNc.jacDet * weights[iIntPts]);
#ifdef NDEBUG
    this->aMat_.Mult_Blas(this->bMat_, result, true, false, this->factor_*fac, 1.0);
#else
    LOG_DBG3(mortarIntMA) << "Integration point #" << iIntPts+1
        << "\n\tglobal coordinates: " << globIntPoint.ToString()
        << "\n\tlocal coordinates in primary: " << ipPrimary.coord.ToString()
        << "\n\tlocal coordinares in secondary: " << ipSecondary.coord.ToString();
    assert(esmPrimary->CoordIsInsideElem(ipPrimary.coord, EPS));
    assert(esmSecondary->CoordIsInsideElem(ipSecondary.coord, EPS));
    result += Transpose(this->aMat_) * this->bMat_ * this->factor_*fac;
    LOG_DBG3(mortarIntMA) << "Element matrix (not transposed): \n" << result;
#endif
  }
  if (doTranspose_) {
    result.Transpose(elemMat);
  } else {
    elemMat = result;
  }
}

template< class COEF_DATA_TYPE, class B_DATA_TYPE>
  void SurfaceMortarABIntMA<COEF_DATA_TYPE, B_DATA_TYPE>::
    SetFeSpace( shared_ptr<FeSpace> feSpace1, shared_ptr<FeSpace> feSpace2 ) {

  BBInt<COEF_DATA_TYPE, B_DATA_TYPE>::SetFeSpace(feSpace1, feSpace2);

  SolutionType resType1 = feSpace1->GetFeFunction().lock()->GetResultInfo()->resultType;
  SolutionType resType2 = feSpace2->GetFeFunction().lock()->GetResultInfo()->resultType;

  if (resType1 == MECH_DISPLACEMENT || resType1 ==  FLUIDMECH_VELOCITY) {
    ptFeSpacePrimary_ = feSpace1;
    ptPrimaryOp_ = this->aOperator_;
  } else if (resType2 == MECH_DISPLACEMENT || resType2 ==  FLUIDMECH_VELOCITY) {
    ptFeSpacePrimary_ = feSpace2;
    ptPrimaryOp_ = this->bOperator_;
    this->name_ = "SurfaceMortarABIntMAtrans";
    doTranspose_ = true;
  } else {
    EXCEPTION("Could not find FeSpace of mechanic region.");
  }
  if (resType1 == ACOU_POTENTIAL || resType1 == ACOU_PRESSURE) {
    ptFeSpaceSecondary_ = feSpace1;
    ptSecondaryOp_ = this->aOperator_;
  } else if (resType2 == ACOU_POTENTIAL || resType2 == ACOU_PRESSURE) {
    ptFeSpaceSecondary_ = feSpace2;
    ptSecondaryOp_ = this->bOperator_;
  } else {
    EXCEPTION("Could not find FeSpace of secondary region.");
  }
  assert(ptFeSpacePrimary_ != nullptr);
  assert(ptFeSpaceSecondary_ != nullptr);
  assert(ptPrimaryOp_ != nullptr);
  assert(ptSecondaryOp_ != nullptr);
}

// Explicit template instantiation
template class ABInt<Double,Double>;
template class ABInt<Double,Complex>;
template class ABInt<Complex,Double>;
template class ABInt<Complex,Complex>;
template class ABIntLem<Double,Double>;
template class ABIntLem<Double,Complex>;
template class ABIntLem<Complex,Double>;
template class ABIntLem<Complex,Complex>;
template class SurfaceABInt<Double,Double>;
template class SurfaceABInt<Double,Complex>;
template class SurfaceABInt<Complex,Double>;
template class SurfaceABInt<Complex,Complex>;
template class SurfaceMortarABInt<Double,Double>;
template class SurfaceMortarABInt<Double,Complex>;
template class SurfaceMortarABInt<Complex,Double>;
template class SurfaceMortarABInt<Complex,Complex>;
template class SurfaceNitscheABInt<Double,Double>;
template class SurfaceNitscheABInt<Double,Complex>;
template class SurfaceNitscheABInt<Complex,Double>;
template class SurfaceNitscheABInt<Complex,Complex>;
template class SurfaceMortarABIntMA<Double,Double>;
template class SurfaceMortarABIntMA<Double,Complex>;
template class SurfaceMortarABIntMA<Complex,Double>;
template class SurfaceMortarABIntMA<Complex,Complex>;

} // namespace CoupledField