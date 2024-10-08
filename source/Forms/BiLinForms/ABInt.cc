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
ABInt<COEF_DATA_TYPE, B_DATA_TYPE>
::ABInt(BaseBOperator * aOp, BaseBOperator * bOp, 
        PtrCoefFct scalCoef, MAT_DATA_TYPE factor,
        bool coordUpdate )
  : BBInt<COEF_DATA_TYPE, B_DATA_TYPE>(bOp, scalCoef, factor, coordUpdate )
{
  this->type_ = BiLinearForm::AB_INT;
  this->name_ = "ABInt"; // no idea, why this doesn't compile: type.ToString(type_); do
  this->aOperator_ = aOp;
  this->solDependent_ = false;
  this->isTimeFrequencyDependent_ = scalCoef->IsTimeFrequencyDependent();

  // Note: In general the AB-Integrator is not symmetric, as is should 
  // get only used, if the A- and B differential operators are distinct.
  this->isSymmetric_ = false;
}

template< class COEF_DATA_TYPE, class B_DATA_TYPE>
void ABInt<COEF_DATA_TYPE, B_DATA_TYPE>
::CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                     EntityIterator& ent1,
                     EntityIterator& ent2 )
{
  // Extract physical element
  const Elem* ptElem = ent1.GetElem();

  MAT_DATA_TYPE fac(0.0);

  // Obtain FE element from feSpace and integration scheme
  IntegOrder order1, order2;
  IntScheme::IntegMethod method1, method2;
  BaseFE* ptFeA = this->ptFeSpace1_->GetFe( ent1, method1, order1 );
  BaseFE* ptFeB = this->ptFeSpace2_->GetFe( ent1, method2, order2 );

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
  LocPointMapped lp;
  const UInt numIntPts = intPoints.GetSize();
  for( UInt i = 0; i < numIntPts; ++i  ) {

    // Calculate for each integration point the LocPointMapped
    lp.Set( intPoints[i], esm, weights[i] );

    // Calculate A-matrix (first differential operator)
    this->aOperator_->CalcOpMat( aMat_, lp, ptFeA );

    // Calculate B-matrix (second differential operator)
    this->bOperator_->CalcOpMat( this->bMat_, lp, ptFeB );

    // Calculate scalar factor
    this->coefScalar_->GetScalar(fac, lp);
    fac *= MAT_DATA_TYPE(lp.jacDet * weights[i]); 

#ifdef NDEBUG
    aMat_.Mult_Blas(this->bMat_,elemMat,true,false,this->factor_*fac,1.0);
#else
    elemMat += Transpose(aMat_) * this->bMat_ * this->factor_*fac;
#endif
  }
}

template< class COEF_DATA_TYPE, class B_DATA_TYPE>
SurfaceABInt<COEF_DATA_TYPE, B_DATA_TYPE>
::SurfaceABInt( BaseBOperator * aOp, BaseBOperator * bOp,
                PtrCoefFct scalCoef, MAT_DATA_TYPE factor,
                const std::set<RegionIdType>& volRegions,
                bool coordUpdate )
  : ABInt<COEF_DATA_TYPE,B_DATA_TYPE>(aOp, bOp, scalCoef, factor, coordUpdate)
{
  this->name_ = "SurfaceABInt";
  volRegions_ = volRegions;
  this->isSymmetric_ = false;
}

template< class COEF_DATA_TYPE, class B_DATA_TYPE>
SurfaceABInt<COEF_DATA_TYPE, B_DATA_TYPE>
::SurfaceABInt( BaseBOperator * aOp, BaseBOperator * bOp,
                const std::map< RegionIdType, PtrCoefFct >& regionCoefs,
                MAT_DATA_TYPE factor,
                const std::set<RegionIdType>& volRegions,
                bool coordUpdate)
  : ABInt<COEF_DATA_TYPE,B_DATA_TYPE>(aOp, bOp,regionCoefs.begin()->second, 
                                      factor, coordUpdate )
{
  this->name_ = "SurfaceABInt";
  volRegions_ = volRegions;
  this->isSymmetric_ = false;
  regionCoefs_ = regionCoefs;
}

template< class COEF_DATA_TYPE, class B_DATA_TYPE>
void SurfaceABInt<COEF_DATA_TYPE, B_DATA_TYPE>
::CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                     EntityIterator& ent1,
                     EntityIterator& ent2 )
{
  // Extract physical element
  const Elem* ptElem1 = ent1.GetElem();
  const Elem* ptElem2 = ent2.GetElem();

  MAT_DATA_TYPE fac(0.0);

  // Obtain FE element from feSpace and integration scheme

  const SurfElem *sElem1 = NULL;
  const SurfElem *sElem2 = NULL;
  const Elem *ptVolElem1 = NULL;
  const Elem *ptVolElem2 = NULL;

  BaseFE *ptFeA = NULL;
  BaseFE *ptFeB = NULL;
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
  LocPointMapped lp1,lp2;
  const UInt numIntPts = intPoints.GetSize();
  for( UInt i = 0; i < numIntPts; ++i ) {

    // Calculate for each integration point the LocPointMapped
    lp1.SetWithSurface( intPoints[i], esm1, volRegions_, weights[i] );
    lp2.SetWithSurface( intPoints[i], esm2, volRegions_, weights[i] );

    // Calculate A-matrix (first differential operator)
    this->aOperator_->CalcOpMat( this->aMat_, lp1, ptFeA );

    // Calculate B-matrix (second differential operator)
    this->bOperator_->CalcOpMat( this->bMat_, lp2, ptFeB );

    // Calculate scalar factor
    if(!regionCoefs_.empty()) {
      regionCoefs_[lp1.lpmVol->ptEl->regionId]->GetScalar(fac, lp1);
    } else {
      this->coefScalar_->GetScalar(fac, lp1);
    }

    fac *= MAT_DATA_TYPE(lp1.jacDet * weights[i]);

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
SurfaceNitscheABInt<COEF_DATA_TYPE, B_DATA_TYPE>
::SurfaceNitscheABInt( BaseBOperator * aOp, BaseBOperator * bOp,
                       PtrCoefFct scalCoef, MAT_DATA_TYPE factor,
                       BiLinearForm::CouplingDirection cplDir,
                       bool coordUpdate, bool isSym, bool isPenalty )
  : ABInt<COEF_DATA_TYPE,B_DATA_TYPE>(aOp, bOp, scalCoef, factor, coordUpdate)
{
  this->name_ = "SurfaceNitscheABInt";

  this->myDirection_ = cplDir;

  this->isSymmetric_ = isSym;
  this->isPenalty_ = isPenalty;
}


template< class COEF_DATA_TYPE, class B_DATA_TYPE>
void SurfaceNitscheABInt<COEF_DATA_TYPE, B_DATA_TYPE>
::CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                     EntityIterator& ent1,
                     EntityIterator& ent2 )
{
  //do now the following>
  //1. Perform some security checks
  //2. fill a volregion set according to the coupling direction
  //   in order to determine a consistent normal direction
  //3. obtain pointers to associated volume elements
  //4. loop over int points of surface element and do a standard integration
  //5. Hamideh: ptFeSpace1_ and ptFeSpace2_ points to Master and Slave unknowns.
  // The true feSpace for test function and  unknow will be taken from coupling direction
  // For example in the case of LinFlowMechNitsche coupling:
  // Master --> Fespace1_(mechanics displacement)
  // Slave --> FeSpace2_(LinFlow Velocity)

  const SurfElem* sElem1 = ent1.GetSurfElem();
  const SurfElem* sElem2 = ent2.GetSurfElem();
  //just to be sure, this context requires it1 and it2 to be identical
  if(sElem1->elemNum != sElem2->elemNum){
    EXCEPTION("SurfaceBiLinFormContext requires identical iterators")
  }

  //now lets try to downcast to MortarNcSurfElem
  const MortarNcSurfElem * mSe = dynamic_cast<const MortarNcSurfElem*>(sElem1);

  //this is more or less a switch case statement which terms to consider and
  // where to acquire the equation numbers
  bool useMaster1 = false,useMaster2 = false;
  this->GetVolFromSurfElem(useMaster1,useMaster2);
  const Elem* ptVolElem1 = (useMaster1) ? mSe->ptMaster->ptVolElems[0] : mSe->ptSlave->ptVolElems[0] ; // take pointer volume element of test function
  const Elem* ptVolElem2 = (useMaster2) ? mSe->ptMaster->ptVolElems[0] : mSe->ptSlave->ptVolElems[0] ; // take pointer volume element of unknown

  Matrix<MAT_DATA_TYPE> aMat, bMat;
  MAT_DATA_TYPE fac(0.0);

  // Integration order must be the same for all Nitsche integrators,
  // regardless of coupling direction.
  IntegOrder order1, order2;
  IntScheme::IntegMethod method1, method2;
  BaseFE* ptFeMaster = this->ptFeSpace1_->GetFe(mSe->ptMaster->ptVolElems[0]->elemNum);
  BaseFE* ptFeSlave =  this->ptFeSpace2_->GetFe(mSe->ptSlave->ptVolElems[0]->elemNum);
  this->ptFeSpace1_->GetIntegration(ptFeMaster, mSe->ptMaster->ptVolElems[0]->regionId, method1, order1);
  this->ptFeSpace2_->GetIntegration(ptFeSlave, mSe->ptSlave->ptVolElems[0]->regionId, method2, order2);

  // Obtain FE element from FeSpace for operators according to the coupling direction
  BaseFE* ptFeA = NULL;
  BaseFE* ptFeB = NULL;
  if (useMaster1 && useMaster2){
    ptFeA = this->ptFeSpace1_->GetFe(ptVolElem1->elemNum); // test function
    ptFeB = this->ptFeSpace1_->GetFe(ptVolElem2->elemNum); // unknown
  }else if(useMaster1 && !useMaster2){
    ptFeA = this->ptFeSpace1_->GetFe(ptVolElem1->elemNum); // test function
    ptFeB = this->ptFeSpace2_->GetFe(ptVolElem2->elemNum); // unknown
  } else if(!useMaster1 && useMaster2){
    ptFeA = this->ptFeSpace2_->GetFe(ptVolElem1->elemNum); // test function
    ptFeB = this->ptFeSpace1_->GetFe(ptVolElem2->elemNum); // unknown
  }else{
    ptFeA = this->ptFeSpace2_->GetFe(ptVolElem1->elemNum); // test function
    ptFeB = this->ptFeSpace2_->GetFe(ptVolElem2->elemNum); // unknown
  }

  const UInt nrFncsA = ptFeA->GetNumFncs();
  const UInt nrFncsB = ptFeB->GetNumFncs();

  // Get shape map from grid
  shared_ptr<ElemShapeMap> esm = ent1.GetGrid()->GetElemShapeMap( mSe, this->coordUpdate_ );

  // Get integration points
  StdVector<LocPoint> intPoints;
  StdVector<Double> weights;
  this->intScheme_->GetIntPoints( Elem::GetShapeType(mSe->type),
                                  method1, order1, method2, order2,
                                  intPoints, weights );

  elemMat.Resize( nrFncsA * this->aOperator_->GetDimDof(), nrFncsB * this->bOperator_->GetDimDof() );
  elemMat.Init();

  MAT_DATA_TYPE myFactor = this->factor_;
  if (this->isPenalty_) {
    shared_ptr<ElemShapeMap> esm1T = ent1.GetGrid()->GetElemShapeMap( mSe->ptMaster, this->coordUpdate_ );
    shared_ptr<ElemShapeMap> esm2T = ent1.GetGrid()->GetElemShapeMap( mSe->ptSlave, this->coordUpdate_ );

    //obtain pointer to basis functions
    BaseFE* SFe1 = this->ptFeSpace1_->GetFe(mSe->ptMaster->elemNum);
    BaseFE* SFe2 = this->ptFeSpace2_->GetFe(mSe->ptSlave->elemNum);

    MAT_DATA_TYPE tmp(2.0);
    Double min,max;
    esm1T->GetMaxMinEdgeLength(min,max);
    UInt order1 = SFe1->GetIsoOrder();

    if(order1 == 0){
      StdVector<UInt> aIsoOrder;
      SFe1->GetAnisoOrder(aIsoOrder);
      for(UInt i=0;i<aIsoOrder.GetSize();i++)
        order1 = (aIsoOrder[i]>order1)? aIsoOrder[i] : order1;

      if(order1 == 0 && this->ptFeSpace1_->GetSpaceType() == this->ptFeSpace1_->HCURL){
        // that's a bit dirty because zero order Nedelec elements are neither zero nor first order...
        // but here we set it to first, otherwise we divide by zero
        order1 = 1;
      }
    }
    MAT_DATA_TYPE surface1(min/(Double)order1);

    esm2T->GetMaxMinEdgeLength(min,max);
    UInt order2 = SFe2->GetIsoOrder();

    if(order2 == 0){
      StdVector<UInt> aIsoOrder;
      SFe2->GetAnisoOrder(aIsoOrder);
      for(UInt i=0;i<aIsoOrder.GetSize();i++)
        order2 = (aIsoOrder[i]>order2)? aIsoOrder[i] : order2;

      if(order2 == 0 && this->ptFeSpace1_->GetSpaceType() == this->ptFeSpace1_->HCURL){
        // that's a bit dirty because zero order Nedelec elements are neither zero nor first order...
        // but here we set it to first, otherwise we divide by zero
        order2 = 1;
      }
    }
    MAT_DATA_TYPE surface2(min/(Double)order2);

    myFactor *= tmp/(surface1+surface2);
  }

  // Loop over all integration points
  LocPointMapped lp1, lp2;
  const UInt numIntPts = intPoints.GetSize();
  for( UInt i = 0; i < numIntPts; ++i ) {

    // Calculate for each integration point the LocPointMapped
    lp1.SetMortar( intPoints[i], esm, weights[i],useMaster1 );
    lp2.SetMortar( intPoints[i], esm, weights[i],useMaster2 );

    // Calculate A-matrix (first differential operator)
    this->aOperator_->CalcOpMat( this->aMat_, lp1, ptFeA );

    // Calculate B-matrix (second differential operator)
    this->bOperator_->CalcOpMat( this->bMat_, lp2, ptFeB );

    // Calculate scalar factor
    this->coefScalar_->GetScalar(fac, lp1);

    LOG_DBG(nitscheInt) << "fac in region with ID " << lp1.ptEl->regionId <<" = "<<fac;
    LOG_DBG(nitscheInt) << "myFactor without weights = "<<myFactor;

    fac *= MAT_DATA_TYPE(lp1.jacDet * weights[i]);


#ifdef NDEBUG
    this->aMat_.Mult_Blas(this->bMat_, elemMat, true, false,
        myFactor*fac, 1.0);
#else
    elemMat += Transpose(this->aMat_) * this->bMat_ * myFactor *fac;
#endif
  }

}

template< class COEF_DATA_TYPE, class B_DATA_TYPE>
void SurfaceNitscheABInt<COEF_DATA_TYPE, B_DATA_TYPE>
::GetVolFromSurfElem(bool & uMaster1, bool & uMaster2)
{
  switch(this->myDirection_){
  case BiLinearForm::MASTER_MASTER:
    uMaster1 = true;
    uMaster2 = true;
    break;
  case BiLinearForm::SLAVE_SLAVE:
    uMaster1 = false;
    uMaster2 = false;
    break;
  case BiLinearForm::MASTER_SLAVE:
    uMaster1 = true;
    uMaster2 = false;
    break;
  case BiLinearForm::SLAVE_MASTER:
    uMaster1 = false;
    uMaster2 = true;
    break;
  default:
    EXCEPTION("Undefined coupling direction");
    break;
  }
}

DEFINE_LOG(mortarInt, "Forms.MortarInt")

template< class COEF_DATA_TYPE, class B_DATA_TYPE>
SurfaceMortarABInt<COEF_DATA_TYPE, B_DATA_TYPE>
::SurfaceMortarABInt( BaseBOperator * aOp, BaseBOperator * bOp,
                      PtrCoefFct scalCoef, MAT_DATA_TYPE factor,
                      RegionIdType masterVolRegion,
                      RegionIdType slaveVolRegion,
                      bool coplanar,
                      bool coordUpdate,
                      BiLinearForm::CouplingDirection cplDirection)
: ABInt<COEF_DATA_TYPE, B_DATA_TYPE>( aOp, bOp, scalCoef, factor, coordUpdate),
  isCoplanar_(coplanar), cplDirection_(cplDirection)
{
  this->name_ = "SurfaceMortarABInt";
  this->isSymmetric_ = true;
  
  this->masterVolRegion_ = masterVolRegion;
  this->slaveVolRegion_ = slaveVolRegion;
  
  // create set of volume regions
  volRegions_.insert(masterVolRegion_);
  volRegions_.insert(slaveVolRegion_);
}

template< class COEF_DATA_TYPE, class B_DATA_TYPE>
void SurfaceMortarABInt<COEF_DATA_TYPE, B_DATA_TYPE>
::CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                     EntityIterator& ent1,
                     EntityIterator& ent2 )
{
  assert( ent1.GetElem() == ent2.GetElem() );

  // Extract MortarNcSurfElem element
  const NcSurfElem* ptNcElem = ent1.GetNcSurfElem();
  const MortarNcSurfElem * ptMortarElem =
      dynamic_cast<const MortarNcSurfElem*>(ptNcElem);
  assert( ptMortarElem );

  const SurfElem *ptSurfMaster = ptMortarElem->ptMaster;
  const SurfElem *ptSurfSlave = ptMortarElem->ptSlave;
  assert(ptSurfMaster);
  assert(ptSurfSlave);

  LOG_DBG(mortarInt) << "NcSurfElem #" << ptMortarElem->elemNum
                     << ", type " << Elem::feType.ToString(ptMortarElem->type)
                     << "\n\tmaster parent #" << ptSurfMaster->elemNum
                     << "\n\tslave parent #" << ptSurfSlave->elemNum;
  LOG_DBG2(mortarInt) << "\tconnectivity: " << ptMortarElem->connect.ToString();
  
  MAT_DATA_TYPE fac(0.0);

  // Obtain FE element from feSpace and integration scheme
  BaseFE* ptFeMaster = this->ptFeSpaceField_->GetFe(ptSurfMaster->elemNum);
  BaseFE* ptFeSlave = this->ptFeSpaceLM_->GetFe(ptSurfSlave->elemNum);
  assert(ptFeMaster);
  assert(ptFeSlave);

  IntegOrder orderMaster, orderSlave;
  IntScheme::IntegMethod methodMaster, methodSlave;
  this->ptFeSpaceField_->GetIntegration( ptFeMaster, masterVolRegion_, methodMaster, orderMaster );
  this->ptFeSpaceField_->GetIntegration( ptFeSlave, slaveVolRegion_, methodSlave, orderSlave );

  const UInt nrFncsMaster = ptFeMaster->GetNumFncs();
  const UInt nrFncsSlave = ptFeSlave->GetNumFncs();

  // Get shape map from grid
  shared_ptr<ElemShapeMap> esmNc = ent1.GetGrid()->GetElemShapeMap( ptMortarElem, this->coordUpdate_ );
  shared_ptr<ElemShapeMap> esmMaster = ent1.GetGrid()->GetElemShapeMap( ptSurfMaster, this->coordUpdate_ );
  shared_ptr<ElemShapeMap> esmSlave = ent1.GetGrid()->GetElemShapeMap( ptSurfSlave, this->coordUpdate_ );

  // Get integration points
  StdVector<LocPoint> intPoints;
  StdVector<Double> weights;
  this->intScheme_->GetIntPoints( Elem::GetShapeType(ptMortarElem->type), 
      methodMaster, orderMaster, methodSlave, orderSlave, intPoints, weights );

  // initialize element matrix
  elemMat.Resize( nrFncsMaster * this->aOperator_->GetDimDof(), nrFncsSlave * this->bOperator_->GetDimDof() );
  elemMat.Init();

  // Loop over all integration points
  Vector<Double> globIntPoint;
  LocPoint ipMaster, ipSlave;
  LocPointMapped lpmNc, lpmMaster, lpmSlave;

  const UInt numIntPts = intPoints.GetSize();
  for( UInt i = 0; i < numIntPts; ++i ) {
    // Calculate global coordinates of integration point
    esmNc->Local2Global(globIntPoint, intPoints[i]);

    // Calculate local coordinates of integration point in slave element
    esmSlave->Global2Local(ipSlave.coord, globIntPoint);
    assert( esmSlave->CoordIsInsideElem(ipSlave.coord) );

    // Projection for curved interfaces
    if ( !isCoplanar_ ) {
      /* Instead of projecting the integration point back into the master
       * element, we use a better approach: Orthogonal projection is a linear
       * operation, so any point in the projected master element has the same
       * local coordinates as its orthogonal projection in the original master
       * element. So we perform the Global2Local mapping in the projected
       * element, which has been computed already by the intersection
       * algorithm. Sometimes the intersection is the projected master element
       * itself, in which case we need no coordinate mapping at all.  
       */ 
      if ( ptMortarElem->projectedMaster ) {
        shared_ptr<ElemShapeMap> esmProj = ent1.GetGrid()->GetElemShapeMap(
            ptMortarElem->projectedMaster.get(), this->coordUpdate_ );
        esmProj->Global2Local(ipMaster.coord, globIntPoint);
      } else { // projectedMaster == NULL means it is the MortarElem itself
        ipMaster.coord = intPoints[i].coord;
      }
    } else {
      // Calculate local coordinates of integration point in master element
      esmMaster->Global2Local(ipMaster.coord, globIntPoint);
    }
    assert( esmMaster->CoordIsInsideElem(ipMaster.coord) );

    LOG_DBG3(mortarInt) << "Integration point #" << i+1
        << "\n\tglobal coordinates: " << globIntPoint.ToString()
        << "\n\tlocal coordinates in master: " << ipMaster.coord.ToString()
        << "\n\tlocal coordinares in slave: " << ipSlave.coord.ToString();

    // Calculate for each integration point the LocPointMapped
    lpmNc.Set( intPoints[i], esmNc, weights[i] );
    lpmMaster.SetWithSurface( ipMaster, esmMaster, this->volRegions_, weights[i] );
    lpmSlave.SetWithSurface( ipSlave, esmSlave, this->volRegions_, weights[i] );

    // Calculate A-matrix (first differential operator)
    // Calculate B-matrix (second differential operator)
    if (cplDirection_ == BiLinearForm::MASTER_SLAVE)
    {
      this->aOperator_->CalcOpMat( this->aMat_, lpmMaster, ptFeMaster );
      this->bOperator_->CalcOpMat( this->bMat_, lpmSlave, ptFeSlave );
    }
    else if (cplDirection_ == BiLinearForm::SLAVE_MASTER)
    {
      this->aOperator_->CalcOpMat( this->aMat_, lpmSlave, ptFeSlave );
      this->bOperator_->CalcOpMat( this->bMat_, lpmMaster, ptFeMaster );
    }
    else
    {
      EXCEPTION("Cannot calculate operator matrices for the current coupling direction: " + boost::lexical_cast<std::string>(cplDirection_));
    }

    // Calculate scalar factor
    this->coefScalar_->GetScalar(fac, lpmSlave);

    fac *= MAT_DATA_TYPE(lpmNc.jacDet * weights[i]);

#ifdef NDEBUG
    this->aMat_.Mult_Blas(this->bMat_, elemMat, true, false,
        this->factor_*fac, 1.0);
#else
    elemMat += Transpose(this->aMat_) * this->bMat_ * this->factor_*fac;
#endif
  }
  
  LOG_DBG2(mortarInt) << "Element matrix of NcSurfElem #"
                      << ptMortarElem->elemNum << ":\n" << elemMat;
}

template< class COEF_DATA_TYPE, class B_DATA_TYPE>
void SurfaceMortarABInt<COEF_DATA_TYPE, B_DATA_TYPE>
::SetFeSpace( shared_ptr<FeSpace> feSpace1, shared_ptr<FeSpace> feSpace2 ) {
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
}

template< class COEF_DATA_TYPE, class B_DATA_TYPE>
SurfaceMortarABIntMA<COEF_DATA_TYPE, B_DATA_TYPE>
::SurfaceMortarABIntMA( BaseBOperator * aOp, BaseBOperator * bOp,
                      PtrCoefFct scalCoef, MAT_DATA_TYPE factor,
					  bool coplanar,
                      bool coordUpdate)
: SurfaceABInt<COEF_DATA_TYPE, B_DATA_TYPE>( aOp, bOp, scalCoef, factor,
		                                     std::set<RegionIdType>(), coordUpdate)
//  isCoplanar_(coplanar), cplDirection_(cplDirection)
{
  this->isSymmetric_ = false;
  this->isCoplanar_ = coplanar;
  this->ptMasterOp_ = NULL;
  this->ptSlaveOp_ = NULL;
  doTranspose_ = false;
}

//mechanical-acoustic coupling on non-conforming grids!
template<class COEF_DATA_TYPE, class B_DATA_TYPE>
SurfaceMortarABIntMA<COEF_DATA_TYPE, B_DATA_TYPE>
::SurfaceMortarABIntMA( BaseBOperator * aOp, BaseBOperator * bOp,
                      const std::map< RegionIdType, PtrCoefFct >& regionCoefs,
                      MAT_DATA_TYPE factor,
                      bool coplanar,
                      bool coordUpdate)
: SurfaceABInt<COEF_DATA_TYPE, B_DATA_TYPE>(aOp, bOp, regionCoefs, factor,
                                     std::set<RegionIdType>(), coordUpdate)
{
  this->name_ = "SurfaceMortarABIntMA";
  this->isSymmetric_ = false;
  this->regionCoefs_ = regionCoefs;
  this->isCoplanar_  = coplanar;
  this->ptMasterOp_  = NULL;
  this->ptSlaveOp_   = NULL;
  doTranspose_ = false;
}


template< class COEF_DATA_TYPE, class B_DATA_TYPE>
void SurfaceMortarABIntMA<COEF_DATA_TYPE, B_DATA_TYPE>
::CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                     EntityIterator& ent1,
                     EntityIterator& ent2 ) {

  // Make sure both entities are the same
  assert( ent1.GetElem() == ent2.GetElem() );

  // Make sure both master and slave FeSpaces have been assigned
  assert(ptFeSpaceMaster_);
  assert(ptFeSpaceSlave_);

  // Make sure both master and slave operators have been assigned
  assert(ptMasterOp_);
  assert(ptSlaveOp_);

  // Extract MortarNcSurfElem element
  const NcSurfElem* ptNcElem = ent1.GetNcSurfElem();
  const MortarNcSurfElem * ptMortarElem = dynamic_cast<const MortarNcSurfElem*>(ptNcElem);
  assert( ptMortarElem );


  // Obtain master and slave parents of Mortar element
  const SurfElem *ptSurfMaster = ptMortarElem->ptMaster;
  const SurfElem *ptSurfSlave = ptMortarElem->ptSlave;
  assert(ptSurfMaster);
  assert(ptSurfSlave);

  assert(ptSurfMaster->ptVolElems[0]);
  assert(ptSurfSlave->ptVolElems[0]);
  RegionIdType masterVolRegion = ptSurfMaster->ptVolElems[0]->regionId;
  RegionIdType slaveVolRegion = ptSurfSlave->ptVolElems[0]->regionId;
  this->volRegions_.insert(masterVolRegion);
  this->volRegions_.insert(slaveVolRegion);

  LOG_DBG(mortarInt) << "NcSurfElem #" << ptMortarElem->elemNum
                     << ", type " << Elem::feType.ToString(ptMortarElem->type)
                     << "\n\tmaster parent #" << ptSurfMaster->elemNum
                     << "\n\tslave parent #" << ptSurfSlave->elemNum;
  LOG_DBG2(mortarInt) << "\tconnectivity: " << ptMortarElem->connect.ToString();

  MAT_DATA_TYPE fac(0.0);

  // Obtain FE element from feSpace and integration scheme
  BaseFE *ptFeMaster, *ptFeSlave;
  ptFeMaster = this->ptFeSpaceMaster_->GetFe(ptSurfMaster->elemNum);
  ptFeSlave = this->ptFeSpaceSlave_->GetFe(ptSurfSlave->elemNum);

  assert(ptFeMaster);
  assert(ptFeSlave);

  // Obtain integration schemes
  IntegOrder orderMaster, orderSlave;
  IntScheme::IntegMethod methodMaster, methodSlave;
  this->ptFeSpaceMaster_->GetIntegration( ptFeMaster, masterVolRegion, methodMaster, orderMaster );
  this->ptFeSpaceSlave_->GetIntegration( ptFeSlave, slaveVolRegion, methodSlave, orderSlave );

  const UInt nrFncsMaster = ptFeMaster->GetNumFncs();
  const UInt nrFncsSlave = ptFeSlave->GetNumFncs();

  // Get shape map from grid
  shared_ptr<ElemShapeMap> esmNc = ent1.GetGrid()->GetElemShapeMap( ptMortarElem, this->coordUpdate_);
  shared_ptr<ElemShapeMap> esmMaster = ent1.GetGrid()->GetElemShapeMap( ptSurfMaster, this->coordUpdate_);
  shared_ptr<ElemShapeMap> esmSlave = ent1.GetGrid()->GetElemShapeMap( ptSurfSlave, this->coordUpdate_ );

  // Get integration points
  StdVector<LocPoint> intPoints;
  StdVector<Double> weights;

  this->intScheme_->GetIntPoints( Elem::GetShapeType(ptMortarElem->type),
      methodMaster, orderMaster, methodSlave, orderSlave, intPoints, weights );

  // initialize element matrix
  Matrix<MAT_DATA_TYPE> result;
  result.Resize( nrFncsMaster * this->ptMasterOp_->GetDimDof(), nrFncsSlave * this->ptSlaveOp_->GetDimDof() );
  result.Init();

  // Loop over all integration points
  Vector<Double> globIntPoint;
  LocPoint ipMaster, ipSlave;
  LocPointMapped lpmNc, lpmMaster, lpmSlave;

  const UInt numIntPts = intPoints.GetSize();
  for( UInt i = 0; i < numIntPts; ++i ) {
    // Calculate global coordinates of integration point
    esmNc->Local2Global(globIntPoint, intPoints[i]);

    // Calculate local coordinates of integration point in slave element
    esmSlave->Global2Local(ipSlave.coord, globIntPoint);
    assert( esmSlave->CoordIsInsideElem(ipSlave.coord) );

    // Projection for curved interfaces
    if ( !isCoplanar_ ) {
      /* Instead of projecting the integration point back into the master
       * element, we use a better approach: Orthogonal projection is a linear
       * operation, so any point in the projected master element has the same
       * local coordinates as its orthogonal projection in the original master
       * element. So we perform the Global2Local mapping in the projected
       * element, which has been computed already by the intersection
       * algorithm. Sometimes the intersection is the projected master element
       * itself, in which case we need no coordinate mapping at all.
       */
      if ( ptMortarElem->projectedMaster ) {
        shared_ptr<ElemShapeMap> esmProj = ent1.GetGrid()->GetElemShapeMap(
            ptMortarElem->projectedMaster.get(), this->coordUpdate_ );
        esmProj->Global2Local(ipMaster.coord, globIntPoint);
      } else { // projectedMaster == NULL means it is the MortarElem itself
        ipMaster.coord = intPoints[i].coord;
      }
    } else {
      // Calculate local coordinates of integration point in master element
      esmMaster->Global2Local(ipMaster.coord, globIntPoint);
    }
    assert( esmMaster->CoordIsInsideElem(ipMaster.coord) );

    LOG_DBG3(mortarInt) << "Integration point #" << i+1
        << "\n\tglobal coordinates: " << globIntPoint.ToString()
        << "\n\tlocal coordinates in master: " << ipMaster.coord.ToString()
        << "\n\tlocal coordinares in slave: " << ipSlave.coord.ToString();

    // Calculate for each integration point the LocPointMapped
    lpmNc.Set( intPoints[i], esmNc, weights[i] );
    lpmMaster.SetWithSurface( ipMaster, esmMaster, this->volRegions_, weights[i] );
    lpmSlave.SetWithSurface( ipSlave, esmSlave, this->volRegions_, weights[i] );

    // Calculate A-matrix (first differential operator)
    this->ptMasterOp_->CalcOpMat( this->aMat_, lpmMaster, ptFeMaster );

    // Calculate B-matrix (second differential operator)
    this->ptSlaveOp_->CalcOpMat( this->bMat_, lpmSlave, ptFeSlave );

    // Calculate scalar factor
    if (!this->regionCoefs_.empty()) {
      if (this->regionCoefs_.find(slaveVolRegion) != this->regionCoefs_.end()) {
        this->regionCoefs_[slaveVolRegion]->GetScalar(fac, lpmSlave);
      } else if (this->regionCoefs_.find(ALL_REGIONS) != this->regionCoefs_.end()) {
        this->regionCoefs_[ALL_REGIONS]->GetScalar(fac, lpmSlave);
      } else {
        EXCEPTION("Could not find coefficient function for region "
                  << ent1.GetGrid()->GetRegion().ToString(slaveVolRegion));
      }
    } else {
      this->coefScalar_->GetScalar(fac, lpmSlave);
    }

    fac *= MAT_DATA_TYPE(lpmNc.jacDet * weights[i]);

#ifdef NDEBUG
    this->aMat_.Mult_Blas(this->bMat_, result, true, false,
        this->factor_*fac, 1.0);
#else
    result += Transpose(this->aMat_) * this->bMat_ * this->factor_*fac;
#endif
  }
  if ( doTranspose_ ) {
	  result.Transpose(elemMat);
  }
  else {
	  elemMat = result;
//  if ( ptMasterOp_ != this->bOperator_) {
//	  elemMat = result;
//  } else {
//	  result.Transpose(elemMat);
  }

  LOG_DBG2(mortarInt) << "Element matrix of NcSurfElem #"
                      << ptMortarElem->elemNum << ":\n" << elemMat;
}

template< class COEF_DATA_TYPE, class B_DATA_TYPE>
void SurfaceMortarABIntMA<COEF_DATA_TYPE, B_DATA_TYPE>
::SetFeSpace( shared_ptr<FeSpace> feSpace1, shared_ptr<FeSpace> feSpace2 ) {

  BBInt<COEF_DATA_TYPE, B_DATA_TYPE>::SetFeSpace(feSpace1, feSpace2);

  SolutionType resType1 = feSpace1->GetFeFunction().lock()->GetResultInfo()->resultType;
  SolutionType resType2 = feSpace2->GetFeFunction().lock()->GetResultInfo()->resultType;

  if (resType1 == MECH_DISPLACEMENT || resType1 ==  FLUIDMECH_VELOCITY) {
	  ptFeSpaceMaster_ = feSpace1;
	  ptMasterOp_ = this->aOperator_;
  } else if (resType2 == MECH_DISPLACEMENT || resType2 ==  FLUIDMECH_VELOCITY) {
	  ptFeSpaceMaster_ = feSpace2;
	  ptMasterOp_ = this->bOperator_;
	  this->name_ = "SurfaceMortarABIntMAtrans";
	  doTranspose_ = true;
  } else {
	  EXCEPTION("Could not find FeSpace of mechanic region.");
  }
  if (resType1 == ACOU_POTENTIAL || resType1 == ACOU_PRESSURE) {
	  ptFeSpaceSlave_ = feSpace1;
	  ptSlaveOp_ = this->aOperator_;
  } else if (resType2 == ACOU_POTENTIAL || resType2 == ACOU_PRESSURE) {
	  ptFeSpaceSlave_ = feSpace2;
	  ptSlaveOp_ = this->bOperator_;
  } else {
	  EXCEPTION("Could not find FeSpace of slave region.");
  }
  assert(ptFeSpaceMaster_);
  assert(ptFeSpaceSlave_);
}

// Explicit template instantiation
template class ABInt<Double,Double>;
template class ABInt<Double,Complex>;
template class ABInt<Complex,Double>;
template class ABInt<Complex,Complex>;
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
