#include "ADBInt.hh"

namespace CoupledField{

 template< class COEF_DATA_TYPE, class B_DATA_TYPE >
  ADBInt<COEF_DATA_TYPE, B_DATA_TYPE>::
  ADBInt(BaseBOperator * aOp, BaseBOperator * bOp,
         PtrCoefFct dData, MAT_DATA_TYPE factor,
         bool coordUpdate )
  : BDBInt<COEF_DATA_TYPE, B_DATA_TYPE>(bOp, dData, factor, coordUpdate  )  {
      this->type_ = BiLinearForm::ADB_INT;
      this->name_ = "ADBInt";
      this->isSymmetric_ = false; // in general the ADB
      this->isTimeFrequencyDependent_ = dData->IsTimeFrequencyDependent();
      assert(dData->GetDimType() == CoefFunction::TENSOR);
#ifndef NDEBUG
      if(dData->GetDimType() != CoefFunction::TENSOR){
        Exception("BDB integrator expects the coefficient function to be tensorial");
      }
#endif
      this->aOperator_ = aOp;
      this->factor_ = factor;
      
  }

  //! Destructor
 template< class COEF_DATA_TYPE, class B_DATA_TYPE >
   ADBInt<COEF_DATA_TYPE, B_DATA_TYPE>::
  ~ADBInt(){
   delete aOperator_;
  }

 template< class COEF_DATA_TYPE, class B_DATA_TYPE >
  void  ADBInt<COEF_DATA_TYPE, B_DATA_TYPE>::
  CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                     EntityIterator& ent1,
                     EntityIterator& ent2) {

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
    shared_ptr<ElemShapeMap> esm =
        ent1.GetGrid()->GetElemShapeMap( ptElem, this->coordUpdate_ );

    // Get integration points
    StdVector<LocPoint> intPoints;
    StdVector<Double> weights;
    this->intScheme_->GetIntPoints( Elem::GetShapeType(ptElem->type), 
                                    method1, order1, method2, order2,
                                    intPoints, weights );

    elemMat.Resize( nrFncsA * aOperator_->GetDimDof(), 
                    nrFncsB * this->bOperator_->GetDimDof() );
    elemMat.Init();
    
    // Loop over all integration points
    LocPointMapped lp;
    const UInt numIntPts = intPoints.GetSize();
    for( UInt i = 0; i < numIntPts; ++i) {

      // Calculate for each integration point the LocPointMapped
      lp.Set( intPoints[i], esm, weights[i] );

      // Calculate the A-matrix
      aOperator_->CalcOpMat( aMat_, lp, ptFeA);
      
      // Calculate the B-matrix
      this->bOperator_->CalcOpMat( this->bMat_, lp, ptFeB);

      // Calculate D-Mat
      this->dData_->GetTensor(this->dMat_,lp);
      
      fac = MAT_DATA_TYPE(lp.jacDet * weights[i]);
      
      this->dbMat_.Resize(this->dMat_.GetNumRows(), 
                          nrFncsB * this->bOperator_->GetDimDof() );

#ifdef NDEBUG
      this->dMat_.Mult_Blas(this->bMat_,this->dbMat_,false,false,1.0,0);
      aMat_.Mult_Blas(this->dbMat_,elemMat,true,false,this->factor_*fac,1.0);
#else
      this->dbMat_ = (this->dMat_ * this->bMat_) * fac;
      elemMat += Transpose(this->aMat_) * this->dbMat_ * this->factor_;
#endif
    }
  }
 
 // ===========================
 //  Apply transposed A Matrix
 // ============================
 // 1) General template for double-case -> not implemented
 template< class COEF_DATA_TYPE, class B_DATA_TYPE >
 void ADBInt<COEF_DATA_TYPE, B_DATA_TYPE>::
 ApplyATransMat( Vector<Double>&ret, 
                 const Vector<Double>& sol,
                 const LocPointMapped& lpm ) {
   EXCEPTION( "Not implemented" );
 }
 // 2) Special double-only case
 template<>
 void ADBInt<Double, Double>::
 ApplyATransMat( Vector<Double>&ret, 
                 const Vector<Double>& sol,
                 const LocPointMapped& lpm ) {
   BaseFE* ptFe = this->ptFeSpace1_->GetFe( lpm.ptEl->elemNum );
   aOperator_->CalcOpMat(aMat_, lpm, ptFe);
   ret = Transpose(aMat_) * sol;
 }

 // 3) Complex case is always possible
 template< class COEF_DATA_TYPE, class B_DATA_TYPE >
 void ADBInt<COEF_DATA_TYPE, B_DATA_TYPE>::
 ApplyATransMat( Vector<Complex>&ret, 
                 const Vector<Complex>& sol,
                 const LocPointMapped& lpm ) {
   BaseFE* ptFe = this->ptFeSpace1_->GetFe( lpm.ptEl->elemNum );
   aOperator_->CalcOpMat(aMat_, lpm, ptFe);
   ret = Transpose(aMat_) * sol;
 }
 
 // ============================
 //  Apply transposed dA Matrix
 // =============================
 // 1) General template for double-case -> not implemented
 template< class COEF_DATA_TYPE, class B_DATA_TYPE >
 void ADBInt<COEF_DATA_TYPE, B_DATA_TYPE>::
 ApplydATransMat( Vector<Double>&ret, 
                  const Vector<Double>& sol,
                  const LocPointMapped& lpm ) {
   EXCEPTION( "Not implemented" );
 }
 // 2) Special double-only case
  template<>
 void ADBInt<Double, Double>::
 ApplydATransMat( Vector<Double>&ret, 
                  const Vector<Double>& sol,
                  const LocPointMapped& lpm ) {
   Matrix<MAT_DATA_TYPE> aOp, dMat, dAMat;
   BaseFE* ptFe = this->ptFeSpace1_->GetFe( lpm.ptEl->elemNum );
   aOperator_->CalcOpMat(aMat_, lpm, ptFe);
   this->dData_->GetTensor(this->dMat_,lpm);
   dAMat_.Resize(this->dMat_.GetNumCols(), aMat_.GetNumCols()  );
   this->dMat_.Mult_Blas(aMat_,dAMat_,true,false,1.0,0);
   ret = dAMat_* sol;
 }
 
  // 3) Complex case is always possible
 template< class COEF_DATA_TYPE, class B_DATA_TYPE >
 void ADBInt<COEF_DATA_TYPE, B_DATA_TYPE>::
 ApplydATransMat( Vector<Complex>&ret, 
                  const Vector<Complex>& sol,
                  const LocPointMapped& lpm ) {
   Matrix<MAT_DATA_TYPE> aOp, dMat, dAMat;
   BaseFE* ptFe = this->ptFeSpace1_->GetFe( lpm.ptEl->elemNum );
   aOperator_->CalcOpMat(aMat_, lpm, ptFe);
   this->dData_->GetTensor(this->dMat_,lpm);
   dAMat_.Resize(this->dMat_.GetNumCols(), aMat_.GetNumCols()  );
   this->dMat_.Mult_Blas(aMat_,dAMat_,true,false,1.0,0);
   ret = dAMat_* sol;
 }
 
 
  //! Calculate the integration kernel
 template< class COEF_DATA_TYPE, class B_DATA_TYPE >
 void ADBInt<COEF_DATA_TYPE, B_DATA_TYPE>::
 CalcKernel( Matrix<MAT_DATA_TYPE>& kernel, 
             const LocPointMapped& lpm ) {

    // Obtain FE element from feSpace and integration scheme
    BaseFE* ptFeA = this->ptFeSpace1_->GetFe( lpm.ptEl->elemNum );
    BaseFE* ptFeB = this->ptFeSpace2_->GetFe( lpm.ptEl->elemNum );

    const UInt nrFncsA = ptFeA->GetNumFncs();
    const UInt nrFncsB = ptFeB->GetNumFncs();
    kernel.Resize( nrFncsA * this->aOperator_->GetDimDof(), 
                   nrFncsB * this->bOperator_->GetDimDof() );
    kernel.Init();


    // Calculate the A-matrix
    aOperator_->CalcOpMat( aMat_, lpm, ptFeA);
    
    // Calculate the B-matrix
    this->bOperator_->CalcOpMat( this->bMat_, lpm, ptFeB);

    // Calculate D-Mat
    this->dData_->GetTensor(this->dMat_,lpm);
    this->dbMat_.Resize(this->dMat_.GetNumRows(),
                        nrFncsB * this->bOperator_->GetDimDof() );

#ifdef NDEBUG
    this->dMat_.Mult_Blas(this->bMat_, this->dbMat_, false, false, 1.0, 0);
    aMat_.Mult_Blas(this->dbMat_, kernel, true, false, this->factor_,1.0);
#else
    this->dbMat_ = (this->dMat_ * this->bMat_);
    kernel += Transpose(this->aMat_) * this->dbMat_ * this->factor_;
#endif

  }
  
 // Explicit template instantiation
  template class ADBInt<Double,Double>;
  template class ADBInt<Double,Complex>;
  template class ADBInt<Complex,Double>;
  template class ADBInt<Complex,Complex>;
}
