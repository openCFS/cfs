// =====================================================================================
// 
//       Filename:  bdbInt.cc
// 
//    Description:  Implementation file for bdbIntegrators
//                  TAKE care:
//                  This file is just for inclusion in the header file!
// 
//        Version:  1.0
//        Created:  11/03/2011 07:36:28 PM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Andreas Hueppe (AHU), andreas.hueppe@uni-klu.ac.at
//        Company:  Universitaet Klagenfurt
// 
// =====================================================================================
#include "BDBInt.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

DECLARE_LOG(bdbint)
DEFINE_LOG(bdbint, "bdbint")

namespace CoupledField{

  template< class COEF_DATA_TYPE,
            class B_DATA_TYPE>
  BDBInt<COEF_DATA_TYPE,B_DATA_TYPE>::
  BDBInt(BaseBOperator* bOp, PtrCoefFct dData, MAT_DATA_TYPE factor, bool coordUpdate )
  : BaseBDBInt(coordUpdate) 
    {
      name_ = "BDBInt";
      isSymmetric_ = true;

      //assert(dData->GetDimType() == CoefFunction::TENSOR);
#ifndef NDEBUG
      if(dData->GetDimType() != CoefFunction::TENSOR){
        Exception("BDB integrator expects the coefficient function to be tensorial");
      }
#endif
      bOperator_ = bOp;
      dData_ = dData;
      factor_ = factor;
  }

  //! Destructor
  template< class COEF_DATA_TYPE,
            class B_DATA_TYPE>
  BDBInt<COEF_DATA_TYPE,B_DATA_TYPE>::~BDBInt(){
    delete this->bOperator_;
  }

  //! Calculate the element matrix
  template< class COEF_DATA_TYPE,
            class B_DATA_TYPE>
  void BDBInt<COEF_DATA_TYPE,B_DATA_TYPE>::CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat, EntityIterator& ent1, EntityIterator& ent2) {

    // Extract physical element
    const Elem* ptElem = ent1.GetElem();

    MAT_DATA_TYPE fac(0.0);

    // Obtain FE element from feSpace and integration scheme
    IntegOrder order;
    IntScheme::IntegMethod method;
    BaseFE* ptFe = ptFeSpace1_->GetFe( ent1, method, order );

    const UInt nrFncs = ptFe->GetNumFncs();

    // Get shape map from grid
    shared_ptr<ElemShapeMap> esm =
        ent1.GetGrid()->GetElemShapeMap( ptElem, this->coordUpdate_ );

    // Get integration points
    StdVector<LocPoint> intPoints;
    StdVector<Double> weights;
    intScheme_->GetIntPoints( Elem::GetShapeType(ptElem->type), method, order, 
                              intPoints, weights );

    elemMat.Resize( nrFncs * bOperator_->GetDimDof());
    elemMat.Init();
    

    // Loop over all integration points
    LocPointMapped lp;
    const UInt numIntPts = intPoints.GetSize();
    for( UInt i = 0; i < numIntPts; i++  ) {

      // Calculate for each integration point the LocPointMapped
      lp.Set( intPoints[i], esm, weights[i] );

      // Call the CalcBMat()-method
      bOperator_->CalcOpMat( bMat_, lp, ptFe);

      // LOG_DBG3(bdbint) << "CEM e1=" << ptElem->elemNum << " i=" << i << " bMat=" << bMat_.ToString();


      // Calculate D-Mat
      dData_->GetTensor(dMat_,lp);
      
      fac = MAT_DATA_TYPE(lp.jacDet * weights[i]);

      dbMat_.Resize(dMat_.GetNumRows(),nrFncs * bOperator_->GetDimDof());

#ifdef NDEBUG
      dMat_.Mult_Blas(bMat_,dbMat_,false,false,1.0,0);
      bMat_.Mult_Blas(dbMat_,elemMat,true,false,factor_*fac,1.0);
#else
      dbMat_ = (dMat_ * bMat_) * fac;
      elemMat += Transpose(bMat_) * dbMat_ * factor_;
#endif

    }

    // LOG_DBG3(bdbint) << "CEM e1=" << ptElem->elemNum << " dbMat=" << dbMat_.ToString();

  }

  // ===============
  //  Apply ElemMat
  // ===============
  // 1) General template for double-case -> not implemented
  template< class COEF_DATA_TYPE, class B_DATA_TYPE>
  void BDBInt<COEF_DATA_TYPE,B_DATA_TYPE>::
  ApplyElemMat( Vector<Double>&ret, 
                const Vector<Double>& sol,
                EntityIterator& ent1,
                EntityIterator& ent2 ) {
    EXCEPTION("Not implemented")
  }  
  
  // 2) Special double-only case
  template<>
  void BDBInt<Double, Double>::
  ApplyElemMat( Vector<Double>&ret, 
                const Vector<Double>& sol,
                EntityIterator& ent1,
                EntityIterator& ent2 ) {
    Matrix<MAT_DATA_TYPE> elemMat;
    CalcElementMatrix(elemMat, ent1, ent2);
    ret = elemMat * sol;
  }  
  
  // 3) Complex case is always possible
  template< class COEF_DATA_TYPE, class B_DATA_TYPE>
  void BDBInt<COEF_DATA_TYPE,B_DATA_TYPE>::
  ApplyElemMat( Vector<Complex>&ret, 
                const Vector<Complex>& sol,
                EntityIterator& ent1,
                EntityIterator& ent2 ) {
    Matrix<MAT_DATA_TYPE> elemMat;
    CalcElementMatrix(elemMat, ent1, ent2);
    ret = elemMat * sol;
  }  
  
  // ===============
  //  Apply BMatrix
  // ===============
  // 1) General template for double-case -> not implemented
  template< class COEF_DATA_TYPE, class B_DATA_TYPE>
  void BDBInt<COEF_DATA_TYPE,B_DATA_TYPE>::
  ApplyBMat( Vector<Double>&ret, 
                        const Vector<Double>& sol,
                        const LocPointMapped& lpm ) {
    EXCEPTION("Not implemented");
  }

  // 2) Special double-only case
  template<>
  void BDBInt<Double, Double>::
  ApplyBMat( Vector<Double>&ret, 
             const Vector<Double>& sol,
             const LocPointMapped& lpm ) {
    BaseFE* ptFe = ptFeSpace1_->GetFe( lpm.ptEl->elemNum );
    bOperator_->CalcOpMat(bMat_, lpm, ptFe);
    ret = bMat_ * sol;
  }
  
  // 3) Complex case is always possible
  template< class COEF_DATA_TYPE, class B_DATA_TYPE>
  void BDBInt<COEF_DATA_TYPE,B_DATA_TYPE>::
  ApplyBMat( Vector<Complex>&ret, 
             const Vector<Complex>& sol,
             const LocPointMapped& lpm ) {
    BaseFE* ptFe = ptFeSpace1_->GetFe( lpm.ptEl->elemNum );
    bOperator_->CalcOpMat(bMat_, lpm, ptFe);
    ret = bMat_ * sol;
  }
  

  // ===============
  //  Apply dBMatrix
  // ===============
  // 1) General template for double-case -> not implemented
  template< class COEF_DATA_TYPE, class B_DATA_TYPE>
  void BDBInt<COEF_DATA_TYPE,B_DATA_TYPE>::
  ApplydBMat( Vector<Double>&ret, 
              const Vector<Double>& sol,
              const LocPointMapped& lpm ) {
    EXCEPTION("Not implemented")
  }
  
  // 2) Special double-only case
  template<>
  void BDBInt<Double, Double>::
  ApplydBMat( Vector<Double>&ret, 
              const Vector<Double>& sol,
              const LocPointMapped& lpm ) {
    BaseFE* ptFe = ptFeSpace1_->GetFe( lpm.ptEl->elemNum );
    bOperator_->CalcOpMat(bMat_, lpm, ptFe);
    dData_->GetTensor(dMat_,lpm);
    dbMat_.Resize(dMat_.GetNumRows(), bMat_.GetNumCols());
    dMat_.Mult_Blas(bMat_,dbMat_,false,false,1.0,0);
    ret = dbMat_ * sol;
  }
  
  // 3) Complex case is always possible
  template< class COEF_DATA_TYPE, class B_DATA_TYPE>
  void BDBInt<COEF_DATA_TYPE,B_DATA_TYPE>::
  ApplydBMat( Vector<Complex>&ret, 
              const Vector<Complex>& sol,
              const LocPointMapped& lpm ) {
    BaseFE* ptFe = ptFeSpace1_->GetFe( lpm.ptEl->elemNum );
    bOperator_->CalcOpMat(bMat_, lpm, ptFe);
    dData_->GetTensor(dMat_,lpm);
    dbMat_.Resize(dMat_.GetNumRows(), bMat_.GetNumCols());
    dMat_.Mult_Blas(bMat_,dbMat_,false,false,1.0,0);
    ret = dbMat_ * sol;
  }
  
  
  // ===========================
  //  Apply transposed A Matrix
  // ============================
  // 1) General template for double-case -> not implemented
  template< class COEF_DATA_TYPE, class B_DATA_TYPE>
  void BDBInt<COEF_DATA_TYPE,B_DATA_TYPE>::
  ApplyATransMat( Vector<Double>&ret, 
                  const Vector<Double>& sol,
                  const LocPointMapped& lpm ) {
    EXCEPTION("Not implemented");
  }
  
  // 2) Special double-only case
  template<>
  void BDBInt<Double, Double>::
  ApplyATransMat( Vector<Double>&ret, 
                  const Vector<Double>& sol,
                  const LocPointMapped& lpm ) {
    BaseFE* ptFe = ptFeSpace1_->GetFe( lpm.ptEl->elemNum );
    bOperator_->CalcOpMat(bMat_, lpm, ptFe);
    ret = bMat_ * sol;
  }

  // 3) Complex case is always possible
  template< class COEF_DATA_TYPE, class B_DATA_TYPE>
  void BDBInt<COEF_DATA_TYPE,B_DATA_TYPE>::
  ApplyATransMat( Vector<Complex>&ret, 
                  const Vector<Complex>& sol,
                  const LocPointMapped& lpm ) {
    BaseFE* ptFe = ptFeSpace1_->GetFe( lpm.ptEl->elemNum );
    bOperator_->CalcOpMat(bMat_, lpm, ptFe);
    ret = bMat_ * sol;
  }
  
  // ============================
  //  Apply transposed dA Matrix
  // =============================
  
  template< class COEF_DATA_TYPE, class B_DATA_TYPE>
  void BDBInt<COEF_DATA_TYPE,B_DATA_TYPE>::
  ApplydATransMat( Vector<Double>&ret, 
                   const Vector<Double>& sol,
                   const LocPointMapped& lpm ) {
    EXCEPTION("Not implemented");
  }
  
  // 2) Special double-only case
  template<>
  void BDBInt<Double, Double>::
  ApplydATransMat( Vector<Double>&ret, 
                   const Vector<Double>& sol,
                   const LocPointMapped& lpm ) {
    BaseFE* ptFe = ptFeSpace1_->GetFe( lpm.ptEl->elemNum );
    bOperator_->CalcOpMat(bMat_, lpm, ptFe);
    dData_->GetTensor(dMat_,lpm);
    dbMat_.Resize(dMat_.GetNumRows(), bMat_.GetNumCols());
    dMat_.Mult_Blas(bMat_,dbMat_,false,false,1.0,0);
    ret = dbMat_* sol;
  }
  
  // 3) Complex case is always possible
  template< class COEF_DATA_TYPE, class B_DATA_TYPE>
  void BDBInt<COEF_DATA_TYPE,B_DATA_TYPE>::
  ApplydATransMat( Vector<Complex>&ret, 
                   const Vector<Complex>& sol,
                   const LocPointMapped& lpm ) {
    BaseFE* ptFe = ptFeSpace1_->GetFe( lpm.ptEl->elemNum );
    bOperator_->CalcOpMat(bMat_, lpm, ptFe);
    dData_->GetTensor(dMat_,lpm);
    dbMat_.Resize(dMat_.GetNumRows(), bMat_.GetNumCols());
    dMat_.Mult_Blas(bMat_,dbMat_,false,false,1.0,0);
    ret = dbMat_* sol;
  }
  
  //! Calculate the integration kernel
  template< class COEF_DATA_TYPE, class B_DATA_TYPE>
  void BDBInt<COEF_DATA_TYPE,B_DATA_TYPE>::CalcKernel( Matrix<MAT_DATA_TYPE>& kernel,const LocPointMapped& lpm ) {

    // Obtain FE element from feSpace and integration scheme
    BaseFE* ptFe = ptFeSpace1_->GetFe( lpm.ptEl->elemNum );

    UInt nrFncs = ptFe->GetNumFncs();

    kernel.Resize( nrFncs * bOperator_->GetDimDof());
    kernel.Init();


    // Call the CalcBMat()-method
    bOperator_->CalcOpMat( bMat_, lpm, ptFe);

    // LOG_DBG3(bdbint) << "CK e=" << lpm.ptEl->elemNum << " bMat=" << bMat_.ToString();

    // Calculate D-Mat
    dData_->GetTensor(dMat_,lpm);
    dbMat_.Resize(dMat_.GetNumRows(),nrFncs* bOperator_->GetDimDof());

#ifdef NDEBUG
    dMat_.Mult_Blas(bMat_,dbMat_,false,false,1.0,0);
    bMat_.Mult_Blas(dbMat_,kernel,true,false,factor_,1.0, true); // conjugate complex
#else
    dbMat_ = (dMat_ * bMat_);
    kernel += TransposeConjugate(bMat_) * dbMat_ * factor_;
#endif

  }
  
  template< class COEF_DATA_TYPE,
            class B_DATA_TYPE>
   void BDBInt<COEF_DATA_TYPE,B_DATA_TYPE>::__Instantiate() {
    EXCEPTION("This method may never be called.");
    
    Vector<Double> rVec;
    Vector<Complex> cVec;
    LocPointMapped lpm;
    ApplyBMat(rVec, rVec, lpm);
    ApplyBMat(cVec, cVec, lpm);
  }
  
  
// Explicit template instantiation
  template class BDBInt<Double,Double>;
  template class BDBInt<Double,Complex>;
  template class BDBInt<Complex,Double>;
  template class BDBInt<Complex,Complex>;
  
} // namespace
