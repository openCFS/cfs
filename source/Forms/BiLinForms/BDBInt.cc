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


namespace CoupledField{

  template< class B_OP,
            class MAT_DATA_TYPE,
            class COEF_DATA_TYPE>
  BDBInt<B_OP,MAT_DATA_TYPE,COEF_DATA_TYPE>::
  BDBInt(shared_ptr<CoefFunction> dData, MAT_DATA_TYPE factor, bool coordUpdate )
  : BaseBDBInt(coordUpdate) 
  {
      name_ = "BDBInt";
      isSymmetric_ = true;
      assert(dData->GetDimType() == CoefFunction::TENSOR);
#ifndef NDEBUG
      if(dData->GetDimType() != CoefFunction::TENSOR){
        Exception("BDB integrator expects the coefficient function to be tensorial");
      }
#endif
      dData_ = dData;
      factor_ = factor;
  }

  //! Destructor
  template< class B_OP,
              class MAT_DATA_TYPE,
              class COEF_DATA_TYPE>
  BDBInt<B_OP,MAT_DATA_TYPE,COEF_DATA_TYPE>::~BDBInt(){

  }

  //! Calculate the element matrix
  template< class B_OP,
              class MAT_DATA_TYPE,
              class COEF_DATA_TYPE>
  void BDBInt<B_OP,MAT_DATA_TYPE,COEF_DATA_TYPE>::
  CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                     EntityIterator& ent1,
                     EntityIterator& ent2) {

    // Extract physical element
    const Elem* ptElem = ent1.GetElem();

    Matrix<MAT_DATA_TYPE> bMat, dbMat;
    Matrix<COEF_DATA_TYPE> dMat;
    MAT_DATA_TYPE fac = 0.0;

    // Obtain FE element from feSpace and integration scheme
    IntegOrder order;
    IntScheme::IntegMethod method;
    BaseFE* ptFe = ptFeSpace1_->GetFe( ent1, method, order );

    
    UInt nrFncs = ptFe->GetNumFncs();

    // Get shape map from grid
    shared_ptr<ElemShapeMap> esm = 
        domain->GetGrid()->GetElemShapeMap( ptElem, this->coordUpdate_ );

    // Get integration points
    StdVector<LocPoint> intPoints;
    StdVector<Double> weights;
    intScheme_->GetIntPoints( Elem::GetShapeType(ptElem->type), method, order, 
                              intPoints, weights );

    elemMat.Resize( nrFncs * B_OP::DIM_DOF);
    elemMat.Init();
    
#define USE_BLAS_VERSION

    // Loop over all integration points
    LocPointMapped lp;
    for( UInt i = 0; i < intPoints.GetSize(); i++  ) {

      // Calculate for each integration point the LocPointMapped
      lp.Set( intPoints[i], esm );

      // Call the CalcBMat()-method
      bOperator_.CalcOpMat( bMat, lp, ptFe);

      // Calculate D-Mat
      dData_->GetTensor(dMat,lp);

      fac = MAT_DATA_TYPE(lp.jacDet * weights[i]);
      bOperator_.TransformJacDet(fac,lp,ptFe);

      dbMat.Resize(dMat.GetNumRows(),nrFncs * B_OP::DIM_DOF);

#ifdef USE_BLAS_VERSION
      dMat.Mult_Blas(bMat,dbMat,false,false,1.0,0);
      bMat.Mult_Blas(dbMat,elemMat,true,false,factor_*fac,1.0);
#else
      dbMat = (dMat * bMat) * fac;
      elemMat += Transpose(bMat) * dbMat * factor_;
#endif

    }

  }

  //! Multiply element matrix with vector
  template< class B_OP,
  class MAT_DATA_TYPE,
  class COEF_DATA_TYPE>
  void BDBInt<B_OP,MAT_DATA_TYPE,COEF_DATA_TYPE>::
  ApplyElemMat( Vector<MAT_DATA_TYPE>&ret, 
                const Vector<Double>& sol,
                EntityIterator& ent1,
                EntityIterator& ent2 ) {
    Matrix<MAT_DATA_TYPE> elemMat;
    CalcElementMatrix(elemMat, ent1, ent2);
    ret = elemMat * sol;
  }

  //! Apply B-operator on vector
  template< class B_OP,
  class MAT_DATA_TYPE,
  class COEF_DATA_TYPE>
  void BDBInt<B_OP,MAT_DATA_TYPE,COEF_DATA_TYPE>::
  ApplyBMat( Vector<MAT_DATA_TYPE>&ret, 
             const Vector<MAT_DATA_TYPE>& sol,
             const LocPointMapped& lpm ) {
    Matrix<MAT_DATA_TYPE> bOp;
    BaseFE* ptFe = ptFeSpace1_->GetFe( lpm.ptEl->elemNum );
    bOperator_.CalcOpMat(bOp, lpm, ptFe);
    ret = bOp * sol;
  }

  //! Apply dB-operator on vector
  //template<class VEC_DATA_TYPE> 
  template< class B_OP,
  class MAT_DATA_TYPE,
  class COEF_DATA_TYPE>
  void BDBInt<B_OP,MAT_DATA_TYPE,COEF_DATA_TYPE>::
  ApplydBMat( Vector<MAT_DATA_TYPE>&ret, 
              const Vector<MAT_DATA_TYPE>& sol,
              const LocPointMapped& lpm ) {
    Matrix<MAT_DATA_TYPE> bOp, dMat, dbMat;
    BaseFE* ptFe = ptFeSpace1_->GetFe( lpm.ptEl->elemNum );
    bOperator_.CalcOpMat(bOp, lpm, ptFe);
    dData_->GetTensor(dMat,lpm);
    dbMat.Resize(dMat.GetNumRows(), bOp.GetNumCols());
    dMat.Mult_Blas(bOp,dbMat,false,false,1.0,0);
    ret = dbMat* sol;
  }

  //! Calculate the integration kernel
  template< class B_OP,
            class MAT_DATA_TYPE,
            class COEF_DATA_TYPE>
  void BDBInt<B_OP,MAT_DATA_TYPE,COEF_DATA_TYPE>::
  CalcKernel( Matrix<MAT_DATA_TYPE>& kernel, 
              const LocPointMapped& lpm ) {

    Matrix<MAT_DATA_TYPE> bMat, dbMat;
    Matrix<COEF_DATA_TYPE> dMat;
    MAT_DATA_TYPE fac = 0.0;

    // Obtain FE element from feSpace and integration scheme
    BaseFE* ptFe = ptFeSpace1_->GetFe( lpm.ptEl->elemNum );

    UInt nrFncs = ptFe->GetNumFncs();

    kernel.Resize( nrFncs * B_OP::DIM_DOF);
    kernel.Init();

#define USE_BLAS_VERSION

    // Call the CalcBMat()-method
    bOperator_.CalcOpMat( bMat, lpm, ptFe);

    // Calculate D-Mat
    dData_->GetTensor(dMat,lpm);
    bOperator_.TransformJacDet(fac,lpm,ptFe);
    dbMat.Resize(dMat.GetNumRows(),nrFncs* B_OP::DIM_DOF);

#ifdef USE_BLAS_VERSION
    dMat.Mult_Blas(bMat,dbMat,false,false,1.0,0);
    bMat.Mult_Blas(dbMat,kernel,true,false,factor_,1.0);
#else
    dbMat = (dMat * bMat) * fac;
    kernel += Transpose(bMat) * dbMat * factor_;
#endif

  }
  
  template< class B_OP,
             class MAT_DATA_TYPE,
             class COEF_DATA_TYPE>
   void BDBInt<B_OP,MAT_DATA_TYPE,COEF_DATA_TYPE>::__Instantiate() {
    EXCEPTION("This method may never be called.");
    
    Vector<Double> rVec;
    Vector<Complex> cVec;
    LocPointMapped lpm;
    ApplyBMat(rVec, rVec, lpm);
    ApplyBMat(cVec, cVec, lpm);
  }
  
} // namespace
