namespace CoupledField{

 template< class A_OP,
           class B_OP,
           class MAT_DATA_TYPE,
           class COEF_DATA_TYPE >
  ADBInt<A_OP, B_OP, MAT_DATA_TYPE, COEF_DATA_TYPE>::
  ADBInt(shared_ptr<CoefFunction> dData, 
         MAT_DATA_TYPE factor,
         bool coordUpdate )
  : BDBInt<B_OP, MAT_DATA_TYPE, COEF_DATA_TYPE>( dData, factor, coordUpdate  ) 
  {
      this->name_ = "ADBInt";
      this->isSymmetric_ = false; // in general the ADB
      assert(dData->GetDimType() == CoefFunction::TENSOR);
#ifndef NDEBUG
      if(dData->GetDimType() != CoefFunction::TENSOR){
        Exception("BDB integrator expects the coefficient function to be tensorial");
      }
#endif
      this->factor_ = factor;
      
      WARN("We need a clean way to determine the highest integration order if two "
          << "spaces are involved" ); 
  }

  //! Destructor
 template< class A_OP,
           class B_OP,
           class MAT_DATA_TYPE,
           class COEF_DATA_TYPE >
  ADBInt<A_OP, B_OP, MAT_DATA_TYPE, COEF_DATA_TYPE>::
  ~ADBInt(){

  }

  //! Calculate the element matrix
 template< class A_OP,
           class B_OP,
           class MAT_DATA_TYPE,
           class COEF_DATA_TYPE >
  void ADBInt<A_OP, B_OP, MAT_DATA_TYPE, COEF_DATA_TYPE>::
  CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                     EntityIterator& ent1,
                     EntityIterator& ent2) {

    // Extract physical element
    const Elem* ptElem = ent1.GetElem();

    Matrix<MAT_DATA_TYPE> aMat, bMat, dbMat;
    Matrix<COEF_DATA_TYPE> dMat;
    MAT_DATA_TYPE fac = 0.0;

    // Obtain FE element from feSpace and integration scheme
    shared_ptr<IntScheme> intScheme;
    BaseFE* ptFeA = this->ptFeSpace1_->GetFe( ent1, intScheme );
    BaseFE* ptFeB = this->ptFeSpace2_->GetFe( ent1, intScheme );

    UInt nrFncsA = ptFeA->GetNumFncs();
    UInt nrFncsB = ptFeB->GetNumFncs();

    // Get shape map from grid
    shared_ptr<ElemShapeMap> esm = 
        domain->GetGrid()->GetElemShapeMap( ptElem, this->coordUpdate_ );

    // Get integration points
    StdVector<LocPoint> intPoints;
    StdVector<Double> weights;
    intScheme->GetIntPoints( Elem::GetShapeType(ptElem->type), intPoints, weights );

    elemMat.Resize( nrFncsA * A_OP::DIM_DOF, 
                    nrFncsB * B_OP::DIM_DOF );
    elemMat.Init();
    
#define USE_BLAS_VERSION

    // Loop over all integration points
    LocPointMapped lp;
    for( UInt i = 0; i < intPoints.GetSize(); i++  ) {

      // Calculate for each integration point the LocPointMapped
      lp.Set( intPoints[i], esm );

      // Calculate the A-matrix
      aOperator_.CalcOpMat( aMat, lp, ptFeA);
      
      // Calculate the B-matrix
      this->bOperator_.CalcOpMat( bMat, lp, ptFeB);

      // Calculate D-Mat
      this->dData_->GetTensor(dMat,lp);

      fac = MAT_DATA_TYPE(lp.jacDet * weights[i]);
      this->aOperator_.TransformJacDet(fac,lp,ptFeA);
      this->bOperator_.TransformJacDet(fac,lp,ptFeB );

      dbMat.Resize(dMat.GetNumRows(), nrFncsB * B_OP::DIM_DOF);

#ifdef USE_BLAS_VERSION
      dMat.Mult_Blas(bMat,dbMat,false,false,1.0,0);
      aMat.Mult_Blas(dbMat,elemMat,true,false,this->factor_*fac,1.0);
#else
      dbMat = (dMat * bMat) * fac;
      elemMat += Transpose(aMat) * dbMat * this->factor_;
#endif

    }
  }
 
  //! Calculate the integration kernel
 template< class A_OP,
           class B_OP,
           class MAT_DATA_TYPE,
           class COEF_DATA_TYPE >
 void ADBInt<A_OP, B_OP, MAT_DATA_TYPE, COEF_DATA_TYPE>::
  CalcKernel( Matrix<MAT_DATA_TYPE>& kernel, 
              const LocPointMapped& lpm ) {

    Matrix<MAT_DATA_TYPE> aMat, bMat, dbMat;
    Matrix<COEF_DATA_TYPE> dMat;
    MAT_DATA_TYPE fac = 0.0;

    // Obtain FE element from feSpace and integration scheme
    BaseFE* ptFeA = this->ptFeSpace1_->GetFe( lpm.ptEl->elemNum );
    BaseFE* ptFeB = this->ptFeSpace2_->GetFe( lpm.ptEl->elemNum );

    UInt nrFncsA = ptFeA->GetNumFncs();
    UInt nrFncsB = ptFeB->GetNumFncs();

    kernel.Resize( nrFncsA * A_OP::DIM_DOF, 
                   nrFncsB * B_OP::DIM_DOF );
    kernel.Init();

#define USE_BLAS_VERSION

    // Calculate the A-matrix
    aOperator_.CalcOpMat( aMat, lpm, ptFeA);
    
    // Calculate the B-matrix
    this->bOperator_.CalcOpMat( bMat, lpm, ptFeB);

    // Calculate D-Mat
    this->dData_->GetTensor(dMat,lpm);
    this->aOperator_.TransformJacDet(fac,lpm,ptFeA);
    this->bOperator_.TransformJacDet(fac,lpm,ptFeB);
    dbMat.Resize(dMat.GetNumRows(),nrFncsB* B_OP::DIM_DOF);

#ifdef USE_BLAS_VERSION
    dMat.Mult_Blas(bMat,dbMat,false,false,1.0,0);
    aMat.Mult_Blas(dbMat,kernel,true,false,this->factor_,1.0);
#else
    dbMat = (dMat * bMat) * fac;
    kernel += Transpose(aMat) * dbMat * this->factor_;
#endif

  }
  
}
