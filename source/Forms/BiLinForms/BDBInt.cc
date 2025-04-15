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
#include "Domain/Domain.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Domain/CoefFunction/CoefFunctionOpt.hh"

DEFINE_LOG(bdbint, "bdbint")

namespace CoupledField{

  template< class COEF_DATA_TYPE,
            class B_DATA_TYPE>
  BDBInt<COEF_DATA_TYPE,B_DATA_TYPE>::
  BDBInt(BaseBOperator* bOp, PtrCoefFct dData, MAT_DATA_TYPE factor, bool coordUpdate )
  : BaseBDBInt(coordUpdate) 
    {
      this->type_ = BDB_INT;
      this->name_ = type.ToString(type_);
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
      isSolDependent_ = dData->GetDependency() == CoefFunction::SOLUTION;
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
  void BDBInt<COEF_DATA_TYPE,B_DATA_TYPE>::CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat, EntityIterator& ent1, EntityIterator& ent2)
  {
    // Extract physical element
    const Elem* ptElem = ent1.GetElem();

    // for optimization we cache the local element matrices, e.g. to speed up assembly and apply
    // ersatz material design parametrization directly on the elements
    if(domain->HasDesign())
      if(domain->GetDesign()->ApplyPhysicalDesignElementMatrix<MAT_DATA_TYPE>(this, elemMat, ptElem))
        return;

    MAT_DATA_TYPE fac(0.0);

    // Obtain FE element from feSpace and integration scheme
    IntegOrder order;
    IntScheme::IntegMethod method;
    BaseFE* ptFe = ptFeSpace1_->GetFe(ent1, method, order );

    const UInt nrFncs = ptFe->GetNumFncs();

    // Get shape map from grid
    shared_ptr<ElemShapeMap> esm = ent1.GetGrid()->GetElemShapeMap( ptElem, this->coordUpdate_ );
    // Get integration points
    StdVector<LocPoint> intPoints;
    StdVector<Double> weights;
    intScheme_->GetIntPoints( Elem::GetShapeType(ptElem->type), method, order, intPoints, weights );

    elemMat.Resize(nrFncs * bOperator_->GetDimDof());
    elemMat.Init();
    
    // Loop over all integration points
    LocPointMapped lp;
    // if MSFEM get Element Matrix from material catalog
    if(domain->HasDesign() && domain->GetDesign()->DoMSFEM())
    {
      lp.Set( intPoints[0], esm, weights[0]);
      // it would be nicer to have a separate MS-FEM call
      if(domain->GetDesign()->ApplyPhysicalDesign<MAT_DATA_TYPE>(dynamic_cast<const CoefFunctionOpt*>(dData_.get()), elemMat, &lp))
      {
        //LOG_DBG3(bdbint) << "BDBInt: MS-FEM elem=" << ptElem->elemNum << " == "<< lp.ptEl->elemNum << " elemMat=" << elemMat.ToString();
        return;
      }
    }
    const UInt numIntPts = intPoints.GetSize();
    for( UInt i = 0; i < numIntPts; i++  )
    {
      // Calculate for each integration point the LocPointMapped
      lp.Set( intPoints[i], esm, weights[i] );

      // Call the CalcBMat()-method
      bOperator_->CalcOpMat( bMat_, lp, ptFe);

      LOG_DBG3(bdbint) << "CEM e=" << ptElem->elemNum << " ip=" << i << " bMat=" << bMat_.ToString();

      // Calculate D-Mat
      dData_->GetTensor(dMat_,lp);
      assert(dMat_.IsSymmetric(1e-8));
      LOG_DBG3(bdbint) << "CEM e=" << ptElem->elemNum << " ip=" << i << " dMat=" << dMat_.ToString();

      fac = MAT_DATA_TYPE(lp.jacDet * weights[i]);

      LOG_DBG3(bdbint) << "CEM e=" << ptElem->elemNum << " ip=" << i << " factor_=" << factor_ << " jacDet=" << lp.jacDet << " weight=" << weights[i];

      dbMat_.Resize(dMat_.GetNumRows(),nrFncs * bOperator_->GetDimDof());

#ifdef NDEBUG
      dMat_.Mult_Blas(bMat_,dbMat_,false,false,1.0,0); // dbMat_ = 1.0 * dMat_ * bMat_ + 0.0 * dbMat_
      bMat_.Mult_Blas(dbMat_,elemMat,true,false,factor_*fac,1.0, true); // conjugate complex; elemMat = factor_*fac * bMat_^H * dbMat_ + 1.0 * elemMat
#else
      dbMat_ = (dMat_ * bMat_) * fac;
      elemMat += TransposeConjugate(bMat_) * dbMat_ * factor_;

      //LOG_DBG3(bdbint) << "bMat: " << bMat_.ToString() << " transpose:"  << TransposeConjugate(bMat_).ToString();
      LOG_DBG3(bdbint) << "CEM e=" << ptElem->elemNum << " ip=" << i << " fac=" << fac << " factor_=" << factor_ << " dMat= " << dMat_.ToString() << " bmat=" << bMat_.ToString();
      //LOG_DBG3(bdbint) << "CEM e=" << ptElem->elemNum << " dBMat=" << dbMat_.ToString() << " -> K_" << i << "=" << elemMat.ToString();
#endif

    }

    //LOG_DBG3(bdbint) << "CEM e=" << ptElem->elemNum << " elemMat=" << elemMat.ToString();
  }


  template< class COEF_DATA_TYPE, class B_DATA_TYPE> 
  void BDBInt<COEF_DATA_TYPE, B_DATA_TYPE>::
  CalcElementMatrixLpm( Matrix<MAT_DATA_TYPE>& elemMat,
                      BaseFE* ptFe, 
                      const LocPointMapped& lp, 
                      bool overrideIsSurfOpt ) {
    
    MAT_DATA_TYPE fac = 0.0;

    const UInt nrFncs = ptFe->GetNumFncs();

    elemMat.Resize( nrFncs * bOperator_->GetDimDof() );
    elemMat.Init();

    // if we are dealing with network coupling, we might evaluate a surface within an
    // volume element integrator
    // therefore, we have the ability to manually override the integrator and tell it
    // to be a surface integrator
    
    this->bOperator_->OverrideIsSurfOperator(overrideIsSurfOpt);

    // Call the CalcBMat()-method
    this->bOperator_->CalcOpMat( bMat_, lp, ptFe);
    LOG_DBG3(bdbint) << "point " << lp.lp.coord << " bMat= " << bMat_;

    // reset it
    this->bOperator_->OverrideIsSurfOperator(false);

    // Calculate D-Mat
    dData_->GetTensor(dMat_,lp);
    assert(dMat_.IsSymmetric(1e-8));
    LOG_DBG3(bdbint) << "point " << lp.lp.coord << " dMat=" << dMat_.ToString();
    
    fac = MAT_DATA_TYPE(lp.jacDet *lp.weight); 

    LOG_DBG3(bdbint) << "point " << lp.lp.coord << " factor_=" << factor_ << " jacDet=" << lp.jacDet << " weight=" << lp.weight;

    dbMat_.Resize(dMat_.GetNumRows(),nrFncs * bOperator_->GetDimDof());

#ifdef NDEBUG
    dMat_.Mult_Blas(bMat_,dbMat_,false,false,1.0,0); // dbMat_ = 1.0 * dMat_ * bMat_ + 0.0 * dbMat_
    bMat_.Mult_Blas(dbMat_,elemMat,true,false,factor_*fac,1.0, true); // conjugate complex; elemMat = factor_*fac * bMat_^H * dbMat_ + 1.0 * elemMat
#else
    dbMat_ = (dMat_ * bMat_) * fac;
    elemMat += TransposeConjugate(bMat_) * dbMat_ * factor_;
    LOG_DBG3(bdbint) << "CEM point " << lp.lp.coord << " fac=" << fac << " factor_=" << factor_ << " bmat=" << bMat_.ToString();
    LOG_DBG3(bdbint) << "CEM point " << lp.lp.coord << " -> K_" << "=" << elemMat.ToString();
#endif
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

    // this calculates K = B^* C B with B^* is conjugate complex. The effect is only seen form Bloch mode analysis
    // where we have a complex B matrix.

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
