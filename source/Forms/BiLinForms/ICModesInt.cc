
#include "ICModesInt.hh"
#include "FeBasis/H1/H1Elems.hh"

namespace CoupledField{


template< class COEF_DATA_TYPE, class B_DATA_TYPE>
ICModesInt<COEF_DATA_TYPE, B_DATA_TYPE>::
ICModesInt(BaseBOperator * bOp, 
           BaseBOperator * gOp,
           PtrCoefFct dData, MAT_DATA_TYPE factor,
           bool coordUpdate )
           : BDBInt<COEF_DATA_TYPE, B_DATA_TYPE>(bOp, dData, factor, coordUpdate ) {
  
  this->type_ = BiLinearForm::IC_MODES_INT;
  this->name_ = "ICMOdesInt";
  this->gOperator_ = gOp;
  this->isSymmetric_ = true;
}

template< class COEF_DATA_TYPE, class B_DATA_TYPE>
void ICModesInt<COEF_DATA_TYPE, B_DATA_TYPE>::
CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                   EntityIterator& ent1,
                   EntityIterator& ent2 ) {
  

  // Extract physical element
  const Elem* ptElem = ent1.GetElem();
  
  MAT_DATA_TYPE fac(0.0);

  // Obtain FE element from feSpace and integration scheme
  IntegOrder order;
  IntScheme::IntegMethod method;
  FeH1 * ptFe = static_cast<FeH1*>(this->ptFeSpace1_->GetFe( ent1, method, order ));

  const UInt nrFncs = ptFe->GetNumFncs();
  const UInt nrIcModes = ptFe->GetNumICModes();
  const UInt nrDofs = this->bOperator_->GetDimDof();

  // ensure, that incompatible modes are present at all
  if( !nrIcModes ) {
    EXCEPTION( "Element does not seem to support incompatible modes softening!\n"
        << "Please check that only 1st order quadrilateral and hexahedral "
        << "elements are present!");
  }
  // Get shape map from grid
  // shared_ptr<ElemShapeMap> esm =
  //     ent1.GetGrid()->GetElemShapeMap( ptElem, this->coordUpdate_ );
  // shared_ptr<ElemShapeMap> esm(ptElem->ptrShapeMap);
  shared_ptr<ElemShapeMap> esm = (ptElem)->GetElemShapeMap(ent1.GetGrid(), this->coordUpdate_);

  // Get integration points
  StdVector<LocPoint> intPoints;
  StdVector<Double> weights;
  this->intScheme_->GetIntPoints( Elem::GetShapeType(ptElem->type), method, order, 
                            intPoints, weights );

  elemMat.Resize(nrFncs * nrDofs);
  elemMat.Init();
  
  // helper matrices
  Matrix<MAT_DATA_TYPE> gMat, dgMat;
  Matrix<MAT_DATA_TYPE> elemMat12; // mixed stiffness matrix
  Matrix<MAT_DATA_TYPE> elemMat22; // incompatible modes matrix
  elemMat12.Resize(nrIcModes * nrDofs, nrFncs * nrDofs);
  elemMat12.Init();
  elemMat22.Resize(nrIcModes * nrDofs);
  elemMat22.Init();
  
  
  // Loop over all integration points
  LocPointMapped lp;
  const UInt numIntPts = intPoints.GetSize();
  for( UInt i = 0; i < numIntPts; i++  ) {

    // Calculate for each integration point the LocPointMapped
    lp.Set( intPoints[i], esm, weights[i] );

    // Call the CalcBMat()-method
    this->bOperator_->CalcOpMat( this->bMat_, lp, ptFe);
    gOperator_->CalcOpMat( gMat_, lp, ptFe);

    // Calculate D-Mat
    this->dData_->GetTensor(this->dMat_,lp);

    fac = MAT_DATA_TYPE(lp.jacDet * weights[i]);

    // Note: 
    // The following section seems obsolete at least for the
    // plane case, as in the trunk version, the G-Matrix was
    // already scaled by * jacDet, so it cancels out.
    // However, for the axisymmetric case this is not verified
    // yet.
    
    //    // Scale g-matrix appropriately
    //    Double jacDetInv = 1.0 / lp.jacDet;
    //    if ( lp.shapeMap->IsAxi() ) {
    //      for ( UInt i=0; i<gMat_.GetNumCols(); i++)
    //        for ( UInt j=0; j<gMat_.GetNumRows()-1; j++ )
    //          gMat_[j][i] *=  jacDetInv;
    //    }
    //    else
    //      gMat_ *=  jacDetInv;


    this->dbMat_.Resize(this->dMat_.GetNumRows(), nrDofs * nrFncs);
    this->dgMat_.Resize(this->dMat_.GetNumRows(), nrDofs * nrIcModes);

#ifdef NDEBUG
    this->dMat_.Mult_Blas(this->bMat_, this->dbMat_, false, false, 1.0, 0);
    this->bMat_.Mult_Blas(this->dbMat_,elemMat,true,false,
                          this->factor_ * fac,1.0);
    
    this->dMat_.Mult_Blas(this->gMat_, dgMat_, false, false, 1.0, 0.0);
    this->gMat_.Mult_Blas(this->dgMat_,elemMat22, true, false, 
                          this->factor_ * fac, 1.0);

    this->gMat_.Mult_Blas(this->dbMat_, elemMat12, true, false, 
                         this->factor_ * fac, 1.0);
    
#else 
    // We now compute B^T * D * B and scale it by the determinant
    // of the Jacobian and the weight of the current integration
    // point. The result is added to the element matrix.
    this->dbMat_ = (this->dMat_ * this->bMat_) * fac;
    elemMat += Transpose(this->bMat_) * this->dbMat_ * this->factor_;

    // We now compute G^T * D * G and scale it by the determinant
    // of the Jacobian and the weight of the current integration
    // point. The result is added to the element matrix.
    dgMat_ = (this->dMat_ * gMat_) * fac;
    elemMat22 += Transpose(this->gMat_) * this->dgMat_ * this->factor_;

    // We now compute G^T * D * B and scale it by the determinant
    // of the Jacobian and the weight of the current integration
    // point. The result is added to the element matrix.
    elemMat12 += Transpose(this->gMat_) * this->dbMat_ * this->factor_;
#endif
    
  }
  
  // Invert K22 matrix
  Matrix<MAT_DATA_TYPE> invElemMat22;
  elemMat22.Invert(invElemMat22);

  Matrix<MAT_DATA_TYPE> part1 (nrIcModes * nrDofs, nrFncs * nrDofs);
#ifdef NDEBUG
  // We now compute elemMat - k12^T * k22^(-1) * k12
  invElemMat22.Mult_Blas(elemMat12, part1, false, false, 1.0, 0);
  elemMat12.Mult_Blas(part1, elemMat, true, false, -1.0, 1.0);
#else
  part1 = invElemMat22 * elemMat12;
  elemMat -= Transpose(elemMat12)*part1;
#endif
}

  // Explicit template instantiation
  template class ICModesInt<Double,Double>;
  template class ICModesInt<Complex,Double>;

} // namespace
