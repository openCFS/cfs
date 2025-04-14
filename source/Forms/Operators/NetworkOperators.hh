// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     NetworkOperators.hh
 *       \brief    Collection of Operators used for network coupling
 *
 *       \date     Jan 7, 2015
 *       \author   dmayrhof
 */
//================================================================================================

#ifndef NETWORKOPERATORS_HH
#define NETWORKOPERATORS_HH

#include "BaseBOperator.hh"

namespace CoupledField{


template<class FE, UInt D = 1, UInt D_DOF = 1, class TYPE = Double>
class SurfaceMagneticStiffnessMatrixH1 : public BaseBOperator{

public:

  // ------------------
  //  STATIC CONSTANTS
  // ------------------
  //@{
  //! \name Static constants

  //! Order of differentiation
  static const UInt ORDER_DIFF = 1;

  //! Number of components of the problem (scalar, vector)
  static const UInt DIM_DOF = D_DOF;

  //! Dimension of the underlying domain / space
  static const UInt DIM_SPACE = D;

  //! Dimension of the finite element
  static const UInt DIM_ELEM = D;

  //! Dimension of the related material
  static const UInt DIM_D_MAT = 1;
  //@}

  
  SurfaceMagneticStiffnessMatrixH1(){
    return;
  }

  SurfaceMagneticStiffnessMatrixH1(const SurfaceMagneticStiffnessMatrixH1 & other)
    : BaseBOperator(other){
  }

  virtual SurfaceMagneticStiffnessMatrixH1 * Clone(){
    return new SurfaceMagneticStiffnessMatrixH1(*this);
  }

  virtual ~SurfaceMagneticStiffnessMatrixH1(){
    return;
  }

  virtual void CalcOpMat(Matrix<Double> & bMat,
                          const LocPointMapped& lp, BaseFE* ptFe );

  virtual void CalcOpMatTransposed(Matrix<Double> & bMat,
                                    const LocPointMapped& lp, BaseFE* ptFe );

  //avoid reimplementation of complex operator by making the bas class function
  //available
  using BaseBOperator::CalcOpMat;

  using BaseBOperator::CalcOpMatTransposed;

  // ===============
  //  QUERY METHODS
  // ===============
  //@{ \name Query Methods
  //! \copydoc BaseBOperator::GetDiffOrder
  virtual UInt GetDiffOrder() const {
    return ORDER_DIFF;
  }

  //! \copydoc BaseBOperator::GetDimDof()
  virtual UInt GetDimDof() const {
    return DIM_DOF;
  }

  //! \copydoc BaseBOperator::GetDimSpace()
  virtual UInt GetDimSpace() const {
    return DIM_SPACE;
  }

  //! \copydoc BaseBOperator::GetDimElem()
  virtual UInt GetDimElem() const {
    return DIM_ELEM;
  }

  //! \copydoc BaseBOperator::GetDimDMat()
  virtual UInt GetDimDMat() const {
    return DIM_D_MAT;
  }
  //@}


protected:

};

template<class FE, UInt D, UInt D_DOF, class TYPE>
void SurfaceMagneticStiffnessMatrixH1<FE, D, D_DOF, TYPE>::CalcOpMat(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
{
  //check if lp is surface and ptFe is volume
  assert(lp.isSurface);
  assert(D == ptFe->shape_.dim);
  assert(D==D_DOF);
  assert(D==2);

  const UInt numFncs = ptFe->GetNumFncs();
      
  // Set correct size of matrix B and initialize with zeros
  Matrix<Double>& bMatIntermed;
  bMatIntermed.Resize( 2, numFncs );

  // Get derivatives of local shape functions with respect to global
  // coords (format: nrNodes x spaceDim)
  Matrix<Double> xiDx;
  FeH1 *fe = (static_cast<FeH1*>(ptFe));
  if (this->isSurfOpt_)
    fe->GetGlobDerivShFnc(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
  else
    fe->GetGlobDerivShFnc( xiDx, lp, lp.shapeMap->GetElem() , 1 );

  UInt iFunc = 0;
  for( iFunc = 0; iFunc < numFncs; ++iFunc ) {
    bMatIntermed[0][iFunc] =  xiDx[iFunc][1];
  }
  for( iFunc = 0; iFunc < numFncs; ++iFunc ) {
    bMatIntermed[1][iFunc] = -xiDx[iFunc][0];
  }

  EXCEPTION("SurfaceMagneticStiffnessMatrixH1 not fully implmeneted yet");
  //bMat = bMatIntermed.Transpose * this->coef_ * bMatIntermed
}

template<class FE, UInt D, UInt D_DOF, class TYPE>
void SurfaceMagneticStiffnessMatrixH1<FE, D, D_DOF, TYPE>::CalcOpMatTransposed(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
{
  Matrix<Double> dummyMat;
  CalcOpMat(dummyMat, lp, ptFe);

  dummyMat.Transpose(bMat);
}






template<class FE, UInt D = 1, UInt D_DOF = 1, class TYPE = Double>
class FemLemAllocationOperator : public BaseBOperator{

public:

  // ------------------
  //  STATIC CONSTANTS
  // ------------------
  //@{
  //! \name Static constants

  //! Order of differentiation
  static const UInt ORDER_DIFF = 1;

  //! Number of components of the problem (scalar, vector)
  static const UInt DIM_DOF = D_DOF;

  //! Dimension of the underlying domain / space
  static const UInt DIM_SPACE = D;

  //! Dimension of the finite element
  static const UInt DIM_ELEM = D;

  //! Dimension of the related material
  static const UInt DIM_D_MAT = 1;
  //@}

  
  FemLemAllocationOperator(){
    return;
  }

  FemLemAllocationOperator(const FemLemAllocationOperator & other)
    : BaseBOperator(other){
  }

  virtual FemLemAllocationOperator * Clone(){
    return new FemLemAllocationOperator(*this);
  }

  virtual ~FemLemAllocationOperator(){
    return;
  }

  virtual void CalcOpMat(Matrix<Double> & bMat,
                          const LocPointMapped& lp, BaseFE* ptFe );

  virtual void CalcOpMatTransposed(Matrix<Double> & bMat,
                                    const LocPointMapped& lp, BaseFE* ptFe );

  //avoid reimplementation of complex operator by making the bas class function
  //available
  using BaseBOperator::CalcOpMat;

  using BaseBOperator::CalcOpMatTransposed;

  // ===============
  //  QUERY METHODS
  // ===============
  //@{ \name Query Methods
  //! \copydoc BaseBOperator::GetDiffOrder
  virtual UInt GetDiffOrder() const {
    return ORDER_DIFF;
  }

  //! \copydoc BaseBOperator::GetDimDof()
  virtual UInt GetDimDof() const {
    return DIM_DOF;
  }

  //! \copydoc BaseBOperator::GetDimSpace()
  virtual UInt GetDimSpace() const {
    return DIM_SPACE;
  }

  //! \copydoc BaseBOperator::GetDimElem()
  virtual UInt GetDimElem() const {
    return 0; // this one is only used for points
  }

  //! \copydoc BaseBOperator::GetDimDMat()
  virtual UInt GetDimDMat() const {
    return DIM_D_MAT;
  }
  //@}


protected:

};

template<class FE, UInt D, UInt D_DOF, class TYPE>
void FemLemAllocationOperator<FE, D, D_DOF, TYPE>::CalcOpMat(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
{
  //check if lp is surface and ptFe is volume
  //assert(lp.isSurface); TODO
  //assert(D == ptFe->shape_.dim);
  //assert(D==D_DOF);
  //assert(D==2); TODO

  const UInt numFncs = ptFe->GetNumFncs();
      
  // Set correct size of matrix B and initialize with zeros
  bMat.Resize( 1, numFncs );

  // the allocation vector only consists of ones since we only consier elements that actually have a contribution
  // the sign is governed by the integrator based on the network element definition.
  UInt iFunc = 0;
  for( iFunc = 0; iFunc < numFncs; ++iFunc ) {
    bMat[0][iFunc] =  1;
  }
}

template<class FE, UInt D, UInt D_DOF, class TYPE>
void FemLemAllocationOperator<FE, D, D_DOF, TYPE>::CalcOpMatTransposed(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
{
  Matrix<Double> dummyMat;
  CalcOpMat(dummyMat, lp, ptFe);

  dummyMat.Transpose(bMat);
}



}
#endif
