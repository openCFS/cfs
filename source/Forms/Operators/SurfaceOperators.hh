// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     SurfaceOperators.hh
 *       \brief    Collection of Operators acting on the surfaces of volume elements
 *                 These operators assume always a local point defined on a surface
 *                 and a volume BaseFE. If Not, they will throw an exception
 *
 *       \date     Feb 2, 2012
 *       \author   ahueppe
 */
//================================================================================================

#ifndef SURFACEOPERATORS_HH
#define SURFACEOPERATORS_HH

#include "BaseBOperator.hh"

namespace CoupledField{

template<class FE, UInt D = 1, UInt D_DOF = 1, class TYPE = Double>
class SurfaceIdentityOperator : public BaseBOperator{

public:

   // ------------------
   //  STATIC CONSTANTS
   // ------------------
   //@{
   //! \name Static constants

   //! Order of differentiation
   static const UInt ORDER_DIFF = 0;

   //! Number of components of the problem (scalar, vector)
   static const UInt DIM_DOF = D_DOF;

   //! Dimension of the underlying domain / space
   static const UInt DIM_SPACE = D;

   //! Dimension of the finite element
   static const UInt DIM_ELEM = D;

   //! Dimension of the related material
   static const UInt DIM_D_MAT = D_DOF;
   //@}

   SurfaceIdentityOperator(){
      return;
    }

   SurfaceIdentityOperator(const SurfaceIdentityOperator & other)
    : BaseBOperator(other){
   }

   virtual SurfaceIdentityOperator * Clone(){
     return new SurfaceIdentityOperator(*this);
   }

    virtual ~SurfaceIdentityOperator(){
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

template<class FE,  UInt D, UInt D_DOF, class TYPE>
void SurfaceIdentityOperator<FE,D,D_DOF,TYPE>::
  CalcOpMat(Matrix<Double> & bMat,const LocPointMapped& lp, BaseFE* ptFe){

  //check if lp is surface and ptFe is volume
  assert(lp.isSurface);
  assert(D == ptFe->shape_.dim);

  UInt numFncs = ptFe->GetNumFncs();

  FE *fe = (static_cast<FE*>(ptFe));

  if(fe->GetFeSpaceName() == "HCurl"){
    // Set correct size of matrix B and initialise with zeros
    bMat.Resize( 3, numFncs);
    bMat.Init();

    Matrix<Double> v;
    FeHCurl *feHC = (static_cast<FeHCurl*>(ptFe));
    feHC->GetShFnc( bMat, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem());
  }else{
    // Set correct size of matrix B and initialize with zeros
    bMat.Resize( DIM_DOF, numFncs * DIM_DOF );
    bMat.Init();

    Vector<Double> s;
    for(UInt d = 0; d < DIM_DOF ; d ++){
      fe->GetShFnc( s, lp.lpmVol->lp, lp.lpmVol->shapeMap->GetElem() , d );
      for(UInt sh = 0; sh < numFncs; sh ++){
        bMat[d][sh*DIM_DOF + d] = s[sh];
      }
    }
  }

}


template<class FE,  UInt D, UInt D_DOF, class TYPE>
void SurfaceIdentityOperator<FE,D,D_DOF,TYPE>::
  CalcOpMatTransposed(Matrix<Double> & bMat,
                    const LocPointMapped& lp, BaseFE* ptFe){
  //check if lp is surface and ptFe is volume
  assert(lp.isSurface);
  assert(D == ptFe->shape_.dim);

  UInt numFncs = ptFe->GetNumFncs();
  // Set correct size of matrix B and initialize with zeros
  bMat.Resize( numFncs * DIM_DOF , DIM_DOF );
  bMat.Init();

  // Get derivatives of local shape functions with respect to global
  // coords (format: nrNodes x spaceDim)
  FE *fe = (static_cast<FE*>(ptFe));

  if(fe->GetFeSpaceName() == "HCurl"){
    EXCEPTION("SurfaceIdentityOperator Transposed not implemented for edge elements...just c&p");
  }else{
    Vector<Double> s;
    for(UInt d = 0; d < DIM_DOF ; d ++){
      fe->GetShFnc( s, lp.lpmVol->lp, lp.lpmVol->shapeMap->GetElem() , d );
      for(UInt sh = 0; sh < numFncs; sh ++){
        bMat[sh*DIM_DOF + d][d] = s[sh];
      }
    }
  }
}

template<class FE, UInt D = 1, UInt D_DOF = 1, class TYPE = Double>
class SurfaceMultiIdOp : public BaseBOperator{

public:

  // ------------------
  //  STATIC CONSTANTS
  // ------------------
  //@{
  //! \name Static constants

  //! Order of differentiation
  static const UInt ORDER_DIFF = 0;

  //! Number of components of the problem (scalar, vector)
  static const UInt DIM_DOF = D_DOF;

  //! Dimension of the underlying domain / space
  static const UInt DIM_SPACE = D;

  //! Dimension of the finite element
  static const UInt DIM_ELEM = D;

  //! Dimension of the related material
  static const UInt DIM_D_MAT = 1;
  //@}


  SurfaceMultiIdOp(){
    return;
  }

  SurfaceMultiIdOp(const SurfaceMultiIdOp & other)
   : BaseBOperator(other){
  }

  virtual SurfaceMultiIdOp * Clone(){
    return new SurfaceMultiIdOp(*this);
  }

  virtual ~SurfaceMultiIdOp(){
    return;
  }

  virtual void CalcOpMat(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe);

  virtual void CalcOpMatTransposed(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe);

  //avoid reimplementation of complex operator by making the base class function available
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

template<class FE,  UInt D, UInt D_DOF, class TYPE>
void SurfaceMultiIdOp<FE, D, D_DOF, TYPE>::CalcOpMat(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
{
  //check if lp is surface and ptFe is volume
  assert(lp.isSurface);
  assert(D == ptFe->shape_.dim);

  const UInt numFncs = ptFe->GetNumFncs();

  // Set correct size of matrix B and initialize with zeros
  bMat.Resize(DIM_SPACE, numFncs*DIM_DOF);

  Vector<Double> s;
  FE *fe = (static_cast<FE*>(ptFe));
  for(UInt d = 0; d < DIM_SPACE ; d ++)
  {
    fe->GetShFnc(s, lp.lpmVol->lp, lp.lpmVol->shapeMap->GetElem(), d);
    for(UInt sh = 0; sh < numFncs; sh ++)
    {
      bMat[d][sh] = s[sh];
    }
  }
}

template<class FE,  UInt D, UInt D_DOF, class TYPE>
void SurfaceMultiIdOp<FE, D, D_DOF, TYPE>::CalcOpMatTransposed(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
{
  Matrix<Double> dummyMat;
  CalcOpMat(dummyMat, lp, ptFe);

  dummyMat.Transpose(bMat);
}

template<class FE, UInt D = 1, UInt D_DOF = 1, class TYPE = Double>
class SurfaceNormalDerivOperator : public BaseBOperator{

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
   static const UInt DIM_D_MAT = D_DOF;
   //@}

   SurfaceNormalDerivOperator(){
      return;
    }

   SurfaceNormalDerivOperator(const SurfaceNormalDerivOperator & other)
    : BaseBOperator(other){
   }

   virtual SurfaceNormalDerivOperator * Clone(){
     return new SurfaceNormalDerivOperator(*this);
   }


    virtual ~SurfaceNormalDerivOperator(){
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

template<class FE,  UInt D, UInt D_DOF, class TYPE>
void SurfaceNormalDerivOperator<FE,D,D_DOF,TYPE>::CalcOpMat(Matrix<Double> & bMat,
                                              const LocPointMapped& lp,
                                              BaseFE* ptFe ){
  //check if lp is surface and ptFe is volume
  assert(lp.isSurface);
  assert(D == ptFe->shape_.dim);

  UInt numFncs = ptFe->GetNumFncs();
  // Set correct size of matrix B and initialise with zeros
  bMat.Resize( D_DOF, numFncs );
  bMat.InitValue(0.0);

  // Get derivatives of local shape functions with respect to global
  // coords (format: nrNodes x spaceDim)
  Matrix<Double> xiDx;
  FE *fe = (static_cast<FE*>(ptFe));
  fe->GetGlobDerivShFnc( xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem() , 1 );

  //perform scalar mult with surface normal
  for(UInt d = 0; d < DIM_DOF ; ++d){
    for(UInt d1 = 0; d1 < DIM_SPACE ; ++d1){
      for(UInt sh = 0; sh < numFncs; ++sh){
        bMat[d][sh*DIM_DOF + d] += xiDx[sh][d1] * lp.normal[d1];
      }
    }
  }
}

template<class FE,  UInt D, UInt D_DOF, class TYPE>
void SurfaceNormalDerivOperator<FE,D,D_DOF,TYPE>::CalcOpMatTransposed(Matrix<Double> & bMat,
                                                        const LocPointMapped& lp,
                                                        BaseFE* ptFe ){
  //check if lp is surface and ptFe is volume
  assert(lp.isSurface);
  assert(D == ptFe->shape_.dim);
  Matrix<Double> tmpMat;
  this->CalcOpMat(tmpMat,lp,ptFe);
  bMat = Transpose(tmpMat);
}


//! Calculate the divergence of the vector shape functions scaled with the normal vector
//!    / N_1x*n_x  N_1y*n_x N_1z*n_x   ...\
//! b =| N_1x*n_y  N_1y*n_y N_1z*n_y   ...|
//!    \ N_1x*n_z  N_1y*n_z N_1z*n_z .../
//!  here N_1x denotes the x-derivative of the first
//!  shape function at a given local point and n_x
//!  the x-component of the normal vector
template<class FE, UInt D = 1, class TYPE = Double>
class SurfaceNormalDivOperator : public BaseBOperator{

public:

   // ------------------
   //  STATIC CONSTANTS
   // ------------------
   //@{
   //! \name Static constants

   //! Order of differentiation
   static const UInt ORDER_DIFF = 1;

   //! Number of components of the problem (scalar, vector)
   static const UInt DIM_DOF = D;

   //! Dimension of the underlying domain / space
   static const UInt DIM_SPACE = D;

   //! Dimension of the finite element
   static const UInt DIM_ELEM = D;

   //! Dimension of the related material
   static const UInt DIM_D_MAT = D;
   //@}

   SurfaceNormalDivOperator(){
      return;
    }

   SurfaceNormalDivOperator(const SurfaceNormalDivOperator & other)
    : BaseBOperator(other){
   }

   virtual SurfaceNormalDivOperator * Clone(){
     return new SurfaceNormalDivOperator(*this);
   }


    virtual ~SurfaceNormalDivOperator(){
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

template<class FE,  UInt D, class TYPE>
void SurfaceNormalDivOperator<FE,D,TYPE>::CalcOpMat(Matrix<Double> & bMat,
                                              const LocPointMapped& lp,
                                              BaseFE* ptFe ){
  //check if lp is surface and ptFe is volume
  assert(lp.isSurface);
  assert(D == ptFe->shape_.dim);

  UInt numFncs = ptFe->GetNumFncs();
  // Set correct size of matrix B and initialise with zeros
  bMat.Resize( DIM_D_MAT, numFncs * DIM_SPACE );
  bMat.InitValue(0.0);

  // Check if the normal vector is parallel to the coordinate system
  // In the loop below we check if the sum of all components is equal to 1 in order to check if the normal vector lines up with one coordinate axis
  // At the moment this check is necessary due to some problems with the doNothing BC
  // It seems that the already implemented SurfaceNormalStressOperator is not working properly since results differ when rotating a testcase e.g. 30° in plane
  Double normalSum = 0.0;

  // Get derivatives of local shape functions with respect to global
  // coords (format: nrNodes x spaceDim)
  Matrix<Double> xiDx;
  FE *fe = (static_cast<FE*>(ptFe));
  fe->GetGlobDerivShFnc( xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem() , 1 );

  for(UInt iDim = 0; iDim < DIM_SPACE; ++iDim) {
    for( UInt i = 0; i < numFncs; ++i ) {
      bMat[iDim][i*DIM_SPACE + iDim] += xiDx[i][iDim] * lp.normal[iDim];
      bMat[DIM_SPACE-1-iDim][i*DIM_SPACE + iDim] += xiDx[i][iDim] * lp.normal[DIM_SPACE-1-iDim];
    }
    normalSum += lp.normal[iDim];
  }

  if( !(abs(abs(normalSum)-1)<1e-12) ) {
    //EXCEPTION("Please line up your BCs with the coordinate system (normal vector in direction of coordinate axis) since atm this BC misbehaves for arbitrary angles!");
  }
}

template<class FE,  UInt D, class TYPE>
void SurfaceNormalDivOperator<FE,D,TYPE>::CalcOpMatTransposed(Matrix<Double> & bMat,
                                                        const LocPointMapped& lp,
                                                        BaseFE* ptFe ){
  Matrix<Double> tmpMat;
  this->CalcOpMat(tmpMat,lp,ptFe);
  bMat = Transpose(tmpMat);
}


template<class FE, UInt D = 1, UInt D_DOF = 1, class TYPE = Double>
class SurfaceIdentityOperatorScaledBySurface : public BaseBOperator{

public:

   // ------------------
   //  STATIC CONSTANTS
   // ------------------
   //@{
   //! \name Static constants

   //! Order of differentiation
   static const UInt ORDER_DIFF = 0;

   //! Number of components of the problem (scalar, vector)
   static const UInt DIM_DOF = D_DOF;

   //! Dimension of the underlying domain / space
   static const UInt DIM_SPACE = D;

   //! Dimension of the finite element
   static const UInt DIM_ELEM = D;

   //! Dimension of the related material
   static const UInt DIM_D_MAT = D_DOF;
   //@}

   SurfaceIdentityOperatorScaledBySurface(){
      return;
    }

   SurfaceIdentityOperatorScaledBySurface(const SurfaceIdentityOperatorScaledBySurface & other)
    : BaseBOperator(other){
   }

   virtual SurfaceIdentityOperatorScaledBySurface * Clone(){
     return new SurfaceIdentityOperatorScaledBySurface(*this);
   }

    virtual ~SurfaceIdentityOperatorScaledBySurface(){
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

template<class FE,  UInt D, UInt D_DOF, class TYPE>
void SurfaceIdentityOperatorScaledBySurface<FE,D,D_DOF,TYPE>::
  CalcOpMat(Matrix<Double> & bMat,const LocPointMapped& lp, BaseFE* ptFe){

  //check if lp is surface and ptFe is volume
  assert(lp.isSurface);
  assert(D == ptFe->shape_.dim);

  //Double factor2 = lp.shapeMap->CalcVolume();
  const NcSurfElem* sElem = dynamic_cast<const NcSurfElem*>(lp.ptEl);

  shared_ptr<ElemShapeMap> esm1 = lp.shapeMap->GetGrid()->GetElemShapeMap(sElem->neighbors[0].get(),true);
  shared_ptr<ElemShapeMap> esm2 = lp.shapeMap->GetGrid()->GetElemShapeMap(sElem->neighbors[1].get(),true);

  /*
   * Enable depth scaling for 2d plane case here, as it is done in NACS
   */
  Double v1 = esm1->CalcVolume(true);
  Double v2 = esm1->CalcVolume(true);

  Double factor2 = 2/(v1+v2);

  UInt numFncs = ptFe->GetNumFncs();

  // Set correct size of matrix B and initialize with zeros
  bMat.Resize( DIM_DOF, numFncs * DIM_DOF );
  bMat.Init();

  Vector<Double> s;
  FE *fe = (static_cast<FE*>(ptFe));
  for(UInt d = 0; d < DIM_DOF ; d ++){
    fe->GetShFnc( s, lp.lpmVol->lp, lp.lpmVol->shapeMap->GetElem() , d );
    for(UInt sh = 0; sh < numFncs; sh ++){
      bMat[d][sh*DIM_DOF + d] = s[sh] / factor2;
    }
  }
}


template<class FE,  UInt D, UInt D_DOF, class TYPE>
void SurfaceIdentityOperatorScaledBySurface<FE,D,D_DOF,TYPE>::
  CalcOpMatTransposed(Matrix<Double> & bMat,
                    const LocPointMapped& lp, BaseFE* ptFe){
  //check if lp is surface and ptFe is volume
  assert(lp.isSurface);
  assert(D == ptFe->shape_.dim);

  /*
   * Enable depth scaling for 2d plane case here, as it is done in NACS
   */
  Double factor2 = lp.shapeMap->CalcVolume(true);

  UInt numFncs = ptFe->GetNumFncs();
  // Set correct size of matrix B and initialize with zeros
  bMat.Resize( numFncs * DIM_DOF , DIM_DOF );
  bMat.Init();

  // Get derivatives of local shape functions with respect to global
  // coords (format: nrNodes x spaceDim)
  Vector<Double> s;
  FE *fe = (static_cast<FE*>(ptFe));
  for(UInt d = 0; d < DIM_DOF ; d ++){
    fe->GetShFnc( s, lp.lpmVol->lp, lp.lpmVol->shapeMap->GetElem() , d );
    for(UInt sh = 0; sh < numFncs; sh ++){
      bMat[sh*DIM_DOF + d][d] = s[sh] / factor2;
    }
  }
}


template<class FE, UInt D = 1, UInt D_DOF = 1, class TYPE = Double>
class SurfaceCurlOperator : public BaseBOperator{

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
   static const UInt DIM_D_MAT = D_DOF;
   //@}

   SurfaceCurlOperator(){
      return;
    }

   SurfaceCurlOperator(const SurfaceCurlOperator & other)
    : BaseBOperator(other){
   }

   virtual SurfaceCurlOperator * Clone(){
     return new SurfaceCurlOperator(*this);
   }


    virtual ~SurfaceCurlOperator(){
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



template<class FE,  UInt D, UInt D_DOF, class TYPE>
void SurfaceCurlOperator<FE,D,D_DOF,TYPE>::CalcOpMat(Matrix<Double> & bMat,
                                              const LocPointMapped& lp,
                                              BaseFE* ptFe ){
  //check if lp is surface and ptFe is volume
  assert(lp.isSurface);
  assert(D == ptFe->shape_.dim);

  UInt numFncs = ptFe->GetNumFncs();
  // Set correct size of matrix B and initialise with zeros
  bMat.Resize( 3, numFncs);
  bMat.Init();

  Matrix<Double> v;
  FE *fe = (static_cast<FE*>(ptFe));

  if(fe->GetFeSpaceName() != "HCurl"){
    EXCEPTION("SurfaceCurlOperator only for HCurl function space");
  }

  FeHCurl *feHC = (static_cast<FeHCurl*>(ptFe));
  feHC->GetCurlShFnc( bMat, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem());

}

template<class FE,  UInt D, UInt D_DOF, class TYPE>
void SurfaceCurlOperator<FE,D,D_DOF,TYPE>::CalcOpMatTransposed(Matrix<Double> & bMat,
                                                        const LocPointMapped& lp,
                                                        BaseFE* ptFe ){
  //check if lp is surface and ptFe is volume
  assert(lp.isSurface);
  assert(D == ptFe->shape_.dim);
  Matrix<Double> tmpMat;
  this->CalcOpMat(tmpMat,lp,ptFe);
  bMat = Transpose(tmpMat);
}



template<class FE, UInt D = 1, UInt D_DOF = 1, class TYPE = Double>
class SurfaceCurlNormalOperator : public BaseBOperator{

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
   static const UInt DIM_D_MAT = D_DOF;
   //@}

   SurfaceCurlNormalOperator(){
      return;
    }

   SurfaceCurlNormalOperator(const SurfaceCurlNormalOperator & other)
    : BaseBOperator(other){
   }

   virtual SurfaceCurlNormalOperator * Clone(){
     return new SurfaceCurlNormalOperator(*this);
   }


    virtual ~SurfaceCurlNormalOperator(){
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



template<class FE,  UInt D, UInt D_DOF, class TYPE>
void SurfaceCurlNormalOperator<FE,D,D_DOF,TYPE>::CalcOpMat(Matrix<Double> & bMat,
                                              const LocPointMapped& lp,
                                              BaseFE* ptFe ){
  //check if lp is surface and ptFe is volume
  assert(lp.isSurface);
  assert(D == ptFe->shape_.dim);

  UInt numFncs = ptFe->GetNumFncs();
  // Set correct size of matrix B and initialise with zeros
  bMat.Resize( 3, numFncs);
  bMat.Init();

  Matrix<Double> curl;
  FE *fe = (static_cast<FE*>(ptFe));

  if(fe->GetFeSpaceName() != "HCurl"){
    EXCEPTION("SurfaceCurlOperator only for HCurl function space");
  }

  FeHCurl *feHC = (static_cast<FeHCurl*>(ptFe));
  feHC->GetCurlShFnc( curl, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(),1);

  // normal vector
  Vector<Double> n = lp.normal;
  Vector<Double> tmp;
  Vector<Double> curlXn ;
  for(UInt sh = 0; sh < numFncs ; sh ++){
    curl.GetCol(tmp, sh);

    tmp.CrossProduct( n, curlXn);
    for(UInt d = 0; d < 3; d ++){
      bMat[d][sh] = curlXn[d];
    }
  }
}

template<class FE,  UInt D, UInt D_DOF, class TYPE>
void SurfaceCurlNormalOperator<FE,D,D_DOF,TYPE>::CalcOpMatTransposed(Matrix<Double> & bMat,
                                                        const LocPointMapped& lp,
                                                        BaseFE* ptFe ){
  //check if lp is surface and ptFe is volume
  assert(lp.isSurface);
  assert(D == ptFe->shape_.dim);
  Matrix<Double> tmpMat;
  this->CalcOpMat(tmpMat,lp,ptFe);
  bMat = Transpose(tmpMat);
}




template<class FE, UInt D = 1, UInt D_DOF = 1, class TYPE = Double>
class SurfaceNormalOperator : public BaseBOperator{

public:

   // ------------------
   //  STATIC CONSTANTS
   // ------------------
   //@{
   //! \name Static constants

   //! Order of differentiation
   static const UInt ORDER_DIFF = 0;

   //! Number of components of the problem (scalar, vector)
   static const UInt DIM_DOF = D_DOF;

   //! Dimension of the underlying domain / space
   static const UInt DIM_SPACE = D;

   //! Dimension of the finite element
   static const UInt DIM_ELEM = D;

   //! Dimension of the related material
   static const UInt DIM_D_MAT = D_DOF;
   //@}

   SurfaceNormalOperator(){
      return;
    }

   SurfaceNormalOperator(const SurfaceNormalOperator & other)
    : BaseBOperator(other){
   }

   virtual SurfaceNormalOperator * Clone(){
     return new SurfaceNormalOperator(*this);
   }


    virtual ~SurfaceNormalOperator(){
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


template<class FE,  UInt D, UInt D_DOF, class TYPE>
void SurfaceNormalOperator<FE,D,D_DOF,TYPE>::CalcOpMat(Matrix<Double> & bMat,
                                              const LocPointMapped& lp,
                                              BaseFE* ptFe ){
  //check if lp is surface and ptFe is volume
  assert(lp.isSurface);
  assert(D == ptFe->shape_.dim);

  UInt numFncs = ptFe->GetNumFncs();
  // Set correct size of matrix B and initialise with zeros
  // Set correct size of matrix B and initialise with zeros
  bMat.Resize( 3, numFncs);
  bMat.Init();

  Matrix<Double> v;
  FE *fe = (static_cast<FE*>(ptFe));

  if(fe->GetFeSpaceName() != "HCurl"){
    EXCEPTION("SurfaceNormalOperator only for HCurl function space");
  }

  // normal vector
  Vector<Double> n = lp.normal;

  FeHCurl *feHC = (static_cast<FeHCurl*>(ptFe));
  feHC->GetShFnc( v, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);



  Vector<Double> tmp, curl;
  Vector<Double> transNormal ;
  for(UInt sh = 0; sh < numFncs ; sh ++){
    v.GetCol(tmp, sh);
    transNormal = lp.normal;
    tmp.CrossProduct( transNormal ,curl);
    for(UInt d = 0; d < 3; d ++){
      bMat[d][sh] = curl[d];
    }
  }

}

template<class FE,  UInt D, UInt D_DOF, class TYPE>
void SurfaceNormalOperator<FE,D,D_DOF,TYPE>::CalcOpMatTransposed(Matrix<Double> & bMat,
                                                        const LocPointMapped& lp,
                                                        BaseFE* ptFe ){
  //check if lp is surface and ptFe is volume
  assert(lp.isSurface);
  assert(D == ptFe->shape_.dim);
  Matrix<Double> tmpMat;
  this->CalcOpMat(tmpMat,lp,ptFe);
  bMat = Transpose(tmpMat);
}



template<class FE, UInt D = 1, UInt D_DOF = 1, class TYPE = Double>
class SurfaceTangentialIncompressibleStrainOperator2D : public BaseBOperator{

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

  
  SurfaceTangentialIncompressibleStrainOperator2D(){
    return;
  }

  SurfaceTangentialIncompressibleStrainOperator2D(const SurfaceTangentialIncompressibleStrainOperator2D & other)
    : BaseBOperator(other){
  }

  virtual SurfaceTangentialIncompressibleStrainOperator2D * Clone(){
    return new SurfaceTangentialIncompressibleStrainOperator2D(*this);
  }

  virtual ~SurfaceTangentialIncompressibleStrainOperator2D(){
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
void SurfaceTangentialIncompressibleStrainOperator2D<FE, D, D_DOF, TYPE>::CalcOpMat(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
{
  //check if lp is surface and ptFe is volume
  assert(lp.isSurface);
  assert(D == ptFe->shape_.dim);
  assert(D==D_DOF);
  assert(D==2);

  UInt numFncs = ptFe->GetNumFncs();

  // here we calculate (I - n x n) \cdot ((\nabla u + (\nabla u)^T) \cdot n)
  // This is essentially a projection into the tangential plane of the normal derivative of the surface velocity plus its transposed part

  // E.g. for x entry of 3D version
  /*for(UInt i = 0; i < numFncs; ++i) {
    // x entry scaled with u_x
    bMat[0,i*D_DOF] = (1-nVec[0]*nVec[0])*(2*xiDx[i][0]*nVec[0]+xiDx[i][1]*nVec[1]+xiDx[i][2]*nVec[2])-xiDx[i][1]*nVec[0]*nVec[0]*nVec[1]-xiDx[i][2]*nVec[0]*nVec[0]*nVec[2];
    // x entry scaled with u_y
    bMat[0,i*D_DOF+1] = (1-nVec[0]*nVec[0])*xiDx[i][0]*nVec[1]-nVec[0]*nVec[1]*(xiDx[i][0]*nVec[0]+2*xiDx[i][1]*nVec[1]+xiDx[i][2]*nVec[2])-xiDx[i][2]*nVec[0]*nVec[1]*nVec[2];
    // x entry scaled with u_z
    bMat[0,i*D_DOF+2] = (1-nVec[0]*nVec[0])*xiDx[i][0]*nVec[2]-nVec[0]*nVec[2]*(xiDx[i][0]*nVec[0]+xiDx[i][1]*nVec[1]+2*xiDx[i][2]*nVec[2])-xiDx[i][1]*nVec[0]*nVec[1]*nVec[2]
  }*/

  Matrix<Double> xiDx;
  FE *fe = (static_cast<FE*>(ptFe));
  bMat.Resize(D, D_DOF*numFncs);
  fe->GetGlobDerivShFnc(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
  Vector<Double> nVec = lp.normal;

  // Check if the normal vector is parallel to the coordinate system
  // In the loop below we check if the sum of all components is equal to 1 in order to check if the normal vector lines up with one coordinate axis
  // At the moment this check is necessary due to some problems with the doNothing BC
  // It seems that the already implemented SurfaceNormalStressOperator is not working properly since results differ when rotating a testcase e.g. 30° in plane
  Double normalSum = 0.0;
  
  //Checking the normal vector befor Operator construction
  for(UInt iDim = 0; iDim < DIM_SPACE; ++iDim) {
    normalSum += lp.normal[iDim];
  }

  if( !(abs(abs(normalSum)-1)<1e-12) ) {
   //EXCEPTION("Please line up your BCs with the coordinate system (normal vector in direction of coordinate axis) since atm this BC misbehaves for arbitrary angles!");
  }

  for(UInt i = 0; i < numFncs; ++i) {
    // x entry scaled with u_x
    bMat[0][i*D_DOF] = (1-nVec[0]*nVec[0])*(2*xiDx[i][0]*nVec[0]+xiDx[i][1]*nVec[1])-xiDx[i][1]*nVec[0]*nVec[0]*nVec[1];
    // x entry scaled with u_y
    bMat[0][i*D_DOF+1] = (1-nVec[0]*nVec[0])*xiDx[i][0]*nVec[1]-nVec[0]*nVec[1]*(xiDx[i][0]*nVec[0]+2*xiDx[i][1]*nVec[1]);
    // y entry scaled with u_x
    bMat[1][i*D_DOF] = (1-nVec[1]*nVec[1])*xiDx[i][1]*nVec[0]-nVec[0]*nVec[1]*(2*xiDx[i][0]*nVec[0]+xiDx[i][1]*nVec[1]);
    // y entry scaled with u_y
    bMat[1][i*D_DOF+1] = (1-nVec[1]*nVec[1])*(xiDx[i][0]*nVec[0]+2*xiDx[i][1]*nVec[1])-xiDx[i][0]*nVec[0]*nVec[1]*nVec[1];
  }
}

template<class FE, UInt D, UInt D_DOF, class TYPE>
void SurfaceTangentialIncompressibleStrainOperator2D<FE, D, D_DOF, TYPE>::CalcOpMatTransposed(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
{
  Matrix<Double> dummyMat;
  CalcOpMat(dummyMat, lp, ptFe);

  dummyMat.Transpose(bMat);
}

template<class FE, UInt D = 1, UInt D_DOF = 1, class TYPE = Double>
class SurfaceTangentialCompressibleStrainOperator2D : public BaseBOperator{

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

  
  SurfaceTangentialCompressibleStrainOperator2D(){
    return;
  }

  SurfaceTangentialCompressibleStrainOperator2D(const SurfaceTangentialCompressibleStrainOperator2D & other)
    : BaseBOperator(other){
  }

  virtual SurfaceTangentialCompressibleStrainOperator2D * Clone(){
    return new SurfaceTangentialCompressibleStrainOperator2D(*this);
  }

  virtual ~SurfaceTangentialCompressibleStrainOperator2D(){
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
void SurfaceTangentialCompressibleStrainOperator2D<FE, D, D_DOF, TYPE>::CalcOpMat(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe){
  //check if lp is surface and ptFe is volume
  assert(lp.isSurface);
  assert(D == ptFe->shape_.dim);
  assert(D==D_DOF);
  assert(D==2);

  UInt numFncs = ptFe->GetNumFncs();
 
  Matrix<Double> xiDx;
  FE *fe = (static_cast<FE*>(ptFe));
  bMat.Resize(D, D_DOF*numFncs);
  fe->GetGlobDerivShFnc(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
  Vector<Double> nVec = lp.normal;

  // Check if the normal vector is parallel to the coordinate system
  // In the loop below we check if the sum of all components is equal to 1 in order to check if the normal vector lines up with one coordinate axis
  // At the moment this check is necessary due to some problems with the doNothing BC
  // It seems that the already implemented SurfaceNormalStressOperator is not working properly since results differ when rotating a testcase e.g. 30° in plane
  Double normalSum = 0.0;
  
  //Checking the normal vector befor Operator construction
  for(UInt iDim = 0; iDim < DIM_SPACE; ++iDim) {
    normalSum += lp.normal[iDim];
  }

  if( !(abs(abs(normalSum)-1)<1e-12) ) {
    //EXCEPTION("Please line up your BCs with the coordinate system (normal vector in direction of coordinate axis) since atm this BC misbehaves for arbitrary angles!");
  }

  Double coefLambdaOverMu;
  this->coef_->GetScalar(coefLambdaOverMu,lp);

  // here we calculate (I - n x n) \cdot ((\nabla u + (\nabla u)^T) \cdot n)+(kappa-2/3mu) div(u)I
  // This is essentially a projection into the tangential plane of the normal derivative of the surface velocity plus its transposed part
  for(UInt i = 0; i < numFncs; ++i) {
    // x entry scaled with u_x
    bMat[0][i*D_DOF] = (1-nVec[0]*nVec[0])*(2*xiDx[i][0]*nVec[0]+xiDx[i][1]*nVec[1])-xiDx[i][1]*nVec[0]*nVec[0]*nVec[1]+(coefLambdaOverMu-2/3)*(1-nVec[0]*nVec[0]-nVec[1]*nVec[1])*nVec[0]*xiDx[i][0];
    // x entry scaled with u_y
    bMat[0][i*D_DOF+1] = (1-nVec[0]*nVec[0])*xiDx[i][0]*nVec[1]-nVec[0]*nVec[1]*(xiDx[i][0]*nVec[0]+2*xiDx[i][1]*nVec[1])+(coefLambdaOverMu-2/3)*(1-nVec[0]*nVec[0]-nVec[1]*nVec[1])*nVec[0]*xiDx[i][1];
    // y entry scaled with u_x
    bMat[1][i*D_DOF] = (1-nVec[1]*nVec[1])*xiDx[i][1]*nVec[0]-nVec[0]*nVec[1]*(2*xiDx[i][0]*nVec[0]+xiDx[i][1]*nVec[1])+(coefLambdaOverMu-2/3)*(1-nVec[0]*nVec[0]-nVec[1]*nVec[1])*nVec[1]*xiDx[i][0];
    // y entry scaled with u_y
    bMat[1][i*D_DOF+1] = (1-nVec[1]*nVec[1])*(xiDx[i][0]*nVec[0]+2*xiDx[i][1]*nVec[1])-xiDx[i][0]*nVec[0]*nVec[1]*nVec[1]+(coefLambdaOverMu-2/3)*(1-nVec[0]*nVec[0]-nVec[1]*nVec[1])*nVec[1]*xiDx[i][1];
  }
}

template<class FE, UInt D, UInt D_DOF, class TYPE>
void SurfaceTangentialCompressibleStrainOperator2D<FE, D, D_DOF, TYPE>::CalcOpMatTransposed(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
{
  Matrix<Double> dummyMat;
  CalcOpMat(dummyMat, lp, ptFe);

  dummyMat.Transpose(bMat);
}

template<class FE, UInt D = 1, UInt D_DOF = 1, class TYPE = Double>
class SurfaceTangentialHollowIncompressibleStrainOperator2D : public BaseBOperator{

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

  
  SurfaceTangentialHollowIncompressibleStrainOperator2D(){
    return;
  }

  SurfaceTangentialHollowIncompressibleStrainOperator2D(const SurfaceTangentialHollowIncompressibleStrainOperator2D & other)
    : BaseBOperator(other){
  }

  virtual SurfaceTangentialHollowIncompressibleStrainOperator2D * Clone(){
    return new SurfaceTangentialHollowIncompressibleStrainOperator2D(*this);
  }

  virtual ~SurfaceTangentialHollowIncompressibleStrainOperator2D(){
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
void SurfaceTangentialHollowIncompressibleStrainOperator2D<FE, D, D_DOF, TYPE>::CalcOpMat(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
{
  //check if lp is surface and ptFe is volume
  assert(lp.isSurface);
  assert(D == ptFe->shape_.dim);
  assert(D==D_DOF);
  assert(D==2);

  UInt numFncs = ptFe->GetNumFncs();

  // here we calculate (I - n x n) \cdot (((\nabla u + (\nabla u)^T) - diag(\nabla u +  (\nabla u)^T)) \cdot n)
  // This is essentially a projection into the tangential plane of the normal derivative of the surface velocity plus its transposed part without the diagonal part (--> hollow matrix)

  //-ny*(2*nx^2 - 1)*(Nx*uy + Ny*ux)
  //-nx*(2*ny^2 - 1)*(Nx*uy + Ny*ux)

  Matrix<Double> xiDx;
  FE *fe = (static_cast<FE*>(ptFe));
  bMat.Resize(D, D_DOF*numFncs);
  fe->GetGlobDerivShFnc(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
  Vector<Double> nVec = lp.normal;

  // Check if the normal vector is parallel to the coordinate system
  // In the loop below we check if the sum of all components is equal to 1 in order to check if the normal vector lines up with one coordinate axis
  // At the moment this check is necessary due to some problems with the doNothing BC
  // It seems that the already implemented SurfaceNormalStressOperator is not working properly since results differ when rotating a testcase e.g. 30° in plane
  Double normalSum = 0.0;
  
  //Checking the normal vector befor Operator construction
  for(UInt iDim = 0; iDim < DIM_SPACE; ++iDim) {
    normalSum += lp.normal[iDim];
  }

  if( !(abs(abs(normalSum)-1)<1e-12) ) {
   //EXCEPTION("Please line up your BCs with the coordinate system (normal vector in direction of coordinate axis) since atm this BC misbehaves for arbitrary angles!");
  }

  for(UInt i = 0; i < numFncs; ++i) {
    // x entry scaled with u_x
    //-ny*(2*nx^2 - 1)*(Ny*ux)
    bMat[0][i*D_DOF] = -nVec[1]*(2*nVec[0]*nVec[0]-1)*xiDx[i][1];
    // x entry scaled with u_y
    //-ny*(2*nx^2 - 1)*(Nx*uy)
    bMat[0][i*D_DOF+1] = -nVec[1]*(2*nVec[0]*nVec[0]-1)*xiDx[i][0];
    // y entry scaled with u_x
    //-nx*(2*ny^2 - 1)*(Ny*ux)
    bMat[1][i*D_DOF] = -nVec[0]*(2*nVec[1]*nVec[1]-1)*xiDx[i][1];
    // y entry scaled with u_y
    //-nx*(2*ny^2 - 1)*(Nx*uy)
    bMat[1][i*D_DOF+1] = -nVec[0]*(2*nVec[1]*nVec[1]-1)*xiDx[i][0];
  }
}

template<class FE, UInt D, UInt D_DOF, class TYPE>
void SurfaceTangentialHollowIncompressibleStrainOperator2D<FE, D, D_DOF, TYPE>::CalcOpMatTransposed(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
{
  Matrix<Double> dummyMat;
  CalcOpMat(dummyMat, lp, ptFe);

  dummyMat.Transpose(bMat);
}


//--------------------------------------------------------
//Implementation of 3D Surface Tangential Strain Operator
//--------------------------------------------------------



template<class FE, UInt D = 1, UInt D_DOF = 1, class TYPE = Double>
class SurfaceTangentialIncompressibleStrainOperator3D : public BaseBOperator{

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

  
  SurfaceTangentialIncompressibleStrainOperator3D(){
    return;
  }

  SurfaceTangentialIncompressibleStrainOperator3D(const SurfaceTangentialIncompressibleStrainOperator3D & other)
    : BaseBOperator(other){
  }

  virtual SurfaceTangentialIncompressibleStrainOperator3D * Clone(){
    return new SurfaceTangentialIncompressibleStrainOperator3D(*this);
  }

  virtual ~SurfaceTangentialIncompressibleStrainOperator3D(){
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
void SurfaceTangentialIncompressibleStrainOperator3D<FE, D, D_DOF, TYPE>::CalcOpMat(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
{
  //check if lp is surface and ptFe is volume
  assert(lp.isSurface);
  assert(D == ptFe->shape_.dim);
  assert(D==D_DOF);
  assert(D==3);

  UInt numFncs = ptFe->GetNumFncs();

  // here we calculate (I - n x n) \cdot ((\nabla u + (\nabla u)^T) \cdot n)
  // This is essentially a projection into the tangential plane of the normal derivative of the surface velocity plus its transposed part

  Matrix<Double> xiDx;
  FE *fe = (static_cast<FE*>(ptFe));
  bMat.Resize(D, D_DOF*numFncs);
  fe->GetGlobDerivShFnc(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
  Vector<Double> nVec = lp.normal;

  // Check if the normal vector is parallel to the coordinate system
  // In the loop below we check if the sum of all components is equal to 1 in order to check if the normal vector lines up with one coordinate axis
  // At the moment this check is necessary due to some problems with the doNothing BC
  // It seems that the already implemented SurfaceNormalStressOperator is not working properly since results differ when rotating a testcase e.g. 30° in plane
  Double normalSum = 0.0;
  
  //Checking the normal vector befor Operator construction
  for(UInt iDim = 0; iDim < DIM_SPACE; ++iDim) {
    normalSum += lp.normal[iDim];
  }

  if( !(abs(abs(normalSum)-1)<1e-12) ) {
    //EXCEPTION("Please line up your BCs with the coordinate system (normal vector in direction of coordinate axis) since atm this BC misbehaves for arbitrary angles!");
  }

  for(UInt i = 0; i < numFncs; ++i) {
    // x entry scaled with u_x
    bMat[0][i*D_DOF] = (- (nVec[0]*nVec[0] - 1)*(2*xiDx[i][0]*nVec[0] + xiDx[i][1]*nVec[1] + xiDx[i][2]*nVec[2]) - xiDx[i][1]*nVec[0]*nVec[0]*nVec[1] - xiDx[i][2]*nVec[0]*nVec[0]*nVec[2]);
    // x entry scaled with u_y
    bMat[0][i*D_DOF+1] = (- xiDx[i][0]*nVec[1]*(nVec[0]*nVec[0] - 1) - nVec[0]*nVec[1]*(xiDx[i][0]*nVec[0] + 2*xiDx[i][1]*nVec[1] + xiDx[i][2]*nVec[2]) - xiDx[i][2]*nVec[0]*nVec[1]*nVec[2]);
    // x entry scaled with u_z
    bMat[0][i*D_DOF+2] = (- xiDx[i][0]*nVec[2]*(nVec[0]*nVec[0] - 1) - nVec[0]*nVec[2]*(xiDx[i][0]*nVec[0] + xiDx[i][1]*nVec[1] + 2*xiDx[i][2]*nVec[2]) - xiDx[i][1]*nVec[0]*nVec[1]*nVec[2]);
    // y entry scaled with u_x
    bMat[1][i*D_DOF] =  (- xiDx[i][1]*nVec[0]*(nVec[0]*nVec[0] - 1) - nVec[0]*nVec[1]*(2*xiDx[i][0]*nVec[0] + xiDx[i][1]*nVec[1] + xiDx[i][2]*nVec[2]) - xiDx[i][2]*nVec[0]*nVec[1]*nVec[2]);
    // y entry scaled with u_y
    bMat[1][i*D_DOF+1] = (- (nVec[1]*nVec[1] - 1)*(xiDx[i][0]*nVec[0] + 2*xiDx[i][1]*nVec[1] + xiDx[i][2]*nVec[2]) - xiDx[i][0]*nVec[0]*nVec[1]*nVec[1] - xiDx[i][2]*nVec[1]*nVec[1]*nVec[2]);
    // y entry scaled with u_z
    bMat[1][i*D_DOF+2] = (- xiDx[i][1]*nVec[2]*(nVec[1]*nVec[1] - 1) - nVec[1]*nVec[2]*(xiDx[i][0]*nVec[0] + xiDx[i][1]*nVec[1] + 2*xiDx[i][2]*nVec[2]) - xiDx[i][0]*nVec[0]*nVec[1]*nVec[2]);
    // z entry scaled with u_x
    bMat[2][i*D_DOF] = (- xiDx[i][2]*nVec[0]*(nVec[2]*nVec[2] - 1) - nVec[0]*nVec[2]*(2*xiDx[i][0]*nVec[0] + xiDx[i][1]*nVec[1] + xiDx[i][2]*nVec[2]) - xiDx[i][1]*nVec[0]*nVec[1]*nVec[2]);
    // z entry scaled with u_y
    bMat[2][i*D_DOF+1] = (- xiDx[i][2]*nVec[1]*(nVec[2]*nVec[2] - 1) - nVec[1]*nVec[2]*(xiDx[i][0]*nVec[0] + 2*xiDx[i][1]*nVec[1] + xiDx[i][2]*nVec[2]) - xiDx[i][0]*nVec[0]*nVec[1]*nVec[2]);
    // z entry scaled with u_z
    bMat[2][i*D_DOF+2] = (- (nVec[2]*nVec[2] - 1)*(xiDx[i][0]*nVec[0] + xiDx[i][1]*nVec[1] + 2*xiDx[i][2]*nVec[2]) - xiDx[i][0]*nVec[0]*nVec[2]*nVec[2] - xiDx[i][1]*nVec[1]*nVec[2]*nVec[2]);
  }
}

template<class FE, UInt D, UInt D_DOF, class TYPE>
void SurfaceTangentialIncompressibleStrainOperator3D<FE, D, D_DOF, TYPE>::CalcOpMatTransposed(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
{
  Matrix<Double> dummyMat;
  CalcOpMat(dummyMat, lp, ptFe);

  dummyMat.Transpose(bMat);
}


template<class FE, UInt D = 1, UInt D_DOF = 1, class TYPE = Double>
class SurfaceTangentialCompressibleStrainOperator3D : public BaseBOperator{

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

  
  SurfaceTangentialCompressibleStrainOperator3D(){
    return;
  }

  SurfaceTangentialCompressibleStrainOperator3D(const SurfaceTangentialCompressibleStrainOperator3D & other)
    : BaseBOperator(other){
  }

  virtual SurfaceTangentialCompressibleStrainOperator3D * Clone(){
    return new SurfaceTangentialCompressibleStrainOperator3D(*this);
  }

  virtual ~SurfaceTangentialCompressibleStrainOperator3D(){
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
void SurfaceTangentialCompressibleStrainOperator3D<FE, D, D_DOF, TYPE>::CalcOpMat(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
{
  //check if lp is surface and ptFe is volume
  assert(lp.isSurface);
  assert(D == ptFe->shape_.dim);
  assert(D==D_DOF);
  assert(D==3);

  UInt numFncs = ptFe->GetNumFncs();

  // here we calculate (I - n x n) \cdot ((\nabla u + (\nabla u)^T) \cdot n)
  // This is essentially a projection into the tangential plane of the normal derivative of the surface velocity plus its transposed part

  Matrix<Double> xiDx;
  FE *fe = (static_cast<FE*>(ptFe));
  bMat.Resize(D, D_DOF*numFncs);
  fe->GetGlobDerivShFnc(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
  Vector<Double> nVec = lp.normal;

  // Check if the normal vector is parallel to the coordinate system
  // In the loop below we check if the sum of all components is equal to 1 in order to check if the normal vector lines up with one coordinate axis
  // At the moment this check is necessary due to some problems with the doNothing BC
  // It seems that the already implemented SurfaceNormalStressOperator is not working properly since results differ when rotating a testcase e.g. 30° in plane
  Double normalSum = 0.0;
  
  //Checking the normal vector befor Operator construction
  for(UInt iDim = 0; iDim < DIM_SPACE; ++iDim) {
    normalSum += lp.normal[iDim];
  }

  if( !(abs(abs(normalSum)-1)<1e-12) ) {
    //EXCEPTION("Please line up your BCs with the coordinate system (normal vector in direction of coordinate axis) since atm this BC misbehaves for arbitrary angles!");
  }

  Double coefLambdaOverMu;
  this->coef_->GetScalar(coefLambdaOverMu,lp);

  for(UInt i = 0; i < numFncs; ++i) {
    // x entry scaled with u_x
    bMat[0][i*D_DOF] = (- (nVec[0]*nVec[0] - 1)*(2*xiDx[i][0]*nVec[0] + xiDx[i][1]*nVec[1] + xiDx[i][2]*nVec[2]) - xiDx[i][1]*nVec[0]*nVec[0]*nVec[1] - xiDx[i][2]*nVec[0]*nVec[0]*nVec[2]) + (nVec[0]*nVec[0] + nVec[1]*nVec[1] + nVec[2]*nVec[2]- 1) * (coefLambdaOverMu-2/3) * nVec[0] * (xiDx[i][0]);
    // x entry scaled with u_y
    bMat[0][i*D_DOF+1] = (- xiDx[i][0]*nVec[1]*(nVec[0]*nVec[0] - 1) - nVec[0]*nVec[1]*(xiDx[i][0]*nVec[0] + 2*xiDx[i][1]*nVec[1] + xiDx[i][2]*nVec[2]) - xiDx[i][2]*nVec[0]*nVec[1]*nVec[2]) + (nVec[0]*nVec[0] + nVec[1]*nVec[1] + nVec[2]*nVec[2]- 1) * (coefLambdaOverMu-2/3) * nVec[0] * (xiDx[i][1]);
    // x entry scaled with u_z
    bMat[0][i*D_DOF+2] = (- xiDx[i][0]*nVec[2]*(nVec[0]*nVec[0] - 1) - nVec[0]*nVec[2]*(xiDx[i][0]*nVec[0] + xiDx[i][1]*nVec[1] + 2*xiDx[i][2]*nVec[2]) - xiDx[i][1]*nVec[0]*nVec[1]*nVec[2]) + (nVec[0]*nVec[0] + nVec[1]*nVec[1] + nVec[2]*nVec[2]- 1) * (coefLambdaOverMu-2/3) * nVec[0] * (xiDx[i][2]);
    // y entry scaled with u_x
    bMat[1][i*D_DOF] =  (- xiDx[i][1]*nVec[0]*(nVec[0]*nVec[0] - 1) - nVec[0]*nVec[1]*(2*xiDx[i][0]*nVec[0] + xiDx[i][1]*nVec[1] + xiDx[i][2]*nVec[2]) - xiDx[i][2]*nVec[0]*nVec[1]*nVec[2]) + (nVec[0]*nVec[0] + nVec[1]*nVec[1] + nVec[2]*nVec[2]- 1) * (coefLambdaOverMu-2/3) * nVec[1] * (xiDx[i][0]);
    // y entry scaled with u_y
    bMat[1][i*D_DOF+1] = (- (nVec[1]*nVec[1] - 1)*(xiDx[i][0]*nVec[0] + 2*xiDx[i][1]*nVec[1] + xiDx[i][2]*nVec[2]) - xiDx[i][0]*nVec[0]*nVec[1]*nVec[1] - xiDx[i][2]*nVec[1]*nVec[1]*nVec[2]) + (nVec[0]*nVec[0] + nVec[1]*nVec[1] + nVec[2]*nVec[2]- 1) * (coefLambdaOverMu-2/3) * nVec[1] * (xiDx[i][1]);
    // y entry scaled with u_z
    bMat[1][i*D_DOF+2] = (- xiDx[i][1]*nVec[2]*(nVec[1]*nVec[1] - 1) - nVec[1]*nVec[2]*(xiDx[i][0]*nVec[0] + xiDx[i][1]*nVec[1] + 2*xiDx[i][2]*nVec[2]) - xiDx[i][0]*nVec[0]*nVec[1]*nVec[2]) + (nVec[0]*nVec[0] + nVec[1]*nVec[1] + nVec[2]*nVec[2]- 1) * (coefLambdaOverMu-2/3) * nVec[1] * (xiDx[i][2]);
    // z entry scaled with u_x
    bMat[2][i*D_DOF] = (- xiDx[i][2]*nVec[0]*(nVec[2]*nVec[2] - 1) - nVec[0]*nVec[2]*(2*xiDx[i][0]*nVec[0] + xiDx[i][1]*nVec[1] + xiDx[i][2]*nVec[2]) - xiDx[i][1]*nVec[0]*nVec[1]*nVec[2]) + (nVec[0]*nVec[0] + nVec[1]*nVec[1] + nVec[2]*nVec[2]- 1) * (coefLambdaOverMu-2/3) * nVec[2] * (xiDx[i][0]);
    // z entry scaled with u_y
    bMat[2][i*D_DOF+1] = (- xiDx[i][2]*nVec[1]*(nVec[2]*nVec[2] - 1) - nVec[1]*nVec[2]*(xiDx[i][0]*nVec[0] + 2*xiDx[i][1]*nVec[1] + xiDx[i][2]*nVec[2]) - xiDx[i][0]*nVec[0]*nVec[1]*nVec[2]) + (nVec[0]*nVec[0] + nVec[1]*nVec[1] + nVec[2]*nVec[2]- 1) * (coefLambdaOverMu-2/3) * nVec[2] * (xiDx[i][1]);
    // z entry scaled with u_z
    bMat[2][i*D_DOF+2] = (- (nVec[2]*nVec[2] - 1)*(xiDx[i][0]*nVec[0] + xiDx[i][1]*nVec[1] + 2*xiDx[i][2]*nVec[2]) - xiDx[i][0]*nVec[0]*nVec[2]*nVec[2] - xiDx[i][1]*nVec[1]*nVec[2]*nVec[2]) + (nVec[0]*nVec[0] + nVec[1]*nVec[1] + nVec[2]*nVec[2]- 1) * (coefLambdaOverMu-2/3) * nVec[2] * (xiDx[i][2]);
  }
}

template<class FE, UInt D, UInt D_DOF, class TYPE>
void SurfaceTangentialCompressibleStrainOperator3D<FE, D, D_DOF, TYPE>::CalcOpMatTransposed(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
{
  Matrix<Double> dummyMat;
  CalcOpMat(dummyMat, lp, ptFe);

  dummyMat.Transpose(bMat);
}




template<class FE, UInt D = 1, UInt D_DOF = 1, class TYPE = Double>
class SurfaceTangentialHollowIncompressibleStrainOperator3D : public BaseBOperator{

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

  
  SurfaceTangentialHollowIncompressibleStrainOperator3D(){
    return;
  }

  SurfaceTangentialHollowIncompressibleStrainOperator3D(const SurfaceTangentialHollowIncompressibleStrainOperator3D & other)
    : BaseBOperator(other){
  }

  virtual SurfaceTangentialHollowIncompressibleStrainOperator3D * Clone(){
    return new SurfaceTangentialHollowIncompressibleStrainOperator3D(*this);
  }

  virtual ~SurfaceTangentialHollowIncompressibleStrainOperator3D(){
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
void SurfaceTangentialHollowIncompressibleStrainOperator3D<FE, D, D_DOF, TYPE>::CalcOpMat(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
{
  //check if lp is surface and ptFe is volume
  assert(lp.isSurface);
  assert(D == ptFe->shape_.dim);
  assert(D==D_DOF);
  assert(D==3);

  UInt numFncs = ptFe->GetNumFncs();

  // here we calculate (I - n x n) \cdot ((\nabla u + (\nabla u)^T) \cdot n)
  // This is essentially a projection into the tangential plane of the normal derivative of the surface velocity plus its transposed part

  Matrix<Double> xiDx;
  FE *fe = (static_cast<FE*>(ptFe));
  bMat.Resize(D, D_DOF*numFncs);
  fe->GetGlobDerivShFnc(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
  Vector<Double> nVec = lp.normal;

  // Check if the normal vector is parallel to the coordinate system
  // In the loop below we check if the sum of all components is equal to 1 in order to check if the normal vector lines up with one coordinate axis
  // At the moment this check is necessary due to some problems with the doNothing BC
  // It seems that the already implemented SurfaceNormalStressOperator is not working properly since results differ when rotating a testcase e.g. 30° in plane
  Double normalSum = 0.0;
  
  //Checking the normal vector befor Operator construction
  for(UInt iDim = 0; iDim < DIM_SPACE; ++iDim) {
    normalSum += lp.normal[iDim];
  }

  if( !(abs(abs(normalSum)-1)<1e-12) ) {
    //EXCEPTION("Please line up your BCs with the coordinate system (normal vector in direction of coordinate axis) since atm this BC misbehaves for arbitrary angles!");
  }

  for(UInt i = 0; i < numFncs; ++i) {
    // x entry scaled with u_x
    bMat[0][i*D_DOF] = (- (nVec[0]*nVec[0] - 1)*(xiDx[i][1]*nVec[1] + xiDx[i][2]*nVec[2]) - xiDx[i][1]*nVec[0]*nVec[0]*nVec[1] - xiDx[i][2]*nVec[0]*nVec[0]*nVec[2]);
    // x entry scaled with u_y
    bMat[0][i*D_DOF+1] = (- xiDx[i][0]*nVec[1]*(nVec[0]*nVec[0] - 1) - nVec[0]*nVec[1]*(xiDx[i][0]*nVec[0] + xiDx[i][2]*nVec[2]) - xiDx[i][2]*nVec[0]*nVec[1]*nVec[2]);
    // x entry scaled with u_z
    bMat[0][i*D_DOF+2] = (- xiDx[i][0]*nVec[2]*(nVec[0]*nVec[0] - 1) - nVec[0]*nVec[2]*(xiDx[i][0]*nVec[0] + xiDx[i][1]*nVec[1]) - xiDx[i][1]*nVec[0]*nVec[1]*nVec[2]);
    // y entry scaled with u_x
    bMat[1][i*D_DOF] =  (- xiDx[i][1]*nVec[0]*(nVec[0]*nVec[0] - 1) - nVec[0]*nVec[1]*(xiDx[i][1]*nVec[1] + xiDx[i][2]*nVec[2]) - xiDx[i][2]*nVec[0]*nVec[1]*nVec[2]);
    // y entry scaled with u_y
    bMat[1][i*D_DOF+1] = (- (nVec[1]*nVec[1] - 1)*(xiDx[i][0]*nVec[0] + xiDx[i][2]*nVec[2]) - xiDx[i][0]*nVec[0]*nVec[1]*nVec[1] - xiDx[i][2]*nVec[1]*nVec[1]*nVec[2]);
    // y entry scaled with u_z
    bMat[1][i*D_DOF+2] = (- xiDx[i][1]*nVec[2]*(nVec[1]*nVec[1] - 1) - nVec[1]*nVec[2]*(xiDx[i][0]*nVec[0] + xiDx[i][1]*nVec[1]) - xiDx[i][0]*nVec[0]*nVec[1]*nVec[2]);
    // z entry scaled with u_x
    bMat[2][i*D_DOF] = (- xiDx[i][2]*nVec[0]*(nVec[2]*nVec[2] - 1) - nVec[0]*nVec[2]*( xiDx[i][1]*nVec[1] + xiDx[i][2]*nVec[2]) - xiDx[i][1]*nVec[0]*nVec[1]*nVec[2]);
    // z entry scaled with u_y
    bMat[2][i*D_DOF+1] = (- xiDx[i][2]*nVec[1]*(nVec[2]*nVec[2] - 1) - nVec[1]*nVec[2]*(xiDx[i][0]*nVec[0]  + xiDx[i][2]*nVec[2]) - xiDx[i][0]*nVec[0]*nVec[1]*nVec[2]);
    // z entry scaled with u_z
    bMat[2][i*D_DOF+2] = (- (nVec[2]*nVec[2] - 1)*(xiDx[i][0]*nVec[0] + xiDx[i][1]*nVec[1]) - xiDx[i][0]*nVec[0]*nVec[2]*nVec[2] - xiDx[i][1]*nVec[1]*nVec[2]*nVec[2]);
  }
}

template<class FE, UInt D, UInt D_DOF, class TYPE>
void SurfaceTangentialHollowIncompressibleStrainOperator3D<FE, D, D_DOF, TYPE>::CalcOpMatTransposed(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
{
  Matrix<Double> dummyMat;
  CalcOpMat(dummyMat, lp, ptFe);

  dummyMat.Transpose(bMat);
}





}
#endif
