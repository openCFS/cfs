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

  // Set correct size of matrix B and initialize with zeros
  bMat.Resize( DIM_DOF, numFncs * DIM_DOF );
  bMat.Init();

  FE *fe = (static_cast<FE*>(ptFe));

  if(fe->GetFeSpaceName() == "HCurl"){
    EXCEPTION("SurfaceIdentityOperator cannot be used for edge elements");
  }else{
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
    EXCEPTION("SurfaceIdentityOperator cannot be used for edge elements");
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

  Double v1 = esm1->CalcVolume();
  Double v2 = esm1->CalcVolume();

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

  Double factor2 = lp.shapeMap->CalcVolume();

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


/*
 * For the penalty term in HCurl
 */
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
  bMat.InitValue(0.0);

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


/*
 * For the consistency resp. symmetrization term in HCurl
 */
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
  bMat.InitValue(0.0);

  Matrix<Double> v;
  FE *fe = (static_cast<FE*>(ptFe));

  if(fe->GetFeSpaceName() != "HCurl"){
    EXCEPTION("SurfaceNormalOperator only for HCurl function space");
  }

  // normal vector
  Vector<Double> n = lp.normal;

  FeHCurl *feHC = (static_cast<FeHCurl*>(ptFe));
  feHC->GetShFnc( v, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem());


  Vector<Double> tmp, curl;
  for(UInt sh = 0; sh < numFncs ; sh ++){
    Vector<Double> s;
    v.GetCol(tmp, sh);
    tmp.CrossProduct( lp.normal ,curl);
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


}
#endif
