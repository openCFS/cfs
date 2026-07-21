// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     preStressOperator.hh
 *       \brief    Definition of the derivative operator used
 *                 for mechanical prestressing.
 *
 *       \date     Oct 11, 2014
 *       \author   ahueppe
 */
//================================================================================================

#ifndef PRESTRESSOPERATOR_HH_
#define PRESTRESSOPERATOR_HH_

#include "BaseBOperator.hh"

namespace CoupledField{
  //! Calculate the operator for the computation of Piola stress tensor
  //! used in geometric nonlinearity and prestressing
  //!    / N_1x  0    0   N_2x   0   ...\
  //! b =| N_1y  0    0   N_2y   0   ...|
  //!    | N_1z  0    0   N_2z   0   ...|
  //!    |  0   N_1x  0    0    N_2x ...|
  //!    |  0   N_1y  0    0    N_2y ...|
  //!    |  0   N_1z  0    0    N_2z ...|
  //!    |  0    0   N_1x  0     0   ...|
  //!    |  0    0   N_1y  0     0   ...|
  //!    \  0    0   N_1z  0     0   .../
  //!  here N_1x denotes the x-derivative of the first
  //!  shape function at a given local point
  //!  this operator is only valid for vectorial unknowns
  //!  Valid for cartesian coordinate systems.
  //! \tparam FE Type of Finite Element used
  //! \tparam D     Dimension of the problem space and unknown
  //! \tparam TYPE Data type (DOUBLE, COMPLEX)
  template<class FE, UInt D, class TYPE = Double>
  class PiolaStressOperator : public BaseBOperator{

  public:
    // ------------------
    //  STATIC CONSTANTS
    // ------------------
    //@{
    //! \name Static constants

    //! Order of differentiation
    static const UInt ORDER_DIFF = 1;

    //! Number of components of the problem (scalar, vector)
    static const UInt DIM_DOF = D; // m=1 (Skalar), m=2,3 Vektor

    //! Dimension of the underlying domain / space
    static const UInt DIM_SPACE = D; // n=2,3

    //! Dimension of the finite element
    static const UInt DIM_ELEM = D; // Dimension von Referenzelement

    //! Dimension of the related material
    static const UInt DIM_D_MAT = D*D;

    //! Constructor
    //! \param useIcModes Use incompatible modes shape functions
    PiolaStressOperator( bool useICModes = false )
    :  useICModes_(useICModes) {
      this->name_ = "preStressOp";
    }

    //! Copy constructor
    PiolaStressOperator(const PiolaStressOperator & other)
       : BaseBOperator(other),
         useICModes_(other.useICModes_){
      this->name_ = other.name_;
    }

    //! \copydoc BaseBOperator::Clone()
    virtual PiolaStressOperator * Clone(){
      return new PiolaStressOperator(*this);
    }

    //! Destructor
    virtual ~PiolaStressOperator(){

    }

    //! Calculate operator matrix
    virtual void CalcOpMat(Matrix<Double> & bMat,
                           const LocPointMapped& lp,
                           BaseFE* ptFe );

    //! Calculate transposed operator matrix
    virtual void CalcOpMatTransposed(Matrix<Double> & bMat,
                                     const LocPointMapped& lp,
                                     BaseFE* ptFe );

//    //! Calculate operator matrix
//    virtual void CalcOpMat(Matrix<Complex> & bMat,
//                           const LocPointMapped& lp,
//                           BaseFE* ptFe );
//
//    //! Calculate transposed operator matrix
//    virtual void CalcOpMatTransposed(Matrix<Complex> & bMat,
//                                     const LocPointMapped& lp,
//                                     BaseFE* ptFe );

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

    bool useICModes_;

  };

  template<class FE, UInt D, class TYPE>
  void PiolaStressOperator<FE,D,TYPE>::CalcOpMat(Matrix<Double> & bMat,
                                               const LocPointMapped& lp,
                                               BaseFE* ptFe ){
    // Get derivatives of local shape functions with respect to global
    // coords (format: nrNodes x spaceDim)
    Matrix<Double> xiDx;
    FE *fe = (static_cast<FE*>(ptFe));
    if( useICModes_ ) {
      if ( isSurfOpt_ )
        fe->GetGlobDerivShFncICModes(  xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem() , 1 );
      else
        fe->GetGlobDerivShFncICModes( xiDx, lp, lp.shapeMap->GetElem() , 1 );

    } else {
      if ( isSurfOpt_ )
        fe->GetGlobDerivShFnc(  xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem() , 1 );
      else
        fe->GetGlobDerivShFnc( xiDx, lp, lp.shapeMap->GetElem() , 1 );
    }

    const UInt numFncs = xiDx.GetNumRows();
    // Set correct size of matrix B and initialize with zeros
    bMat.Resize( DIM_D_MAT, numFncs * DIM_DOF );
    bMat.Init();

    UInt iFunc = 0;
    for( iFunc = 0; iFunc < numFncs; iFunc++ ) {
      for(UInt iDim1 = 0; iDim1 < DIM_DOF; iDim1++ ) {
        for(UInt iDim2 = 0; iDim2 < DIM_DOF; iDim2++ ) {
          bMat[iDim1*DIM_DOF+iDim2][iFunc*DIM_DOF+iDim1] = xiDx[iFunc][iDim2];
        }
      }
    }
  }

  template<class FE, UInt D, class TYPE>
  void PiolaStressOperator<FE,D,TYPE>::CalcOpMatTransposed(Matrix<Double> & bMat,
                                                 const LocPointMapped& lp,
                                                 BaseFE* ptFe ){
    Matrix<Double> tmpMat;
    this->CalcOpMat(tmpMat,lp,ptFe);
    tmpMat.Transpose(bMat);
  }

  /*****************************************************************************
   * PiolaStressOperatorAxi
   ****************************************************************************/
  
  //! Calculate the operator for the computation of prestressing
  //! in the axi-symmetric case.
  //!    / N_1rr   0   N_2rr   0   ...\
  //! b =| N_1zz   0   N_2zz   0   ...|
  //!    | N_1rz   0   N_2rz   0   ...|
  //!    | N_1pp   0   N_2pp   0   ...|
  //!    |  0    N_1rr  0    N_2rr ...|
  //!    |  0    N_1zz  0    N_2zz ...|
  //!    |  0    N_1rz  0    N_2rz ...|
  //!    \  0    N_1pp  0    N_2pp .../
  //! here N_1rr denotes the r-derivative of the r-component of the first shape
  //! function N_1r at a given local point. This operator is only valid for
  //! vectorial unknowns.
  //! \tparam FE Type of Finite Element used
  //! \tparam TYPE Data type (DOUBLE, COMPLEX)
  template<class FE, class TYPE = Double>
  class PiolaStressOperatorAxi : public BaseBOperator{
  
  public:
    // ------------------
    //  STATIC CONSTANTS
    // ------------------
    //@{
    //! \name Static constants
  
    //! Order of differentiation
    static const UInt ORDER_DIFF = 1;
  
    //! Number of components of the problem (scalar, vector)
    static const UInt DIM_DOF = 2;
  
    //! Dimension of the underlying domain / space
    static const UInt DIM_SPACE = 2;
  
    //! Dimension of the finite element
    static const UInt DIM_ELEM = 2;
  
    //! Dimension of the related material
    static const UInt DIM_D_MAT = 5;
  
    //! Constructor
    //! \param useIcModes Use incompatible modes shape functions
    PiolaStressOperatorAxi( bool useICModes = false )
    :  useICModes_(useICModes)
    {
      this->name_ = "preStressOpAxi";
    }
    
    //! Copy constructor
    PiolaStressOperatorAxi(const PiolaStressOperatorAxi & other)
        : BaseBOperator(other),
          useICModes_(other.useICModes_)
    {
      this->name_ = other.name_;
    }
    
    //! Destructor
    virtual ~PiolaStressOperatorAxi(){}
  
    //! \copydoc BaseBOperator::Clone()
    virtual PiolaStressOperatorAxi * Clone(){
      return new PiolaStressOperatorAxi(*this);
    }

    //! Calculate operator matrix
    virtual void CalcOpMat(Matrix<Double> & bMat,
                           const LocPointMapped& lpm,
                           BaseFE* ptFe );
  
    //! Calculate transposed operator matrix
    virtual void CalcOpMatTransposed(Matrix<Double> & bMat,
                                     const LocPointMapped& lpm,
                                     BaseFE* ptFe );
  
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
    
  private:
  
    //! Flag, if incompatible modes are used
    bool useICModes_;
  
    //! Cached entry for partial derivatives
    Matrix<Double> xiDx_;
    
    //! Cached vector for shape functions
    Vector<Double> shape_;
    
    //! Cached coordinates of integration point
    Vector<Double> ipCoord_;
  };
  
  template<class FE, class TYPE>
  void PiolaStressOperatorAxi<FE,TYPE>::CalcOpMat(Matrix<Double> & bMat,
                                                  const LocPointMapped& lpm,
                                                  BaseFE* ptFe )
  {
    // Get derivatives of local shape functions with respect to global
    // coords (format: nrNodes x spaceDim)
    FE *fe = (static_cast<FE*>(ptFe));
    fe->GetGlobDerivShFnc( xiDx_, lpm, lpm.shapeMap->GetElem(), 1 );
  
    // Calculate phi-phi component
    if (useICModes_) {
      fe->GetShFncICModes( shape_, lpm.lp, lpm.shapeMap->GetElem() );
    } else {
      fe->GetShFnc( shape_, lpm.lp, lpm.shapeMap->GetElem() );
    }
    
    // get radius of integration piint
    if (lpm.lp.coord.GetSize() > 0) {
      ipCoord_ = lpm.lp.coord;
    }
    else {
      lpm.shapeMap->Local2Global(ipCoord_, lpm.lp);
    }
    const Double oneOverR = 1.0 / ipCoord_[0];
    
    const UInt numFncs = xiDx_.GetNumRows();
    // Set correct size of matrix B and initialize with zeros
    bMat.Resize( DIM_D_MAT, numFncs * DIM_DOF );
    bMat.Init();
  
    UInt iFunc = 0;
    for (iFunc = 0; iFunc < numFncs; ++iFunc) {
      for (UInt iDim1 = 0; iDim1 < DIM_DOF; ++iDim1) {
        for (UInt iDim2 = 0; iDim2 < DIM_DOF; ++iDim2) {
          bMat[iDim1*DIM_DOF+iDim2][iFunc*DIM_DOF+iDim1] = xiDx_[iFunc][iDim2];
        }
      }
      bMat[4][iFunc*DIM_DOF] = shape_[iFunc] * oneOverR;
    }
  }
  
  template<class FE, class TYPE>
  void PiolaStressOperatorAxi<FE,TYPE>::CalcOpMatTransposed(Matrix<Double> & bMat,
                                                            const LocPointMapped& lpm,
                                                            BaseFE* ptFe )
  {
    Matrix<Double> tmpMat;
    this->CalcOpMat(tmpMat, lpm, ptFe);
    tmpMat.Transpose(bMat);
  }

  //! Calculate the operator for the computation of prestressing in 2.5D case
  //!    / N_1x  0    0   N_2x   0   ...\
  //! b =| N_1y  0    0   N_2y   0   ...|
  //!    |  0    0    0    0     0   ...|
  //!    |  0   N_1x  0    0    N_2x ...|
  //!    |  0   N_1y  0    0    N_2y ...|
  //!    |  0    0    0    0     0   ...|
  //!    |  0    0   N_1x  0     0   ...|
  //!    |  0    0   N_1y  0     0   ...|
  //!    \  0    0    0    0     0   .../
  //!  here N_1x denotes the x-derivative of the first
  //!  shape function at a given local point
  //!  this operator is only valid for vectorial unknowns
  //!  Valid for cartesian coordinate systems. Not sure about Axi...
  //! \tparam FE Type of Finite Element used
  //! \tparam D     Dimension of the problem space
  //! \tparam D_DOF Dimension of the unknown
  //! \tparam TYPE Data type (DOUBLE, COMPLEX)
  template<class FE, UInt D = 2, UInt D_DOF = 3, class TYPE = Double>
  class PreStressOperator2p5D : public BaseBOperator{

  public:
    // ------------------
    //  STATIC CONSTANTS
    // ------------------
    //@{
    //! \name Static constants

    //! Order of differentiation
    static const UInt ORDER_DIFF = 1;

    //! Number of components of the problem (scalar, vector)
    static const UInt DIM_DOF = D_DOF; // m=1 (Skalar), m=2,3 Vektor

    //! Dimension of the underlying domain / space
    static const UInt DIM_SPACE = D; // n=2,3

    //! Dimension of the finite element
    static const UInt DIM_ELEM = D; // Dimension von Referenzelement

    //! Dimension of the related material
    static const UInt DIM_D_MAT = D_DOF*D_DOF;

    //! Constructor
    //! \param useIcModes Use incompatible modes shape functions
    PreStressOperator2p5D( bool useICModes = false )
    :  useICModes_(useICModes) {
      this->name_ = "preStressOp2p5D";
    }

    //! Copy constructor
    PreStressOperator2p5D(const PreStressOperator2p5D & other)
       : BaseBOperator(other),
         useICModes_(other.useICModes_){
      this->name_ = other.name_;
    }

    //! \copydoc BaseBOperator::Clone()
    virtual PreStressOperator2p5D * Clone(){
      return new PreStressOperator2p5D(*this);
    }

    //! Destructor
    virtual ~PreStressOperator2p5D(){

    }

    //! Calculate operator matrix
    virtual void CalcOpMat(Matrix<Double> & bMat,
                           const LocPointMapped& lp,
                           BaseFE* ptFe );

    //! Calculate transposed operator matrix
    virtual void CalcOpMatTransposed(Matrix<Double> & bMat,
                                     const LocPointMapped& lp,
                                     BaseFE* ptFe );

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

    bool useICModes_;

  };

  template<class FE, UInt D, UInt D_DOF, class TYPE>
  void PreStressOperator2p5D<FE,D,D_DOF,TYPE>::CalcOpMat(Matrix<Double> & bMat,
                                               const LocPointMapped& lp,
                                               BaseFE* ptFe ){
    // Get derivatives of local shape functions with respect to global
    // coords (format: nrNodes x spaceDim)
    Matrix<Double> xiDx;
    FE *fe = (static_cast<FE*>(ptFe));
    if( useICModes_ ) {
      if ( isSurfOpt_ )
        fe->GetGlobDerivShFncICModes(  xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem() , 1 );
      else
        fe->GetGlobDerivShFncICModes( xiDx, lp, lp.shapeMap->GetElem() , 1 );

    } else {
      if ( isSurfOpt_ )
        fe->GetGlobDerivShFnc(  xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem() , 1 );
      else
        fe->GetGlobDerivShFnc( xiDx, lp, lp.shapeMap->GetElem() , 1 );
    }

    const UInt numFncs = xiDx.GetNumRows();
    // Set correct size of matrix B and initialize with zeros
    bMat.Resize( DIM_D_MAT, numFncs * DIM_DOF );
    bMat.Init();

    UInt iFunc = 0;
    for( iFunc = 0; iFunc < numFncs; iFunc++ ) {
      for(UInt iDim1 = 0; iDim1 < DIM_DOF; iDim1++ ) {
        for(UInt iDim2 = 0; iDim2 < DIM_SPACE; iDim2++ ) {
          bMat[iDim1*DIM_DOF+iDim2][iFunc*DIM_DOF+iDim1] = xiDx[iFunc][iDim2];
        }
      }
    }
  }

  template<class FE, UInt D, UInt D_DOF, class TYPE>
  void PreStressOperator2p5D<FE,D,D_DOF,TYPE>::CalcOpMatTransposed(Matrix<Double> & bMat,
                                                 const LocPointMapped& lp,
                                                 BaseFE* ptFe ){
    Matrix<Double> tmpMat;
    this->CalcOpMat(tmpMat,lp,ptFe);
    tmpMat.Transpose(bMat);
  }

  //! Calculate a scaled operator for the computation of prestressing
  //!    / sx*N_1x     0       0   sx*N_2x     0   ...\
  //! b =| sy*N_1y     0       0   sy*N_2y     0   ...|
  //!    | sz*N_1z     0       0   sz*N_2z     0   ...|
  //!    |     0   sx*N_1x     0      0    sx*N_2x ...|
  //!    |     0   sy*N_1y     0      0    sy*N_2y ...|
  //!    |     0   sz*N_1z     0      0    sz*N_2z ...|
  //!    |     0       0   sx*N_1x    0        0   ...|
  //!    |     0       0   sy*N_1y    0        0   ...|
  //!    \     0       0   sz*N_1z    0        0   .../
  //!  here N_1x denotes the x-derivative of the first
  //!  shape function at a given local point
  //!  this operator is only valid for vectorial unknowns
  //!  Valid for cartesian coordinate systems. Not sure about Axi...
  //! \tparam FE Type of Finite Element used
  //! \tparam D     Dimension of the problem space and unknown
  //! \tparam TYPE Data type (DOUBLE, COMPLEX)
  template<class FE, UInt D, class TYPE>
  class ScaledPreStressOperator : public PiolaStressOperator<FE, D, TYPE>
  {

  public:
    // ------------------
    //  STATIC CONSTANTS
    // ------------------
    //@{
    //! \name Static constants

    //! Order of differentiation
    static const UInt ORDER_DIFF = 1;

    //! Number of components of the problem (scalar, vector)
    static const UInt DIM_DOF = D; // m=1 (Skalar), m=2,3 Vektor

    //! Dimension of the underlying domain / space
    static const UInt DIM_SPACE = D; // n=2,3

    //! Dimension of the finite element
    static const UInt DIM_ELEM = D; // Dimension von Referenzelement

    //! Dimension of the related material
    static const UInt DIM_D_MAT = D*D;

    //! Constructor
    //! \param useIcModes Use incompatible modes shape functions
    ScaledPreStressOperator(bool useICModes = false) : PiolaStressOperator<FE, D, TYPE>(useICModes)
    {
      this->name_ = "ScaledPreStressOp";
    }

    //! Copy constructor
    ScaledPreStressOperator(const ScaledPreStressOperator & other) : PiolaStressOperator<FE, D, TYPE>(other)
    {
      this->name_ = other.name_;
    }

    //! \copydoc BaseBOperator::Clone()
    virtual ScaledPreStressOperator * Clone(){
      return new ScaledPreStressOperator(*this);
    }

    //! Destructor
    virtual ~ScaledPreStressOperator() { }

    //! Calculate operator matrix
    virtual void CalcOpMat(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe);

    //! Calculate operator matrix
    virtual void CalcOpMat(Matrix<Complex>& bMat, const LocPointMapped& lp, BaseFE* ptFe);

  };

  template<class FE, UInt D, class TYPE>
  void ScaledPreStressOperator<FE,D,TYPE>::CalcOpMat(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
  {
    assert(this->coef_ != NULL);
    Vector<Double> coefs;
    this->coef_->GetVector(coefs, lp);

    // Get derivatives of local shape functions with respect to global
    // coords (format: nrNodes x spaceDim)
    Matrix<Double> xiDx;
    FE *fe = (static_cast<FE*>(ptFe));
    if (this->useICModes_)
      if (this->isSurfOpt_)
        fe->GetGlobDerivShFncICModes(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
      else
        fe->GetGlobDerivShFncICModes(xiDx, lp, lp.shapeMap->GetElem(), 1);
    else
      if (this->isSurfOpt_)
        fe->GetGlobDerivShFnc(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
      else
        fe->GetGlobDerivShFnc(xiDx, lp, lp.shapeMap->GetElem(), 1);

    const UInt numFncs = xiDx.GetNumRows();
    // Set correct size of matrix B and initialize with zeros
    bMat.Resize(DIM_D_MAT, numFncs*DIM_DOF);
    bMat.Init();

    UInt iFunc = 0;
    for (iFunc = 0; iFunc < numFncs; iFunc++)
    {
      for (UInt iDim1 = 0; iDim1 < DIM_DOF; iDim1++)
      {
        for (UInt iDim2 = 0; iDim2 < DIM_DOF; iDim2++)
          bMat[iDim1*DIM_DOF + iDim2][iFunc*DIM_DOF + iDim1] = xiDx[iFunc][iDim2]*coefs[iDim2];
      }
    }
  }

  template<class FE, UInt D, class TYPE>
  void ScaledPreStressOperator<FE,D,TYPE>::CalcOpMat(Matrix<Complex>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
  {
    assert(this->coef_ != NULL);
    Vector<Complex> coefs;
    this->coef_->GetVector(coefs, lp);

    // Get derivatives of local shape functions with respect to global
    // coords (format: nrNodes x spaceDim)
    Matrix<Double> xiDx, xiDxTmp, rotMat;
    FE *fe = (static_cast<FE*>(ptFe));
    if (this->coef_->GetCoordinateSystem())
    {
      if (this->useICModes_)
        if (this->isSurfOpt_)
          fe->GetGlobDerivShFncICModes(xiDxTmp, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
        else
          fe->GetGlobDerivShFncICModes(xiDxTmp, lp, lp.shapeMap->GetElem(), 1);
      else
        if (this->isSurfOpt_)
          fe->GetGlobDerivShFnc(xiDxTmp, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
        else
          fe->GetGlobDerivShFnc(xiDxTmp, lp, lp.shapeMap->GetElem(), 1);

      // If coordinate system is set at the coefficient function, rotate B-matrix
      Vector<Double> globPoint;
      lp.shapeMap->Local2Global(globPoint, lp.lp.coord);
      this->coef_->GetCoordinateSystem()->GetGlobRotationMatrix(rotMat, globPoint);
      xiDx = xiDxTmp*rotMat;
    }
    else
    {
      if (this->useICModes_)
        if (this->isSurfOpt_)
          fe->GetGlobDerivShFncICModes(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
        else
          fe->GetGlobDerivShFncICModes(xiDx, lp, lp.shapeMap->GetElem(), 1);
      else
        if (this->isSurfOpt_)
          fe->GetGlobDerivShFnc(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
        else
          fe->GetGlobDerivShFnc(xiDx, lp, lp.shapeMap->GetElem(), 1);
    }

    const UInt numFncs = xiDx.GetNumRows();
    // Set correct size of matrix B and initialize with zeros
    bMat.Resize(DIM_D_MAT, numFncs*DIM_DOF);
    bMat.Init();

    UInt iFunc = 0;
    for (iFunc = 0; iFunc < numFncs; iFunc++)
    {
      for (UInt iDim1 = 0; iDim1 < DIM_DOF; iDim1++)
      {
        for (UInt iDim2 = 0; iDim2 < DIM_DOF; iDim2++)
          bMat[iDim1*DIM_DOF + iDim2][iFunc*DIM_DOF + iDim1] = xiDx[iFunc][iDim2]*coefs[iDim2];
      }
    }
  }

  //! Calculate a scaled operator for the computation of prestressing in 2.5D case
  //!    / sx*N_1x     0       0   sx*N_2x     0   ...\
  //! b =| sy*N_1y     0       0   sy*N_2y     0   ...|
  //!    |     0       0       0      0        0   ...|
  //!    |     0   sx*N_1x     0      0    sx*N_2x ...|
  //!    |     0   sy*N_1y     0      0    sy*N_2y ...|
  //!    |     0       0       0      0        0   ...|
  //!    |     0       0   sx*N_1x    0        0   ...|
  //!    |     0       0   sy*N_1y    0        0   ...|
  //!    \     0       0       0      0        0   .../
  //!  here N_1x denotes the x-derivative of the first
  //!  shape function at a given local point
  //!  this operator is only valid for vectorial unknowns
  //!  Valid for cartesian coordinate systems. Not sure about Axi...
  //! \tparam FE Type of Finite Element used
  //! \tparam D     Dimension of the problem space
  //! \tparam D_DOF Dimension of the unknown
  //! \tparam TYPE Data type (DOUBLE, COMPLEX)
  template<class FE, UInt D = 2, UInt D_DOF = 3, class TYPE = Double>
  class ScaledPreStressOperator2p5D : public PreStressOperator2p5D<FE, D, D_DOF, TYPE>
  {
  public:
    // ------------------
    //  STATIC CONSTANTS
    // ------------------
    //@{
    //! \name Static constants

    //! Order of differentiation
    static const UInt ORDER_DIFF = 1;

    //! Number of components of the problem (scalar, vector)
    static const UInt DIM_DOF = D_DOF; // m=1 (Skalar), m=2,3 Vektor

    //! Dimension of the underlying domain / space
    static const UInt DIM_SPACE = D; // n=2,3

    //! Dimension of the finite element
    static const UInt DIM_ELEM = D; // Dimension von Referenzelement

    //! Dimension of the related material
    static const UInt DIM_D_MAT = D_DOF*D_DOF;

    //! Constructor
    //! \param useIcModes Use incompatible modes shape functions
    ScaledPreStressOperator2p5D(bool useICModes = false)
        : PreStressOperator2p5D<FE, D, D_DOF, TYPE>(useICModes)
    {
      this->name_ = "ScaledPreStressOp2p5D";
    }

    //! Copy constructor
    ScaledPreStressOperator2p5D(const ScaledPreStressOperator2p5D& other)
        : PreStressOperator2p5D<FE, D, D_DOF, TYPE>(other)
    {
      this->name_ = other.name_;
    }

    //! \copydoc BaseBOperator::Clone()
    virtual ScaledPreStressOperator2p5D * Clone(){
      return new ScaledPreStressOperator2p5D(*this);
    }

    //! Destructor
    virtual ~ScaledPreStressOperator2p5D() { }

    //! Calculate operator matrix
    virtual void CalcOpMat(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe);

    //! Calculate operator matrix
    virtual void CalcOpMat(Matrix<Complex>& bMat, const LocPointMapped& lp, BaseFE* ptFe);

  };

  template<class FE, UInt D, UInt D_DOF, class TYPE>
  void ScaledPreStressOperator2p5D<FE, D, D_DOF, TYPE>::CalcOpMat(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
  {
    assert(this->coef_ != NULL);
    Vector<Double> coefs;
    this->coef_->GetVector(coefs, lp);

    // Get derivatives of local shape functions with respect to global
    // coords (format: nrNodes x spaceDim)
    Matrix<Double> xiDx;
    FE *fe = (static_cast<FE*>(ptFe));
    if (this->useICModes_)
      if (this->isSurfOpt_)
        fe->GetGlobDerivShFncICModes(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
      else
        fe->GetGlobDerivShFncICModes(xiDx, lp, lp.shapeMap->GetElem(), 1);
    else
      if (this->isSurfOpt_)
        fe->GetGlobDerivShFnc(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
      else
        fe->GetGlobDerivShFnc(xiDx, lp, lp.shapeMap->GetElem(), 1);

    const UInt numFncs = xiDx.GetNumRows();
    // Set correct size of matrix B and initialize with zeros
    bMat.Resize(DIM_D_MAT, numFncs*DIM_DOF);
    bMat.Init();

    UInt iFunc = 0;
    for (iFunc = 0; iFunc < numFncs; iFunc++)
    {
      for (UInt iDim1 = 0; iDim1 < DIM_DOF; iDim1++)
      {
        for (UInt iDim2 = 0; iDim2 < DIM_SPACE; iDim2++)
          bMat[iDim1*DIM_DOF + iDim2][iFunc*DIM_DOF + iDim1] = xiDx[iFunc][iDim2]*coefs[iDim2];
      }
    }
  }

  template<class FE, UInt D, UInt D_DOF, class TYPE>
  void ScaledPreStressOperator2p5D<FE, D, D_DOF, TYPE>::CalcOpMat(Matrix<Complex>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
  {
    assert(this->coef_ != NULL);
    Vector<Complex> coefs;
    this->coef_->GetVector(coefs, lp);

    // Get derivatives of local shape functions with respect to global
    // coords (format: nrNodes x spaceDim)
    Matrix<Double> xiDx, xiDxTmp, rotMat;
    FE *fe = (static_cast<FE*>(ptFe));
    if (this->coef_->GetCoordinateSystem())
    {
      if (this->useICModes_)
        if (this->isSurfOpt_)
          fe->GetGlobDerivShFncICModes(xiDxTmp, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
        else
          fe->GetGlobDerivShFncICModes(xiDxTmp, lp, lp.shapeMap->GetElem(), 1);
      else
        if (this->isSurfOpt_)
          fe->GetGlobDerivShFnc(xiDxTmp, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
        else
          fe->GetGlobDerivShFnc(xiDxTmp, lp, lp.shapeMap->GetElem(), 1);

      // If coordinate system is set at the coefficient function, rotate B-matrix
      Vector<Double> globPoint;
      lp.shapeMap->Local2Global(globPoint, lp.lp.coord);
      this->coef_->GetCoordinateSystem()->GetGlobRotationMatrix(rotMat, globPoint);
      xiDx = xiDxTmp*rotMat;
    }
    else
    {
      if (this->useICModes_)
        if (this->isSurfOpt_)
          fe->GetGlobDerivShFncICModes(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
        else
          fe->GetGlobDerivShFncICModes(xiDx, lp, lp.shapeMap->GetElem(), 1);
      else
        if (this->isSurfOpt_)
          fe->GetGlobDerivShFnc(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
        else
          fe->GetGlobDerivShFnc(xiDx, lp, lp.shapeMap->GetElem(), 1);
    }

    const UInt numFncs = xiDx.GetNumRows();
    // Set correct size of matrix B and initialize with zeros
    bMat.Resize(DIM_D_MAT, numFncs*DIM_DOF);
    bMat.Init();

    UInt iFunc = 0;
    for (iFunc = 0; iFunc < numFncs; iFunc++)
    {
      for (UInt iDim1 = 0; iDim1 < DIM_DOF; iDim1++)
      {
        for (UInt iDim2 = 0; iDim2 < DIM_SPACE; iDim2++)
          bMat[iDim1*DIM_DOF + iDim2][iFunc*DIM_DOF + iDim1] = xiDx[iFunc][iDim2]*coefs[iDim2];
      }
    }
  }


  template<class FE, UInt D, UInt D_DOF, class TYPE = Double>
  class SurfaceNormalPreStressOperator : public BaseBOperator
  {
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

    SurfaceNormalPreStressOperator(const std::string& subType, bool useICModes = false)
    {
      if(subType != "axi" )
      {
        if (subType == "2.5d")
          preStressOp_ = new PreStressOperator2p5D<FeH1, D, D_DOF, TYPE>(useICModes);
        else
          preStressOp_ = new PiolaStressOperator<FeH1, D, TYPE>(useICModes);
      }
      else
      {
        EXCEPTION("Unknown subtype '" << subType << "' in SurfaceNormalPreStressOperator");
      }

      this->name_ = "surfNormPreStressOp";
      preStressOp_->SetOperator2SurfOperator();
    }

    SurfaceNormalPreStressOperator(const std::string& subType, PtrCoefFct baseOpCoef, bool useICModes = false)
    {
      if(subType != "axi" )
      {
        if (subType == "2.5d")
          preStressOp_ = new ScaledPreStressOperator2p5D<FeH1, D, D_DOF, TYPE>(useICModes);
        else
          preStressOp_ = new ScaledPreStressOperator<FeH1, D, TYPE>(useICModes);
      }
      else
      {
        EXCEPTION("Unknown subtype '" << subType << "' in SurfaceNormalPreStressOperator");
      }

      this->name_ = "surfNormPreStressOp";
      preStressOp_->SetOperator2SurfOperator();
      preStressOp_->SetCoefFunction(baseOpCoef);
    }

    //! Copy constructor
    SurfaceNormalPreStressOperator(const SurfaceNormalPreStressOperator & other)
       : BaseBOperator(other) {
      this->name_ = other.name_;
      this->preStressOp_ = other.preStressOp_->Clone();
      this->preStressOp_->SetOperator2SurfOperator();
    }

    //! \copydoc BaseBOperator::Clone()
    virtual SurfaceNormalPreStressOperator * Clone(){
      return new SurfaceNormalPreStressOperator(*this);
    }

    virtual ~SurfaceNormalPreStressOperator() { return; }

    virtual void CalcOpMat(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe);

    virtual void CalcOpMatTransposed(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe);

    virtual void CalcOpMat(Matrix<Complex>& bMat, const LocPointMapped& lp, BaseFE* ptFe);

    virtual void CalcOpMatTransposed(Matrix<Complex>& bMat, const LocPointMapped& lp, BaseFE* ptFe);

    // ===============
    //  QUERY METHODS
    // ===============
    //@{ \name Query Methods
    //! \copydoc BaseBOperator::GetDiffOrder
    virtual UInt GetDiffOrder() const
    {
      return ORDER_DIFF;
    }

    //! \copydoc BaseBOperator::GetDimDof()
    virtual UInt GetDimDof() const
    {
      return DIM_DOF;
    }

    //! \copydoc BaseBOperator::GetDimSpace()
    virtual UInt GetDimSpace() const
    {
      return DIM_SPACE;
    }

    //! \copydoc BaseBOperator::GetDimElem()
    virtual UInt GetDimElem() const
    {
      return DIM_ELEM;
    }

    //! \copydoc BaseBOperator::GetDimDMat()
    virtual UInt GetDimDMat() const
    {
      return preStressOp_->GetDimDMat();
    }
    //@}

  protected:

    //! computes the normal operator
    virtual void ComputeNormalOperator(Matrix<Double>& normalOp, const LocPointMapped& lp)
    {
      UInt dimD = preStressOp_->GetDimDMat();
      normalOp.Resize(dimD, DIM_DOF);
      normalOp.Init();

      if (DIM_SPACE == 2)
      {
        //2D and 2.5D cases
        for (UInt i = 0; i < DIM_SPACE; ++i)
        {
          normalOp[i*DIM_SPACE][i] = lp.normal[0];
          normalOp[i*DIM_SPACE + 1][i] = lp.normal[1];
        }
      }
      else if (DIM_SPACE == 3)
      {
        // 3D case
        for (UInt i = 0; i < DIM_SPACE; ++i)
        {
          normalOp[i*DIM_SPACE][i] = lp.normal[0];
          normalOp[i*DIM_SPACE + 1][i] = lp.normal[1];
          normalOp[i*DIM_SPACE + 2][i] = lp.normal[2];
        }
      }
      else
      {
        EXCEPTION("Unknown formulation for 'DIM_SPACE' = " << DIM_SPACE);
      }
    }

    //! computes the normal operator
    virtual void ComputeNormalOperator(Matrix<Complex>& normalOp, const LocPointMapped& lp)
    {
      UInt dimD = preStressOp_->GetDimDMat();
      normalOp.Resize(dimD, DIM_DOF);
      normalOp.Init();

      Vector<Complex> coeffs;
      PtrCoefFct preStrOpCoef = preStressOp_->GetCoefFunction();
      if (preStrOpCoef == NULL)
      {
        coeffs.Resize(3, Complex(1.0, 0.0));
      }
      else
      {
        preStrOpCoef->GetVector(coeffs, lp);
      }

      if (DIM_SPACE == 2)
      {
        //2D and 2.5D cases
        for (UInt i = 0; i < DIM_SPACE; ++i)
        {
          normalOp[i*DIM_SPACE][i] = lp.normal[0]*coeffs[0];
          normalOp[i*DIM_SPACE + 1][i] = lp.normal[1]*coeffs[1];
        }
      }
      else if (DIM_SPACE == 3)
      {
        // 3D case
        for (UInt i = 0; i < DIM_SPACE; ++i)
        {
          normalOp[i*DIM_SPACE][i] = lp.normal[0]*coeffs[0];
          normalOp[i*DIM_SPACE + 1][i] = lp.normal[1]*coeffs[1];
          normalOp[i*DIM_SPACE + 2][i] = lp.normal[2]*coeffs[2];
        }
      }
      else
      {
        EXCEPTION("Unknown formulation for 'DIM_SPACE' = " << DIM_SPACE);
      }
    }

    //! prestressing operator
    BaseBOperator* preStressOp_;

  };

  template<class FE, UInt D, UInt D_DOF, class TYPE>
  void SurfaceNormalPreStressOperator<FE, D, D_DOF, TYPE>::CalcOpMat(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
  {
    //check if lp is surface and ptFe is volume
    assert(lp.isSurface);
    assert(D == ptFe->shape_.dim);

    Matrix<Double> tmpMat;
    this->CalcOpMatTransposed(tmpMat, lp, ptFe);
    bMat = Transpose(tmpMat);
  }

  template<class FE, UInt D, UInt D_DOF, class TYPE>
  void SurfaceNormalPreStressOperator<FE, D, D_DOF, TYPE>::CalcOpMatTransposed(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
  {
    //check if lp is surface and ptFe is volume
    assert(lp.isSurface);
    assert(D == ptFe->shape_.dim);

    UInt nFunc = ptFe->GetNumFncs();
    bMat.Resize(nFunc*D_DOF, D_DOF);
    bMat.Init();

    // get the tensor of the initial stress
    Matrix<Double> initStressTens;
//    this->coef_->GetTensor(initStressTens, lp);
    // TODO kirill:
    // see the implementation for Matrix<Complex>
    this->coef_->GetTensor(initStressTens, *(lp.lpmVol.get()));

    // get the prestressing B-operator matrix
    Matrix<Double> preStrOpMat;
    FE* fe = static_cast<FE*>(ptFe);
    preStressOp_->CalcOpMatTransposed(preStrOpMat, lp, fe);

    bMat = preStrOpMat*initStressTens;

    // get normal operator
    Matrix<Double> normOp;
    ComputeNormalOperator(normOp, lp);

    bMat *= normOp;
  }

  template<class FE, UInt D, UInt D_DOF, class TYPE>
  void SurfaceNormalPreStressOperator<FE, D, D_DOF, TYPE>::CalcOpMat(Matrix<Complex>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
  {
    //check if lp is surface and ptFe is volume
    assert(lp.isSurface);
    assert(D == ptFe->shape_.dim);

    Matrix<Complex> tmpMat;
    this->CalcOpMatTransposed(tmpMat, lp, ptFe);
    bMat = Transpose(tmpMat);
  }

  template<class FE, UInt D, UInt D_DOF, class TYPE>
  void SurfaceNormalPreStressOperator<FE, D, D_DOF, TYPE>::CalcOpMatTransposed(Matrix<Complex>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
  {
    //check if lp is surface and ptFe is volume
    assert(lp.isSurface);
    assert(D == ptFe->shape_.dim);

    UInt nFunc = ptFe->GetNumFncs();
    bMat.Resize(nFunc*D_DOF, D_DOF);
    bMat.Init();

    // get the tensor of the initial stress
    Matrix<Complex> initStressTens;
//    this->coef_->GetTensor(initStressTens, lp);
    // TODO kirill:
    // there're problems trying to get the value of the coefficient function
    // if it is a coefficient function compound and 'lp' is a local point of
    // an NC-surface-element; that's why the 'lp' of the corresponding volume
    // element is taken instead.
    this->coef_->GetTensor(initStressTens, *(lp.lpmVol.get()));

    // get the prestressing B-operator matrix
    Matrix<Complex> preStrOpMat;
    FE* fe = static_cast<FE*>(ptFe);
    preStressOp_->CalcOpMatTransposed(preStrOpMat, lp, fe);

    bMat = preStrOpMat*initStressTens;

    // get normal operator
    Matrix<Complex> normOp;
    ComputeNormalOperator(normOp, lp);

    bMat *= normOp;
  }
}
#endif /* PRESTRESSOPERATOR_HH_ */
