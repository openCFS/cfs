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
  //! Calculate the operator for the computation of prestressing
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
  //!  Valid for cartesian coordinate systems. Not sure about Axi...
  //! \tparam FE Type of Finite Element used
  //! \tparam D     Dimension of the problem space and unknown
  //! \tparam TYPE Data type (DOUBLE, COMPLEX)
  template<class FE, UInt D, class TYPE = Double>
  class PreStressOperator : public BaseBOperator{

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
    PreStressOperator( bool useICModes = false )
    :  useICModes_(useICModes) {
      this->name_ = "preStressOp";
    }
    //! Destructor
    virtual ~PreStressOperator(){

    }

    //! Calculate operator matrix
    virtual void CalcOpMat(Matrix<Double> & bMat,
                           const LocPointMapped& lp,
                           BaseFE* ptFe );

    //! Calculate transposed operator matrix
    virtual void CalcOpMatTransposed(Matrix<Double> & bMat,
                                     const LocPointMapped& lp,
                                     BaseFE* ptFe );

    //! Calculate operator matrix
    virtual void CalcOpMat(Matrix<Complex> & bMat,
                           const LocPointMapped& lp,
                           BaseFE* ptFe );

    //! Calculate transposed operator matrix
    virtual void CalcOpMatTransposed(Matrix<Complex> & bMat,
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
  private:

    bool useICModes_;

  };

  template<class FE, UInt D, class TYPE>
  void PreStressOperator<FE,D,TYPE>::CalcOpMat(Matrix<Double> & bMat,
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
  void PreStressOperator<FE,D,TYPE>::CalcOpMatTransposed(Matrix<Double> & bMat,
                                                 const LocPointMapped& lp,
                                                 BaseFE* ptFe ){
    Matrix<Double> tmpMat;
    this->CalcOpMat(tmpMat,lp,ptFe);
    tmpMat.Transpose(bMat);
  }

  template<class FE, UInt D, class TYPE>
  void PreStressOperator<FE,D,TYPE>::CalcOpMat(Matrix<Complex> & bMat,
                                               const LocPointMapped& lp,
                                               BaseFE* ptFe ){
    if (coef_ == NULL)
    {
      BaseBOperator::CalcOpMat(bMat, lp, ptFe);
    }
    else
    {
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
      Vector<TYPE> coeffs;
      coef_->GetVector(coeffs, lp);

      const UInt numFncs = xiDx.GetNumRows();
      // Set correct size of matrix B and initialize with zeros
      bMat.Resize( DIM_D_MAT, numFncs * DIM_DOF );
      bMat.Init();

      UInt iFunc = 0;
      for( iFunc = 0; iFunc < numFncs; iFunc++ ) {
        for(UInt iDim1 = 0; iDim1 < DIM_DOF; iDim1++ ) {
          for(UInt iDim2 = 0; iDim2 < DIM_DOF; iDim2++ ) {
            bMat[iDim1*DIM_DOF+iDim2][iFunc*DIM_DOF+iDim1] = xiDx[iFunc][iDim2]*coeffs[iDim2];
          }
        }
      }
    }
  }

  template<class FE, UInt D, class TYPE>
  void PreStressOperator<FE,D,TYPE>::CalcOpMatTransposed(Matrix<Complex> & bMat,
                                                 const LocPointMapped& lp,
                                                 BaseFE* ptFe ){
    Matrix<Complex> tmpMat;
    this->CalcOpMat(tmpMat,lp,ptFe);
    tmpMat.Transpose(bMat);
  }

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

    //! Calculate operator matrix
    virtual void CalcOpMat(Matrix<Complex> & bMat,
                           const LocPointMapped& lp,
                           BaseFE* ptFe );

    //! Calculate transposed operator matrix
    virtual void CalcOpMatTransposed(Matrix<Complex> & bMat,
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
  private:

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

  template<class FE, UInt D, UInt D_DOF, class TYPE>
  void PreStressOperator2p5D<FE,D,D_DOF,TYPE>::CalcOpMat(Matrix<Complex> & bMat,
                                               const LocPointMapped& lp,
                                               BaseFE* ptFe ){
    if (coef_ == NULL)
    {
      BaseBOperator::CalcOpMat(bMat, lp, ptFe);
    }
    else
    {
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
      Vector<TYPE> coeffs;
      coef_->GetVector(coeffs, lp);

      const UInt numFncs = xiDx.GetNumRows();
      // Set correct size of matrix B and initialize with zeros
      bMat.Resize( DIM_D_MAT, numFncs * DIM_DOF );
      bMat.Init();

      UInt iFunc = 0;
      for( iFunc = 0; iFunc < numFncs; iFunc++ ) {
        for(UInt iDim1 = 0; iDim1 < DIM_DOF; iDim1++ ) {
          for(UInt iDim2 = 0; iDim2 < DIM_SPACE; iDim2++ ) {
            bMat[iDim1*DIM_DOF+iDim2][iFunc*DIM_DOF+iDim1] = xiDx[iFunc][iDim2]*coeffs[iDim2];
          }
        }
      }
    }
  }

  template<class FE, UInt D, UInt D_DOF, class TYPE>
  void PreStressOperator2p5D<FE, D, D_DOF, TYPE>::CalcOpMatTransposed(Matrix<Complex> & bMat,
                                                 const LocPointMapped& lp,
                                                 BaseFE* ptFe ){
    Matrix<Complex> tmpMat;
    this->CalcOpMat(tmpMat,lp,ptFe);
    tmpMat.Transpose(bMat);
  }


  template<class FE, UInt D = 1, UInt D_DOF = 1, class TYPE = Double>
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

    SurfaceNormalPreStressOperator(std::string subType, bool useICModes = false)
    {
      if(subType != "axi" )
      {
        if (subType == "2.5d")
          preStressOp_ = new PreStressOperator2p5D<FeH1, D, D_DOF, TYPE>();
        else
          preStressOp_ = new PreStressOperator<FeH1, D, TYPE>();
      }
      else
      {
        EXCEPTION("Unknown subtype '" << subType << "' in SurfaceNormalPreStressOperator");
      }

      this->name_ = "surfNormPreStressOp";
      preStressOp_->SetOperator2SurfOperator();
      useICModes_ = useICModes;
    }

    SurfaceNormalPreStressOperator(std::string subType, PtrCoefFct baseOpCoef, bool useICModes = false)
    {
      if(subType != "axi" )
      {
        if (subType == "2.5d")
          preStressOp_ = new PreStressOperator2p5D<FeH1, D, D_DOF, TYPE>();
        else
          preStressOp_ = new PreStressOperator<FeH1, D, TYPE>();
      }
      else
      {
        EXCEPTION("Unknown subtype '" << subType << "' in SurfaceNormalPreStressOperator");
      }

      this->name_ = "surfNormPreStressOp";
      preStressOp_->SetOperator2SurfOperator();
      useICModes_ = useICModes;
      preStressOp_->SetCoefFunction(baseOpCoef);
    }

    virtual ~SurfaceNormalPreStressOperator() { return; }

    virtual void CalcOpMat(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe);

    virtual void CalcOpMatTransposed(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe);

    virtual void CalcOpMat(Matrix<Complex>& bMat, const LocPointMapped& lp, BaseFE* ptFe);

    virtual void CalcOpMatTransposed(Matrix<Complex>& bMat, const LocPointMapped& lp, BaseFE* ptFe);

    //avoid reimplementation of complex operator by making the base class function
    //available
    using BaseBOperator::CalcOpMat;

    using BaseBOperator::CalcOpMatTransposed;

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
      PtrCoefFct gradOpCoef = preStressOp_->GetCoefFunction();
      if (gradOpCoef == NULL)
      {
        coeffs.Resize(3, Complex(1.0, 0.0));
      }
      else
      {
        gradOpCoef->GetVector(coeffs, lp);
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

    //! Flag, if incompatible modes are used
    bool useICModes_;
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
