// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     NormalFluxOperators.hh
 *       \brief    All operators for computing fluxes in normal direction
 *
 *       \date     Feb 22, 2012
 *       \author   ahueppe
 */
//================================================================================================

#ifndef NORMALFLUXOPERATORS_HH_
#define NORMALFLUXOPERATORS_HH_

#include "IdentityOperator.hh"

namespace CoupledField{


  template<class FE, UInt D = 1, UInt D_DOF = 1, class TYPE = Double>
  class NormalFluxOperator : public IdentityOperator<FE,D,D_DOF,TYPE>{
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

    NormalFluxOperator(){
      return;
    }

    //! Copy constructor
    NormalFluxOperator(const NormalFluxOperator & other)
       : IdentityOperator<FE,D,D_DOF,TYPE>(other){
    }

    //! \copydoc BaseBOperator::Clone()
    virtual NormalFluxOperator * Clone(){
      return new NormalFluxOperator(*this);
    }

    virtual ~NormalFluxOperator(){
      return;
    }

    virtual void CalcOpMat(Matrix<Double> & bMat,
                           const LocPointMapped& lp, BaseFE* ptFe ){
      assert(this->coef_ != NULL);
      if(!lp.isSurface){
        EXCEPTION("NormalFluxOperator only suited for surface elements");
      }
      Vector<Double> myVec;
      this->coef_->GetVector(myVec,lp);
      Double normalScal = myVec * lp.normal;
      IdentityOperator<FE,D,D_DOF,TYPE>::CalcOpMat(bMat,lp,ptFe);
      bMat *= normalScal;

    }

    virtual void CalcOpMatTransposed(Matrix<Double> & bMat,
                                     const LocPointMapped& lp, BaseFE* ptFe ){
      assert(this->coef_ != NULL);
      if(!lp.isSurface){
        EXCEPTION("NormalFluxOperator only suited for surface elements");
      }
      Vector<Double> myVec;
      this->coef_->GetVector(myVec,lp);
      Double normalScal = myVec * lp.normal;
      IdentityOperator<FE,D,D_DOF,TYPE>::CalcOpMatTransposed(bMat,lp,ptFe);
      bMat *= normalScal;
    }

    //avoid reimplementation of complex operator by making the bas class function
    //available
    using IdentityOperator<FE,D,D_DOF,TYPE>::CalcOpMat;

    using IdentityOperator<FE,D,D_DOF,TYPE>::CalcOpMatTransposed;
    
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
  };

  template<class FE, UInt D = 1, UInt D_DOF = 1, class TYPE = Double>
  class NormalFluxOperatorPiola : public IdentityOperatorPiola<FE,D,D_DOF,TYPE>{
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

    NormalFluxOperatorPiola(){
      return;
    }

    //! Copy constructor
    NormalFluxOperatorPiola(const NormalFluxOperatorPiola & other)
       : IdentityOperatorPiola<FE,D,D_DOF,TYPE>(other){
    }

    //! \copydoc BaseBOperator::Clone()
    virtual NormalFluxOperatorPiola * Clone(){
      return new NormalFluxOperatorPiola(*this);
    }

    virtual ~NormalFluxOperatorPiola(){
      return;
    }

    virtual void CalcOpMat(Matrix<Double> & bMat,
                           const LocPointMapped& lp, BaseFE* ptFe ){
      assert(this->coef_ != NULL);
      if(!lp.isSurface){
        EXCEPTION("NormalFluxOperator only suited for surface elements");
      }
      Vector<Double> myVec;
      this->coef_->GetVector(myVec,lp);
      Double normalScal = myVec * lp.normal;
      IdentityOperatorPiola<FE,D,D_DOF,TYPE>::CalcOpMat(bMat,lp,ptFe);
      bMat *= normalScal;
    }

    virtual void CalcOpMatTransposed(Matrix<Double> & bMat,
                                     const LocPointMapped& lp, BaseFE* ptFe ){
      assert(this->coef_ != NULL);
      if(!lp.isSurface){
        EXCEPTION("NormalFluxOperator only suited for surface elements");
      }
      Vector<Double> myVec;
      this->coef_->GetVector(myVec,lp);
      Double normalScal = myVec * lp.normal;
      IdentityOperatorPiola<FE,D,D_DOF,TYPE>::CalcOpMatTransposed(bMat,lp,ptFe);
      bMat *= normalScal;
    }

    //avoid reimplementation of complex operator by making the bas class function
    //available
    using IdentityOperatorPiola<FE,D,D_DOF,TYPE>::CalcOpMat;

    using IdentityOperatorPiola<FE,D,D_DOF,TYPE>::CalcOpMatTransposed;

  };

}

#endif /* NORMALFLUXOPERATORS_HH_ */
