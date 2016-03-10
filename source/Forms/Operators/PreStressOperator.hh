// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     PreStressOperator.hh
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

}
#endif /* PRESTRESSOPERATOR_HH_ */
