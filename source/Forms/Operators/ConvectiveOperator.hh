// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*- 
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     ConvectiveOperator.hh
 *       \brief    <ENTER BRIEF DESCRIPTION HERE>
 *
 *       \date     Feb 12, 2012
 *       \author   ahueppe
 */ 
//================================================================================================

#ifndef CONVECTIVEOPERATOR_HH_
#define CONVECTIVEOPERATOR_HH_


#include "BaseBOperator.hh"

namespace CoupledField{

  template<class FE, UInt D, UInt D_DOF, class TYPE = Double>
  class ConvectiveOperator : public BaseBOperator<TYPE>{
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
    static const UInt DIM_D_MAT = D;
    //@}

    ConvectiveOperator(){
      this->name_ = "DivOperator";
    }

    virtual ~ConvectiveOperator(){

    }

    virtual void CalcOpMat(Matrix<Double> & bMat,
                           const LocPointMapped& lp,
                           BaseFE* ptFe );

    virtual void CalcOpMatTransposed(Matrix<Double> & bMat,
                                     const LocPointMapped& lp,
                                     BaseFE* ptFe );

    virtual void CalcOpMat(Matrix<Complex> & bMat,
                           const LocPointMapped& lp,
                           BaseFE* ptFe ){
      EXCEPTION("ConvectiveOperator::CalcOpMat not implement for Complex");
    }

    virtual void CalcOpMatTransposed(Matrix<Complex> & bMat,
                                     const LocPointMapped& lp,
                                     BaseFE* ptFe ){
      EXCEPTION("ConvectiveOperator::CalcOpMatTransposed not implement for Complex");
    }


    protected:

  };


  template<class FE, UInt D, UInt D_DOF, class TYPE>
  void ConvectiveOperator<FE,D,D_DOF,TYPE>::CalcOpMat(Matrix<Double> & bMat,
                                         const LocPointMapped& lp,
                                         BaseFE* ptFe ){

    //obtain external field
    Vector<TYPE> myVec;
    this->coef_->GetVector(myVec,lp);

    UInt numFncs = ptFe->GetNumFncs();
    // Set correct size of matrix B and initialize with zeros
    bMat.Resize( DIM_DOF, numFncs * DIM_DOF );
    bMat.Init();

    // Get derivatives of local shape functions with respect to global
    // coords (format: nrNodes x spaceDim)
    Matrix<Double> xiDx;
    FE *fe = (static_cast<FE*>(ptFe));
    fe->GetGlobDerivShFnc( xiDx, lp, lp.shapeMap->GetElem() , 1 );
    for(UInt iDimDof = 0; iDimDof < DIM_DOF; ++iDimDof) {
      for(UInt iDim = 0; iDim < DIM_SPACE; ++iDim) {
        for( UInt i = 0; i < numFncs; ++i ) {
          bMat[iDimDof][i*DIM_DOF + iDimDof] += xiDx[i][iDim] * myVec[iDim];
        }
      }
    }
  }

  template<class FE, UInt D, UInt D_DOF, class TYPE>
  void ConvectiveOperator<FE,D,D_DOF,TYPE>::CalcOpMatTransposed(Matrix<Double> & bMat,
                                                   const LocPointMapped& lp,
                                                   BaseFE* ptFe ){
    //obtain external field
    Vector<TYPE> myVec;
    this->coef_->GetVector(myVec,lp);

    UInt numFncs = ptFe->GetNumFncs();
    // Set correct size of matrix B and initialize with zeros
    bMat.Resize( numFncs * DIM_DOF , DIM_DOF);
    bMat.Init();

    // Get derivatives of local shape functions with respect to global
    // coords (format: nrNodes x spaceDim)
    Matrix<Double> xiDx;
    FE *fe = (static_cast<FE*>(ptFe));
    fe->GetGlobDerivShFnc( xiDx, lp, lp.shapeMap->GetElem() , 1 );
    for(UInt iDimDof = 0; iDimDof < DIM_DOF; ++iDimDof) {
      for(UInt iDim = 0; iDim < DIM_SPACE; ++iDim) {
        for( UInt i = 0; i < numFncs; ++i ) {
          bMat[i*DIM_DOF + iDimDof][iDimDof] += xiDx[i][iDim] * myVec[iDim];
        }
      }
    }
  }

  template<class FE, UInt D, UInt D_DOF, class TYPE = Double>
   class ConvectiveOperatorPiola : public ConvectiveOperator<FE,D,D_DOF,TYPE>{
     public:

     ConvectiveOperatorPiola(){
       this->name_ = "DivOperator";
     }

     virtual ~ConvectiveOperatorPiola(){

     }

     virtual void CalcOpMat(Matrix<Double> & bMat,
                            const LocPointMapped& lp,
                            BaseFE* ptFe );

     virtual void CalcOpMatTransposed(Matrix<Double> & bMat,
                                      const LocPointMapped& lp,
                                      BaseFE* ptFe );

     virtual void CalcOpMat(Matrix<Complex> & bMat,
                            const LocPointMapped& lp,
                            BaseFE* ptFe ){
       EXCEPTION("ConvectiveOperator::CalcOpMat not implement for Complex");
     }

     virtual void CalcOpMatTransposed(Matrix<Complex> & bMat,
                                      const LocPointMapped& lp,
                                      BaseFE* ptFe ){
       EXCEPTION("ConvectiveOperator::CalcOpMatTransposed not implement for Complex");
     }


     protected:

   };


   template<class FE, UInt D, UInt D_DOF, class TYPE>
   void ConvectiveOperatorPiola<FE,D,D_DOF,TYPE>::CalcOpMat(Matrix<Double> & bMat,
                                          const LocPointMapped& lp,
                                          BaseFE* ptFe ){

     Matrix<Double> bMatTmp;
     ConvectiveOperator<FE,D,D_DOF,TYPE>::CalcOpMat(bMatTmp,lp,ptFe);

     //now apply piola transform
     bMat = lp.jac * bMatTmp;
     bMat /= lp.jacDet;
     //std::cout << lp.jac << std::endl;
     //std::cout << bMatTmp << std::endl;

   }

   template<class FE, UInt D, UInt D_DOF, class TYPE>
   void ConvectiveOperatorPiola<FE,D,D_DOF,TYPE>::CalcOpMatTransposed(Matrix<Double> & bMat,
                                                    const LocPointMapped& lp,
                                                    BaseFE* ptFe ){

     Matrix<Double> bMatTmp;
     ConvectiveOperator<FE,D,D_DOF,TYPE>::CalcOpMatTransposed(bMatTmp,lp,ptFe);
     //now apply piola transform
     bMat =  bMatTmp * lp.jac;
     bMat /= lp.jacDet;
   }

}

#endif /* CONVECTIVEOPERATOR_HH_ */
