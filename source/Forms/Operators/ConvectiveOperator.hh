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
  class ConvectiveOperator : public BaseBOperator{
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

    ConvectiveOperator(){
      this->name_ = "ConvectiveOperator";
    }

    ConvectiveOperator(const ConvectiveOperator & other)
     : BaseBOperator(other){
    }

    virtual ConvectiveOperator * Clone(){
      return new ConvectiveOperator(*this);
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
                           BaseFE* ptFe );
    

    virtual void CalcOpMatTransposed(Matrix<Complex> & bMat,
                                     const LocPointMapped& lp,
                                     BaseFE* ptFe );
    

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
  void ConvectiveOperator<FE,D,D_DOF,TYPE>::CalcOpMat(Matrix<Double> & bMat,
                                         const LocPointMapped& lp,
                                         BaseFE* ptFe ){

    //obtain external field
    Vector<Double> myVec;
    this->coef_->GetVector(myVec,lp);

    //std::cout << "Velocity at IP" << std::endl << myVec << std::endl;

    const UInt numFncs = ptFe->GetNumFncs();
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
    Vector<Double> myVec;
    this->coef_->GetVector(myVec,lp);

    const UInt numFncs = ptFe->GetNumFncs();
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

  template<class FE, UInt D, UInt D_DOF, class TYPE>
  void ConvectiveOperator<FE,D,D_DOF,TYPE>::CalcOpMat(Matrix<Complex> & bMat,
                                         const LocPointMapped& lp,
                                         BaseFE* ptFe ){

    //obtain external field
	Vector<Complex> myVec;
    this->coef_->GetVector(myVec,lp);

    //std::cout << "Velocity at IP" << std::endl << myVec << std::endl;

    const UInt numFncs = ptFe->GetNumFncs();
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
  void ConvectiveOperator<FE,D,D_DOF,TYPE>::CalcOpMatTransposed(Matrix<Complex> & bMat,
                                                   const LocPointMapped& lp,
                                                   BaseFE* ptFe ){
    //obtain external field
    Vector<Complex> myVec;
    this->coef_->GetVector(myVec,lp);

    const UInt numFncs = ptFe->GetNumFncs();
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
       this->name_ = "ConvectiveOperatorPiola";
     }

     ConvectiveOperatorPiola(const ConvectiveOperatorPiola & other)
      : ConvectiveOperator<FE,D,D_DOF,TYPE>(other){
     }

     virtual ConvectiveOperatorPiola * Clone(){
       return new ConvectiveOperatorPiola(*this);
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

     Matrix<Double> bMatInitial;
     ConvectiveOperator<FE,D,D_DOF,TYPE>::CalcOpMat(bMatInitial,lp,ptFe);

     bMat.Resize(bMatInitial.GetNumRows(),bMatInitial.GetNumCols());
     bMat.Init();
     //now apply piola transform
     //in case of NC_SURF_ELEMs we evaluate the piola matrix on the volume element

     if(lp.isSurface){
#ifdef NDEBUG
       Double jacDetInv = (1.0/lp.lpmVol->jacDet);
       lp.lpmVol->jac.Mult_Blas(bMatInitial,bMat,false,false,jacDetInv,0.0);
#else
       bMat = lp.lpmVol->jac * bMatInitial;
       bMat *= (1.0/lp.lpmVol->jacDet);
#endif
     }else{
#ifdef NDEBUG
       Double jacDetInv = (1.0/lp.jacDet);
       lp.jac.Mult_Blas(bMatInitial,bMat,false,false,jacDetInv,0.0);
#else
       bMat = lp.jac * bMatInitial;
       bMat *= (1.0/lp.jacDet);
#endif
     }


   }

   template<class FE, UInt D, UInt D_DOF, class TYPE>
   void ConvectiveOperatorPiola<FE,D,D_DOF,TYPE>::CalcOpMatTransposed(Matrix<Double> & bMat,
                                                    const LocPointMapped& lp,
                                                    BaseFE* ptFe ){

     Matrix<Double> bMatInitial;
     ConvectiveOperator<FE,D,D_DOF,TYPE>::CalcOpMatTransposed(bMatInitial,lp,ptFe);
     bMat.Resize(bMatInitial.GetNumRows(),bMatInitial.GetNumCols());
     bMat.Init();
     //now apply piola transform
     //in case of NC_SURF_ELEMs we evaluate the piola matrix on the volume element

     if(lp.isSurface){
#ifdef NDEBUG
       Double jacDetInv = (1.0/lp.lpmVol->jacDet);
       bMatInitial.Mult_Blas(lp.lpmVol->jac,bMat,false,true,jacDetInv,0.0);
#else
       Matrix<Double> jacTmp;
       lp.lpmVol->jac.Transpose(jacTmp);
       bMat = bMatInitial * jacTmp;
       bMat *= (1.0/lp.lpmVol->jacDet);
#endif
     }else{
#ifdef NDEBUG
       Double jacDetInv = (1.0/lp.jacDet);
       bMatInitial.Mult_Blas(lp.jac,bMat,false,true,jacDetInv,0.0);
#else
       Matrix<Double> jacTmp;
       lp.jac.Transpose(jacTmp);
       bMat = bMatInitial * jacTmp;
       bMat *= (1.0/lp.jacDet);
#endif
     }
   }

}

#endif /* CONVECTIVEOPERATOR_HH_ */
