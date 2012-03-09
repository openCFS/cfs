#ifndef DIV_OPERATOR_HH
#define DIV_OPERATOR_HH

#include "BaseBOperator.hh"


namespace CoupledField{

//! Calculate the divergence of the vector shape functions
//!    / N_1x  0        ...\
//! b =| 0    N_1y      ...|
//!    \ 0     0   N_1z .../
//!  here N_1x denotes the x-derivative of the first
//!  shape function at a given local point
//! \tparam FE Type of Finite Element used
//! \tparam D Dimension of the problem space
//! \tparam TYPE Data type (DOUBLE, COMPLEX)
template<class FE, UInt D, class TYPE = Double>
class DivOperator : public BaseBOperator<TYPE>{
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


    DivOperator(){
      this->name_ = "DivOperator";
    }

    ~DivOperator(){

    }

    virtual void CalcOpMat(Matrix<Double> & bMat,
                           const LocPointMapped& lp,
                           BaseFE* ptFe );

    virtual void CalcOpMatTransposed(Matrix<Double> & bMat,
                                     const LocPointMapped& lp, 
                                     BaseFE* ptFe );

    using BaseBOperator<TYPE>::CalcOpMat;

    using BaseBOperator<TYPE>::CalcOpMatTransposed;


    protected:

  };

  template<class FE, UInt D, class TYPE>
  void DivOperator<FE,D,TYPE>::CalcOpMat(Matrix<Double> & bMat,
                                         const LocPointMapped& lp, 
                                         BaseFE* ptFe ){
    UInt numFncs = ptFe->GetNumFncs();
    // Set correct size of matrix B and initialize with zeros
    bMat.Resize( DIM_D_MAT, numFncs * DIM_SPACE );
    bMat.Init();
    
    // Get derivatives of local shape functions with respect to global
    // coords (format: nrNodes x spaceDim)
    Matrix<Double> xiDx;
    FE *fe = (static_cast<FE*>(ptFe));
    fe->GetGlobDerivShFnc( xiDx, lp, lp.shapeMap->GetElem() , 1 );
    for(UInt iDim = 0; iDim < DIM_SPACE; ++iDim) {
      for( UInt i = 0; i < numFncs; ++i ) {
        bMat[iDim][i*DIM_SPACE + iDim] = xiDx[i][iDim];
      }
    }
  }

  template<class FE, UInt D, class TYPE>
  void DivOperator<FE,D,TYPE>::CalcOpMatTransposed(Matrix<Double> & bMat,
                                                   const LocPointMapped& lp, 
                                                   BaseFE* ptFe ){
    UInt numFncs = ptFe->GetNumFncs();
    // Set correct size of matrix B and initialize with zeros
    bMat.Resize( numFncs * DIM_SPACE, DIM_D_MAT );
    bMat.Init();
    
    // Get derivatives of local shape functions with respect to global
    // coords (format: nrNodes x spaceDim)
    Matrix<Double> xiDx;
    FE *fe = (static_cast<FE*>(ptFe));
    fe->GetGlobDerivShFnc( xiDx, lp, lp.shapeMap->GetElem() , 1 );
    for( UInt i = 0; i < numFncs; ++i ) {
      for(UInt iDim = 0; iDim < DIM_SPACE; ++iDim ) {
        bMat[i+iDim][iDim] = xiDx[i][iDim];
      }
    }
  }

} // namespace
#endif
