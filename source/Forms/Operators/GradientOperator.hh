#ifndef GRADIENTOP_HH
#define GRADIENTOP_HH

#include "BaseBOperator.hh"

//! Calculate the gradient of the shape functions
//!    / N_1z N_2z ...\
//! b =| N_1y N_2y ...|
//!    \ N_1z N_2z .../
//!  here N_1x denotes the x-derivative of the first
//!  shape function at a given local point
//! \tparam FE Type of Finite Element used
//! \tparam D Dimension of the problem space
//! \tparam TYPE Data type (DOUBLE, COMPLEX)
namespace CoupledField{
  template<class FE, UInt D, class TYPE = Double>
  class GradientOperator : public BaseBOperator<FE,TYPE>{
    public:

    // ------------------
    //  STATIC CONSTANTS 
    // ------------------
    //@{ 
    //! \name Static constants

    //! Order of differentiation
    static const UInt ORDER_DIFF = 1;

    //! Number of components of the problem (scalar, vector)
    static const UInt DIM_DOF = 1;

    //! Dimension of the underlying domain / space
    static const UInt DIM_SPACE = D;

    //! Dimension of the finite element
    static const UInt DIM_ELEM = D;
    
    //! Dimension of the related material 
    static const UInt DIM_D_MAT = D; 
    //@}


    GradientOperator(){
      this->name_ = "GradientOperator";
    }

    ~GradientOperator(){

    }

    virtual void CalcOpMat(Matrix<Double> & bMat,
                           const LocPointMapped& lp,
                           BaseFE* ptFe );

    virtual void CalcOpMatTransposed(Matrix<Double> & bMat,
                                     const LocPointMapped& lp, 
                                     BaseFE* ptFe );

    using BaseBOperator<FE,TYPE>::CalcOpMat;

    using BaseBOperator<FE,TYPE>::CalcOpMatTransposed;


    protected:

  };

  template<class FE, UInt D, class TYPE>
  void GradientOperator<FE,D,TYPE>::CalcOpMat(Matrix<Double> & bMat,
                                                const LocPointMapped& lp, 
                                                BaseFE* ptFe ){
    UInt numFncs = ptFe->GetNumFncs();
    // Set correct size of matrix B and initialise with zeros
    UInt eDim = ptFe->shape_.dim;
    bMat.Resize( eDim, numFncs );

    // Get derivatives of local shape functions with respect to global
    // coords (format: nrNodes x spaceDim)
    Matrix<Double> xiDx;
    FE *fe = (static_cast<FE*>(ptFe));
    fe->GetGlobDerivShFnc( xiDx, lp, lp.shapeMap->GetElem() , 1 );
    bMat= Transpose(xiDx);
  }

  template<class FE, UInt D, class TYPE>
  void GradientOperator<FE,D,TYPE>::CalcOpMatTransposed(Matrix<Double> & bMat,
                                                          const LocPointMapped& lp, 
                                                          BaseFE* ptFe ){
    const UInt numFncs = ptFe->GetNumFncs();
    // Set correct size of matrix B and initialise with zeros
    const UInt eDim = ptFe->shape_.dim;
    bMat.Resize(numFncs , eDim );

    // Get derivatives of local shape functions with respect to global
    // coords (format: spaceDim x nrNodes )
    FE *fe = (static_cast<FE*>(ptFe));
    fe->GetGlobDerivShFnc( bMat, lp, lp.shapeMap->GetElem() , 1 );


  }
}
#endif
