#ifndef LAPLOP_HH
#define LAPLOP_HH

#include "BaseBOperator.hh"

namespace CoupledField{
  //! Calculate the gradient of the shape functions
  //!    / N_1x N_2x ...\
  //! b =| N_1y N_2y ...|
  //!    \ N_1z N_2z .../
  //!  here N_1x denotes the x-derivative of the first
  //!  shape function at a given local point
  //! \tparam FE Type of Finite Element used
  //! \tparam D Dimension of the problem space
  //! \tparam TYPE Data type (DOUBLE, COMPLEX)
  template<class FE, UInt D, class TYPE = Double>
  class LaplOperator : public BaseBOperator<TYPE>{
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


    LaplOperator(){
      this->name_ = "LaplOperator";
    }

    virtual ~LaplOperator(){

    }

    //!
    //! Provide detailed desciption of this function
    //! @param parmeter1 Describe this parameter
    //!
    //! Here is an example of inserting a mathematical formula into the text:
    //! The distance is computed as /f$\sqrt{ (x_2-x_1)^2 + (y_2 - y_1)^2 }/f$
    //! If we wanted to insert the formula centered on a separate line:
    //! /f[
    //! \sqrt{ (x_2-x_1)^2 + (y_2 - y_1)^2 }
    //! /f]
    //! Please note that all formulas must be valid LaTeX math-mode commands. 
    //! Additionally, to be processed by Doxygen, the machine used must have 
    //! LaTeX installed. Please see the Doxygen manual for more information 
    //! about installing LaTeX locally.
    //!
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
  void LaplOperator<FE,D,TYPE>::CalcOpMat(Matrix<Double> & bMat,
                                                const LocPointMapped& lp, 
                                                BaseFE* ptFe ){
    UInt numFncs = ptFe->GetNumFncs();
    // Set correct size of matrix B and initialise with zeros
    bMat.Resize( DIM_SPACE * DIM_SPACE, numFncs * DIM_SPACE );
    bMat.Init();
    

    // Get derivatives of local shape functions with respect to global
    // coords (format: nrNodes x spaceDim)
    Matrix<Double> xiDx;
    FE *fe = (static_cast<FE*>(ptFe));
    fe->GetGlobDerivShFnc( xiDx, lp, lp.shapeMap->GetElem() , 1 );

    for(UInt k=0; k<numFncs; k++) 
    {
      for(UInt i=0, n=DIM_SPACE; i<n; i++) 
      {
        for(UInt j=0, m=DIM_SPACE; j<m; j++) 
        {
          bMat[i*DIM_SPACE+j][k*DIM_SPACE+j] = xiDx[k][i];
        }        
      }      
    }
  }

  template<class FE, UInt D, class TYPE>
  void LaplOperator<FE,D,TYPE>::CalcOpMatTransposed(Matrix<Double> & bMat,
                                                          const LocPointMapped& lp, 
                                                          BaseFE* ptFe ){
    const UInt numFncs = ptFe->GetNumFncs();
    // Set correct size of matrix B and initialise with zeros
    bMat.Resize( numFncs * DIM_SPACE, DIM_SPACE * DIM_SPACE );
    bMat.Init();
    

    // Get derivatives of local shape functions with respect to global
    // coords (format: nrNodes x spaceDim)
    Matrix<Double> xiDx;
    FE *fe = (static_cast<FE*>(ptFe));
    fe->GetGlobDerivShFnc( xiDx, lp, lp.shapeMap->GetElem() , 1 );

    for(UInt k=0; k<numFncs; k++) 
    {
      for(UInt i=0, n=DIM_SPACE; i<n; i++) 
      {
        for(UInt j=0, m=DIM_SPACE; j<m; j++) 
        {
          bMat[k*DIM_SPACE+j][i*DIM_SPACE+j] = xiDx[k][i];
        }        
      }      
    }
  }
}
#endif
