// =====================================================================================
//
//       Filename:  gradientOp.hh
//
//    Description:  This class implements the gradient operator
//                  It returns and applies matrices of the form
//                     / N_1x N_2x ...\
//                  b =| N_1y N_2y ...|
//                     \ N_1z N_2z .../
//                  Where N_1x denotes the x-derivative of the first
//                  shape function at a given local point
//
//        Version:  1.0
//        Created:  10/03/2011 11:00:04 AM
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Andreas Hueppe (AHU), andreas.hueppe@uni-klu.ac.at
//        Company:  Universitaet Klagenfurt
//
// =====================================================================================

#ifndef GRADIENTOP_HH
#define GRADIENTOP_HH

#include "BaseBOperator.hh"

namespace CoupledField{
  template<class FE, class TYPE>
  class GradientOperator : public BaseBOperator<FE,TYPE>{
    public:
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

  template<class FE, class TYPE>
  void GradientOperator<FE,TYPE>::CalcOpMat(Matrix<Double> & bMat,
                                            const LocPointMapped& lp, 
                                            BaseFE* ptFe ){
    UInt numFncs = ptFe->GetNumFncs();
    // Set correct size of matrix B and initialise with zeros
    UInt eDim = ptFe->shape_.dim;
    bMat.Resize( eDim, numFncs );
    bMat.Init();

    // Get derivatives of local shape functions with respect to global
    // coords (format: nrNodes x spaceDim)
    Matrix<Double> xiDx;
    FE *fe = (static_cast<FE*>(ptFe));
    fe->GetGlobDerivShFnc( xiDx, lp, lp.shapeMap->GetElem() , 1 );
    xiDx.Transpose(bMat);

  }

  template<class FE, class TYPE>
  void GradientOperator<FE,TYPE>::CalcOpMatTransposed(Matrix<Double> & bMat,
                                                      const LocPointMapped& lp, 
                                                      BaseFE* ptFe ){
    UInt numFncs = ptFe->GetNumFncs();
    // Set correct size of matrix B and initialise with zeros
    UInt eDim = ptFe->shape_.dim;
    bMat.Resize(numFncs , eDim );
    bMat.Init();

    // Get derivatives of local shape functions with respect to global
    // coords (format: spaceDim x nrNodes )
    FE *fe = (static_cast<FE*>(ptFe));
    fe->GetGlobDerivShFnc( bMat, lp, lp.shapeMap->GetElem() , 1 );


  }
}
#endif
