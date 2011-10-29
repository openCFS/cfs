// =====================================================================================
// 
//       Filename:  coefFunctionAnalytic.hh
// 
//    Description:  Provides analytical parameters to operators and bilinear forms
//                  Right now this is just a wrapping for the BaseMaterial
//
// 
//        Version:  1.0
//        Created:  10/19/2011 11:14:50 AM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Andreas Hueppe (AHU), andreas.hueppe@uni-klu.ac.at
//        Company:  Universitaet Klagenfurt
// 
// =====================================================================================

#ifndef COEFFUNCTIONANALYTIC_HH
#define COEFFUNCTIONANALYTIC_HH

#include "coefFunction.hh"

#include "General/environment.hh"
#include "MatVec/matrix.hh"
#include "Domain/elem.hh"
#include "Elements/elemShapeMap.hh"

namespace CoupledField{


class CoefFunctionAnalytic : public CoefFunction{
  public:
    CoefFunctionAnalytic(){
      ;
    }

    ~CoefFunctionAnalytic(){
      ;
    }

    void GetMatrix(Matrix<Double>& CoefMat, LocPointMapped lp, const Elem* elem){
      CoefMat =  constCoefMat_;
    }

    void GetMatrix(Matrix<Complex>& CoefMat, LocPointMapped lp, const Elem* elem){
      CoefMat =  constCoefCMat_;
    }

    void SetConstMatrix(Matrix<Double>& CoefMat){
      this->constCoefMat_ = CoefMat;
    }

    void SetConstMatrix(Matrix<Complex>& CoefMat){
      this->constCoefCMat_ = CoefMat;
    }

  protected:
    Matrix<Double> constCoefMat_;

    Matrix<Complex> constCoefCMat_;

};

}
#endif
