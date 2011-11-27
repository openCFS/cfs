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

#ifndef COEFFUNCTIONCONST_HH
#define COEFFUNCTIONCONST_HH

#include "CoefFunction.hh"

namespace CoupledField{

template<class DATA_TYPE>
class CoefFunctionConst : public CoefFunction{
  public:
    CoefFunctionConst(){
      ;
    }

    ~CoefFunctionConst(){
      ;
    }

    void GetTensor(Matrix<DATA_TYPE>& CoefMat, 
                   const LocPointMapped& lpm ) {
      assert(this->dimType_ == TENSOR);
      CoefMat =  constCoefMat_;
    }

    void GetVector(Vector<DATA_TYPE>& coefVec, 
                   const LocPointMapped& lpm ) {
      assert(this->dimType_ == VECTOR);
      coefVec =  coefVec_;
    }

    void GetScalar(DATA_TYPE& coefScalar, 
                   const LocPointMapped& lpm ) {
      assert(this->dimType_ == SCALAR);
      coefScalar =  coefScalar_;
    }

    void SetTensor(Matrix<DATA_TYPE>& CoefMat){
      assert((this->dimType_ == NO_DIM) || (this->dimType_ == TENSOR) );
      this->dimType_ = TENSOR;
      this->constCoefMat_ = CoefMat;
    }

    void SetScalar(DATA_TYPE val){
      assert((this->dimType_ == NO_DIM) || (this->dimType_ == SCALAR) );
      this->coefScalar_ = val;
      this->dimType_ = SCALAR;
    }

    void SetVector(Vector<DATA_TYPE> val){
      assert((this->dimType_ == NO_DIM) || (this->dimType_ == VECTOR) );
      this->coefVec_ = val;
      this->dimType_ = VECTOR;
    }

    std::string ToString() {
      switch( dimType_ ) {
        case NO_DIM:
          return "";
          break;
        case SCALAR:
          return lexical_cast<std::string>(coefScalar_);
          break;
        case VECTOR:
          return coefVec_.ToString();
          break;
        case TENSOR:
          return constCoefMat_.ToString();
          break;
        default:
          EXCEPTION("Missing case");
      }
    }

  protected:
    Matrix<DATA_TYPE> constCoefMat_;

    Vector<DATA_TYPE> coefVec_;

    DATA_TYPE coefScalar_;

};

}
#endif
