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

#include "coefFunction.hh"

#include "General/environment.hh"
#include "MatVec/matrix.hh"
#include "Domain/elem.hh"
#include "Elements/elemShapeMap.hh"

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

    void GetTensor(Matrix<DATA_TYPE>& CoefMat, LocPointMapped lp, const Elem* elem){
      assert(this->mType_ == TENSOR);
      CoefMat =  constCoefMat_;
    }

    void GetVector(Vector<DATA_TYPE>& coefVec, LocPointMapped lp, const Elem* elem){
      assert(this->mType_ == VECTOR);
      coefVec =  coefVec_;
    }

    void GetScalar(Matrix<DATA_TYPE>& coefScalar, LocPointMapped lp, const Elem* elem){
      assert(this->mType_ == SCALAR);
      coefScalar =  coefScalar_;
    }

    void SetTensor(Matrix<DATA_TYPE>& CoefMat){
      assert((this->mType_ == UNDEF) || (this->mType_ == TENSOR) );
      this->mType_ = TENSOR;
      this->constCoefMat_ = CoefMat;
    }

    void SetScalar(DATA_TYPE val){
      assert((this->mType_ == UNDEF) || (this->mType_ == SCALAR) );
      this->coefScalar_ = val;
      this->mType_ = SCALAR;
    }

    void SetVector(Vector<DATA_TYPE> val){
      assert((this->mType_ == UNDEF) || (this->mType_ == VECTOR) );
      this->coefVec_ = val;
      this->mType_ = VECTOR;
    }
  protected:
    Matrix<DATA_TYPE> constCoefMat_;

    Vector<DATA_TYPE> coefVec_;

    DATA_TYPE coefScalar_;

};

}
#endif
