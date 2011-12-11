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
#include "Domain/CoordinateSystems/CoordSystem.hh"
#include "CoefFunction.hh"

namespace CoupledField{

template<class DATA_TYPE>
class CoefFunctionConst : public CoefFunction{
  public:
    CoefFunctionConst(){
      ;
    }

    virtual ~CoefFunctionConst(){
      ;
    }

    void GetTensor(Matrix<DATA_TYPE>& CoefMat, 
                   const LocPointMapped& lpm ) {
      assert(this->dimType_ == TENSOR);
      // if no coordinate system is set, just
      // use internal vector
      if( !coordSys_ ) {
        CoefMat =  constCoefMat_;
      } else {
        EXCEPTION(
            "The rotation is not fully finished ':-(\n" << 
            "Here we have to add a call to the method BaseMaterial::PerformRotation "
            "This method should be moved to the base class of the CoefFunction"
            "In addition the initial rotation of the material must be incorporated"
            "somewehre in string-notation, as we are generally dealing with string"
            "parameters."
            "Thus we should treat the case, where rotation angles are multiples of "
            "90 degree separately, where the entries are just interchanged");
      }
    }

    void GetVector(Vector<DATA_TYPE>& coefVec, 
                   const LocPointMapped& lpm ) {
      assert(this->dimType_ == VECTOR);

      // if no coordinate system is set, just
      // use internal vector
      if( !coordSys_ ) {
        coefVec = coefVec_;

      } else {
        // otherwise, perform local -> global mapping
        Vector<Double> pointCoord;
        lpm.shapeMap->Local2Global(pointCoord,lpm.lp);
        // Afterwards rotate the local vector back to global coordinates
        this->coordSys_->Local2GlobalVector( coefVec, coefVec_, pointCoord );
      }
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

    std::string ToString() const {
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
