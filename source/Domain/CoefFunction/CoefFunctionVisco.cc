/*
 * CoefFunctionVisco.cc
 *
 *  Created on: 22 Feb 2016
 *      Author: INTERN\ftoth
 */

#include "CoefFunctionVisco.hh"

namespace CoupledField {
  CoefFunctionVisco::CoefFunctionVisco(StdVector<PtrCoefFct> relaxationTensors, StdVector<PtrCoefFct> pronyFunctions) {
    dimType_ = VECTOR; // should be determined from input
    dependType_ = SOLUTION; // should be determined from input
    isAnalytic_ = false;
    isComplex_  = false;
    // check length
    if (relaxationTensors.GetSize() == pronyFunctions.GetSize()){
        N_ = relaxationTensors.GetSize();
    }
    else {
        EXCEPTION("Vectors must be of the same length");
    }
    // get result size by checking first one
    UInt aM,aN,xN;
    relaxationTensors[0]->GetTensorSize(aM,aN);
    xN = pronyFunctions[0]->GetVecSize();
    if (aN == xN) {
        vecSize_ = aN;
    }
    else {
        EXCEPTION( "Size must be compatible" )
    }
    // check all
    for (UInt i=0; i<N_;i++) {
        // type
        if (relaxationTensors[i]->GetDimType() != TENSOR) {
            EXCEPTION("All 'relaxationTensors' must be of type 'TENSOR', but "<<(i+1)
                      <<"th element is "<<relaxationTensors[i]->GetDimType());
        }
        if (pronyFunctions[i]->GetDimType() != VECTOR) {
            EXCEPTION("All 'pronyFunctions' must be of type 'VECTOR', but "<<(i+1)
                      <<"th element is "<<relaxationTensors[i]->GetDimType());
        }
        // size
        relaxationTensors[i]->GetTensorSize(aM,aN);
        xN = pronyFunctions[i]->GetVecSize();
        if (aN != xN) {
            EXCEPTION( "Size must be compatible" )
        }
        if (aN != vecSize_) {
            EXCEPTION( "Return size must be constant for all prony terms" )
        }
    }
    // check size
    Cs_ = relaxationTensors;
    ps_ = pronyFunctions;
  }

  CoefFunctionVisco::~CoefFunctionVisco() {
    ;
  }

  void CoefFunctionVisco::GetVector(Vector<Double>& vec, const LocPointMapped& lpm ) {
    vec.Resize(vecSize_,0.0);
    Matrix<Double> C;
    Vector<Double> p;
    // loop over all and sum
    for (UInt i=0; i<N_;i++) {
        Cs_[i]->GetTensor(C,lpm);
        ps_[i]->GetVector(p,lpm);
        vec = vec + C*p;
    }
  }

  UInt CoefFunctionVisco::GetVecSize() const{
    return vecSize_;
  }

}
