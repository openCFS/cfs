// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
#include "CoefFunctionAccumulator.hh"

namespace CoupledField  {
CoefFunctionAccumulator::CoefFunctionAccumulator(PtrCoefFct fct, 
                                                 bool integrate )
: CoefFunction(), integrate_(integrate) {
  fct_ = fct;
}

CoefFunctionAccumulator::~CoefFunctionAccumulator(){

}

void CoefFunctionAccumulator::GetTensor(Matrix<Complex>& coefMat,
                                        const LocPointMapped& lpm ) {
  
}

void CoefFunctionAccumulator::GetVector(Vector<Complex>& coefVec,
                                        const LocPointMapped& lpm ){
  
}

void CoefFunctionAccumulator::GetScalar(Complex& coef,
                                        const LocPointMapped& lpm ){
  
}

void CoefFunctionAccumulator::GetTensor(Matrix<Double>& coefMat,
                                        const LocPointMapped& lpm ){
  
}
void CoefFunctionAccumulator::GetVector(Vector<Double>& coefVec,
                                        const LocPointMapped& lpm ){
  
}

void CoefFunctionAccumulator::GetScalar(Double& coef,
                                        const LocPointMapped& lpm ){
  
}

} // end of namespace

