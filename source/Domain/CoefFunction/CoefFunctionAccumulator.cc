// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
#include "CoefFunctionAccumulator.hh"

namespace CoupledField  {
CoefFunctionAccumulator::CoefFunctionAccumulator(PtrCoefFct fct, 
                                                 bool integrate )
: CoefFunction(), integrate_(integrate) {
  squaredSum_ = 0.0;
  fct_ = fct;
  dependType_ = fct->GetDependency(); 
  
//  sum_.Resize(fct_->GetVecSize());
  
}

CoefFunctionAccumulator::~CoefFunctionAccumulator(){

}

void CoefFunctionAccumulator::GetTensor(Matrix<Complex>& coefMat,
                                        const LocPointMapped& lpm ) {
  REFACTOR
}

void CoefFunctionAccumulator::GetVector(Vector<Complex>& coefVec,
                                        const LocPointMapped& lpm ){
  #pragma omp critical (COEFFUNCTIONACCUMULATOR_SUMMATION)
  {
    fct_->GetVector(coefVec, lpm);
    Complex temp;
    for( UInt i = 0; i < coefVec.GetSize(); ++i ) {
      if( integrate_ ) {
        temp += coefVec[i] * coefVec[i] * lpm.weight * lpm.jacDet;
      } else {
        temp += coefVec[i] * coefVec[i];
      }
    }
    squaredSum_ += std::abs(temp);
  }
}

void CoefFunctionAccumulator::GetScalar(Complex& coef, const LocPointMapped& lpm ){
	fct_->GetScalar(coef, lpm);

#pragma omp critical (CoefFunctionAccumulator)
	{
		//Double square = coef.real()*coef.real() - coef.imag()*coef.imag();
		if( integrate_ )
			//squaredSum_ = square * lpm.weight * lpm.jacDet;
			squaredSum_ += std::abs(coef*coef) * lpm.weight * lpm.jacDet;
		else
			squaredSum_ += std::abs(coef*coef);
	}

}

void CoefFunctionAccumulator::GetTensor(Matrix<Double>& coefMat,
                                        const LocPointMapped& lpm ){
  REFACTOR
  
}
void CoefFunctionAccumulator::GetVector(Vector<Double>& coefVec,
                                        const LocPointMapped& lpm ){

#pragma omp critical (COEFFUNCTIONACCUMULATOR_SUMMATION)
{
  fct_->GetVector(coefVec, lpm);
  for( UInt i = 0; i < coefVec.GetSize(); ++i ) {
    if( integrate_ ) {
      squaredSum_ += coefVec[i] * coefVec[i] * lpm.weight * lpm.jacDet;
    } else {
      squaredSum_ += coefVec[i] * coefVec[i];
    }
  }
}
  // code for vector norm
  //  if( integrate_ ) {
  //    sum_ += (coefVec * lpm.weight * lpm.jacDet);
  //  } else {
  //    sum_ += coefVec;
  //  }
}

void CoefFunctionAccumulator::GetScalar(Double& coef,
                                        const LocPointMapped& lpm ){

	fct_->GetScalar(coef, lpm);

#pragma omp critical (CoefFunctionAccumulator)
{
	if( integrate_ )
		squaredSum_ += coef * coef * lpm.weight * lpm.jacDet;
	else
		squaredSum_ += coef * coef;
}

}

} // end of namespace

