// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef COEF_FUNCTION_ACCUMULATOR_HH
#define COEF_FUNCTION_ACCUMULATOR_HH

#include "CoefFunction.hh"

namespace CoupledField  {

//! Meta-CoefFunction (proxy) which collects norm information

//! This class is used to encapsulate a normal scalar/tensor-valued 
//! CoefFunction and compute the norm of accessed values. 
//! Most calls are just forwarded to the internal CoefFunction instance.
class CoefFunctionAccumulator : public CoefFunction {
public:
  //! Constructor
  CoefFunctionAccumulator(PtrCoefFct fct, bool integrate );
  
  //! Destructor 
  virtual ~CoefFunctionAccumulator();
  
  std::string GetName() const { return "CoefFunctionAccumulator"; }

  // ========================
   //  ACCUMULATOR METHODS
   // ========================
   //@{ \name Statistic methods
  void ResetSampling() {
    squaredSum_ = 0.0;
//    sum_.Init();
  }
   
  //! Return norm of accessed values since last reset
  Double GetNorm() {
    return std::sqrt(squaredSum_);
  }
  
//  //! Get vector norm
//  Double GetVecNorm() {
//    return sum_.NormL2();
//  }
  //@}
  
  virtual CoefDimType GetDimType() const{
	  return fct_->GetDimType();
  }

  void SetIntegrateFlag( bool integrate ) {
    integrate_ = integrate;
  }

  
  // ========================
  //  ACCESS METHODS
  // ========================
  //@{ \name Access Methods
  //! \copydoc CoefFunction::GetTensor
  virtual void GetTensor(Matrix<Complex>& coefMat,
                         const LocPointMapped& lpm );
  
  //! \copydoc CoefFunction::GetVector
  virtual void GetVector(Vector<Complex>& coefVec,
                         const LocPointMapped& lpm );

  //! \copydoc CoefFunction::GetScalar
  virtual void GetScalar(Complex& coef,
                         const LocPointMapped& lpm );
  
  //! \copydoc CoefFunction::GetTensor
  virtual void GetTensor(Matrix<Double>& coefMat,
                         const LocPointMapped& lpm );
  
  //! \copydoc CoefFunction::GetVector
  virtual void GetVector(Vector<Double>& coefVec,
                         const LocPointMapped& lpm );

  //! \copydoc CoefFunction::GetScalar
  virtual void GetScalar(Double& coef,
                         const LocPointMapped& lpm );
  
  //! \copydoc CoefFunction::GetVecSize
  virtual UInt GetVecSize() const {
    return fct_->GetVecSize();
  }

  //! \copydoc CoefFunction::GetTensorSize
  virtual void GetTensorSize( UInt& numRows, UInt& numCols ) const{
    fct_->GetTensorSize(numRows, numCols);
  }
  //@}

  //! \copydoc CoefFunction::ToString
  virtual std::string ToString() const {
    return fct_->ToString();
  }
  
private:
  
  //! Encapsulated CoefFunction
  PtrCoefFct fct_;
  
  //! Flag, if coefFunction should be integrated
  bool integrate_;
  
  //! Squared total sum of accessed entries
  Double squaredSum_;
  
//  //! Total vector sum
//  Vector<Double> sum_;
//  
  
};

} // end of namespace
#endif
