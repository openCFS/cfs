// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef COEF_FUNCTION_COMPLEXTOREAL_HH
#define COEF_FUNCTION_COMPLEXTOREAL_HH

#define _USE_MATH_DEFINES
#include <cmath>

#include <string>
#include <map>

#include "CoefFunction.hh"
#include "FeBasis/BaseFE.hh"
#include "Materials/BaseMaterial.hh"
#include "FeBasis/FeFunctions.hh"

namespace CoupledField  {

class MathParser;

template <class T>
class CoefFunctionComplexToReal : public CoefFunction {

public:
  //! Constructor
  CoefFunctionComplexToReal(PtrCoefFct coefToChange,
                         Grid* ptGrid);

  //! Destructor
  virtual ~CoefFunctionComplexToReal();


  //! \copydoc CoefFunction::GetScalar
  virtual void GetScalar(T& coefScal, const LocPointMapped& lpm );


  //! \copydoc CoefFunction::GetVector
  virtual void GetVector(Vector<T>& coefVec, const LocPointMapped& lpm);

  //! \copydoc CoefFunction::GetTensor
  virtual void GetTensor(Matrix<T>& coefMat, const LocPointMapped& lpm);

  //! \copydoc CoefFunction::ToString
  virtual std::string ToString() const;

  virtual string GetName() const { return "CoefFunctionComplexToReal"; }


  //! \copydoc CoefFunction::GetVecSize
  virtual UInt GetVecSize() const;

private:

  //! Pointer to grid object
  Grid * ptGrid_;

  //! CoefFunction to change it's type
  PtrCoefFct coefToChange_ = NULL;

};
} //end of namespace

#endif
