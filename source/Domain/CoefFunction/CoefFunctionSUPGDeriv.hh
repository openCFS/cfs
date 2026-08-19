/*
 * CoefFunctionSUPGDeriv.hh
 *  Created on: 20.08.2026
 *      Author: Lorenz Klimon
 * 
 */

#ifndef SOURCE_DOMAIN_COEFFUNCTION_COEFFUNCTIONSUPGDERIV_HH
#define SOURCE_DOMAIN_COEFFUNCTION_COEFFUNCTIONSUPGDERIV_HH

#include "CoefFunction.hh"
#include "FeBasis/BaseFE.hh"
#include "FeBasis/FeFunctions.hh"
#include "FeBasis/H1/H1Elems.hh"

namespace CoupledField {

class CoefFunctionSUPGDeriv : public CoefFunction {
public:
  //! Constructor
  CoefFunctionSUPGDeriv(PtrCoefFct velocityField, PtrCoefFct materialCoeff, shared_ptr<BaseFeFunction> feFnc);

  //! Destructor
  virtual ~CoefFunctionSUPGDeriv() { }

  virtual string GetName() const override { return "CoefFunctionSUPGDeriv"; }

  //! \see CoefFunction::GetScalar
  virtual void GetScalar(Double& scal, const LocPointMapped& lpm) override;

protected:
  PtrCoefFct velField_;
  PtrCoefFct matCoeff_;
  shared_ptr<BaseFeFunction> feFct_;
};

} // namespace CoupledField

#endif /* SOURCE_DOMAIN_COEFFUNCTION_COEFFUNCTIONSUPGDERIV_HH */