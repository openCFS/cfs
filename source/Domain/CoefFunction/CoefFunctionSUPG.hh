/*
 * CoefFunctionSUPG.hh
 *
 *  Created on: 25.04.2024
 *      Author: Vladimir Filkin and Lukas Schafferhofer
 */

#ifndef SOURCE_DOMAIN_COEFFUNCTION_COEFFUNCTIONUPCOTH
#define SOURCE_DOMAIN_COEFFUNCTION_COEFFUNCTIONUPCOTH

#include "CoefFunction.hh"
#include "FeBasis/BaseFE.hh"
#include "FeBasis/FeFunctions.hh"
#include "FeBasis/H1/H1Elems.hh"
#include "FeBasis/HCurl/HCurlElems.hh"

namespace CoupledField {
// forward class declaration
// ============================================================================
//  Coef Function Up Coth
// ============================================================================
//! Coefficient function for calculating the stabilization factor for SUPG
//! \note This class only works for real-valued scalar data.

  class CoefFunctionSUPG : public CoefFunction {

  public:

    //! Constructor
    CoefFunctionSUPG(PtrCoefFct velocityField,  PtrCoefFct materialCoeff, shared_ptr<BaseFeFunction> FeFnc);

    //! Destructor
    virtual ~CoefFunctionSUPG(){;}

    virtual string GetName() const { return "CoefFunctionSUPG"; }

    //! \see CoefFunction::GetScalar
    virtual void GetScalar(Double& scal, const LocPointMapped& lpm);

    //! \see CoefFunction::GetVector
    virtual void GetVector(Vector<Double>& vec, const LocPointMapped& lpm );

  protected:

    //! Set feFunction and xml input data
    //shared_ptr<BaseFeFunction> FeFnc;
    shared_ptr<BaseFeFunction>  feFct_;
    
    //! Set feFunction and xml input data
    //shared_ptr<BaseFeFunction> FeFnc;
    PtrCoefFct matCoeff_;

    //! velocity (or weighted with material parameters velocity)
    PtrCoefFct velField_;

    //Coeffunction, which is used to evalute the model
    PtrCoefFct depCoef_ ;

  };

}


#endif /* SOURCE_DOMAIN_COEFFUNCTION_COEFFUNCTIONUPCOTH */
