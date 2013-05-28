//================================================================================================
/*!
 *       \file     CoefFunctionStabParams.hh
 *       \brief    Coefficient function which compute sthe stabilization parameters for SUPG
 *
 *       \date     07/04/2013
 *       \author   Manfred Kaltenbacher
 */
//================================================================================================

#ifndef FILE_COEFFUNCTION_STABPARAMS_HH
#define FILE_COEFFUNCTION_STABPARAMS_HH

#include <boost/tr1/type_traits.hpp>
#include "CoefFunction.hh"
#include "CoefFunctionMulti.hh"
#include "FeBasis/FeFunctions.hh"
#include "Forms/Operators/BaseBOperator.hh"

namespace CoupledField {

// forward class declaration
class BaseBOperator;

// ============================================================================
//  Stabilization parameters for SUPG
// ============================================================================
//! Computes the stabilization parameters for SUPG


class CoefFunctionStabParams : public CoefFunction{
public:

  //! Define stabilization type
	typedef
			enum { NO_TYPE, SUPG, PSPG, LSIC } StabType;

  //! Constructor
  CoefFunctionStabParams(PtrCoefFct density, PtrCoefFct viscosity,
		  	  	  	  	  std::string type);


  //! Constructor for additional mean flow coefFunction
  CoefFunctionStabParams(PtrCoefFct density, PtrCoefFct viscosity,
		  	  	  	  	  PtrCoefFct meanFlow,
		  	  	  	  	  BaseBOperator* opt,
		  	  	  	  	  shared_ptr<BaseFeFunction > feFnc,
		  	  	  	  	  std::string type);

  //! Destructor
  virtual ~CoefFunctionStabParams(){;}

  //! Return complex-valued scalar at integration point
  void GetScalar(Double& val,
                 const LocPointMapped& lpm ) ;

protected:

  //! do initilization
  void PerformInitilization(PtrCoefFct density, PtrCoefFct viscosity,
		  	  	  	  	  	  std::string type);

  //! density
  PtrCoefFct density_;

  //! viscosity
  PtrCoefFct viscosity_;

  //! true, if mean flow is set
  bool IsSetMeanFlow_;

  //! stabilization type
  StabType stabType_;

  //! Coefficient function for the flow field

  //! This coefficient function describes the flow field. As this
  //! is in general different for each region and will most likely
  //! not be given in a close form, it is described by a CoefFunctionMulti.
  PtrCoefFct meanFlowCoef_;

  //! Differential operator to calculate the field value
  BaseBOperator* bOperator_;


  //! Depending FeFunction
  shared_ptr<BaseFeFunction>  feFct_;

};
}

#endif
