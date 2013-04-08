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

namespace CoupledField {

// ============================================================================
//  Stabilization parameters for SUPG
// ============================================================================
//! Computes the stabilization parameters for SUPG


class CoefFunctionStabParams : public CoefFunction{
public:

  //! Constructor
  CoefFunctionStabParams(PtrCoefFct density, PtrCoefFct viscosity);
  //! Destructor
  virtual ~CoefFunctionStabParams(){;}

  //! Return cpmplex-valued scalar at integration point
  void GetScalar(Double& val,
                 const LocPointMapped& lpm ) ;

protected:
  //! density
  PtrCoefFct density_;

  //! viscosity
  PtrCoefFct viscosity_;
};
}

#endif
