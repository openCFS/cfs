
#ifndef FILE_PHYSICAL_CONSTANTS
#define FILE_PHYSICAL_CONSTANTS

#include <boost/math/constants/constants.hpp>

#include "General/defs.hh"

namespace CoupledField {

  //! Physical constants which are shared by several PDEs / material models
  //! Define them here once instead of hard-coding them in every class
  namespace PhysicalConstants {

    //! Vacuum permeability mu_0 in Vs/(Am)
    inline const Double MU_ZERO = 4.0 * boost::math::constants::pi<Double>() * 1e-7;

    //! Vacuum reluctivity nu_0 = 1/mu_0 in Am/(Vs)
    inline const Double NU_ZERO = 1.0 / MU_ZERO;

    //! Vacuum permittivity eps_0 in As/(Vm)
    inline const Double EPS_ZERO = 8.8541878128e-12;

  } // namespace PhysicalConstants

} // namespace CoupledField

#endif
