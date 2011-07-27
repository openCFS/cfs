#ifndef FILTER_HH_
#define FILTER_HH_

#include "General/Enum.hh"

namespace CoupledField
{

class DesignSpace;

/** This is an information container for the Filter which is stored in the DesignElemen!.
 * It is used by DesignStructure to initialize the element filters.
 * Set in DesignStructure, stored in DesignElement */
struct Filter
{
public:
  Filter();

  /** type of of filter */
  typedef enum { NO_FILTERING, SENSITIVITY, DENSITY } Type;

  /** Subtype for senisitivy filtering. w = weight, p is density, f' is cost gradient
   * Sigmund  = sum_i w(x_i) * p_i * f_i' / p_e * sum_i w(x_i)
   * sharp Sigmund  = sum_i (i=e ? 1:0 : w(x_i)) * p_i * f_i' / p_e * sum_i w(x_i) + "bug" in normalized weighting
   * Borrvall = sum_i w(x_i) * p_i * f_i' / sum_i p_i * w(x_i)
   * plain    = sum_i w(x_i) * f_i' / sum_i w(x_i)
   * sharp plain = plain but the "bug" in normalized weighting as in sharp Sigmund */

  typedef enum { PLAIN, SHARP_PLAIN, SIGMUND, SHARP_SIGMUND, BORRVALL } Sensitivity;

  /** Subtype of design filter
   * See Sigmund; Morphology based black and white filters for topology optimization; 2007
   * Standard: the plain filter
   * modified heaviside, Sigmund (29), postproc stanard to 1 or max
   * tanh: Variant of the Xu-Filter as in Wang,Lazarow,Sigmund; On projection methods, ...;2010 but simpler implementation!*/
  typedef enum { STANDARD, HEAVISIDE, MOD_HEAVISIDE, TANH } Density;

  /** Handled in DesignElement.cc */
  static Enum<Type>        type;
  static Enum<Sensitivity> sensitivity;
  static Enum<Density>     density;

  bool Enabled() const { return type_ != NO_FILTERING; }

  double GetBeta() const { return beta_; }

  /** Setting beta triggers the calculation of heaviside_corr.
   * The other settings need to be set! */
  void SetBeta(double beta, const DesignSpace* space);

  Type        type_;
  Sensitivity sensitivity_;
  Density     density_;

  /** this is the correction parameter for the (not modified) heaviside filter,
   * such that H(rho_min) = rho_min -> otherwise H(0.001) -> 1 for beta -> inf
   * the correction value found by bisection on any new beta value */
  double    heaviside_corr;

  /** switching parameter for tanh */
  double eta;

private:
  /** this is the beta parameter for the (modified) heaviside or tanh design filter.
   * Private such that we force setting heaviside_corr on setting beta */
   double     beta_;

};

}

#endif /* FILTER_HH_ */
