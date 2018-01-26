#ifndef FILTER_HH_
#define FILTER_HH_

#include "General/Enum.hh"

namespace CoupledField
{

class DesignSpace;
class DesignElement;

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
  typedef enum { STANDARD, SOLID_HEAVISIDE, VOID_HEAVISIDE, TANH } Density;

  /** Handled in DesignElement.cc */
  static Enum<Type>        type;
  static Enum<Sensitivity> sensitivity;
  static Enum<Density>     density;

  bool Enabled() const { return type_ != NO_FILTERING; }

  /** Convenience function. Gives a lower bound. The explicit filter bound if given, otherwise
   * from the design element */
  Double GetLowerBound(const DesignElement* de) const;

  /** Set an explicit lower bound to overwrite the design's lower bound */
  void SetLowerBound(Double value) { explicit_lower_bound_ = value; }

  /** Set non_lin_scale and non_lin_offset. considers explicit_lower_bound_
   * @peram de reference design element for bounds */
  void SetNonLinCorrection(const DesignElement* de, unsigned int filter_idx);

  void SetType(Type t) { type_ = t; }

  Type GetType() const { return type_; }

  /** Sums up the weights of the neighbors and optionally the own element */
  Double CalcWeightSum(bool include_this) const;

  void Dump() const;

  /** Neighborhood is element and pre-calculated distance. It is stored in the Filter! */
  struct NeighbourElement
  {
  public:
    /** read the variable */
    DesignElement* neighbour;

    /** pre-calculated weight: radius - distanance and >= 0 */
    Double        weight;

    /** the distance in domain dimensions! */
    Double        distance;
  };


  /** pre-calculated weight sum */
  Double weight_sum;

  /** The weight of THIS element which is radius */
  Double weight;

  /** The neighbors if filter otherwise empty.
   * The element itself is NOT part of the neighborhood!
   * @see DesignStructure::DesignStructure() */
  StdVector<NeighbourElement> neighborhood;

  Sensitivity sensitivity_;
  Density     density_;

  /** to have F(rho_max) = rho_max and F(rho_min) = rho_min we need a scaling and a offset.
   * This applies for Heaviside and tanh.
   * With solid_heaviside for a large beta F(x > 0) -> 1
   * With tanh for a small beta F(0) >> 0 and F(1) << 1. With eta != 0.5 this is unsymmetric */
  Double non_lin_scale;
  Double non_lin_offset;

  /** this is the beta parameter for the heaviside filters or tanh design filter. */
  Double     beta;

  /** switching parameter for tanh */
  Double eta;

  /** to check where we belong to */
  RegionIdType region;

private:

  Type type_;

  /** Holds the optional "force_lower_bound" attribute to overwrite the design lower bound
   * for Heaviside and tanh type filters. Necessary for mixed design scenarios where one has
   * a fixed design with equal lower and upper bound. */
   Double explicit_lower_bound_;

};

}

#endif /* FILTER_HH_ */
