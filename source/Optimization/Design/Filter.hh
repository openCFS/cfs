#ifndef FILTER_HH_
#define FILTER_HH_

#include "General/Enum.hh"
#include "Optimization/Tune.hh"
#include "Optimization/TransferFunction.hh"
#include "Optimization/Design/BaseDesignElement.hh"

namespace CoupledField
{

class DesignSpace;
class DesignElement;
class MathParser;
struct GlobalFilter;

/** This is an information container for the Filter which is stored in the DesignElemen!.
 * The Filter information is split between local element data and a connected GlobalFilter.
 * Some "global" functions are just forwarded to global filter.
 * It is used by DesignStructure to initialize the element filters.
 * Set in DesignStructure, stored in DesignElement */
struct Filter
{
public:
  /** type of of filter */
  typedef enum { NO_FILTERING, SENSITIVITY, DENSITY } Type;

  /** Subtype for senisitivy filtering. w = weight, p is density, f' is cost gradient
   * Sigmund  = sum_i w(x_i) * p_i * f_i' / p_e * sum_i w(x_i)
   * Sigmund_Safe: use Sigmund, but check that denominator is not close to zero otherwise don't filter
   * sharp Sigmund  = sum_i (i=e ? 1:0 : w(x_i)) * p_i * f_i' / p_e * sum_i w(x_i) + "bug" in normalized weighting
   * Borrvall = sum_i w(x_i) * p_i * f_i' / sum_i p_i * w(x_i)
   * Borrvall_Safe: use Borrvall, but check that denominator is not close to zero otherwise don't filter
   * plain    = sum_i w(x_i) * f_i' / sum_i w(x_i)
   * sharp plain = plain but the "bug" in normalized weighting as in sharp Sigmund */

  typedef enum { PLAIN, SHARP_PLAIN, SIGMUND,SIGMUND_TRACE, SHARP_SIGMUND, BORRVALL } Sensitivity;

  /** Subtype of design filter
   * See Sigmund; Morphology based black and white filters for topology optimization; 2007
   * Standard: the plain filter
   * modified heaviside, Sigmund (29), postproc stanard to 1 or max
   * tanh: Variant of the Xu-Filter as in Wang,Lazarow,Sigmund; On projection methods, ...;2010 but simpler implementation!
   * expression: math parser exression with child element 'exression'
   * material: from Lukas: rho_material(x)=mp(rho(x)) * ms(density_filter(mf(x))) with material filter transfer functions phase, scale, and filter
   * material_part: for computational use only: only density_filter(mf(x)) without transfer functions phase and scale */
  typedef enum { STANDARD, SOLID_HEAVISIDE, VOID_HEAVISIDE, TANH, EXPRESSION, MATERIAL, MATERIAL_PART } Density;

  /** the way of the weighting in the filter. CONSTANT e.g. for MAX filter */
  typedef enum { NO_CONTRIBUTION, LINEAR, CONSTANT } Contribution;

  /** Filter types we have
   * <ul>
   *   <li>RADIUS: this is the implementation following Sigmund in the 99lines paper.
   *               The drawback is the discretization dependency.</li>
   *   <li>VOLUME_RADIUS: The radius is *value* times square/cube edge length where the
   *               square/cube has the volume of the element</li>
   *   <li>MAX_EDGE: The largest edge size, discretization independent and preferable</li>
   * </ul>
   * This does not tell if we have design or sensitivity filtering! */
  typedef enum { NO_FILTER = -1, RADIUS, VOLUME_RADIUS, MAX_EDGE } FilterSpace;

  /** Also to be used in Function::Local */
  static Enum<FilterSpace> filterSpace;

  /** Handled in DesignElement.cc */
  static Enum<Type>        type;
  static Enum<Sensitivity> sensitivity;
  static Enum<Density>     density;

  bool Enabled() const;
  
  void SetType(Type t);

  Type GetType() const;

  double GetBeta() const;

  double GetEta() const;

  /** Sums up the weights of the neighbors and optionally the own element */
  double CalcWeightSum(bool include_this) const
  {
	  double res = 0.0;

	  for(unsigned int i = 0, n = neighborhood.GetSize(); i < n; i++)
		  res += neighborhood[i].weight;

	  if(include_this)
		  res += this->weight;

	  return res;
  }

  void Dump() const;

  /** Neighborhood is element and pre-calculated distance. It is stored in the Filter! */
  struct NeighbourElement
  {
  public:
    /** read the variable */
    DesignElement* neighbour = NULL;

    /** pre-calculated weight: radius - distanance and >= 0 */
    double        weight;

    /** the distance in domain dimensions! */
    double        distance;
  };

  /** pre-calculated weight sum including this-weight. Set in DesignStructure and always available */
  double weight_sum = -1.0;

  /** The weight of THIS element which is radius */
  double weight = 1.0;

  /** The neighbors if filter otherwise empty.
   * The element itself is NOT part of the neighborhood!
   * @see DesignStructure::DesignStructure() */
  StdVector<NeighbourElement> neighborhood;

  /** here we story common data for a region/design/excitation unit.
   * Object is in DesignSpace */
  GlobalFilter* global = NULL;
};


/** a global filter unit is per region/design/excitation and stored in DesignSpace */
struct GlobalFilter
{
public:
  GlobalFilter();
  ~GlobalFilter();

  /** the copy constructor cannot handle the math parser and Tune properly.
   * ReInitParser() creates a new math parser handle and does not release the old.
   * We need to call this every time we create GF instance! */
  void PostCopy(bool register_tune);

  /** math parser expression */
  void InitParser(const string& func_expr, const string& deriv_expr);

  /** the default copy constructor kills the math parser, it is registered with
   * the double reference of the original class.
   * This overwrites the parser but does not release the handles as the id's are from
   * the original class */
  void ReInitParser();

  bool Enabled() const { return type != Filter::NO_FILTERING; }

  /** Convenience function. Gives a lower bound. The explicit filter bound if given, otherwise
   * from the design element */
  double GetLowerBound(const DesignElement* de) const;

  /** Set non_lin_scale and non_lin_offset. Harmless for not nonlinear filters
  * @param de reference design element for bounds */
  void SetNonLinCorrection(const DesignElement* ref);
  
  void SetType(Filter::Type t) {type = t;}

  Filter::Type GetType() const { return type; }

  /** evaluate the math parser expression with optional applying nonlinear scaling and offset */
  double EvalExpressionFunction(double x, bool scale) const;
  double EvalExpressionDerivative(double x, bool scale) const;

  /** some debug information */
  std::string ToString() const;

  /** write data of the filter to the given info node */
  void ToInfo(PtrParamNode info);

  Filter::Type type = Filter::NO_FILTERING;

  /** this is the beta parameter for the heaviside filters or tanh design filter.
   * the mathparser uses it as reference variable to allow it updated */
  double beta = -1.0;

  /** switching parameter for tanh */
  double eta = -1.0;

  /** the filter value */
  double value = -1.0;

  /** to check where we belong to */
  RegionIdType region = NO_REGION_ID;

  BaseDesignElement::Type design = BaseDesignElement::NO_TYPE;

  /** robust excitation index */
  int robust = -1;

  /** Parameter for filter */
  Filter::Contribution contribution = Filter::NO_CONTRIBUTION;

  Filter::FilterSpace filterspace = Filter::NO_FILTER;

  Filter::Sensitivity sensitivity = Filter::PLAIN;
  Filter::Density     density = Filter::STANDARD;

  /** number of elements */
  int elements = 0;

  /** average radius */
  double avg_radius = 0.0;

  /** average neighbor size */
  double avg_neigbor = 0.0;

  /** to have F(rho_max) = rho_max and F(rho_min) = rho_min we need a scaling and a offset.
   * This applies for Heaviside and tanh.
   * With solid_heaviside for a large beta F(x > 0) -> 1
   * With tanh for a small beta F(0) >> 0 and F(1) << 1. With eta != 0.5 this is unsymmetric */
  double non_lin_scale = 1.0;
  double non_lin_offset = 0.0;

  /** Material transfer function to scale the rho within the density filter */
  TransferFunction mat_filter;

  /** transfer function to scale the filtered density for physical correction */
  TransferFunction mat_scale;

  /** transfer function as function of plain rho be multiplied with scaled filter */
  TransferFunction mat_phase;

  /** optional automatic beta setter. Has Init() and a late post init Register() because of the GlobalFilter copy constructor issue */
  Tune tune;

  /** in the robust case, only one filter has a own beta Tune, the other filters register directly there */
  Tune* ext_tune = nullptr;

  /** when we want to set ext_tune but another tune is not registered yet */
  static const int UNSET_EXT_BETA_TUNE = 0x12345;

  /** for expression */
  MathParser* parser_ = nullptr;

  /** math parser handle */
  unsigned int function_handle_ = 0;
  unsigned int derivative_handle_ = 0;
private:
  /** the math parser input (density filtered value) - requires ReInitiParser() */
  mutable double expression_x_ = -1.2;

};

/** needs to be here due to forward declaration :( */
inline bool Filter::Enabled() const
{
  return global->GetType() != NO_FILTERING;
}

inline void Filter::SetType(Type t)
{
  global->SetType(t);
}

inline Filter::Type Filter::GetType() const
{
  return global->GetType();
}

inline double Filter::GetBeta() const
{
  return global->beta;
}

inline double Filter::GetEta() const
{
  return global->eta;
}

} // end of namespace

#endif /* FILTER_HH_ */
