#ifndef DESIGNELEMENT_HH_
#define DESIGNELEMENT_HH_

#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <limits>
#include <string>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/Enum.hh"
#include "General/environment.hh"
#include "MatVec/vector.hh"
#include "Optimization/Design/Filter.hh"
#include "Utils/StdVector.hh"

namespace CoupledField {
class Point;
}  // namespace CoupledField

namespace CoupledField
{
class Condition;
class DesignElement;
class DesignSpace;
class DesignStructure;
class Function;
class LevelSetElement;
class Objective;
class ResultDescription;
class SIMPElement;
class SinglePDE;
class TransferFunction;
class DesignSpace;
class MultiMaterial;
struct Elem;

/** This DesignElement package provides information about the direct neighbours for uniform cartesian
 * quadrilateral or hexahedral meshes. */
class VicinityElement
{
public:
  /** resizes design and set entries to NULL */
  VicinityElement();

  /** Initializes all elements, ignores multiple design elements!
   * Adds VicinityElement packages to the design elements.
   * Make use of of the element neighbors and works only for regular grids!
   * Makes it only once! If the first element has vicinity we return silent!
   * @param structure in the periodic case this is a helper */
  static void Init(DesignSpace* design, DesignStructure* structure);

  /** Gives the number of not NULL entries. */
  int GetNumberOfEntries() const;

  /** This are indices for the entries to design. */
  enum Neighbour {X_P = 0, X_N = 1, Y_P = 2, Y_N = 3, Z_P = 4, Z_N = 5, NONE = -1 };

  /** Helper which translates X_P and X_N to 0, Y_P and Y_N to 1, Z_P and Z_N to 2. Nothing else! */
  static int ToMainAxis(Neighbour neigh)
  {
    if(neigh == X_P || neigh == X_N) return 0;
    if(neigh == Y_P || neigh == Y_N) return 1;
    assert(neigh == Z_P || neigh == Z_N);
    return 2;
  }

  /** @axis = 0, 1, 2
   * @param dir = -1 or 1 */
  static Neighbour ToNeighbour(int axis, int dir)
  {
    assert(axis >= 0 && axis <= 2);
    assert(dir == -1 || dir == 1);
    return (Neighbour) (2 * axis + (dir == 1 ? 0 : 1));
  }

  /** Gives the neighbor elements */
  DesignElement* GetNeighbour(Neighbour idx) { return design[idx]; }

  /** Gives far away neighbors.
   * @param n = 1  is next neighbor */
  static DesignElement* GetNeighbour(DesignElement* base, Neighbour idx, int n, bool throw_exception = true);

  /** convenience to check if the neighbor is NULL */
  bool HasNeighbor(Neighbour idx) const { return design[idx] != NULL; }

  /** @see GetNeighbour(Neighbour idx, int n) */
  static bool HasNeighbor(DesignElement* base, Neighbour idx, int n);


  /** dump method for logging */
  std::string ToString() const;

  /** Contains the next neighbors (only +/- x,y(,z) and not diagonal.
   * Ordered as +x, -x, +y, -y (+z, -z). As the elements are DesignElements only within this region.
   * If there is no neighbor (e.g. z in 2D or on the boundary) the value is NULL.
   * Periodic BC could be treated as 'real' element. Please check and clear this comment!*/
  StdVector<DesignElement*> design;

  /** Has the element neighbors by periodic boundary conditions? */
  bool periodic;

private:

  /** Compares the node coordinates of two elements and decides which orientation we have.
   * Only face/line (3D/2D) neighbors are assumed. We are happy with the first nodes
   * coordinates to perform this comparison. Note that we regular!
   * @param reference from 'reference' element the first node's coordinates
   * @param other for a 'close' neighbor element the first node's coordinate
   * @param spacing to identify periodic b.c. correctly
   * @return the index within VicinitiyElement::design */
  static Neighbour FindRelativeNeighborLocation(Point& referenence, Point& other, StdVector<double> spacing);
};


/** holds the value of the topology gradient */
class TopGradElement
{
public:
  /** standard ctor */
  explicit TopGradElement(const double val = 0.0) : value(val) {}
  
  /** the value of the topgrad on this element */
  double value;
};

/** This class is introduced as a lightweight basic subset of the DesignElement
 * it does not yet correspond to a finite element
 * it is merely a class used for an optimizable variable
 * only with basic properties as designvalue, bounds, derivatives
 * no filtering or alike at all */
class BaseDesignElement
{
public:

  /** types for GetValue() and for the optimization results in the xml file */
  typedef enum { DESIGN, COST_GRADIENT, CONSTRAINT_GRADIENT, WEIGHT, OBJECTIVE, NUM_NEIGHBOURS,
    LEVEL_SET_VALUE, LEVEL_SET_STATE, TOPGRAD_VALUE, SHAPEGRAD_VALUE, SHAPEGRAD_NODE_VALUE,
    MAX_SLOPE, /* the max(abs()) of the 2 * dim slope constraints for each element */
    MAX_CHECKERBOARD, /* the max value per element */
    MAX_MOLE, /* the max mole value */
    MAX_OSCILLATION, /* the max value per element */
    MAX_JUMP, /* weak greyness constraint formulation */
    PENALIZED_STRESS, /* stess with own transfer function */
    DESIGN_TRACKING, /* (rho-rho^*)^2 but without 1/N */
    PROJECTION, /* local value from projection || nu(rho_i) - H_eta_beta(rho_i) ||^2 */
    LEVEL_SET_GRAD_XP, LEVEL_SET_GRAD_XN, LEVEL_SET_GRAD_YP, LEVEL_SET_GRAD_YN, LEVEL_SET_GRAD_ZP, LEVEL_SET_GRAD_ZN } ValueSpecifier;

    /** The type of this design element, influences the Get*Bound() methods.
     * By definition the design elements are stored in the ordering of the type!!
     * make sure, that ALL_DESIGNS is the last with the highest number!!! */
    typedef enum { UNITY = -10, NO_DERIVATIVE = -9, NO_MULTIMATERIAL = -8, TENSOR_TRACE = -7, ELAST_ALL = -6, DIELEC_TRACE = -5, DIELEC_ALL = -4, PIEZO_ALL = -3, DEFAULT = -2, NO_TYPE = -1, DENSITY = 0,
                   POLARIZATION = 1, ACOU_DENSITY = 2, EMODUL, POISSON, LAMELAMBDA, LAMEMU, EMODULISO, POISSONISO,
                   GMODUL, MASS, DAMPINGALPHA, DAMPINGBETA, TENSOR11, TENSOR22, TENSOR33, TENSOR23, TENSOR13, TENSOR12, SLACK,
                   DIELEC_11, DIELEC_12, DIELEC_22, PIEZO_11, PIEZO_12, PIEZO_13, PIEZO_21, PIEZO_22, PIEZO_23,
                   ROTANGLE, STIFF1, STIFF2,STIFF3, MULTIMATERIAL, ALL_DESIGNS} Type;

  BaseDesignElement(Type type = NO_TYPE);
  virtual ~BaseDesignElement() {};

  Type GetType() const { return type_; }

  /** Checks if test matches super. To be used in Function::SetElements()
   * if super == test it is compatible
   * @param super shall TENSOR_TRACE, ELAST_ALL, ... DEFAULT
   * @param test ege. DIELEC_11, ... */
  static bool IsCompatible(Type super, Type test);

  /** Allows to set the design element. */
  void SetDesign(double value) { this->design = value; }

  /** Return the design value.
   * In the derived DesignElement() the instance is overloaded and invalidated! */
  virtual double GetDesign() const { return(this->design); }

  /** returns the type */
  virtual std::string ToString() const;

  /** Get the gradient values for either objective or constraint.
   * if neither f nor g is given the objective gradient sum is returned */
  double GetPlainGradient(const Objective* c, const Condition* g) const;

  double GetPlainGradient(const Function* f) const;

  /** Sum app the old value (get and set together) */
  void AddGradient(const Objective* c, const Condition* g, double value);

  void AddGradient(const Function* f, double value);

  /** Reset either gradients of the class
   * @param vs either COST_GRADIENT or CONSTRAINT_GRADIENT 
   * @param g this should preferably be a Function*, but it didn't work and
   *  it is currently only needed for Condition anyways */
  void Reset(ValueSpecifier vs, Function* f = NULL);

  /**  Gets the lower bound of the design variable -
   * up to now this are defaults by type */
  double GetLowerBound() const { return lower_; }

  /** The upper bound of the design variable for the optimizer */
  double GetUpperBound() const { return upper_; }

  /** Set the lower bound of the design variable */
  void SetLowerBound(const double v) { lower_ = v; }

  /** Set the upper bound of the design variable */
  void SetUpperBound(const double v) { upper_ = v; }

  /** adjusts length of the gradient vectors possibly not known during creation */
  void PostInit(int objectives, int constraints);

  static Enum<Type> type;

protected:

  /** The scalar value. Public access only via getter to handle filtering. */
  double design;

  /** Sums up the costGradient values (they include penalty) */
  double SumObjectiveGradient() const;

  /** <p>The gradient of the constraint function w.r.t. this element of the design space.
   * Every constraint contains an unique index attribute (which is the order in the xml file)
   * only for the purpose to index this vector.</p>
   * <p>Therefore this vector has to be initialized on runtime</p> */
  StdVector<double> constraintGradient;

  /** For multiple objective functions. It already includes penalty!
   * @see constraintGradient */
  StdVector<double> costGradient;

  /** The lower bound of this design variable. Redundant but such faster than look it up */
  double lower_;

  double upper_;

  /** what is our design type */
  Type type_;

};


/** This is the base design element, storing all information related to 
 * elements. For ersatz material this is a 1:1 relation to finite elements.
 * Method specific information (SIMP, Level-Set, FreeMat, ...) is in add-ons. */
class DesignElement : public BaseDesignElement
{     
public:


  /** This defines how to access variables (design, objective_gradient, ...),
   *  PLAIN is the value and SMART does a filtering if enabled otherwise also as PLAIN */
  typedef enum { PLAIN, SMART } Access;

  /** The empty constructor is the StdVector and for ghost elements */
  DesignElement();

  /** This sets the DesignElement with the values from the XML file.
   * Is slow as it does the same evaluation often but is only O(n)
   * @param space to output 'penalizedDesign' the pointer is needed to find the transfer function
   * @param index location within the design space */
  DesignElement(Type dt, double lower, double upper, Elem* elem, unsigned int index, MultiMaterial* mm);

  /** Dummy elements for Function */
  DesignElement(Elem* elem, Type type, unsigned int index, int pseudoElementIndex);

  virtual ~DesignElement();

   /** We might need the transfer functions! */
  static void SetDesignSpace(DesignSpace* space)
  {
    space_ = space;
  }

  static DesignSpace* GetDesignSpace()
  {
    return space_;
  }

  /** Default mapping from the PDE */
  static Type Default(const SinglePDE* pde);

    /** This specifies result details for various ValueSpecifier/Detail combinations:
     * OBJECTIVE/SYMMETRY (check!)
     * COST_GRADIENT/FINITE_DIFF_COST_GRADIENT
     * COST_GRADIENT/ERROR_COST_GRADIENT
     * For the PiezoSIMP case:
     * COST_GRADIENT/MECH_MECH ... as ErsatzMaterial::CalcU1KU2() Parameters. QUAD is calcMode
     * or it does the very special purpose symmetry for the ValueSpecifier */
    typedef enum { NONE, SYMMETRY, FINITE_DIFF_COST_GRADIENT, ERROR_COST_GRADIENT,
      MECH_MECH, ELEC_ELEC, ELEC_ELEC_QUAD, ELEC_MECH, MECH_ELEC,
      COMPLIANCE, VOLUME, PENALIZED_VOLUME, GAP, TRACKING, HOMOGENIZATION_TRACKING,
      POISSONS_RATIO, YOUNGS_MODULUS, YOUNGS_MODULUS_E1, YOUNGS_MODULUS_E2,
      TYCHONOFF, GREYNESS, REALVOLUME,
      GLOBAL_SLOPE, GLOBAL_CHECKERBOARD, STRESS,
      /*!< only for the projection function. This is the element wise fake filter part */
      PROJECTION_FILTER } Detail;

    /** Gets the design element
     * @param access if plain the rho value if SMART and filtering is enabled the filtered value */
    double GetDesign(Access access) const;

    /** Gives the physical design, which is penalized and filtered if we have density filtering.
     * Therefore there is no access as we are implicit SMART */
    double GetPhysicalDesign(bool densForElec = false) const;
    
    /** Return whether physical design is reasonable for this DesignElement::Type */
    bool HasPhysicalDesign() const;

    /** Overloads the original BaseDesignElement() method and invalidates it to force the access version */
    double GetDesign() const;

    /** Checks out specialResult[]!
     * @param dofs 1 for scalar values */
    void GetValue(ResultDescription& rd, StdVector<double>& out, unsigned int dofs) const;

    /** This method decides if either GetFilteredValue() or GetPlainValue() is to be returned.
     * @param g mandatory for vs = CONSTRAINT_GRADIENT only */
    double GetValue(ValueSpecifier vs, Access access, Condition* g = NULL) const;

    /** internal helper to get the value by type
     * @param g for sp = CONSTRAINT_GRADIENT only */
    double GetPlainValue(ValueSpecifier valueSpecifier, Condition* g = NULL) const;

    /** This is only for the Heaviside Filter!! as is so often called there that it makes a real difference! */
    double GetPlainDesignValue() const { return design; }

    /** Initilize the Enum. Currently called by Optimization::CreateInstance() */
    void static SetEnums();
    
    /** Write key values as attributes
     * @param tf if given prints the physical lower bound */
    void ToInfo(PtrParamNode in, TransferFunction* tf) const;

    /** @see BaseDesignElement::ToString() */
    std::string ToString() const { return ToString(this); }

    /** makes a short dump, handles NULL */
    static std::string ToString(const DesignElement* de);
    
    /** helper for LOG output */
    static std::string ToString(const StdVector<DesignElement*>& vec, bool print_type = false);
    static std::string ToString(const StdVector<DesignElement>& vec, bool print_val = false, bool print_type = false);

    /** Calculates the volume of the element, used static helpers.
     * caches the result, hence cheap to query again */
    double CalcVolume();

    /** to make the class polymorphi and we can dynamic_cast<> it */
    /** Pointer to the element of the region, parameter for integration, ... */
    Elem*  elem;

    /** The index of this element within the design space - 0 based */
    unsigned int GetIndex() const { assert(index_ != std::numeric_limits<unsigned int>::max()); return index_; }

    /** In case we are a pseudo design element which is not within the design domain but from the
     *  non-design region of a function (e.g. stress) this index stores the index within the element storage.
     *  It is -1 if we are a standard design element, otherwise it is >= the number of standard design elements.
     *  This element is not stored within data but with the registered pseudoDesigns_ and might be non-unique */
    int GetPseudoElementIndex() const { return pseudoElementIndex_; }

    /** Gets the index within the element vector solutions store in the solutions space.
     * Handles multiple design elements and pesudo elements */
    unsigned int GetElementSolutionIndex() const;

    /** extracts the OPT_RESULT_x or -1 if some other solution type */
    static int GetOptResultIndex(SolutionType st);

    /** This stores special values during calculation for visualization via
     * result description in XML: See DesignSpace::GetSpecialResultIndex()
     * <result id="optResult_2" design="density" access="plain" value="costGradient" />
     * <result id="optResult_3" design="density" access="plain" value="objective" /> */
    Vector<double> specialResult;

    static Enum<ValueSpecifier> valueSpecifier;

    static Enum<Access> access;

    static Enum<Detail> detail;

    /** This are our add-ons, not all are initialized! */
    SIMPElement* simp;

    /** The vicinity (structured grids, e.g. for Level-Set */
    VicinityElement* vicinity;

    /** The level-set element, will be destroyed by LevelSet */
    LevelSetElement* lse_;
    
    /** The topgrad element, will be destroyed by TopGrad */
    TopGradElement *tge;

    /** if we are a multimaterial this is our material
     * and the index there is our own index*/
    MultiMaterial* multimaterial;

    /** calculates the location on request and stores it */
    Point* GetLocation();

private:

  /** Here we share the base initialization for the constructors */
  void Init();

  /** returns the non-scalar values */
  void GetValue(ValueSpecifier sp, StdVector<double>& out) const;  
  
  /** the barycenter of this element only set on request. */
  Point* location_;

  /** @see GetIndex() */
  unsigned int index_;

  /** @see GetPseudoElementIndex() */
  int pseudoElementIndex_;

  /** the element volume calculated on request by CalcVolume() */
  double elemVol_;

  /** up to now only needed to extract 'penalizedDesign'. Make it protected
   * if you need it. */
  static DesignSpace* space_;
};


/** This is a add-on to DesignElement and contains all data specific for SIMP,
 * this is especiall the filtering stuff */
class SIMPElement
{
public:
  /** Base init of element */
  SIMPElement(DesignElement* base);

  /** Does sensitvity filtering
   * @param g @see GetPlainValue() */
  double GetSensitivityFilteredValue(DesignElement::ValueSpecifier valueSpecifier, Condition* g) const;

  /** Does design filtering. */
  double GetDensityFilteredValue(DesignElement::ValueSpecifier sp, Filter::Density fd) const;

  /** Helper for GetDensityFilteredValue() */
  double CalcHeaviside(double input_value) const;

  /** Calculates the tanh function. This is a variant of the Xu-Filter, see also Wang/Laraow/Sigmund;2010 */
  double CalcTanh(double input_value) const;

  /** only for sensitivities for density filtering.
   * See Sigmund; Morpology-based black and white filters for topology optimization; 2007; (35) and (36)
   * @param sp COST_GRADIENT, CONSTRAINT_GRADIENT or DENSITY for PROJECTION only */
  double GetDensityFilteredGradient(DesignElement::ValueSpecifier sp, Condition* g) const;

  /** Sums up the weights of the neighbors and optionally the own element */
  double CalcWeightSum(bool include_this) const;

  /** Neighborhood is element and pre-calculated distance */
  struct NeighbourElement
  {
  public:
    /** read the variable */
    DesignElement* neighbour;

    /** pre-calculated weight: radius - distanance and >= 0 */
    double        weight;

    /** the distance in domain dimensions! */
    double        distance;
  };

  /** The weight of THIS element which is radius */
  double weight;

  /** The neighbors if filter otherwise empty.
   * The element itself is NOT part of the neighborhood!
   * @see DesignStructure::DesignStructure() */
  StdVector<NeighbourElement> neighborhood;

  /** string representation for logging, includes neighborhood.
   * @param level 0 is elements, 1 is with weighting and distance */
  std::string ToString(int level = 0) const;

  /** for debugging. Sums the weights of all neighbors, ... */
  void Dump();

  /** The complete filter settings */
  Filter filter;

private:

  /** We need our base design element to do the filtering */
  DesignElement* de_;
};

/** requires in DesignSpace and DesignMaterial */
class DesignID
{
public:
   DesignID(DesignElement::Type design = DesignElement::NO_TYPE, MultiMaterial* mm = NULL)
   {
     this->design = design;
     this->multimaterial = mm;
   }

   DesignElement::Type design;
   /** index. -1 for non-multimaterial */
   MultiMaterial*      multimaterial;
};


/** implemented in StdVector.cc, there we need it */
std::ostream & operator << ( std::ostream & out, const DesignID& id);


/** <p>A result description holds the result element in the xml file which describes what data from
 * a DesignElement is to be written to the cfs output. The following parameters have to be given
 * in the optimization elment: <result id="optResult_1" design="density" access="plain" value="costGradient"/>
 * This is referenced by the solution type optResult_1/2/3.</p> */
class ResultDescription
{
public:
  /** empty destructor for StdVector() and default in DesignSpace::ExtractResult().
   * Does NOT initialize all elements! */
  ResultDescription();

  /** reads itsef from the xml element.
   * @param pn our data */
  ResultDescription(PtrParamNode pn);

  SolutionType solutionType;

  /** Finds the proper design element by element number */
  DesignElement::Type design;

  /** optionally filtered or plain */
  DesignElement::Access access;

  /** Which of the values? */
  DesignElement::ValueSpecifier value;

  /** An optional detail for values COST_GRADIENT and OBJECTIVE in PiezoSIMP case */
  DesignElement::Detail detail;

  /** An optional excitation label */
  std::string excitation;
};

inline
double SIMPElement::CalcHeaviside(double input_value) const
{
  const Filter* f = &de_->simp->filter;
  assert(f->type_ == Filter::DENSITY);
  assert(f->density_ == Filter::HEAVISIDE || f->density_ == Filter::MOD_HEAVISIDE);

  double result;

  double b = f->GetBeta();
  assert(b >= 0.0 && b < 2000);

  if(f->density_ == Filter::HEAVISIDE)
  {
    // we apply the correction factor in a way that H(rho_min) = rho_min and H(1) = 1
    double corr = (1.0 - (1.0 - input_value) * f->heaviside_corr) * input_value;
    result = 1.0 - std::exp(-1.0 * b * corr) + corr * std::exp(-1.0 * b);

    // no LOG_DBG possible due to inline
    // std::cout << "CH: de=" << de_->elem->elemNum << " f=" << f->density.ToString(f->density_)
    ///          << " hc=" << f->heaviside_corr << " corr=" << corr << " iv=" << input_value << " -> " << result << std::endl;
  }
  else // if(f->density_ == Filter::MOD_HEAVISIDE)
  {
    // make sure we are within the bounds
    double ub = this->de_->GetUpperBound();
    double lb = f->GetLowerBound(this->de_); // might be force_lower_bound from the filter setting

    double first    = std::exp(-1.0 * b * (1.0 - input_value));
    double second   = -1.0 * (1.0 - input_value) * std::exp(-1.0 * b);

    assert((ub-lb) > 1e-2); // if not you probably forgot to set force_lower_bound in the filter definition

    result = (ub-lb) * (first + second) + lb;

    // std::cout << "CH: el=" << de_->elem->elemNum << " iv=" << input_value << " b=" << b << " lb=" << lb << " (ub-lb)=" << (ub-lb)
    //           << " +1st=" << first << " +2nd=" << second << " -> " << result << std::endl;
  }

  return result;
}

} // end of namespace

#endif /*DESIGNELEMENT_HH_*/
