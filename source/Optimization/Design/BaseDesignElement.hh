#ifndef BASEDESIGNELEMENT_HH_
#define BASEDESIGNELEMENT_HH_

#include "General/Enum.hh"
#include "General/Environment.hh"
#include "Utils/StdVector.hh"
#include "MatVec/Vector.hh"


/** We try to make this include maximal lightweight to prevent cyclic inclusions.
 * The implementation is in DesignElement.cc */

namespace CoupledField
{
class Function;
class Condition;
class Objective;

/** This class is introduced as a lightweight basic subset of the DesignElement
 * it does not yet correspond to a finite element
 * it is merely a class used for an optimizable variable
 * only with basic properties as designvalue, bounds, derivatives
 * no filtering or alike at all */
class BaseDesignElement
{
public:

  /** types for GetValue() and for the optimization results in the xml file.
   * The numbered values can also be applied to the python exported function opt_get_value */
  typedef enum { DESIGN = 0, COST_GRADIENT = 1, CONSTRAINT_GRADIENT = 2, WEIGHT = 3, NUM_NEIGHBOURS = 4,
    FUNCTION_GRADIENT, /* combines COST_GRADIENT and CONSTRAINT_GRADIENT */
    OBJECTIVE, LEVEL_SET_VALUE, LEVEL_SET_STATE, TOPGRAD_VALUE, SHAPEGRAD_VALUE, SHAPEGRAD_NODE_VALUE,
    MAX_SLOPE, /* the max(abs()) of the 2 * dim slope constraints for each element */
    MAX_CHECKERBOARD, /* the max value per element */
    MAX_MOLE, /* the max mole value */
    MAX_OSCILLATION, /* the max value per element */
    MAX_JUMP, /* weak greyness constraint formulation */
    QUADRATIC_VM_STRESS, /* element wise <Bu, M Bu> with transfer function for stress. Note that the stess constraints are evaluated at the integration points an the cfs pde results at the center only */
    LOCAL_LOAD_FACTOR, /* local buckling load factor (e.g. for two scale structures) */
    DESIGN_TRACKING, /* (rho-rho^*)^2 but without 1/N */
    HEAT_NODAL_TRACK_VAL, /* for each node: (stateSol - trackVal)^2 */
    TEMP_AT_INTERFACE, /* temperature at interface between solid and void */
    PROJECTION, /* local value from projection || nu(rho_i) - H_eta_beta(rho_i) ||^2 */
    LEVEL_SET_GRAD_XP, LEVEL_SET_GRAD_XN, LEVEL_SET_GRAD_YP, LEVEL_SET_GRAD_YN, LEVEL_SET_GRAD_ZP, LEVEL_SET_GRAD_ZN,
    SHAPE_MAP_GRAD, /* the sum of all dtanh_da over all ip for a rho element for shape mapping */
    SHAPE_MAP_ORDER, /* the number of integration points for this element */
    SHAPE_MAP_CORNER, /* the difference between the minimal and maximal corner values (min and max for all shapes) Makes only sense for 1 shape!*/
    MMA_ASYMPTOTE, /* element and function wise MMA asymptotes for own MMA implementation */
    MMA_LOWER_VAL ,/* own MMA implementation lower bound for design variable */
    MMA_UPPER_VAL, /* own MMA implementation lower bound for design variable */
    MMA_OBJ_GRADIANT, /* own MMA implementation gradiant of objective */
    MMA_CON_GRADIANT_1, /* own MMA implementation gradiant of constraint */
    MMA_CON_GRADIANT_2, /* own MMA implementation gradiant of constraint */
    SPLINE_BOX_GRAD_X, /* gradient (sum over all ip w.r.t. control point movement in x direction) for a rho element for spline box mapping */
    SPLINE_BOX_GRAD_Y, /* gradient (sum over all ip w.r.t. control point movement in y direction) for a rho element for spline box mapping */
    SPLINE_BOX_GRAD_Z, /* gradient (sum over all ip w.r.t. control point movement in z direction) for a rho element for spline box mapping */
    SPLINE_BOX_INT_ORDER, /* the number of integration points for this element */
    SPLINE_BOX_INT_CORNER, /* the difference between the minimal and maximal corner values*/
    GENERIC_ELEM, /** with attribute generic with content from .info.xml, e.g. for Python spaghetti */
    FILTERED_DESIGN,
    DIFF_FILTERED_DESIGN, /** x - F*x, where F is filter matrix and x is the design vector */
  } ValueSpecifier;

  /** The type of this design element, influences the Get*Bound() methods.
   * By definition the design elements are stored in the ordering of the type!!
   * make sure, that ALL_DESIGNS is the last with the highest number!!! */
  typedef enum { UNITY = -20, // unused stuff from ShapeOpt
                 // wrapper design types representing a class of design types
                 SHAPE_MAP, // NODE or PROFILE
                 SPAGHETTI, // NODE, PROFILE, NORMAL
                 SPLINE_BOX, // e.g. CP
                 MECH_TRACE, // MECH_11, MECH_22, MECH_33
                 MECH_ALL, DIELEC_TRACE, DIELEC_ALL, PIEZO_ALL,
                 // special types
                 MULTIMATERIAL, NO_MULTIMATERIAL, INTERPOLATION,
                 NO_DERIVATIVE = -3, DEFAULT = -2, NO_TYPE = -1,
                 // real design types
                 DENSITY = 0, POLARIZATION, EMODUL, POISSON, LAMELAMBDA, LAMEMU, EMODULISO, POISSONISO,
                 GMODUL, MASS, DAMPINGALPHA, DAMPINGBETA = 12,
                 MECH_11, MECH_12, MECH_13, MECH_14, MECH_15, MECH_16,
                 MECH_22, MECH_23, MECH_24, MECH_25, MECH_26,
                 MECH_33, MECH_34, MECH_35, MECH_36,
                 MECH_44, MECH_45, MECH_46,
                 MECH_55, MECH_56,
                 MECH_66,
                 RHS_DENSITY, // for mag opt, e.g. coil modelling (scaling current)
                 DIELEC_11, DIELEC_12, DIELEC_22, PIEZO_11, PIEZO_12, PIEZO_13, PIEZO_21, PIEZO_22, PIEZO_23,
                 ROTANGLE, SHEAR1, STIFF1, STIFF2, STIFF3, LOWER_EIG_BOUND, ROTANGLEFIRST, ROTANGLESECOND, ROTANGLETHIRD,
                 SLACK, ALPHA,  // slack variables or spaghetti scaling
                 NODE, PROFILE, // shape mapping and spaghetti
                 NORMAL, RADIUS, // spaghetti height and radius
                 CP,            // spline box control point
                 FEATURE_MAPPING_PX, FEATURE_MAPPING_PY, FEATURE_MAPPING_QX, FEATURE_MAPPING_QY, FEATURE_MAPPING_P, // Pill::GradDistance() 
                 ALL_DESIGNS } Type; // ALL_DESIGNS needs to be last

  /** This defines how to access variables (design, objective_gradient, ...),
   *  PLAIN is the value and SMART does a filtering if enabled otherwise also as PLAIN */
  typedef enum { PLAIN = 0, SMART = 1} Access; // not used here but needed for virtual method GetDesign(Access)

  BaseDesignElement(Type type = NO_TYPE);
  virtual ~BaseDesignElement() {};

  Type GetType() const { return type_; }

  /** Checks if test matches super. To be used in Function::SetElements()
   * if super == test it is compatible
   * @param super shall TENSOR_TRACE, ELAST_ALL, ... DEFAULT
   * @param test ege. DIELEC_11, ... */
  static bool IsCompatible(Type super, Type test);

  /** checks if the type is a shape mapping type. Then there is a counterpart with the same name in ShapeMapDesign::Type
   * @see ShapeMapDesign::Convert() */
  static bool IsShapeMapType(Type type) { return type == NODE || type == PROFILE || type == SHAPE_MAP; }

  static bool IsSpaghettiType(Type type) { return type == NODE || type == PROFILE || type == NORMAL || type == SPAGHETTI; }

  /** checks if the type is a splinebox feature mapping type. */
  static bool IsSplineBoxType(Type type) { return type == SPLINE_BOX || type == CP; }

  static bool IsFeatureMappingType(Type type) { return IsShapeMapType(type) || IsSpaghettiType(type) || IsSplineBoxType(type); }

  /** maps from SolutionType to DesignElement::Type. In the PHYSICAL_* case to the standard case */
  static Type MapSolutionType(SolutionType soltype, bool throw_exception = true);

  /** checks if this is a physical solution type where we filter and penalize */
  static bool IsPhysical(SolutionType soltype);

  /** Allows to set the design element. */
  void SetDesign(double value) { this->design = value; }

  /** Return the design value.
   * In the derived DesignElement() the instance is overloaded and invalidated! */
  virtual double GetDesign() const { return(this->design); }

  /** Overwritten in DesignElement.
   * Note the difference between BaseDesignElement::Access and Function::Access.
   * Due to C++ and cyclic inclusions we cannot have the Function::Access option here :(
   * @see DesignElement::GetDesign(Function*) and DesignElement::GetPhysicalDesign() */
  virtual double GetDesign(BaseDesignElement::Access access) const { return GetDesign(); };

  /** The index of this element within the design space - 0 based.
   * @see GetOptIndex()*/
  unsigned int GetIndex() const { assert(index_ != std::numeric_limits<unsigned int>::max()); return index_; }

  /** When only part of the design space is really for optimization (e.g. ShapeMapDesign with symmetry) the design exported
   * to the optimizers need their own consecutive index returned by GetOptIndex() -> overloaded in ShapeParamElement */
  virtual unsigned int GetOptIndex() const { return GetIndex(); }

  /** returns the type */
  virtual std::string ToString() const;

  /** Set porosity per element, currently only implemented for two-scale volume */
  virtual void SetElemPorosity (double vol) {EXCEPTION("Not implemented");};

  /** Get porosity per element, currently only implemented for two-scale volume */
  virtual double GetElemPorosity () {EXCEPTION("Not implemented");};

  /** This is only for the Heaviside Filter!! as is so often called there that it makes a real difference! */
  double GetPlainDesignValue() const { return design; }

  /** for f = objective gives the value by index. For multiple objective use SumCostGradient() */
  double GetPlainGradient(const Function* f) const;

  /** return from the function. for multi objective use SumCostGradient() */
  double GetPlainGradient(const Objective* c) const;
  double GetPlainGradient(const Condition* g) const;

  /** Sum app the old value (get and set together) */
  void AddGradient(const Objective* c, const Condition* g, double value);

  void AddGradient(const Function* f, double value);

  /** Reset either gradients of the class
   * @param vs either COST_GRADIEN, CONSTRAINT_GRADIENT or FUNCTION_GRADIENT (needs f given!)
   * @param g this should preferably be a Function*, but it didn't work and
   *  it is currently only needed for Condition anyways */
  void Reset(ValueSpecifier vs, Function* f = nullptr);

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
  void PostInit(int objectives, int constraints)
  {
    // resize and init with 0.0 so constraint, which only act on one design variable,
    // do not have to set all others explicitly to zero
    costGradient.Resize(objectives, 0.0);
    constraintGradient.Resize(constraints, 0.0);
  }

  /** helper for LOG output */
  static std::string ToString(const StdVector<BaseDesignElement*>& vec, bool print_type = false);

  static Enum<Type> type;

  /** <p>The gradient of the constraint function w.r.t. this element of the design space.
   * Every constraint contains an unique index attribute (which is the order in the xml file)
   * only for the purpose to index this vector.</p>
   * <p>Therefore this vector has to be initialized on runtime</p> */
  StdVector<double> constraintGradient;

  /** For multiple objective functions. It already includes penalty!
   * Don't access directly with the exception of ShapeMapDesign!
   * @see constraintGradient */
  StdVector<double> costGradient;

  /** Sums up the costGradient values (they include penalty) */
  double SumObjectiveGradient() const;

  /** for each node: grad = 4* ds/drho_filt * drho_filt/drho * (1-2*s)
   * s is the interpolated density at a node, e.g. s = 1/4*(rho_1+rho_2+rho_3+rho_4)
   * Need to store this in order to calculate the right derivative with filtering and penalization*/
  Vector<double> interfaceDrivenLoadGrad_;

protected:

  /** The scalar value. Public access only via getter to handle filtering. */
  double design;

  /** The lower bound of this design variable. Redundant but faster than look it up */
  double lower_;

  double upper_;

  /** what is our design type */
  Type type_;

protected:
  /** @see GetIndex() */
  unsigned int index_;

};

} // end of namespace

#endif /*BASEDESIGNELEMENT_HH_*/
