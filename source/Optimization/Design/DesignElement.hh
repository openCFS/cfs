#ifndef DESIGNELEMENT_HH_
#define DESIGNELEMENT_HH_

#include <assert.h>
#include <cmath>
#include <stddef.h>
#include <limits>
#include <string>

#include "Optimization/Design/BaseDesignElement.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/Environment.hh"
#include "Optimization/Design/Filter.hh"

namespace CoupledField
{
class Point;
class Excitation;
class DesignElement;
class DesignSpace;
class DesignStructure;
class LevelSetElement;
class Objective;
class ResultDescription;
class SIMPElement;
class SinglePDE;
class TransferFunction;
class DesignSpace;
struct MultiMaterial;
class ErsatzMaterial;
class Context;
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
  enum Neighbour {X_P = 0, X_N = 1, Y_P = 2, Y_N = 3, Z_P = 4, Z_N = 5, XY_P = 6, XY_N = 7, NONE = -1 };

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

/** Bastian's shape optimization */
class ShapeDesignElement : public BaseDesignElement
{
public:
  /** ShapeDesignElement have an index, needed for sparse gradients, i.e. shape constraints */
  ShapeDesignElement(unsigned int index);
};

/** for FeaturedDesign. Holds a shape parameter. E.g. node in shape mapping with dof=x ny+1 times with ny is number of elements */
class ShapeParamElement : public BaseDesignElement
{
public:
  ShapeParamElement(Type type = NO_TYPE, unsigned int index = std::numeric_limits<unsigned int>::max());

  virtual ~ShapeParamElement() {};

  /** at the moment we don't have filter for shape mapping but the slope constraints are implemented using GetDesign(SMART).
   * So just ignore it! */
  virtual double GetDesign(BaseDesignElement::Access access) const { return design; }

/** to handle exporting only parts of the design (due to symmetry) to the external optimizers */
  virtual unsigned int GetOptIndex() const { assert(opt_index_ != std::numeric_limits<unsigned int>::max()); return opt_index_; }

  void SetOptIndex(unsigned int idx) { this->opt_index_ = idx; }

  /** overwrite to add opt_idx */
  virtual std::string ToString() const;

  /** gives a human readable label without values, e.g. for WriteGradFile() */
  virtual std::string GetLabel() const;

  /** The dof for shape elements. This is the design variable, the other coordinates are implicitly given be the mesh.*/
  typedef enum { NOT_SET = -1, X=0, Y=1, Z=2 } Dof; // X=0 to Z=2 must not be changed, it is index to ShapeMapDesign::n_

  static Enum<Dof> dof;

  /** for node which dof BaseDesignElement::value is for. value correspond to the missing entry in coord and idx*/
  Dof dof_ = NOT_SET;

private:
  /** see BaseDesignElement::GetOptIndex() */
  unsigned int opt_index_ = std::numeric_limits<unsigned int>::max();
};


/** Base Variable for FeatureMappingDesign and SpaghettiDesign. The later extends */
class FeatureVariable : public ShapeParamElement
{   
public:
  typedef enum { NO_TIP = -1, START=0, END=1, INNER=2 } Tip; // INNER only for SpaghettiDesign or BAR for Feature::Distance()

  static Enum<Tip> tip_enum;

  /** optional late constructor - only SpaghettiDesign */
  void InitInnerNode(int noodle_idx, Dof dof, double val, double lower, double upper);

  static std::string ToString(const StdVector<FeatureVariable>& vec, bool show_fixed = true);

  static Vector<double> AsVector(const StdVector<FeatureVariable>& vec);

  bool IsVariable() const { return !fixed && map == ""; }

  /** check all components of point if we are fixed. Mapping is OK */
  static bool IsFixed(const StdVector<FeatureVariable>& point);

  /** Count number of real variables. Exclude fixed and mapped */
  static int CountRealVariables(const StdVector<FeatureVariable>& point);

  static void CompareToInfoHelper(const FeatureVariable* v0, const FeatureVariable* test, std::string& fixed, std::string& lower, std::string& upper);

  void Parse(PtrParamNode pn, int feature, double interpolate_value = -12.34);

  /** all other BaseDesignElement children do this in the constructor.
   * @see SetOptIndex() */
  void SetIndex(int idx) { index_ = idx; }

  void ToInfo(PtrParamNode in) const;

  /** for debug info */
  std::string ToString() const override;

  std::string GetLabel() const override;

  /** does not apply for all */
  Tip tip = NO_TIP;

  bool fixed = false;

  /** Pill has mapping and this is the reference, map refers to key   */
  std::string key; 

  /** Pill has mapping and we point to a key. Check and housekeeping is not done by FeatureVariable
   * Consider to implement formulas (-key, key+5, ...) */
  std::string map; 

  FeatureVariable* map_to = nullptr; // only for FeatureMappingDesign. Make sure to point to non-resizable vector to avoid dangling!

  /** pill/noodle index */
  int feature = -1;
};


/** This is the base design element, storing all information related to
 * elements. For ersatz material this is a 1:1 relation to finite elements.
 * Method specific information (SIMP, Level-Set, FreeMat, ...) is in add-ons. */
class DesignElement : public BaseDesignElement
{
public:

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
  static void SetDesignSpace(DesignSpace* space) { space_ = space; }

  static DesignSpace* GetDesignSpace() { return space_; }

  /** Default mapping from the PDE */
  static Type Default(const Context* ctxt);

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
      GLOBAL_SLOPE, GLOBAL_DESIGN, GLOBAL_CHECKERBOARD, STRESS,
      /*!< only for the projection function. This is the element wise fake filter part */
      PROJECTION_FILTER,
      SM_NODE, /*!< for 2D shape map */
      SM_NODE_A, /*!< for 3D shape map */
      SM_NODE_B, /*!< for 3D shape map */
      SM_PROFILE,  /*!< for shape map */
      SP_CP, /*!< for spline box mapping */
      BUCKLING_LOAD_FACTOR, /*!< macroscopic load factor */
      LOCAL_BUCKLING_LOAD_FACTOR, /*!< microscopic load factor */
      GLOBAL_BUCKLING_LOAD_FACTOR, /*!< globalized microscopic load factor */
      SIN,  /*!< sine value, e.g. for filteredDesign: F*sin(x)*/
      COS,  /*!< cosine value, e.g. for filteredDesign: F*cos(x)*/
    } Detail;

    /** Gets the design element
     * @param access if plain the rho value if SMART and filtering is enabled the filtered value */
    double GetDesign(BaseDesignElement::Access access) const;

    /** Returns the design based on Function::access. This needs to be PHYSICAL, PLAIN or FILTERED.
     * @param tf only required if PHYSICAL, if PHYSICAL and not given tries GetPyhsicalDesign() */
    double GetDesign(const Function* f, const TransferFunction* tf = NULL) const;

    /** Gives the physical design, which is penalized and filtered if we have density filtering.
     * Therefore there is no access as we are implicit SMART */
    double GetPhysicalDesign(const Context* ctxt = NULL) const;

    /** Return whether physical design is reasonable for this DesignElement::Type */
    bool HasPhysicalDesign() const;

    /** Overloads the original BaseDesignElement() method and invalidates it to force the access version */
    double GetDesign() const;

    /** Checks out specialResult[]!
     * @param dofs 1 for scalar values */
    void GetValue(ResultDescription& rd, StdVector<double>& out, unsigned int dofs) const;

    /** This method decides if either GetFilteredValue() or GetPlainValue() is to be returned.
     * @param f mandatory for vs = CONSTRAINT_GRADIENT and to determine if we are filtered! */
    double GetValue(ValueSpecifier vs, Access access, Function* f = nullptr) const;

    /** internal helper to get the value by type
     * @param g for sp = CONSTRAINT_GRADIENT only.
     * We are not virtual but we need the keyword to make the inline implementation link */
    virtual double GetPlainValue(ValueSpecifier valueSpecifier, Condition* g = nullptr) const;

    /** internal helper to get the cost gradient
     *  handles single and multi objective */
    double GetPlainCostGradient() const;

    /** returns the sum of all design elements compatible to MECH_TRACE */
    double GetPlainMechTrace() const;

    /** Initilize the Enum. Currently called by Optimization::CreateInstance() */
    void static SetEnums();

    /** Write key values as attributes
     * @param tf if given prints the physical lower bound */
    void ToInfo(PtrParamNode in, TransferFunction* tf, ErsatzMaterial* em) const;

    /** @see BaseDesignElement::ToString() */
    std::string ToString(bool barycenter = false) const { return ToString(this, barycenter); }

    /** makes a short dump, handles NULL */
    static std::string ToString(const DesignElement* de, bool barycenter = false);

    /** helper for LOG output */
    static std::string ToString(const StdVector<DesignElement*>& vec, bool print_type = false);
    static std::string ToString(const StdVector<DesignElement>& vec, bool print_val = false, bool print_type = false);

    /** Calculates the volume of the element, used static helpers.
     * caches the result, hence cheap to query again */
    double CalcVolume();

    /** Sets porosity per element, currently only implemented for two-scale volume*/
    virtual void SetElemPorosity(double vol);

    /** Get porosity per element, currently only implemented for two-scale volume*/
    virtual double GetElemPorosity();

    /** to make the class polymorphi and we can dynamic_cast<> it */
    /** Pointer to the element of the region, parameter for integration, ... */
    Elem* elem = NULL;

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
    SIMPElement* simp = NULL;

    /** The vicinity (structured grids, e.g. for Level-Set */
    VicinityElement* vicinity = NULL;

    /** The level-set element, will be destroyed by LevelSet */
    LevelSetElement* lse_ = NULL;

    /** The topgrad element, will be destroyed by TopGrad */
    TopGradElement *tge = NULL;

    /** if we are a multimaterial this is our material
     * and the index there is our own index*/
    MultiMaterial* multimaterial = NULL;

    /** calculates the location on request and stores it */
    const Point* GetLocation();

private:

  /** Here we share the base initialization for the constructors */
  void Init();

  /** returns the non-scalar values */
  void GetValue(ValueSpecifier sp, StdVector<double>& out) const;

  /** the barycenter of this element only set on request. */
  Point* location_;

  /** @see GetPseudoElementIndex() */
  int pseudoElementIndex_;

  /** the element volume calculated on request by CalcVolume() */
  double elemVol_;

  /** the element porosity calculated by eg. CalcTwoScaleVolume() */
  double elemPorosity_;

  /** up to now only needed to extract 'penalizedDesign'. Make it protected
   * if you need it. */
  static DesignSpace* space_;
};


/** This is a add-on to DesignElement and contains all data specific for SIMP,
 * this is especially the filtering stuff */
class SIMPElement
{
public:
  /** Base init of element */
  SIMPElement(DesignElement* base);

  /** Does sensitvity filtering
   * @param g @see GetPlainValue() */
  double GetSensitivityFilteredValue(DesignElement::ValueSpecifier valueSpecifier, Function* g) const;

  /** Does design filtering. */
  double GetDensityFilteredValue(DesignElement::ValueSpecifier sp, Filter::Density fd) const;

  /** Helper for GetDensityFilteredValue()
   * @param filter_idx to handle robust. 0 shall work in the non-robust case */
  static double CalcHeaviside(double input_value, const GlobalFilter* global);

  /** Calculates the tanh function. This is a variant of the Xu-Filter, see also Wang/Laraow/Sigmund;2010
   * @param filter_idx see CalcHeaviside()*/
  static double CalcTanh(double input_value, const GlobalFilter* global);

  /** Calculates the material filter of Lukas. We have function mainly for CalcNonLinFilter for GlobalFilter::SetNonLinCorrection() */
  static double CalcMaterial(double filter_value, double plain_value, const GlobalFilter* global);

  /** this is a convenience cover function for CalcHeaviside(), CalcTanh() and CalcMaterial() for SetNonLinCorrection
   * @param plain_value only for CalcMaterial*/
  static double CalcNonLinFilter(double input_value, const GlobalFilter* global, double plain_value=-1);

  /** only for sensitivities for density filtering.
   * See Sigmund; Morpology-based black and white filters for topology optimization; 2007; (35) and (36)
   * @param sp COST_GRADIENT, CONSTRAINT_GRADIENT or DENSITY for PROJECTION only */
  double GetDensityFilteredGradient(DesignElement::ValueSpecifier sp, Function* func) const;

  /** gives the filter index from DesignSpace::Context::Excitation.
   * Currently assumes multiple filter because of robust optimization.
   * Does not handle the case for multiple designs, for that one would need the design type as argument.
   * A combination multiple designs and robust would be tricky, but the question is, if this would make sense at all */
  unsigned int DetermineFilterIndex() const;

  /** Of the inlined DetermineFilterIndex() does not link */
  unsigned int DetermineFilterIndexNonInlined() const;

  /** Gives the appropriate filter based on DetermineFilterIndex() */
  const Filter& DetermineFilter() const;

  /** string representation for logging, includes neighborhood.
   * @param level 0 is elements, 1 is with weighting and distance */
  std::string ToString(int level = 0) const;

  /** for debugging. Sums the weights of all neighbors, ... */
  void Dump();

  /** The complete filter settings. There is a lot of copy and paste data to other elements.
   * But as we might be design and region dependent, might be robust and robust might also
   * store the neighborhood it is simply better to have the expensive instances.
   * Without regularization or with slopes/ perimeter this is empty. Normally we have one entry, with robust it is more  */
  StdVector<Filter> filter;

private:

  /** We need our base design element to do the filtering */
  DesignElement* de_;
};


/** required in DesignSpace and DesignMaterial */
class DesignID
{
public:
   DesignID(DesignElement::Type design = DesignElement::NO_TYPE, MultiMaterial* mm = NULL, double rb = -1.0, bool enforce_bounds = false)
   {
     this->design = design;
     this->multimaterial = mm;
     this->relative_bound = rb;
     this->enforce_bounds = enforce_bounds;
   }

   DesignElement::Type design;
   /** index. -1 for non-multimaterial */
   MultiMaterial* multimaterial;

   /** relative bounds for design, negative if not applicable, size by number of design types */
   double relative_bound;

   /** shall the bounds be enforced when loading external designs */
   bool enforce_bounds;
};

/**required in Design space**/
class TensorElement
{
public:
   TensorElement(DesignElement::Type design = DesignElement::NO_TYPE, MultiMaterial* mm = NULL)
   {
     this->design = design;
     this->multimaterial = mm;
   }

   DesignElement::Type design;
   /** index. -1 for non-multimaterial */
   MultiMaterial* multimaterial;
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

  /** debug output */
  std::string ToString();

  SolutionType solutionType;

  /** Finds the proper design element by element number. DEFAULT if not given */
  DesignElement::Type design;

  /** optionally filtered or plain */
  DesignElement::Access access;

  /** Which of the values? */
  DesignElement::ValueSpecifier value;

  /** An optional detail for values COST_GRADIENT and OBJECTIVE in PiezoSIMP case */
  DesignElement::Detail detail;

  /** for value = generic only */
  std::string generic;

  /** An optional excitation index. Negative for not set*/
  int excitation;
};

inline
double SIMPElement::CalcHeaviside(double input_value, const GlobalFilter* f)
{
  assert(f->GetType() == Filter::DENSITY);
  assert(f->density == Filter::SOLID_HEAVISIDE || f->density == Filter::VOID_HEAVISIDE);

  double result;

  double b = f->beta;
  assert(b >= 0.0 && b < 2000);

  if(f->density == Filter::SOLID_HEAVISIDE)
  {
    double tmp = 1.0 - std::exp(-b * input_value) + input_value * std::exp(-b);
    result = f->non_lin_scale * tmp + f->non_lin_offset;

    // no LOG_DBG possible due to inline
  }
  else // if(f.density_ == Filter::VOID_HEAVISIDE)
  {

    double first    = std::exp(-1.0 * b * (1.0 - input_value));
    double second   = -1.0 * (1.0 - input_value) * std::exp(-1.0 * b);

    assert(f->non_lin_scale > 1e-2);

    result = f->non_lin_scale * (first + second) + f->non_lin_offset;
    // std::cout << "CH: el=" << de_->elem->elemNum << " iv=" << input_value << " b=" << b << " lb=" << lb << " (ub-lb)=" << (ub-lb)
    //           << " +1st=" << first << " +2nd=" << second << " -> " << result << std::endl;
  }

  return result;
}

inline
double SIMPElement::CalcTanh(double input_value, const GlobalFilter* f)
{
  assert(f->GetType() == Filter::DENSITY);
  assert(f->density == Filter::TANH);

  double b = f->beta;
  double e = f->eta;

  assert(b >= 0.0 && b < 2000);

  // 1 - 1/(exp(2*beta*(x-param)) + 1)
  double func = 1.0 - 1.0/(std::exp(2.0 * b * (input_value - e)) + 1.0);
  double result = f->non_lin_scale * (func) + f->non_lin_offset;

  //std::cout << "CT: b=" << b << " e=" << e << " s=" << f->non_lin_scale << " o=" << f->non_lin_offset << " iv=" << input_value << " func=" << func << " -> " << result << std::endl;
  return result;
}

inline
double SIMPElement::CalcMaterial(double filter_value, double plain_value, const GlobalFilter* f)
{
  assert(f->GetType() == Filter::DENSITY);
  assert(f->density == Filter::MATERIAL);
  assert(plain_value != -1);

  double func = f->mat_phase.Transform(plain_value) * f->mat_scale.Transform(filter_value);
  double result = f->non_lin_scale * (func) + f->non_lin_offset;

  return result;
}


} // end of namespace

#endif /*DESIGNELEMENT_HH_*/
