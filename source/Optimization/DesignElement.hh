#ifndef DESIGNELEMENT_HH_
#define DESIGNELEMENT_HH_

#include "General/Enum.hh"
#include "Domain/elem.hh"
#include "Utils/tools.hh"
#include "MatVec/vector.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField
{
class Elem;
class ParamNode;
class ParamNode;
class DesignSpace;
class DesignElement;
class ResultDescription;
class SIMPElement;
class VicinityElement;
class LevelSetElement;
class Objective;
class Condition;

/** This DesignElement package provides information about the direct neighbours for uniform cartesian
 * quadrilateral or hexahedral meshes. */
class VicinityElement
{
public:
  /** resizes design and set entries to NULL */
  VicinityElement();

  /** Initializes all elements, ignores multiple design elements!
   * Adds VicinityElement packages to the design elements*/
  static void Init(DesignSpace* design);

  /** Gives the number of not NULL entries. */
  int GetNumberOfEntries() const;

  /** This are indices for the entries to design */
  enum Neighbour {X_P = 0, X_N = 1, Y_P = 2, Y_N = 3, Z_P = 4, Z_N = 5, NONE = -1};

  /** Gives the neighbour elements */
  DesignElement* GetNeighbour(Neighbour idx) { return design[idx]; }

  /** dump method for logging */
  std::string ToString();

  /** Contains the next neighbours (only +/- x,y(,z) and not diagonal.
   * Ordered as +x, -x, +y, -y (+z, -z). As the elements are DesignElements only within this region.
   * If there is no neighbour (e.g. z in 2D or on the boundary) the value is NULL */
  StdVector<DesignElement*> design;

private:

  /** Helper method that identifies the neighborhood of two elements given by their nodal coordinates.
   * TODO: make use of the element neighbors containing the number of common nodes! */
  static bool IdentifyNeighbor(Matrix<Double>& reference, Matrix<Double>& other, int& dimension, bool& positive);
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

  /** types for GetFilteredValue() */
  typedef enum { DESIGN, DESIGN_COST_GRADIENT, COST_GRADIENT, CONSTRAINT_GRADIENT, WEIGHT, OBJECTIVE, NUM_NEIGHBOURS,
    LEVEL_SET_VALUE, LEVEL_SET_STATE, TOPGRAD_VALUE, SHAPEGRAD_VALUE, SHAPEGRAD_NODE_VALUE,
    LEVEL_SET_GRAD_XP, LEVEL_SET_GRAD_XN, LEVEL_SET_GRAD_YP, LEVEL_SET_GRAD_YN, LEVEL_SET_GRAD_ZP, LEVEL_SET_GRAD_ZN } ValueSpecifier;

  BaseDesignElement();

  /** Allows to set the design element. */
  void SetDesign(double value) { this->design = value; }

  /** Return the design value */
  double GetDesign() const { return(this->design); }

  /** Get the gradient values for either objective or constraint.
   * @param f if given the the objective inclusive penalty. If NULL and g is NULL then the sum over the objectives
   * @param g if given the constraint
   * @param f, g or sum f */
  double GetGradient(const Objective* f, const Condition* g) const;

  /** Get the gradien values with summed objective gradients
   * @param objectives if given then the values of all gradiens with the penalties from objective are returned
   * @param g if given as the other GetGradient() */
  double GetGradient(const StdVector<Objective*> objectives, const Condition* g) const;

  /** Sum app the old value (get and set together) */
  void AddGradient(const Objective* f, const Condition* g, double value);

  /** Reset either gradients of the class
   * @param vs either COST_GRADIENT or CONSTRAINT_GRADIENT */
  void Reset(ValueSpecifier vs);

  /**  Gets the lower bound of the desing variable -
   * up to now this are defaults by type */
  double GetLowerBound() const { return lower_; }

  /** The upper bound of the desing variable for the optimizer */
  double GetUpperBound() const { return upper_; }

  /** Set the lower bound of the design variable */
  void SetLowerBound(const double v) { lower_ = v; }

  /** Set the upper bound of the design variable */
  void SetUpperBound(const double v) { upper_ = v; }

  /** adjusts length of the gradient vectors possibly not known during creation */
  void PostInit(int objectives, int constraints);

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
   * Is slow as it does the same evaluation often but is only O(n) */
  DesignElement(PtrParamNode pn, Elem* elem);

  ~DesignElement();


  /** This defines how to acces variables (design, objective_gradient, ...),
   *  PLAIN is the value and SMART does a filtering if enabled otherwise also as PLAIN */
  typedef enum { PLAIN, SMART } Access;

  /** Filter types we have
   * <ul>
   *   <li>RADIUS: this is the implementation following Sigmund in the 99lines paper.
   *               The drawback is the discretization dependency.</li>
   *   <li>VOLUME_RADIUS: The radius is *value* times square/cube edge length where the
   *               square/cube has the volume of the element</li>
   * </ul> */
  typedef enum { RADIUS, VOLUME_RADIUS, MAX_EDGE } Filter;

  /** The type of this design element, influences the Get*Bound() methods.
   * By definition the design elements are stored in the ordering of the type!! */
  typedef enum { UNITY = -5, NO_DERIVATIVE = -4, TENSOR_TRACE = -3, DEFAULT = -2, NO_TYPE = -1, DENSITY = 0, POLARIZATION = 1, EMODUL, POISSON, LAMELAMBDA, LAMEMU, EMODULISO, POISSONISO, GMODUL} Type;


    /** This specifies result details for various ValueSpecifier/Detail combinations:
     * OBJECTIVE/SYMMETRY (check!)
     * COST_GRADIENT/FINITE_DIFF_COST_GRADIENT
     * COST_GRADIENT/ERROR_COST_GRADIENT
     * For the PiezoSIMP case:
     * COST_GRADIENT/MECH_MECH ... as ErsatzMaterial::CalcU1KU2() Parameters. QUAD is calcMode
     * or it does the very special purpose symmetry for the ValueSpecifier */
    typedef enum { NONE, SYMMETRY, FINITE_DIFF_COST_GRADIENT, ERROR_COST_GRADIENT, MECH_MECH, ELEC_ELEC, ELEC_ELEC_QUAD, ELEC_MECH, MECH_ELEC } Detail;

    /** Gets the design element
     * @param access if plain the rho value if SMART and filtering is enabled the filtered value */
    double GetDesign(Access access) const;

    /** Checks out specialResult[]!
     * @param dofs 1 for scalar values */
    void GetValue(ResultDescription& rd, StdVector<double>& out, unsigned int dofs) const;

    /** This is a generic getter. <p>
     * <p>Note, that GetObjectiveGradient(SMART) gives different values than
     * GetValue(COST_GRADIENT, SMART)!!! because of an necessary multiplication by
     * the densities of the element in the filter.!</p>
     * </p>Not that this works no to CONSTRAINT_GRADIENT as there is also and index required! */
    double GetValue(ValueSpecifier vs, Access access) const;

    /** internal helper to get the value by type */
    double GetValue(ValueSpecifier valueSpecifier) const;

    /** Gets the gradient of the cost function w.r.t. the design element.
     * <p>Note, that in the SMART case, the filtering is done with a multiplication by the desing value!
     * In other words, speaking of GetValue():
     * <ul><li>GetObjectiveGradient(PLAIN) = GetValue(COST_GRADIENT, PLAIN)</li>
     *     <li>GetObjectiveGradient(SMART) = GetValue(DESIGN_COST_GRADIENT, SMART)</li></ul>
     * </p>
     * @param access if plain the rho value if SMART and filtering is enabled the filtered value */
    double GetObjectiveGradient(Access access) const;

    /** Initilize the Enum. Currently called by Optimization::CreateInstance() */
    void static SetEnums();

    Type GetType() const { return type_; }
    
    /** Write key values as attributes */
    void ToInfo(PtrParamNode in) const;

    std::string ToString() { return ToString(this); }

    /** makes a short dump, handles NULL */
    static std::string ToString(DesignElement* de);
    
    /** to make the class polymorphic and we can dynamic_cast<> it */
    /** Pointer to the element of the region, paramter for integration, ... */
    Elem*  elem;

    /** This stores special values during calculation for visualization via
     * result description in XML: See DesignSpace::GetSpecialResultIndex()
     * <result id="optResult_2" design="density" access="plain" value="costGradient" />
     * <result id="optResult_3" design="density" access="plain" value="objective" /> */
    double specialResult[6];

    static Enum<Filter> filter;

    static Enum<Type> type;

    static Enum<ValueSpecifier> valueSpecifier;

    static Enum<Access> access;

    static Enum<Detail> detail;

    /** This are our add-ons, not all are initialized! */
    SIMPElement* simp;

    /** The vicinity (structured grids, e.g. for Level-Set */
    VicinityElement* vicinity_;

    /** The level-set element, will be destroyed by LevelSet */
    LevelSetElement* lse_;
    
    /** The topgrad element, will be destroyed by TopGrad */
    TopGradElement *tge;

    /** calculates the location on request and stores it */
    Point* GetLocation();

private:
  /** Here we share the base initialization for the constructors */
  void Init();

  /** returns the non-scalar values */
  void GetValue(ValueSpecifier sp, StdVector<double>& out) const;  
  
  /** the barycenter of this element only set on request. */
  Point* location_;

  /** what is our design type */
  Type type_;

};


/** This is a add-on to DesignElement and contains all data specific for SIMP,
 * this is especiall the filtering stuff */
class SIMPElement
{
public:
  /** Base init of element */
  SIMPElement(DesignElement* base);

  /** does the core for GetDesign(). GetObjectiveGradient(), ... */
  double GetFilteredValue(DesignElement::ValueSpecifier valueSpecifier, bool design_weighted) const;

  /** Neigborhood is element and precalculated distance */
  struct NeighbourElement
  {
  public:
    /** read the variable */
    DesignElement* neighbour;

    /** precalculated weight by distance - to be multiplied with the density! */
    double        weight;

    /** the distance in domain dimensions! */
    double        distance;
  };

  /** The weight of THIS element to be summed to 1.0 with all neighbor weights */
  double weight;

  /** The neighbors if filter otherwise empty. Set by InitFilter().
   * The element itself is NOT part of the neighborhood! */
  StdVector<NeighbourElement> neighborhood;

  /** for debugging. Sums the weights of all neighbors, ... */
  void Dump();

  /** Do filtering? */
  bool filter;

private:

  /** We need our base design element to do the filtering */
  DesignElement* de_;
};


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
};


} // end of namespace

#endif /*DESIGNELEMENT_HH_*/
