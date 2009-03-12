#ifndef DESIGNELEMENT_HH_
#define DESIGNELEMENT_HH_

#include "General/Enum.hh"
#include "Domain/elem.hh"
#include "Utils/tools.hh"
#include "MatVec/vector.hh"

namespace CoupledField
{
class Elem;
class ParamNode;
class InfoNode;
class DesignSpace;
class DesignElement;
class Condition;
class ResultDescription;
class SIMPElement;
class VicinityElement;
class LevelSetNode;

/** This DesignElement package provides information about the direct neighbours for uniform cartesian
 * quadrilateral or hexahedral meshes. */
class VicinityElement
{
public:
  /** resises design and set entries to NULL */
  VicinityElement();

  /** Initializes all elements, ignores multiple design elements!
   * Adds VicinityElement packages to the design elements*/
  static void Init(DesignSpace* design);

  /** Gives the number of not NULL entries. This is size of design for enterior elements
   * as long as there are no Ghost elements set (LevelSet!!) */
  int GetNumberOfEntries();

  /** This are indices for the entries to design */
  enum Neighbour {X_P = 0, X_N = 1, Y_P = 2, Y_N = 3, Z_P = 4, Z_N = 5, NONE = -1};

  /** Gives the Neighbour elements */
  DesignElement* GetNeighbour(Neighbour idx) { return design[idx]; }

  /** Gives a chained neighbour, breaks and returns NULL on not existence */
  DesignElement* GetNeighbour(Neighbour first, Neighbour second);

  /** Sets the full neighbourhood, which includes also edge and corner elements. Is
   * as such larges than the sorted design neighbourhood
   * @return out from X_P clockwise. NULL at positions where there is no neighbour */
  void GetFullNeighbourhood(StdVector<DesignElement*>& out);

  /** dump method for logging */
  std::string ToString();

  /** Contains the next neighbours (only +/- x,y(,z) and not diagonal.
   * Ordered as +x, -x, +y, -y (+z, -z). As the elements are DesignElements only within this region.
   * If there is no neighbour (e.g. z in 3D or on the boundary the value is NULL */
  StdVector<DesignElement*> design;

private:

  /** Helper method that identifies the neighbourhood of two elements given by their nodal coordinates. */
  static bool IdentifyNeighbor(Matrix<Double>& reference, Matrix<Double>& other, int& dimension, bool& positive);

};


/** This is the base design element, storing all information related to
 * elements. For ersatz material this is a 1:1 relation to finite elements.
 * Method specific information (SIMP, Level-Set, FreeMat, ...) is in add-ons. */
class DesignElement
{
public:
  /** The empty constructor is the StdVector and for ghost elements */
  DesignElement();

  /** This sets the DesignElement with the values from the XML file.
   * Is slow as it does the same evaluation often but is only O(n) */
  DesignElement(ParamNode* pn, Elem* elem);

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
  typedef enum { NO_DERIVATIVE = -4, TENSOR_TRACE = -3, DEFAULT = -2, NO_TYPE = -1, DENSITY = 0, POLARIZATION = 1, EMODUL, POISSON, LAMELAMBDA, LAMEMU} Type;

  /** types for GetFilteredValue() */
  typedef enum { DESIGN, DESIGN_COST_GRADIENT, COST_GRADIENT, CONSTRAINT_GRADIENT, WEIGHT, OBJECTIVE, NUM_NEIGHBOURS,
    LEVEL_SET_VALUE, LEVEL_SET_NORMAL, LEVEL_SET_CURVATURE, LEVEL_SET_FIRST_GRAD, LEVEL_SET_SECOND_GRAD } ValueSpecifier;

    /** This specifies result! details for value = COST_GRADIENT or OBJECTICE in the PiezoSIMP case
     * or it does the very special purpose symmetry for the ValueSpecifier */
    typedef enum { NONE, SYMMETRY } Detail;


    /** Allows to set the design element. Does it PLAIN */
    void SetDesign(double value) { this->design = value; }

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

    /** Allows to set the gradient of the cost element. Does it PLAIN */
    void SetObjectiveGradient(double value) { this->cost_gradient = value; }

    /** Gets the gradient of the cost function w.r.t. the design element.
     * <p>Note, that in the SMART case, the filtering is done with a multiplication by the desing value!
     * In other words, speaking of GetValue():
     * <ul><li>GetObjectiveGradient(PLAIN) = GetValue(COST_GRADIENT, PLAIN)</li>
     *     <li>GetObjectiveGradient(SMART) = GetValue(DESIGN_COST_GRADIENT, SMART)</li></ul>
     * </p>
     * @param access if plain the rho value if SMART and filtering is enabled the filtered value */
    double GetObjectiveGradient(Access access) const;

    /** @param condition contains a unique index! */
    void SetConstraintGradient(const Condition* condition, double value);

    /** @see SetConstraintGradient() */
    double GetConstraintGradient(const Condition* condition) const;

    /** Initilize the Enum. Currently called by Optimization::CreateInstance() */
    void static SetEnums();

    Type GetType() const { return type_; }

    /**  Gets the lower bound of the desing variable -
     * up to now this are defaults by type */
    double GetLowerBound() const { return lower_; }

    /** The upper bound of the desing variable for the optimizer */
    double GetUpperBound() const { return upper_; }

    /** Set the lower bound of the design variable */
    void SetLowerBound(const double v) { lower_ = v; }

    /** Set the upper bound of the design variable */
    void SetUpperBound(const double v) { upper_ = v; }

    /** Writes key values as attributes */
    void ToInfo(InfoNode* in) const;

    std::string ToString() { return ToString(this); }

    /** makes a short dump, handles NULL */
    static std::string ToString(DesignElement* de);


    /** <p>The gradient of the constraint function w.r.t. this element of the design space.
     * Every constraint contains an unique index attribute (which is the order in the xml file)
     * only for the purpose to index this vector.</p>
     * <p>Therefor this vector has to be initalized on runtime</p> */
    StdVector<double> constraintGradient;

    /** to make the class polymorphic and we can dynamic_cast<> it */
    /** Pointer to the element of the region, paramter for integration, ... */
    Elem*  elem;

    /** This stores special values during calculation for visualization via
     * result description in XML:<pre>. See DesignSpace::GetSpecialResultIndex()
     * <result id="optResult_2" design="density" access="plain" value="costGradient" />
     * <result id="optResult_3" design="density" access="plain" value="objective" /></pre> */
    double specialResult[3];

    static Enum<Filter> filter;

    static Enum<Type> type;

    static Enum<ValueSpecifier> valueSpecifier;

    static Enum<Access> access;

    static Enum<Detail> detail;

    /** This are our add-ons, not all are initialized! */
    SIMPElement* simp;

    /** The vicinity (structured grids, e.g. for Level-Set */
    VicinityElement* vicinity;

    /** The level-set element */
    LevelSetNode* lsn;

    /** calculates the location on request and stores it */
    Point* GetLocation();

    /** Set the location indirect by a reference element. Assumes uniform elements!
     * Used for the level set method for the location of ghost elements
     * @param ref the position of the reference element
     * @param pos1 if X_P, this->location = ref + length(ref, x). if X_N it is ref - length(ref, x)
     * @param pos2 if not NONE allows a diagonal position*/
    void SetLocation(DesignElement* ref, VicinityElement::Neighbour pos1, VicinityElement::Neighbour pos2 = VicinityElement::NONE);

private:
  /** Here we share the base initialization for the constructors */
  void Init();

  /** returns the non-scalar values */
  void GetValue(ValueSpecifier sp, StdVector<double>& out) const;

  /** The scalar value. Public access only via getter to handle filtering. */
  double design;

  /** The gradient of the cost function w.r.t. this element of the design space */
  double cost_gradient;

  /** The lower bound of this design variable. Redundant but such faster than look it up */
  double lower_;

  double upper_;

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

  /** Initialized filtering. Sets the locations and neighbours and enables the flag
   * @param data the region, we substitue a seperate container class by this static method
   * @param pn our parameter section */
  static void InitFilter(StdVector<DesignElement>& data, ParamNode* pn);

  /** does the core for GetDesign(). GetObjectiveGradient(), ... */
  double GetFilteredValue(DesignElement::ValueSpecifier valueSpecifier, bool design_weighted) const;

  /** Neigbourhood is element and precalculated distance */
  struct NeighbourElement
  {
  public:
    /** read the variable */
    DesignElement* neighbour;

    /** precalculated weight by distance - to be multiplid with the density! */
    double          weight;

    /** the distance in domain dimensions! */
    double          distance;
  };

  /** The weight of THIS element to be summed to 1.0 with all neighbour weights */
  double weight;

  /** The neighbours if filter otherwise emtpy. Set by InitFilter().
   * The element itself is NOT part of the neighbourhood! */
  StdVector<NeighbourElement> neighbourhood;

  /** for debugging. Summs the weights of all neighbours, ... */
  void Dump();

private:

  /** Do filtering? */
  bool filter_;

  /** We need our base design element to do the filtering */
  DesignElement* de_;
};



/** The LevelSetNodes are the part of LevelSet::LevelSetElements that correspond to design elements.
 * The level set consists of elements, but the properties are actually stored on nodes. This
 * makes it more convenient to do a bilinear interpolation for the continious level set value.
 * The level set nodes coincide with design elements. Ghost (design) elements/ (level set) nodes
 * are stored in the design vicinity for sides and diagonal in this element in cornder/edge */
class LevelSetNode
{
public:

  /* create all level set nodes before calling stuff like UpdateCurvature(), ...
   * @param value the initial level set value -> see LevelSet.cc*/
  LevelSetNode(DesignElement* de, double value);

  /** We have to delete the ghost elements manually! */
  ~LevelSetNode();

  /** updates the design (density) value by an interpretation of the current level set value */
  void UpdateDesign();

  /** creates ghost elements and add them to vicinity where we need such elements.
   * Identified by missing entries to vicinity. The ghost elements are deleted by
   * the LeveSetLevelSetNode destructor. See code for convention to handle "edge ghost elements".
   * The ghost elements have only "this" elemens in the vicinity set.
   * @param value the inital level set value
   * @return the number of created ghost elements (3 in the corner)
   * @see SetGhostVicinity() */
  int CreateGhostElements(double value);

  /** This sets the vicinity for ghost elements, the problem is to integrate the corner elements! */
  void SetGhostVicinity();

  /** This is the base for many further methods, like the calculation of curvatuere and normal.
   * updates first_order_grad and second_order_grad */
  void CalcGradients();

  /** calculates the mean curvature, assumes design to be up to date. Assumes CalcGradients() to be called already */
  void UpdateCurvature();

  /** calculate the normalized (length = 1) outer normal vector n = grad phi/|| grad phi|| where phi is value.
   * Assumes CalcGradients() to be called already */
  void UpdateNormal();

  /** Calculates a first order space convex finite difference evolution of the Hamilton Jacobi equation.
   * This takes the cost_gradient from the underlying design element as velocity
   * @param delta_t shall confitm the CFL condition such that we don't move futher than one element */
  void CalcNextFirstOrderHamiltonJacobiStep(double delta_t);

  /** This are indices for the entries to design. This is poor! It's a clone from VicinityElement :( */
  enum Neighbour {X_P = 0, X_N = 1, Y_P = 2, Y_N = 3, Z_P = 4, Z_N = 5};

  /** Assumes, the vicinity of the underlying design element to be properly set! */
  DesignElement* GetNeighbour(Neighbour idx);

  /** Gives a chained neighbour, breaks and returns NULL on not existence. Handles the corner!! */
  DesignElement* GetNeighbour(Neighbour first, Neighbour second);

  /** set the level set value, mirrors to the ghost nodes */
  void SetValue(double value, bool mirror_to_ghosts = true);

  /** the level set value, usually denoted phi */
  double GetValue() { return value_ ; }

  /** Returns topological gradient */
  double GetTopGradient() { return top_gradient; }

  /** Allows to set the topological gradient */
  void SetTopGradient(double const value) { this->top_gradient = value; }

  /** Are we a ghost element? This is an additional layer of element around the border, they are
   * not mapped to "traditional" design elements and need to be constructed and deleted therefore
   * explicitly */
  bool IsGhost() const { return de_->elem == NULL; }

  /** Check the given element for ghost, if NULL then false */
  static bool IsGhost(DesignElement* de);

  /** Is this node representing solid material */
  bool IsSolid() const { return value_ < 0; }

  /** Sets this value to represent a hole.
   * Does *NOT* update the underlying design value via UpdateDesign() but set ghosts! */
  void SetHole() { SetValue(1.0); }

  /** This is a shortcut to ask the vicinity if they are ghosts */
  StdVector<DesignElement*>& GetGhosts() { return ghosts_; }

  /** The mean curvature value */
  double curvature;

  /** In this vector element we store the normalized outer surface normal. */
  Vector<double> normal;

  /** Here we store the gradient, which is the first spacial derivative.
   * To be evaluated by CalcGradients() */
  Vector<double> first_order_grad;

  /** Here we store the second order spacial derivative.
   * To be evaluated by CalcGradients() */
  Vector<double> second_order_grad;

  /** if not within the domain we don't have all neighbours */
  bool domain_boundary;

  DesignElement* de_;

private:
  /** trivial inline function, for power 2 of doubles */
  static double square(const double in) { return in * in; }

  /** the leve set value is private because we mirror the values to the ghosts! */
  double value_;

  /** The value of the current topological gradient for this element */
  double top_gradient;

  /** Here we store the real ghost elements we create, this includes corners and edges (3d).
   * These elements are linked via the vicinity of the underlying design element. This reference here is just
   * for easy and clean deletion */
  StdVector<DesignElement*> ghosts_;
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
  ResultDescription(ParamNode* pn);

  /** killme -> use the typedef version when the Enums eliminated environment.hh */
  int solutionType;

  /** Finds the proper design element by element number */
  DesignElement::Type design;

  /** optionally filtered or plain */
  DesignElement::Access access;

  /** Which of the values? */
  DesignElement::ValueSpecifier value;

  /** An optional detail for values COST_GRADIENT and OBJECTIVE in PiezoSIMP case */
  DesignElement::Detail detail;

  /** As long as solutionType is an int this gives the string representation */
  std::string ToString();
};


} // end of namespace

#endif /*DESIGNELEMENT_HH_*/
