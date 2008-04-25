#ifndef DESIGNELEMENT_HH_
#define DESIGNELEMENT_HH_

#include "General/Enum.hh"
#include "Domain/elem.hh"

namespace CoupledField
{
class Elem;
class ParamNode;
class Condition;
class ResultDescription;
class SIMPElement;

/** This is the base design element, storing all information related to 
 * elements. For ersatz material this is a 1:1 relation to finite elements.
 * Method specific information (SIMP, Level-Set, FreeMat, ...) is in add-ons. */
class DesignElement
{     
public:
  /** The empty constructor is only for the StdVector */
  DesignElement();

  /** This sets the DesignElement with the values from the XML file.
   * Is slow as it does the same evaluation often but is only O(n) */
  DesignElement(ParamNode* pn);

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
  typedef enum { DEFAULT = -2, NO_TYPE = -1, DENSITY = 0, POLARIZATION = 1} Type;

  /** types for GetFilteredValue() */
  typedef enum { DESIGN, DESIGN_COST_GRADIENT, COST_GRADIENT, CONSTRAINT_GRADIENT, 
    WEIGHT, OBJECTIVE, NUM_NEIGHBOURS} ValueSpecifier;

    /** This specifies result! details for value = COST_GRADIENT or OBJECTICE in the PiezoSIMP case
     * or it does the very special purpose symmetry for the ValueSpecifier */
    typedef enum { NONE, SYMMETRY } Detail;


    /** Allows to set the design element. Does it PLAIN */
    void SetDesign(double value) { this->design = value; } 

    /** Gets the design element
     * @param access if plain the rho value if SMART and filtering is enabled the filtered value */
    double GetDesign(Access access) const;

    /** Checks out specialResult[]! */
    double GetValue(ResultDescription* rd) const;

    /** This is a generic getter. <p>
     * <p>Note, that GetObjectiveGradient(SMART) gives different values than 
     * GetValue(COST_GRADIENT, SMART)!!! because of an necessary multiplication by
     * the densities of the element in the filter.!</p>
     * </p>Not that this works no to CONSTRAINT_GRADIENT as there is also and index required! */
    double GetValue(ValueSpecifier vs, Access access) const;   

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

    /** internal helper to get the value by type */
    double GetValue(ValueSpecifier valueSpecifier) const;
    
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
    
private:
  

  /** The scalar value. Public access only via getter to handle filtering. */
  double design;

  /** The gradient of the cost function w.r.t. this element of the design space */
  double cost_gradient;

  /** The lower bound of this design variable. Redundant but such faster than look it up */
  double lower_;

  double upper_;

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
  /** the barycenter of this element, set only when filter. In 2D the last
   * coordinate is zero -> it is too complicated to make a template for DesignElement */
  Point location_;

  /** Do filtering? */
  bool filter_;     

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
