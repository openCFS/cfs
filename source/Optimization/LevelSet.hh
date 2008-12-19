#ifndef LEVELSET_HH_
#define LEVELSET_HH_
#include "Optimization/BaseOptimizer.hh"
#include "Optimization/DesignElement.hh"

namespace CoupledField
{
  class Optimization;
  class ParamNode;
  class DesignSpace;

  /** A LevelSetElement has links to its nodes, where the actual information
   * is stored. This LevelSetNode elements the correspond to our classical design elements
   * of type density as we directly correlate the level set value with the element density.
   * There are ghost (design) elements stored in the LevelSetNode elements stored in the
   * DesignElement units.
   * This is clearly not optimized for (cache) speed, but for a good matching into the
   * (simp optimized) optimization framework */
  class LevelSetElement
  {
  public:

    /** for StdVector only */
    LevelSetElement() {};

    /** resises nodes */
    LevelSetElement(double subLevelSpacing);

    /** Gives the center value by averaging */
    double GetCenterValue();

    /** do we contain the front? */
    bool ContainsFront();

    /** Finds the closest distance from the node to the front within this element.
     * Does a subdivision of this element and works based on this discretization.
     * This element needs to contain the front! Check first! */
    double CalcDistance(LevelSetNode* lsn);

    /** Debug output for a single element */
    std::string ToString();

    /** Debug output for the complete space */
    static std::string ToString(StdVector<LevelSetElement>& space);

    /** in 2D 4 entries (starting with the upper right corner and then counter clock wise,
     * in 3D 8 entries. The LevelSetNode has a link to the DesignElement which has
     * the vicinity. Some nodes might be properties of ghost elements! */
    StdVector<LevelSetNode*> nodes;

    /** This is the node number of the FE design space that coincides with the LSE center */
    unsigned int nodeNumber;

  private:

    /** This is a helper data-structure for CalcDistance */
    struct DistanceHelper
    {
      Point  sub_location;
      double value;    // by bilinear interopolation
      double distance; // cache
    };

    /** This is a static helper reuse attribute for CalcDistance to save memory services */
    static StdVector<DistanceHelper> distance_tmp_;

    /** Our spacing for CalcDistance, we subdivide the level set element with this spacing
     * to finde a discrete discante value from the implicit front */
    double subLevelSpacing_;
  };


  /** The level set solver sets the vicinity element in the design elements and the level set package */
  class LevelSet : public BaseOptimizer
  {
  public:
    LevelSet(Optimization* optimization, ParamNode* pn);

    ~LevelSet() {}

    void SolveProblem();

    /** Search for level set element by the node number */
    LevelSetElement* Find(unsigned int nodeNumber);

    /** This defines a action, as the level set method is naturally a combination of multiple actions */
    class Action
    {
    public:
      /** empty constructor just for StdVector() */
      Action() {};

      /** read ourself from XML and does basic plausibility checks */
      Action(ParamNode* pn);

      /** the action types we have: */
      enum Type { SIGNED_DISTANCE_FIELD  /*!< A real, expensive signed distance field */,
                  FIRST_ORDER_FD         /*!< Do first oder finite differences evolution of the Hamilton-Jacobi equation */,
                  TOP_GRAD               /*!< Calcutate the topology gradient and create one hole */,
                  TRIVIAL_HOLE           /*!< Simply searches for the first possible one-element hole to create */ };

      static Enum<Type> type;

      /** our actual action */
      Type type_;

      /** this value gives, how often to trigger this action */
      int modulus;

      /** how often to perform this action within this action/iteration. Default = 1*/
      int perform;

    };

  private:

    /** based on the current level set value sets the front, updated the front bit in the level set elements */
    void FindFront();

    /** Calls for all level set nodes CalcGradients() and then UpdateNormal() and UpdateCurvature() */
    void UpdateLevelSetNodes();

    /** outsourced constructor stuff, adds the LevelSetNodes to the Optimization::design DesignSpace and
     * adds and sets ghost (design) elements to the design vicinity. Also sets the nodes array. */
    void SetupDesign();

    /** outsourced constructor stuff, based on the "extended" design greates a dual design space of level
     * set elements */
    void SetupSpace(double subLevelSpace);

    /** calculates a signed distance field, density is updated.
     * Is based on current boundary information in the level set elements, call FindFront() first! */
    void CalcSignedDistanceField();

    /** Calculates a single finite differences hamilton jacobi step */
    void EvolveFirstOrderFiniteDifferences();

    /** calculates from the larges abs(velocity=costGradient) and the min edge size (assuming uniform
     * design domain) by the CFL condition the appropriate delta_t.
     * @return this delta_t will not cross more than one elmenet for all velocities in the worst case */
    double FindMaxTimeStep();

    /** This traverses the design elements and makes a hole at the first element which has a
     * complete solid full neighbourhood */
    void MakeTrivialHole();

    /** appends another level set element to space, there as much level set elements as
     * nodes of the (finite element) design elements */
    void AddLevelSetElement(DesignElement* upper_right, double subLevelSpace);

    /** This queries if a specific action is to be executed. Evaluates the action list
     * @return NULL if not triggered or a pointer to the relevant action (to evaluate the perform attribute) */
    Action* IsTriggered(Action::Type action, int iteration);


    /** The level set element design space */
    StdVector<LevelSetElement> space;

    /** The original design spaces holds the actual LevelSetNode elements in its DesignElement elements.
     * Ghost elements are only accessible via this standard elements (ghost_ list or vicinity) */
    StdVector<DesignElement>* design;

    /** This is a list of all level set nodes, this is design + ghosts. Note, that via the de_ attribute
     * one gets the design elements as with design. Consider the data as unsorted! */
    StdVector<LevelSetNode*> nodes;

    /** This are the defined actions */
    StdVector<Action> action_;

    /** to speed up find */
    unsigned int last_node_index_;
  };


} // end of namespace

#endif /*LEVELSET_HH_*/
