#ifndef LEVELSET_HH_
#define LEVELSET_HH_
#include <assert.h>
#include <stddef.h>
#include <algorithm> // for for_each
#include <functional> // for function objects
#include <map>
#include <string>
#include <utility> // for pair
#include <vector>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/Enum.hh"
#include "MatVec/vector.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Utils/StdVector.hh"

namespace CoupledField
{
  class LevelSetNode;
  class Optimization;
  class ShapeGrad;
  class linElastInt;
  
  /** typedef for convenience */
  typedef std::pair<const LevelSetNode *, const LevelSetNode *> lsnodepair;
  
  using std::vector;
  
  /** When using shape optimization with a levelset we have to calculate the
   *  implicitly given isolines/surfaces where the levelset is 0. 
   *  On these objects, we then have to perform an integration to get the shape gradient. */
  struct IntersectionObject
  {
    /** contains all the points the object is made of (2 in 2D, 3 - 6(?) in 3D) */
    vector<Vector<double> > points;
    vector<Vector<double> > displacements; // displacements in the points
    vector<unsigned int> indices; // indices of the points
    double integralValue;
    unsigned int elemNum; // maybe do not need this
  };

  /** helper function to determine if two point coordinates match */
  inline bool IsClose(const double d1, const double d2);
  
  /** takes two pointers to levelset nodes and tests if the levelset-0-line
   *  intersects the element */
  bool HasIntersection(const LevelSetNode *node1, const LevelSetNode *node2);
  bool HasIntersection(const lsnodepair &p);
  
  /** The LevelSetNodes are the part of LevelSet::LevelSetElements that correspond to design elements.
   * The level set consists of elements, but the properties are actually stored on nodes. This
   * makes it more convenient to do a bilinear interpolation for the continuous level set value.
   * The level set nodes coincide with grid nodes. */
  class LevelSetNode
  {
  public:

    /** create all level set nodes
      * @param value the initial level set value -> see LevelSet.cc*/
    explicit LevelSetNode(const double value = -1.0, const unsigned int cfs_number = 0);
    ~LevelSetNode() {}
    
    /** state of the node for fast marching algorithm */
    enum State {FAR = -3, CLOSE = -2, ACCEPTED = -1, NONE = 0, ACCEPTED_POS = 1, CLOSE_POS = 2, FAR_POS = 3};

    /** Assumes, the vicinity of the underlying design element to be properly set!
     * sorting is according to VicinityElement::Neighbour: 0 = X_P, 1 = X_N, 2 = Y_P, 3 = Y_N, 4 = Z_P, 5 = Z_N */
    LevelSetNode* GetNeighbour(const VicinityElement::Neighbour idx) const { return neighbours_[idx]; }

    /** check all neighbours to see if it is accepted. this is used for the fast marching algorithm. */
    bool HasAcceptedNeighbour(const State state_accepted) const;

    /** get the cfs node number */
    unsigned int GetNumber() const { return global_node_number ; }

    /** iterates over all neighbours and checks for a sign change; thus, a boundary nodes
        is a node with that property. we have to be careful with this when determining the
        boundary after the initialization step where several nodes have 0.0 as levelset value */
    bool IsBoundary() const;
    
    /** current state of the node w. r. t. fast marching */
    State state;
    
    /** data of the levelset node */
    double value; // levelset value
    double phi_temp;
    double intersection_length;
    double f_ext; // extension velocity
    double shapegrad; // value of shape gradient

    /** a link to all the neighbours, sorted as in VicinityElement::Neighbour: 0 = X_P, 1 = X_N, 2 = Y_P, 3 = Y_N, 4 = Z_P, 5 = Z_N */
    vector<LevelSetNode*> neighbours_;
    
  private:
    
    unsigned int global_node_number;
  };
  
  /** Helper struct for SetNodeStates */
  struct SetStateToNone : public std::unary_function<LevelSetNode, void>
  {
    void operator()(LevelSetNode &node) 
    { 
      node.state = LevelSetNode::NONE;
    }
  };


  /** A LevelSetElement has links to its nodes, where the actual information
   * is stored. It can be asked for intersections with the zero levelset isosurface
   * and can update the simp density using an interpolated value from the attached
   * levelset nodes.
   * This is clearly not optimized for (cache) speed, but for a good matching into the
   * (simp optimized) optimization framework */
  class LevelSetElement
  {
  public:

    /** for StdVector only */
    explicit LevelSetElement(DesignElement* de = NULL, const double val = 0.0);
    ~LevelSetElement() {}

    /** updates the design (SIMP density) value by an interpretation of the current level set value 
      * so that the optimizer weighs elements less if they
      * are only partially covered by the computational domain */
    void UpdateDesign();

    /** do we contain the front (= the levelset-0-line/surface)? */
    bool ContainsFront() const;
    
    /** returns the edge nodes corresponding to the given index, where the edges are
     *  numbered as follows:
     *  
     *        ____6_____
     *       /|       /|
     *    11/ |7   10/ |5 
     *     /__|_2___/  |
     *    |   |_____|__|
     *   3|  /   4  |  /
     *    | /8     1| /9
     *    |/________|/        
     *         0    
     *  
     *  for 2D compatibility, the front face edges get numbers 0 - 3, after that the
     *  back face edges 4 - 7, last the remaining edges 8 - 11;
     *  @return a pair of two LevelSetNode pointers */
    const lsnodepair GetEdgeNodes(const unsigned int idx) const;

    /** Calculates the shape grad on the element and writes it to shapeGrad_;
     *  if ContainsFront is false, we unset shapeGrad_ */
    void CalcShapeGrad(const Vector<double>& elem_sol, linElastInt* bdb_form);
    
    /** pointer to the underlying designelement */
    DesignElement *de_;

    /** in 2D 4 entries (starting with the lower left corner and then counter clock wise,
      * in 3D 8 entries. The LevelSetNode has a link to the DesignElement which has the vicinity. */
    vector<LevelSetNode*> nodes_;
    
    /** calculate one shape gradient value for the whole element by averaging the
     *  calculations on the intersection objects */
    double shapeGradValue;

  private:

    /** calculates the intersection objects of the element. Check for intersection
     *  of levelset-0-isosurface before calling this function! After computation is done,
     *  the intersection objects are saved in the member vector */
    void CalcIntersectionObjects(const Vector<double> &elem_sol);
    
    /** Given two levelset nodes, this function calculates the barycenter point according
     *  to the levelset values in those nodes. It does not check for sign changes, so it
     *  should only be called if it is already clear that it is needed. The resulting point
     *  is put into the out Vector.
     *  @return the barycentrical weight, for later use */
    double CalcBarycenterNode(Vector<double> &outpoint, const lsnodepair &p) const;
    
    void CalcBarycenterDisplacement(Vector<double>& u_out, const double bary_weight,
                                    const Vector<double>& elem_sol, const int idx) const;
    
    double IntegrateIntersectionObject(IntersectionObject &o, linElastInt *bdb_form);
    
    vector<IntersectionObject> objects;
  };
  
  /** Debug output for a single element */
  const std::string ToString(const LevelSetElement &elem);

  /** helper class for sorting a vector of LevelSetNode-pointers by their levelset value */
  struct SortLevelSetNodesUsingLevelSetValue :
      public std::binary_function<const LevelSetNode*, const LevelSetNode*, bool>
  {
    public:
      bool operator()(const LevelSetNode *lhs, const LevelSetNode *rhs) const
      {
          return (lhs->value < rhs->value);
      }
  };

  /** helper class for sorting a vector of LevelSetNode-pointers by their node number */
  struct SortLevelSetNodesUsingNodeNumber :
      public std::binary_function<const LevelSetNode*, const LevelSetNode*, bool>
  {
    public:
      bool operator()(const LevelSetNode *lhs, const LevelSetNode *rhs) const
      {
        return (lhs->GetNumber() < rhs->GetNumber());
      }
  };

  /** This class holds the levelset needed for shape optimization
      in here, we assume the use of 2D rectangular or 3D cube elements at several places */
  class LevelSet
  {
  public:
    explicit LevelSet(Optimization* opt = 0, PtrParamNode pn = PtrParamNode());
    ~LevelSet() {}
    
    /** called from shapeoptimizer to solve the problem */
    void SolveProblem(const unsigned int iter);
    
    /** uses fast marching to build up the signed distance levelset function
     *  @param really_dump flag for debug output */
    void BuildSignedDistanceFunction(const bool really_dump = false);
    
    /** create a hole at node trivial_hole_node_nr_ */
    void MakeTrivialHole(const unsigned int node_nr);
    
    /** This function is called from DesignSpace::GetNodalValue(); */
    LevelSetNode* GetNodePointer(const unsigned int nodeNumber) const;
    
    /** calculates the gradient at a levelset node
     * @param node_nr the node number
     * @param dir the direction of the gradient:
     * 0 = X_P, 1 = X_N, 2 = Y_P, 3 = Y_N, 4 = Z_P, 5 = Z_N */
    double GetGradientAtNode(const unsigned int node_nr, const unsigned int dir) const;
    
    /** calculates the gradient at a levelset node
     *  @param node_nr the node number */
    double GetNormGradientAtNode(const unsigned int node_nr) const;
    
    /** Shortcut to update the SIMP-densities for all LevelSetElements; e. g. called by TopGrad */
    void UpdateDesignForAllLevelSetElements()
    {
      std::for_each(space_.begin(), space_.end(), std::mem_fun_ref(&LevelSetElement::UpdateDesign));
    }
    
    /** This defines an action, as the level set method is naturally a combination of multiple actions */
    class Action
    {
    public:
      /** read ourself from XML and does basic plausibility checks */
      explicit Action(PtrParamNode const pn);
      ~Action() {}

      /** the action types we have: */
      enum Type { SIGNED_DISTANCE_FIELD,  /*!< Calculate signed distance field */
                  TRIVIAL_HOLE,           /*!< Takes nodes from the xml node list and creates a hole */
                  DO_SHAPE_STEP           /*!< Perform a shape step */
                };

      static Enum<Type> type;

      int GetModulus() const { return modulus_; }
      int GetPerform() const { return perform_; }
      Type GetType()   const { return type_; }
      bool GetFirst()  const { return first_; }

    private:
      /** standard constructor is not to be used! */
      explicit Action() {}

      /** this value gives, how often to trigger this action */
      int modulus_;

      /** how often to perform this action within this action/iteration. Default = 1*/
      int perform_;
      
      /** perform the action in the first iteration? */
      bool first_;

      /** our actual action */
      Type type_;
    };

  private:
		/** do not copy the levelset! */
		LevelSet(const LevelSet &other);
		/** do not assign the levelset! */
		LevelSet& operator=(const LevelSet &other);

    /** outsourced constructor stuff; adds the LevelSetNodes to the Optimization::design DesignSpace. 
        In this function we attach a levelset element to each designelement. Also, for the construction of
        the new levelset element we compute the correct ordering of the levelset nodes so that we can
        retrieve the neighbourhood of the levelset nodes afterwards. We sort the levelset nodes counterclockwise
        starting with the lower left corner (in 3D the lower left corner of the front face = where the z-coordinate
        is smaller) 
        The main computation of the correct order of the points is done in 
        AddOrderedLevelsetNodesToLevelsetElement() below */
    void SetupSpace();

    /** this function executes the fast marching algorithm for the negative and the positive levelset values
     *  @param negative mode of fast marching, once for the negative, once the positive values
     *  @param really_dump controls dumping of the fast marching steps to the output file */
    void ExecuteFastMarching(const bool negative, const bool really_dump);

    /** helper function for signed distance calculation; walks over all levelset nodes and 
        pushes them into one of the sets below for the fast marching algorithm 
        @return signals if at least one boundary node could be found */
    bool SetSignedDistanceNodeStates();
    
    /** in the first iteration of the fast marching we need the values for the points in close. 
        this is done here using the information in the accepted sets */
    void FindNewNodeValueFromNeighbours(LevelSetNode *nodeptr);

    /** Here we figure out the correct order of the LevelsetNodes (counterclockwise, starting with lower left corner)
        and attach them to a LevelsetElement. We use the point coordinates to do this, which is not very nice.
        this function will be called for every levelset element */
    void AddOrderedLevelsetNodesToLevelsetElement(LevelSetElement &lse);

    /** function to build the neighbourhood information for the levelset nodes */
    void BuildNeighbourhoodOfLevelsetNodes(const DesignElement* de);
    
    /** helper function for BuildNeighbourhoodOfLevelsetNodes, handles upper left nodes (which have
     *  index 0 and 4 respectively */
    void SetNeighboursOfLowLeftNode(const int startIdx, const DesignElement *de);
    
    /** helper function for BuildNeighbourhoodOfLevelsetNodes, handles upper left nodes (which have
         *  index 1 and 5 respectively */
    void SetNeighboursOfLowRightNode(const int startIdx, const DesignElement *de);
    
    /** helper function for BuildNeighbourhoodOfLevelsetNodes, handles upper left nodes (which have
     *  index 3 and 7 respectively */
    void SetNeighboursOfUpLeftNode(const int startIdx, const DesignElement *de);
    
    /** helper function for BuildNeighbourhoodOfLevelsetNodes, handles upper right nodes (which have
     *  index 2 and 6 respectively */
    void SetNeighboursOfUpRightNode(const int startIdx, const DesignElement *de);

    /** walks over all levelset elements and calculates the shape gradient, writes it back 
      * to lse->shapeGrad_ */
    void CalcShapeGradientOnAllElements();
    
    /** walk over all elements and transport the levelset using the shape gradient value
     *  @param dt the time step length for this transport */
    void TransportLevelSet(const double dt);
    
    /** helper function for loops; returns for a given index the corresponding axis index */
    unsigned int GetAxisFromIndex(const unsigned int in) const
    {
      // we need to get 0 for 0,1, 1 for 2,3 and 2 for 4,5, so we use integer division 
      assert(in < 6);
      return in/2;
    }

    /** This queries if a specific action is to be executed. Evaluates the action list
     * @return NULL if not triggered or a pointer to the relevant action (to evaluate the perform attribute) */
    Action* IsTriggered(const Action::Type action, const int iteration);

    /** The level set element design space */
    vector<LevelSetElement> space_;

    /** This is a list of all level set nodes. Note, that via the de_ attribute
     * one gets the design elements as with design. Consider the data as unsorted! */
    vector<LevelSetNode> nodes_;

    /** This are the defined actions */
    vector<Action> action_;

    /** holds the edge lengths of the grid, [0] = dx, [1] = dy (, [2] = dz) */
    StdVector<double> edge_length_;

    /** map of cfs node number and levelset node pointer, has to be done because of regions */
    std::map<unsigned int, LevelSetNode*> lsnode_map_;

    /** to speed up find */
    unsigned int last_node_index_;
    
    /** read from xml how long the time step should be */
    double time_step_;

    /** read from xml if we dump every step of the fast marching to the output file */
    bool dump_fast_marching_;
    
    /** read from xml a node number for which we can dump all neighbours
     *  after construction of the vicinity elements in the ctor */
    unsigned int dump_lselement_;

    /** contains node numbers for the creation of trivial holes */
    vector<unsigned int> trivial_holes_;

    /** The original design spaces holds the actual LevelSetNode elements in its DesignElement elements */
    StdVector<DesignElement> *design_;

    /** \var ShapeGrad* optimization
     * \brief Reference to our optimization problem */
    ShapeGrad *optimization;
  };

  /** Debug output for the complete space */
  const std::string ToString(const vector<LevelSetElement> &space);
  
} // end of namespace

#endif /*LEVELSET_HH_*/
