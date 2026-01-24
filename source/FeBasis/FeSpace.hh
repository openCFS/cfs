#ifndef FILE_CFS_FESPACE_HH
#define FILE_CFS_FESPACE_HH

#include "BaseFE.hh"
#include "FeFunctions.hh"

#include "General/Environment.hh"
#include "Domain/ElemMapping/EntityLists.hh"
#include "Domain/ElemMapping/Elem.hh"
#include "Domain/Domain.hh"
#include "Domain/Results/ResultInfo.hh"

#include "General/Exception.hh"
#include "MatVec/Matrix.hh"
#include "Forms/IntScheme.hh"
#include "OLAS/algsys/AlgebraicSys.hh"
#include <boost/unordered_map.hpp>

namespace CoupledField {


// forward class declarations
class BaseFe;
class BaseFeFunction;
class SolStrategy;

// ========================================================================
 //  TYPE DEFINTION FOR POLYNOMIAL ORDER
 // ========================================================================
 //! This class describes the approximation order of an FeFuncion
 
 //! This class encapsulates the polynomial order, used to approximate
 //! a given function. It can be isotropic (like for Lagrangian Space)
 //! or anisotropic (e.g. hierarchical approximation).
class ApproxOrder {
public:

  //! Define mapping type 
  typedef enum { UNDEF, GRID_TYPE, SERENDIPITY_TYPE, TENSOR_TYPE  } PolyMapType;
  
  //! Constructor
  ApproxOrder();

  //! Constructor
  //! \param dim Spatial dimension
  ApproxOrder(UInt dim);

  //! Set isotropic order
  void SetIsoOrder( UInt isoOrder );

  //! Set anisotropic order
  void SetAnisoOrder( const Matrix<UInt>& anisoOrder );

  //! Specify the type of polynomial mapping
  void SetPolyMapping(PolyMapType type);

  //! Return, if approximation order is isotropic
  bool IsIsotropic() const;

  //! Return the polynomial mapping type
  PolyMapType GetPolyMapType() const;

  //! Return isotropic polynomial degree
  UInt GetIsoOrder() const;

  //! Return anisotropic polynomial degree
  void GetAnisoOrder( Matrix<UInt>& order ) const;

  //! Return maximum overall polynomial degree
  UInt GetMaxOrder() const;

  //! Return maximum polynomial degree per local direction
  const StdVector<UInt>& GetMaxOrderLocDir( ) const;

  //! Return polynomial order for a given dof (spatial ansisotropic)
  StdVector<UInt> GetDofOrder( UInt dof ) const;

  //! Print to string
  std::string ToString() const;


private:

  //! Spatial dimension
  UInt dim_;

  //! Isotropic order
  UInt isoOrder_;

  //! Anisotropic order
  Matrix<UInt> anisoOrder_;

  //! Maximum order w.r.t. to local coordinate direction
  StdVector<UInt> maxOrderLocDir_;

  //! Maximum order (independent of isotropy)
  UInt maxOrder_;

  //! Flag, if order is isotropic
  bool isIsotropic_;

  //! Type of polynonmial mapping
  PolyMapType polyMapType_;
};

  


 // ========================================================================
 //  Finite Element Space Definition
 // ========================================================================
 
//!  Base class for the Finite Element Space (FeSpace) 
/*!
  The Finite Element Space (FeSpace) represents the mathematical discrete 
  function space in which a unknown function is defined.
  
  It has the following functionality:
  - Collect pointers to reference elements of related functional space
  - Perform equation numbering
  
  Each unknown function is associated with an appropriate FeSpace. 
  The spaces differ in the following senses
  
    - Continuity of derivatives (H1, HCurl, HDiv, L2, Constants)
    - Order of shape functions
      -* Lower Order: The approximation order is explicitly encoded in the 
                      element itself (i.e. Lagrange Elements)
      -* Higher Order: The functions of the elements are composed of 1D 
                       polynomials (e.g. Legendre, Jacobi, Gegenbauer, Dubiner).
  
  The FeSpaceConst plays a special role: It can be used in conjunction with
  discrete unknowns, e.g. additional network equations.
*/

class FeSpace {
public:

  //! Enumeration type for element space
  //@{
  typedef enum { UNDEF_SPACE, CONSTANT, H1, HCURL, HDIV, L2 } SpaceType;
  static Enum<SpaceType> SpaceTypeEnum;
  //@}

  //! Enumeration type for polynomial shape functions
  //@{
  typedef enum { UNDEF_POLY, LAGRANGE, LEGENDRE, JACOBI } PolyType;
  static Enum<PolyType> PolyTypeEnum;
  //@}
  
  //! Enumeration type for boundary conditions
  //@{
  typedef enum{NOBC, HDBC = -2, IDBC = -3, CS = -4} BcType;
  static Enum<BcType> BcTypeEnum;
  //@}

  //! Flag for Grid determined mapping of element DOFs or polynomial based mapping 
  typedef enum {GRID,POLYNOMIAL} MappingType; 

  typedef enum {INTEG_MODE_ABSOLUTE,INTEG_MODE_RELATIVE} IntegOrderMode;
  static Enum<IntegOrderMode> IntegOrderModeEnum;

  struct IntegDefinition{
      IntScheme::IntegMethod method;
      IntegOrder order;
      IntegOrderMode mode;
  };

 
  //! struct for the single equationMap i.e. Nodal Eqns
  //! for the BC map we code the following
  //! <ul> 
  //!   <li>No BC</li>
  //!   <li>HomDirichlet</li>
  //!   <li>InHomDirichlet</li>
  //!   <li>Constraint</li>
  //! </ul>
  struct SingleEqnMap{

    //! Map for every node (key to map) its equations (values)
    std::unordered_map< Integer, StdVector<Integer> > eqns;
    
    //! Map for every node (key to map) its boundary conditions types
    std::unordered_map< Integer,StdVector<BcType> > BcKeys;
    
    
    //! Map for storing constraint information
    
    //! Special treatment of constraints:  we store for every slave 
    //! (node, dof)-pair the master (node, dof)-pair
    std::map< std::pair<Integer,Integer>, 
              std::pair<Integer,Integer>  > constraintNodes;
    
    //! Convenience operator for accessing the equations of a (virtual) node
    StdVector<Integer>& operator[](Integer key){
      return eqns[key];
    }
  };
  
  // ========================================================================
  
  //! Constructor
  FeSpace(PtrParamNode paramNode, PtrParamNode infoNode, Grid* ptGrid );

  //! Destructor
  virtual ~FeSpace();
  
  //! Read parameters from ParamNode and initialize
  
  //! Reads the parameters from a ParamNode and initializes itself,
  //! incorporating the information from a SolStrategy object.
  //! \param solStrategy Solution strategy object, which encapsulates
  //!                    additional information (splitting, number of
  //!                    solution steps  etc.)
  virtual void Init( shared_ptr<SolStrategy> solStrat ) {};

  //! Return type of FeSpace
  SpaceType GetSpaceType() { return type_;}
  
  //! Generate instance of specific element space type
  static shared_ptr<FeSpace> CreateInstance(PtrParamNode aNode, 
                                            PtrParamNode infoNode,
                                            SpaceType reqType,
                                            Grid* ptGrid,
                                            bool isAVExc = false );
  
  // ========================================================================
  //  INITIALIZATION 
  // ========================================================================
  //@{ \name Initialization

  //! Sets an offset to element orders. Right now this is only needed for
  //! Mixed TaylorHood approximation
  void SetElementOrderOffset(UInt offset);

  //! Return, if the space is hierarchical
  bool IsHierarchical() { return isHierarchical_;};
  
  //@}
  

  // ========================================================================
  // Multiple Region handling
  // ========================================================================
  //@{ \name Region Handling

  //! Get the Integration method and its order for the current region
  void GetIntegration(BaseFE * fe, RegionIdType region, 
                      IntScheme::IntegMethod & method,
                      IntegOrder & order);

  //! Return integration scheme
  shared_ptr<IntScheme> GetIntScheme() { return intScheme_;}
  
  //! Set the approximation for the given region according to the passed poly and integ Ids
  //! ToDO if integration and polynomial nodes are merged, of course only an orderId has to be passed
  void SetRegionApproximation(RegionIdType region, std::string polyId, std::string integId);

  //! Sets a default approximation for all regions
  void SetDefaultRegionApproximation();

  virtual void SetUseGradients(RegionIdType region){
    EXCEPTION("Not implemented here in baseclass!");
  }

  //@}
  // ========================================================================
  //  ELEMENT HANDLING
  // ========================================================================
  //@{ \name Element Handling

  //! Return pointer to reference element
  virtual BaseFE* GetFe( const EntityIterator ent ) = 0;
  
  //! Return pointer to reference element and adjusted integrationScheme
  
  //! This method returns the pointer to the reference according to the
  //! physical element. In addition it returns the integration scheme,
  //! which is already initialized to the order depending on the definition
  //! of the region of the element and the order of the element (if a 
  //! relative integration order was specified). 
  virtual BaseFE* GetFe( const EntityIterator ent ,
                         IntScheme::IntegMethod& method,
                         IntegOrder & order ) = 0;
  
  //! Return pointer to reference element (by element number)
  virtual BaseFE* GetFe( UInt elemNum ) = 0;

  //! Return number of functions for given entity
  virtual UInt GetNumFunctions( const EntityIterator ent );

  //! Obtain a pointer with interior surface elements for a given volume region
  //! \param (in) region The volume region Id
  //! \param (out) surfElems entity list with surface elements
  //! \param (out) opposingElems entitylist ordered such that the elements are
  //!                  direct neighbors of the first entity list
  virtual void GetInteriorSurfaceElems(RegionIdType region,
                               shared_ptr<EntityList> & surfElems,
                               shared_ptr<EntityList> & opposingElems){
    EXCEPTION("This FeSpace does not feature the dynamic generation of surface elements");
  }

  //! Obtain a pointer with interior surface elements for a given volume region
  //! \param (in) region The volume region Id
  //! \param (out) surfElems entity list with surface elements
  virtual void GetInteriorSurfaceElems(RegionIdType region,
                               shared_ptr<EntityList> & surfElems){
    EXCEPTION("This FeSpace does not feature the dynamic generation of surface elements");
  }

  //! Obtain a pointer with exterior surface elements for a given volume region
  //! NOTE: Be very careful, this entity list will contain the exterior surfaces
  //!       of all regions the space has created interoir surfaces for
  //! TODO: THINK OF ANOTHER NAME WHIHC REFLECTS THIS FACT
  //! \param (in) region The volume region Id
  //! \param (out) surfElems entity list with surface elements
  virtual void GetExteriorSurfaceElems(RegionIdType region, shared_ptr<EntityList> & surfElems){
    EXCEPTION("This FeSpace does not feature the dynamic generation of surface elements");
  }

  //! Return for a given lower dimensional element a valid volume element

  //! This method takes a lower dimensional element as argument and returns a 
  //! valid volume neighbor.
  //! A valid neighbor is considered to be a volume element
  //! which belongs to one of the regions the space / function is defined on.
  //! If more than one elements match, simply the first one is returned (which
  //! should be okay for continuous spaces), if none matches, an exception is
  //! thrown.
  //! \param[in] lowDimElem Pointer to surface element, which neighbor is seeked
  //! \return Pointer to neighboring volume element
  virtual const Elem* GetVolElem( const Elem* lowDimElem ) const;

  //@}
  
  // ========================================================================
  //  EQUATION HANDLING
  // ========================================================================
  //@{ \name Equation Handling
  
  //! Return equation numbers
  virtual void GetEqns( StdVector<Integer>& eqns, const EntityIterator ent ); 

  //! Return equation numbers for a specific dof
  virtual void GetEqns( StdVector<Integer>& eqns, const EntityIterator ent,
                        UInt dof );
  
  //! Return equation numbers for a specific dof and entitytype
  virtual void GetEqns( StdVector<Integer>& eqns, const EntityIterator ent,
                        UInt dof, BaseFE::EntityType ); 
  
  //! Return equation numbers for a specific entitytype
  virtual void GetEqns( StdVector<Integer>& eqns, const EntityIterator ent,
                        BaseFE::EntityType ); 

  //! Get Equation numbers for a specific element
  virtual void GetElemEqns(StdVector<Integer>& eqns, const Elem* elem);

  //! Get Equation numbers for a specific element
  virtual void GetElemEqns(StdVector<Integer>& eqns, const Elem* elem, UInt dof);

  //! Get equations for a given entityList
  virtual void GetEntityListEqns(StdVector<Integer>& eqns, 
                                 shared_ptr<EntityList> ent,
                                 BaseFE::EntityType = BaseFE::ALL );
  
  //! Get equations for a given entityList for a given dof
  virtual void GetEntityListEqns(StdVector<Integer>& eqns, 
                                 shared_ptr<EntityList> ent, UInt dof,
                                 BaseFE::EntityType = BaseFE::ALL );
  
  //! Return SBM-block and Matrix-SubBlock definition according to strategy
  
  //! This methods returns a list of SBMBlock definitions and minorBlock
  //! definitions, according to the solution strategy provided.
  //! \param sbmBlocks   Vector (length: numBlocks) with all equations in 
  //!                    this block(see AlgebraicSys)
  //! \param minorBlocks Map (key: sbmBlockIndex) with Vector (index:
  //!                    subBlockIndex) of set of equations for each minorBlock
  virtual void GetOlasMappings( StdVector<AlgebraicSys::SBMBlockDef>& sbmBlocks,
                                std::map<UInt,StdVector<std::set<Integer> > >&
                                minorBlocks );
  
  //! Add result
  virtual void AddFeFunction( shared_ptr<BaseFeFunction> fct );

  //! Precalculate integration points
  virtual void PreCalcShapeFncs() {};

  //! Get number of equations which are not fixed by BCs this space
  
  //! This method returns the number of equations, which are free,
  //! i.e. not fixed by a (non)homogeneous Dirichlet boundary
  //! condition. As the free equations are numbered first internally
  //! we can use this method to determine, if an equation number belongs
  //! to a free equations (<= numFreeEquations) or if it is a Dirichlet
  //! equations (>numFreeEquations)
  virtual UInt GetNumFreeEquations(){
    return numFreeEquations_;
  }

  //! Get number of equations this space has assigned
  virtual UInt GetNumEquations(){
    return numEqns_;
  }

  //! Returns the MapType of the equation numbering
  MappingType GetMapType(RegionIdType region){
    return mapType_;
  }

  UInt GetNumHomDirichletBc(){
    return bcCounter_[HDBC];
  }

  UInt GetNumInhomDirichletBc(){
    return bcCounter_[IDBC];
  }

  UInt GetNumConstraints(){
    return bcCounter_[CS];
  }
  
  //! Return number of unknowns per virtual "node"

  //! This method returns the number of components / unknowns per virtual 
  //! "node", which differs depending on the type of space: If a vectorial 
  //! function is approximated in a H1 space, every scalar component gets
  //! approximated by a scalar basis function. In case of a HCurl or HDiv 
  //! element, the shape functions are already vectorial, so only one
  //! coefficient is needed to describe the function value.
  virtual UInt GetNumDofs() const;

  weak_ptr<BaseFeFunction> GetFeFunction(){
    return feFunction_;
  }
    
  const SingleEqnMap& GetNodeMap() const { return nodeMap_; }

  //! Map equations i.e. intialize object
  virtual void Finalize() = 0;

  /** Dump information of the equation map to the console
   * @param file if NULL std::cout is used, otherwise give an fil ostream */
  virtual void PrintEqnMap(std::ostream* file = NULL);
  
  //! Update the FeSpace in case of a multistep solution strategy
  
  //! This method updates the internal data structure to reflect the a new step
  //! in a multistep simulation. E.g. for a twoLevel solution strategy, the
  //! active set of equations is switched between lowest order equations only 
  //! and the full set.
  virtual void UpdateToSolStrategy() {;}

  virtual void InsertElemsToCoilList(shared_ptr<ElemList> eL, shared_ptr<CoilList> cL){}

 //@}

  void SetLagrSurfSpace() 
  {
    lagrangeSurfSpace_ = true;
  }

  //! Map a general coefficient function onto a given FeFunction
  
  //! This method can be used to map a general coefficient function
  //! to the current finite element space. It returns a map, containing the 
  //! equations numbers and the corresponding coefficients.
  //! \param support Entitylists on which the function is defined
  //! \param coefFct Coefficient function to be mapped 
  //! \param feFct FeFunction onto which the coeffunction get mapped.
  //!              The feFunction must point to the same FeSpace.
  //! \param vals Map containing the equations numbers (key) and the
  //!             coefficient values (value)
  //! \param cache Flag, if mapping should be cached (e.g. for boundary
  //!              conditions, depending on frequency, time)
  //! \param comp Set containing the components, which should get mapped.
  //!             If empty, all components of the (vector-valued) function
  //!             get mapped
  //! \param updatedGeo Flag, if the coefFct lives on the updated geometry
  virtual void MapCoefFctToSpace(StdVector<shared_ptr<EntityList> > support, 
                                 shared_ptr<CoefFunction> coefFct,
                                 shared_ptr<BaseFeFunction> feFct,
                                 std::map<Integer, Double>& vals,
                                 bool cache,
                                 const std::set<UInt>& comp = std::set<UInt>(),
                                 bool updatedGeo = true)=0;

  virtual void MapCoefFctToSpace(StdVector<shared_ptr<EntityList> > support, 
                                 shared_ptr<CoefFunction> coefFct,
                                 shared_ptr<BaseFeFunction> feFct,
                                 std::map<Integer, Complex>& vals,
                                 bool cache,
                                 const std::set<UInt>& comp = std::set<UInt>(),
                                 bool updatedGeo = true)=0;
  
  //! Check if entity approximation is the same in other space
  
  //! This method checks, if the approximation  of an entity list
  //! (e.g. elements of a region) have the same approximation in this space
  //! compared to another one.
  //! In this case the coefficients in the related FeFunction correspond
  //! and can be copied.
  virtual bool IsSameEntityApproximation( shared_ptr<EntityList> list,
                                          shared_ptr<FeSpace> space ) = 0;


  //! Create a map of equation-index-geometry for AMG-solver/preconditioner
  //! ONLY for single-PDE's with lowest order Lagrange- or edge-elements
  void CreateEquIndGeomMap(std::unordered_map< Integer, BaseFeFunction::EqNodeGeom>&,
                          UInt&, UInt&);

protected:
  
  bool lagrangeSurfSpace_;  

  //! Parameter node
  PtrParamNode myParam_;
  
  //! Info node
  PtrParamNode infoNode_;
 
  //! Pointer to grid
  Grid* ptGrid_;
  
  //! Type of element space
  SpaceType type_;
  
  //! Stores if the space is in polynomial or grid mode
  //! Therefore, we set this variable to POLYNOMIAL as soon as
  //! some region needs e.g. higher order or edge elements
  MappingType mapType_;

  //! store the Polynomial of the space
  //! not very clean. in the near future we should make the spaces
  //! capable to store any kind of polynomial
  PolyType polyType_;

  //! Stores an offset to the element orders created
  //! is ignored in case of grid mapping
  UInt orderOffset_;

  //!Map for the parameter nodes for polyList for easier access
  std::map< std::string, PtrParamNode > integNodes_;

  //!Map for the parameter nodes for integrationList for easier access
  std::map< std::string, PtrParamNode > polyNodes_;

  //! Storing the FeFunctions associated with this space
  
  //! Note: Here we need a weak pointer to break the cyclic reference
  //! FeSpace <-> FeFunction
  weak_ptr<BaseFeFunction> feFunction_;
  
  //! Set with all regions the space is defined on
  std::set<RegionIdType> regions_;
  
  //! Pointer to solution strategy object
  shared_ptr<SolStrategy> solStrat_;

  //! Solution step (in case of multistep solution)
  UInt solStep_;
  
  //! Flag indicating if the FeSpace is already initialized
  bool isFinalized_;
  
  //! Flag indicating use of discontinuous (L2) functions.
  
  //! If false, we perform the numbering of nodes/edges/faces 
  //! element local.
  bool isContinuous_;

  //! Flag indicating use of hierarchical polynomials
  bool isHierarchical_;
  
  //! Number of equations administrated by this space
  UInt numEqns_;

  //! Number of equations administrated by this space not fixed by BCs
  UInt numFreeEquations_;

  //! map for storing the number of different boundary conditions
  std::map< BcType, UInt> bcCounter_;
  
  //! Pointer to integration scheme
  shared_ptr<IntScheme> intScheme_;

  // =====================================================
  // REGION SPECIFIC DATA
  // =====================================================

  //! Stores the Integration for each region
  //! if the user does not specify anything we just fill a default
  //! in the anisotrpopic case,the order matrix represents in each row, the component
  //! of (vectorial) unknowns and in each column the order for each space direction
  //! to distinguish cleanly one could define an isIsoOrderInt flag....
  std::map<RegionIdType, IntegDefinition > regionIntegration_;

  //! Maps polyregionIDs to integIds
  //! needed only for precalulation of shape function in FeSpaceH1Nodal right now.
  std::map<RegionIdType, std::set<RegionIdType> > polyToIntegMap;

  //! Set the order and mapping type of a specific region
  virtual void SetRegionElements( RegionIdType region, MappingType mType,
                                  const ApproxOrder& order,
                                  PtrParamNode infoNode )=0;

  //! read in integration data and set defaults
  virtual void SetRegionIntegration(RegionIdType region, 
                                    IntScheme::IntegMethod method,
                                    const IntegOrder& order,
                                    IntegOrderMode mode,
                                    PtrParamNode infoNode );
  
  // ====================================================================
  // INTERNAL INITIALIZATION
  // ====================================================================
  //! Here the spaces have the possibility to check if user definitions makes sense
  //! e.g. if the chosen integration is correct or the element order is nice
  //! here one could e.g. adjust the integration order according to the element order
  virtual void CheckConsistency() = 0;

  //! sets the default integration scheme and order
  virtual void SetDefaultIntegration(PtrParamNode infoNode ) = 0;

  //! Create default finite elements to be used if nothing else is requested
  virtual void SetDefaultElements(PtrParamNode infoNode) = 0;

  //! reads in special options for the space under consideration
  virtual void ReadCustomAttributes(PtrParamNode pNode,RegionIdType region){;}

  // =====================================================
  //  Nodal section
  // =====================================================
  
  //! A Function Creating the virtual Node Array, this method adds virtual
   //! node numbers to each element, thus making the equation numbering more
   //! convenient
   virtual void CreateVirtualNodes();

   virtual void CreatePolynomialNodes();

   virtual void CreateGridNodes();
   
   //! Get all Node Numbers according to the mapping GRID based or
   //! POLYNOMIAL Based according to the Type requested
   virtual void GetNodesOfEntities( StdVector<UInt>& nodes,
                                    shared_ptr<EntityList> ent,
                                    BaseFE::EntityType entType = BaseFE::ALL);

   //! Get all Node Numbers of a specific element according to the mapping GRID based or
   //! POLYNOMIAL Based according to the Type requested
   virtual void GetNodesOfElement( StdVector<UInt>& nodes,
                                   const Elem* ptElem,
                                   BaseFE::EntityType entType = BaseFE::ALL);
   
  
  //! Associates GridNodeNumbers to virtual node numbers
  //! Needed e.g. for quadratic Grids in combination with NodalBCs and linear,
  //! Polynomial approximation
  std::map<UInt,StdVector<UInt> > gridToVirtualNodes_;

  //! Map for every node the type geometric entitytype it belongs to (V/E/F/I)
  typedef boost::unordered_map<UInt, BaseFE::EntityType> NodeTypeMap;
  NodeTypeMap nodesType_;

  
  //! Auxilliary map for assigning an enumerable entity a list of nodes
  
  //! This map type can be used to assign an enumerable geometric 
  //! entity (e.g. node/edge/face/element number, used as key in the map)
  //! a list of virtual nodes (value of map).
  typedef boost::unordered_map<UInt, StdVector<UInt> > EntityNodesType;
  
  //! Global map for continuous virtual node numbers
  
  //! This map contains all continuously numbered virtual nodes. For each
  //! type of entity (key, e.g. VERTEX, EDGE, FACE, INTERIOR) it
  //! contains a map (see EntityNodesType) which globally holds for each
  //! unique entity number (e.g. node/edge/face/interior number) a 
  //! list of virtual Nodes
  boost::unordered_map<BaseFE::EntityType, EntityNodesType> vNodesCont_;
  
  // ====================================================================
  // CACHE ACCESS TO LAST USED EQUATIONS
  // ====================================================================
  //@{
  //! Last accessed element number (only for element access)
  CfsTLS<UInt> lastElemNum_;

  //! Last accessed equations numbers
  CfsTLS< StdVector<Integer> > lastEqns_;

  //! Temporary vector for equation numbers
  CfsTLS< StdVector<UInt> > eqnNodes_;

  //@}

  // ====================================================================
  // Equation Map
  // ====================================================================
  //! Map Nodal BC Equation NUmbers
  virtual void MapNodalBCs() = 0;
  
  //! Map Nodal Equation Numbers
  virtual void MapNodalEqns(UInt phase);
  
  //! Nodal Equation Map
  SingleEqnMap nodeMap_;
  
  // =========================================================
  // READ USER INPUT
  // =========================================================

  //! Read the contents of the given parameter node and call the SetRegionIntegration Function
  virtual void ReadIntegNode(PtrParamNode node,  IntScheme::IntegMethod & iMeth,
                             IntegOrder& order, IntegOrderMode & mode);

  //! Here we pass a fePolynomial node such that the feSpace can extract the information
  //! which is important for the specific space
  virtual void ReadPolyNode(PtrParamNode node, MappingType & mapType,
                            ApproxOrder & order);

  //! Read the IntegrationSchemeList from the xml and put the nodes in the map
  //! according to their id's
  virtual void ReadIntegList();

  //! Read the fePolynomialList from the xml and put the nodes in the map
  //! according to their id's
  //! This makes it easier later on to assign the correct elements
  virtual void ReadPolyList();

  // =========================================================
  // MISCELLANEOUS
  // =========================================================
  
  //! Get set of all elements of this space
  virtual void GetAllElems(std::set<const Elem*>& allElems );
};


}


  // namespace CoupledField
#endif
