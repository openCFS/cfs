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

namespace CoupledField {


// forward class declarations
class BaseFe;
class BaseFeFunction;
class SolStrategy;

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
  typedef enum { UNDEF_SPACE, CONST, H1, HCURL, HDIV, L2 } SpaceType;
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

  
  //! Struct containing all virtual nodes for on entity type (Vertex, Face
  struct EntityTypeNodes {
    //! Nodes for all numbers of entity (edgeNodes, faceNodes, innerNodes)
    StdVector<UInt> vNodes;
    //! Offset to vNodes array
    StdVector<UInt> offset;
  };
  
  //! Enum which stores the (Virtual) Nodes of an element according to
    //! their definition on vertices,edges,faces and interior
  typedef std::map< BaseFE::EntityType , EntityTypeNodes> ElemVirtualNodes;

  typedef enum {ABSOLUTE,RELATIVE} IntegOrderMode;
  static Enum<IntegOrderMode> IntegOrderModeEnum;

  struct IntegDefinition{
      IntScheme::IntegMethod method;
      Matrix<Integer> order;
      IntegOrderMode mode;
  };

  //! Constructor
  FeSpace(PtrParamNode paramNode, PtrParamNode infoNode);

  //! Destructor
  ~FeSpace();
  
  //! Read parameters from ParamNode and initialize
  virtual void Init() {};

  //! Return type of FeSpace
  SpaceType GetSpaceType() { return type_;}
  
  //! Generate instance of specific element space type
  static shared_ptr<FeSpace> CreateInstance(PtrParamNode aNode, 
                                            PtrParamNode infoNode,
                                            SpaceType reqType  );
  
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
                      IntScheme::IntegMethod & method,Matrix<Integer> & order);

  //! Set the approximation for the given region according to the passed poly and integ Ids
  //! ToDO if integration and polynomial nodes are merged, of course only an orderId has to be passed
  void SetRegionApproximation(RegionIdType region, std::string polyId, std::string integId);

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
                         shared_ptr<IntScheme>& intScheme ) = 0;
  
  //! Return pointer to reference element (by element number)
  virtual BaseFE* GetFe( UInt elemNum ) = 0;

  //! Return number of functions for given entity
  virtual UInt GetNumFunctions( const EntityIterator ent ) = 0;
  //@}
  
  // ========================================================================
  //  EQUATION HANDLING
  // ========================================================================
  //@{ \name Equation Handling
  
  //! Return equation numbers
  virtual void GetEqns( StdVector<Integer>& eqns, const EntityIterator ent ) = 0; 

  //! Return equation numbers for a specific dof
  virtual void GetEqns( StdVector<Integer>& eqns, const EntityIterator ent
                        , UInt dof ) = 0;
  
  //! Return equation numbers for a specific dof and entitytype
  virtual void GetEqns( StdVector<Integer>& eqns, const EntityIterator ent
                        , UInt dof, BaseFE::EntityType ) = 0; 

  //! Get a Nodal Equation number
  virtual UInt GetNodeEqn(UInt nodeNr, UInt dof) = 0;

  //! Get Equation numbers for a specific element
  virtual void GetElemEqns(StdVector<Integer>& eqns,const Elem* elem) = 0;

  //! Get Equation numbers for a specific element
  virtual void GetElemEqns(StdVector<Integer>& eqns,const Elem* elem, UInt dof) = 0;

  //! Return SBM-block and Matrix-SubBlock definition according to strategy
  
  //! This methods returns a list of SBMBlock definitions and minorBlock
  //! definitions, according to the solution strategy provided.
  //! \param solStrat    Solution strategy object, which defines splitting
  //! \param sbmBlocks   Vector (length: numBlocks) with all equations in 
  //!                    this block(see AlgebraicSys)
  //! \param minorBlocks Map (key: sbmBlockIndex) with Vector (index:
  //!                    subBlockIndex) of set of equations for each minorBlock
  virtual void GetOlasMappings( shared_ptr<SolStrategy> solStrat, 
                                StdVector<AlgebraicSys::SBMBlockDef>& sbmBlocks,
                                std::map<UInt,StdVector<std::set<Integer> > >&
                                minorBlocks ) = 0;
  
  //! Add result
  virtual void AddFeFunction( shared_ptr<BaseFeFunction> fct ) = 0;

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

  shared_ptr<BaseFeFunction> GetFeFunction(){
    return feFunction_;
  }

  //! Get polynomial order per entity per dof
  
  //! This method returns the number of unknowns per entitytype (NODE / EDGE / FACE / INTERIOR)
  //! per component of the unknown if its a vectorial unkown. If the order is isotropic
  //! (e.g. Lagrangian space), this method returns always the same number. 
  //! In the case of anisotropic polynomial orders (e.g. hierarchical FE) this is more involved
  virtual UInt GetEntityOrder( UInt elemNum, BaseFE::EntityType type, 
                               UInt entityNum, UInt comp = 1 ) = 0;
  
  //! Get maximum polynomial orde per entity (maximum over all components) 
  virtual UInt GetMaxEntityOrder( UInt elemNum, BaseFE::EntityType type, 
                                  UInt entityNum ) = 0;
    
  //! Map equations i.e. intialize object
  virtual void Finalize() = 0;

  //! Dump information of the equation map to the console
  virtual void PrintEqnMap() = 0;
  
  //@}
  
protected:

  //! Parameter node
  PtrParamNode myParam_;
  
  //! Infor node
  PtrParamNode infoNode_;
  
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
  shared_ptr<BaseFeFunction> feFunction_;

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
                                  const Matrix<Integer>& order,
                                  PtrParamNode infoNode )=0;

  //! read in integration data and set defaults
  virtual void SetRegionIntegration(RegionIdType region, 
                                    IntScheme::IntegMethod method,
                                    Matrix<Integer> order,
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

  //! Stores every assigned virtual node
  StdVector<UInt> nodes_;
  
  //! Map for every node the type geometric entitytype it belongs to (V/E/F/I)
  std::map<UInt, BaseFE::EntityType> nodesType_;

  //! This is the virtual node Map for standard element it just contains
  //! the connectivity of the element, for higher order elements it contains also 
  //! the virtual node numbers in the correct ordering
  //! This Variable could be extended to store also the coordinates of all nodes
  //! created
  std::map< UInt, ElemVirtualNodes > virtualNodes_;

  // =========================================================
  // READ USER INPUT
  // =========================================================

  //! Read the contents of the given parameter node and call the SetRegionIntegration Function
  virtual void ReadIntegNode(PtrParamNode node,  IntScheme::IntegMethod & iMeth,
                             Matrix<Integer> &orderMat, IntegOrderMode & mode);

  //! Here we pass a fePolynomial node such that the feSpace can extract the information
  //! which is important for the specific space
  virtual void ReadPolyNode(PtrParamNode node, MappingType & mapType,
                            Matrix<Integer> & order);

  //! Read the IntegrationSchemeList from the xml and put the nodes in the map
  //! according to their id's
  virtual void ReadIntegList();

  //! Read the fePolynomialList from the xml and put the nodes in the map
  //! according to their id's
  //! This makes it easier later on to assign the correct elements
  virtual void ReadPolyList();

};


}


  // namespace CoupledField
#endif
