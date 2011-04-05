// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
#ifndef FILE_CFS_FESPACE_HH
#define FILE_CFS_FESPACE_HH

#include "General/environment.hh"
#include "Domain/entityList.hh"
#include "Domain/Composite.hh"
#include "Domain/elem.hh"
#include "General/exception.hh"
#include "Domain/domain.hh"
#include "Elements/basefe.hh"
//#include "Elements/2D/quad1fe.hh"
//#include "Elements/3D/hexa1FE.hh"
//#include "Elements/1D/line1fe.hh"
#include "Domain/resultInfo.hh"
#include "Elements/fefunction.hh"
#include "MatVec/matrix.hh"
#include "Elements/integrationScheme.hh"

namespace CoupledField {


// forward class declarations
class BaseFe;
class BaseFeFunction;

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
                       polynomials (e.g. Legendre, Jacobi, Gegenbauer).
  
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

  //! Enum which stores the (Virtual) Nodes of an element according to
  //! thier definition on vertices,edges,faces and interor
  typedef std::map< BaseFE::EntityType , StdVector<UInt> > ElemVirtualNodes;

  //! Constructor
  FeSpace(ParamNode * paramNode);

  //! Destructor
  ~FeSpace();
  
  //! Read parameters from ParamNode and initialize
  virtual void Init() {};

  //! Return type of FeSpace
  SpaceType GetSpaceType() { return type_;}
  
  //! Generate instance of specific element space type
  static shared_ptr<FeSpace> CreateInstance( ParamNode* aNode  );
  
  // ========================================================================
  //  INITIALIZATION 
  // ========================================================================
  //@{ \name Initialization
  
    //! Set isotropic order of functions to be used
  virtual void SetIsoOrder( UInt order );
  
  //! Return, if the space is hierarchical
  bool IsHierarchical() { return isHierarchical_;};
  
  //! Set Solution strategy and solution step
  virtual void SetStrategy( SolStrategyType strategy, UInt step ) {}; 
  
  //@}
  
  // ========================================================================
  //  ELEMENT HANDLING
  // ========================================================================
  //@{ \name Element Handling

  //! Return pointer to reference element
  virtual BaseFE* GetFe( const EntityIterator ent ) = 0;
  
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

  //! Reorder the equation Map (just for comptibility)
  virtual void ReorderEqnMap( StdVector<UInt> newOrder ) = 0;

  //! Add result
  virtual void AddFeFunction( shared_ptr<BaseFeFunction> fct ) = 0;

  //! Precalculate integration points
  virtual void PreCalcShapeFncs() {};

  //! Get number of equaitons thich are not fixed by BCs this space has assinged
  virtual UInt GetNumFreeEquations(){
    // in this approach we assume the penalty approach towards
    // solution. otherwise we would need to return only the free equations
    return numEqns_;
    //return numFreeEquations_;
  }

  
  //! Get number of equaitons this space has assinged
  virtual UInt GetNumEquations(){
    return numEqns_;
  }

  //! Sets the MapType of the equation numbering to GridBased or PolynomialBased Mapping
  virtual void SetMapType(MappingType type){
    mapType_ = type;
  }

  //! Returns the MapType of the equation numbering
  MappingType GetMapType(){
    return mapType_;
  }

  //! Get number of (vectorial) unknowns this space has assinged
  //! Including those, of which one Dof might be fixed by BCs
  //! Deprecated. Just for compatibility with NodeStoreSolution Class
  UInt GetNumUnknowns(){
    return numUnknowns_;
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

  //! Get polynomial order per entity per dof
  
  //! This method returns the number of unknowns per entitytype (NODE / EDGE / FACE / INTERIOR)
  //! per component of the unknown if its a vectorial unkown. If the order is isotropic
  //! (e.g. Lagrangian space), this method returns always the same number. 
  //! In the case of anisotropic polynomial orders (e.g. hierarchical FE) this is more involved
  virtual UInt GetEntityOrder( UInt elemNum, BaseFE::EntityType type, 
                               UInt entityNum, UInt comp = 1 ) = 0;
  
  //! Get maximum polynomial orde per enntity (maximum over all components) 
  virtual UInt GetMaxEntityOrder( UInt elemNum, BaseFE::EntityType type, 
                                  UInt entityNum ) = 0;
    
  //! Map equations i.e. intialize object
  virtual void Finalize() = 0;

  //! Dump information of the equation map to the console
  virtual void PrintEqnMap() = 0;

  //@}
  
protected:

  //! Parameter node
  ParamNode * myParam_;
  
  //! Type of element space
  SpaceType type_;
  
  //! Map for reference elmenets
  std::map<Elem::FEType, BaseFE* > refElems_;

  //! Storing the FeFunctions associated with this space
  shared_ptr<BaseFeFunction> feFunction_;

  //! Solution strategy to be used
  SolStrategyType solStrategy_;
  
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
  
  //! Isotropic polynomial approximation order
  UInt isoOrder_;
  
  //! Number of equations administrated by this space
  UInt numEqns_;

  //! Number of equations administrated by this space not fixed by BCs
  UInt numFreeEquations_;

  //! Number of (vectorial) unknowns administrated by this space
  //! Including those, of which one Dof might be fixed by BCs
  //TODO: Perhaps think of another name. Unknowns is a little misleading
  UInt numUnknowns_;
  
  //! map for storing the number of different boundary conditions
  std::map< BcType, UInt> bcCounter_;

  //! sotres if equation numbering is grid based or order based
  MappingType mapType_;

  // =====================================================
  //  Nodal section
  // =====================================================
  
  //! A Function Creating the virtual Node Array, this method adds virtual
   //! node numbers to each element, thus making the equation numbering more
   //! convenient
   virtual void CreateVirtualNodes();

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
  std::map<UInt,UInt> gridToVirtualNodes_;

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
};


}


  // namespace CoupledField
#endif
