// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
//=================================
/*
 * \file   FeSpaceConst.cc
 * \brief  This FeSpace provides an interface for adding single equations to the system,
 *         e.g. equations of an equivalent circuit of lumped elements.
 *
 * \date   May 6, 2014
 * \author dperchto
 */
//=================================

#ifndef FESPACECONST_HH
#define FESPACECONST_HH

#include "FeBasis/FeSpace.hh"

namespace CoupledField {

class FeSpaceConst : public FeSpace {
public:
  
  //! Constructor
  FeSpaceConst(PtrParamNode paramNode, PtrParamNode infoNode, Grid* ptGrid );

  //! Destructor
  virtual ~FeSpaceConst();
  
  //! Read parameters from ParamNode and initialize
  
  //! Reads the parameters from a ParamNode and initializes itself,
  //! incorporating the information from a SolStrategy object.
  //! \param solStrategy Solution strategy object, which encapsulates
  //!                    additional information (splitting, number of
  //!                    solution steps  etc.)
  virtual void Init( shared_ptr<SolStrategy> solStrat );

  //! Return type of FeSpace
  //SpaceType GetSpaceType() { return type_;}
  
  //! Generate instance of specific element space type
  /*static shared_ptr<FeSpace> CreateInstance(PtrParamNode aNode,
                                            PtrParamNode infoNode,
                                            SpaceType reqType,
                                            Grid* ptGrid );*/
  
  // ========================================================================
  //  INITIALIZATION 
  // ========================================================================
  //@{ \name Initialization

  //! Sets an offset to element orders. Right now this is only needed for
  //! Mixed TaylorHood approximation
  //void SetElementOrderOffset(UInt offset); // according to Ctrl+Shift+T this is not needed anywhere

  //! Return, if the space is hierarchical
  //bool IsHierarchical() { return isHierarchical_;};
  
  //@}
  

  // ========================================================================
  // Multiple Region handling
  // ========================================================================
  //@{ \name Region Handling

  //! Get the Integration method and its order for the current region
  /*void GetIntegration(BaseFE * fe, RegionIdType region,
                      IntScheme::IntegMethod & method,
                      IntegOrder & order);*/ // maybe the PDE should set the same integration as the main unknown

  //! Return integration scheme
  //shared_ptr<IntScheme> GetIntScheme() { return intScheme_;}
  
  //! Set the approximation for the given region according to the passed poly and integ Ids
  //! ToDO if integration and polynomial nodes are merged, of course only an orderId has to be passed
  //void SetRegionApproximation(RegionIdType region, std::string polyId, std::string integId);

  //! Sets a default approximation for all regions
  //void SetDefaultRegionApproximation();

  //@}
  // ========================================================================
  //  ELEMENT HANDLING
  // ========================================================================
  //@{ \name Element Handling

  //! Return pointer to reference element
  virtual BaseFE* GetFe( const EntityIterator ent );
  
  //! Return pointer to reference element and adjusted integrationScheme
  
  //! This method returns the pointer to the reference according to the
  //! physical element. In addition it returns the integration scheme,
  //! which is already initialized to the order depending on the definition
  //! of the region of the element and the order of the element (if a 
  //! relative integration order was specified). 
  virtual BaseFE* GetFe( const EntityIterator ent ,
                         IntScheme::IntegMethod& method,
                         IntegOrder & order );
  
  //! Return pointer to reference element (by element number)
  virtual BaseFE* GetFe( UInt elemNum );

  //! Return number of functions for given entity
  virtual UInt GetNumFunctions( const EntityIterator ent );

  //! Obtain a pointer with interior surface elements for a given volume region
  //! \param (in) region The volume region Id
  //! \param (out) surfElems entity list with surface elements
  //! \param (out) opposingElems entitylist ordered such that the elements are
  //!                  direct neighbors of the first entity list
  /*virtual void GetInteriorSurfaceElems(RegionIdType region,
                               shared_ptr<EntityList> & surfElems,
                               shared_ptr<EntityList> & opposingElems){
    EXCEPTION("This FeSpace does not feature the dynamic generation of surface elements");
  }*/

  //! Obtain a pointer with interior surface elements for a given volume region
  //! \param (in) region The volume region Id
  //! \param (out) surfElems entity list with surface elements
  /*virtual void GetInteriorSurfaceElems(RegionIdType region,
                               shared_ptr<EntityList> & surfElems){
    EXCEPTION("This FeSpace does not feature the dynamic generation of surface elements");
  }*/

  //! Obtain a pointer with exterior surface elements for a given volume region
  //! NOTE: Be very careful, this entity list will contain the exterior surfaces
  //!       of all regions the space has created interoir surfaces for
  //! TODO: THINK OF ANOTHER NAME WHIHC REFLECTS THIS FACT
  //! \param (in) region The volume region Id
  //! \param (out) surfElems entity list with surface elements
  /*virtual void GetExteriorSurfaceElems(RegionIdType region, shared_ptr<EntityList> & surfElems){
    EXCEPTION("This FeSpace does not feature the dynamic generation of surface elements");
  }*/

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
  //virtual const Elem* GetVolElem( const Elem* lowDimElem ) const;

  //@}
  
  // ========================================================================
  //  EQUATION HANDLING
  // ========================================================================
  //@{ \name Equation Handling
  
  //! Return equation numbers
  virtual void GetEqns( StdVector<Integer>& eqns, const EntityIterator ent ); 

  //! Return equation numbers for a specific dof
  virtual void GetEqns( StdVector<Integer>& eqns, const EntityIterator ent
                        , UInt dof );
  
  //! Return equation numbers for a specific dof and entitytype
  virtual void GetEqns( StdVector<Integer>& eqns, const EntityIterator ent,
                        UInt dof, BaseFE::EntityType ); 
  
  //! Return equation numbers for a specific entitytype
  virtual void GetEqns( StdVector<Integer>& eqns, const EntityIterator ent,
                        BaseFE::EntityType ); 

  //! Get Equation numbers for a specific element
  virtual void GetElemEqns(StdVector<Integer>& eqns,const Elem* elem);

  //! Get Equation numbers for a specific element
  virtual void GetElemEqns(StdVector<Integer>& eqns,const Elem* elem, UInt dof);

  // Following 2 methods perhaps not needed (if base class calls the GetEqns of FeSpaceConst)

  //! Get equations for a given entityList
  /*virtual void GetEntityListEqns(StdVector<Integer>& eqns,
                                 shared_ptr<EntityList> ent,
                                 BaseFE::EntityType = BaseFE::ALL );
  
  //! Get equations for a given entityList for a given dof
  virtual void GetEntityListEqns(StdVector<Integer>& eqns, 
                                 shared_ptr<EntityList> ent, UInt dof,
                                 BaseFE::EntityType = BaseFE::ALL );*/
  
  //! Return SBM-block and Matrix-SubBlock definition according to strategy
  
  //! This methods returns a list of SBMBlock definitions and minorBlock
  //! definitions, according to the solution strategy provided.
  //! \param sbmBlocks   Vector (length: numBlocks) with all equations in 
  //!                    this block(see AlgebraicSys)
  //! \param minorBlocks Map (key: sbmBlockIndex) with Vector (index:
  //!                    subBlockIndex) of set of equations for each minorBlock
  /*virtual void GetOlasMappings( StdVector<AlgebraicSys::SBMBlockDef>& sbmBlocks,
                                std::map<UInt,StdVector<std::set<Integer> > >&
                                minorBlocks );*/
  
  //! Add result
  //virtual void AddFeFunction( shared_ptr<BaseFeFunction> fct );

  //! Precalculate integration points
  //virtual void PreCalcShapeFncs() {};

  //! Get number of equations which are not fixed by BCs this space
  
  //! This method returns the number of equations, which are free,
  //! i.e. not fixed by a (non)homogeneous Dirichlet boundary
  //! condition. As the free equations are numbered first internally
  //! we can use this method to determine, if an equation number belongs
  //! to a free equations (<= numFreeEquations) or if it is a Dirichlet
  //! equations (>numFreeEquations)
  /*virtual UInt GetNumFreeEquations(){
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
  }*/
  
  //! Return number of unknowns per virtual "node"

  //! This method returns the number of components / unknowns per virtual 
  //! "node", which differs depending on the type of space: If a vectorial 
  //! function is approximated in a H1 space, every scalar component gets
  //! approximated by a scalar basis function. In case of a HCurl or HDiv 
  //! element, the shape functions are already vectorial, so only one
  //! coefficient is needed to describe the function value.
  //virtual UInt GetNumDofs() const;

  /*weak_ptr<BaseFeFunction> GetFeFunction(){
    return feFunction_;
  }*/
    
  //! Map equations i.e. intialize object
  virtual void Finalize();

  //! Dump information of the equation map to the console
  //virtual void PrintEqnMap();
  
  //! Update the FeSpace in case of a multistep solution strategy
  
  //! This method updates the internal data structure to reflect the a new step
  //! in a multistep simulation. E.g. for a twoLevel solution strategy, the
  //! active set of equations is switched between lowest order equations only 
  //! and the full set.
  //virtual void UpdateToSolStrategy() {;}
 //@}

  /*void SetLagrSurfSpace()
  {
    lagrangeSurfSpace_ = true;
  }*/


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
  virtual void MapCoefFctToSpace(StdVector<shared_ptr<EntityList> > support, 
                                 shared_ptr<CoefFunction> coefFct,
                                 shared_ptr<BaseFeFunction> feFct,
                                 std::map<Integer, Double>& vals,
                                 bool cache,
                                 const std::set<UInt>& comp = std::set<UInt>() );

  virtual void MapCoefFctToSpace(StdVector<shared_ptr<EntityList> > support, 
                                 shared_ptr<CoefFunction> coefFct,
                                 shared_ptr<BaseFeFunction> feFct,
                                 std::map<Integer, Complex>& vals,
                                 bool cache,
                                 const std::set<UInt>& comp = std::set<UInt>() );
  
  //! Check if entity approximation is the same in other space
  
  //! This method checks, if the approximation  of an entity list
  //! (e.g. elements of a region) have the same approximation in this space
  //! compared to another one.
  //! In this case the coefficients in the related FeFunction correspond
  //! and can be copied.
  virtual bool IsSameEntityApproximation( shared_ptr<EntityList> list,
                                          shared_ptr<FeSpace> space );

protected:
  
  //! can be removed if this space is generalized
  std::set<EntityList::ListType> allowedEntities_;

  /*bool lagrangeSurfSpace_;

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
  shared_ptr<IntScheme> intScheme_;*/

  // =====================================================
  // REGION SPECIFIC DATA
  // =====================================================

  //! Stores the Integration for each region
  //! if the user does not specify anything we just fill a default
  //! in the anisotrpopic case,the order matrix represents in each row, the component
  //! of (vectorial) unknowns and in each column the order for each space direction
  //! to distinguish cleanly one could define an isIsoOrderInt flag....
  //std::map<RegionIdType, IntegDefinition > regionIntegration_;

  //! Maps polyregionIDs to integIds
  //! needed only for precalulation of shape function in FeSpaceH1Nodal right now.
  //std::map<RegionIdType, std::set<RegionIdType> > polyToIntegMap;

  //! Set the order and mapping type of a specific region
  virtual void SetRegionElements( RegionIdType region, MappingType mType,
                                  const ApproxOrder& order,
                                  PtrParamNode infoNode );

  //! read in integration data and set defaults
  /*virtual void SetRegionIntegration(RegionIdType region,
                                    IntScheme::IntegMethod method,
                                    const IntegOrder& order,
                                    IntegOrderMode mode,
                                    PtrParamNode infoNode );*/
  
  // ====================================================================
  // INTERNAL INITIALIZATION
  // ====================================================================
  //! Here the spaces have the possibility to check if user definitions makes sense
  //! e.g. if the chosen integration is correct or the element order is nice
  //! here one could e.g. adjust the integration order according to the element order
  virtual void CheckConsistency();

  //! sets the default integration scheme and order
  virtual void SetDefaultIntegration(PtrParamNode infoNode );

  //! Create default finite elements to be used if nothing else is requested
  virtual void SetDefaultElements(PtrParamNode infoNode);

  //! reads in special options for the space under consideration
  //virtual void ReadCustomAttributes(PtrParamNode pNode,RegionIdType region){;}

  // =====================================================
  //  Nodal section
  // =====================================================
  
  //! A Function Creating the virtual Node Array, this method adds virtual
   //! node numbers to each element, thus making the equation numbering more
   //! convenient
   /*virtual void CreateVirtualNodes();

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
                                   BaseFE::EntityType entType = BaseFE::ALL);*/
   
  
  //! Associates GridNodeNumbers to virtual node numbers
  //! Needed e.g. for quadratic Grids in combination with NodalBCs and linear,
  //! Polynomial approximation
  /*std::map<UInt,StdVector<UInt> > gridToVirtualNodes_;

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
  boost::unordered_map<BaseFE::EntityType, EntityNodesType> vNodesCont_;*/
  
  // ====================================================================
  // Equation Map
  // ====================================================================
  //! Map Nodal BC Equation NUmbers
  virtual void MapNodalBCs();
  
  //! Map Nodal Equation Numbers
  //virtual void MapNodalEqns(UInt phase);
  
  //! Nodal Equation Map
  //SingleEqnMap nodeMap_;
  
  // =========================================================
  // READ USER INPUT
  // =========================================================

  //! Read the contents of the given parameter node and call the SetRegionIntegration Function
  /*virtual void ReadIntegNode(PtrParamNode node,  IntScheme::IntegMethod & iMeth,
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
  virtual void ReadPolyList();*/

  // =========================================================
  // MISCELLANEOUS
  // =========================================================
  
  //! Get set of all elements of this space
  //virtual void GetAllElems(std::set<const Elem*>& allElems );

private:

  //! currently the space is restricted to certain entities
  void CheckEntityType(const EntityIterator ent) const;

  //! maps equation number to entity id provided by EntityIterator::GetIdString()
  boost::unordered_map<std::string,Integer> equationMap_;

};

}

#endif
