#ifndef FILE_CFS_EQNMAP_HH
#define FILE_CFS_EQNMAP_HH


#include "Utils/StdVector.hh"
#include "Domain/resultInfo.hh"
#include "Domain/bcs.hh"

namespace CoupledField {

  // Forward class declarations
  class EntityList;
  class EntityIterator;
  class ResultInfo;
  
  //! Class for mapping entities and ansatz functions to equation numbers
  class EqnMap {
    
  public:
    
    //! Typedef for mapping results to entitylists
    typedef std::map<ResultInfo, StdVector<shared_ptr<EntityList> > > ResultEntityMap;
    typedef std::map<ResultInfo, StdVector<shared_ptr<ElemList> > > ResultElemListMap;
    typedef std::map<ResultInfo, StdVector<shared_ptr<NodeList> > > ResultNodeListMap;
    typedef std::map<ResultInfo, Matrix<Integer> > EqnMapType;
    typedef std::map<ResultInfo, StdVector<Vector<Integer> > > VecEqnMapType;
    typedef std::map<ResultInfo, HdBcList> ResultHdBcMap;
    typedef std::map<ResultInfo, IdBcList> ResultIdBcMap;
    typedef std::map<ResultInfo, ConstraintList> ResultConstraintMap;
    
    //! Constructor
    EqnMap( Grid* ptGrid, PdeIdType, bool sortEqns);
    
    //! Destructor
    ~EqnMap();
    
    // ======================================================================
    // SET AND INITIALIZATION METHODS
    // ======================================================================
    
    //@{
    //! \name Set And Initialization Methods
    
    //! Add a new result
    void AddResult( ResultInfo& result, 
                    shared_ptr<EntityList> list );
    
    //! Set the homogeneous boundary conditions
    void SetHomoDirichletBCs( HdBcList& hdBcs );
        
    //! Set the in-homogeneous boundary conditions
    void SetInhomDirichletBCs( IdBcList& idBcs );

    //! Set the constraint conditions
    void SetConstraints( ConstraintList& constraints );
    
    //! Finalize setup and calculate mapping
    void Finalize();

    //! Maps the equation numbers according to the reordering

    //! If the algebraic system has re-ordered the equations, e.g. to improve
    //! the performance of a sparse direct solver, this method must be used to
    //! incorporate the resulting change into the dof to eqn map.
    //! \param order permutation vector containing the re-ordering, i.e.
    //!              order[k-1] is the new equation number for the dof that
    //!              was previously associated to the k-th equation number.
    //! \note
    //! - It is safe to pass a NULL pointer to the method. In this case
    //!   the method assumes that no re-ordering was performed by the
    //!   algebraic system.
    //! - Since the re-ordering is private property of the equation data
    //!   object, we re-set order to NULL in this method.
    //! - For the sake of a low memory foot-print we delete the permutation
    //!   buffer once the reordering was incorporated
    void ReorderMapping( Integer **order );
    
    //@}
    // ======================================================================
    // GENERAL QUERY METHODS
    // ======================================================================

    //@{
    //! \name General Query Methods
    //! Return true if eqns are assigned
    inline bool IsFinalized() const {return isFinalized_;}

    //! Return the total number of equations
    inline UInt GetNumEqns() const { return numEqns_; }

    //! Return the equation number of the last unfixed degree of freedom

    //! This method returns the equation number of the last free degree of
    //! freedom. We consider a degree of freedom to be free, if it is not
    //! fixed by either
    //! - homogeneous Dirichlet boundary conditions
    //! - inhomogeneous Dirichlet boundary conditions
    //! - constraints (master - slave dof relations)
    //!
    //! The return value actually depends on the status of the sortEqns_ flag.
    //! If the latter is true, then we return numRealEqns_. If it is false,
    //! then inhomogeneous Dirichlet boundary conditions are treated by the
    //! penalty method and must be considered free dofs as well. Thus, in this
    //! case we return numEqns_.
    inline UInt GetNumLastFreeDof() const {
      return ( ( sortEqns_ == true ) ? numRealEqns_ : numEqns_ );
    }
                              
    //! Return number of real inhomogeneous Dirichlet boundary conditions

    //! This method returns the number of 'real' boundary conditions, i.e.
    //! douplets of nodes are already removed, so this number represents
    //! the number of equations, which are not 'free'
    inline UInt GetNumInHomDirichletEqns() const {
      return numIdBcs_;
    }
                              
    //! Return id of associated PDE
    PdeIdType GetPdeId() const {
      return pdeId_;
    }
    
    //@}

    // ======================================================================
    // EQUATION MAPPING
    // ======================================================================
    
    //@{
    //! \name Equation Mapping
    //! Map result and entities to equation numbers
    void GetEqns( StdVector<Integer>& eqns,
                  const ResultInfo& result, const EntityIterator& it ) const;
    
    //! Mpa combination of result, entity and dof to single equation number
    Integer GetEqn( const ResultInfo& result, 
                    const EntityIterator& it, UInt dof ) const;

    //! Name mapping function for obtaining a nodal equation
    Integer GetNodeEqn( const ResultInfo& result, UInt nodeNr, UInt dof );

    //! Name mapping function for obtaining a nodal equation
    Integer GetNodeEqn( UInt nodeNr, UInt dof );

    void GetNodeEqn( const StdVector<UInt>& nodeNrs, 
                     StdVector<Integer>& eqnNrs );

    //@}
    
    // ======================================================================
    // LOCAL/GLOBAL MAPPING OF MESH ENTITIES
    // ======================================================================
    
    //@{ 
    //! \name Local/Global Mapping of Mesh Entities

    //! Get number of local nodes
    inline UInt GetNumLocalNodes() const {return numLocNodes_;}
    
    //! Get number of local elements
    inline UInt GetNumLocalElems() const {return numLocElems_;}

    //! Get number of local edges
    inline UInt GetNumLocalEdges() const {return numLocEdges_;}

    //! Get number of local faces
    inline UInt GetNumLocalFaces() const {return numLocFaces_;}

    //! Map global to local node numbers
    //! (needed for nodal displacement of grid)
    void Mesh2PdeNode(StdVector<UInt> & PdeNodes,
                      const StdVector<UInt> & MeshNodes) const;

    //! Map global to local node number
    //! (needed for nodal displacement of grid)
    UInt Mesh2PdeNode(const UInt meshNode) const;
    
    //! Map local to global node number
    void Pde2MeshNode(StdVector<UInt> & meshNodes,
                      const StdVector<UInt> & pdeNodes) const;
    
    //! Map local to global node number
    UInt Pde2MeshNode(const UInt pdeNode) const;

    //! Map global to local elem number
    UInt Mesh2PdeElem(const UInt elemNumGlob) const;

    //! Map local to global elem number
    UInt Pde2MeshElem(const UInt elemNumLoc) const;

    //@}

    // ======================================================================
    // MISCELLANEOUS
    // ======================================================================
    //@{
    //! \name Miscellaneous

    //! Print the mapping nodes <-> Eqns
    void Print(std::ostream & out) const;

    //@}
    
  private:

    // ======================================================================
    // PRIVATE HELPER METHODS
    // ======================================================================

    //! Trigger local/global mapping of element/nodal numbers
    void CalcNodeElemMapping();

    //! Trigger local/global mapping of edge numbers
    void CalcEdgeMapping();

    //! Trigger local/global mapping of face numbers
    void CalcFaceMapping();

    //! Calculate equation numbers for nodes

    //! Triggers the mapping for nodal equations
    //! \param phase Parameter indicating if numbering of free equations (1)
    //!              has to be performed or if the fixed equations have to
    //!              be mapped (2)
    void CalcNodalEquations( UInt phase );

    //! Calculate equation numbers for elements interior (pfem)

    //! Triggers the mapping for element equations
    //! \param phase Parameter indicating if numbering of free equations (1)
    //!              has to be performed or if the fixed equations have to
    //!              be mapped (2)
    void CalcElemInteriorEquations( UInt phase );

    //! Calculate equation numbers for elements (constants)

    //! Triggers the mapping for element equations
    //! \param phase Parameter indicating if numbering of free equations (1)
    //!              has to be performed or if the fixed equations have to
    //!              be mapped (2)
    void CalcElemConstEquations( UInt phase );

    //! Calculate equation numbers for edges

    //! Triggers the mapping for edge equations
    //! \param phase Parameter indicating if numbering of free equations (1)
    //!              has to be performed or if the fixed equations have to
    //!              be mapped (2)
    void CalcEdgeEquations( UInt phase );

    //! Calculate equation numbers for faces

    //! Triggers the mapping for face equations
    //! \param phase Parameter indicating if numbering of free equations (1)
    //!              has to be performed or if the fixed equations have to
    //!              be mapped (2)
    void CalcFaceEquations( UInt phase );

    //! Extract node numbers from given entitylist
    void GetNodesOfEntities( StdVector<UInt>& nodeNr, 
                             shared_ptr<EntityList> ent );


    // ======================================================================
    // FLAGS
    // ======================================================================
    
    //! flag indicating initialization
    bool isFinalized_;

    //! Flag for turning on ordering of equation numbers

    //! If this flag is true, the equation numbers are ordered such that
    //! the set of highest equation numbers belongs to "unknowns" whose
    //! values are fixed by inhomogeneous Dirichlet boundary conditions.
    bool sortEqns_;
    
    // ======================================================================
    // EQUATION DATA
    // ======================================================================
    
    //! Pde Id of associated pde
    PdeIdType pdeId_;

    //! Number of equations
    UInt numEqns_;

    //! Number of "real" equations in Pde

    //! This attribute stores the number of 'real' equation numbers. By the
    //! latter we understand equation numbers for degrees of freedom that
    //! are not fixed by either
    //! - homogeneous Dirichlet boundary conditions
    //! - inhomogeneous Dirichlet boundary conditions
    //! - constraints (master - slave dof relations)
    UInt numRealEqns_;

    //! Number equations with inhomogeneous Dirichelt boundary condition
    UInt numIdBcs_;
    
    //! Store all results (before mapping )
    ResultEntityMap resEntMap_;

    //! Store all nodal-mapped element lists
    ResultEntityMap nodeMappedList_;

    //! Store all edge-mapped element list
    ResultEntityMap edgeMappedList_;
    
    //! Store all face-mapped entity lists
    ResultEntityMap faceMappedList_;

    //! Store all interio element-mapped lists
    ResultEntityMap elemIntMappedList_;

    //! Store all constant element-mapped entity lists
    ResultEntityMap elemConstMappedList_;
    
    //! Store all free entity lists
    ResultEntityMap freeMappedList_;
    
    //! Map of results and their nodal equations (local node -> eqn)
    EqnMapType nodeEqns_;
    
    //! Map of results and their nodal-edge equations (local edge -> eqn)
    VecEqnMapType edgeEqns_;

    //! Map of results and their nodal-face equations (local face -> eqn)
    VecEqnMapType faceEqns_;
    
    //! Map of results and their element equations (local element -> eqn)
    VecEqnMapType elemEqns_;

    //! Map of results and their free equations
    EqnMapType freeEqns_;
    
    // ======================================================================
    // GLOBAL/LOCAL MAPPING
    // ======================================================================

    //! Number of local nodes
    UInt numLocNodes_;

    //! Number of local elements
    UInt numLocElems_;
    
    //! Number of local edges
    UInt numLocEdges_;
 
    //! Number of local faces
    UInt numLocFaces_;

    //! Elements which are locally element/nodal/edge/face mapped
    StdVector<shared_ptr<ElemList> > locElems_;

    //! Mapping Local -> Global node numbering
    StdVector<UInt> pde2MeshNode_;

    //! Mapping Global -> Local node numbering
    StdVector<Integer> mesh2PdeNode_;
 
    //! Element mapping from local->global
    StdVector<UInt> pde2MeshElem_;
  
    //! Element mapping from global->local
    StdVector<Integer> mesh2PdeElem_;

    //! Edge mapping from global->local
    StdVector<Integer> mesh2PdeEdge_;

    //! Face mapping from global->local
    StdVector<Integer> mesh2PdeFace_;

    // ======================================================================
    // BOUNDARY CONDITIONS
    // ======================================================================

    //! List of homogeneous Dirichlet boundary conditions
    ResultHdBcMap hdBcs_;
    
    //! List of inhomogeneous Dirichlet boundary conditions
    ResultIdBcMap idBcs_;
    
    //! List of constraints
    ResultConstraintMap constraints_;
    
    //! ======================================================================
    //! POINTERS
    //! ======================================================================

    //! Pointer to grid
    Grid * ptGrid_;
     
  };
  
  
#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class EqnMap
  //! 
  //! \purpose
  //!
  //! This class is the base class which defines the interface of all
  //! objects that deal with the assignment of equation numbers to degrees
  //! of freedom. The following gives an overview on the meaning of the
  //! different equation numbers
  //! <center><img src="../AddDoc/splitting.png"></center>
  //!
  //! \collab 
  //!
  //! \implement 
  //!
  //! \status In use
  //!
  //! \unused 
  //!
  //! \improve
  //! 

#endif
}

#endif
