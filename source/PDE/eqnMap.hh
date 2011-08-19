// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_EQNMAP_HH
#define FILE_CFS_EQNMAP_HH


#include "Utils/StdVector.hh"
#include "Domain/resultInfo.hh"
#include "Domain/bcs.hh"
#include "MatVec/vector.hh"

namespace CoupledField {

  // Forward class declarations
  class EntityList;
  class NodeList;
  class ElemList;
  class EntityIterator;
  struct ResultInfo;
  class ParamNode; 
  
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
    typedef std::map<ResultInfo, IdFileBcList> ResultIdFileBcMap;
    typedef std::map<ResultInfo, IdFileBcList> ResultIdFiBcMap;
    typedef std::map<ResultInfo, ConstraintList> ResultConstraintMap;
    
    //! Constructor
    EqnMap( Grid* ptGrid, FeFctIdType, bool usePenalty);
    
    //! Destructor
    virtual ~EqnMap();
    
    // ======================================================================
    // SET AND INITIALIZATION METHODS
    // ======================================================================
    
    //@{
    //! \name Set And Initialization Methods
    
    //! Add a new result
    virtual void AddResult( ResultInfo& result, 
                    shared_ptr<EntityList> list );
    
    //! Set the homogeneous boundary conditions
    virtual void SetHomoDirichletBCs( HdBcList& hdBcs );
        
    //! Set the in-homogeneous boundary conditions
    virtual void SetInhomDirichletBCs( IdBcList& idBcs );
        
    //! Set the in-homogeneous boundary conditions
    virtual void SetInhomDirichFileBCs( IdFileBcList& idFiBcs );

    //! Set the constraint conditions
    virtual void SetConstraints( ConstraintList& constraints );
    
    //! Finalize setup and calculate mapping
    virtual void Finalize();

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
    virtual void ReorderMapping( const StdVector<UInt>& order );
    
    //@}
    // ======================================================================
    // GENERAL QUERY METHODS
    // ======================================================================

    //@{
    //! \name General Query Methods
    //! Return true if eqns are assigned
    virtual inline bool IsFinalized() const {return isFinalized_;}

    //! Return the total number of equations
    virtual inline UInt GetNumEqns() const { return numEqns_; }

    //! Return the equation number of the last unfixed degree of freedom

    //! This method returns the equation number of the last free degree of
    //! freedom. We consider a degree of freedom to be free, if it is not
    //! fixed by either
    //! - homogeneous Dirichlet boundary conditions
    //! - inhomogeneous Dirichlet boundary conditions
    //! - constraints (master - slave dof relations)
    //!
    //! The return value actually depends on the status of the usePenalty_ flag.
    //! If it is true,
    //! then inhomogeneous Dirichlet boundary conditions are treated by the
    //! penalty method and must be considered free dofs as well. Thus, in this
    //! case we return the total number of equations. Otherwise we return only
    //! the number of equations withou a dirichlet boundary condition
    virtual UInt GetNumLastFreeDof() const;
                              
    //! Return number of real inhomogeneous Dirichlet boundary conditions

    //! This method returns the number of 'real' boundary conditions, i.e.
    //! douplets of nodes are already removed, so this number represents
    //! the number of equations, which are not 'free'
    virtual inline UInt GetNumInHomDirichletEqns() const {
      return numIdBcs_;
    }
    /** Return number of real inhomogeneous Dirichlet boundary conditions which
     * are read from a file
     * This method returns the number of 'real' boundary conditions, i.e.
     * douplets of nodes are already removed, so this number represents
     * the number of equations, which are not 'free'
     * @return the number of inhomogeneous Dirichlet boundary conditions which
     * are read from a file
     */
    virtual inline UInt GetNumInHomDirichletFileEqns() const {
      return numIdFiBcs_;
    }

    //! Return number of constraint slave equations.

    //! This method returns the number of equations, which are fixed
    //! by a slave constraint condition.
    virtual inline UInt GetNumConstraintSlaveEqns() const {
      return numCs_;
    }
                              
    //! Return id of associated PDE
    virtual FeFctIdType GetPdeId() const {
      return pdeId_;
    }
    
    //@}

    // ======================================================================
    // EQUATION MAPPING
    // ======================================================================
    
    //@{
    //! \name Equation Mapping
    //! Map result and entities to equation numbers
    virtual void GetEqns( StdVector<Integer>& eqns,
                  const ResultInfo& result, const EntityIterator& it ) const;

    //! Map result, entities and dof numbers to equation numbers
    virtual void GetEqns( StdVector<Integer>& eqns,
                  const ResultInfo& result, const EntityIterator& it,
                  UInt dof ) const;
    
    //! Map combination of result, entity and dof to single equation number
    virtual Integer GetEqn( const ResultInfo& result, 
                    const EntityIterator& it, UInt dof ) const;

    //! Name mapping function for obtaining a nodal equation
    virtual Integer GetNodeEqn( const ResultInfo& result, UInt nodeNr, UInt dof );

    /** Name mapping function for obtaining a nodal equation
     @param dof 1 based! */
    virtual Integer GetNodeEqn( UInt nodeNr, UInt dof );

    /** An expensive reverse of GetNodeEqn().
     * So expensive that it should be only used for debugging */
    Integer SearchNode(Integer eqn, UInt dof);

    virtual void GetNodeEqn( const StdVector<UInt>& nodeNrs, 
                     StdVector<Integer>& eqnNrs );

    //@}
    
    // ======================================================================
    // LOCAL/GLOBAL MAPPING OF MESH ENTITIES
    // ======================================================================
    
    //@{ 
    //! \name Local/Global Mapping of Mesh Entities

    //! Get number of local nodes
    virtual inline UInt GetNumLocalNodes() const {return numLocNodes_;}
    
    //! Get number of local elements
    virtual inline UInt GetNumLocalElems() const {return numLocElems_;}

    //! Get number of local edges
    virtual inline UInt GetNumLocalEdges() const {return numLocEdges_;}

    //! Get number of local faces
    virtual inline UInt GetNumLocalFaces() const {return numLocFaces_;}

    //! Map global to local node numbers
    //! (needed for nodal displacement of grid)
    virtual void Mesh2PdeNode(StdVector<UInt> & PdeNodes,
                      const StdVector<UInt> & MeshNodes) const;

    //! Map global to local node number
    //! (needed for nodal displacement of grid)
    virtual UInt Mesh2PdeNode(const UInt meshNode) const;
    
    //! Map local to global node number
    virtual void Pde2MeshNode(StdVector<UInt> & meshNodes,
                      const StdVector<UInt> & pdeNodes) const;
    
    //! Map local to global node number
    virtual UInt Pde2MeshNode(const UInt pdeNode) const;

    //! Map global to local elem number
    virtual UInt Mesh2PdeElem(const UInt elemNumGlob) const;

    //! Map local to global elem number
    virtual UInt Pde2MeshElem(const UInt elemNumLoc) const;

    //@}

    /** Give all details of the mapping to the info.xml file it triggered 
     * such in the comman line. One can gain similar data via the loggin 
     * (in the debug version) */ 
    virtual void ToInfo(PtrParamNode in) const;

    //! Reads the comst significant infos from the given Map and copys them
    //! used e.g. in the mixed equation map to determine the number of unknowns already assigned
    virtual void CopyMapInfo(EqnMap* map);
    
  protected:

    // ======================================================================
    // PRIVATE HELPER METHODS
    // ======================================================================

    //! Trigger local/global mapping of element/nodal numbers
    virtual void CalcNodeElemMapping();

    //! Trigger local/global mapping of edge numbers
    virtual void CalcEdgeMapping();

    //! Trigger local/global mapping of face numbers
    virtual void CalcFaceMapping();

    //! Calculate equation numbers for nodes

    //! Triggers the mapping for nodal equations
    //! \param phase Parameter indicating if numbering of free equations (1)
    //!              has to be performed or if the fixed equations have to
    //!              be mapped (2)
    virtual void CalcNodalEquations( UInt phase );

    //! Calculate equation numbers for elements interior (pfem)

    //! Triggers the mapping for element equations
    //! \param phase Parameter indicating if numbering of free equations (1)
    //!              has to be performed or if the fixed equations have to
    //!              be mapped (2)
    virtual void CalcElemInteriorEquations( UInt phase );

    //! Calculate equation numbers for elements (constants)

    //! Triggers the mapping for element equations
    //! \param phase Parameter indicating if numbering of free equations (1)
    //!              has to be performed or if the fixed equations have to
    //!              be mapped (2)
    virtual void CalcElemConstEquations( UInt phase );

    //! Calculate equation numbers for edges

    //! Triggers the mapping for edge equations
    //! \param phase Parameter indicating if numbering of free equations (1)
    //!              has to be performed or if the fixed equations have to
    //!              be mapped (2)
    virtual void CalcEdgeEquations( UInt phase );

    //! Calculate equation numbers for faces

    //! Triggers the mapping for face equations
    //! \param phase Parameter indicating if numbering of free equations (1)
    //!              has to be performed or if the fixed equations have to
    //!              be mapped (2)
    virtual void CalcFaceEquations( UInt phase );

    //! Extract node numbers from given entitylist
    virtual void GetNodesOfEntities( StdVector<UInt>& nodeNr, 
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
    bool usePenalty_;
    
    // ======================================================================
    // EQUATION DATA
    // ======================================================================
    
    //! Pde Id of associated pde
    FeFctIdType pdeId_;

    //! Number of equations
    UInt numEqns_;

    //! Number equations with inhomogeneous Dirichelt boundary condition
    UInt numIdBcs_;

    //! Number equations with inhomogeneous Dirichelt boundary condition
    UInt numIdFiBcs_;

    //! Number equations with slave constraint condition
    UInt numCs_;
        
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

    //! Entities which are locally element/nodal/edge/face mapped
    StdVector<shared_ptr<EntityList> > locEntities_;

    //! Mapping Local -> Global node numbering
    StdVector<UInt> pde2MeshNode_;

    //! Element mapping from local->global
    StdVector<UInt> pde2MeshElem_;
  
    //! Element mapping from global->local
    StdVector<Integer> mesh2PdeElem_;

    // ======================================================================
    // BOUNDARY CONDITIONS
    // ======================================================================

    //! List of homogeneous Dirichlet boundary conditions
    ResultHdBcMap hdBcs_;
    
    //! List of inhomogeneous Dirichlet boundary conditions
    ResultIdBcMap idBcs_;
    
    //! List of inhomogeneous Dirichlet boundary conditions
    ResultIdFileBcMap idFiBcs_;
    
    //! List of constraints
    ResultConstraintMap constraints_;
    
    //! ======================================================================
    //! POINTERS
    //! ======================================================================

    //! Pointer to grid
    Grid * ptGrid_;
     

private:
    //! Edge mapping from global->local
    StdVector<Integer> mesh2PdeEdge_;

    //! Face mapping from global->local
    StdVector<Integer> mesh2PdeFace_;

    //! Mapping Global -> Local node numbering
    StdVector<Integer> mesh2PdeNode_;
 
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
