// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef GRAPH_MANAGER_HH
#define GRAPH_MANAGER_HH

#include <string>
#include "General/environment.hh"
#include "OLAS/graph/baseordering.hh"

namespace CoupledField {


// forward class declarations
class BaseGraph;
class IDBC_Graph;

  //! This is the base class for all graph manager classes. A graph manager is
  //! an object that is repsonsible for administering the graphs representing
  //! the connectivity information related to either a single feFcuntion or e.g. in
  //! the case of direct coupled field simulations that of a collection of
  //! feFunctions.
  //! The graph manager functions as an intermediary between the graph
  //! object(s) and the feFunctions.
  //!
  //! The following piece of pseudo-code demonstrates the sequence of method
  //! calls that must be observed for setting up the graphs:
  //!
  //! \code
  //! /* The graph manager must be told how many PDEs lie in its
  //!    responsibility */
  //! SetupInit( numPDEs );
  //!
  //! /* Now each PDE must register itself with the graph manager using
  //!    a unique identifier */
  //! for ( i = 1; i <= numPDEs; i++ ) {
  //!    RegisterPDE( ... );
  //! }
  //!
  //! /* One PDE after the other can now pass its connectivity information to
  //!    the graph manager */  
  //! for ( i = 1; i <= numPDEs; i++ ) {
  //!
  //!    /* Initialise assembly of a new PDE sub-graph */
  //!    AssembleInit( ... );
  //!
  //!    /* Set information for each Finite Element */
  //!    for ( k = 1; k <= numElemInThisPDE; k++ ) {
  //!    SetElementPos( ... );
  //!    }
  //!
  //!    /* Finish assembly of the PDE sub-graph */
  //!    AssembleDone( ... );
  //! }
  //!
  //! /* After the PDEs the coupling objects (if any) pass their connectivity
  //!    information to the graph manager */  
  //! for ( j = 1; j <= numCouplings; j++ ) {
  //!
  //!    /* Initialise assembly of a new coupling sub-graph */
  //!    AssembleInit( ... );
  //!
  //!    /* Set information for each Finite Element belonging to both PDEs
  //!       that are coupled by the coupling object */
  //!    for ( k = 1; k <= numElemInThisCoupling; k++ ) {
  //!    SetElementPos( ... );
  //!    }
  //!
  //!    /* Finish assembly of the coupling sub-graph */
  //!    AssembleDone( ... );
  //! }
  //!
  //! /* By calling this method we inform the graph manager that all coupling
  //!    objects have passed their connectivity and that the setup phase now
  //!    is completed */
  //! SetupDone();
  //! \endcode
  class GraphManager {

  public:

    
    //! Struct containing information for SBM-system

    //! Lightweight struct for storing information per SBM Matrix block
    typedef struct {

      //! Map (fctId,eqnNr) to matrix / vector index
      StdVector<std::map<UInt, UInt> > eqnToIndex;

      //! Define mappings for subBlocks
      //! 1st index: blockIndex
      //! pair: [begin,end] of one index block
      StdVector<std::pair<UInt,UInt> > indexBlocks;

      //! Total number of equations in this block
      UInt size;
      
      //! Number of last free index in this block
      UInt numLastFreeIndex;

      //! Flag if minor blocks are defined

      //! If the SBM matrix has minor blocks, no reordering can
      //! be applied
      bool hasSubBlocks;

    } SBMBlockInfo;
    
    // =======================================================================
    // SETUP
    // =======================================================================

    //@{ \name Methods for construction, destruction and setup

    //! Default Constructor
    GraphManager();

    //! Destructor
    virtual ~GraphManager();

    //! Prepares graph manager for generation of sub-graphs

    //! This method must be called after construction of the graph manager
    //! object and before any other method. It allows the graph manager to
    //! prepare itself for the registration and generation of the sub-graphs.
    //! \param numBlocks number of SBM blocks whose connectivity graphs the
    //!                          graph manager will administer.
    //! \para, useDistinctGraphs if true, every matrixType 
    void SetupInit( UInt numBlocks, bool useDistinctGraphs );

    //! Finalises setup of the graph manager

    //! This method must be called after the assembly of all sub-graphs was
    //! done, i.e. once all blocks have conveyed their
    //! connectivity information to the graph manager.
    void SetupDone();

    //! Register block with graph manager
    
    //! Before a matrix block can set its sparsity pattern, the matrix
    //! block has to be registered and the number of unknowns and 
    //! non free equations has to be set.
    //! \param blockNum       block number of matrix to be defined 
    //! \param numEqns        number of unknowns in the linear system
    //!                       associated with the PDE (i.e. number of
    //!                       equation numbers from CFS++)
    //! \param numLastFreeEqn specifies the equation number of the last
    //!                       unfixed degree of freedom (if inhomogeneous
    //!                       Dirichlet boundary conditions are treated by
    //!                       the penalty approach numLastFreeDof
    //!                       \f$\stackrel{!}{=}\f$ numEqns)
    //! \param reorder        Re-ordering type for this PDE (optional)
    void RegisterBlock( const UInt blockNum,
                        SBMBlockInfo* blockInfo,
                        const BaseOrdering::ReorderingType reorder =
                        BaseOrdering::NOREORDERING );

    //@}

    // =======================================================================
    // SETUP OF SUB_GRAPH
    // =======================================================================

    //@{ \name Methods for generation of sub-graphs
    
    //! Insert connectivity of a finite element into the matrix graph

    //! This method is used to insert the connectivity of the unknowns
    //! attached to a finite element into the sub-graph representing the
    //! structure of the matrix associated to a single PDE or the coupling
    //! between two PDEs. In the case that the sub-graph belongs to a PDE only
    //! the corresponding identifier (obtained from RegisterPDE) needs to be
    //! supplied. In the case that the graph describes the coupling between
    //! two PDEs both related identifiers must be specified.
    //! \param identifierPDE1 identifier for first PDE related to sub-graph
    //! \param connect1       array containing equation numbers for the
    //!                       degrees of freedom attachted to the element
    //!                       (unknowns at nodes or edges, whatsoever), must
    //!                       be a one based pointer; for first PDE
    //! \param elemSize1      length of the connect array; for first PDE
    //! \param identifierPDE2 identifier for first PDE related to sub-graph
    //! \param connect2       array containing equation numbers for the
    //!                       degrees of freedom attachted to the element
    //!                       (unknowns at nodes or edges, whatsoever), must
    //!                       be a one based pointer; for second PDE
    //! \param elemSize2      length of the connect array; for second PDE
    //! \param setCounterPart flag signalling whether the counterpart, i.e.
    //!                       the connectivity with the PDE identifiers and
    //!                       equation numbers reversed, should also be
    //!                       inserted into the graph
    virtual void SetElementPos( const StdVector<UInt>& rowBlock,
                                const StdVector<UInt>& rowIndices,
                                const StdVector<UInt>& colBlock,
                                const StdVector<UInt>& colIndices,
                                FEMatrixType matrixType,
                                bool setCounterPart );
    //@}

    // =======================================================================
    // QUERY METHODS
    // =======================================================================

    //@{ \name Methods for querying information

    //! Get a specified (sub-)graph

    //! After the assembly of the graph is completed, the specified
    //! (sub-)graph can be accessed to communicate directly with it.
    //! This is needed for example to create a sparse matrix, where the
    //! matrix graph has to be supported.
    //! \param identifierPDE1 identifier for first PDE related to sub-graph
    //! \param identifierPDE2 identifier for second PDE related to sub-graph
    //! \return a BaseGraph object according to the specified identifier
    //! combination
    //! \note The method refuses to return a pointer to a non-existant
    //!       graph and will instead report an error
    BaseGraph* GetGraph( UInt rowNum = 0,
                         UInt colNum = 0 );

    //! Get a specified (sub-)graph for inhom. Dirichlet BC

    //! After the assembly of the graphs storing the connectivity between
    //! real dofs and dofs fixed by inhomogeneous Dirichlet boundary
    //! conditions is completed, the specified IDBC (sub-)graph for a certain
    //! PDE pair can be accessed to communicate directly with it.
    //! \param pdeID1 identifier for first PDE of pair
    //! \param pdeID2 identifier for second PDE of pair
    //! \return a BaseGraph object according to the specified identifier
    BaseGraph* GetIDBCGraph( UInt rowNum = 0,
                             UInt colNum = 0 ) const;

    //! Obtain the permutation/re-ordering vector for a graph

    //! The PDE attached to the GraphManagerSBMMat can call this method to
    //! obtain the re-ordering vector for its unknowns, once assembly of the
    //! graph was completed by a call to AssembleDone().
    //! \param identifier unique identifier for the PDE registered with the
    //!                   graph manager
    //! \return an integer array containing the new equation numbers after
    //!         the re-ordering or a NULL pointer, if no re-ordering was
    //!         performed.
    //! \note While memory for the permutation vector is allocated in this
    //!       method, it is the caller's responsibility to dispose of that
    //!       memory once it no longer needs the array.
    void GetReordering( UInt blockNum, StdVector<UInt>& order ) ;

    //! Print some statistics on the graph manager
    void PrintStats();

    //! Query whether a certain sub-graph exists

    //! This method can be used to query whether a sub-graph with the index
    //! pair (idPDE1, idPDE2) was generated in the Setup phase.
    //! \return true if graph exists, false otherwise
    bool SubGraphExists( UInt rowNum = 0,
                         UInt colNum = 0 )  {
      return graph_[ ComputeIndex( rowNum, colNum ) ] != NULL;
    }

    //! Query whether a certain IDBC sub-graph exists

    //! This method can be used to query whether a sub-graph of IDBC type
    //! was generated for the specified PDE in the Setup phase.
    //! \return true if IDBC graph exists, false otherwise
    bool IDBCGraphExists( UInt rowNum = 0,
                          UInt colNum = 0 ) const {
      return graphIDBC_[ ComputeIndex( rowNum, colNum ) ] != NULL;
    }

    //@}


  private:

    //! Computes the index into the graph pointer matrix

    //! The class stores the pointers to the graphs associated with the PDEs
    //! and the coupling objects in a matrix. This matrix is stored as 1D
    //! array. This auxilliary method computes for a given pair of PDE
    //! identifiers the 1D index of the corresponding matrix entry.
    inline UInt ComputeIndex( UInt rowNum,
                              UInt colNum ) const {
      return numBlocks_*rowNum + colNum;
    };

    //! Auxilliary method for generation of IDBC graph objects

    //! This method will check for a pair of PDEs (pdeID1, pdeID2) whether
    //! it is necessary to generate of an associated IDBC graph. This will
    //! be the case, if the second PDE of the pair contains degrees of
    //! freedom that are fixed by inhomogeneous Dirichlet boundary
    //! conditions.
    void GenerateIDBCGraph( UInt rowNum, UInt colNum );

    //! Auxilliary method for generation of graphs of coupling objects

    //! This method can be called to generate the graph of a coupling object
    //! for the coupling between the PDEs in the specified pair
    //! (pdeId1, pdeID2). It is called by AssembleInit().
    void GenerateCouplingGraph( UInt rowNum, UInt colNum );

    //! Matrix storing pointers to the graph objects

    //! We use a flattened 2D structure to store for each PDE and coupling
    //! object a pointer to the associated graph object. In analogy to the
    //! resulting %SBM_Matrix we use a 2D matrix for storing the pointers.
    //! The graph belonging to the object with identifiers id1 and id2 is
    //! stored as matrix entry (id1,id2). We use C-style mapping to transform
    //! the matrix to a 1D storage layout. Entries in the matrix may be NULL
    //! pointers, if no corresponding coupling object exists.
    StdVector<BaseGraph *> graph_;

    //! Matrix storing pointers to the IDBC graph objects

    //! This matrix stores for each PDE pair a pointer to an associated graph
    //! for handling the connection between free degrees of freedom and
    //! those fixed by inhomogeneous Dirichlet boundary conditions.
    StdVector<IDBC_Graph *> graphIDBC_;

    //! Matrix storing coupling flags for PDE pairs

    //! We use a flattened 2D structure to store for each PDE pair
    //! (pde1,pde2) a flag that indicates whether there is a coupling object
    //! for this pair. In analogy to the resulting %SBM_Matrix we use a 2D
    //! matrix for storing these flags. The flag for PDE pair (pde1,pde2) is
    //! stored as matrix entry (pde1,pde2). We use C-style mapping to
    //! transform the matrix to a 1D storage layout. See the ComputeIndex()
    //! method.
    StdVector<bool> isCoupled_;

    //! Array containing pointers to permutation vectors

    //! We allocate memory for the permutation vector of each PDE
    //! in this class, but do not assume that we are responsible for their
    //! de-allocation. This is the task of the respective PDE. However,
    //! if the PDE does not claim its permutation vector by calling
    //! GetReordering, we will de-allocate the memory in the destructor
    //! of this class. An entry of this array may be a NULL pointer in one
    //! of two cases
    //! - no reordering is performed for the PDE
    //! - the permutation vector was already claimed via GetReordering()
    StdVector< StdVector<UInt> > newOrdering_;

    //! Vector storing for each block meta informaation
    StdVector<SBMBlockInfo*> blockInfo_;

    //! Attribute to store the number of PDEs that belong to the manager
    UInt numBlocks_;

    //! Attribute to store the number of PDEs that finished assembly

    //! We use this attribute to keep track of the number of PDEs that have
    //! already assembled their sub-graphs. This is useful for some
    //! consistency checks. We can e.g. avoid that coupling objects try
    //! to start their assembly before all PDEs have finished theirs, which
    //! could lead to problems due to the re-ordering of the PDEs.
    UInt numAssembledBlocks_;

    //! Attribute to keep track of number of PDEs that were registered
    UInt numRegisteredBlocks_;

    //! Attribute for testing whether it is safe to query the graph re-ordering

    //! Since there is a fundamental difference between the different graph
    //! manager implementations w.r.t. to the place at which the graph manager
    //! instructs the graph(s) to re-order themselves, we use this attribute
    //! for checking in GetReordering() whether it is safe to return the
    //! computed re-ordering.
    //! \note In the case of the %GraphManagerSBMMat class the re-ordering of
    //! the graph associated with each PDE triggered once the assembly of
    //! that PDE is finalised by calling AssembleDone().
    bool reorderingDone_;

    //! Flag indicating whether registration of PDEs and Couplings is done

    //! This attribute is used as flag which indicates whether the
    //! registration of PDEs and coupling objects is finished. Since there
    //! is currently no method which could be called in this case, we
    //! simply assume that the registration phase is over when AssembleInit()
    //! is called for the first time and set the flag to true then.
    bool registrationDone_;

    //@{
    //! Auxilliary vector for passing connectivities to the graph objects

    //! This auxilliary vector is employed by the SetElementPos() method to
    //! pass the connectivity information that was transformed by a call
    //! to the AdaptConnects() method to the graph objects.
    std::vector<std::vector<UInt> > vertexList1_;
    std::vector<std::vector<UInt> > vertexList2_;
    std::vector<std::vector<UInt> > edgeList1_;
    std::vector<std::vector<UInt> > edgeList2_;
    //@}

  };

}

#endif
