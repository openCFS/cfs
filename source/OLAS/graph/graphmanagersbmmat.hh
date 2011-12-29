// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef GRAPH_MANAGER_SBMMAT
#define GRAPH_MANAGER_SBMMAT

#include <string>
#include <vector>

#include "General/defs.hh"
#include "General/environment.hh"
#include "OLAS/graph/basegraphmanager.hh"
#include "OLAS/graph/baseordering.hh"
#include "Utils/StdVector.hh"

namespace CoupledField {

  //! Graph manager for the case of several PDEs and a SBM_Matrix

  //! This class implements the graph manager concept for the case that
  //! several PDEs are directly coupled and the corresponding matrices are
  //! to be assembled into a super-block matrix (SBM_Matrix) structure.
class BaseGraph;
class IDBC_Graph;

  class GraphManagerSBMMat : public BaseGraphManager {

  public:

    // =======================================================================
    // SETUP
    // =======================================================================

    //@{ \name Methods for construction, destruction and setup

    //! Default Constructor
    GraphManagerSBMMat();

    //! Destructor
    virtual ~GraphManagerSBMMat();

    //! Prepares graph manager for generation of sub-graphs

    //! This method must be called after construction of the graph manager
    //! object and before any other method. It allows the graph manager to
    //! prepare itself for the registration and generation of the sub-graphs.
    //! In the case of this class ...
    //! \param numPDEs number of PDE objects whose connectivity graphs the
    //!                          graph manager will administer.
    void SetupInit( UInt numPDEs );

    //! Finalises setup of the graph manager

    //! This method must be called after the assembly of all sub-graphs was
    //! done, i.e. once all PDEs a coupling objects have conveyed their
    //! connectivity information to the graph manager.
    //! \note In the case of the %GraphManagerStdMat case this method triggers
    //! the re-ordering of the unique graph. As a consequence the permutation
    //! vector resulting from the re-ordering may only be queried after
    //! SetupDone() was called.
    void SetupDone();

    //! Register PDE with graph manager
    
    //! A PDE must call this method to register itself with the Graph manager
    //! and to give its number of unkonwns. In the case of CFS++ this method 
    //! is not called  directly, but via the algebraic system interface.
    //! \param identifierPDE  identifier of PDE related to sub-graph
    //! \param numEqns        number of unknowns in the linear system
    //!                       associated with the PDE (i.e. number of
    //!                       equation numbers from CFS++)
    //! \param numLastFreeDof specifies the equation number of the last
    //!                       unfixed degree of freedom (if inhomogeneous
    //!                       Dirichlet boundary conditions are treated by
    //!                       the penalty approach numLastFreeDof
    //!                       \f$\stackrel{!}{=}\f$ numEqns)
    //! \param reorder        Re-ordering type for this PDE (optional)
    void RegisterPDE( const PdeIdType identifierPDE,
                      const UInt numEqns,
                      const UInt numLastFreeDof,
                      const BaseOrdering::ReorderingType reorder =
                      BaseOrdering::NOREORDERING );

    //@}

    // =======================================================================
    // SETUP OF SUB_GRAPH
    // =======================================================================

    //@{ \name Methods for generation of sub-graphs
    
    //! Prepare graph manager for assembling a sub-graph

    //! Before the assembly phase of a sub-graph associated with a PDE or the
    //! coupling between two PDEs can be started this method must be called
    //! in order for the graph manager to do administrative stuff. In the case
    //! that the sub-graph belongs to a PDE only the corresponding identifier
    //! (obtained from RegisterPDE) needs to be supplied. In the case that the
    //! graph describes the coupling between two PDEs both related identifiers
    //! must be specified.
    //! \param identifierPDE1 identifier for first PDE related to sub-graph
    //! \param identifierPDE2 identifier for second PDE related to sub-graph
    //! \param assemblingTranspose indicates whether the graph of the
    //!                            transpose coupling object will be assembled
    //!                            together with the coupling object.
    void AssembleInit( const PdeIdType identifierPDE1,
                       const PdeIdType identifierPDE2,
                       bool assemblingTranspose );

    //! Finalise assembly of a sub-graph

    //! Once the assembly phase of a sub-graph associated with a PDE or the
    //! coupling between two PDEs is done, this method must be called in order
    //! to inform the graph about this. The call allows the graph manager to
    //! immediately trigger associated tasks, like e.g. graph conversion or
    //! re-ordering. In the case that the sub-graph belongs to a PDE only the
    //! corresponding identifier (obtained from RegisterPDE) needs to be
    //! supplied. In the case that the graph describes the coupling between
    //! two PDEs both related identifiers must be specified.
    //! \param identifierPDE1 identifier for first PDE related to sub-graph
    //! \param identifierPDE2 identifier for second PDE related to sub-graph
    //! \param assemblingTranspose indicates whether the graph of the
    //!                            transpose coupling object was assembled
    //!                            together with the coupling object.
    void AssembleDone( const PdeIdType identifierPDE1,
                       const PdeIdType identifierPDE2,
                       bool assemblingTranspose );

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
    virtual void SetElementPos( const PdeIdType identifierPDE1,
                                const StdVector<Integer>& eqnNrs1,
                                const PdeIdType identifierPDE2,
                                const StdVector<Integer>& eqnNrs2,
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
    BaseGraph* GetGraph( const PdeIdType identifierPDE1 = NO_PDE_ID,
                         const PdeIdType identifierPDE2 = NO_PDE_ID );

    //! Get a specified (sub-)graph for inhom. Dirichlet BC

    //! After the assembly of the graphs storing the connectivity between
    //! real dofs and dofs fixed by inhomogeneous Dirichlet boundary
    //! conditions is completed, the specified IDBC (sub-)graph for a certain
    //! PDE pair can be accessed to communicate directly with it.
    //! \param pdeID1 identifier for first PDE of pair
    //! \param pdeID2 identifier for second PDE of pair
    //! \return a BaseGraph object according to the specified identifier
    BaseGraph* GetIDBCGraph( const PdeIdType pdeID1,
                             const PdeIdType pdeID2 ) const;

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
    void GetReordering( const PdeIdType identifier, StdVector<UInt>& order ) ;

    //! Print some statistics on the graph manager
    void PrintStats( std::ostream *log );

    //! Query whether a certain sub-graph exists

    //! This method can be used to query whether a sub-graph with the index
    //! pair (idPDE1, idPDE2) was generated in the Setup phase.
    //! \return true if graph exists, false otherwise
    bool SubGraphExists( const PdeIdType idPDE1,
                         const PdeIdType idPDE2 ) const;

    //! Query whether a certain IDBC sub-graph exists

    //! This method can be used to query whether a sub-graph of IDBC type
    //! was generated for the specified PDE in the Setup phase.
    //! \return true if IDBC graph exists, false otherwise
    bool IDBCGraphExists( const PdeIdType idPDE1,
                          const PdeIdType idPDE2 ) const;

    //@}


  private:

    //! Computes the index into the graph pointer matrix

    //! The class stores the pointers to the graphs associated with the PDEs
    //! and the coupling objects in a matrix. This matrix is stored as 1D
    //! array. This auxilliary method computes for a given pair of PDE
    //! identifiers the 1D index of the corresponding matrix entry.
    inline UInt ComputeIndex( const PdeIdType identifierPDE1,
                              const PdeIdType identifierPDE2 ) const {
      UInt i = identifierPDE1;
      UInt j = identifierPDE2 == NO_PDE_ID ? identifierPDE1 : identifierPDE2;
      return numPDEs_*(i-1) + j;
    };

    //! Auxilliary method for generation of IDBC graph objects

    //! This method will check for a pair of PDEs (pdeID1, pdeID2) whether
    //! it is necessary to generate of an associated IDBC graph. This will
    //! be the case, if the second PDE of the pair contains degrees of
    //! freedom that are fixed by inhomogeneous Dirichlet boundary
    //! conditions.
    void GenerateIDBCGraph( PdeIdType pdeID1, PdeIdType pdeID2 );

    //! Auxilliary method for generation of graphs of coupling objects

    //! This method can be called to generate the graph of a coupling object
    //! for the coupling between the PDEs in the specified pair
    //! (pdeId1, pdeID2). It is called by AssembleInit().
    void GenerateCouplingGraph( PdeIdType pdeID1, PdeIdType pdeID2 );

    //! Matrix storing pointers to the graph objects

    //! We use a flattened 2D structure to store for each PDE and coupling
    //! object a pointer to the associated graph object. In analogy to the
    //! resulting %SBM_Matrix we use a 2D matrix for storing the pointers.
    //! The graph belonging to the object with identifiers id1 and id2 is
    //! stored as matrix entry (id1,id2). We use C-style mapping to transform
    //! the matrix to a 1D storage layout. Entries in the matrix may be NULL
    //! pointers, if no corresponding coupling object exists.
    BaseGraph **graph_;

    //! Matrix storing pointers to the IDBC graph objects

    //! This matrix stores for each PDE pair a pointer to an associated graph
    //! for handling the connection between free degrees of freedom and
    //! those fixed by inhomogeneous Dirichlet boundary conditions.
    IDBC_Graph **graphIDBC_;

    //! Matrix storing coupling flags for PDE pairs

    //! We use a flattened 2D structure to store for each PDE pair
    //! (pde1,pde2) a flag that indicates whether there is a coupling object
    //! for this pair. In analogy to the resulting %SBM_Matrix we use a 2D
    //! matrix for storing these flags. The flag for PDE pair (pde1,pde2) is
    //! stored as matrix entry (pde1,pde2). We use C-style mapping to
    //! transform the matrix to a 1D storage layout. See the ComputeIndex()
    //! method.
    bool *isCoupled_;

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

    //! Vector storing for each PDE the last free equation number

    //! Vector storing for each PDE the last free equation number, i.e. the
    //! last equation number for a degree of freedom that is not fixed by
    //! an inhomogeneous Dirichlet boundary condition.
    UInt *numLastFreeDof_;

    //! Vector storing for each PDE the total number of equation numbers

    //! Vector storing for each PDE the total number of equation numbers,
    //! i.e. the sum of the equation numbers for degrees of freedom that
    //! are fixed by inhomogeneous Dirichlet boundary conditions and those
    //! that are not.
    UInt *numEqn_;

    //! Attribute to store the number of PDEs that belong to the manager
    UInt numPDEs_;

    //! Attribute to store the number of PDEs that finished assembly

    //! We use this attribute to keep track of the number of PDEs that have
    //! already assembled their sub-graphs. This is useful for some
    //! consistency checks. We can e.g. avoid that coupling objects try
    //! to start their assembly before all PDEs have finished theirs, which
    //! could lead to problems due to the re-ordering of the PDEs.
    UInt numAssembledPDEs_;

    //! Attribute to keep track of number of PDEs that were registered
    UInt numRegisteredPDEs_;

    //! Attribute for testing wether it is safe to query the graph re-ordering

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
    std::vector<UInt> vertexList1_;
    std::vector<UInt> vertexList2_;
    std::vector<UInt> edgeList1_;
    std::vector<UInt> edgeList2_;
    //@}

  };

}

#endif
