#ifndef GRAPH_MANAGER_STDMAT
#define GRAPH_MANAGER_STDMAT

#include <map>
#include <string>
#include <vector>
#include <iostream>

#include "graph/basegraph.hh"
#include "graph/idbcgraph.hh"
#include "graph/basegraphmanager.hh"

namespace OLAS {

  //! Graph manager for the case of several PDEs and a StdMatrix

  //! This class implements the graph manager concept for the case that
  //! several PDEs are directly coupled and the corresponding matrices are
  //! to be assembled into a single (StdMatrix) matrix structure.
  //! \note In the case that all PDEs are assembled into a single matrix,
  //!       the following must be kept in mind. If the user wants to re-order
  //!       the matrix, new equation numbers are computed for each PDE. This
  //!       is handled by the graph object itself that delegates the work to
  //!       a re-ordering object. In the case that the user does not want to
  //!       perform a re-ordering, however, the class still must provide a
  //!       sort of re-ordering for the following reason. All PDEs, but the
  //!       one that was registered first, require on the CFS++ level equation
  //!       numbers that include an offset, so that the equation number points
  //!       to the respective sub-matrix. The class takes care of this, but
  //!       for efficiency reasons does not compute an identity permutation
  //!       for the first registered PDE, unless the elimination approach
  //!       for IDBC handling is used and there are fixed dofs in that PDE.
  //!       In the former case GetReordering() will return a NULL pointer.
  //! \todo Change the concept of the reordering another time. If no real
  //!       reordering is to be performed, it would be more efficient on the
  //!       CFS++ side, if we would pass back the offset via RegisterPDE() and
  //!       the equation data object just adds the offset, if the permutation
  //!       vector obtained from GetReordering() equals NULL.
  class GraphManagerStdMat : public BaseGraphManager {

  public:

    // =======================================================================
    // SETUP
    // =======================================================================

    //@{ \name Methods for construction, destruction and setup

    //! Default constructor
    GraphManagerStdMat();

    //! Destructor
    virtual ~GraphManagerStdMat();

    //! Prepares graph manager for generation of sub-graphs

    //! This method must be called after construction of the graph manager
    //! object and before any other method. It allows the graph manager to
    //! prepare itself for the registration and generation of the sub-graphs.
    //! In the case of this class this only includes minor tasks like, e.g.
    //! memory allocation for the vector storing the offsets and so on.
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
		      const ReorderingType reorder 
		      = NOREORDERING );

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
    //! \note Currently this method does not do anything in the special
    //!       case of this class, but must be provided for consistency.
    //! \param identifierPDE1 identifier for first PDE related to sub-graph
    //! \param identifierPDE2 identifier for second PDE related to sub-graph
    //! \param assemblingTranspose indicates whether the graph of the
    //!                            transpose coupling object was assembled
    //!                            together with the coupling object.
    void AssembleDone( const PdeIdType identifierPDE1,
		       const PdeIdType identifierPDE2,
                       bool assemblingTranspose ) {
      ENTER_FCN( "GraphManagerStdMat::AssembleDone" );
    }

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
    //! \note
    //!     - For performance reasons the offsets that are added to the
    //!       entries of the connect array are computed by SetupInit(). Thus,
    //!       it is crucial to not inter-mix the assembly of different
    //!       sub-graphs. The latter will only be detected in debug mode,
    //!       since only in this mode the identifiers passed to this method
    //!       are validated against those provided to SetupInit().
    //!     - In the case that the sub-graph belongs to a PDE and not to
    //!       a coupling object, it is safe to pass identical identifiers
    //!       and connect arrays to the method. The offset will only be added
    //!       once. 
    void SetElementPos( const PdeIdType identifierPDE1,
			Integer *connect1,
			Integer elemSize1,
			const PdeIdType identifierPDE2,
			Integer *connect2,
			Integer elemSize2,
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
    BaseGraph* GetGraph( const PdeIdType identifierPDE1 = NO_PDE_ID,
			 const PdeIdType identifierPDE2 = NO_PDE_ID );

    //! Get a specified (sub-)graph for inhom. Dirichlet BC

    //! After the assembly of the graph storing the connectivity between
    //! real dofs and dofs fixed by inhomogeneous Dirichlet boundary
    //! conditions is completed, the IDBC (sub-)graph for the specified
    //! PDE pair can be accessed to communicate directly with it.
    //! \param pdeID1 identifier for the first PDE of the pair
    //! \param pdeID2 identifier for the second PDE of the pair
    //! \return a BaseGraph object according to the specified identifier or
    //!         NULL, if no IDBC graph was generated
    BaseGraph* GetIDBCGraph( const PdeIdType pdeID1 = NO_PDE_ID,
                             const PdeIdType pdeID2 = NO_PDE_ID ) const;

    //! Obtain the permutation/re-ordering vector for a PDE

    //! A PDE attached to the GraphManagerStdMat can call this method to
    //! obtain the re-ordering vector for its unknowns, once assembly of the
    //! graph was completed by a call to SetupDone().
    //! \param identifier unique identifier for the PDE registered with the
    //!                   graph manager
    //! \return an integer array containing the new equation numbers after
    //!         the re-ordering or a NULL pointer, if no re-ordering was
    //!         performed.
    //! \note While memory for the permutation vector is allocated in this
    //!       method, it is the caller's responsibility to dispose of that
    //!       memory once it no longer needs the array.
    Integer *GetReordering( const PdeIdType identifier );

    //! Print some statistics on the graph manager
    void PrintStats( std::ostream *log );

    //@}


  private:

    //! Permutation vector stored as collection of one-based arrays

    //! We allocate memory for the permutation vectors of each PDE
    //! in this class, but do not assume that we are responsible for their
    //! de-allocation. This is the task of the respective PDE. However,
    //! if the PDE does not claim its permutation vector by calling
    //! GetReordering, we will de-allocate the memory in the destructor
    //! of this class.
    Integer **newOrdering_;

    //! Type of re-ordering to be applied to the graph

    //! This attribute stores the type of re-ordering that will be applied to
    //! the graph once the latter was completely assembled.
    //! \note Since the connectivity of all PDEs and all coupling objects is
    //! assembled into a single graph by this graph manager class, there can
    //! only be one type of re-ordering. The graph manager will stick to the
    //! reordering type specified by the PDE that registers first.
    ReorderingType reorderType_;

    //! Vector containing the offsets

    //! This vector stores for each PDE the offset that must be applied to
    //! the PDE's vertex numbers when inserting the connectivity information
    //! of the PDE and the coupling objects related to it. The offsets are
    //! stored in the number of registration.
    //! \note In order to be able to compute from the offsets the number of
    //!       vertices in the PDE sub-graph we add a virtual offset for an
    //!       (numPDEs_ + 1)st PDE.
    UInt *offset_;

    //! Vector containing offsets for graphIDBC_

    //! This vector stores for each PDE the offset that must be applied to
    //! the PDE's vertex numbers for fixed degrees of freedom when inserting
    //! the connectivity information of the PDE and the coupling objects
    //! related to it into the auxilliary IDBC graph. The offsets are
    //! stored in the number of registration.
    //! \note In order to be able to compute from the offsets the number of
    //!       vertices in the PDE sub-graph we add a virtual offset for an
    //!       (numPDEs_ + 1)st PDE.
    Integer *offsetIDBC_;

    //! Offset to be applied to vertexList1_ in SetElementPos

    //! This is the offset that is added in SetElementPos() to the vertex
    //! numbers in vertexList1_ for the equation numbers coming from the
    //! first PDE . The offset is determined when AssembleInit() is called.
    UInt offset1_;

    //! Offset to be applied to edgeList1_ in SetElementPos

    //! This is the offset that is added in SetElementPos() to the edge
    //! numbers in edgeList1_ for the equation numbers coming from the
    //! second PDE. The offset is determined when AssembleInit() is called.
    UInt offset2_;

    //! Offset to be applied to vertexList2_ in SetElementPos

    //! This is the offset that is added in SetElementPos() to the vertex
    //! numbers in vertexList2_ for the equation numbers coming from the
    //! first PDE . The offset is determined when AssembleInit() is called.
    Integer offsetIDBC1_;

    //! Offset to be applied to edgeList2_ in SetElementPos

    //! This is the offset that is added in SetElementPos() to the edge
    //! numbers in edgeList2_ for the equation numbers coming from the
    //! second PDE. The offset is determined when AssembleInit() is called.
    Integer offsetIDBC2_;

    //! Identifier for first PDE as passed to AssembleInit

    //! We use this attribute to store the identifier of the first PDE that
    //! was supplied to AssembleInit() for the insert of the connectivity of
    //! the PDE or an associated coupling object. We can use this identifier
    //! for performing consistency checks in SetElementPos() when debugging
    //! is enabled.
    PdeIdType idPDE1_;

    //! Identifier for second PDE as passed to AssembleInit

    //! We use this attribute to store the identifier of the second PDE that
    //! was supplied to AssembleInit() for the insert of the connectivity of
    //! an associated coupling object. We can use this identifier for
    //! performing consistency checks in SetElementPos() when debugging is
    //! enabled.
    PdeIdType idPDE2_;

    //! Pointer to the graph object for real/un-fixed dofs

    //! This graph object is used for storing the connectivity information
    //! for equation numbers for free dofs of the PDE.
    BaseGraph *graph_;

    //! Pointer to the graph object for dofs fixed by inhom. Dirichlet BC

    //! This graph object is used for storing the connectivity information
    //! between the equation numbers for free dofs of the PDE and those for
    //! dofs which are fixed by inhomogeneous Dirichlet boundary conditions.
    IDBC_Graph *graphIDBC_;

    //! Attribute to store the number of PDEs that belong to the manager
    UInt numPDEs_;

    //! Attribute to keep track of number of PDEs that were registered
    UInt numRegisteredPDEs_;

    //! Attribute to store the number of fixed degrees of freedom
    UInt numFixedDofs_;

    //! Attribute to store the number of free degrees of freedom
    UInt numFreeDofs_;

    //! Vector storing for each PDE the last free equation number
    UInt *numLastFreeDof_;

    //! Vector storing for each PDE the total number of equation numbers
    UInt *numEqn_;

    //! Attribute for testing wether it is safe to query the graph re-ordering

    //! Since there is a fundamental difference between the different graph
    //! manager implementations w.r.t. to the place at which the graph manager
    //! instructs the graph(s) to re-order themselves, we use this attribute
    //! for checking in GetReordering() whether it is safe to return the
    //! computed re-ordering.
    //! \note In the case of the %GraphManagerStdMat class the re-ordering of
    //! the graph is triggered by the call to SetupDone().
    bool reorderingDone_;

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
