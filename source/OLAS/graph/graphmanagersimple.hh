// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef GRAPH_MANAGER_SIMPLE
#define GRAPH_MANAGER_SIMPLE

#include <iostream>
#include <string>
#include <vector>

#include "General/defs.hh"
#include "General/environment.hh"
#include "OLAS/graph/basegraphmanager.hh"
#include "OLAS/graph/baseordering.hh"
#include "Utils/StdVector.hh"


namespace CoupledField {


  //! Simple graph manager for the case of a single PDE

  //! This class implements the graph manager concept for the simple case that
  //! only a single PDE and only one associated graph exists.
class BaseGraph;
class IDBC_Graph;

  class GraphManagerSimple : public BaseGraphManager {

  public:

    // =======================================================================
    // SETUP
    // =======================================================================

    //@{ \name Methods for construction, destruction and setup

    //! Default constructor
    GraphManagerSimple();

    //! Destructor
    ~GraphManagerSimple();

    //! Prepares graph manager for generation of sub-graphs

    //! This method must be called after construction of the graph manager
    //! object and before any other method. It allows the graph manager to
    //! prepare itself for the registration and generation of the sub-graphs.
    //! In the case of this class there is nothing to be done. We only check
    //! that numPDEs equals 1 in order to ensure inappropriate use of the
    //! class.
    //! \param numPDEs number of PDE objects whose connectivity graphs the
    //!                          graph manager will administer.
    void SetupInit( UInt numPDEs );

    //! Finalises setup of the graph manager

    //! This method must be called after the assembly of all sub-graphs was
    //! done, i.e. once all PDEs a coupling objects have conveyed their
    //! connectivity information to the graph manager.
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
                      const BaseOrdering::ReorderingType reorder 
                      = BaseOrdering::NOREORDERING );

    //@}


    // =======================================================================
    // GENERATION OF SUB_GRAPHS
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

    //! Finialise assembly of a sub-graph

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
    //!                       inserted into the graph (ignored by this class)
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

    //! Obtain the permutation/re-ordering vector for a graph

    //! The PDE attached to the GraphManagerSimple can call this method to
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
    //! \return a BaseGraph object according to the specified identifier
    BaseGraph* GetIDBCGraph( const PdeIdType pdeID1 = NO_PDE_ID,
                             const PdeIdType pdeID2 = NO_PDE_ID ) const;

    //! Print some statistics on the graph manager
    void PrintStats( std::ostream *log );

    //@}


  private:

    //! Permutation vector stored as one-based array
    StdVector<UInt> newOrdering_;

    //! This attribute keeps track of memory de-allocation for newOrdering_

    //! We allocate the memory for the permutation vector, but do not assume
    //! that we are responsible for its de-allocation. This is the task of the
    //! PDE. However, if the PDE does not query the permutation vector, we
    //! will de-allocate the memory in this class. The attribute is used for
    //! the respective book-keeping
    bool newOrderingPassedOn_;

    //! Identifier of attached PDE

    //! This attribute stores the identifier of the PDE attached to the
    //! simple graph manager. Since there is only one PDE the identifier can
    //! be used to check for inappropriate use of the class.
    PdeIdType myPDEIdent_;

    //! Pointer to the graph object for real/un-fixed dofs

    //! This graph object is used for storing the connectivity information
    //! for equation numbers for free dofs of the PDE.
    BaseGraph *graph_;

    //! Pointer to the graph object for dofs fixed by inhom. Dirichlet BC

    //! This graph object is used for storing the connectivity information
    //! between the equation numbers for free dofs of the PDE and those for
    //! dofs which are fixed by inhomogeneous Dirichlet boundary conditions.
    IDBC_Graph *graphIDBC_;

    //! Type of re-ordering to be applied to the graph

    //! This attribute stores the type of re-ordering that will be applied to
    //! the graph once the latter was completely assembled.
    BaseOrdering::ReorderingType reorderType_;

    //! Attribute for testing wether it is safe to query the graph re-ordering

    //! Since there is a fundamental difference between the different graph
    //! manager implementations w.r.t. to the place at which the graph manager
    //! instructs the graph(s) to re-order themselves, we use this attribute
    //! for checking in GetReordering() whether it is safe to return the
    //! computed re-ordering.
    //! \note In the case of the %GraphManagerSimple class the re-ordering of
    //! the graph is triggered by the call to SetupDone().
    bool reorderingDone_;

    //! Auxilliary method for consistency checks in debug mode

    //! This auxilliary method is called by some methods in debug mode in
    //! order to ensure that the passed pde identifiers are valid for this
    //! class and that a graph object exists.
    //! \param  idPDE1  identifier for first PDE related to sub-graph
    //! \param  idPDE2  identifier for second PDE related to sub-graph
    //! \param  caller  string identifying the calling method for the
    //!                 generation of a proper error message
    void CheckConsistency( const PdeIdType idPDE1, const PdeIdType idPDE2,
                           std::string caller ) const;

    //! Total number of dofs (real and fixed by inhom. Dirichlet BC)

    //! This attribute stores the sum of the number of equation numbers
    //! that represent either real degrees of freedom or degrees of freedom
    //! that are fixed by inhomogeneous Dirichlet boundary conditions.
    UInt numEqns_;

    //! Equation number of last real / unfixed dofs

    //! This attribute stores the highest equation numbers that represent
    //! a real degree of freedom that is not fixed by constraints or
    //! Dirichlet boundary conditions.
    UInt numLastFreeDof_;

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

    //! Auxilliary method for assembling the permutation vector
    void BuildPermutationVector();
  };

}

#endif
