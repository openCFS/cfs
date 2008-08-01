// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef BASE_GRAPH_MANAGER
#define BASE_GRAPH_MANAGER

#include <string>
#include <vector>
#include <iostream>

#include "utils/utils.hh"

namespace OLAS {

  // forward class declarations
  class BaseGraph;

  //! This is the base class for all graph manager classes. A graph manager is
  //! an object that is repsonsible for administering the graphs representing
  //! the connectivity information related to either a single PDE or e.g. in
  //! the case of direct coupled field simulations that of a collection of
  //! PDEs and the related coupling objects.
  //! The graph manager functions as an intermediary between the graph
  //! object(s) and the PDE and coupling objects.
  //!
  //! A typical example of a direct coupled simulation is the static
  //! piezo-electric problem. In this case the system matrix \f$A\f$
  //! can be written in the following form:
  //! \f[
  //! A = \left(\begin{array}{cc} M & P^T \\ P & -E \end{array}\right)
  //! \f]
  //! Here \f$M\f$ represents the FEM (stiffness) matrix derived from the
  //! mechanical PDE, \f$E\f$ that derived from the electrostatic PDE, while
  //! \f$P\f$ is the FEM discretisation of the (direct) PDE coupling.\n
  //! Depending on the solution process one wants to apply it is appropriate
  //! to setup \f$A\f$ either as an StdMatrix or in block form as an
  //! SBM_Matrix. In the firsrt case the matrix object expects a single graph
  //! object for setting up its sparsity pattern, while in the second case,
  //! each sub-matrix wants to talk to its individual subgraph.
  //!
  //! The distinction between these cases is handled by the different graph
  //! manager classes who manage the dissemination of connectivity information
  //! between the different sub-graphs. Thus CFS++ is mostly relieved from
  //! taking care of the two different approaches for representing the system
  //! matrix and instead can deal with a uniform interface provided by the
  //! methods of the %BaseGraphManager class.
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
  //! \note Since the PDE identifier is no longer assigned by the
  //!       GraphManager the derived classes must rely on the caller of
  //!       the public methods to satisfy the following requirements:
  //! - The PDE identifiers are of integral type
  //! - The PDE identifiers start from 1 and are assigned with out holes
  //!   in their sequence
  class BaseGraphManager {

  public:


    // =======================================================================
    // SETUP
    // =======================================================================

    //@{ \name Methods for construction, destruction and setup

    //! Default constructor
    BaseGraphManager(){};

    //! Destructor
    virtual ~BaseGraphManager(){
    };

    //! Prepares graph manager for generation of sub-graphs

    //! This method must be called after construction of the graph manager
    //! object and before any other method. It allows the graph manager to
    //! prepare itself for the registration and generation of the sub-graphs.
    //! \param numPDEs number of PDE objects whose connectivity graphs the
    //!                          graph manager will administer.
    virtual void SetupInit( UInt numPDEs ) = 0;

    //! Finalises setup of the graph manager

    //! This method must be called after the assembly of all sub-graphs was
    //! done, i.e. once all PDEs a coupling objects have conveyed their
    //! connectivity information to the graph manager.
    virtual void SetupDone() = 0;

    //! Register PDE with graph manager

    //! A PDE must call this method to register itself with the Graph manager
    //! and provide information on its number of unknowns. In the case of
    //! CFS++ this method  is not called directly, but via the algebraic
    //! system interface.
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
    virtual void RegisterPDE( const PdeIdType identifierPDE,
                              const UInt numEqns,
                              const UInt numLastFreeDof,
                              const ReorderingType reorder 
                              = NOREORDERING ) = 0;

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
    virtual void AssembleInit( const PdeIdType identifierPDE1,
                               const PdeIdType identifierPDE2,
                               bool assemblingTranspose ) = 0;

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
    virtual void AssembleDone( const PdeIdType identifierPDE1,
                               const PdeIdType identifierPDE2,
                               bool assemblingTranspose ) = 0;

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
                                Integer *connect1,
                                Integer elemSize1,
                                const PdeIdType identifierPDE2,
                                Integer *connect2,
                                Integer elemSize2,
                                bool setCounterPart ) = 0;

    //@}


    // =======================================================================
    // QUERY METHODS
    // =======================================================================

    //@{ \name Methods for querying information

    //! Obtain the permutation/re-ordering vector for a graph

    //! A PDE can call this method to obtain the re-ordering vector for its
    //! unknowns. Depending of the type of graph manager class. This method
    //! can either be called, once assembly of the sub-graph for the PDE was
    //! completed by a call to AssembleDone() or only after all sub-graphs for
    //! PDEs and their coupling have completely been assembled.
    //! \param identifier unique identifier for a PDE registered with the
    //!                   graph manager
    //! \return an integer array containing the new equation numbers after
    //!         the re-ordering or a NULL pointer, if no re-ordering was
    //!         performed.
    //! \note While memory for the permutation vector is allocated in this
    //!       method, it is the caller's responsibility to dispose of that
    //!       memory once it no longer needs the array.
    virtual Integer *GetReordering( const PdeIdType identifier ) = 0;

    //! Get a specified (sub-)graph

    //! After the assembly of the graph is completed, the specified
    //! (sub-)graph can be accessed to communicate directly with it.
    //! This is needed for example to create a sparse matrix, where the
    //! matrix graph has to be supported.
    //! \param identifierPDE1 identifier for first PDE related to sub-graph
    //! \param identifierPDE2 identifier for second PDE related to sub-graph
    //! \return a BaseGraph object according to the specified identifier
    //! combination
    virtual BaseGraph* GetGraph( const PdeIdType identifierPDE1 
                                 = NO_PDE_ID,
                                 const PdeIdType identifierPDE2
                                 = NO_PDE_ID) = 0;

    //! Get a specified (sub-)graph for inhom. Dirichlet BC

    //! After the assembly of the graph storing the connectivity between
    //! real dofs and dofs fixed by inhomogeneous Dirichlet boundary
    //! conditions is completed, the IDBC (sub-)graph for the specified
    //! PDE pair can be accessed to communicate directly with it.
    //! \param pdeID1 identifier for the first PDE of the pair
    //! \param pdeID2 identifier for the second PDE of the pair
    //! \return a BaseGraph object according to the specified identifier
    virtual
    BaseGraph* GetIDBCGraph( const PdeIdType pdeID1 = NO_PDE_ID,
                             const PdeIdType pdeID2 = NO_PDE_ID ) const = 0;

    //! Print some statistics on the graph manager
    virtual void PrintStats( std::ostream *log ) = 0;

    //@}


    // =======================================================================
    // AUXILLIARY METHODS
    // =======================================================================

    //! \name Auxilliary methods for manipulating equation numbers from CFS++
    //! The connectivity arrays passed to the various graph manager classes
    //! via the algebraic system interface (see: BaseSystem) contain equation
    //! numbers assigned by CFS++ to the degrees of freedom of the problem.
    //! The mapping of dof to equation number results in the following layout
    //! of equation numbers:\n\n
    //! <center><img src="../AddDoc/splitting.png"></center>\n
    //! The following methods are used to manipulate the entries in the
    //! connect arrays and adapt them to the needs of the graph manager class.
    //@{

    //! Generate vertex and edge lists from two connect arrays splitting info

    //! This auxilliary method takes as input two connect arrays and uses
    //! these to generate a vertex list and two edge lists that can be passed
    //! on to the graphs of a PDE. The first edge list is intended for the
    //! graph object managing the connectivity of the real dofs, while the
    //! second list is intended for the graph object that deals with the
    //! connectivity of the real dofs to those fixed by inhomogeneous
    //! Dirichlet boundary conditions.
    //!
    //! Both connect arrays contain equation numbers that in graph speak
    //! correspond to vertices. However, the first set of vertex numbers
    //! translates to row indices in the matrix graph, while the second one
    //! translates to column indices. Thus, connect1 contains equation numbers
    //! that identify vertices in the graph for which connect2 contains their
    //! neighbouring vertices.
    //!
    //! The method first builts from connect1 the vertexList by changing the
    //! sign of all equation numbers for dofs fixed by constraints and by
    //! then dropping all equation numbers that correspond to dofs fixed by
    //! homogeneous or inhomogeneous Dirichlet boundary conditions.
    //!
    //! It then generates the two edge lists from connect2 by changing the
    //! sign of all equation numbers for dofs fixed by constraints, dropping
    //! all dofs fixed by homogeneous Dirichlet boundary conditions and
    //! then inserting all equation numbers for free dofs into edgeList1
    //! and those for dofs fixed by inhomogeneous Dirichlet boundary
    //! conditions into edgeList2.
    //! \param connect1   equation numbers for vertices (i.e. row indices)
    //! \param length1    length of connect1 array
    //! \param connect2   equation numbers for  vertices (i.e. column indices)
    //! \param length2    length of connect2 array
    //! \param limit      marks the equation number of the last free dof
    //! \param vertexList on return contains the list of vertices for
    //!                   both graphs
    //! \param edgeList1  on return contains the edges (i.e. connectivity)
    //!                   for the vertices in the vertexList with respect to
    //!                   vertices for free dofs
    //! \param edgeList2  on return contains the edges (i.e. connectivity)
    //!                   for the vertices in the vertexList with respect to
    //!                   vertices for dofs fixed by inhomogeneous Dirichlet
    //!                   boundary conditions
    //! \param offsetVertexList  offset to add to the equation numbers
    //!                          in vertexList
    //! \param offsetEdgeList1   offset to add to the equation numbers
    //!                          in edgeList1
    //! \param offsetEdgeList2   offset to add to the equation numbers
    //!                          in edgeList2 (may be negative)
    inline void AdaptConnects( Integer *connect1, const UInt length1,
                               Integer *connect2, const UInt length2,
                               const UInt limit,
                               std::vector<UInt> &vertexList,
                               std::vector<UInt> &edgeList1,
                               std::vector<UInt> &edgeList2,
                               UInt offsetVertexList,
                               UInt offsetEdgeList1,
                               Integer offsetEdgeList2 ) {


      UInt aux;

      // Clear the arrays
      vertexList.clear();
      edgeList1.clear();
      edgeList2.clear();

      // STEP 1: Generate vertex list from first connect array, dropping
      //         equation numbers for dofs fixed by inhomogeneous Dirichlet
      //         boundary conditions and changing the sign of those fixed by
      //         constraints
      for ( UInt i = 1; i <= length1; i++ ) {
        aux = std::abs( (double) connect1[i] );
        if ( aux <= limit && aux > 0 ) {
          vertexList.push_back( aux + offsetVertexList );
        }
      }

      // STEP 2: Split the second connect array into two edge lists, one for
      //         the graph and one for the IDBCgraph (which handles the dofs
      //         fixed by inhomogeneous Dirichlet boundary conditions)
      for ( UInt i = 1; i <= length2; i++ ) {
        aux = std::abs( (double) connect2[i] );
        if ( aux > 0 ) {
          if ( aux > limit ) {
            edgeList2.push_back( UInt(aux + offsetEdgeList2) );
          }
          else {
            edgeList1.push_back( aux + offsetEdgeList1 );
          }
        }
      }

#ifdef DEBUG_BASEGRAPHMANAGER
      // output original connectivity
      (*debug) << "\n BaseGraphManager::AdaptConnects\n"
               << " numLastFreeDof   = " << limit            << '\n'
               << " offsetVertexList = " << offsetVertexList << '\n'
               << " offsetEdgeList1  = " << offsetEdgeList1  << '\n'
               << " offsetEdgeList2  = " << offsetEdgeList2  << '\n';
      (*debug) << " connect1 ";
      for ( UInt i = 1; i <= length1; i++ ) {
        (*debug) << connect1[i] << " ";
      }
      (*debug) << std::endl;
      (*debug) << " connect2 ";
      for ( UInt i = 1; i <= length2; i++ ) {
        (*debug) << connect2[i] << " ";
      }
      (*debug) << std::endl;

    // output original connectivity
      (*debug) << " vertexList: ";
      for ( UInt i = 0; i < vertexList.size(); i++ ) {
        (*debug) << vertexList[i] << " ";
      }
      (*debug) << std::endl;
      (*debug) << " edgeList1: ";
      for ( UInt i = 0; i < edgeList1.size(); i++ ) {
        (*debug) << edgeList1[i] << " ";
      }
      (*debug) << std::endl;
      (*debug) << " edgeList2: ";
      for ( UInt i = 0; i < edgeList2.size(); i++ ) {
        (*debug) << edgeList2[i] << " ";
      }
      (*debug) << std::endl;
#endif

    }

    //@}

  };

}

#endif
