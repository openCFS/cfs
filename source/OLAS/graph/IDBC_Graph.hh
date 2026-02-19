#ifndef OLAS_IDBC_GRAPH_HH
#define OLAS_IDBC_GRAPH_HH

#include <iostream>
#include <vector>

#include "BaseGraph.hh"

namespace CoupledField {

  //! Special graph class for handling IDBC graphs

  //! This class implements a special graph type that is used for the
  //! handling of graphs that describe the coupling between free degrees of
  //! freedom and those that are fixed by inhomogeneous Dirichlet boundary
  //! conditions. However, there are only few (currently one) methods that
  //! must be re-implemented, since the IDBC case requires special
  //! treatment.
  class IDBC_Graph : public BaseGraph {

  public:


    // =======================================================================
    // CONSTRUCTORS
    // =======================================================================
  
    //! \name Methods for construction, assembly and destruction
    //@{

    //! Constructor

    //! The constructor is a quite simple affair, since all it does is call
    //! the constructor of the base class setting NOREORDERING as strategy
    //! for re-ordering.
    //! \param numFreeDofs  number of vertices related to free degrees of
    //!                     freedom; this number is used as number of graph
    //!                     vertices and number of rows of associated matrix
    //! \param numFixedDofs number of vertices related to fixed degrees of
    //!                    freedom; this number is used as number of columns
    //!                    of associated matrix
    IDBC_Graph( UInt numFreeDofs, UInt numFixedDofs, unsigned int estimated_row_size) :
      BaseGraph( numFreeDofs, numFixedDofs, estimated_row_size) {
    }

    //! Default destructor
    virtual ~IDBC_Graph() {
    }

    //! Method to be called after all vertices were inserted into the graph

    //! This method must be called after all vertices were inserted into the
    //! graph. This finalises the assembly phase of the graph and triggers
    //! the re-ordering and the conversion into CRS storage format. The
    //! IDBC_Graph re-implements this method for the following reason.
    //! Once the graph object containing the coupling between all free
    //! degrees of freedom has re-ordered itself, the equation numbers for
    //! the free dofs have changed. The free dofs correspond to the vertices
    //! of the graph for which we store neighbourhood relations. Thus, these
    //! new equation numbers must be mapped to these vertex numbers.
    //! \note If newEqn = NULL we assume that the equation numbers have not
    //!       changed and that there is no need to re-number the vertices.
    //! \param newEqn one-based array containing the new equation numbers
    //!               of the free degrees of freedom
    void FinaliseAssembly( StdVector<UInt>* newEqn );

    private:

    //! Default constructor is not allowed
    IDBC_Graph( UInt numNodes, BaseOrdering::ReorderingType reorder ) {
      EXCEPTION( "Call to forbidden default constructor of "
               << "IDBC_Graph class" );
    }

    //! Copy constructor is not allowed
    IDBC_Graph( const IDBC_Graph &idbcGraph ) {
      EXCEPTION( "Call to forbidden copy constructor of "
               << "IDBC_Graph class");
    }

  };

}

#endif
