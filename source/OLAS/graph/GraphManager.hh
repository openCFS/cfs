#ifndef GRAPH_MANAGER_HH
#define GRAPH_MANAGER_HH

#include <string>
#include "General/Environment.hh"
#include "OLAS/graph/BaseOrdering.hh"
#include <unordered_map>

namespace CoupledField {


// forward class declarations
class BaseGraph;
class IDBC_Graph;

   //! Class for administering the connectivity information of SBMblocks

  //! A graph manager is
  //! an object that is responsible for administering the graphs representing
  //! the connectivity information related to either a single block or e.g. in
  //! the case of direct coupled field simulations that of a collection of
  //! feFunctions.
  //! The graph manager serves as an intermediary between the graph
  //! object(s) and the feFunctions.
  //! Internally it stores for each SBM-Block a base-graph, which holds
  //! only the connectivity information of purely free equations, i.e. equations
  //! not fixed by a Dirichlet boundary condition. Each of these graph objects
  //! has a square shape.
  //! For equations fixed by an inhomogeneous Dirichlet boundary condition,
  //! the GraphManager holds a rectangular IDBCGraph object, which handles
  //! the connectivity of free equations to fixed equations.
  //! This coupling graph will be passed in the end to the IDBCHandler.
  //!
  //! The following piece of pseudo-code demonstrates the sequence of method
  //! calls that must be observed for setting up the graphs:
  //!
  //! \code
  //! /* The graph manager must be told how many SBM blocks lie in its
  //!    responsibility */
  //! SetupInit( numBlocks );
  //!
  //! /* Now each Block must register itself with the graph manager */
  //! for ( i = 1; i <= numBlocks; i++ ) {
  //!    RegisterBlock( ... );
  //! }
  //!
  //! /* Pass connectivity information to the GraphManager for all
  //!    elements */
  //! for ( i = 1; i <= numElems; i++ ) {
  //!
  //!    SetElementPos( ... );
  //! }
  //!
  //! /* By calling this method we inform the graph manager that all SBMBlocks
  //!    have defined their connectivity and that the setup phase now
  //!    is completed */
  //! SetupDone();
  //! \endcode
  class GraphManager {

  public:
    
    //! Struct containing information for SBM-system

    //! Lightweight struct for storing information per SBM Matrix block.
    typedef struct {

      //! Map (fctId,eqnNr) to matrix / vector index
      StdVector<std::unordered_map<UInt, UInt> > eqnToIndex;

      //! Define mappings for subBlocks
      //! 1st index: subBlockIndex
      //! 2nd vector contains indices for each subBlock
      StdVector<StdVector<UInt> > indexBlocks;

      //! Total number of equations in this block
      UInt size;
      
      //! Number of last free index in this block
      UInt numLastFreeIndex;

      //! Flag if minor/subBlocks are defined

      //! If the SBM matrix has minor blocks, no reordering can
      //! be applied!
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
    //! \param useDistinctGraphs if true, every matrixType (MASS, STIFFNESS)
    //!                          can have a different sparsity and subBlock
    //!                          pattern.
    //! \param isMultHarm        True if we perform a multiharmonic analysis
    //! \param N, M              Number of harmonics
    //! \param size              Number of considered frequencies (num of sbm blocks)
    void SetupInit( UInt numBlocks, bool useDistinctGraphs,
                    bool isMultHarm = false,
                    UInt N = 0,
                    UInt M = 0,
                    UInt size = 0,
                    bool isFullSys = false);

    //! Finalises setup of the graph manager

    //! This method must be called after the assembly of all sub-graphs was
    //! done, i.e. once all blocks have conveyed their
    //! connectivity information to the graph manager.
    //! \param reorder  re-ordering type for each block
    void SetupDone(const StdVector<BaseOrdering::ReorderingType>& reorder);

    //! Finalises setup of the graph manager in the multiharmonic case

    //! This method must be called after the assembly of all sub-graphs was
    //! done, i.e. once all blocks have conveyed their
    //! connectivity information to the graph manager.
    //! \param reorder  re-ordering type for each block
    //! \param N        Number of harmonics for solution
    //! \param M        Number of harmonics for nonlinearity
    //! \param isFullSys Boolean if we are considering all (even & odd) harmonics
    void SetupDoneMH(const StdVector<BaseOrdering::ReorderingType>& reorder,
                     const UInt N,
                     const UInt M);


    //! Register block with graph manager
    
    //! Before a matrix block can set its sparsity pattern, the matrix
    //! block has to be registered and the number of unknowns and 
    //! non free equations has to be set.
    //! \param blockNum       block number of matrix to be defined 
    //! \param blockInfo      pointer to structure containing the information
    //!                       how the SBMBlock is defined, i.e. a mapping
    //!                       of all equations defining this block and their
    //!                       col/row index, an (optional) definition of
    //!                       sub-blocks, as well as a flag, if the block
    //!                       can be reordered.
    //!
    //! \note A SBMBlock can only be reordered, if it has no subblocks defined,
    //!       as otherwise we can not guarantee for continuous indices within
    //!       one block.
    void RegisterBlock( const UInt blockNum,
                        SBMBlockInfo* blockInfo );


    //! Register block with graph manager for multiharmonic analysis

    //! Before a matrix block can set its sparsity pattern, the matrix
    //! block has to be registered and the number of unknowns and
    //! non free equations has to be set.
    //! \param blockInfo      pointer to structure containing the information
    //!                       how the SBMBlock is defined, i.e. a mapping
    //!                       of all equations defining this block and their
    //!                       col/row index, an (optional) definition of
    //!                       sub-blocks, as well as a flag, if the block
    //!                       can be reordered.
    //!
    //! \note A SBMBlock can only be reordered, if it has no subblocks defined,
    //!       as otherwise we can not guarantee for continuous indices within
    //!       one block.
    void RegisterBlockMultHarm( SBMBlockInfo* blockInfo);




    //@}

    // =======================================================================
    // SETUP OF SUB_GRAPH
    // =======================================================================

    //@{ \name Methods for generation of sub-graphs
    
    //! Insert connectivity of a finite element into the matrix graph

    //! This method is used to insert the connectivity of the unknowns
    //! attached to a finite element into the sub-graph representing the
    //! structure of the matrix associated to a single FeFunction or the 
    //! coupling between two FeFunctions. Thus the AlgebraicSystem
    //! has already performed the mapping from (feFctId,eqn) to row/column
    //! block numbers and indices therein.
    //! \param rowBlock       vector containing for every rowIndex the
    //!                       corresponding SBMMatrix block 
    //! \param rowIndices     vector containing all row indices for the 
    //!                       element
    //! \param colBlock       vector containing for every colIndex the
    //!                       corresponding SBMMatrix block
    //! \param colIndices     vector containing all column indices for the 
    //!                       element
    //! \param matrixType     matrix type (MASS, STIFFNESS) for which the 
    //!                       connectivity is set. Only meaningful, if 
    //!                       GraphManager was created with option
    //!                       distinctGraphs = true.
    //! \param setCounterPart if true, also the transposed connectivity 
    //!                       information is set
    virtual void SetElementPos( const StdVector<UInt>& rowBlock,
                                const StdVector<UInt>& rowIndices,
                                const StdVector<UInt>& colBlock,
                                const StdVector<UInt>& colIndices,
                                FEMatrixType matrixType,
                                bool setCounterPart);


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
    //! \param rowNum row number of the SBMBlock for the requested graph
    //! \param colNum column number of the SBMBlock for the requested graph
    //! \return a BaseGraph object according to the specified SBMBlock
    //! 
    //! \note The method refuses to return a pointer to a non-existant
    //!       graph and will report an error instead
    BaseGraph* GetGraph( UInt rowNum = 0,
                         UInt colNum = 0 );

    //! Get a specified (sub-)graph for inhom. Dirichlet BC

    //! After the assembly of the graphs storing the connectivity between
    //! real unknown entries and entries fixed by inhomogeneous Dirichlet 
    //! boundary conditions is completed, the specified IDBC (sub-)graph for
    //! a certain SBMBlock can be accessed to communicate directly with it.
    //! \param rowNum row number of the SBMBlock for the requested IDBCGraph
    //! \param colNum column number of the SBMBlock for the requested IDBCGraph
    //! \return a IDBCGraph object according to the specified SBMBlock
    BaseGraph* GetIDBCGraph( UInt rowNum = 0,
                             UInt colNum = 0 ) const;

    //! Obtain the permutation/re-ordering vector for a graph

    //! This method returns the permutation array of all unknowns of a 
    //! SBMBlock, once assembly of the graph was completed by a call 
    //! to AssembleDone().
    //! \param blockNum row/colNum of the SBMBlock for which the 
    //!                 reordering permutation vector is requested
    //! \param order  an integer array containing the new equation numbers after
    //!               the re-ordering 
    void GetReordering( UInt blockNum, StdVector<UInt>& order ) ;

    //! Print some statistics on the graph manager
    void PrintStats();

    //! Query whether a certain sub-graph exists

    //! This method can be used to query whether a sub-graph with the index
    //! pair (rowNum, colNum) was generated in the Setup phase.
    //! \return true if graph exists, false otherwise
    bool SubGraphExists( UInt rowNum = 0,
                         UInt colNum = 0 )  {
      return graph_[ ComputeIndex( rowNum, colNum ) ] != NULL;
    }

    //! Query whether a certain IDBC sub-graph exists

    //! This method can be used to query whether a sub-graph of IDBC type
    //! was generated for the specified SBMBlock in the setup phase.
    //! \return true if IDBC graph exists, false otherwise
    bool IDBCGraphExists( UInt rowNum = 0,
                          UInt colNum = 0 ) const {
      return graphIDBC_[ ComputeIndex( rowNum, colNum ) ] != NULL;
    }

    //@}


    //! Computes the index into the graph pointer matrix

    //! The class stores the pointers to the graphs associated with the 
    //! SBMBlocks and the coupling objects in a matrix. This matrix is stored
    //! as 1D array. This auxiliary method computes for a given pair of row 
    //! and column index the 1D index of the corresponding matrix entry.
    inline UInt ComputeIndex( UInt rowNum,
                              UInt colNum ) const {
      return numBlocks_*rowNum + colNum;
    }

  private:
    //! Auxiliary method for generation of IDBC graph objects

    //! This method will check for a SBMBlock defined by rowNum and colNum
    //! whether it is necessary to generate an associated IDBC graph. 
    //! This will be the case, if the second column block contains degrees of 
    //! freedom that are fixed by inhomogeneous Dirichlet boundary
    //! conditions.
    void GenerateIDBCGraph( UInt rowNum, UInt colNum);

    //! Auxiliary method for generation of graphs for coupled SBMBlocks

    //! This method can be called to generate the graph of a coupling object
    //! for the coupling between the SBMBlocks in the specified pair
    //! (rowNum, colNum). It is called by AssembleInit().
    void GenerateCouplingGraph( UInt rowNum, UInt colNum );

    //! Matrix storing pointers to the graph objects

    //! We use a flattened 2D structure to store for each diagonal SBMBlock 
    //! and their coupling a pointer to the associated graph object. 
    //! In analogy to the resulting %SBM_Matrix we use a 2D matrix for 
    //! storing the pointers.
    //! The graph belonging to the object with identifiers rowNum and colNum is
    //! stored as matrix entry (rowNum,colNum). We use C-style mapping to 
    //! transform the matrix to a 1D storage layout. Entries in the matrix may 
    //! be NULL pointers, if no corresponding coupling object exists.
    StdVector<BaseGraph *> graph_;

    //! Matrix storing pointers to the IDBC graph objects

    //! This matrix stores for each SBMBlock pair a pointer to an associated 
    //! graph for handling the connection between free degrees of freedom and
    //! those fixed by inhomogeneous Dirichlet boundary conditions.
    StdVector<IDBC_Graph *> graphIDBC_;

    //! Array containing pointers to permutation vectors

    //! We allocate memory for the permutation vector of each diagonal SBMBlock
    //! in this class.
    StdVector< StdVector<UInt> > newOrdering_;

    //! Vector storing for each block meta information
    StdVector<SBMBlockInfo*> blockInfo_;

    //! In multiharmonic analysis only one blockInfo is present
    SBMBlockInfo* blockInfoMH_;


    //! Attribute to store the number of SBMBlocks that belong to the manager
    UInt numBlocks_;

    //! Attribute to keep track of number of SBMBlocks that were registered
    UInt numRegisteredBlocks_;

    //! Attribute for testing whether it is safe to query the graph re-ordering

    //! This flag is only set to true, if all subGraphs are reordered and
    //! only then it is allowed to query the reordering array.
    //! This will be case after the GraphManager has been finalised by a call
    //! to AssembleDone().
    bool reorderingDone_;

    // =======================================================================
    // MULITHARMONIC VARIABLES
    // =======================================================================

    //@{ \name Variables for multiharmonic analysis

    //! Number of harmonics for solution quantity
    UInt N_;
    //! Number of harmonics for reluctivity
    UInt M_;
    //! Number of considered frequencies (positive and negative)
    //! which is the number of row/col blocks in the SBM matrix
    UInt sizeMH_;
    //! Flag wether we are considering all harmonics (even and odd) or the optimized version
    bool isFullSys_;
    //@}

    //! Flag indicating whether registration of all SBMBlocks is done

    //! This attribute is used as flag which indicates whether the
    //! registration of all SBMBlocks is finished.
    bool registrationDone_;

    //! Flag if multiharmonic analyis is performed
    bool isMultHarm_;

    //@{
    //! Auxiliary vector for passing connectivities to the graph objects

    //! This auxiliary vector is employed by the SetElementPos() method to
    //! pass the connectivity informations to the graph objects.
    std::vector<std::vector<UInt> > vertexList1_;
    std::vector<std::vector<UInt> > vertexList2_;
    std::vector<std::vector<UInt> > edgeList1_;
    std::vector<std::vector<UInt> > edgeList2_;
    //@}

  };

}

#endif
