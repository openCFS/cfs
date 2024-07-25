#ifndef OLAS_ALGEBRAIC_SYS_HH
#define OLAS_ALGEBRAIC_SYS_HH

#include <map>
#include <set>

#include "General/defs.hh"
#include "General/Environment.hh"
#include "MatVec/Matrix.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/Exception.hh"
#include "OLAS/graph/GraphManager.hh"
#include "Utils/ThreadLocalStorage.hh"
//#include "OLAS/solver/BaseEigenSolver.hh"

namespace CoupledField {

  // Forward declarations
  class BaseSBMPrecond;
  class BaseSolver;
  class BaseEigenSolver;
  class BasePrecond;
  class BaseIDBC_Handler;
  class SBM_Matrix;
  class SBM_Vector;
  class BaseMatrix;
  class PatternPool;
  class StdMatrix;
  class BaseVector;
  class SingleVector;
  class SolStrategy;
  template <class TYPE> class StdVector;
  template <class TYPE> class Vector;


  //! Class responsible for solving a linear system of the form \f$\mathbf{A}x=b\f$l
  
  //! The purpose of this class is to administrate and assemble a linear
  //! algebraic matrix-vector system for FE simulation. Thus, it offers methods
  //! to assemble the element-wise contributions onto matrices or onto vectors.
  //! 
  //! Furthermore, it is responsible for solving the arising systems using
  //! (preconditioned) solvers. To allow for efficient solution strategies,
  //! it maintains internally a flexible mapping and separates thus the 
  //! FeFunction-oriented view of CFS from the matrix-oriented view of the 
  //! object oriented algebraic system (OLAS).
  //!
  //! OLAS internally uses so-called super-block-matrices (SBM, \ref SBM_Matrix) 
  //! which can be seen as dense matrices whose entries are again sparse
  //! matrices. Thus, the internal way of addressing entries is defined by a 
  //! combination of (sbmBlock, index).
  //! CFS on the other side uses a function/physic oriented view, i.e. unknowns
  //! are identified by a combination of FeFunctionId and a equation number
  //! associated with this FeFunction / FeSpace.
  //! The corresponding mapping
  //! \f[ \mathrm{(FeFunctionId, eqnNr) \rightarrow (sbmBlock, index)} \f]
  //! can an be externally defined. This allows different ways of splittings
  //! for different solution strategies. In addition, single SBM blocks
  //! can be re-ordered to minimize the bandwidth. 
  //!
  //! Furthermore we can define a sub-matrix structure, so called sub-blocks.
  //! which cann be used e.g. for block-preconditioners. They are only
  //! stored as additional indices within one matrix block.
  //! 
  //! \note The indices for FeFunctionId and SBM-blockindices are 0-based,
  //!       whereas equation numbers and row/col indices are 1-based, in
  //!       order to represent homogeneous values (Dirichlet values)
  //!       by a 0 equation number / index.
  class AlgebraicSys {
  public:

    //@{
    //! \name Public typedefs and helper classes
    
    //! Define constant for addressing all sbmBlocks
    static const Integer ALL_BLOCKS = -1;

    //! Triple identifying a sub-matrix of a Finite-Element matrix
    typedef struct {
      FEMatrixType feType;
      Integer rowInd;
      Integer colInd;
    } FeSubMatrixID;

    //! Tuple storing row and column index of an Finite-Element sub-matrix
    typedef struct {
      Integer rowInd;
      Integer colInd;
    } SubMatrixID;
    
    //! Definition of SBM matrix blocks by combination of (FeFctId,eqnNr)
    typedef std::map<FeFctIdType, std::set<Integer> > SBMBlockDef;

    //! For convenience and to help g++ parsing we define this shortcut
    typedef std::map<FEMatrixType,Double> factorMap;
    //@}

    //! Default Constructor
    AlgebraicSys( PtrParamNode param, PtrParamNode info, bool isSolutionComplex, bool isMultiHarm = false );

    //! Default Destructor
    virtual ~AlgebraicSys();

    //! Return struct with solution strategy
    shared_ptr<SolStrategy> GetSolStrategy() {return solStrat_;}
    
    //! Update algebraic system in multistep solution strategy
    void UpdateToSolStrategy();

    /** return current list of FeFunction names. We assume them to be ordered and 0-based */
    StdVector<std::string> GetFeFunctionsNames();
    
    /** return total number of equations.
     * If empty, the function is called too early */
    const StdVector<unsigned int>& GetFeFunctionsTotalEquations() const { return numEqnsPerFct_; }

    // ***********************************************************************
    //   Methods for creating and initializing the algebraic system
    // ***********************************************************************

    //@{
    //! \name Creation and Initialization

    //! Generate objects for the linear system(s)

    //! The method will generate the matrix and vector objects required for
    //! the different system matrices and the solution and right hand side
    //! vector.
    void CreateLinSys();

    //! Generate preconditioner
    void CreatePrecond();

    //! Generate solver object

    //! Calling this method triggers the generation of a solver object to
    //! solve a linear system of type Ax=b.
    //! Before this object can be used the SetupSolver() method should be
    //! called and also the generation of a preconditioner via the
    //! CreatePrecond() and SetupPrecond() methods should be finished, if
    //! a preconditioner is to be used.
    //! \note This method must not be called if an eigenfrequency analysis
    //! is performed, since this method creates only a solver to solve
    //! a system Ax=b.
    void CreateSolver();

    //! Generate EigenSolver object

    //! Calling this method triggers the generation of an EigenSolver object
    //! which ist used to calculate a generalized eigenvalue problem of the
    //! form M(du^2/dt^2)=Ku or of the quadratic eigenvalue problem, if
    //! a damping matrix is present.
    //! Before this object can be used the SetupEigenSolver() method should
    //! be called.
    //! \note If an Eigenfrequency analysis is performed, the methods
    //! SetupPrecond() and SetupSolver() must not be called!
    void CreateEigenSolver();

    //! Trigger setup of preconditioner.
    
    //! Calling this method will trigger the setup phase of the
    //! preconditioner. The setup is performed using the system matrix of the
    //! linear system.
    //! \note This method must not be called if an eigenfrequency analysis
    //! is performed, since this method creates only a preconditioner
    //! to solver which solves a system Ax=b. */
    void SetupPrecond();

    //! Checks if a preconditioner has been defined
    bool HasPrecond();

    //! Trigger setup of solution method.
    
    //! Calling this method will trigger the setup phase of the solver. This
    //! is especially important for direct solvers, where typically the
    //! factorisation of the problem matrix will be performed at this stage.
    //! The setup is performed using the system matrix of the linear system.*/
    void SetupSolver();

    //! Trigger setup of eigenvalue solver

    //! Calling this method will trigger the setup phase of the eigensolver.
    //! The setup is performed using the system matrix of the linear system.
    //! \param numFreq Number of eigenfrequencies to be calculated
    //! \param shift Frequency shift applied to the system
    //! \param quadratic Flag indicating if a quadratic eigenvalue problem
    //!        (true) or a generalized problem (false) is to be solved
    //! \param sort shall the evs be sorted?
    //! \param bloch mode problems are complex but not quadratic
    void SetupEigenSolver(UInt numFreq, Double shift, bool quadratic, bool sort, bool bloch);

    //! Solve the linear system.
    //! Calling this method triggers the solution of the linear system
    //! defined by the Finite-Element system matrix, i.e. sysMat_[SYSTEM],
    //! and the given right-hand side vector, i.e. #rhs_. The computed
    //! (approximate) solution is stored in #sol_.
    //! The method used for solving the system depends on the solver and
    //! preconditioner constructed and the parameters in myParams_.
    //! For iterative solvers an initial guess is created by inserting the
    //! Dirichlet values in the correct positions of the solution vector in
    //! case of the penalty formulation.
    //! \param setIDBC true: considers inhomog. Dirichlet b.c. (standard)
    //!                false: is needed in case of nonlinear PDE, incremental formulation!
    //! \note This method must not be called if an eigenfrequency analysis
    //! is performed, since this method is only used to solve a system of the
    //! form Ax=b.*/
    //! if deltaIDBC = true, instead of the current IDBC values idbc_cur,
    //!                         the delta between current and old values idbc_old
    //!                         is set
    void Solve(bool setIDBC = true, bool deltaIDBC = false );

    //! Calculate eigenfrequencies of a generalized eigenvalue problem

    //! Calling this method triggers the solution of a generalized eigenvalue
    //! problem, i.e. \f[ \mathbf{M} \ddot{x} = \mathbf{K} x \f].
    //! The number of eigenfrequencies is specified by setting the parameter
    //! "numEigenFreqs" of the OLAS_Params object.
    //! The resulting eigenfrequencies are stored into an array, as well as the
    //! resulting error.
    //! \param frequencies Read-only buffer which contains the eigenfrequencies
    //!                    of the generalized eigenvalue problem. The eigenvalues are (2*pi*f)^2 !!
    //!                    The size of the vector is the number of converged ev
    //! \param err Reached error norm for each eigenvalue
    //! \note This method may only be calaled if SetupEigenfrequencySolver()
    //!       was called previously.
    void CalcEigenFrequencies(Vector<Double>& frequencies, Vector<Double>& err);
 
    //! Calculate eigenfrequencies of a quadratic eigenvalue problem

    //! Calling this method triggers the solution of a quadratic eigenvalue
    //! problem, where the resulting frequencies are always complex.
    //! The number of eigenfrequencies is specified by setting the parameter
    //! "numEigenFreqs" of the OLAS_Params object.
    //! The resulting eigenfrequencies are stored into an array, as well as the
    //! resulting error.
    //! \param frequencies Read-only buffer which contains the eigenfrequencies
    //!                    of the generalized eigenvalue problem
    //! \param err Reached error norm for each eigenvalue
    //! \return Number of converged eigenvalues
    //! \note This method may only be called if SetupEigenfrequencySolver()
    //!       was called previously.
    void CalcEigenFrequencies(Vector<Complex>& frequencies, Vector<Double>& err);
    void CalcEigenValues(BaseVector &sol, BaseVector &err, Double minVal, Double maxVal);
    void CalcEigenValues(BaseVector &sol, BaseVector &err);

    //! Calculate eigenmodes of a generalized eigenvalue problem

    //! Calling this method triggers the calculation of the n-th eigenmode
    //! corresponding to the n-th eigenvalue.
    //! This method may only be called after CalcEigenfrequencies(), as
    //! the eigenmodes are only a postprocessing result.
    //! The calculated mode can be accessed by GetSolutionVal().
    //! \param numMode Mode number corresponding to the numMode-th
    //!                eigenfrequency
    //! \note This method may only ba called if SetupEigenfrequencySolver()
    //!       and CalcEigenfrequencies was called previously.
    void GetEigenMode(UInt numMode, bool right=true );

    //@}

    // ***********************************************************************
    //   Methods for creating and assembling the matrix graph
    // ***********************************************************************

    //@{
    //! \name Matrix Graph Methods

    //! This method must be called after construction of the graph manager
    //! object and before any other method. It allows the graph manager to
    //! prepare itself for the registration and generation of the sub-graphs.
    //! \param numFcts number of FeFunction objects 
    //! \param useDistinctGraphs flag, indicating if different matrix types
    //!                          (STIFFNESS, MASS) should have different
    //!                          sparsity patterns (e.g. pure (block)diagonal,
    //!                          CRS etc.)
    //!
    void GraphSetupInit( UInt numFcts, 
                         bool useDistinctGraphs );

    //! Finalises setup of the graph manager

    //! This method must be called after the assembly of all sub-graphs was
    //! done, i.e. once all functions or coupling objects have conveyed their
    //! connectivity information to the graph manager.
    void GraphSetupDone();

    //! Finalises setup of the graph manager

    //! Special graph setup for multiharmonic case!
    //! This method must be called after the assembly of all sub-graphs was
    //! done, i.e. once all functions or coupling objects have conveyed their
    //! connectivity information to the graph manager.
    void GraphSetupDoneMH();


    //! Obtain a unique FeFunction identifier

    //! Return a unique FeFunction identifier for each FeFunction, which will be used
    //! for all subsequent communication with the algebraic system.
    //! \param fctString string which identifies the feFunction
    //! \return unique identifier for the feFunction
    FeFctIdType ObtainFctId( const std::string& fctString );

    //! Register FeFunction with graph manager

    //! A feFunction must call this method to register itself with the graph 
    //! manager. 
    //! \param fctId       unique identifier for the feFunction
    //! \param numEqns     number of unknowns in the linear system associated
    //!                    with the feFunction (i.e. number of equation numbers 
    //!                    from openCFS)
    //! \param numLastFreeEqn specifies the equation number of the last
    //!                       unfixed degree equation (if inhomogeneous
    //!                       Dirichlet boundary conditions are treated by the
    //!                       penalty method numLastFreeEqn
    //!                       \f$  \stackrel{!}{=}  \f$  numEqns)
    void RegisterFct( const FeFctIdType fctId, UInt const numEqns,
                      UInt const numLastFreeEqn );

    //! Define entries of one SBM matrix block 

    //! A SBMBlock is physically a sparse matrix, which can be composed
    //! of entries of several feFunctions. With this function we explicitly
    //! define the entries of one of those blocks.
    //! \param eqns contains the set of equations for one super matrix block
    //! \param isInnerBlock flag, if block contains inner functions which can
    //!                     get eliminated using static condensation
    //! \return index of the matrix within the SuperBlockMatrix
    //! 
    //! \note Even in case of just one SBM matrix block, a call to this method
    //!       is mandatory.
    //! \note If a matrix block contains no equations at all, the return value
    //!       of this method is -1, as OLAS internally does not allow empty
    //!       matrices  
    Integer DefineSBMMatrixBlock( const std::map<FeFctIdType, 
                                  std::set<Integer> >& eqns,
                                  bool isInnerBlock );

    //! Register submatrix blocks within a sparse matrix

    //! Within a sparse matrix we allow additional blocks (e.g. for block-
    //! preconditioners / solvers). The restriction on those blocks is, that
    //! they have to be sequentially ordered.
    //! \param sbmIndex index of the super block matrix
    //! \param numMinorBlocks number of minor blocks
    void RegisterSubMatrixBlocks( UInt sbmIndex, UInt numMinorBlocks );


    //! Define minor blocks within a sparse matrix

    //! Define the range of one matrix sub-block within a sparse matrix
    void DefineSubMatrixBlocks( UInt sbmIndex, UInt blockIndex,
                                const StdVector<FeFctIdType>& fctIds,
                                const StdVector<Integer>& eqns );
    
    // Alternative new interface, which is easier
    //void DefineSubMatrixBlocks( UInt sbmIndex, UInt blockIndex,
    //                            const StdVector<std::pair<FeFctIdType,Integer> >& fctIds );
    
    
    //! Finalize registration and definition of SBMMatrix / Submatrix blocks
    void FinishRegistration( );

    //! Insert connectivity of a finite element into the matrix graph

    //! This method is used to insert the connectivity of the unknowns
    //! attached to a finite element into the sub-graph representing the
    //! structure of the matrix associated to a single Fct or the coupling
    //! between two Fcts. In the case that the sub-graph belongs to a Fct only
    //! the corresponding identifier (obtained from RegisterFct) needs to be
    //! supplied. In the case that the graph describes the coupling between
    //! two Fcts both related identifiers must be specified.
    //! \param identifierFct1 identifier for first Fct related to sub-graph
    //! \param eqnNrs1       array containing equation numbers for the
    //!                       degrees of freedom attached to the element
    //!                       (unknowns at nodes or edges, whatsoever), must
    //!                       be a one based pointer; for first Fct
    //! \param identifierFct2 identifier for second Fct related to sub-graph
    //! \param eqnNrs2        array containing equation numbers for the
    //!                       degrees of freedom attached to the element
    //!                       (unknowns at nodes or edges, whatsoever), must
    //!                       be a one based pointer; for second Fct
    //! \param matrixType     type of FE matrix (SYSTEM, STIFFNESS, ..), if
    //!                       we have different sparsity patterns for different
    //!                       matrix types
    //! \param setCounterPart flag signaling whether the counterpart, i.e.
    //!                       the connectivity with the Fct identifiers and
    //!                       equation numbers reversed, should also be
    //!                       inserted into the graph
    void SetElementPos( const FeFctIdType identifierFct1,
                        const StdVector<Integer>& eqnNrs1,
                        const FeFctIdType identifierFct2,
                        const StdVector<Integer>& eqnNrs2,
                        FEMatrixType matrixType,
                        bool setCounterPart );

    //! Helper method for mapping (fctId, eqnNr) to (blockNum,index)
    
    //! This method maps the CFS-oriented/physical oriented tuple (fctId,eqnNr)
    //! to the OLAS-block oriented (blockNum, row/colIndex). This allows a 
    //! transparent re-mapping and re-ordering of the block structure
    //! independent of the physical view
    //! \param fctId identifier of the feFunction to be mapped
    //! \param eqns vector containing the equations numbers related to the 
    //!             feFunction
    //! \param blockNum contains for each equation the block number
    //! \param indices contains for each equation the row/col index in the 
    //!                related \blockNum block
    //! \note equation numbers and indices within OLAS are 1-based!
    void MapFctIdEqnToIndex( const FeFctIdType fctId,
                             const StdVector<Integer>& eqns,
                             StdVector<UInt>& blockNums,
                             StdVector<UInt>& indices );

    //! @see MapFctIdEqnToIndex()
    void MapFctIdEqnToIndex_MultHarm(const FeFctIdType fctId,
                                     const StdVector<Integer>& eqns,
                                     StdVector<UInt>& blockNums,
                                     StdVector<UInt>& indices,
                                     const StdVector<UInt>& sbmIndices );

    //! Helper method for inserting an element method into the correct block
    //! in a reordered multiharmonic analysis
    void MapFctIdEqnToIndex_MultHarm( const FeFctIdType fctId,
                                      const StdVector<Integer>& eqns,
                                      StdVector<UInt>& blockNums,
                                      StdVector<UInt>& indices );


    //! @see MapFctIdEqnToIndex()
    void MapFctIdEqnToIndex( const StdVector<FeFctIdType>& fctId,
                             const StdVector<Integer>& eqns,
                             StdVector<UInt>& blockNums,
                             StdVector<UInt>& indices );

    //! @see MapFctIdEqnToIndex()
    //! We only map one specific (fctId,eqn) to (blockNum,index)
    void MapFctIdEqnToIndex( const FeFctIdType fctId,
                             const Integer eqn,
                             UInt& blockNum,
                             UInt& index );

    //! @see MapFctIdEqnToIndex()
    //! Return all (blockNum,index)-pairs for a given fctId 
    void MapCompleteFctIdToIndex( const FeFctIdType fctId,
                                  StdVector<UInt>& blockNums,
                                  StdVector<UInt>& indices );
    
    //! @see MapFctIdEqnToIndex()
    //! Return all for each SBM block all indices for a given fctId.
    //! In case of the elimination approach to handle IDBCs, the second
    //! set contains all indices corresponding to the auxiliary matrices
    //! stored in the IDBC handler. 
    //! This method is used to assemble e.g. the effective matrix, where
    //! just a sub-set of the indices of a sparse matrix is addressed.
    //! \param fctId identifier of the feFunction to be mapped
    //! \param freeIndPerBlock contains for each SBM-block (key to first map)
    //!                        a set of indices (1-based) of rows / cols 
    //!                        corresponding to the non-fixed equations of a
    //!                        feFunction.
    //! \param fixedIndPerBlock in case of the IDBC elimination approach, it 
    //!                         contains for each SBM-block (key to first map)
    //!                         of the IDBC auxiliary matrix a set of indices 
    //!                         (0-based) of rows corresponding to 
    //!                         the fixed equations of a feFunction.
    //! \param indexZeroBased flag indicating if the indices should be numbered
    //!                       0-based (= are dircetly used as indices in a 
    //!                       matrix) or if they are used in conjunction with 
    //!                       the graph manager / SetElementMatrix method 
    //!                       (1-based indices) 
    void MapCompleteFctIdToIndex( const FeFctIdType fctId,
                                  std::map<UInt, 
                                  std::set<UInt> >& freeIndPerBlock,
                                  std::map<UInt, 
                                  std::set<UInt> >& fixedIndPerBlock,
                                  bool indexZeroBased );

    
    // ***********************************************************************
    //   Methods for assembling the matrices, rhs and bc
    // ***********************************************************************

    //@{
    //! \name Assembling Routines

    //! Initialize specified matrix with zeros

    //! Calling this method initializes the Finite-Element matrix designated
    //! by the provided FEMatrixType identifier with zeros.
    //! After initialization all matrix entries in the matrix' sparsity
    //! pattern are zero. If no FEMatrixType identifier is specified all
    //! matrices in the internal set (see #matrixTypes_) will be initialized.
    //! In the case of an SBM_System the FeFctIdType identifier can be used to
    //! initialize only a sub-matrix of the Finite Element matrix. If no
    //! FeFctIdType identifier is specified the complete matrix will be
    //! initialized.
    //! \param matrixType type of Finite Element matrix
    //!                   (STIFFNESS, MASS, ...)
    //! \param fctId unique identifier obtained from the ObtainFctId()
    //!                      method
    void InitMatrix( FEMatrixType matrixType = NOTYPE,
                     const FeFctIdType fctId = NO_FCT_ID );

    //! Set right hand side vector to zero

    //! This method sets the right hand side vector associated with the
    //! specified Fct to zero. If no Fct identifier is given, the complete
    //! right-hand side vector is zeroed.
    //! \param fctId unique identifier obtained from the ObtainFctId()
    //!              method
    void InitRHS( const FeFctIdType fctId = NO_FCT_ID );

    //! Set the global RHS to given values

    //! Set the right hand side vector to the values specified in the 
    //! SBM-vector \a newRHS. The SBM-vector contains for each FeFunction
    //! a sub-vector. 
    //! \param newRHS SBM-Vector with RHS for all FeFunctions of the system
    //! \note The vector \a newRHS has the same layout as the one returned
    //!       by GetRHSVal(SBM_Vector)
    void InitRHS( const SBM_Vector& newRHS );

    //! Set the solution vector of the specified Fct to zero

    //! This method sets the part of the solution vector associated with the
    //! specified Fct to zero. If no Fct identifier is given, the complete
    //! solution vector is zeroed.
    //! \param fctId unique identifier obtained from the ObtainFctId()
    //!              method
    void InitSol( const FeFctIdType fctId = NO_FCT_ID );

    //! Set the solution vector to the values specified in the 
    //! SBM-vector \a newRHS. The SBM-vector contains for each FeFunction
    //! a sub-vector. 
    //! \param newSol SBM-Vector with solution for all FeFunctions of the 
    //!               system
    //! \note The vector \a newRHS has the same layout as the one returned
    //!       by GetSolutionVal(SBM_Vector)
    void InitSol( const SBM_Vector& newSol );

    //! Assemble an element matrix into the global one

    //! This methods assembles the given element matrix into a specified
    //! global one (MASS, STIFFNESS, etc.), i.e. adds the entries to the.
    //! global matrix. For element matrices which are associated
    //! only with one Fct, only one identifier, the eqn numbers of the matrix
    //! and the number of eqnNrs have to be specified.
    //! For matrices which are associated with two different Fct identifiers,
    //! the matrix will be assembled into the upper off-diagonal matrix-block
    //! and the lower transposed one.
    //! \param matrixType type of finite element destination matrix
    //!                 (STIFFNESS, MASS, ...)
    //! \param elemmat entries of the element matrix 
    //! \param fctId1 identifier for first Fct related to sub-graph
    //! \param eqnNrs1 equation numbers (1-based) of the element matrix
    //!                w.r.t. sub-graph associated with identifierFct1
    //! \param numEqn1 number of equations related to sub-graph of
    //!                identifierFct1
    //! \param fctId2 identifier for first Fct related to sub-graph
    //! \param eqnNrs2 equation numbers (1-based) of the element matrix
    //!                w.r.t. sub-graph associated with identifierFct2
    //! \param numEqn2 number of equations related to sub-graph of
    //!                identifierFct2
    //! \param setCounterPart if this flag is true, then the method will
    //!                not only insert the element matrix    \f$  E  \f$ , but also
    //!                its transpose    \f$  E^T  \f$ . In doing so also the
    //!                row and column indices derived from the equation
    //!                numbers are interchanged. Note that this is only
    //!                supported for off-diagonal blocks, i.e. for cases
    //!                with different Fct identifiers.
    //! \param noStaticCond If set to false, no static condensation will be
    //!                     applied to this element matrix. This is needed e.g.
    //!                     for matrices not being assembled to the system matrix.
    //! \param isDiagonal if this flag is true, then fctId1 == fctId2
    //!                   and eqnNrs1 == eqnNrs2. In this case the
    //!                   assembly get speed up.
    template<typename T>
    void SetElementMatrix( FEMatrixType matrixType, 
                           Matrix<T>& elemmat,
                           FeFctIdType fctId1,
                           const StdVector<Integer>& eqnNrs1,
                           FeFctIdType fctId2,
                           const StdVector<Integer>& eqnNrs2,
                           bool setCounterPart, 
                           bool noStaticCond,
                           bool isDiagonal );


    //! Assemble an element matrix into the global one for multiharmonic analysis

    //! This methods assembles the given element matrix into a specified
    //! global one (MASS, STIFFNESS, etc.), i.e. adds the entries to the.
    //! global matrix. For element matrices which are associated
    //! only with one Fct, only one identifier, the eqn numbers of the matrix
    //! and the number of eqnNrs have to be specified.
    //! For matrices which are associated with two different Fct identifiers,
    //! the matrix will be assembled into the upper off-diagonal matrix-block
    //! and the lower transposed one.
    //! \param matrixType type of finite element destination matrix
    //!                 (STIFFNESS, MASS, ...)
    //! \param elemmat entries of the element matrix
    //! \param fctId1 identifier for first Fct related to sub-graph
    //! \param eqnNrs1 equation numbers (1-based) of the element matrix
    //!                w.r.t. sub-graph associated with identifierFct1
    //! \param numEqn1 number of equations related to sub-graph of
    //!                identifierFct1
    //! \param fctId2 identifier for first Fct related to sub-graph
    //! \param eqnNrs2 equation numbers (1-based) of the element matrix
    //!                w.r.t. sub-graph associated with identifierFct2
    //! \param numEqn2 number of equations related to sub-graph of
    //!                identifierFct2
    //! \param setCounterPart if this flag is true, then the method will
    //!                not only insert the element matrix    \f$  E  \f$ , but also
    //!                its transpose    \f$  E^T  \f$ . In doing so also the
    //!                row and column indices derived from the equation
    //!                numbers are interchanged. Note that this is only
    //!                supported for off-diagonal blocks, i.e. for cases
    //!                with different Fct identifiers.
    //! \param sbmIndices SBM-indices of the sbm-blocks,
    //!                   which have to be assembled
    template<typename T>
    void SetElementMatrix_MultHarm( FEMatrixType matrixType,
                                   Matrix<T>& elemmat,
                                   FeFctIdType fctId1,
                                   const StdVector<Integer>& eqnNrs1,
                                   FeFctIdType fctId2,
                                   const StdVector<Integer>& eqnNrs2,
                                   bool setCounterPart,
                                   const StdVector<UInt>& sbmIndices );

    //! Assemble the local rhs vector to the global one

    //! This method adds the entries of the element right hand side to the
    //! global ones, according the FctIdentifier and the associated equation
    //! numbers.
    //! \param elemRHS entries of the element right hand side in a sequential
    //!                array
    //! \param fctId   identifier of the Fct related to sub-graph
    //! \param eqnNrs  equation numbers (1-based) of the element rhs
    //!                w.r.t. sub-graph associated with idFct
    //! \param harm  for multiharmonic analysis, we need the info
    //!                which harmonic we are currently considerng
    template<typename T>
    void SetElementRHS( const Vector<T>& elemRHS, 
                        const FeFctIdType fctId,
                        StdVector<Integer>& eqnNrs,
                        UInt& harm);


    //! Adds a value to a given global rhs entry

    //! This method adds a given value to the global right hand side vector.
    //! Therefore a Fct identifier, the related equation number and an
    //! equation number have to be specified.
    //! \param val value to be added
    //! \param fctId identifier of the Fct related to sub-graph
    //! \param eqnNr equation number of the node to be set
    template<typename T>
    void SetNodeRHS(T val, FeFctIdType fctId,
                    Integer eqnNr);

    //! Adds an feFunction to the right hand side
    template<typename T>
    void SetFncRHS( const Vector<T>& fncRHS, FeFctIdType fctId);

    //! Performs a matrix-vector multiplication and adds the vector to the rhs

    //! This method multiplies a specified global matrix (STIFFNES, MASS,
    //! etc.) with a given vector and adds the result vector to the global
    //! rhs.
    //! \f[ rhs = rhs + \mathbf A_{matrixType} \cdot fup \f]
    //! This method is mainly used in a time step scheme to adapt the rhs.
    //! \param matrixType type of finite element matrix (STIFFNESS, MASS, ...)
    //!                 which gets multiplied
    //! \param fup array with vector entries, which get multiplied
    //! \param SysMatUpdated indicates if we need to allocate new memory for the tmpRHS_ vector
    void UpdateRHS(FEMatrixType matrixType, const SBM_Vector& fup,bool SysMatUpdated, bool useTransposed=false);

    void ComputeSysMatTimesVector(FEMatrixType matrixType, SBM_Vector& inputVec, SBM_Vector& outputVec, bool transpose);

    //! Performs a matrix-vector multiplication and adds the vector to the rhs

    //! This method multiplies a specified global matrix (STIFFNES, MASS,
    //! etc.) with a given vector and adds the result vector to the global
    //! rhs.
    //! \f[ rhs = rhs + \mathbf A_{matrixType} \cdot fup \f]
    //! This method is currently only used to adapt the rhs for nonlinear multiharmonic analysis.
    //! \param matrixType type of finite element matrix (STIFFNESS, MASS, ...)
    //!                 which gets multiplied
    //! \param fup array with vector entries, which get multiplied
    //! \param SysMatUpdated indicates if we need to allocate new memory for the tmpRHS_ vector
    void UpdateRHS_MultHarm(FEMatrixType matrixType, const SBM_Vector& fup,bool SysMatUpdated);

    //! Add a value to a diagonal matrix entry

    //! This method allows to directly add a value to the value of an entry
    //! on the main diagonal of a certain of the different problem matrices.
    //! \param matrixType specifies which matrix is to be manipulated
    //! \param fctId    identifier for Fct related to sub-graph
    //! \param eqnNum   determines the diagonal entry to be manipulated, i.e.
    //!                 the parameter gives the equation number of the unknown
    //!                 coupled to the row/column of the entry; e.g. for
    //!                    \f$  \mbox{eqnNum}=k  \f$  it is    \f$  a_{kk}  \f$  that will be
    //!                 altered
    //! \param val      zero-based array containing the numerical value that
    //!                 is to be added to the matrix entry; in the case of
    //!                 real entries, the first array value is used, in the
    //!                 case of complex entries the first entry is used as
    //!                 real and the second as imaginary part
    template<typename T>
    void AddToDiagMatrixEntry( FEMatrixType matrixType,
                               const FeFctIdType fctId,
                               Integer eqnNum,
                               T val );
    
    //! Get value of specific matrix entry

    //! This method allows to directly get the value of a specific matrix
    //! entry of a given problem matrix. If the given index pair is not
    //! contained in the matrix, an error is thrown.
    //! \param matrixType  specifies which matrix is to be queried
    //! \param rowFctId  identifier for Fct related to row sub-graph
    //! \param rowEqnNum equation number of the row index to get
    //! \param colFctId  identifier for Fct related to column sub-graph
    //! \param colEqnNum equation number of the column index to get
    //! \param val       return value
    //! \note For sparse matrices this method may be very costly and slow,
    //!       as the given index pair has to be searched for and may not
    //!       be accessible directly!
    template<typename T>
    void GetMatrixEntry( FEMatrixType matrixType,
                         const FeFctIdType rowFctId,
                         Integer eqnNum1,
                         const FeFctIdType colFctId,
                         Integer eqnNum2,
                         T & val );

    //! Set value of specific matrix entry

    //! This method allows to directly set the value of a specific matrix
    //! entry of a given problem matrix. If the given index pair is not
    //! contained in the matrix, an error is thrown.
    //! \param matrixType  specifies which matrix is to be manipulated
    //! \param rowFctId  identifier for Fct related to row sub-graph
    //! \param rowEqnNum equation number of the row index to be set
    //! \param colFctId  identifier for Fct related to column sub-graph
    //! \param colEqnNum equation number of the column index to be set
    //! \param val       value to be set
    //! \param setCounterPart if this flag is true, then the method will
    //!                not only set the value on position (eqnNum1, eqnNum2)
    //!                but also on its transposed position (eqnNum2, eqnNum1).
    //!                If the entry is located on an off-diagonal entry (i.e.
    //!                different FctIDs) the counterpart in the related
    //!                transposed block is inserted.
    //! \note For sparse matrices this method may be very costly and slow,
    //!       as the given index pair has to be searched for and may not
    //!       be accessible directly!
    template<typename T>
    void SetMatrixEntry( FEMatrixType matrixType,
                         const FeFctIdType rowFctId,
                         Integer rowEqnNum,
                         const FeFctIdType colFctId,
                         Integer colEqnNum,
                         T val,
                         bool setCounterPart );

    //! Constructs the effective system matrix related to a given fctId

    //! Constructs the effective system matrix for a given functionId by 
    //! multiplying the individual matrices (except AUXILIARY!) with 
    //! factors and adding them to one system matrix
    //! \param fctId identifier for function related to sub-graph
    //! \param matFactors a map which contains for each matrixtype (DAMPING,
    //!                   STIFFNESS,...) the according factor
    //! \param isMultHarm flag if multiharmonic analysis
    void ConstructEffectiveMatrix( 
        const FeFctIdType fctId,
        const std::map<FEMatrixType, Double>& matFactors,
        const bool isMultHarm = false);

    //! Pass a Dirichlet value to %OLAS

    //! This method passes the value of a given inhomogeneous Dirichlet
    //! equation number to %OLAS.
    //! \param fctId  identifier for Fct related to sub-graph
    //! \param eqnNr  equation number of the dirichlet restraint
    //! \param val    value of the dirichlet restraint
    //!
    //! \note This method is only used to give the information of the current
    //! values to %OLAS. This method does NOT assemble them into the matrix.
    //! This is done by the method BuildInDirichlet.
    template<typename T>
    void SetDirichlet( const FeFctIdType fctId,
                       Integer eqnNr, const T &val );

    template<typename T>
    void SetDirichletMH( const FeFctIdType fctId, Integer eqnNr, const T &val , UInt &harmInt);

    //! Assemble the previously defined Dirichlet values into the algsys

    //! This method triggers the assembling of the dirichlet values into
    //! the algebraic system. The values itself must have been specified
    //! before via the SetDirichlet() method.
    //! \note After assembling the dirichlet values, the preconditioner and
    //! the solver have to be set up again.
    void BuildInDirichlet();
    
    template <typename T>
    void GetRidOfZeros(Double tol);

    void SetSysMatCopy( SBM_Matrix* actMat );

    SBM_Matrix* GetSysMatCopy(){return sysMat_[BACKUP];};

    void RestoreSysMat();

    //! correct RHS according to inhomogeneous Dirichlet bcs
    void AddIDBCToRHS(bool deltaIDBC = false);

    void SetOldDirichletValues();

    //! Return complete solution vector

    
    //! This method returns the complete current solution as 
    //! SBM-vector. The blocks are split per FctIdType,
    //! i.e. every FeFunction solution resides in its own SBM-block.
    //! \param sbmSolVec solution vector for all fctIds (FeFctId can be used
    //!                  as SBM-index)
    //! \param setIDBC true: considers inhomog. Dirichlet b.c. (standard)
    //!                false: is needed in case of nonlinear PDE, incremental formulation!
    void GetSolutionVal( SBM_Vector& sbmSolVec, bool setIDBC=true, bool deltaIDBC=false );
    

    void GetSolutionVal( const UInt& h, SBM_Vector& sbmSolVec, bool setIDBC=true, bool deltaIDBC=false );


    //! Return solution vector for one single FeFct

    //! This method returns the solution vector associated with one specific
    //! FeFunction. The entries are numbered according to the equations numbers
    //! of this FeFunction.
    //! \param solVec solution vector for specified FeFcunction
    //! \param fctId identifier for function related to sub-graph
    void GetSolutionVal( SingleVector& solVec,
                         const FeFctIdType fctId,
                         bool setIDBC, bool deltaIDBC=false  );

    //! Return solution vector for a given block in multiharmonic analysis

    //! This method returns the solution vector associated with one specific
    //! SBM-block. The entries are numbered according to the equations numbers
    //! \param solVec solution vector for specified FeFcunction
    //! \param ident (has no function) only purpose is to call the correct method
    void GetSolutionVal( SingleVector& ptSol,
                        const UInt& block,
                        bool setIDBC,
                        bool deltaIDBC,
                        const bool ident);

    //! Return solution vector for a given block in multiharmonic analysis

    //! This method returns the solution vector associated with one specific
    //! SBM-block. The entries are numbered according to the equations numbers
    //! \param solVec solution vector for specified FeFcunction
    void GetFullMultiHarmSolutionVal(SBM_Vector& solVec, bool setIDBC, bool deltaIDBC = false);



    //! Helper function introduced for non-linear-transient stepping
    //! Background:
    //!   During transient, non-linear stepping, the first solution step
    //!   always should call GetSolutionVal with setIDBC=true in
    //!   order to obtain the current IDBC values in the solution vector.
    //!   As long as we obtain a full solution of the system, where we
    //!   can completely overwrite the previous solution, this is no problem.
    //!   In case that we only want to compute and update to the solution of
    //!   the previous time step, we cannot simply add the returned solution
    //!   vector from GetSolutionVal, as this would also add up the IDBC nodes
    //!   instead of replacing them.
    //! Idea:
    //!   Before adding the updated solution that features the correct IDBC value of
    //!   the current time step, we clear the old IDBC values from the solution of
    //!   the last time step.
    //!   -> Function goes over solVec and sets all entries which correspond to
    //!       to IDBC nodes to 0.
    void ClearIDBCFromSolutionVal( SBM_Vector& solVec );
    void ClearIDBCFromSolutionVal( SingleVector& solVec,const FeFctIdType fctId);

    //! Return complete right-hand-side (RHS) vector

    //! This method returns the complete right-hand-side as 
    //! SBM-vector. The blocks are split per FctIdType,
    //! i.e. every FeFunction solution resides in its own SBM-block.
    //! \param sbmRhsVec RHS vector for all fctIds (FeFctId can be used
    //!                  as SBM-index)
    void GetRHSVal( SBM_Vector& sbmRhsVec );
    
    //! Return block right-hand-side (RHS) vector from multiharmonic anlysis

    //! This method returns the specified right-hand-side block (h) as
    //! SBM-vector
    void GetRHSVal( const UInt& h, SBM_Vector& sbmRhsVec );


    //! Return full right-hand-side (RHS) vector from multiharmonic anlysis

    //! This method returns the full multiharmonic right-hand-side
    void GetFullMultiHarmRHSVal(SBM_Vector& rhsVec );

    //! Return right-hand-side (RHS) vector for one single FeFct

    //! This method returns the RHS vector associated with one specific
    //! FeFunction. The entries are numbered according to the equations numbers
    //! of this FeFunction.
    //! \param rhsVec RHS vector for specified FeFcunction
    //! \param fctId identifier for function related to sub-graph
    void GetRHSVal( SingleVector &rhsVec,
                    const FeFctIdType fctId );

    //! Return right-hand-side (RHS) vector for one SBM Block in multiharmonic analysis

    //! This method returns the RHS vector associated with one specific
    //! SBM Block. The entries are numbered according to the equations numbers
    //! of the blockInfo_.
    //! \param rhsVec RHS vector for specified FeFcunction
    //! \param blockVec SBM block number
    //! \param ident just an identifier, not to mix it up with a call to the above
    //!              method, where the FeFctIdType is an argument (which is also an Integer)
    void GetRHSVal( SingleVector &rhsVec,
                    const UInt& blockVec,
                    const bool ident);


    // ***********************************************************************
    //   Methods for passing general data
    // ***********************************************************************


    //@{
    //! \name General Data Management

    //! Set type of needed matrix (STIFFNESS, MASS, ...)

    //! This method must be called to inform the algebraic system that a
    //! certain Finite-Element matrix (STIFFNESS, MASS, DAMPING, ...) is
    //! required in the simulation. In the case of an SBM_System we use
    //! the supplied Fct identifiers to determine which sub-matrices of
    //! the specified Finite-Element matrix will be needed.
    //! Thi method gets called once per bilinearform.
    //! \note The final problem matrix of type SYSTEM is always generated
    //!       by default.
    //! \param matrixType type of finite element matrix (STIFFNESS, MASS, ...)
    //! \param isSymmetric true, if matrix is symmetric (only possible for 
    //!                    fctId1 = fctid2).
    //! \param isComplex true, if matrix entries are complex
    //! \param fctId1  unique identifier for a Fct registered with the
    //!                graph manager
    //! \param fctId2  unique identifier for a Fct registered with the
    //!                graph manager
    void SetFEMatrixType( const FEMatrixType matrixType,
                          const bool isSymmetric,
                          const bool isComplex,
                          const FeFctIdType fctId1,
                          const FeFctIdType fctId2 = NO_FCT_ID );

    //! Get types of used matrices (STIFFNESS, MASS, ...)

    //! Get the types of the matrices, which were created by OLAS,
    //! for example STIFFNESS, MASS...
    //! \param matrixTypes set of enums of global matrices
    void GetFEMatrixTypes( std::set<FEMatrixType> &matrixTypes ) const;

    /** This is meant mainly for debug purpose */
//    SBM_Matrix* GetMatrix(FEMatrixType type) { assert(sysMat_.find(type) != sysMat_.end()); return sysMat_.find(type)->second; }
    SBM_Matrix* GetMatrix(FEMatrixType type) { assert(sysMat_.find(type) != sysMat_.end()); return sysMat_[type]; }
    
    SBM_Matrix* GetKeff(){return effMat_;};

    void PrintKeff();
    void PrintMatrixPart(FEMatrixType matrixPart);

    //! Return, if a non-zero static condensation block is present
    bool UseStaticCondensation() {
      return statCond_;
    }
    //@}



    // ***********************************************************************
    //   Methods for Algebraic Multigrid (AMG)
    // ***********************************************************************

    //@{
    //! \name Geometry Inclusion Routines

    //! Struct for one edge, it contains the edge-nodes (correct order)
    //! and its length
    struct EdgeGeom{
      StdVector<Integer> eNodes; // edge nodes
      Double length;  // length of edge
    };

    //! Set geometry information - nodal version

    //! Create a mapping between indices in the system matrix and
    //! coordinates of dof-positions (nodes)
    //! \note Fills geomInd_ variable
    //! \param indGeomMap Vector of coordinates of the specific points
    //! \param dim Dimension of the problem
    void SetGeomIndexMap(const StdVector< Vector<Double> >& indGeomMap,
                         const UInt& dim);

    //! Set geometry information - edge version

    //! Create a mapping between indices in the system matrix and
    //! edge-nodes respectively edge-lengths
    //! \note Fills geomIndEdge_ variable
    //! \param lengths Map where system matrix-index is the key and the length of the
    //!                corresponding edge is the value
    //! \param eNodes  Map where system matrix-index is the key and the node-numbers
    //!                of the corresponding edge are the values
    void SetEdgeIndexMap(boost::unordered_map<Integer, Double>& lengths,
                         boost::unordered_map< Integer, StdVector<Integer> >& eNodes);

    //! Build auxiliary matrix for algebraic multigrid solver/preconditioner

    //! In this method, we only differentiate between nodal- (Lagrangian-)
    //! and edge- (Nédéléc-) version
    void BuildAMGAuxMatrix();

    //! Build the auxiliary matrix for the nodal- (Lagrangian-) version of AMG
    void BuildAMGLagrangeAuxMatrix();

    //! Build the auxiliary matrix for the edge- (Nédéléc-) version of AMG
    void BuildAMGEdgeAuxMatrix();

    AMGType GetAMGType(){ return amgType_;}

    //! Return if AMG is even used
    bool UseAMG(){
      return useAMG_;
    }

    //@}


    // ***********************************************************************
    //   Various methods
    // ***********************************************************************

    //@{
    //! \name Miscellaneous
//
//    //! Export matrix to a file in MatrixMarket Format
//
//    //! The method will export the matrix to an ascii file according to the
//    //! MatrixMarket specifications. For details of the specification see
//    //! http://math.nist.gov/MatrixMarket
//    //! \param matrixType specifies matrix to be exported
//    //! \param filename name of output file
//    //! \param comment  string to be inserted into file header (optional)
//    virtual void Export( FEMatrixType matrixType, char *filename,
//                         char *comment = NULL ) const;
//
//    //! Print the specified matrix into the .las-file
//
//    //! The method will print the specified in the .las-file. The format
//    //! depends on the storage type of the matrix, which is normally an
//    //! indexed format (for crs and scrs). This mainly used for debugging
//    //! purpose.
//    //! \param matrixType specifies matrix to be exported
//    virtual void Print( FEMatrixType matrixType ) const;

    //! Print information about registered Fcts

    //! The method prints the name of the Fcts, their associated identifiers
    //! and the number of unkonwns per Fct into a specified ostream.
    void PrintRegistrationInfo() const;
//
//    /** This method gives back global matrix. For almost all cases you want
//     * the SYSTEM matrix */
//    virtual StdMatrix* GetSysMat(FEMatrixType type = SYSTEM) const { return NULL;}
//
//    //! returns system matrix (if stored in scrs!) as a vector
//    //! containing all upper triangle entries
//    virtual void GetFullSystemMatrixAsVec(Complex* &sysVec) const {;};

    //! Removes the IDBC Information from algebraic system

    //! This method is needed in nonlinear computations
    //! while computing solution updates
    void RemoveIDBCInfoFromMatrix() const {;};
    //@}

    /** Handle export linear system at the different phases. Checks by itself what needs to be done if anything.
     * Shall be called for each phase, set each time exactly one parameter to true */
    void ExportLinSys(bool setup, bool pre_solve, bool post_solve);

    void ExportMHSys(int step);

    //! In multiharmonic analysis, set the nonzero sbm-blocks
    inline void SetNnzSBMInd(const StdVector<UInt>& sbmInd){ nnzSBMInd_ = sbmInd;};

    BaseEigenSolver* GetEigenSolver(){ return eigenSolver_; };

    BaseSolver* GetSolver(){ return solver_; };

    //! Return if it is a multiharmonic analysis
    bool IsMultHarm(){return isMultHarm_; };

    PtrParamNode GetExportLinSysParam();

    bool IsMatrixComplex(){return isMatrixComplex_;};

    /*
     * function for StdSolveStepHyst
     *
     * background: when solving the nonlinear system of equations via newtons method,
     *             we need an approximation of the Jacobian; this approximation is delivered
     *             from CoefFunctionHyst and is assembled to the stiffness matrix;
     *             for the computation of the residual, we need the linear stiffness matrix, i.e.
     *             without considering hysteresis; as this matrix should be constant, we need not
     *             to reassemble it each time; instead we could keep a deep copy of the system and
     *             reuse it during residual computation
     */
    void BackupCurrentSystemMatrix(FEMatrixType storageLocation);
    void LoadBackupToCurrentSystemMatrix(FEMatrixType storageLocation);

    // copys one matrix in the algsys to another storage location. useful for backups
    void CopyMatrixToOther(FEMatrixType matrix, FEMatrixType other, bool add = false);

  protected:

    //! Auxiliary method for logging information on matrix patterns

    //! This method can be used to print a summary of the Finite-Element
    //! matrices of the algebraic system and their sub-matrix pattern to
    //! a stream. It is mostly intended as an auxiliary method for
    //! CreateLinSys().
    void PrintFeMatrixInfo();

    //! Check consistency of setup (solver, preconditioner, matrix format)
    
    //! This method performs consistency checks on the type of solver,
    //! the preconditioner and the matrix to be created.
    //! It checks for example, if the symmetry type of the matrix matches 
    //! the solver and matrix specifications.
    //! \note This method basically performs the former CFSOLASParams 
    //!       functionality.
    void CheckConsistency();
    
    /** helper for CheckConsistency(), determines reordering setting based on
     * xml and compile time option. If set to default, compile settings are used
     * @param can_change true if obtained by "default" in xml or CFS_REORDERING and false in precise setting
     * @return metis, sloan or noreodering */
    BaseOrdering::ReorderingType ReorderingDefault(PtrParamNode pn, bool& can_change);

    //! Generate  SBM matrix according to graph information
    
    //! This method generates the SBM matrix for a given entry type. In
    //! addition, the flag sharePattern denotes, if  the matrix pattern
    //! can be shared.
    SBM_Matrix* GenerateSBM_Matrix( FEMatrixType matType,
                                    BaseMatrix::EntryType entryType,
                                    bool sharePattern );
    
    //! Pointer to the graph manager object
    GraphManager *graphManager_;

    //! Pointer to solver
    BaseSolver *solver_;

    //! Pointer to the preconditioner object
    BasePrecond *precond_;
    
    //! Pointer to eigenvalue solver
    BaseEigenSolver *eigenSolver_;
    
    //! Pointer to struct with solution strategy information
    shared_ptr<SolStrategy> solStrat_; 

    //@{ \name Attributes storing information on individual functions / blocks

    //! Number of FeFunctions

    //! Holds the number of FeFunctions which is equivalent to the number of
    //! sbm blocks in the algebraic system matrix.
    UInt numFcts_;

    //! Number of super block matrices
    UInt numBlocks_;

    //! Number of registered FeFunctions so far
    UInt numRegFcts_;
    
    //! Flag, if registration is already finished
    bool registrationFinished_;
    
    //! Flag, if system is already created
    bool systemCreated_;

    //! Total number of equations per FeFunction
    StdVector<UInt> numEqnsPerFct_;
    
    //! Last free equation number per fct
    StdVector<UInt> lastFreeEqnPerFct_;

    /** STL map for associating Fct-identifiers with FeFunction name
     * TODO: why is this a map and no a vector as we start with 0?! */
    std::map<FeFctIdType,std::string> fctNames_;
    
    //! Store for each FeFct, if the assoicated system matrix is symmetric
    std::map<FeFctIdType,bool> matIsSymm_;

    //@}

    //! Size of the matrix (i.e. sum of sizes of all blocks)

    //! Since we are talking of square matrices the size equals the order
    //! of the matrix, which corresponds to the number of rows resp. columns.
    UInt size_;
    
    // =======================================================================
    // BLOCK HANDLING / RE-NUMBERING
    // =======================================================================

    //! Map (fctId,eqnNr) to SBM matrix block (only diagonals)
    
    //! This map is used to translate the CFS-oriented view (fctId, eqnNr)
    //! into the OLAS-oriented view, i.e. SuperMatrixBlocks (SBMBlocks).
    //! For each (fctId,eqnNr) it stores the destination SBMBlock. To get
    //! the corresponding row/col index, one has to use the #eqnToIndex array
    //! of the corresponding SBMBlockInfo, stored in #blockInfo_.
    StdVector<StdVector<UInt> > eqnToSBMBlock_; 
      
    //! Store for each SBMBlock the block information
    
    //! This structure stores for each SBMBlock the definition, i.e.
    //!  - a map (fctId,eqnNr) to (index)
    //!  - a definition of (optional) subblocks (e.g. for block preconditioners)
    //!  - a flag, if the block can be reordered
    StdVector<GraphManager::SBMBlockInfo*> blockInfo_;
    
    //! Store for each fctId in which sbmBlocks it occurs
    
    //! When interfacing with CFS, the algebraic systems gets requests on a
    //! fctId-Level (GetSolution(), InitMatrix() ...). Inside OLAS we have
    //! to know, which blocks are affected by this operation. Thus
    //! this structure gets filled during the definition of SBM blocks.
    std::map<FeFctIdType, std::set<UInt> > fctIdsInBlocks_;
    
    // =======================================================================
    // BOUNDARY CONDITIONS / CONSTRAINTS
    // =======================================================================

    //@{ \name Handling of Boundary Conditions and Constraints

    //! Number of inhomogeneous Dirichlet per feFunction 
    StdVector<UInt> numDirichletValuesPerBlock_;
    
    //! Object for handling inhomogeneous Dirichlet boundary conditions
    BaseIDBC_Handler *idbcHandler_;

    //! Status of system matrix w.r.t. penalty terms

    //! This flag is used to keep track on whether we must insert penalty
    //! terms into the system matrix or not. The latter can be the case if
    //! either the penalty formulation is not employed, or the insertion
    //! was already performed.
    bool assembleDirichletToSysMat_;

    //@}

    //! Binary predicate used as comparison operator for SubMatrixID sets

    //! This class provides a binary function predicate that is used as
    //! comparison operator in STL sets containing elements of the
    //! SubMatrixID type. The predicate returns true for the comparison
    //! sID1 < sID2, if either the row index of sID1 is smaller than that
    //! of sID2, or if both row indices are equal and the column index of
    //! sID1 is smaller than that of sID2.
    //! \note For the above definition of the predicate equivalance of
    //! sID1 and sID2, i.e. sID1 < sID2 false and sID2 < sID1 false,
    //! is the same as identity.
    class SortSubMatrixID {
    public:
      bool operator() ( const SubMatrixID &sID1,
                        const SubMatrixID &sID2 ) const {
        bool retVal = false;
        if ( sID1.rowInd < sID2.rowInd ) {
          retVal = true;
        }
        else if ( sID1.rowInd == sID2.rowInd ) {
          retVal = sID1.colInd < sID2.colInd;
        }
        return retVal;
      }
    };
    
    //! Helper definition for set of sub-matrices
    typedef std::set<SubMatrixID,SortSubMatrixID> SubMatrixSet;
    
    //! Vector with IDs of sub-matrices of Finite-Element matrices by fctIds

    //! This map contains for each possible Finite-Element matrix managed
    //! by the algebraic system object a set containing the index tuples for
    //! the non-zero sub-matrices of the Finite-Element matrix. These sets
    //! are filled during the setup phase by calls to SetFEMatrixType() and
    //! used in the CreateLinSys() method for generation of the SBM_Matrix
    //! objects.
    std::map<FEMatrixType,SubMatrixSet> feSubMatricesByFctId_;
    
    //! Vector with IDs of sub-matrices of Finite-Element matrices by numBlocks
    std::map<FEMatrixType, SubMatrixSet> feSubMatricesByBlocks_;
    
    //! Set defining which type of matrices (stiffness, mass,...) are used
    std::set<FEMatrixType> matrixTypes_;

    //! Pointer to <system> ParamNode
    PtrParamNode myParam_;
    
    //! Pointer to information ParamNode
    PtrParamNode myInfo_;
    
    //! Flag indicating use of idbc elimination via Penalty approach
    bool usingPenalty_;
    
    //! Flag indicating use of static condensation
    bool statCond_;
    
    //! Flag indicating use of multiharmonic analysis
    bool isMultHarm_;

    StdVector<UInt> nnzSBMInd_;

    //! Flag indicating, if system matrix is complex
    bool isMatrixComplex_;

    //! Flag indicating, if solution vector is complex
    bool isSolutionComplex_;
  
  // =======================================================================
  // MATRICES / VECTORS
  // =======================================================================

    //@{ \name Matrices and Vector
    //! Pointers to Finite-Element matrices

    //! This attribute points to an array containing pointers to the different
    //! Finite-Element matrices (SYSTEM, MASS, STIFFNESS, ... ) managed by the
    //! algebraic system. In the case of the %SBM_System class these are of
    //! course SBM_Matrix objects.
    std::map<FEMatrixType, SBM_Matrix*> sysMat_;

    //! Effective system matrix (without condensation block)

    //! This is the "effective" matrix to be solved, i.e. the one passed to
    //! the solver. This a shallow copy of the #sysMat and thus contains only
    //! pointers into the global system matrix.
    //! In case we use static condensation, the effective matrix does not
    //! contain the interior block and the related coupling matrices.
    SBM_Matrix *effMat_;

    //! Vector containing right-hand side of the linear system (see #effMat_)
    SBM_Vector *rhs_;

    //! Effective rhs vector
    SBM_Vector *effRhs_;

    //! temporary rhs
    SBM_Vector* tmpRHS_;

    //! Vector containing (approximate) solution of the linear system
    SBM_Vector *sol_;

    //! Effective solution vector
    SBM_Vector *effSol_;

    //! Buffer for storing the eigenvalues of the system
    SingleVector *eigenValues_;

    //! Buffer for storing the error bounds of the eigenvalues
    SingleVector *eigenValError_;

    //! Store for each diagonal SBM-Block if it is symmetric.
    StdVector<bool> isDiagBlockSymm_;

    //! Pointer to a pool for sharing sparsity patterns between matrices
    PatternPool *patternPool_;

    //! Store for each SBM block the pattern Id (if applicable)
    std::map<SubMatrixID, PatternIdType, SortSubMatrixID> sbmPatternIds_;

    //! Flag, if matrix pattern can be shared
    bool sharedPatternPossible_;

    //! Flag signaling symmetry of system matrix
    bool sbmSymm_;

    //! Flag signaling use of only one matrix block

    //! This flag denotes, if the complete system matrix just consists of one
    //! single matrix block. In this case we can simply access the (1,1) entry
    //! as StdMatrix and have the complete system.
    bool onlyOneMatrixBlock_;

    //! Flag if we have distinct matrix graphs for different matric types
    bool distinctMatGraphs_;
    //@}


    // =======================================================================
    // AMG SECTION
    // =======================================================================

    //@{ \name Algebraic Multigrid variables

    //! Flag indicating use of multigrid methods
    bool useAMG_;

    //! Auxiliary matrix, needed for special (geometric-)algebraic multigrid methods
    StdMatrix *auxMatAMG_;

    //! Store coordinate of each index geomInd_[i] = coordinate of index i
    StdVector< Vector<Double> > geomInd_;

    //! Special structure for edge-AMG geomIndEdge_[i] contains
    //! the coordinates of the two nodes of edge i (already in correct
    //! orientation)
    boost::unordered_map< Integer, EdgeGeom> geomIndEdge_;

    //! Special structure for edge-AMG geomIndEdge_[i] contains
    //! the indices in the auxiliary-matrix of the two nodes
    // of edge i (already in correct orientation)
    StdVector< StdVector< Integer> > edgeIndNode_;

    //! Map for nodeNum <-> index in auxiliary matrix
    boost::unordered_map< Integer, UInt> indexNodeNum_;

    //! vector for index <-> nodeNum in auxiliary matrix
    StdVector<Integer> nodeNumIndex_;

    //! only valid for AMG: dimension of singlefield-problem
    UInt dim_;

    //! Flag, needed for the specialized AMG-methods
    AMGType  amgType_;

    //! Switch if edge element discretization is used
    bool edge_;
    //@}


    // =======================================================================
    // CACHE SECTION
    // =======================================================================

    //@{ \name Cached vectors / matrices for fast access
    //! Index vectors for element matrix assembly
    CfsTLS<StdVector< StdVector<UInt> > > rowIndList1_;
    CfsTLS<StdVector< StdVector<UInt> > > rowList1_;
    CfsTLS<StdVector< StdVector<UInt> > > rowIndList2_;
    CfsTLS<StdVector< StdVector<UInt> > > rowList2_;
    CfsTLS<StdVector< StdVector<UInt> > > colIndList1_;
    CfsTLS<StdVector< StdVector<UInt> > > colList1_;
    CfsTLS<StdVector< StdVector<UInt> > > colIndList2_;
    CfsTLS<StdVector< StdVector<UInt> > > colList2_;

    //! Index vector for element position
    CfsTLS<StdVector<UInt> > rowBlocks_, colBlocks_, rowNums_, colNums_;
    //@}
  
  };

  //extern template void AlgebraicSys::GetRidOfZeros<double>(double tol);

} // namespace


#endif // OLAS_ALGEBRAIC_SYS_HH
