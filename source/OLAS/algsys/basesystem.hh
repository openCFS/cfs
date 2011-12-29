// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef OLAS_BASESYSTEM_HH
#define OLAS_BASESYSTEM_HH

#include <stddef.h>
#include <iosfwd>
#include <map>
#include <set>
#include <string>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/defs.hh"
#include "General/environment.hh"

namespace CoupledField {

  class BaseEigenSolver;
  class BaseGraphManager;
  class BaseIDBC_Handler;
  // Forward declarations
  class BaseSolver;
  class BaseVector;
  class SingleVector;
  class StdMatrix;
  struct BaseEntryManipulator;
template <class TYPE> class Matrix;
  template <class TYPE> class StdVector;
  template <class TYPE> class Vector;

  //! Base class for linear algebraic system
  //! Only used as interface.

  //! The BaseSystem object (resp. its derived classes) reads the following
  //! parameters from the myParams_ object
  //! - MatrixEntryType = Type of entries of matrices (DOUBLE,COMPLEX)
  //! - MatrixStorageType = Storage type of matrices (SPARSE_NONSYM,
  //!   SPARSE_SYM, ...)
  //! - NumDof = Degrees of freedoms per node
  //! - MaxIter = maximal number of iterations
  //! - %Parallel = should parallel matrices/vectors/algorithms be used
  //!
  //! and writes the following keys into the myReport_ object
  //! - numIter = number of iterations performed
  //! - finalPrecondResNorm = the Euclidean norm of the final preconditoned
  //!   residual
  class BaseSystem {

  public:

    //! Default Constructor
    BaseSystem(PtrParamNode pn = PtrParamNode());

    //! Default Destructor
    virtual ~BaseSystem();

    /**  killme! check if really used in SIMP optimization */
    BaseEntryManipulator* GetEntryManipulator()
    {
      return assemble_;
    }

    //! Obtain the dimension of the linear system

    //! This method can be used to query the dimension of the linear system.
    //! By dimension we understand the length of the vectors  \f$ x \f$ and
    //!   \f$ b  \f$  in the linear system    \f$  Ax=b  \f$  and the number of rows and
    //! columns in the system matrix    \f$  A  \f$ .
    //!
    //! \note
    //! - The number is understood on a meta level, i.e. the number of
    //!   entries on the scalar level is larger when dof-blocking is employed.
    //! - The dimension coincides with the highest equation number that
    //!   still belongs to a free degree of freedom.
    UInt GetDimension() {
      return (UInt)size_;
    }

    /** Obtain our IDBC-Handler */
    BaseIDBC_Handler* GetIdbcHandler() {
      return idbcHandler_;
     }

    // ***********************************************************************
    //   Methods for creating and initalizing the algebraic system
    // ***********************************************************************

    //@{
    //! \name Creation and Initialization

    //! Generate objects for the linear system(s)

    //! The method will generate the matrix and vector objects required for
    //! the different system matrices and the solution and right hand side
    //! vector.
    virtual void CreateLinSys() = 0;

    //! Generate preconditioner
    virtual void CreatePrecond() = 0;

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
    virtual void CreateSolver() = 0;

    //! Generate EigenSolver object

    //! Calling this method triggers the generation of an EigenSolver object
    //! which ist used to calculate a generalized eigenvalue problem of the
    //! form M(du^2/dt^2)=Ku or of the quadratic eigenvalue problem, if
    //! a damping matrix is present.
    //! Before this object can be used the SetupEigenSolver() method should
    //! be called.
    //! \note If an Eigenfrequency analysis is performed, the methods
    //! SetupPrecond() and SetupSolver() must not be called!
    virtual void CreateEigenSolver() = 0;

    /** Trigger setup of preconditioner.
     * Calling this method will trigger the setup phase of the
     * preconditioner. The setup is performed using the system matrix of the
     * linear system.
     * @param note This method must not be called if an eigenfrequency analysis
     * is performed, since this method creates only a preconditioner
     * to solver which solves a system Ax=b.
     * @param analysis_id contains the reference to the analysis step for info.xml output */
    virtual void SetupPrecond() = 0;

    /** Trigger setup of solution method.
     * Calling this method will trigger the setup phase of the solver. This
     * is especially important for direct solvers, where typically the
     * factorisation of the problem matrix will be performed at this stage.
     * The setup is performed using the system matrix of the linear system.
     * @param analysis_id @see SetupPrecond() */
    virtual void SetupSolver(PtrParamNode analysis_id) = 0;

    //! Trigger setup of eigenvalue solver

    //! Calling this method will trigger the setup phase of the eigensolver.
    //! The setup is performed using the system matrix of the linear system.
    //! \param numFreq Number of eigenfrequencies to be calculated
    //! \param shift Frequency shift applied to the system
    //! \param quadratic Flag indicating if a quadratic eigenvalue problem
    //! (true) or a generalized problem (false) is to be solved
    virtual void SetupEigenSolver( UInt numFreq, Double shift,
                                   bool quadratic ) = 0;


    /** Solve the linear system.
     * Calling this method triggers the solution of the linear system
     * defined by the Finite-Element system matrix, i.e. sysMat_[SYSTEM],
     * and the given right-hand side vector, i.e. #rhs_. The computed
     * (approximat) solution is stored in #sol_.
     * The method used for solving the system depends on the solver and
     * preconditioner constructed and the parameters in myParams_.
     * For iterative solvers an initial guess is created by inserting the
     * Dirichlet values in the correct positions of the solution vector in
     * case of the penalty formulation.
     * @note This method must not be called if an eigenfrequency analysis
     * is performed, since this method is only used to solve a system of the
     * form Ax=b.
     * @param analysis_id identifies the analysis step.
     *        When the linear system is exported via file, the comment is used. */
    virtual void Solve(PtrParamNode analysis_id) = 0;

    //! Calculate eigenfrequencies of a generalized eigenvalue problem

    //! Calling this method triggers the solution of a generalized eigenvalue
    //! problem, i.e. \f[ \mathbf{M} \ddot{x} = \mathbf{K} x \f].
    //! The number of eigenfrequencies is specified by setting the parameter
    //! "numEigenFreqs" of the OLAS_Params object.
    //! The resulting eigenfrequencies are stored into an array, as well as the
    //! resulting error.
    //! \param frequencies Read-only buffer which contains the eigenfrequencies
    //!                    of the generalized eigenvalue problem
    //! \param err Reached error norm for each eigenvalue
    //! \return Number of converged eigenvalues
    //! \note This method may only be calaled if SetupEigenfrequencySolver()
    //!       was called previously.
    virtual void CalcEigenFrequencies( Vector<Double>& frequencies,
                                       Vector<Double>& err ) = 0;
 
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
    virtual void CalcEigenFrequencies( Vector<Complex>& frequencies,
                                       Vector<Double>& err ) = 0;

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
    virtual void CalcEigenMode( UInt numMode ) = 0;

    //@}

    // ***********************************************************************
    //   Methods for creating and assembling the matrix graph
    // ***********************************************************************

    //@{
    //! \name Matrix Graph Methods

    //! Prepares graph manager for generation of sub-graphs

    //! This method must be called after construction of the graph manager
    //! object and before any other method. It allows the graph manager to
    //! prepare itself for the registration and generation of the sub-graphs.
    //! \param numPDEs number of PDE objects whose connectivity graphs the
    //!                graph manager will administer.
    void GraphSetupInit( UInt numPDEs );

    //! Finalises setup of the graph manager

    //! This method must be called after the assembly of all sub-graphs was
    //! done, i.e. once all PDEs a coupling objects have conveyed their
    //! connectivity information to the graph manager.
    void GraphSetupDone();

    //! Obtain a unique PDE identifier

    //! Return a unique PDE identifier for each PDE, which will be used
    //! for all subsequent communication with the algebraic system.
    //! \param pdeType     string giving the pde type, e.g. electrostatic
    //! \return unique identifier for the PDE
    PdeIdType ObtainPDEId( const std::string &pdeType );

    //! Register PDE with graph manager

    //! A PDE must call this method to register itself with the PDE manager.
    //! In the case of CFS++ this method is not called
    //! directly, but via the algebraic system interface.

    //! \param pdeId       unique identifier for the PDE
    //! \param numEqns     number of unknowns in the linear system associated
    //!                    with the PDE (i.e. number of equation numbers from
    //!                    CFS++)
    //! \param numLastFreeDof specifies the equation number of the last
    //!                       unfixed degree of freedom (if inhomogeneous
    //!                       Dirichlet boundary conditions are treated by the
    //!                       penalty method numLastFreeDof
    //!                          \f$  \stackrel{!}{=}  \f$  numEqns)
    void RegisterPDE( const PdeIdType pdeId, UInt const numEqns,
                      UInt const numLastFreeDof );

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
    //! \param eqnNrs1       array containing equation numbers for the
    //!                       degrees of freedom attachted to the element
    //!                       (unknowns at nodes or edges, whatsoever), must
    //!                       be a one based pointer; for first PDE
    //! \param elemSize1      length of the eqnNrs array; for first PDE
    //! \param identifierPDE2 identifier for second PDE related to sub-graph
    //! \param eqnNrs2       array containing equation numbers for the
    //!                       degrees of freedom attachted to the element
    //!                       (unknowns at nodes or edges, whatsoever), must
    //!                       be a one based pointer; for second PDE
    //! \param elemSize2      length of the eqnNrs array; for second PDE
    //! \param setCounterPart flag signalling whether the counterpart, i.e.
    //!                       the connectivity with the PDE identifiers and
    //!                       equation numbers reversed, should also be
    //!                       inserted into the graph
    void SetElementPos( const PdeIdType identifierPDE1,
                        const StdVector<Integer>& eqnNrs1,
                        const PdeIdType identifierPDE2,
                        const StdVector<Integer>& eqnNrs2,
                        bool setCounterPart );

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
    void GetReordering( const PdeIdType identifier,
                        StdVector<UInt>& newOrder );

    //@}

    // ***********************************************************************
    //   Methods for assembling the matrices, rhs and bc
    // ***********************************************************************

    //@{
    //! \name Assembling Routines

    //! Initialize specified matrix with zeros

    //! Calling this method initialises the Finite-Element matrix designated
    //! by the provided FEMatrixType identifier with zeros.
    //! After initialisation all matrix entries in the matrix' sparsity
    //! pattern are zero. If no FEMatrixType identifier is specified all
    //! matrices in the internal set (see #matrixTypes_) will be initialised.
    //! In the case of an SBM_System the PdeIDType identifier can be used to
    //! initialise only a sub-matrix of the Finite Element matrix. If no
    //! PdeIDType identifier is specified the complete matrix will be
    //! initialised.
    //! \param matrixID      type of Finite Element matrix
    //!                      (STIFFNESS, MASS, ...)
    //! \param identifierPDE unique identifier obtained from the ObtainPDEId()
    //!                      method
    virtual void InitMatrix( FEMatrixType matrixID = NOTYPE,
                             const PdeIdType identifierPDE = NO_PDE_ID ) = 0;

    //! Set right hand side vector to zero

    //! This method sets the right hand side vector associated with the
    //! specified PDE to zero. If no PDE identifier is given, the complete
    //! right-hand side vector is zeroed.
    //! \note In the case of a StandardSystem the PDE identifier is ignored.
    //!       We currently do not support setting only the part of a SingleVector
    //!       related to a single PDE to zero in this case.
    //! \param identifierPDE unique identifier obtained from the ObtainPDEId()
    //!                      method
    virtual void InitRHS( const PdeIdType identifierPDE = NO_PDE_ID ) = 0;

    //! Set the global rhs to given values

    //! Set the right hand side vector to the values specified in the array
    //! newRHS
    //! \param newRHS Pointer to new right-hand-side values
    //! \note The values of newRHS are copied, so the pointer to newRHS
    //! can be changed afterwards!
    virtual void InitRHS( const BaseVector& newRHS ) = 0;

    //! Set the solution vector of the specified PDE to zero

    //! This method sets the part of the solution vector associated with the
    //! specified PDE to zero. If no PDE identifier is given, the complete
    //! solution vector is zeroed.
    //! \note In the case of a StandardSystem the PDE identifier is ignored.
    //!       We currently do not support setting only the part of a SingleVector
    //!       related to a single PDE to zero.
    //! \param identifierPDE unique identifier obtained from the ObtainPDEId()
    //!                      method
    virtual void InitSol( const PdeIdType identifierPDE = NO_PDE_ID ) = 0;

    //! Set the global solution vector to given initial values

    //! Set the solution vector to the values specified in the array
    //! newSol
    //! \param newSol Pointer to new solution values
    //! \note The values of newSol are copied, so the pointer to newSol
    //! can be changed afterwards!
    virtual void InitSol( const BaseVector& newSol ) = 0;


    //@{
    //! Assemble an element matrix into the global one

    //! This methods assembles the given element matrix into a specified
    //! global one (MASS, STIFFNESS, etc.), i.e. adds the entries to the.
    //! global matrix. For element matrices which are associated
    //! only with one PDE, only one identifier, the eqn numbers of the matrix
    //! and the number of eqnNrs have to be specified.
    //! For matrices which are associated with two different pde identifiers,
    //! the matrix will be assembled into the upper off-diagonal matrix-block
    //! and the lower transposed one.
    //! \param matrix_id type of finite element destination matrix
    //!                  (STIFFNESS, MASS, ...)
    //! \param elemmat entries of the element matrix in sequential array
    //! \param identifierPDE1 identifier for first PDE related to sub-graph
    //! \param eqnNrs1 equation numbers (1-based) of the element matrix
    //!                w.r.t. sub-graph associated with identifierPDE1
    //! \param numEqn1 number of equations related to sub-graph of
    //!                identifierPDE1
    //! \param identifierPDE2 identifier for first PDE related to sub-graph
    //! \param eqnNrs2 equation numbers (1-based) of the element matrix
    //!                w.r.t. sub-graph associated with identifierPDE2
    //! \param numEqn2 number of equations related to sub-graph of
    //!                identifierPDE2
    //! \param setCounterPart if this flag is true, then the method will
    //!                not only insert the element matrix    \f$  E  \f$ , but also
    //!                its transpose    \f$  E^T  \f$ . In doing so also the
    //!                row and column indices derived from the equation
    //!                numbers are interchanged. Note that this is only
    //!                supported for off-diagonal blocks, i.e. for cases
    //!                with different PDE identifiers.
    virtual void SetElementMatrix( FEMatrixType matrix_id, 
                                   const Matrix<Double>& elemmat,
                                   PdeIdType identifierPDE1,
                                   const StdVector<Integer>& eqnNrs1,
                                   PdeIdType identifierPDE2,
                                   const StdVector<Integer>& eqnNrs2,
                                   bool setCounterPart ) = 0;

    virtual void SetElementMatrix( FEMatrixType matrix_id, 
                                   const Matrix<Complex>& elemmat,
                                   PdeIdType identifierPDE1,
                                   const StdVector<Integer>& eqnNrs1,
                                   PdeIdType identifierPDE2,
                                   const StdVector<Integer>& eqnNrs2,
                                   bool setCounterPart ) = 0;
    //@}

    //@{
    //! Assemble the local rhs vector to the global one

    //! This method adds the entries of the element right hand side to the
    //! global ones, according the pdeIdentifier and the associated equation
    //! numbers.
    //! \param elemRHS entries of the element right hand side in a sequential
    //!                array
    //! \param idPDE   identifier of the PDE related to sub-graph
    //! \param eqnNrs  equation numbers (1-based) of the element rhs
    //!                w.r.t. sub-graph associated with idPDE
    //! \param numEqn  length of eqnNrs array
    virtual void SetElementRHS( const Vector<Double>& elemRHS, 
                                const PdeIdType idPDE,
                                StdVector<Integer>& eqnNrs ) = 0;

    virtual void SetElementRHS( const Vector<Complex>& elemRHS, 
                                const PdeIdType idPDE,
                                StdVector<Integer>& eqnNrs ) = 0;
    //@}


    //@{
    //! Adds a value to a given global rhs entry

    //! This method adds a given value to the global right hand side vector.
    //! Therefore a pde identifier, the related equation number and a
    //! degree of freedom has to be specified.
    //! \param val value to be added
    //! \param identifierPDE identifier of the PDE related to sub-graph
    //! \param eqnNr equation number of the node to be set
    virtual void SetNodeRHS(Double val, PdeIdType identifierPDE,
                            Integer eqnNr) = 0;

    virtual void SetNodeRHS(Complex val, PdeIdType identifierPDE,
                            Integer eqnNr) = 0;
    //@}

    //! Performs a matrix-vector multiplication and ads the vector to the rhs

    //! This method multiplies a specified global matrix (STIFFNES, MASS,
    //! etc.) with a given vector and adds the result vector to the global
    //! rhs.
    //! \f[ rhs = rhs + \mathbf A_{matrix_id} \cdot fup \f]
    //! This method is mainly used in a time step scheme to adapt the rhs.
    //! \param matrix_id type of finite element matrix (STIFFNESS, MASS, ...)
    //!                  which gets multiplied
    //! \param fup array with vector entries, which get multiplied
    virtual void UpdateRHS(FEMatrixType matrix_id, const BaseVector& fup) = 0;

    //! Add a value to a diagonal matrix entry

    //! This method allows to directly add a value to the value of an entry
    //! on the main diagonal of a certain of the different problem matrices.
    //! \param matrixID specifies which matrix is to be manipulated
    //! \param pdeID    identifier for PDE related to sub-graph
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
    virtual void AddToDiagMatrixEntry( FEMatrixType matrixID,
                                       const PdeIdType pdeID,
                                       Integer eqnNum,
                                       Double val ) = 0;
    virtual void AddToDiagMatrixEntry( FEMatrixType matrixID,
                                       const PdeIdType pdeID,
                                       Integer eqnNum,
                                       Complex val ) = 0;
    //@{
    //! Get value of specific matrix entry

    //! This method allows to directly get the value of a specific matrix
    //! entry of a given problem matrix. If the given index pair is not
    //! contained in the matrix, an error is thrown.
    //! \param matrixID  specifies which matrix is to be queried
    //! \param rowPdeID  identifier for PDE related to row sub-graph
    //! \param rowEqnNum equation number of the row index to get
    //! \param colPdeID  identifier for PDE related to column sub-graph
    //! \param colEqnNum equation number of the column index to get
    //! \param val       return value
    //! \note For sparse matrices this method may be very costly and slow,
    //!       as the given index pair has to be searched for and may not
    //!       be accessible directly!
    virtual void GetMatrixEntry( FEMatrixType matrixID,
                                 const PdeIdType rowPdeID,
                                 Integer eqnNum1,
                                 const PdeIdType colPdeID,
                                 Integer eqnNum2,
                                 Double & val ) = 0;

    virtual void GetMatrixEntry( FEMatrixType matrixID,
                                 const PdeIdType rowPdeID,
                                 Integer rowEqnNum,
                                 const PdeIdType colPdeID,
                                 Integer colEqnNum2,
                                 Complex & val ) = 0;
    //@}

    //@{
  //! Set value of specific matrix entry

    //! This method allows to directly set the value of a specific matrix
    //! entry of a given problem matrix. If the given index pair is not
    //! contained in the matrix, an error is thrown.
    //! \param matrixID  specifies which matrix is to be manipulated
    //! \param rowPdeID  identifier for PDE related to row sub-graph
    //! \param rowEqnNum equation number of the row index to be set
    //! \param colPdeID  identifier for PDE related to column sub-graph
    //! \param colEqnNum equation number of the column index to be set
    //! \param val       value to be set
    //! \param setCounterPart if this flag is true, then the method will
    //!                not only set the value on position (eqnNum1, eqnNum2)
    //!                but also on its transposed position (eqnNum2, eqnNum1).
    //!                If the entry is located on an off-diagonal entry (i.e.
    //!                different pdeIDs) the counterpart in the related
    //!                transposed block is inserted.
    //! \note For sparse matrices this method may be very costly and slow,
    //!       as the given index pair has to be searched for and may not
    //!       be accessible directly!
    virtual void SetMatrixEntry( FEMatrixType matrixID,
                                 const PdeIdType rowPdeID,
                                 Integer rowEqnNum,
                                 const PdeIdType colPdeID,
                                 Integer colEqnNum,
                                 Double val,
                                 bool setCounterPart ) = 0;

    virtual void SetMatrixEntry( FEMatrixType matrixID,
                                 const PdeIdType rowPdeID,
                                 Integer rowEqnNum,
                                 const PdeIdType colPdeID,
                                 Integer colEqnNum,
                                 Complex val,
                                 bool setCounterPart ) = 0;
    //@}

    //! Constructs the effective system matrix

    //! Constructs the effective system matrix by multiplying the individual
    //! matrices (expect AUXILIARY!) with factors and adding them to one system matrix
    //! \param matFactors a map which contains for each matrixtype (DAMPING,
    //!                   STIFFNESS,...) the according factor
    virtual void ConstructEffectiveMatrix( const std::map<FEMatrixType,Double>
                                           &matFactors ) = 0;

    //@{
    //! Pass a Dirichlet value to %OLAS

    //! This method passes the value of a given inhomogeneous Dirichlet
    //! equation number to %OLAS.
    //! \param pdeID  identifier for PDE related to sub-graph
    //! \param eqnNr  equation number of the dirichlet restraint
    //! \param val    value of the dirichlet restraint
    //!
    //! \note This method is only used to give the information of the current
    //! values to %OLAS. This method does NOT assemble them into the matrix.
    //! This is done by the method BuildInDirichlet.
    virtual void SetDirichlet( const PdeIdType pdeID,
                               Integer eqnNr, const Double &val ) = 0;

    virtual void SetDirichlet( const PdeIdType pdeID,
                               Integer eqnNr, const Complex &val ) = 0;
    //@}

    //! Assemble the previously defined Dirichlet values into the algsys

    //! This method triggers the assembling of the dirichlet values into
    //! the algebraic system. The values itself must have been specified
    //! before via the SetDirichlet() method.
    //! \note After assembling the dirichlet values, the preconditioner and
    //! the solver have to be set up again.
    virtual void BuildInDirichlet() = 0;

    //! Return the pointer to the current solution of a PDE

    //! This method passes the current solution of one given PDE
    //! as a pointer to a buffer containing the solution entries.
    //! If no identifier is given, the pointer to the beginning of the
    //! complete solution array is passed.
    //!
    //! \param ptSol pointer (0-based) to the solution buffer
    //! \param identifierPDE identifier for PDE related to sub-graph
    //!
    //! \note The return buffer is guaranteed to retain the current solution
    //! until the next call of this method (after solving the next step)!
    virtual void GetSolutionVal( SingleVector& ptSol,
                                 const PdeIdType identifierPDE
                                 = NO_PDE_ID ) = 0;

    //! Return the pointer to the current rhs value of a PDE

    //! This method passes the current rhs as a pointer to buffer with
    //! right hand side entries.
    //!
    //! \param ptRhs pointer (0-based) to the buffer with rhs values
    //! \param identifierPDE identifier for PDE related to sub-graph
    //!
    //! \note The return buffer is guaranteed to retain the current rhs
    //! until the next call of this method!
    virtual void GetRHSVal( SingleVector &ptRhs,
                            const PdeIdType identifierPDE
                            = NO_PDE_ID ) = 0;


    // ***********************************************************************
    //   Methods for passing general data
    // ***********************************************************************


    //@{
    //! \name General Data Management

    //! Set type of needed matrix (STIFFNESS, MASS, ...)

    //! This method must be called to inform the algebraic system that a
    //! certain Finite-Element matrix (STIFFNESS, MASS, DAMPING, ...) is
    //! required in the simulation. In the case of an SBM_System we use
    //! the supplied PDE identifiers to determine which sub-matrices of
    //! the specified Finite-Element matrix will be needed.
    //! \note The final problem matrix of type SYSTEM is always generated
    //!       by default.
    //! \param matType type of finite element matrix (STIFFNESS, MASS, ...)
    //! \param identifierPDE1  unique identifier for a PDE registered with the
    //!                        graph manager
    //! \param identifierPDE2  unique identifier for a PDE registered with the
    //!                        graph manager
    virtual void SetFEMatrixType( const FEMatrixType matType,
                                  const PdeIdType identifierPDE1,
                                  const PdeIdType identifierPDE2
                                  = NO_PDE_ID) = 0;

    //! Get types of used matrices (STIFFNESS, MASS, ...)

    //! Get the types of the matrices, which were created by OLAS,
    //! for example STIFFNESS, MASS...
    //! \param matTypes set of enums of global matrices
    virtual void GetFEMatrixTypes( std::set<FEMatrixType> &matTypes ) const;

    //! Set the number of inhomogeneous Dirichlet boundary conditions

    //! This method must be called by each PDE in order to inform the
    //! algebraic system on the number of CFS++ equation numbers that
    //! represent dofs (scalar case) or of blocks of dofs containing dofs
    //! (if dof blocking is used) that are fixed by inhomogeneous Dirichlet
    //! boundary conditions. The method stores this information and uses it
    //! to determine the correct entries for the bcOffsets_ array.
    //! \param identifier unique identifier for the PDE as assigned to it by
    //!                   the ObtainPDEId() method
    //! \param numBCs     number of unknowns fixed by an inhomogeneous
    //!                   Dirichlet boundary condition
    void SetNumDirichletBCs( const PdeIdType identifier, const UInt numBCs );

    //! Get number of iterations specified for an iterative solver
    Integer GetNumIter();
    //@}



    // ***********************************************************************
    //   Various methods
    // ***********************************************************************

    //@{
    //! \name Miscellaneous

    //! Export matrix to a file in MatrixMarket Format

    //! The method will export the matrix to an ascii file according to the
    //! MatrixMarket specifications. For details of the specification see
    //! http://math.nist.gov/MatrixMarket
    //! \param matrixID specifies matrix to be exported
    //! \param filename name of output file
    //! \param comment  string to be inserted into file header (optional)
    virtual void Export( FEMatrixType matrixID, char *filename,
                         char *comment = NULL ) const = 0;

    //! Print the specified matrix into the .las-file

    //! The method will print the specified in the .las-file. The format
    //! depends on the storage type of the matrix, which is normally an
    //! indexed format (for crs and scrs). This mainly used for debugging
    //! purpose.
    //! \param matrixID specifies matrix to be exported
    virtual void Print( FEMatrixType matrixID ) const = 0;

    //! Print information about registered PDEs

    //! The method prints the name of the PDEs, their associated identifiers
    //! and the number of unkonwns per PDE into a specified ostream.
    //! \param log any kind of filestrean (normally .las-file)
    virtual void PrintRegistrationInfo( std::ostream *log ) const;

    /** This method gives back global matrix. For almost all cases you want
     * the SYSTEM matrix */
    virtual StdMatrix* GetSysMat(FEMatrixType type = SYSTEM) const { return NULL;}

    //! returns system matrix (if stored in scrs!) as a vector
    //! containing all upper triangle entries
    virtual void GetFullSystemMatrixAsVec(Complex* &sysVec) const {;};

    //! Removes the IDBC Information from algebraic system

    //! This method is needed in nonlinear computations
    //! while computing solution updates
    virtual void RemoveIDBCInfoFromMatrix() const {;};
    //@}

  protected:

    //! Pointer to the graph manager object
    BaseGraphManager *graphManager_;

    //! Pointer to solver
    BaseSolver *solver_;

    //! Pointer to eigenvalue solver
    BaseEigenSolver *eigenSolver_;

    //@{ \name Attributes storing information on individual PDEs

    //! Number of PDEs

    //! Holds the number of PDEs which is equivalent to the number of
    //! supber blocks in the algebraic system matrix.
    UInt numPDEs_;

    //! Number of registered PDEs so far
    Integer numRegPDEs_;

    //! Number of unknowns per PDE
    UInt *sizePerPDE_;

    //! STL map for associating pde-identifiers with PDE names
    std::map<PdeIdType,std::string> pdeNames_;

    //@}

    //! Size of the matrix

    //! Since we are talking of square matrices the size equals the order
    //! of the matrix, which corresponds to the number of rows resp. columns.
    Integer size_;

    //! Maximal number of non-zero entries in the matrices per row
    Integer nne_;

    //! Degree of freedom (= blocksize in blocksystems)
    Integer blockSize_;

    /** collect system data, but bot the individual preconditioner/ solver data*/
    PtrParamNode systemInfo_;

    // =======================================================================
    // BOUNDARY CONDITIONS / CONSTRAINTS
    // =======================================================================

    //@{ \name Handling of Boundary Conditions and Constraints

    //! Total number of inhomogeneous Dirichlet values over all PDEs
    UInt numDirichletValues_;

    //! Object for handling inhomogeneous Dirichlet boundary conditions
    BaseIDBC_Handler *idbcHandler_;

    //! Index of last equation number belonging to a free dof
    UInt *numLastFreeDof_;

    //! Array containing offsets for managing inhomogeneous Dirichlet BCs

    //! This array contains the offsets for the consecutive numbers
    //! for the inhomogenous Dirichlet boundary conditions.  The array is
    //! enlarged by one entry to be able to compute from it the number of
    //! IDBCs for an individual PDE.
    UInt *bcOffsets_;

    //! Status of system matrix w.r.t. penalty terms

    //! This flag is used to keep track on whether we must insert penalty
    //! terms into the system matrix or not. The latter can be the case if
    //! either the penalty formulation is not employed, or the insertion
    //! was already performed.
    bool assembleDirichletToSysMat_;

    //@}

    //! Pointer to the BaseAssemble object
    //! which assembles the element matrices into the
    //! global matrices and handles the Dirichlet BCs
    BaseEntryManipulator *assemble_;

    //! Set defining which type of matrices (stiffness, mass,...) are used
    std::set<FEMatrixType> matrixTypes_;

    //! Attribute identifying type of algebraic system
    AlgSysType algSysType_;

    /** Here we store the complete ParamNode descripton of our liner system
     * - As given in the XML - hence it might be NULL! */
    PtrParamNode xml_;
    PtrParamNode olasInfo_;
    
    bool usingPenalty_;
  };

} // namespace


#endif // OLAS_BASESYSTEM_HH
