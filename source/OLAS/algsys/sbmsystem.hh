// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef OLAS_SBM_SYSTEM_HH
#define OLAS_SBM_SYSTEM_HH

#include <set>
#include <map>

#include "algsys/basesystem.hh"
#include "utils/utils.hh"


namespace OLAS {


  // Forward Declarations of classes
  class BaseSBMPrecond;
  class StdMatrix;
  class SparseVector;
  class SBM_Matrix;
  class SBM_Vector;
  class BaseEntryManipulator;
  class BaseIDBC_Handler;
  class SBM_Matrix;


  //! Linear algebraic system for normal scalar- and blocksystems

  /*! This class implements the default user interface for %OLAS in the case
   *hat the system of linear equations that is to be solved is described
   * by standard matrices (not super block matrices) with scalar or block
   * entries. It is this interface by which %OLAS and CFS++ will communicate
   * in this case.
   */
  class SBM_System : public BaseSystem {


  public:

    //! Triple identifying a sub-matrix of a Finite-Element matrix
    typedef struct {
      FEMatrixType feType;
      PdeIdType rowInd;
      PdeIdType colInd;
    } FeSubMatrixID;

    //! Touple storing row and column index of an Finite-Element sub-matrix
    typedef struct {
      PdeIdType rowInd;
      PdeIdType colInd;
    } SubMatrixID;

    //! For convenience and to help g++ parsing we define this shortcut
    typedef std::map<FEMatrixType,Double> factorMap;

    //! Default constructor

    //! This is the default constructor. It will perform some initial memory
    //! allocation (e.g. for the array of matrices) and set some default
    //! values.
    SBM_System(ParamNode* xml = NULL);

    //! Destructor

    //! The destructor is deep, i.e. it frees dynamically allocated
    //! memory.
    virtual ~SBM_System();


    // ***********************************************************************
    //   Methods for creating and initalizing the algebraic system
    // ***********************************************************************

    //! Generate objects for the linear system

    //! The method will generate the matrix and vector objects required for
    //! the different system matrices and the solution and right hand side
    //! vector.
    void CreateLinSys();

    //! Generate solver object

    //! Calling this method triggers the generation of a solver object.
    //! Before this objects can be used the SetupSolver() method should be
    //! called and also the generation of a preconditioner via the
    //! CreatePrecond() and SetupPrecond() methods should be finished, if
    //! a preconditioner is to be used.
    //! \note We currently always generate a GMRES solver
    void CreateSolver();

    //! Generate preconditioner
    void CreatePrecond() {
      (*warning) << "SBM_System::CreatePrecond not yet implemented!";
      Warning( __FILE__, __LINE__ );
    }

    //! Generate EigenSolver object

    //! Calling this method triggers the generation of an EigenSolver object 
    //! which ist used to calculate a generalized eigenvalue problem of the
    //! form M(du^2/dt^2)=Ku.
    //! Before this objects can be used the SetupEigenSolver() method should 
    //! be called.
    //! \note If an Eigenfrequency analysis is performed, the methods
    //! SetupPrecond() and SetupSolver() must not be called!
    void CreateEigenSolver() {
      (*warning) << "SBM_System::CreateEigenSolver not yet implemented!";
      Warning( __FILE__, __LINE__ );
    } 

    //! Trigger setup of preconditioner

    //! Calling this method will trigger the setup phase of the
    //! preconditioner. The setup is performed using the system
    //! matrix of the linear system.
    //! \note This method must not be called if an eigenfrequency analysis
    //! is performed, since this method creates only a preconditioner 
    //! to solver which solves a system Ax=b.
    void SetupPrecond();

    //! Trigger setup of solution method

    //! Calling this method will trigger the setup phase of the solver. This
    //! is especially important for direct solvers, where typically the
    //! factorisation of the problem matrix will be performed at this stage.
    //! The setup is performed using the system matrix of the linear system.
    void SetupSolver();

        //! Trigger setup of eigenvalue solver

    //! Calling this method will trigger the setup phase of the eigensolver. 
    //! The setup is performed using the system matrix of the linear system.
    //! \param numFreq Number of eigenfrequencies to be calculated
    //! \param shift Frequency shift applied to the system
    //! \param quadratic Flag indicating if a quadratic eigenvalue problem
    //! (true) or a generalized problem (false) is to be solved
    void SetupEigenSolver( UInt numFreq, Double shift, bool quadratic ) {
      (*warning) << "SBM_System::SetupEigenSolver not yet implemented!";
      Warning( __FILE__, __LINE__ );
    } 

    //! Solve the linear system

    //! Calling this method triggers the solution of the linear system
    //! defined by the Finite-Element system matrix, i.e. sysMat_[SYSTEM],
    //! and the given right-hand side vector, i.e. #rhs_. The computed
    //! (approximat) solution is stored in #sol_.
    //! The method used for solving the system depends on the solver and 
    //! preconditioner constructed and the parameters in myParams_.
    //! For iterative solvers an initial guess is created by inserting the
    //! Dirichlet values in the correct positions of the solution vector in
    //! case of the penalty formulation.
    //! \note This method must not be called if an eigenfrequency analysis
    //! is performed, since this method is only used to solve a system of the
    //! form Ax=b.
    void Solve(const std::string& comment = "");

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
    //! \note This method may only be called if SetupEigenfrequencySolver()
    //!       was called previously.
    UInt CalcEigenFrequencies( const Double* &frequencies,
                               const Double* &norm  ) {
      (*warning) << "SBM_System::CalcEigenFrequencies not yet implemented!";
      Warning( __FILE__, __LINE__ );
      return 0;
    } 

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
    UInt CalcEigenFrequencies( const Complex* &frequencies,
                               const Double* &err ) {
    (*warning) << "SBM_System::CalcEigenFrequencies not yet implemented!";
      Warning( __FILE__, __LINE__ );
      return 0;
    } 

    //! Calculate eigenmodes of a generalized/quadratic eigenvalue problem

    //! Calling this method triggers the calculation of the n-th eigenmode  
    //! corresponding to the n-th eigenvalue.
    //! This method may only be called after CalcEigenfrequencies(), as
    //! the eigenmodes are only a postprocessing result.
    //! The calculated mode can be accessed by GetSolutionVal().
    //! \param numMode Mode number corresponding to the numMode-th 
    //!                eigenfrequency
    //! \note This method may only ba called if SetupEigenfrequencySolver()
    //!       and CalcEigenfrequencies was called previously.
    void CalcEigenMode( UInt numMode ) {
      (*warning) << "SBM_System::CalcEigenMode not yet implemented!";
      Warning( __FILE__, __LINE__ );
    } 


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
    //! \param matrixID  type of Finite Element matrix (STIFFNESS, MASS, ...)
    //! \param pdeID     identifier obtained from the ObtainPDEId() method
    void InitMatrix( FEMatrixType matrixID = NOTYPE,
                     const PdeIdType pdeID = NO_PDE_ID );

    //! Set right hand side vector to zero

    //! This method sets the right hand side vector associated with the
    //! specified PDE to zero. If no PDE identifier is given, the complete
    //! right-hand side vector is zeroed.
    //! \note In the case of a StandardSystem the PDE identifier is ignored.
    //!       We currently do not support setting only the part of a SparseVector
    //!       related to a single PDE to zero in this case.
    //! \param identifierPDE unique identifier obtained from the ObtainPDEId()
    //!                      method
    void InitRHS( const PdeIdType identifierPDE = NO_PDE_ID );

    //! Set the global rhs to given values

    //! Set the right hand side vector to the values specified in the array
    //! newRHS
    //! \param newRHS Pointer to new right-hand-side values
    //! \note The method is not yet implemented and with the current interface
    //!       probably also useless. We should probably pass a PDE identifier
    //!       here.
    void InitRHS( const Double *newRHS ) {
      (*error) << "Method not yet implemented!";
      Error( __FILE__, __LINE__ );
    }

    //! Set the solution vector of the specified PDE to zero

    //! This method sets the part of the solution vector associated with the
    //! specified PDE to zero. If no PDE identifier is given, the complete
    //! solution vector is zeroed.
    //! \note In the case of a StandardSystem the PDE identifier is ignored.
    //!       We currently do not support setting only the part of a SparseVector
    //!       related to a single PDE to zero.
    //! \param identifierPDE unique identifier obtained from the ObtainPDEId()
    //!                      method
    void InitSol( const PdeIdType identifierPDE = NO_PDE_ID );

    //! Set the global solution vector to given initial values

    //! Set the solution vector to the values specified in the array
    //! newSol
    //! \param newSol Pointer to new solution values
    //! \note The values of newSol are copied, so the pointer to newSol
    //! can be changed afterwards!
    //! not in use
    void InitSol(const Double *  newSol, const UInt size){};
    
    
    //! Assemble an element matrix into the global one

    //! This methods assembles the given element matrix into a specified
    //! global one (MASS, STIFFNESS, etc.), i.e. adds the entries to the. 
    //! global matrix. For element matrices which are associated
    //! only with one PDE, only one identifier, the eqn numbers of the matrix
    //! and the number of eqnNrs have to be specified. 
    //! For matrices which are associated with two different pde identifiers,
    //! the matrix will be assembled into the upper off-diagonal matrix-block
    //! and the lower transposed one.
    //! \param matrixID type of finite element destination matrix 
    //!                 (STIFFNESS, MASS, ...)
    //! \param elemMat  entries of the element matrix in sequential array
    //! \param pdeID1   identifier for first PDE related to sub-graph
    //! \param eqnNrs1  equation numbers (1-based) of the element matrix
    //!                 w.r.t. sub-graph associated with identifierPDE1
    //! \param numEqn1  number of equations related to sub-graph of 
    //!                 identifierPDE1
    //! \param pdeID2   identifier for first PDE related to sub-graph
    //! \param eqnNrs2  equation numbers (1-based) of the element matrix
    //!                 w.r.t. sub-graph associated with identifierPDE2
    //! \param numEqn2  number of equations related to sub-graph of 
    //!                 identifierPDE2
    //! \param setCounterPart if this flag is true, then the method will
    //!                 not only insert the element matrix \f$E\f$, but also
    //!                 its transpose \f$E^T\f$. In doing so also the
    //!                 row and column indices derived from the equation
    //!                 numbers are interchanged. Note that this is only
    //!                 supported for off-diagonal blocks, i.e. for cases
    //!                 with different PDE identifiers.
    void SetElementMatrix( FEMatrixType matrixID, Double *elemMat, 
                           PdeIdType pdeID1, Integer *eqnNrs1,
                           Integer numEqn1, PdeIdType pdeID2,
                           Integer *eqnNrs2, Integer numEqn2,
                           bool setCounterPart );

    //! Assemble the local rhs vector to the global one

    //! This method adds the entries of the element right hand side to the
    //! global ones, according the pdeIdentifier and the associated equation
    //! numbers.
    //! \param elemRHS entries of the element right hand side in a sequential
    //!                array
    //! \param pdeID   identifier of the PDE related to sub-graph
    //! \param eqnNrs  equation numbers (1-based) of the element rhs
    //!                w.r.t. sub-graph associated with idPDE
    //! \param numEqn  length of eqnNrs array
    void SetElementRHS( Double *elemRHS, const PdeIdType pdeID,
                        Integer *eqnNrs, UInt numEqn );


    //@{
    //! Adds a value to a given global rhs entry

    //! This method adds a given value to the global right hand side vector.
    //! Therefore a pde identifier, the related equation number and a 
    //! degree of freedom has to be specified.
    //! \param val value to be added
    //! \param identifierPDE identifier of the PDE related to sub-graph
    //! \param eqnNr equation number of the node to be set
    void SetNodeRHS( Double val, PdeIdType identifierPDE,
                     Integer eqnNr ) {
      (*error) << "Method not yet implemented!";
      Error( __FILE__, __LINE__ );
    }

    void SetNodeRHS( Complex val, PdeIdType identifierPDE,
                     Integer eqnNr) {
      (*error) << "Method not yet implemented!";
      Error( __FILE__, __LINE__ );
    }
    //@}

    //! Performs a matrix-vector multiplication and adds the vector to the rhs

    //! This method multiplies a specified global matrix (STIFFNES, MASS, 
    //! etc.) with a given vector and adds the result vector to the global 
    //! rhs.
    //! \f[ rhs = rhs + \mathbf A_{matrix_id} \cdot fup \f]
    //! This method is mainly used in a time step scheme to adapt the rhs.
    //! \param matrix_id type of finite element matrix (STIFFNESS, MASS, ...)
    //!                  which gets multiplied
    //! \param fup array with vector entries, which get multiplied
    void UpdateRHS( FEMatrixType matrix_id, Double *fup );

    //! Add a value to a diagonal matrix entry

    //! This method allows to directly add a value to the value of an entry
    //! on the main diagonal of a certain of the different problem matrices.
    //! \param matrixId  specifies which matrix is to be manipulated
    //! \param pdeID     identifier for the corresponding PDE sub-matrix
    //! \param eqnNum    determines the diagonal entry to be manipulated, i.e.
    //!                  the parameter gives the equation number of the
    //!                  unknown coupled to the row/column of the entry;
    //!                  e.g. for \f$\mbox{eqnNum}=k\f$ it is \f$a_{kk}\f$
    //!                  that will be altered
    //! \param val       zero-based array containing the numerical value that
    //!                  is to be added to the matrix entry; in the case of
    //!                  real entries, the first array value is used, in the
    //!                  case of complex entries the first entry is used as
    //!                  real and the second as imaginary part
    void AddToDiagMatrixEntry( FEMatrixType matrixId, const PdeIdType pdeID,
                               Integer eqnNum, Double *val );

    //@{
    //! Get value of specific matrix entry

    //! This method allows to directly get the value of a specific matrix
    //! entry of a given problem matrix. If the given index pair is not 
    //! contained in the matrix, an error is thrown.
    //! \param matrixID  specifies which matrix is to be queried
    //! \param rowPdeID identifier for PDE related to row sub-graph
    //! \param rowEqnNum equation number of the row index to get
    //! \param colPdeID  identifier for PDE related to column sub-graph
    //! \param colEqnNum equation number of the column index to get
    //! \param val       return value
    //! \note For sparse matrices this method may be very costly and slow,
    //!       as the given index pair has to be searched for and may not
    //!       be accessible directly!
    void GetMatrixEntry( FEMatrixType matrixID,
                         const PdeIdType rowPdeID,
                         Integer eqnNum1,
                         const PdeIdType colPdeID,
                         Integer eqnNum2,
                         Double & val );

    void GetMatrixEntry( FEMatrixType matrixID,
                         const PdeIdType rowPdeID,
                         Integer rowEqnNum, 
                         const PdeIdType colPdeID,
                         Integer colEqnNum2,
                         Complex & val );
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
    void SetMatrixEntry( FEMatrixType matrixID,
                         const PdeIdType rowPdeID,
                         Integer rowEqnNum,
                         const PdeIdType colPdeID,
                         Integer colEqnNum,
                         Double val, bool setCounterPart );
    
    void SetMatrixEntry( FEMatrixType matrixID,
                         const PdeIdType rowPdeID,
                         Integer rowEqnNum, 
                         const PdeIdType colPdeID,
                         Integer colEqnNum, 
                         Complex val, bool setCounterPart );
    //@}
    
    //@{
    //! Initialize all entries of a specific row with given value

    //! This method allows to set directly the value of all (potentially 
    //! non-zero) entries of a matrix row of a specified problem matrix.
    //! \param matrixID  specifies which matrix is to be manipulated
    //! \param pdeID     identifier for PDE related to sub-graph
    //! \param eqnNum    equation number of the row index to be set
    //! \param dof       degree of freedom of the row to be set; in the case 
    //!                  of scalar entries this should be 1, only in the case
    //!                  of block entries may this be larger than 1
    //! \param val       value the row gets initialized with
    void SetMatrixRowVals( FEMatrixType matrixID,
                           const PdeIdType pdeID,
                           Integer eqnNum, UInt dof,
                           Double val );
    
    void SetMatrixRowVals( FEMatrixType matrixID,
                           const PdeIdType pdeID,
                           Integer eqnNum, UInt dof,
                           Complex val );
    //@}


    //@{
    //! Initialize all entries of a specific column with given value

    //! This method allows to set directly the value of all (potentially 
    //! non-zero) entries of a matrix column of a specified problem matrix.
    //! \param matrixID  specifies which matrix is to be manipulated
    //! \param pdeID     identifier for PDE related to sub-graph
    //! \param eqnNum    equation number of the column index to be set
    //! \param dof       degree of freedom of the column to be set; in the case 
    //!                  of scalar entries this should be 1, only in the case
    //!                  of block entries may this be larger than 1
    //! \param val       value the column gets initialized with
    void SetMatrixColVals( FEMatrixType matrixID,
                           const PdeIdType pdeID,
                           Integer eqnNum, UInt dof,
                           Double val );

    void SetMatrixColVals( FEMatrixType matrixID,
                           const PdeIdType pdeID,
                           Integer eqnNum, UInt dof,
                           Complex val );
    //@}

    //! Constructs the effective system matrix

    //! Constructs the effective system matrix by multiplying the individual
    //! matrices with factors and adding them to one system matrix
    //! \param matFactors a map which contains for each matrixtype (DAMPING,
    //!                   STIFFNESS,...) the according factor
    void ConstructEffectiveMatrix( const factorMap &matFactors = factorMap());

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
    void SetDirichlet( const PdeIdType pdeID, Integer eqnNr,
                       const Double &val );

    void SetDirichlet( const PdeIdType pdeID, Integer eqnNr,
                       const Complex &val );
    //@}

    //! Assemble the previously defined Dirichlet values into the algsys

    //! This method triggers the assembling of the Dirichlet values into
    //! the algebraic system. The values itself must have been specified
    //! before via the SetDirichlet() method.
    //! \note After assembling the Dirichlet values, the preconditioner and
    //! the solver have to be set up again.
    void BuildInDirichlet();

    //@{
    //! Return the pointer to the current solution of a PDE

    //! This method passes the current solution of one given PDE
    //! as a pointer to a buffer containing the solution entries. 
    //! If no identifier is given, the pointer to the beginning of the
    //! complete solution array is passed.
    //!
    //! \param ptSol pointer (0-based) to the solution buffer
    //! \param identifierPDE identifier for PDE related to sub-graph
    //! \return number of array entries
    //! 
    //! \note The return buffer is guaranteed to retain the current solution
    //! until the next call of this method (after solving the next step)!
    Integer GetSolutionVal( Double* &ptSol, 
                            const PdeIdType identifierPDE 
                            = NO_PDE_ID ) {
      (*error) << "Method not yet implemented!";
      Error( __FILE__, __LINE__ );
      return 0;
    }

    Integer GetSolutionVal( Complex* &ptSol, 
                            const PdeIdType identifierPDE 
                            = NO_PDE_ID ) {
      (*error) << "Method not yet implemented!";
      Error( __FILE__, __LINE__ );
      return 0;
    }
    //@}

    //@{
    //! Return the pointer to the current rhs value of a PDE
    
    //! This method passes the current rhs as a pointer to buffer with
    //! right hand side entries.
    //!
    //! \param ptRhs pointer (0-based) to the buffer with rhs values
    //! \param identifierPDE identifier for PDE related to sub-graph
    //! \return number of array entries
    //!
    //! \note The return buffer is guaranteed to retain the current rhs
    //! until the next call of this method!
    Integer GetRHSVal( Double* &ptRhs, 
                       const PdeIdType identifierPDE 
                       = NO_PDE_ID) {
      (*error) << "Method not yet implemented!";
      Error( __FILE__, __LINE__ );
      return 0;
    }
  
    Integer GetRHSVal( Complex* &ptRhs, 
                       const PdeIdType identifierPDE 
                       = NO_PDE_ID) {
      (*error) << "Method not yet implemented!";
      Error( __FILE__, __LINE__ );
      return 0;
    }
    //@}


    // ***********************************************************************
    //   Methods for passing general data
    // ***********************************************************************

    //@{
    //! \name General Data Management

    //! Set type of needed matrix (STIFFNESS, MASS, ...)

    //! This method must be called to inform the algebraic system that a
    //! certain Finite-Element matrix (STIFFNESS, MASS, DAMPING, ...) is
    //! required in the simulation. In the case of an %SBM_System we use
    //! the supplied PDE identifiers to determine which sub-matrices of
    //! the specified Finite-Element matrix will be needed.
    //! \note The final problem matrix of type SYSTEM is always generated
    //!       by default.
    //! \param matType type of finite element matrix (STIFFNESS, MASS, ...)
    //! \param pdeID1  unique identifier for a PDE as obtained by the
    //!                ObtainPDEId() method
    //! \param pdeID2  unique identifier for a PDE as obtained by the
    //!                ObtainPDEId() method
    void SetFEMatrixType( const FEMatrixType matType, const PdeIdType pdeID1,
                          const PdeIdType pdeID2 = NO_PDE_ID );

    //! Set block size of a matrix entry

    //! The number of degrees of freedom associated with a single equation
    //! number depends on the type of equation numbering in CFS++. In the
    //! simple case of a scalar numbering each equation number represents
    //! precisely one unknown. If dof blocking is used the matrices in %OLAS
    //! will be composed of tiny submatrices and the block size will be
    //! larger than one.
    //! \note Since the SBM_Matrix class only supports scalar entries,
    //!       this method generates an error for bs > 1.
    //! \param identifier unique identifier for a PDE registered with the
    //!                   graph manager
    //! \param bs         number of unknonws per equation number
    //!                   (1 = scalar matrix, > 1 = block matrix)
    void SetBlockSize( const PdeIdType identifier, const UInt bs );

    //@}

  

    // ***********************************************************************
    //   Various methods 
    // ***********************************************************************
    
    //@{
    //! \name Miscellaneous

    //! Export matrix to a file in MatrixMarket Format

    //! The method will export the matrix to an ascii file according to the
    //! MatrixMarket specifications. For details of the specification see
    //! http://math.nist.gov/MatrixMarket. The method makes use of the
    //! SBM_Matrix::Export() method.
    //! \param matrixID specifies matrix to be exported
    //! \param filename name of output file
    //! \param comment  string to be inserted into file header (optional)
    void Export( FEMatrixType matrixID, Char *filename,
                 Char *comment = NULL ) const;

    //! Print the specified matrix into the .las-file

    //! The method will print the specified in the .las-file. The format
    //! depends on the storage type of the matrix, which is normally an
    //! indexed format (for crs and scrs). This mainly used for debugging
    //! purpose.
    //! \param matrixID specifies matrix to be exported
    void Print(FEMatrixType matrixID ) const {
      (*error) << "Method not yet implemented!";
      Error( __FILE__, __LINE__ );
    }
    //@}


  private:

    //! Auxilliary method for generating pass-back vectors for CFS++
    void GenerateCFSTransferBuffers() {
      (*error) << "Method not yet implemented!";
      Error( __FILE__, __LINE__ );
    }

    //! Auxilliary method for logging information on matrix patterns

    //! This method can be used to print a summary of the Finite-Element
    //! matrices of the algebraic system and their sub-matrix pattern to
    //! a stream. It is mostly intended as an auxilliary method for
    //! CreateLinSys().
    void PrintFeMatrixInfo( std::ostream *os );
    SBM_Matrix* GenerateSBM_Matrix( FEMatrixType matType,
                                    MatrixEntryType entryType );

    //! Pointers to Finite-Element matrices

    //! This attribute points to an array containing pointers to the different
    //! Finite-Element matrices (SYSTEM, MASS, STIFFNESS, ... ) managed by the
    //! algebraic system. In the case of the %SBM_System class these are of
    //! course SBM_Matrix objects.
    SBM_Matrix **sysMat_;

    //! Vector containing right-hand side of the linear system
    SBM_Vector *rhs_;

    //! Vector containing (approximate) solution of the linear system
    SBM_Vector *sol_;

    //! Flag signalling symmetry of SBM matrices
    bool sbmSymm_;

    //! Pointer to the preconditioner object
    BaseSBMPrecond *precond_; 

    //! Binary predicate used as comparision operator for SubMatrixID sets

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

    //! Vector with IDs of sub-matrices of Finite-Element matrices

    //! This vector contains for each possible Finite-Element matrix managed
    //! by the algebraic system object a set containing the index tuples for
    //! the non-zero sub-matrices of the Finite-Element matrix. These sets
    //! are filled during the setup phase by calls to SetFEMatrixType() and
    //! used in the CreateLinSys() method for generation of the SBM_Matrix
    //! objects.
    std::set<SubMatrixID,SortSubMatrixID> *feSubMatrices_;

  };

}

#endif
