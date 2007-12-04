// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef OLAS_STANDARDSYSTEM_HH
#define OLAS_STANDARDSYSTEM_HH

#include <set>
#include <map>

#include "algsys/basesystem.hh"

namespace OLAS 
{
  // Forward Declarations of classes
  class BaseStdPrecond;
  class StdMatrix;
  class SparseVector;
  class PatternPool;


  //! Linear algebraic system for normal scalar- and blocksystems

  //! This class implements the default user interface for %OLAS in the case
  //! that the system of linear equations that is to be solved is described
  //! by standard matrices (not super block matrices) with scalar or block
  //! entries. It is this interface by which %OLAS and CFS++ will communicate
  //! in this case.
  class StandardSystem : public BaseSystem {


  public:

    //! For convenience and to help g++ parsing we define this shortcut
    typedef std::map<FEMatrixType,Double> factorMap;

    //! Default constructor

    //! This is the default constructor. It will perform some initial memory
    //! allocation (e.g. for the array of matrices) and set some default
    //! values.
    StandardSystem(ParamNode* xml = NULL);

    //! Destructor

    //! The destructor is deep, i.e. it frees dynamically allocated
    //! memory.
    virtual ~StandardSystem();


    // ***********************************************************************
    //   Methods for creating and initalizing the algebraic system
    // ***********************************************************************
   
    //@{
    //! \name Creation and Initialization

    //! Generate objects for the linear system(s)

    /*! The method will generate the matrix and vector objects required for
      the different system matrices and the solution and right hand side
      vector. In the parallel case this function also performs the
      partitioning of the domain. For this task it reads the following
      parameters from myParams_:
      - bool Parallel should the objects created be parallel or sequential.
      the default is false
      - Integer PartitionOverlap if the domain is partitioned (in parallel 
      computations) this indicates how many layers of overlap should be 
      created. The default value is 0.
    */
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
    //! form M(du^2/dt^2)=Ku.
    //! Before this objects can be used the SetupEigenSolver() method should 
    //! be called.
    //! \note If an Eigenfrequency analysis is performed, the methods
    //! SetupPrecond() and SetupSolver() must not be called!
    void CreateEigenSolver();

    //! Trigger setup of preconditioner

    //! Calling this method will trigger the setup phase of the
    //! preconditioner. The setup is performed using the system matrix of the
    //! linear system.
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
    void SetupEigenSolver( UInt numFreq, Double shift, bool quadratic );

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
                               const Double* &err );

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
                               const Double* &err );
    
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
    void CalcEigenMode( UInt numMode );
    
    //@}


    // ***********************************************************************
    //   Methods for assembling the matrices, rhs and bc
    // ***********************************************************************
 
    //@{
    //! \name Assembling Routines

    //! Initialize specified matrix with zeros

    //! Calling this method initialises the Finite-Element matrix designated
    //! by the provided FEMatrixType identifier with zeros. In the case of
    //! a %StandardSystem the PdeIDType identifier is ignored. After
    //! initialisation all matrix entries in the matrix' sparsity pattern
    //! are zero. If no FEMatrixType identifier is specified all matrices
    //! in the internal set (see #matrixTypes_) will be initialised.
    //! \param matrixID type of Finite Element matrix (STIFFNESS, MASS, ...)
    //! \param pdeID    identifier for a PDE as assigned by the ObtainPDEId()
    //!                 method
    void InitMatrix( FEMatrixType matrixID = NOTYPE,
                     const PdeIdType pdeID = NO_PDE_ID );

    //! Set right hand side vector to zero

    //! This method sets the righ hand side vector of a given PDE to zero.
    //! If no PDE identifier is given, the complete rhs vector is zeroed.
    //! \param identifierPDE  identifier for a PDE as assigned by the
    //!                       ObtainPDEId() method
    void InitRHS( const PdeIdType identifierPDE = NO_PDE_ID );

    //! Set the global rhs to given values

    //! Set the right hand side vector to the values specified in the array
    //! newRHS
    //! \param newRHS Pointer to new right-hand-side values
    //! \note The values of newRHS are copied, so the pointer to newRHS
    //! can be changed afterwards!
    void InitRHS(const Double *  newRHS);
    
    /** @see BaseSystem::InitRHS(const Vector<Complex>*) */
    void InitRHS( const Vector<Complex>* newRHS );
    
    //! Set the solution vector of the specified PDE to zero
    
    //! This method initializes the solution associated with the given
    //! PDE identifier. If no identifier is specified, the global solution
    //! vector is zeroed.
    //! \param identifierPDE unique identifier for a PDE registered with the
    //!                       graph manager
    void InitSol( const PdeIdType identifierPDE = NO_PDE_ID );
    
    //! Set the global solution vector to given initial values

    //! Set the solution vector to the values specified in the array
    //! newSol
    //! \param newSol Pointer to new solution values
    //! \note The values of newSol are copied, so the pointer to newSol
    //! can be changed afterwards!
    void InitSol(const Double *  newSol, const UInt size);
    
    
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
    //!                insert the element matrix \f$E\f$, In doing so the
    //!                row and column indices derived from the equation
    //!                numbers are interchanged. Note that this is only
    //!                supported for off-diagonal blocks, i.e. for cases
    //!                with different PDE identifiers.
    //! \param setTransposeInt if this flag is true, then the method will
    //!                insert the transpose of the element matrix \f$E^T\f$. 
    //! \note We currently use setCounterPart = true as default, since we
    //!       have no clear concept how to administer this information within
    //!       CFS++ up to now. In the case of a diagonal matrix block, we
    //!       re-set this to false before passing the flag on to assemble.
    void SetElementMatrix( FEMatrixType matrix_id, Double * elemmat, 
                           PdeIdType identifierPDE1,
                           Integer *eqnNrs1,
                           Integer numEqn1,
                           PdeIdType identifierPDE2,
                           Integer *eqnNrs2,
                           Integer numEqn2,
                           FEMatrix_Flags pFlags );
    
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
    void SetElementRHS( Double *elemRHS, const PdeIdType idPDE,
                        Integer *eqnNrs, UInt numEqn );

    //@{
    //! Adds a value to a given global rhs entry
    
    //! This method adds a given value to the global right hand side vector.
    //! Therefore a pde identifier, the related equation number and a 
    //! degree of freedom has to be specified.
    //! \param val value to be added
    //! \param identifierPDE identifier of the PDE related to sub-graph
    //! \param eqnNr equation number of the node to be set
    void SetNodeRHS(Double val, PdeIdType identifierPDE,
                    Integer eqnNr);
    
    void SetNodeRHS(Complex val, PdeIdType identifierPDE,
                    Integer eqnNr);
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
    void UpdateRHS(FEMatrixType matrix_id, Double * fup);
    
    //! Add a value to a diagonal matrix entry

    //! This method allows to directly add a value to the value of an entry
    //! on the main diagonal of a certain of the different problem matrices.
    //! \param matrixID specifies which matrix is to be manipulated
    //! \param pdeID    identifier for PDE related to sub-graph
    //! \param eqnNum   determines the diagonal entry to be manipulated, i.e.
    //!                 the parameter gives the equation number of the unknown
    //!                 coupled to the row/column of the entry; e.g. for
    //!                 \f$\mbox{eqnNum}=k\f$ it is \f$a_{kk}\f$ that will be
    //!                 altered
    //! \param val      zero-based array containing the numerical value that
    //!                 is to be added to the matrix entry; in the case of
    //!                 real entries, the first array value is used, in the
    //!                 case of complex entries the first entry is used as
    //!                 real and the second as imaginary part
    void AddToDiagMatrixEntry( FEMatrixType matrixID, const PdeIdType pdeID,
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

    //! Constructs the effective system matrix

    //! Constructs the effective system matrix by multiplying the individual
    //! matrices with factors and adding them to one system matrix
    //! \param matFactors a map which contains for each matrixtype (DAMPING,
    //!                   STIFFNESS,...) the according factor
    //! \note It is allowed to supply no matrix factors only in the case
    //!       that the only existing FE matrix is the system matrix. In this
    //!       case no arithmetic operations are performed. If there is more
    //!       than one matrix we terminate.
    void ConstructEffectiveMatrix( const factorMap &matFactors =
                                   factorMap() );

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

    //! This method triggers the assembling of the dirichlet values into
    //! the algebraic system. The values itself must have been specified
    //! before via the SetDirichlet() method.
    //! \note After assembling the dirichlet values, the preconditioner and
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
                            = NO_PDE_ID );
 
    Integer GetSolutionVal( Complex* &ptSol, 
                            const PdeIdType identifierPDE 
                            = NO_PDE_ID );
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
                       = NO_PDE_ID);
  
    Integer GetRHSVal( Complex* &ptRhs, 
                       const PdeIdType identifierPDE 
                       = NO_PDE_ID);
    //@}


    // ***********************************************************************
    //   Methods for passing general data
    // ***********************************************************************

    //@{
    //! \name General Data Management
    
  
    //! Set type of needed matrix (STIFFNESS, MASS, ...)

    //! Set type of needed matrix (STIFFNESS, MASS, ...)
    //! \param matType type of finite element matrix (STIFFNESS, MASS, ...)
    //! \param identifierPDE1 unique identifier for a PDE registered with the
    //!                       graph manager
    //! \param identifierPDE2 unique identifier for a PDE registered with the
    //!                       graph manager
    void SetFEMatrixType( const FEMatrixType matType,
                          const PdeIdType identifierPDE1,
                          const PdeIdType identifierPDE2 = NO_PDE_ID );

    //! Set block size of a matrix entry

    //! The number of degrees of freedom associated with a single equation
    //! number depends on the type of equation numbering in CFS. In the
    //! simple case of a scalar numbering each equation number represents
    //! precisely one unknown. If dof blocking is used the matrices in OLAS
    //! will be composed of tiny submatrices and the block size will be
    //! larger than one.
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
    //! http://math.nist.gov/MatrixMarket
    //! \param type     specifies matrix to be exported
    //! \param filename name of output file
    //! \param comment  string to be inserted into file header (optional)
    void Export( FEMatrixType type, Char *filename,
                 Char *comment = NULL ) const;

    //! Print the specified matrix into the .las-file

    //! The method will print the specified in the .las-file. The format
    //! depends on the storage type of the matrix, which is normall an indexed
    //! format (for crs and scrs). This mainly used for debugging purpose.
    //! \param matrix_id specifies matrix to be exported
    void Print(FEMatrixType matrix_id) const;

    //! removes IDBS Information from algebraic system
    void RemoveIDBCInfoFromMatrix() const; 

    StdMatrix* GetSysMat(FEMatrixType type = SYSTEM)const {return sysmat_[type];}
    SparseVector *GetSolVec(){return sol_;}
    SparseVector *GetRhsVec(){return rhs_;}

    //@}
  
    
  private:

    //! Auxilliary method for generating pass-back vectors for CFS++
    void GenerateCFSTransferBuffers();

    //! Pointer to the systemmatrix (sysmat[1])
    //! and the partial matrices (mass, stiffnes, ...)
    StdMatrix** sysmat_;

    //! Pointer to solution vector
    SparseVector *sol_;

    //! Pointer to right-hand side vector
    SparseVector *rhs_;

    //! Auxilliary vector for passing solution back to CFS++

    //! This auxilliary vector is used for passing the solution back to
    //! CFS++.
    //! \note This is a temporary solution for the non-penalty approach to
    //!       IDBC handling
    //! \todo Adapt handling of solution vector in StoreSolutions classes
    //!       of CFS++ to administer IDBCs itself
    SparseVector *solComm_;

    //! Buffer for storing the eigenvalues of the system
    SparseVector *eigenValues_;

    //! Buffer for storing the error bounds of the eigenvalues
    SparseVector *eigenValError_;

    //! Pointer to the preconditioner object
    BaseStdPrecond *precond_; 
    
    //! Size of system including fixed dofs
    UInt totalSize_;

    //! Pointer to a pool for sharing sparsity patterns between matrices
    PatternPool *patternPool_;
  };

}

#endif // OLAS_STANDARDSYSTEM_HH
