#ifndef PARDISO_SOLVER_PRIMITIVE_HH
#define PARDISO_SOLVER_PRIMITIVE_HH

#include <string>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Utils/StdVector.hh"
#include "OLAS/solver/BaseSolver.hh"
#include "Utils/Timer.hh"

namespace CoupledField {

  //! This class implements the interface to PARDISO's LU decomposition
  //! and solving of a sparse system of linear equations.
  //! The difference to the class PardisoSolver is that this class
  //! does not get any xml-parameters, since it's ONLY used as a
  //! direct solver for the algebraic-multigrid-preconditioner (AMG)
  //! and the parameters should not get altered via an xml-input
  //! There is no need to initialise it.
  //! After Solve() is called, the BaseMatrix (expected to be either a CRS
  //! or SCRS-matrix) will be reordered using the Nested Dissection or
  //! the Minimum Degree Algorithm, then it will be LU-factorized and
  //! finally solved by backward-forward substitution.

  template<typename T>
  class PardisoSolverPrimitive : public BaseDirectSolver {

  public:

    //! Constructor

    //! The constructor does not do anything but set the pointers to the
    //! internal communciation objects.
    //! Note that the pointer to the report object is optional.
    PardisoSolverPrimitive (PtrParamNode solverNode, PtrParamNode olasInfo);

    //! Default Destructor

    //! The destructor is not used at all.
    ~PardisoSolverPrimitive();

    //! Reordering and factorization of the linear system

    //! After Setup is called the BaseMatrix (expected to be either a CRS
    //! or SCRS-matrix) will be reordered using the Nested Dissection or
    //! the Minimum Degree Algorithm and then it will be LU-factorised.
    void Setup( BaseMatrix &sysmat);

    //! Set the flag that a new matrix pattern shall be used in the setup phase
    void SetNewMatrixPattern();

    //! Direct solution of the linear system

    //! After Solve is called the matrix (which has already to be factorised
    //! by a call of Setup) is finally solved by backward-forward substitution.
    void Solve( const BaseMatrix &sysmat,
                const BaseVector &rhs, BaseVector &sol );

    //! Query type of this solver.

    //! This method can be used to query the type of this solver. The answer
    //! is encoded as a value of the enumeration data type SolverType.
    //! \return PARDISO
    SolverType GetSolverType() {
      return PARDISO_SOLVER;
    }

  private:

    //! Default constructor

    //! The default constructor is private, since we want the solver object to
    //! be initialised with a pointer to a parameter object right at
    //! instantiation.
    PardisoSolverPrimitive();

    //! enum of PARIDO's error codes
    typedef enum {
      PARDISO_NO_ERROR = 0,
      INPUT_INCONSISTENT = -1,
      NOT_ENOUGH_MEMORY = -2,
      REORDERING_PROBLEM = -3,
      ZERO_PIVOT = -4,
      INTERNAL_ERROR = -5,
      PREORDERING_FAILED = -6,
      DIAGONAL_MATRIX = -7,
      INT_OVERFLOW = -8,
      NOT_ENOUGH_OOC_MEM = -9,
      NO_LIC_FILE = -10,
      LIC_EXPIRED = -11,
      WRONG_USER_OR_HOSTNAME = -12,
      MAX_KRYLOV_ITERATIONS = -100,
      INSUFF_KRYLOV_CONVERGENCE = -101,
      KRYLOV_ITERATION_ERROR = -102,
      KRYLOV_BREAKDOWN = -103
    } PardisoError;

    //! Returns a string describing the given error code
    std::string GetErrorString(int err_code);

    //@{
    //! Vectors containing the information about the CRS- or SCRS-matrix

    //! This pointer is used to hold the address of a part of the internal
    //! (S)CRS matrix structures. The related memory segment must not
    //! be altered of deleted. Therefore the pointer is const!
    Integer *rowPtr_;
    Integer *colPtr_;
    const T *datPtr_;
    //@}

    //! Stored information about the storage type and entry type of the matrix
    BaseMatrix::StorageType stype;
    BaseMatrix::EntryType etype;

    //! A working array for the Pardiso-Routines

    //! This is the internal memory address array used by Pardiso. It must
    //! contain 64 entries, which are 4-byte void pointers on a 32-bit
    //! architecture and 8-byte void pointers on a 64-bit architecture.
    StdVector<void*> pt_;

    //! Array for passing information to and from Pardiso

    //! This integer array is used to communicate parameters to and from
    //! Pardiso. Instead of Fortran-/one-based indexing we use C-/zero-
    //! based element indexing.
    StdVector<int> iparm_;

    //! The type of the matrix in a special encoding used by Pardiso
    int mType_;

    //! The type of the solver used by Pardiso. Possible values are zero for 
    //! sparse direct solver (0) or one for multi-recursive iterative solver (1).
    int mSolver_;

    //! This double array is used to communicate parameters to and from
    //! Pardiso. It has been introduced in Pardiso 4.0. Instead of 
    //! Fortran-/one-based indexing we use C-/zero-based element indexing.
    StdVector<double> dparm_;
    
    //! Maximal number of factors with identical nonzero-structure Pardiso
    //! should keep in memory at the same time.
    int maxfct_;

    //! The number of the matrix (out of the maxfct ones) that should be used
    //! for the solution process
    int mnum_;

    //! Dimension of the linear system
    int probDim_;

    //! The number of right hand sides Pardiso should solve the system for
    //! at one pass
    int nrhs_;

    //! Specifies verbosity of Pardiso itself

    //! This attribute stores a number specifying the verbosity of the Pardiso
    //! solver itself. For a value of 0 no output is generated, will setting
    //! the message level to 1 results in statistical information being printed
    //! to the standard output. The attribute's value is set in the Setup()
    //! method corresponding to the PARDISO_stats parameter.
    int msgLvl_;

    //! Integer value of zero for passing to Pardiso
    int zeroINT_;

    //! Floating point value of zero for passing to Pardiso
    Double zeroDBL_;

    //! Array containing entries of problem matrix

    //! This array contains the entries of the problem matrix of the linear
    //! system that is to be solved. We imnplicitely assume that a Complex can
    //! directly be cast into an array of two Doubles.
    Double *theMatrix_;

    //! A flag specifying if Setup is being called for the first time
    bool firstCall_;

    //! Array with identity reordering

    //! This pointer is either NULL or points to a one-based array containing
    //! an identity re-ordering. The latter is used to force Pardiso to use
    //! the original matrix pattern, in the case that NOREORDERING is specified
    //! by the user via the PARDISO_ordering parameter.
    int *idPerm_;

    //! Size of identity permutation array
    int idPermSize_;

    //! number of non zero entries
    UInt nnz_;

    //! Timer objects
    Timer tNumfact_, tSymfact_;

    //! Check if we have a new matrix pattern
    bool newMatrixPattern_ = false;
  };

}


#endif
