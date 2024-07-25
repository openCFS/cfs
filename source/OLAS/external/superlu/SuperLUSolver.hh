#ifndef SUPERLU_SOLVER_HH
#define SUPERLU_SOLVER_HH

#include <string>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Utils/StdVector.hh"
#include "OLAS/solver/BaseSolver.hh"

// SuperLU header
#include "slu_ddefs.h"

// SuperLU_MT header
// #include "pcsp_defs.h"

namespace CoupledField {

  template<typename T>
  class SuperLUSolver : public BaseDirectSolver {

  public:

    //! Constructor

    //! The constructor does not do anything but set the pointers to the
    //! internal communciation objects.
    //! Note that the pointer to the report object is optional.
    SuperLUSolver (PtrParamNode solverNode, PtrParamNode olasInfo);

    //! Default Destructor

    //! The destructor is not used at all.
    ~SuperLUSolver();

    //! Reordering and factorization of the linear system

    //! After Setup is called the BaseMatrix (expected to be either a CRS
    //! or SCRS-matrix) will be reordered using the Nested Dissection or
    //! the Minimum Degree Algorithm and then it will be LU-factorised.
    void Setup( BaseMatrix &sysmat);

    //! Dummy method: Notify the solver that a new matrix pattern has been set
    void SetNewMatrixPattern() {EXCEPTION("SetNewMatrixPattern not implemented for SuperLU. GetRidOfZeros for NCIs will not work.");};

    //! Direct solution of the linear system

    //! After Solve is called the matrix (which has already to be factorised
    //! by a call of Setup) is finally solved by backward-forward substitution.
    void Solve( const BaseMatrix &sysmat,
                const BaseVector &rhs,
                BaseVector &sol);

    //! Query type of this solver.

    //! This method can be used to query the type of this solver. The answer
    //! is encoded as a value of the enumeration data type SolverType.
    //! \return PARDISO
    SolverType GetSolverType() {
      return SUPERLU;
    }

  private:

    //! Default constructor

    //! The default constructor is private, since we want the solver object to
    //! be initialised with a pointer to a parameter object right at
    //! instantiation.
    SuperLUSolver();

    //! Stored information about the storage type and entry type of the matrix
    BaseMatrix::StorageType stype = BaseMatrix::NOSTORAGETYPE;
    BaseMatrix::EntryType etype = BaseMatrix::NOENTRYTYPE;

    //! Dimension of the linear system
    int probDim_ = -1;

    //! The number of right hand sides Pardiso should solve the system for
    //! at one pass
    int nrhs_ = 0;

    //! Array containing entries of problem matrix

    //! A flag specifying if Setup is being called for the first time
    bool firstCall_ = true;

    //! Do we solve a system with complex or double entries?
    bool isComplex_ = false;

    //! number of non zero entries
    UInt nnz_ = 0;

    SuperMatrix    A, L, U;
    SuperMatrix    B, X;

    char           equed[1];
    yes_no_t       equil;
    trans_t        trans;
    NCformat       *Astore = nullptr;
    NCformat       *Ustore;
    SCformat       *Lstore;
    GlobalLU_t     Glu; // facilitate multiple factorizations with SamePattern_SameRowPerm
    double         *a;
    int            *asub;
    int            *xa;
    int            *perm_c; /* column permutation vector */
    int            *perm_r; /* row permutations from partial pivoting */
    int            *etree;
    void           *work;
    int            info, lwork, nrhs, ldx;
    int            i = -1;
    double         *rhsb;
    double         *rhsx;
    double         *xact;
    double         *R;
    double         *C;
    double         *ferr;
    double         *berr;
    double         u, rpg, rcond;
    mem_usage_t    mem_usage;
    superlu_options_t options;
    SuperLUStat_t stat;

  };

}

#endif
