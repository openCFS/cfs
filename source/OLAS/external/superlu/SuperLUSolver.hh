#ifndef SUPERLU_SOLVER_HH
#define SUPERLU_SOLVER_HH

#include <string>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Utils/StdVector.hh"
#include "OLAS/solver/BaseSolver.hh"

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
    void Setup( BaseMatrix &sysmat, PtrParamNode analysis_id);

    //! Direct solution of the linear system

    //! After Solve is called the matrix (which has already to be factorised
    //! by a call of Setup) is finally solved by backward-forward substitution.
    void Solve( const BaseMatrix &sysmat,
                const BaseVector &rhs,
                BaseVector &sol,
                PtrParamNode analysis_id );

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
  };

}

#ifndef EXPLICIT_TEMPLATE_INSTANTIATION
//#include "SuperLUSolver.cc"
#endif

#endif
