#ifndef GINKGOSOLVER_HH_
#define GINKGOSOLVER_HH_

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-local-typedefs"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#pragma GCC diagnostic ignored "-Wsign-compare"
#include "ginkgo/ginkgo.hpp"
#pragma GCC diagnostic pop
#pragma clang diagnostic pop

#include "General/Environment.hh"
#include "OLAS/solver/BaseSolver.hh"
#include "OLAS/solver/BaseSolver.hh"
#include "MatVec/CRS_Matrix.hh"
#include "General/Enum.hh"


namespace CoupledField
{
  class GinkgoSolver : public BaseIterativeSolver
  {
  public:
    /** see CFS_Solvers.xsd or Ginkgo for a description of the types */
    typedef enum {NOSOLVER = 0, ST_JSON, CG, BICGSTAB, GMRES} GinkgoSolverType;
    static Enum<GinkgoSolverType> ginkgoSolverType;

    typedef enum {NOPRECOND = 0, PT_JSON, JACOBI, ILU, IC, AMG} GinkgoPrecondType;
    static Enum<GinkgoPrecondType> ginkgoPrecondType;

    typedef enum {ABSOLUTE = 0, INITIAL_RESNORM, RHS_NORM} TolType;
    static Enum<TolType> tolType;

    GinkgoSolver(PtrParamNode param, PtrParamNode olasInfo, BaseMatrix::EntryType type);

    virtual ~GinkgoSolver() {};

    void Setup(BaseMatrix &sysmat) override;

    /** @param sysmat shall be the one Setup() is called with */
    void Solve( const BaseMatrix &sysmat, const BaseVector &rhs, BaseVector &sol) override;

    BaseSolver::SolverType GetSolverType() override { return BaseSolver::GINKGO; }

    /** to be removed from all solvers! */
    void SetNewMatrixPattern() override {EXCEPTION("SetNewMatrixPattern not implemented for Ginkgo. GetRidOfZeros for NCIs will not work.");};

  private:
    /** Ginkgo cannot make use of only half given symmetric matrices */
    template <class T>
    void Setup(CRS_Matrix<T>* mat);

    template <class T>
    void Solve(const BaseVector &rhs, BaseVector &sol);

    /** either OmpExecutor::create(), ReferenceExecutor::create() or the HIP (CUDA) variant in future */
    std::shared_ptr<const gko::Executor> exec;
    std::shared_ptr<gko::log::Convergence<double>> logger;
    std::shared_ptr<gko::LinOpFactory> precond;
    std::shared_ptr<gko::LinOp> solver;

    GinkgoSolverType solver_type = NOSOLVER;
    GinkgoPrecondType precond_type = NOPRECOND;
    TolType tol_type = ABSOLUTE;
    gko::stop::mode tol_mode = gko::stop::mode::absolute;

    unsigned int max_iter = 10000;
    /** given to the solver */
    double tolerance = 1e-9;
    /** abort only when > max(tolerance, min_tol) */
    double min_tol = 1e-9;

    /** json file name for precond/solver description via external file */
    std::string json;

    bool initial_zero = false;
  };
}
#endif // GINKGOSOLVER_HH_
