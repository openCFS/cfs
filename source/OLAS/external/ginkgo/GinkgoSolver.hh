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
#include <variant>


namespace CoupledField
{
  class GinkgoSolver : public BaseIterativeSolver
  {
  public:
    GinkgoSolver(PtrParamNode param, PtrParamNode olasInfo, BaseMatrix::EntryType type);

    virtual ~GinkgoSolver() {};

    void Setup(BaseMatrix &sysmat) override;

    /** @param sysmat shall be the one Setup() is called with */
    void Solve( const BaseMatrix &sysmat, const BaseVector &rhs, BaseVector &sol) override;

    BaseSolver::SolverType GetSolverType() override { return BaseSolver::GINKGO; }

    /** to be removed from all solvers! */
    void SetNewMatrixPattern() override {EXCEPTION("SetNewMatrixPattern not implemented for Ginkgo. GetRidOfZeros for NCIs will not work.");};

  private:
    /** see CFS_Solvers.xsd or Ginkgo for a description of the types */
    typedef enum {NOSOLVER = 0, ST_JSON, CG, BICGSTAB, GMRES} GinkgoSolverType;
    static Enum<GinkgoSolverType> ginkgoSolverType;

    typedef enum {NOPRECOND = 0, PT_JSON, JACOBI, ILU, IC, AMG} GinkgoPrecondType;
    static Enum<GinkgoPrecondType> ginkgoPrecondType;

    typedef enum {ABSOLUTE = 0, INITIAL_RESNORM, RHS_NORM} TolType;
    static Enum<TolType> tolType;

    typedef enum {ON = 0, OFF, AUTO} CudaChoice;
    static Enum<CudaChoice> cudaChoice;

    /** Ginkgo cannot make use of only half given symmetric matrices.
     * CFS_T is either double or complex<double>  
     * GK_T is either float/complex<float>, or double/complex<double>.
     * GK_REAL_T is either float or double. */
    template <typename CFS_T, typename GK_T, typename GK_REAL_T>
    void Setup(CRS_Matrix<CFS_T>* mat);

    template <typename CFS_T, typename GK_T, typename GK_REAL_T>
    void Solve(const BaseVector &rhs, BaseVector &sol);

    /** either OmpExecutor::create(), ReferenceExecutor::create() or the HIP (CUDA) variant in future */
    shared_ptr<gko::Executor> exec;
    std::variant<shared_ptr<gko::log::Convergence<float>>,shared_ptr<gko::log::Convergence<double>>> logger;
    shared_ptr<gko::LinOpFactory> precond;
    shared_ptr<gko::LinOp> solver;

    /** for single precision we need a csr values array copy and cast from the double precision cfs matrix.
        The double precision is not need for runtime but is there to avoid static asserts from std::get */
    std::variant<std::vector<float>,std::vector<std::complex<float>>,std::vector<double>,std::vector<std::complex<double>>> values;

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

    /** ignore given solution on Solve() as initial guess */
    bool initial_zero = false;

    /** solve single precision system */
    bool fp32 = false;
  };
}
#endif // GINKGOSOLVER_HH_
