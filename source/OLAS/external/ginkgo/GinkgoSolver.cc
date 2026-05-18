#include "GinkgoSolver.hh"
// the header in extensions is not shipped with ginkgo.hpp
#include <ginkgo/extensions/config/json_config.hpp>
#include "Domain/Domain.hh"
#include "Driver/BaseDriver.hh"
#include "MatVec/SparseOLASMatrix.hh"
#include "DataInOut/ProgramOptions.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include <string>
#include <sstream>
#include <stdio.h>
#include <filesystem>
#include <boost/predef.h>

using std::string;
using std::to_string;
using std::complex;
using std::vector;

Enum<GinkgoSolver::GinkgoSolverType>  GinkgoSolver::ginkgoSolverType;
Enum<GinkgoSolver::GinkgoPrecondType> GinkgoSolver::ginkgoPrecondType;
Enum<GinkgoSolver::TolType>           GinkgoSolver::tolType;
Enum<GinkgoSolver::CudaChoice>        GinkgoSolver::cudaChoice;

namespace CoupledField
{
DEFINE_LOG(ginkgo, "ginkgo")

GinkgoSolver::GinkgoSolver(PtrParamNode pn, PtrParamNode olasInfo, BaseMatrix::EntryType type)
{
  ginkgoPrecondType.SetName("GinkgoSolver::GinkgoPrecondType");
  ginkgoPrecondType.Add(NOPRECOND, "noprecond");
  ginkgoPrecondType.Add(AMG, "amg");
  ginkgoPrecondType.Add(JACOBI, "jacobi");
  ginkgoPrecondType.Add(ILU, "ilu");
  ginkgoPrecondType.Add(IC, "ic");
  ginkgoPrecondType.Add(PT_JSON, "json");

  ginkgoSolverType.SetName("GinkgoSolver::GinkgoSolverType");
  ginkgoSolverType.Add(NOSOLVER, "nosolver");
  ginkgoSolverType.Add(ST_AUTO, "auto");
  ginkgoSolverType.Add(CG, "cg");
  ginkgoSolverType.Add(BICGSTAB, "bicgstab");
  ginkgoSolverType.Add(GMRES, "gmres");
  ginkgoSolverType.Add(ST_JSON, "json");

  tolType.SetName("GinkgoSolver::TolType");
  tolType.Add(ABSOLUTE, "absolute");
  tolType.Add(INITIAL_RESNORM, "initial_resnorm");
  tolType.Add(RHS_NORM, "rhs_norm");

  cudaChoice.SetName("GinkgoSolver::CudaChoice");
  cudaChoice.Add(ON, "on");
  cudaChoice.Add(OFF, "off");
  cudaChoice.Add(AUTO, "auto");

  infoNode_ =  olasInfo->Get("ginko");
  xml_ = pn;

  // solver and precond might be overwritten by "json"
  solver_type = ginkgoSolverType.Parse(pn->Get("solver")->As<string>());
  precond_type = ginkgoPrecondType.Parse(pn->Get("precond")->As<string>());
  
  if(solver_type == ST_AUTO)
    solver_type = (type == BaseMatrix::COMPLEX) ? GMRES : CG;

  if(type == BaseMatrix::COMPLEX && solver_type == CG) 
    infoNode_->SetWarning("'cg' solver within ginkgo for complex systems will fail. Try e.g. 'gmres'");

  max_iter = pn->Get("maxIter")->As<double>();
  tolerance = pn->Get("tolerance")->As<double>();
  min_tol = pn->Get("minimalTolerance")->As<double>();
  fp32 = pn->Get("fp32")->As<bool>();

  tol_type = tolType.Parse(pn->Get("toleranceType")->As<string>());
  switch(tol_type)
  {
  case ABSOLUTE:
    tol_mode = gko::stop::mode::absolute;
    break;
  case INITIAL_RESNORM:
    tol_mode = gko::stop::mode::initial_resnorm;
    break;
  case RHS_NORM:
    tol_mode = gko::stop::mode::rhs_norm;
    break;
  }

  initial_zero = pn->Get("initial_zero")->As<bool>();

  // do we use an external json configuration file? Not meant for daily use but e.g. to test what to add in cfs via xml.
  if(pn->Has("json"))
  {
    json = pn->Get("json")->As<string>();
    if(!std::filesystem::exists(json)) // now using std::filesystem per user request
      throw Exception("cannot open json file '" + json + "'");
    infoNode_->Get("json")->SetValue(json);
    solver_type = ST_JSON;
    precond_type = PT_JSON;
  }

  CudaChoice cuda = cudaChoice.Parse(pn->Get("cuda")->As<string>());

  if(UseCuda() && (cuda == ON || cuda == AUTO))
  {
    int cdi = pn->Get("cuda_devide_id")->As<int>();
    exec = gko::CudaExecutor::create(cdi, UseOpenMP() ? gko::OmpExecutor::create() : gko::ReferenceExecutor::create());
    infoNode_->Get("exec")->SetValue("cuda");
  }
  if(cuda == ON && !UseCuda())
    throw Exception("cfs is compiled w/o cuda. Set ginkgo/cudo to 'off' or 'auto'");
  if(!exec && UseOpenMP())
  {
    exec = gko::OmpExecutor::create();
    infoNode_->Get("exec")->SetValue("openmp");
  }
  if(!exec)
  {
    exec = gko::ReferenceExecutor::create();
    infoNode_->Get("exec")->SetValue("reference");
  }
  infoNode_->Get("cuda_available")->SetValue(UseCuda());
  if(pn->Get("dump_profiling")->As<bool>())
  {
#if BOOST_OS_LINUX
    infoNode_->SetWarning("'dump_profiling' currently not possible on Linux"); // we get a vtune linker error
#else
    exec->add_logger(gko::log::ProfilerHook::create_nested_summary());
#endif
  }
}

void GinkgoSolver::Setup(BaseMatrix &sysmat)
{
  LOG_DBG(ginkgo) << "Setup: e=" << MatrixEntryTypeEnum.ToString(sysmat.GetEntryType()) << " s=" << MatrixStorageTypeEnum.ToString(sysmat.GetStorageType())
                  << " r=" << sysmat.GetNumRows() << " c=" << sysmat.GetNumCols();

  if(sysmat.GetStorageType() != BaseMatrix::SPARSE_NONSYM)
    throw Exception("Ginkgo solver can only handle non-symmetric matrix storage");

  if(sysmat.GetEntryType() == BaseMatrix::DOUBLE)
  {
    auto mat = dynamic_cast<CRS_Matrix<double>*>(&sysmat);
    if(fp32) {
      values = vector<float>(); // we use this data to have the same templated csr values data in Setup() and Solve()
      Setup<double, float, float>(mat);
    }
    else
      Setup<double, double, double>(mat);
  }
  else
  {
    auto mat = dynamic_cast<CRS_Matrix<complex<double>>*>(&sysmat);
    if(fp32) {
      values = vector<complex<float>>();
      Setup<complex<double>,complex<float>,float>(mat);
    }
    else
      Setup<complex<double>,complex<double>,double>(mat);
  }
}

template<typename CFS_T, typename GK_T, typename GK_REAL_T>
void GinkgoSolver::Setup(CRS_Matrix<CFS_T>* m)
{
  assert(m != nullptr);

  // Ginkgo knows signed 32 and 64 as idx, but not unsigned
  int nnz  = (int) m->GetNnz();
  int nrow = (int) m->GetNumRows();
  int ncol = (int) m->GetNumCols();

  // nice that we need no temporary data for the matrix, when we are of same type
  GK_T* vdp = nullptr;
  // for 32fp we need a storage of the csr values data which also holds for Solve().
  if(std::is_same<CFS_T, GK_T>::value)
    vdp =  (GK_T*) m->GetDataPointer(); // the cast is only to please the compiler, not relevant in reality
  else
  {
    auto& vdv = std::get<vector<GK_T>>(values); // throws bad_cast exception
    vdv.resize(nnz);
    for(unsigned int i = 0; i < (unsigned int) nnz; i++)
      vdv[i] = (GK_T) m->GetDataPointer()[i]; // from fp62 to fp32
    vdp = (GK_T*) vdv.data();

  }
  auto vv = gko::make_const_array_view(exec, nnz, vdp);
  auto cv = gko::make_const_array_view(exec, nnz, (int*) m->GetColPointer());
  auto rv = gko::make_const_array_view(exec, nrow + 1, (int*) m->GetRowPointer());

  // the std::move makes sure we don't call a copy constructor
  auto csr = gko::share(gko::matrix::Csr<GK_T, int>::create_const(exec, gko::dim<2>(nrow,ncol), std::move(vv), std::move(cv), std::move(rv)));
  logger = gko::share(gko::log::Convergence<GK_T>::create());
  shared_ptr<gko::LinOpFactory> factory;

  // we set this information every time as we handle the options here in detail
  PtrParamNode ip = infoNode_->Get("precond");
  ip->Get("type")->SetValue(ginkgoPrecondType.ToString(precond_type)); // automatically set to json

  PtrParamNode is = infoNode_->Get("solver");
  is->Get("type")->SetValue(ginkgoSolverType.ToString(solver_type)); // json in case

  // a json file contains a arbitrary ginkgo configuration, more versatile than the cfs xml configuration subset
  if(!json.empty())
  {
    // this code is copied from file-config-solver.cpp which also has some sample json config files
    auto config = gko::ext::config::parse_json_file(json);
    auto reg = gko::config::registry();
    auto td = gko::config::make_type_descriptor<GK_T, int>();
    factory = gko::config::parse(config, reg, td).on(exec);
  }
  else
  {
    switch(precond_type)
    {
    case JACOBI:
    {
      unsigned int mbs = 1;
      if(xml_->Has("jacobi"))
        mbs = xml_->Get("jacobi/max_block_size")->As<unsigned int>();
      ip->Get("max_block_size")->SetValue(mbs);
      precond = gko::preconditioner::Jacobi<GK_T>::build().with_max_block_size(mbs).on(exec);
      break;
    }
    case ILU:
      precond = gko::preconditioner::Ilu<gko::solver::LowerTrs<GK_T, int>>::build().on(exec);
      break;
    case IC:
      // there is also ParIct with fill in factor limit, however this is extremely slow. Experiment via json
      // gko::factorization::ParIct<T, int>::build().with_fill_in_limit(limit).with_iterations(iter).on(exec);
      precond = gko::preconditioner::Ic<gko::solver::LowerTrs<GK_T, int>>::build().on(exec);
      break;
    case AMG:
    {
      gko::size_type iter = 1;
      if(xml_->Has("amg"))
        iter  = xml_->Get("amg/iterations")->As<unsigned int>();
      ip->Get("iterations")->SetValue((int) iter);

      precond = gko::solver::Multigrid::build()
                  .with_mg_level(gko::multigrid::Pgm<GK_T, int>::build().with_deterministic(true))
                  .with_criteria(gko::stop::Iteration::build().with_max_iters(iter)).on(exec);
      break;
    }
    case PT_JSON:
    case NOPRECOND:
      assert(false);
      break;
    }

    auto iter_crit = gko::stop::Iteration::build().with_max_iters(max_iter);
    const gko::remove_complex<GK_T> tol = tolerance;
    auto norm_crit =  gko::stop::ResidualNorm<GK_T>::build().with_baseline(tol_mode).with_reduction_factor(tol);
    is->Get("max_iter")->SetValue(max_iter);
    is->Get("tolerance")->SetValue(tolerance);
    is->Get("mode")->SetValue(tolType.ToString(tol_type));

    is->Get("cfs")->SetValue(boost::typeindex::type_id<CFS_T>().pretty_name());
    is->Get("ginkgo")->SetValue(boost::typeindex::type_id<GK_T>().pretty_name());
    is->Get("precision")->SetValue(fp32 ? 32 : 64);

    switch(solver_type)
    {
    case CG:
      factory = gko::solver::Cg<GK_T>::build().with_criteria(iter_crit.on(exec),norm_crit.on(exec)).with_preconditioner(precond).on(exec);
      break;
    case BICGSTAB:
      factory = gko::solver::Bicgstab<GK_T>::build().with_criteria(iter_crit.on(exec),norm_crit.on(exec)).with_preconditioner(precond).on(exec);
      break;
    case GMRES:
      factory = gko::solver::Gmres<GK_T>::build().with_criteria(iter_crit.on(exec),norm_crit.on(exec)).with_preconditioner(precond).on(exec);
      break;
    case NOSOLVER:
    case ST_AUTO:
    case ST_JSON:
      assert(false);
      break;
    }
  }

  solver = factory->generate(csr);
  solver->add_logger(std::get<shared_ptr<gko::log::Convergence<GK_T>>>(logger)); // C++ is so ugly :()
}

void GinkgoSolver::Solve(const BaseMatrix &sysmat, const BaseVector &rhs, BaseVector &sol)
{
  if(sysmat.GetEntryType() == BaseMatrix::DOUBLE)
  {
    if(fp32)
      Solve<double, float, float>(rhs, sol);
    else
      Solve<double, double, double>(rhs, sol);
  }
  else
  {
    if(fp32)
      Solve<complex<double>,complex<float>,float>(rhs, sol);
    else
      Solve<complex<double>,complex<double>,double>(rhs, sol);
  }
}

template <typename CFS_T, typename GK_T, typename GK_REAL_T>
void GinkgoSolver::Solve(const BaseVector &rhs, BaseVector &sol)
{
  // the vectors seem to be possibly based on blocks and therefore are not necessarily just Vector

  // create a continuous vector and fill it
  vector<GK_T> rhs_tmp(rhs.GetSize());
  CFS_T val;
  for(unsigned int i = 0; i < rhs.GetSize(); i++)
  {
    rhs.GetEntry(i, val);
    rhs_tmp[i] = (GK_T) val;
  }

  // make a ginkgo right hand side view
  auto rsv = gko::make_array_view(exec, (int) rhs_tmp.size(), rhs_tmp.data()); // again, ginkgo has signed indices
  // make the ginkgo rhs based on the view
  auto b = gko::matrix::Dense<GK_T>::create(exec, gko::dim<2>((int) rhs.GetSize(), 1), rsv, 1);

  // our solution. cfs provides an initial guess. 
  // in the cuda case, we create the data on the CPU to be able to set it. It seems to be automatically copied to GPU when solving
  auto x = gko::matrix::Dense<GK_T>::create(exec->get_master(), gko::dim<2>(rhs.GetSize(), 1));
  for(unsigned int i = 0; i < rhs.GetSize(); i++)
  {
    sol.GetEntry(i,val);
    x->at(i) = initial_zero ? 0.0 : (GK_T) val;
  }
  double initial_sol_norm = initial_zero ? 0.0 : (double) sol.NormL2();

  // solve S x = b for x
  solver->apply(b, x);

  // store solution slowly
  assert(x->get_num_stored_elements() == sol.GetSize());
  for(unsigned int i = 0; i < sol.GetSize(); i++)
    sol.SetEntry(i, (CFS_T) x->at(i));

  assert(solver->get_loggers().size() == 1);
  auto mylogger = std::get<shared_ptr<gko::log::Convergence<GK_T>>>(logger);
  bool conv = mylogger->has_converged();
  size_t iters = mylogger->get_num_iterations();
  // in the cuda case we live on the GPU and need to copy to CPU
  auto grn = gko::as<gko::matrix::Dense<GK_REAL_T>>(mylogger->get_residual_norm());
  auto grn_host = gko::clone(exec->get_master(), grn);
  double res_norm = (double) ((complex<GK_REAL_T>) grn_host->at(0,0)).real();
  LOG_DBG(ginkgo) << "Solve: converged=" << conv << " iters=" << iters << " residual=" << res_norm << " sol=" << sol.NormL2() << " inital=" << initial_sol_norm;

  ParamNode::ActionType at = progOpts->DoDetailedInfo() ? ParamNode::APPEND : ParamNode::DEFAULT;
  auto in = infoNode_->Get("solve", at);
  in->Get("analysis_id")->SetValue(domain->GetDriver()->GetAnalysisId().ToString());
  in->Get("converged")->SetValue(conv);
  in->Get("iters")->SetValue((unsigned int) iters);
  in->Get("residual")->SetValue(res_norm);
  in->Get("sol")->SetValue(sol.NormL2());
  in->Get("initial")->SetValue(initial_sol_norm);
  if(!json.empty())
    in->Get("json")->SetValue("conv criteria set in json");

  if(!conv)
  {
    string msg = "no convergence within " + to_string(iters) + " iterations with residual " + to_string(res_norm);
    if(!json.empty())
      msg += " (set criteria within you json file)";
    if(res_norm < std::max(tolerance, min_tol))
      infoNode_->SetWarning(msg + " but within minimal tolerance " + to_string(std::max(tolerance, min_tol)));
    else
      throw Exception("ginkgo: " + msg);
  }
}

} // end of namespace


