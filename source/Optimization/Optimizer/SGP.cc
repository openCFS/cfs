#include "SGP.hh"

#include <math.h>
#include <iterator>
#include <iostream>
#include <vector>
#include <algorithm>
#include <memory>

#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ProgramOptions.hh"
#include "Domain/CoefFunction/CoefFunctionOpt.hh"
#include "Driver/Assemble.hh"
#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Forms/BiLinForms/BDBInt.hh"
#include "MatVec/SBM_Matrix.hh"
#include "MatVec/Matrix.hh"
#include "MatVec/Vector.hh"
#include "OLAS/algsys/AlgebraicSys.hh"
#include "OLAS/solver/BaseSolver.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Excitation.hh"
#include "Optimization/Optimization.hh"
#include "Optimization/OptimizationMaterial.hh"
#include "Optimization/TransferFunction.hh"
#include "PDE/StdPDE.hh"
#include "Utils/tools.hh"
#include "Utils/StdVector.hh"

DEFINE_LOG(sgpopt, "sgp")

using namespace CoupledField;

Enum<SGP::Approximation>     SGP::approximation;
Enum<SGP::Filtering>     SGP::filtering;

SGP::SGP(Optimization* optimization, BaseOptimizer* base, PtrParamNode pn, sgp::InitParams initargs) : sgp::Algorithm(initargs)
{
  LOG_TRACE(sgpopt) << "Initialize SGP";

  this->optimization_ = optimization;
  this->base_ = base;
  this->base_pn_ = pn;

  space_ = optimization_->GetDesign();
  grid_ = domain->GetGrid();
  n_elems_ = space_->GetNumberOfElements();
  n_designs_ = space_->GetNumberOfVariables();
  dim_tens_ = (grid_->GetDim() == 3) ? 6 : 3;

  // how many design variables per FE?
  n_vars_ = space_->design.GetSize();

  ErsatzMaterial* em = dynamic_cast<ErsatzMaterial*>(domain->GetOptimization());
  // if filtering, make sure that use_mat_filt is on
  if (em->pn->Has("filters")) {
    if (!space_->is_matrix_filt)
      EXCEPTION("Turn 'use_mat_filt' on for external SGP lib!");
  }

  // this is the input xml
  this->optimizer_pn_ = pn->Get(Optimization::optimizer.ToString(Optimization::SGP_SOLVER), ParamNode::PASS);
  // output
  this->sgp_node_ = optimization->optInfoNode->Get("optimizer/sgp", ParamNode::APPEND);

  // check for optional parameters
  if(optimizer_pn_ != NULL)
  {
    ParamNodeList list;

    list = optimizer_pn_->GetListByVal("option", "type", "string");
    assert(list.GetSize() > 0);
    for(unsigned int i = 0; i < list.GetSize(); i++)
      SetStringValue(list[i]->Get("name")->As<std::string>(), list[i]->Get("value")->As<std::string>());

    list = optimizer_pn_->GetListByVal("option", "type", "integer");
    assert(list.GetSize() > 0);
    for(unsigned int i = 0; i < list.GetSize(); i++)
      SetIntegerValue(list[i]->Get("name")->As<std::string>(), list[i]->Get("value")->As<int>());

    list = optimizer_pn_->GetListByVal("option", "type", "real");
    assert(list.GetSize() > 0);
    for(unsigned int i = 0; i < list.GetSize(); i++)
      SetNumericValue(list[i]->Get("name")->As<std::string>(), list[i]->Get("value")->As<double>());

    gradientCheck_ = optimizer_pn_->Get("gradientCheck", ParamNode::EX)->As<bool>();

    solve_params_.precompute_tensors = optimizer_pn_->Get("precomputeTensors", ParamNode::EX)->As<bool>();
    generatePlotData_ = optimizer_pn_->Get("generatePlotData", ParamNode::EX)->As<bool>();
  }

  double manual_scaling = optimizer_pn_ != NULL && optimizer_pn_->Has("obj_scaling_factor") ?
        optimizer_pn_->Get("obj_scaling_factor")->Get("value")->As<Double>() : 1.0;
  base_->PostInitScale(manual_scaling);

  SGP::SetEnums();
}

SGP::~SGP()
{
}

void SGP::Init()
{
  // gradient check doesn't need this stuff
  if (gradientCheck_)
    return;

  // tell external SGP lib which model to use
  if (approximation_ == ZERO_ASYMPTOTES)
    solve_params_.model = sgp::F_PHYS::SIMPLIFIED;
  else {
    assert(approximation_ == EXACT);
    solve_params_.model = sgp::F_PHYS::GENERALIZED;
  }

  if (!domain->GetGrid()->IsGridRegular())
    solve_params_.n_ip = 3; // TODO: use correct number of integration points

  // make sure that exact approximation is only called for a regular mesh
  // as the number of integration points per element is set to a constant
  if (domain->GetGrid()->IsGridRegular()) {
    BaseBDBInt* bdb = dynamic_cast<BaseBDBInt*>(Optimization::context->pde->GetAssemble()->GetBiLinForm("LinElastInt", space_->GetRegionId(), Optimization::context->pde)->GetIntegrator());
    assert(bdb != NULL);

    BaseFE* ptFe = NULL;
    StdVector<LocPoint> intPoints;
    StdVector<double> weights;
    // take first element in design region
    GetIntPoints(bdb, space_->data[0].elem, &ptFe, intPoints, weights);
    assert(!intPoints.IsEmpty());
    nip_ = intPoints.GetSize();
    solve_params_.n_ip = nip_;
  }

  solve_params_.info           = 1;     // print info in each iteration
  solve_params_.iter_max       = optimization_->GetMaxIterations();     // maximum number of outer iterations
  solve_params_.iter_bisec_max = max_bisections_;     // maximum number of inner iterations
  solve_params_.iter_glob_max = max_globs_;     // maximum number of globalization iterations
  solve_params_.J_stop         = tolerance_; // stopping criterion for cost function
  solve_params_.eps            = volume_tolerance_; // stopping criterion for volume bisection
  solve_params_.p_filt_rho     = p_filt_dens_;   // filter penalty
  solve_params_.p_filt_phi     = p_filt_angle_;   // filter penalty
  solve_params_.tau            = tau_init_;   // globalization penalty
  solve_params_.tau_incr       = tau_factor_;   // increment for tau
  solve_params_.lambda_min     = pmin_vol_; // lower starting value for volume bisection
  solve_params_.lambda_max     = pmax_vol_; // upper starting value for volume bisection
  solve_params_.outputfile     = progOpts->GetSimName() + ".sgp.log";
  solve_params_.n_levels       = levels_;
  assert(levels_ >= 1);

  PtrParamNode pnr = sgp_node_->Get("params", ParamNode::APPEND);
  pnr->Get("max_iter")->SetValue(solve_params_.iter_max);
  pnr->Get("max_bisecs")->SetValue(solve_params_.iter_bisec_max);
  pnr->Get("max_globs")->SetValue(solve_params_.iter_glob_max);
  pnr->Get("J_stop")->SetValue(solve_params_.J_stop, 10);
  pnr->Get("vol_stop")->SetValue(solve_params_.eps, 10);
  pnr->Get("p_filt_density")->SetValue(solve_params_.p_filt_rho);
  pnr->Get("p_filt_angle")->SetValue(solve_params_.p_filt_phi);
  pnr->Get("tau_init")->SetValue(solve_params_.tau);
  pnr->Get("tau_incr")->SetValue(solve_params_.tau_incr);
  pnr->Get("n_ip")->SetValue(solve_params_.n_ip);
  pnr->Get("model")->SetValue(approximation.ToString(approximation_));
  pnr->Get("filtering")->SetValue(filtering.ToString(filtering_));
  pnr->Get("generatePlotData")->SetValue(generatePlotData_);

  // dm->HasParameter is not a valid check as all designs and fixed parameters are stored in 'params_'
  // e.g. setting <param name="rotAngleThird" value="0"/> in xml file will give dm->HasParameter(DesignElement::ROTANGLETHIRD) = true
  // FindDesign() returns -1 if design type cannot be found; second argument of FindDesign() prevents throwing an exception
  hasDensity_ = space_->HasDesign(DesignElement::DENSITY);
  hasRotAngleFirst_ = space_->HasDesign(DesignElement::ROTANGLE) || space_->HasDesign(DesignElement::ROTANGLEFIRST);
  hasRotAngleSecond_ = space_->HasDesign(DesignElement::ROTANGLESECOND);
  hasRotAngleThird_ = grid_->GetDim() == 3 && space_->HasDesign(DesignElement::ROTANGLETHIRD);
  // need to know if filtering for angles should be considered
  hasAngle_ = (hasRotAngleFirst_ && hasRotAngleSecond_) || (hasRotAngleFirst_ && hasRotAngleThird_) || (hasRotAngleSecond_ && hasRotAngleThird_);

  // min. number of samples, e.g. if variable is not used
  int n_samples_default = 1;

  // maximum number of design variables per element
  var_idx_.Reserve(4);
  if (hasDensity_) {
    var_idx_.Push_back(sgp::RHO);
    solve_params_.n_rho = samples_per_level_density_; // number of samples per grid level
  } else
      solve_params_.n_rho = n_samples_default;
  if (hasRotAngleFirst_) {
    var_idx_.Push_back(sgp::PHI);
    solve_params_.n_phi = samples_per_level_angle_; // number of samples per grid level
  } else
      solve_params_.n_phi = n_samples_default;
  if (hasRotAngleSecond_) {
    var_idx_.Push_back(sgp::THETA);
    solve_params_.n_theta = samples_per_level_angle_;
  } else
    solve_params_.n_theta = n_samples_default;
  if (hasRotAngleThird_) {
    var_idx_.Push_back(sgp::PSI);
    solve_params_.n_psi = samples_per_level_angle_;
  } else
    solve_params_.n_psi = n_samples_default;

  std::cout << "samples rho=" << solve_params_.n_rho << " phi=" << solve_params_.n_phi << " theta=" << solve_params_.n_theta << " psi=" << solve_params_.n_psi << std::endl;

  pnr = sgp_node_->Get("hierGrid", ParamNode::APPEND);
  pnr->Get("n_levels")->SetValue(solve_params_.n_levels);
  pnr->Get("samples_rho")->SetValue(solve_params_.n_rho);
  pnr->Get("samples_phi")->SetValue(solve_params_.n_phi);
  pnr->Get("samples_theta")->SetValue(solve_params_.n_theta);
  pnr->Get("samples_psi")->SetValue(solve_params_.n_psi);

  SetInitialDesign();

  // which filtering model: only angle field or combination of sin and cos?
  solve_params_.filt = filtering_ == LINEAR ? sgp::F_REG::ANGLE : sgp::F_REG::SIN_COS;

  SetFilterMatrices();

  solve_params_.p = space_->GetTransferFunction(DesignElement::DENSITY, App::MECH)->GetParam();                    // power law parameter

  SetBounds();

  pnr = sgp_node_->Get("bounds", ParamNode::APPEND);
  // write out lower and upper bound for density
  PtrParamNode pnr_rho =  pnr->Get(DesignElement::type.ToString(DesignElement::DENSITY), ParamNode::APPEND);
  pnr_rho->Get("min")->SetValue(solve_params_.rhoMin);
  pnr_rho->Get("max")->SetValue(solve_params_.rhoMax);
  // write out lower and upper bound for rotAngle
  PtrParamNode pnr_phi =  pnr->Get(DesignElement::type.ToString(DesignElement::ROTANGLEFIRST), ParamNode::APPEND);
  pnr_phi->Get("min")->SetValue(solve_params_.phiMin);
  pnr_phi->Get("max")->SetValue(solve_params_.phiMax);
  // write out lower and upper bound for rotAngleSecond
  PtrParamNode pnr_theta =  pnr->Get(DesignElement::type.ToString(DesignElement::ROTANGLESECOND), ParamNode::APPEND);
  pnr_theta->Get("min")->SetValue(solve_params_.thetaMin);
  pnr_theta->Get("max")->SetValue(solve_params_.thetaMax);
  // write out lower and upper bound for rotAngleThird
  PtrParamNode pnr_psi =  pnr->Get(DesignElement::type.ToString(DesignElement::ROTANGLETHIRD), ParamNode::APPEND);
  pnr_psi->Get("min")->SetValue(solve_params_.psiMin);
  pnr_psi->Get("max")->SetValue(solve_params_.psiMax);

  LOG_DBG2(sgpopt) <<"bounds on rho: min=" << solve_params_.rhoMin << " max=" << solve_params_.rhoMax;
  LOG_DBG2(sgpopt) <<"bounds on phi: min=" << solve_params_.phiMin << " max=" << solve_params_.phiMax;
  LOG_DBG2(sgpopt) <<"bounds on theta: min=" << solve_params_.thetaMin << " max=" << solve_params_.thetaMax;
  LOG_DBG2(sgpopt) <<"bounds on psi: min=" << solve_params_.psiMin << " max=" << solve_params_.psiMax;
}

void SGP::SetFilterMatrices()
{
  if (space_->density_filter.GetSize() > 0)
  {
    if(space_->density_filter.GetSize() > 2)
      WARN("All angles are filtered with the same filter (first one declared in xml)!");
//    if(hasDensity_ && hasAngle_ && space_->density_filter.GetSize() < 2)
//      EXCEPTION("2 filters (one for density + one for angles) are required!");

    CRS_Matrix<double>* filtMatDens = NULL;
    CRS_Matrix<double>* filtMatAngle = NULL;
    for (unsigned int fidx = 0; fidx < space_->density_filter.GetSize(); fidx++) {
      // find correct filter matrix for density and angles
      // assume that same filter is applied to all angles
      DesignElement::Type& dt = space_->density_filter[fidx].designType;
      if (dt == DesignElement::DENSITY)
        filtMatDens = &space_->density_filter[fidx].filter_mat;
      else
        if (dt == DesignElement::ROTANGLE || dt == DesignElement::ROTANGLEFIRST || dt == DesignElement::ROTANGLESECOND || dt == DesignElement::ROTANGLETHIRD)
          filtMatAngle = &space_->density_filter[fidx].filter_mat;
        else
          EXCEPTION("Filter matrix for design '" + DesignElement::type.ToString(dt) + "' unknown! Only density and rotAngles are allowed.");
    }

    // make sure that matrices are set
//    assert(filtMatDens->GetNnz() > 0);
    if (filtMatDens != NULL) {
      // filter matrix for density field
      solve_params_.F_rho = sgp::CSR_Matrix(filtMatDens->GetNumRows(), filtMatDens->GetNumCols(), filtMatDens->GetRowPointer(), filtMatDens->GetColPointer(), filtMatDens->GetDataPointer());
    }

    // filter matrix for local orientations
    if (space_->density_filter.GetSize() > 1) {
      assert(filtMatAngle->GetNnz() > 0);
      solve_params_.F_phi = sgp::CSR_Matrix(filtMatAngle->GetNumRows(), filtMatAngle->GetNumCols(), filtMatAngle->GetRowPointer(), filtMatAngle->GetColPointer(), filtMatAngle->GetDataPointer());
    }
    //    LOG_DBG3(sgpopt) << "filtMat rho: " << filtMatDens.ToString();
    //    LOG_DBG3(sgpopt) << "filtMat angles: " << filtMatAngle.ToString();
    // remove filter cause we don't want cfs to do the filtering
    DisableFilters();
  }
  else
    LOG_DBG3(sgpopt) << "no filtering: filtMat is identity matrix";
}

void SGP::DisableFilters() {
  StdVector<DesignElement>& design = optimization_->GetDesign()->data;
  #pragma omp parallel for num_threads(CFS_NUM_THREADS)
  for(unsigned int i = 0; i < design.GetSize(); i++) {
    DesignElement* de = &design[i];
    unsigned int n_filter = de->simp->filter.GetSize();
    for (unsigned int fidx = 0; fidx < n_filter; fidx++)
      de->simp->filter[fidx].SetType(Filter::NO_FILTERING);
  }
  space_->density_filter.Clear(0);
  space_->is_matrix_filt = false;
}

void SGP::SolveProblem()
{
  PtrParamNode summary = optimization_->optInfoNode->Get(ParamNode::SUMMARY);

  // don't solve subproblems if we only want to validate the gradients
  if (gradientCheck_) {
    GradientCheck();
    return;
  }

  sgp::Algorithm::solve(solve_params_, true, generatePlotData_);

  base_->CommitIteration();
  sgp::Info info = sgp::Algorithm::getStatus();

  ToInfo(summary, info);

//  summary->Get("break/converged")->SetValue(converged ? "yes" : "no");
//  if(iter >= optimization_->GetMaxIterations())
//    summary->Get("reason/msg")->SetValue("maximum iterations reached");
//
//  if(converged)
//    std::cout << "\nSGP converged after " << (iter-1) << " iterations: " << summary->Get("reason/msg")->As<string>() << std::endl;
//  else
//    std::cout << "\nSGP did not converge: " << summary->Get("reason/msg")->As<string>() << std::endl;
}

bool SGP::abort(const sgp::Info& info) const {
  return optimization_->DoStopOptimization();
}

void SGP::store_filt_results(const sgp::Material_Vector& x_filt_sin, const sgp::Material_Vector& x_filt_cos, const sgp::Material_Vector& diff_filt_sin, const sgp::Material_Vector& diff_filt_cos) const
{
  // write to output:
  // - filtered design
  // - difference between filtered and unfiltered design
  // - if non-linear filter is used: Freg = 0.5*p_filt * (||sin(2*phi) - F*sin(2*phi)||^2 + ||cos(2*phi) - F*cos(2*phi)||^2)
  // - else (linear): Freg = 0.5*p_filt * ||phi - F*phi||^2; results are stored in "x_filt_sin" and "diff_filt_sin"
  // - the filter for the density is always of type linear, since sin() and cos() does not make sense there
  //
  int res_idx_filt_dens = space_->GetSpecialResultIndex(DesignElement::DENSITY, DesignElement::FILTERED_DESIGN, DesignElement::NONE, DesignElement::PLAIN);
  int res_idx_diff_dens = space_->GetSpecialResultIndex(DesignElement::DENSITY, DesignElement::DIFF_FILTERED_DESIGN, DesignElement::NONE, DesignElement::PLAIN);

  int res_idx_filt_angle = space_->GetSpecialResultIndex(DesignElement::ROTANGLEFIRST, DesignElement::FILTERED_DESIGN, DesignElement::NONE, DesignElement::PLAIN);
  int res_idx_filt_angle_sin = space_->GetSpecialResultIndex(DesignElement::ROTANGLEFIRST, DesignElement::FILTERED_DESIGN, DesignElement::SIN, DesignElement::PLAIN);
  int res_idx_filt_angle_cos = space_->GetSpecialResultIndex(DesignElement::ROTANGLEFIRST, DesignElement::FILTERED_DESIGN, DesignElement::COS, DesignElement::PLAIN);

  int res_idx_diff_angle = space_->GetSpecialResultIndex(DesignElement::ROTANGLEFIRST, DesignElement::DIFF_FILTERED_DESIGN, DesignElement::NONE, DesignElement::PLAIN);
  int res_idx_diff_angle_sin = space_->GetSpecialResultIndex(DesignElement::ROTANGLEFIRST, DesignElement::DIFF_FILTERED_DESIGN, DesignElement::SIN, DesignElement::PLAIN);
  int res_idx_diff_angle_cos = space_->GetSpecialResultIndex(DesignElement::ROTANGLEFIRST, DesignElement::DIFF_FILTERED_DESIGN, DesignElement::COS, DesignElement::PLAIN);


  assert(n_vars_ >= 1);
  unsigned int idx = 0;
  assert(var_idx_.GetSize() > 0);
  for (int var : var_idx_) {
    for (unsigned int e = 0; e < n_elems_; e++) {
      DesignElement* de = &space_->data[idx];
      idx += 1;
      if (var == sgp::RHO) {
        // x_filt_cos and diff_filt_cos are only used for angle variables and in case of non-linear filter
        if (res_idx_filt_dens > 0)
          de->specialResult[res_idx_filt_dens] = x_filt_sin[var][e];
        if (res_idx_diff_dens > 0)
        de->specialResult[res_idx_diff_dens] = diff_filt_sin[var][e];
      }
      if (var == sgp::PHI){
        if (filtering_ == NON_LINEAR) {
          if (res_idx_filt_angle_sin > 0)
            de->specialResult[res_idx_filt_angle_sin] = x_filt_sin[var][e];
          if (res_idx_diff_angle_sin > 0)
            de->specialResult[res_idx_diff_angle_sin] = diff_filt_sin[var][e];
          if (res_idx_filt_angle_cos > 0)
            de->specialResult[res_idx_filt_angle_cos] = x_filt_cos[var][e];
          if (res_idx_diff_angle_cos > 0)
            de->specialResult[res_idx_diff_angle_cos] = diff_filt_cos[var][e];
        } else {
          if (res_idx_filt_angle > 0)
            de->specialResult[res_idx_filt_angle] = x_filt_sin[var][e];
          if (res_idx_diff_angle > 0)
            de->specialResult[res_idx_diff_angle] = diff_filt_sin[var][e];
        }

      }
      //        if (var == sgp::RHO || var == sgp::PHI)
        //          std::cout <<  "var=" << var << " elem=" << e << " value=" << val << std::endl;
    }
  }
}

double SGP::eval_J(const sgp::Material_Vector& x) const
{
  std::vector<double> design;

  if (hasDensity_)
    std::copy(x[sgp::RHO].begin(), x[sgp::RHO].end(), std::back_inserter(design));

  if (hasRotAngleFirst_)
    // 1st rotation angle (default: about z-axis)
    std::copy(x[sgp::PHI].begin(), x[sgp::PHI].end(), std::back_inserter(design));

  if (grid_->GetDim() == 3) {
    if (hasRotAngleSecond_)
      // 2nd rotation angle (default: about y-axis)
      std::copy(x[sgp::THETA].begin(), x[sgp::THETA].end(), std::back_inserter(design));
    if (hasRotAngleThird_)
      // 3rd rotation angle (default: about x-axis)
      std::copy(x[sgp::PSI].begin(), x[sgp::PSI].end(), std::back_inserter(design));
  }

  LOG_DBG3(sgpopt) << "eval_J: copied design to " << ToString(design.data(), design.size());

  double val = base_->EvalObjective((int)design.size(), design.data(), true);
  StdVector<double> vol(1);
  // only 1 constraint (volume)! all other contraints are treated separately by sgp
  base_->EvalConstraints((int)design.size(), design.data(), hasDensity_ ? 1:0, true, vol.GetPointer(), false);

  return val;
}

void SGP::eval_dJ(const sgp::Material_Vector& x, std::vector<sgp::Matrix>& dJ) const
{
  StdVector<Matrix<double>> deriv;
  StdVector<Vector<double>> tmp;
  CalcTensDeriv(deriv, tmp);

  assert(dJ.size() == n_elems_);
  unsigned int dim = domain->GetGrid()->GetDim();
  // copy res to data structure of external sgp
  for (unsigned int e = 0; e < n_elems_; e++) {
    // expects dimension of grid
    sgp::Matrix dJ_e(sgp::n_VoigtNotation(dim));
    LOG_DBG3(sgpopt) << "pass tensor derivative of size " << dJ_e.rows() << " x " << dJ_e.cols() << " to external sgp lib";
    LOG_DBG3(sgpopt) << "e=" << e << "; dJ/dD=" << deriv[e].ToString();
    for (int i = 0; i < dJ_e.rows(); i++)
      for (int j = 0; j < dJ_e.cols(); j++)
        dJ_e(i,j) = deriv[e][i][j];

    dJ[e] = dJ_e;
  }

  // do it here and not in eval_J because eval_J will be also called in globalization loop
  base_->CommitIteration();
}

void SGP::eval_Gamma_epsilon(const sgp::Material_Vector& x, std::vector<sgp::Matrix>& Gamma, std::vector<sgp::Vector>& epsilon) const
{
  // number of rows of material tensor (depends on dimension of problem)
  unsigned int D_rows = (domain->GetGrid()->GetDim() == 3) ? 6 : 3;

  StdVector<Matrix<double>> Gamma_tmp;
  StdVector<Vector<double>> eps;
  CalcTensDeriv(Gamma_tmp, eps);

  std::vector<sgp::Matrix> GammaOut;
  std::vector<sgp::Vector> epsOut;

  // make sure output is already initialized
  assert(Gamma.size() == n_elems_);

  // copy result to data structure of external sgp
  for (unsigned int e = 0; e < n_elems_; e++) {

    sgp::Matrix Gamma_e(D_rows*nip_);
    sgp::Vector eps_e(D_rows*nip_);
    for (int i = 0; i < Gamma_e.rows(); i++) {
      for (int j = 0; j < Gamma_e.cols(); j++)
        Gamma_e(i,j) = Gamma_tmp[e][i][j];

      eps_e[i] = eps[e][i];
    }
    Gamma[e] = Gamma_e;
    epsilon[e] = eps_e;
  }

}


// expensive - hopefully only called once!
Matrix<double> SGP::GetInverseStiffMat(const Function* f) const
{
  // get global system matrix
  AlgebraicSys* algSys = dynamic_cast<StdSolveStep*>(f->ctxt->pde->GetSolveStep())->GetAlgSys();
  const StdMatrix* stiffMat = algSys->GetMatrix(SYSTEM)->GetPointer(0,0);
  LOG_DBG3(sgpopt) << "GGE: Got stiffMat with size=" << stiffMat->GetNumRows() << " x " << stiffMat->GetNumCols() << " and structure type: " << stiffMat->GetStructureType();
  LOG_DBG3(sgpopt) << "GGE: K=" << stiffMat->ToString() << std::flush;
  unsigned int ndofs = stiffMat->GetNumRows();
  Matrix<double> res(ndofs, ndofs);
  Vector<double> rhs(ndofs);
  Vector<double> col_res(ndofs);

  BaseSolver* solver = algSys->GetSolver();
  LOG_DBG3(sgpopt) << "CG: Solving " << ndofs << " LSEs with solver " << BaseSolver::solverType.ToString(solver->GetSolverType());
  for (unsigned int dof = 0; dof < ndofs; dof++) {
    rhs.Init();
    rhs(dof) = 1.0;
    solver->Solve(*stiffMat, rhs, col_res);
    res.SetCol(col_res, dof);
  }

  return res;
}

void SGP::CalcTensDeriv(StdVector<Matrix<double>>& deriv, StdVector<Vector<double>>& eps) const
{
  // get objective function
  assert(optimization_->objectives.Has(Function::COMPLIANCE));
  Function* f =  optimization_->objectives.Get(Function::COMPLIANCE);

  assert(!f->elements.IsEmpty());

  // for element e: df/dD = -B_e * u * u^T * B_e^T
  // need to find Bmat for all elements: basically a copy of MagSIMP:CalcMagFluxDensity, see also BDBInt::CalcElementMatrix
  // get the form first
  Context* context = Optimization::context;

  // GetRegionId() is only valid if we have 1 region
  assert(space_->GetRegionIds().GetSize() == 1);
  BaseBDBInt* bdb = dynamic_cast<BaseBDBInt*>(context->pde->GetAssemble()->GetBiLinForm("LinElastInt", space_->GetRegionId(), context->pde)->GetIntegrator());
  assert(bdb != NULL);

  ErsatzMaterial* em = dynamic_cast<ErsatzMaterial*>(domain->GetOptimization());
  StateContainer& forward = em->GetForwardStates();

  unsigned int nelems = space_->GetNumberOfElements();
  LOG_DBG3(sgpopt) << "GdfDH: f=" << f->ToString() << " nelems in design space=" << nelems;

  StdVector<LocPoint> intPoints; // Get integration Points
  LocPointMapped lp;
  StdVector<double> weights;

  // do we have multiple excitations?
  MultipleExcitation* me = optimization_->me;
  bool mult_excite = me->IsEnabled();

  // deriv is the result vector with n_elems_ matrices of size dim_tens_ x dim_tens_
  deriv.Resize(n_elems_);
  Matrix<double> invK;
  if (approximation_ == ZERO_ASYMPTOTES)
    for (unsigned int e = 0; e < nelems; e++) {
      deriv[e].Resize(dim_tens_,dim_tens_);
      deriv[e].Init();
    }
  else {
    assert(approximation_ == EXACT);
    eps.Resize(n_elems_);
    invK = GetInverseStiffMat(f);
  }

  // the multiple excitation case is a special case - for all other cases this is executed once
  for(unsigned int ex = 0; ex < (mult_excite ? me->excitations.GetSize() : 1); ex++)
  {
    Excitation* excite = mult_excite ? &(me->excitations[ex]) : f->GetExcitation();
    assert(excite != NULL);
    assert(excite->index < (int) me->excitations.GetSize());
    // the stored element solution vector
    StdVector<SingleVector*>& sol = forward.Get(excite)->elem[App::MECH];
    for(unsigned int e = 0; e < nelems; e++)
    {
      // elemNum is 1-based
      const Elem* elem = grid_->GetElem(e+1);

      Vector<double>& u_e = dynamic_cast<Vector<double>&>(*sol[e]);
      // precompute u*u^T - outer product gives 8 x 8 matrix
      // need to convert u to a matrix and let matrix class do the multiplication
      Matrix<double> umat;
      umat.Assign(u_e, u_e.GetSize(), 1, true);
      // transposed u_e
      Matrix<double> umatT;
      umat.Transpose(umatT);

      umat *= umatT;

      LOG_DBG3(sgpopt) << "GdfDH: elem=" << elem->elemNum << " u_e=" << u_e.ToString();
      LOG_DBG3(sgpopt) << " mat u_e^T=" << umatT.ToString();
      LOG_DBG3(sgpopt) << " u_e*u_e^T=" << umat.ToString();

      BaseFE* ptFe = NULL;

      // get integration points for current element
      GetIntPoints(bdb, elem, &ptFe, intPoints, weights);
      assert(!intPoints.IsEmpty());
      assert(ptFe != NULL);

      // Get shape map from grid
      // shared_ptr<ElemShapeMap> esm = grid_->GetElemShapeMap(elem);
      // grid_->UpdateIndividualShapeMap(elem, ptFe);
      // shared_ptr<ElemShapeMap> esm = elem->ptrShapeMap;
      shared_ptr<ElemShapeMap> esm = (elem)->GetElemShapeMap(grid_);
      // for intermediate steps
      Matrix<double> bMatT;
      Matrix<double> bMat;

      // store bMat for all integration points
      StdVector<Matrix<double>> bMat_e;

      for(unsigned int ip = 0; ip < intPoints.GetSize(); ip++)
      {
        // Calculate for each integration point the LocPointMapped
        lp.Set(intPoints[ip], esm, weights[ip]);

        // Call the CalcBMat()-method
        assert(bdb->GetBOp());
        bdb->GetBOp()->CalcOpMat(bMat, lp, ptFe);
        if (approximation_ == ZERO_ASYMPTOTES) {
          bMat.Transpose(bMatT);

          LOG_DBG3(sgpopt) << "GdfDH: elem=" << elem->elemNum << " ip=" << ip << " bMat=" << bMat.ToString();
          // multiplications for derivative of compliance w. r. t. material tensor: B*u*u^T*B^T
          Matrix<double> uuTBT;
          uuTBT.Resize(u_e.GetSize(),dim_tens_);
          // u*u^T*B^T
          umat.Mult(bMatT,uuTBT);
          LOG_DBG3(sgpopt) << "GdfDH: elem=" << elem->elemNum << " ip=" << ip << " u*u^T*B^T=" << uuTBT.ToString();
          // B*u*u^T*B^T
          bMat *= uuTBT;

          deriv[e] -= weights[ip] * lp.jacDet * bMat;

          LOG_DBG3(sgpopt) << "GdfDH: elem=" << elem->elemNum << " ip=" << ip << " B*u*u^T*B^T=" << bMat.ToString();
          LOG_DBG3(sgpopt) << " tmp=" << deriv[e].ToString();
        }
        else {
          assert(approximation_ == EXACT);
          bMat *= sqrt(weights[ip] * lp.jacDet);
          bMat_e.Push_back(bMat);
        }
      } // integration points
      LOG_DBG3(sgpopt) << "GdfDH: elem=" << elem->elemNum << " dfdH=" << deriv[e].ToString();

      if (approximation_ == EXACT) {
        Vector<double> eps_e;
        Matrix<double> gamma_e;
        ComputeElementGammaEps(gamma_e, eps_e, invK, elem, bMat_e, u_e);
        deriv[e] = gamma_e;
        eps[e] = eps_e;
      }

    } // elements
  } // excitation
}

void SGP::ComputeElementGammaEps(Matrix<Double>& gamma, Vector<Double>& eps, const Matrix<double>& invK, const Elem* elem, const StdVector<Matrix<double>> bMat, Vector<double> u) const
{
  // get eqnVec
  // see ErsatzMaterial::ConstructAdjointRHSBuckling
  Assemble* assemble = Optimization::context->pde->GetAssemble();
  RegionIdType regionId = elem->regionId;
  SinglePDE* pde = Optimization::context->pde;
  BiLinFormContext* linElastIntCtxt = assemble->GetBiLinForm("LinElastInt", regionId, pde, pde, false);
  BDBInt<Double,Double>* linElast_bdb = dynamic_cast<BDBInt<Double,Double>*>(linElastIntCtxt->GetIntegrator());

  // get equation numbers of element
  StdVector<Integer> eqnVec;
  FeSpace* feSpace = linElast_bdb->GetFeSpace1();
  feSpace->GetElemEqns(eqnVec, elem);
  LOG_DBG3(sgpopt) << "CG: elem" << elem->elemNum << " eqn vectors:" << eqnVec.ToString();

  // 2D quadriliteral with order 4 integration: 8 x 12 entries
  unsigned int nip = bMat.GetSize();
  // number of free dofs in element e
  unsigned int ndofe = eqnVec.GetSize();

  assert(ndofe == domain->GetGrid()->GetDim()*nip);
  LOG_DBG3(sgpopt) << "elem " << elem->elemNum << " has " << ndofe << " free dofs";
//  stiffMat->Export("stiffMat", BaseMatrix::MATRIX_MARKET);
  Matrix<double> bMat_glob(3*bMat.GetSize(), ndofe);

  // number of rows of material tensor (depends on dimension of problem)
  unsigned int D_rows = (domain->GetGrid()->GetDim() == 3) ? 6 : 3;
  Matrix<double> Bglob_e(ndofe, D_rows*nip);

  assert(eqnVec.GetSize() == ndofe);
  /* building Bglob_e */
  for (unsigned int ip = 0; ip < nip; ip++) {
    // stack B matrix of all integration point of current element e
    Matrix<double> bMat_ip ;
    bMat[ip].Transpose(bMat_ip);
    LOG_DBG3(sgpopt) << "CG: ip=" << ip << " bMat_ip size=" << bMat_ip.GetNumRows() << " x " << bMat_ip.GetNumCols();
    LOG_DBG3(sgpopt) << "CG: bMat_ip=" << bMat_ip.ToString();
    // at each intergration point, copy each entry at row eqn into Bglob
    Bglob_e.SetSubMatrix(bMat_ip, 0, ip*D_rows);
    LOG_DBG3(sgpopt) << "CG: stacking Bmats: " << " ip=" << ip << " bMat_ip:" << bMat_ip.ToString();
    LOG_DBG3(sgpopt) << "Bglob_e eqnNr=" << eqnVec[ip];
  }

  LOG_DBG3(sgpopt) << "CG: Bglob_e = " << Bglob_e.ToString();
  // eps = B^\top*u
  eps.Resize(D_rows*nip);
  Bglob_e.MultT(u, eps);

  LOG_DBG3(sgpopt) << "CG: epsilon = " << std::setprecision(5) << eps.ToString();
  gamma.Resize(D_rows*nip, D_rows*nip);

  ComputeElementGamma(gamma, eqnVec, invK, Bglob_e);
}

void SGP::ComputeElementGamma(Matrix<Double>& gamma, const StdVector<Integer> eqnVec, const Matrix<double>& invK, Matrix<double> Be) const
{
  unsigned int ndofe = Be.GetNumRows();
  assert(ndofe > 0);
  assert(ndofe == eqnVec.GetSize());

  //\Gamma = B_e K^(-1) B_e^\top for each finite element e

  // temporary storage for submatrix of stiffness matrix with rows respective to current element
  Matrix<double> Kinv_e(ndofe,ndofe);
  for (unsigned int eqn_idx = 0; eqn_idx < ndofe; eqn_idx++) {
    int eqnNum = eqnVec[eqn_idx]-1;
    // skip homogeneous Dirichlet nodes
    if (eqnNum < 0)
      continue;
    /* building Kglob_e */
    // extract row of global stiffness matrix
    for (unsigned int col_idx = 0; col_idx < ndofe; col_idx++) {
      if (eqnVec[col_idx]-1 >= 0) {
        Kinv_e[eqn_idx][col_idx] = invK[eqnNum][eqnVec[col_idx]-1];
        LOG_DBG3(sgpopt) << "Kglob_e set entry (" << eqn_idx << "," << col_idx << ") with val=" << invK[eqnNum][eqnVec[col_idx]-1];
      }
    }
  }
  LOG_DBG3(sgpopt) << "CG: Kglob_e = " << Kinv_e.ToString();

  LOG_DBG3(sgpopt) << "CG: Mglob=" << Be.ToString();

  assert(gamma.GetNumCols() > 0 && gamma.GetNumRows() > 0);
  // Gamma_e = B_e^\top K^(-1) B_e = B_e^\top * M
  Matrix<Double> tmp(Be.GetNumCols(), Be.GetNumRows());
  // tmp = B_e^\top K^(-1)
  Be.MultT(Kinv_e, tmp);
  tmp.Mult(Be, gamma);

  LOG_DBG3(sgpopt) << "CG: gamma=" << gamma.ToString();
}

void SGP::GetIntPoints(const BaseBDBInt* bdb, const Elem* elem, BaseFE** ptFe, StdVector<LocPoint>& intPoints, StdVector<double>& weights) const{
  // annoying entity iterator that holds the elem
  ElemList el(grid_);
  el.SetElement(elem);

  IntegOrder order;
  IntScheme::IntegMethod method = IntScheme::UNDEFINED;
  // also sets method and order
  *ptFe = bdb->GetFeSpace1()->GetFe(el.GetIterator(), method, order );
  bdb->GetIntScheme()->GetIntPoints(Elem::GetShapeType(elem->type), method, order, intPoints, weights);
}

void SGP::GradientCheck() {
  unsigned int dim = grid_->GetDim();
  assert(dim == 2 || dim == 3);
  unsigned int n_var = dim == 2 ? 6 : 9;
  // largest difference between fd approximation and analytical gradient
  double max_err = 1e-12;
  // element where largest difference occurs
  int max_err_elemNum = -1;
  int max_err_idx = -2;

  // calculate finite difference gradient from objective function with central difference quotient
  StdVector<unsigned int> design_idx(n_var);
  design_idx[0] = space_->FindDesign(DesignElement::MECH_11);
  design_idx[1] = space_->FindDesign(DesignElement::MECH_12);
  design_idx[2] = space_->FindDesign(DesignElement::MECH_13);
  design_idx[3] = space_->FindDesign(DesignElement::MECH_22);
  design_idx[4] = space_->FindDesign(DesignElement::MECH_23);
  design_idx[5] = space_->FindDesign(DesignElement::MECH_33);

  if (dim == 3) {
    // assume orthotropic tensor
    design_idx[6] = space_->FindDesign(DesignElement::MECH_44);
    design_idx[7] = space_->FindDesign(DesignElement::MECH_55);
    design_idx[8] = space_->FindDesign(DesignElement::MECH_66);
  }

  // step size for finite difference gradient
  double h = 1e-6;
  StdVector<double> diff_grad(design_idx.GetSize() * n_elems_);
  // symmetric tensor: 6 non-zero entries in 2d and 9 in 3d
  StdVector<Vector<double>> fd(n_elems_);
  for (unsigned int e = 0; e < n_elems_; e++)
    fd[e].Resize(n_var);

  // get initial design
  Vector<double> x_init;
  optimization_->GetDesign()->WriteDesignToExtern(x_init);
//  std::cout << "x_init=" << x_init.ToString() << std::endl;
  double f_x = base_->EvalObjective((int) x_init.GetSize(), x_init.GetPointer(), true);
  base_->CommitIteration();

  StdVector<Matrix<double>> grad_ref;
  StdVector<Vector<double>> eps;
  CalcTensDeriv(grad_ref, eps);

  // write result to file
  std::string out_name = progOpts->GetSimName() + ".gradcheck.sgp";
  std::ofstream out_file(out_name.c_str());
  out_file << std::scientific << std::setprecision(7);
  if (!out_file.is_open())
    EXCEPTION("Could not open file");

  out_file << "### Gradientcheck for tensor derivative with step size " << h << "\n";
  if (dim == 2)
    out_file << "### order of gradient vector: [mech_11, mech_12, mech_13, mech_22, mech_23, mech_33]\n";
  else
    out_file << "### order of gradient vector: [mech_11, mech_12, mech_13, mech_22, mech_23, mech_33, mech_44, mech_55, mech_66]\n";

  // simple forward finite difference check: [f(x + h) - f(x)] / h
  for (unsigned int e = 0; e < n_elems_; e++) {
    // loop over all tensor entries and vary one tensor entry per loop iteration
    for (unsigned int j = 0; j < design_idx.GetSize();j++) {
      Vector<double> x(x_init);
      // e.g. get MECH_11 entry of element e
      unsigned int idx = design_idx[j] * n_elems_ + e;
//      std::cout << "j=" << j << " idx= " << idx << " x[j]=" << x[idx] << std::endl;
      x[idx] += h;
//      std::cout << "x_var=" << x.ToString(2) << std::endl;
      double f_xh = base_->EvalObjective((int) x.GetSize(), x.GetPointer(), true);
      double den = h;
      // for off-diagonal entries:
      // divide by 2*h as cfs automatically mirrors the change +h to ensure a symmetric tensor
      if (j == 1 || j == 2 || j == 4)
        den *= 2;
      fd[e][j] = (f_xh - f_x) / den;
    }

    out_file << "\nElement " << e+1 << ":\n";
    // assume symmetry
    // order 2D: mech_11, mech_12, mech_13, mech_22, mech_23, mech_33
    // order 3D: mech_11, mech_12, mech_13, mech_22, mech_23, mech_33, mech_44, mech_55, mech_66
//    double grad_vec[6] = {grad_ref[e][0][0], grad_ref[e][0][1], grad_ref[e][0][2], grad_ref[e][1][1], grad_ref[e][1][2], grad_ref[e][2][2]};
    double grad_vec[n_var];
    grad_vec[0] = grad_ref[e][0][0];
    grad_vec[1] = grad_ref[e][0][1];
    grad_vec[2] = grad_ref[e][0][2];
    grad_vec[3] = grad_ref[e][1][1];
    grad_vec[4] = grad_ref[e][1][2];
    grad_vec[5] = grad_ref[e][2][2];
    if (dim == 3) {
      // grad_44
      grad_vec[6] = grad_ref[e][3][3];
      // grad_55
      grad_vec[7] = grad_ref[e][4][4];
      // grad_66
      grad_vec[8] = grad_ref[e][5][5];
    }

    // put it into a vector object to get a nice output format
    Vector<double> vec(n_var, grad_vec, false);
    Vector<double> diff;
    diff = vec - fd[e];
    int loc;
    double tmp = std::abs(diff.MaxAbs(loc));
    if (tmp > max_err) {
      max_err = tmp;
      max_err_elemNum = e;
      max_err_idx = loc;
    }
    out_file << "analytical: " << vec.ToString(TS_PYTHON,"",7) << "\n";
    out_file << "fd approxi: " << fd[e].ToString(TS_PYTHON,"",7) << "\n";
  }

  out_file << "\n### largest absolute error: " << max_err << " at element " << max_err_elemNum << " idx=" << max_err_idx << "\n";

  out_file.close();
}

void SGP::SetStringValue(const std::string& key, const std::string& value)
{
  if (key == "approximation")
  { // switch does not support strings
    if (value == "asymptotes") {
      approximation_ = ZERO_ASYMPTOTES;
      L_ = 0;
      LOG_DBG3(sgpopt) << "set 'approximation' to type 'asymptotes'";
    }
    else if (value == "exact") {
      approximation_ = EXACT;
      LOG_DBG3(sgpopt) << "set 'approximation' to type 'exact'";
    }
    else
      EXCEPTION(value + " is not valid for " + key + ". Know only 'asymptotes' or 'exact'.");
    return;
  }

  if (key == "filtering")
  {
    if (value == "linear")
      filtering_ = LINEAR;
    else if (value == "non_linear")
      filtering_ = NON_LINEAR;
    else
      EXCEPTION(value + " is not valid for " + key + ". Know only 'linear' or 'non_linear'.");

    LOG_DBG3(sgpopt) << "set 'filtering' to type " << value;
    return;
  }

  EXCEPTION(key + " does not exist for type 'string'!");
}

void SGP::SetIntegerValue(const std::string& key, int value)
{
  if (key == "levels") {
    levels_ = value;
    assert(levels_ > 0);
    LOG_DBG3(sgpopt) << "SIV: levels of hierarichal grid = " << levels_;
    return;
  }

  if (key == "samples_per_level") {
    // set same value for density and angles
    samples_per_level_density_ = value;
    samples_per_level_angle_ = value;
    LOG_DBG3(sgpopt) << "SIV: samples per level density + angles = " << value;
    return;
  } else {
    // allow different sample numbers for density and angles
    if (key == "samples_per_level_density") {
      samples_per_level_density_ = value;
      LOG_DBG3(sgpopt) << "SIV: samples per level density = " << samples_per_level_density_;
      assert(levels_ > 0);
      return;
    }

    if (key == "samples_per_level_angle") {
      samples_per_level_angle_ = value;
      LOG_DBG3(sgpopt) << "SIV: samples per level angle = " << samples_per_level_angle_;
      assert(levels_ > 0);
      return;
    }
  }

  if (key == "max_bisections") {
    max_bisections_ = value;
    LOG_DBG3(sgpopt) << "SIV: max. number of iterations for vol. bisection = " << max_bisections_;
    return;
  }

  if (key == "max_globs") {
      max_globs_ = value;
      LOG_DBG3(sgpopt) << "SIV: max. number of globalization increments = " << max_globs_;
      return;
    }

  EXCEPTION(key + " unknown for type integer! Choices: levels, samples_per_level, max_bisections");
}

void SGP::SetNumericValue(const std::string& key, double value)
{
  if (key == "tau_init") {
    tau_init_ = value;
    LOG_DBG3(sgpopt) << "SNV: initial globalization tau = " << tau_init_;
    return;
  }
  if (key == "tau_factor"){
    tau_factor_ = value;
    LOG_DBG3(sgpopt) << "SNV: factor for tau update = " << tau_factor_;
    return;
  }
  if (key == "p_filt_density"){
    p_filt_dens_ = value;
    LOG_DBG3(sgpopt) << "SNV: penalty for density filtering term = " << p_filt_dens_;
    return;
  }
  if (key == "p_filt_angle"){
    p_filt_angle_ = value;
    LOG_DBG3(sgpopt) << "SNV: penalty for orientation filtering term = " << p_filt_angle_;
    return;
  }
  if (key == "pmin_vol"){
    pmin_vol_ = value;
    LOG_DBG3(sgpopt) << "SNV: lower bound for volume bisection = " << pmin_vol_;
    return;
  }
  if (key == "pmax_vol"){
    pmax_vol_ = value;
    LOG_DBG3(sgpopt) << "SNV: upper bound for volume bisection = " << pmax_vol_;
    return;
  }
  if (key == "tolerance"){
      tolerance_ = value;
      LOG_DBG3(sgpopt) << "SNV: tolerance for stopping criterion = " << std::scientific << std::setprecision(9) << tolerance_;
      return;
  }
  if (key == "volume_tolerance"){
    volume_tolerance_ = value;
    LOG_DBG3(sgpopt) << "SNV: stopping tolerance for volume bisection = " << std::scientific << std::setprecision(9) <<  volume_tolerance_;
    return;
  }

  EXCEPTION(key + " unknown for type real! Choices: tau_init, tau_factor, penalty_filter, pmin_vol, pmax_vol");
}

void SGP::ToInfo(PtrParamNode summary, sgp::Info info)
{
  if (info.status == sgp::Status::CONVERGED)
    summary->Get("break/converged")->SetValue("yes");
  else
    summary->Get("break/converged")->SetValue("no");

  switch (info.status) {
    case sgp::Status::CONVERGED:
      summary->Get("break/converged")->SetValue("yes");
      return;
    case sgp::Status::MAX_ITERATIONS:
      summary->Get("reason/msg")->SetValue("maximum iterations reached");
      return;
    case sgp::Status::RUNNING:
      summary->Get("reason/msg")->SetValue("running");
      return;
    default:
      summary->Get("break/converged")->SetValue("no");
  }

  assert(info.status != sgp::Status::NOT_AVAILABLE);
}

void SGP::SetBounds() //TODO: avoid repetitions of code to make it less prone to copy-paste errors
{
  unsigned int n_constraints = optimization_->constraints.view->GetNumberOfActiveConstraints();
  StdVector<double> gl, gu;
  gl.Resize(n_constraints);
  gu.Resize(n_constraints);

  if (n_constraints > 0)
  {
    LOG_TRACE(sgpopt) << "SB: volume bound = " << gu[0];
  }

  // we get bounds for each design element type
  StdVector<double> xl, xu;
  double lower = 1e-3;
  double upper = 0;
  double vol_bound = 0;

  if (hasDensity_) {
    space_->WriteBoundsToExtern(xl, xu, DesignElement::DENSITY);
    lower = xl[0]; // minimum density value
    upper = xu[0]; // maximum density value

    // set volume bound
    base_->GetConstraintsBounds(n_constraints, gl.GetPointer(), gu.GetPointer());
    assert(gl.GetSize() == gu.GetSize());
    // so far: only one real constraint for volume
    assert(gu[0] >= -1e-6);
    vol_bound = gu[0];
  } else {
    // assume density is fixed to a constant value and set volume bound and upper density bound to this value
    upper = Optimization::context->dm->GetParameter(DesignElement::DENSITY);
    lower = upper;
    vol_bound = upper;
  }

  solve_params_.v_star = vol_bound;
  solve_params_.rhoMin = lower;
  solve_params_.rhoMax = upper;

  LOG_TRACE(sgpopt) << "SB: bounds on density: lower=" << solve_params_.rhoMin << "  upper=" << solve_params_.rhoMax;

  lower = 0;
  upper = 0;
  if (hasRotAngleFirst_) {
    xl.Clear();
    xu.Clear();
    if (grid_->GetDim() == 2)
      space_->WriteBoundsToExtern(xl, xu, DesignElement::ROTANGLE);
    else
      space_->WriteBoundsToExtern(xl, xu, DesignElement::ROTANGLEFIRST);

    lower = xl[0];
    upper = xu[0];
  }
  solve_params_.phiMin   = lower;  // minimum first angle value
  solve_params_.phiMax   = upper;  // maximum first angle value
  LOG_TRACE(sgpopt) << "SB: bounds on first angle: lower=" << solve_params_.phiMin << "  upper=" << solve_params_.phiMax;

  lower = 0;
  upper = 0;
  if (hasRotAngleSecond_) {
    xl.Clear();
    xu.Clear();
    space_->WriteBoundsToExtern(xl, xu, DesignElement::ROTANGLESECOND);

    lower = xl[0];
    upper = xu[0];
  }
  solve_params_.thetaMin = lower;  // minimum second value
  solve_params_.thetaMax = upper;  // maximum second value
  LOG_TRACE(sgpopt) << "SB: bounds on second angle: lower=" << solve_params_.thetaMin << "  upper=" << solve_params_.thetaMax;

  lower = 0;
  upper = 0;
  if (hasRotAngleThird_) {
    xl.Clear();
    xu.Clear();
    space_->WriteBoundsToExtern(xl, xu, DesignElement::ROTANGLETHIRD);

    lower = xl[0];
    upper = xu[0];
  }
  solve_params_.psiMin   = lower;  // minimum third value
  solve_params_.psiMax   = upper;  // maximum third value
  LOG_TRACE(sgpopt) << "SB: bounds on third angle: lower=" << solve_params_.psiMin << "  upper=" << solve_params_.psiMax;
}

void SGP::SetInitialDesign()
{
  // get initial design
  Vector<double> x_init(n_elems_);

  sgp::Vector rho0(n_elems_);
  sgp::Vector phi0(n_elems_);
  sgp::Vector theta0(n_elems_);
  sgp::Vector psi0(n_elems_);

  if (hasDensity_) {
    optimization_->GetDesign()->WriteDesignToExtern(x_init, false, DesignElement::DENSITY);
  }  else {
    double rho_const = Optimization::context->dm->GetParameter(DesignElement::DENSITY);
    x_init.Init(rho_const);
  }
  rho0.assign(x_init.GetPointer(), x_init.GetPointer() + n_elems_);
  std::ostringstream os;
  std::for_each(rho0.begin(), rho0.end(), [&os](double i)->void {os << i << " ";});
  LOG_DBG2(sgpopt) << "SID: initial densities: " << os.str();
  if (hasRotAngleFirst_) {
    x_init.Clear();
    // check the correct type
    DesignElement::Type type = grid_->GetDim() == 2 ? DesignElement::ROTANGLE : DesignElement::ROTANGLEFIRST;
    optimization_->GetDesign()->WriteDesignToExtern(x_init, false, type);
    phi0.assign(x_init.GetPointer(), x_init.GetPointer() + n_elems_);

    std::ostringstream os;
    std::for_each(phi0.begin(), phi0.end(), [&os](double i)->void {os << i << " ";});
    LOG_DBG2(sgpopt) << "SID: initial first angle: " << os.str();
  }

  if (hasRotAngleSecond_) {
    assert(grid_->GetDim() == 3);
    x_init.Clear();
    optimization_->GetDesign()->WriteDesignToExtern(x_init, false, DesignElement::ROTANGLESECOND);
    theta0.assign(x_init.GetPointer(), x_init.GetPointer() + n_elems_);

    std::ostringstream os;
    std::for_each(theta0.begin(), theta0.end(), [&os](double i)->void {os << i << " ";});
    LOG_DBG2(sgpopt) << "SID: initial second angle: " << os.str();
  }

  if (hasRotAngleThird_) {
    assert(grid_->GetDim() == 3);
    x_init.Clear();
    optimization_->GetDesign()->WriteDesignToExtern(x_init, false, DesignElement::ROTANGLETHIRD);
    psi0.assign(x_init.GetPointer(), x_init.GetPointer() + n_elems_);

    std::ostringstream os;
    std::for_each(psi0.begin(), psi0.end(), [&os](double i)->void {os << i << " ";});
    LOG_DBG2(sgpopt) << "SID: initial third angle: " << os.str();
  }

  solve_params_.rho0   = rho0; // initial material configuration
  solve_params_.phi0   = phi0; // initial material configuration
  solve_params_.theta0 = theta0; // initial material configuration
  solve_params_.psi0   = psi0; // initial material configuration
}

void SGP::SetEnums() {
  approximation.SetName("SGP::Approximation");
  approximation.Add(EXACT, "exact");
  approximation.Add(ZERO_ASYMPTOTES, "zeroAsymptotes");

  filtering.SetName("SGP::Filtering");
  filtering.Add(LINEAR, "linear");
  filtering.Add(NON_LINEAR, "non-linear");
}

