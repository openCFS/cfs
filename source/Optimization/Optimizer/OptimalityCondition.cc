#include <assert.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <boost/math/special_functions/fpclassify.hpp>

#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Domain/Domain.hh"
#include "Domain/ElemMapping/Elem.hh"
#include "Driver/Assemble.hh"
#include "General/defs.hh"
#include "General/Environment.hh"
#include "General/Exception.hh"
#include "Optimization/Condition.hh"
#include "Optimization/Excitation.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Function.hh"
#include "Optimization/Objective.hh"
#include "Optimization/Optimization.hh"
#include "Optimization/Optimizer/OptimalityCondition.hh"
#include "PDE/BasePDE.hh"
#include "Utils/StdVector.hh"
#include "Utils/tools.hh"

using namespace CoupledField;
using std::abs;

DEFINE_LOG(ocm, "ocm")

OptimalityCondition::OptimalityCondition(Optimization* optimization, PtrParamNode pn)
 : BaseOptimizer(optimization, pn, Optimization::OCM_SOLVER)
{
  type.SetName("OptimalityCondition::Type");
  type.Add(FRAMED, "framed");
  type.Add(FUMBLE, "fumble");
  type.Add(TRAJECTORY, "trajectory");
  type.Add(EXTREMIZE, "extremize");
  
  pn = pn->Get(Optimization::optimizer.ToString(Optimization::OCM_SOLVER), ParamNode::PASS);

  // read the xml values
  if(pn != NULL)
  {
    type_             = type.Parse(pn->Get("type")->As<std::string>());
    if(pn->Get("move_limit")->As<string>() == "tuned")
    {
      tuned = new Tuned(pn->Get("tuned", ParamNode::PASS), &move_limit_, 5, 5, this);
      move_limit_ = tuned->max; // set by Tuned
    }
    else
      move_limit_ = pn->Get("move_limit")->As<double>();

    initial_move_limit_ = move_limit_;
    oc_damping_       = pn->Get("damping")->As<double>();
    lambda_min_       = pn->Get("lambda_min")->As<double>();
    feasibility_      = pn->Get("feasibility")->As<double>();
    max_lambda_iters_ = pn->Get("max_lambda_iters")->As<int>();
    constr_grad_      = pn->Get("constr_grad")->As<bool>();

    // it doesn't harm to read the parameters for all types!
    if(pn->Has(type.ToString(FRAMED)))
    {
      PtrParamNode t = pn->Get(type.ToString(FRAMED));
      // there are defaults in XML
      always_enlarge_    = t->Get("alwaysEnlarge")->As<bool>();
      start_lower_       = t->Get("lower")->As<double>();
      start_upper_       = t->Get("upper")->As<double>();
      enlarge_lower_     = t->Get("enlargeLower")->As<double>();
      enlarge_upper_     = t->Get("enlargeUpper")->As<double>();
      check_stalled_err_ = t->Get("checkErr")->As<bool>();
      overcome_deadlock_ = t->Get("overcome_deadlock")->As<bool>();

      if(start_lower_ < 0)
        throw Exception("no negative lower bound frame allowed in current implementation");
      if(enlarge_lower_ > 1.0 || enlarge_lower_ < 0.0 ) 
        throw Exception("the 'frame' value 'enlargeLower' shall be in [0;1]");
      if(enlarge_upper_ < 1.0) 
        throw Exception("the 'frame' value 'enlargeUpper' shall be in >= 1.0");
    }
    if(pn->Has(type.ToString(FUMBLE)))    
    {
      PtrParamNode t = pn->Get(type.ToString(FUMBLE));
      
      step_     = t->Get("step")->As<double>();
      contract_ = t->Get("contract")->As<double>();
      expand_   = t->Get("expand")->As<double>();
      
      if(expand_ <= 1.0)
        throw Exception("expand shall be > 1.0, e.g. 1.99");
      if(contract_ >= 1.0 )
        throw Exception("contract shall be < 1.0, e.g. 0.49");
    }

    if(pn->Has("multiObjective"))
    {
      PtrParamNode multiObj = pn->Get("multiObjective");
      this->weight_ = multiObj->Get("multiObjectiveWeight")->As<double>();
      this->beta_ = multiObj->Get("beta")->As<double>();

      // Get index of "real" constraint
      Function::Type constraint_type = Function::type.Parse(multiObj->Get("constraint")->As<std::string>());
      constraint_idx_ = -1;
      for(int i=0; i < optimization->constraints.view->GetNumberOfActiveConstraints(); ++i)
        if(optimization->constraints.view->Get(i)->GetType() == constraint_type)
            constraint_idx_ = i;
      optimization->constraints.view->Done();

      if(constraint_idx_ == -1)
        EXCEPTION("multiObjective constraint type " << Function::type.ToString(constraint_type) << " is not given as constraint.")

      LOG_DBG(ocm) << "mO: w=" << weight_ << " b=" << beta_ << " c=" << Function::type.ToString(constraint_type) << " constraint_idx_=" << constraint_idx_;
    }
  }

  // some plausibility about optimality condition
  if(type_ == EXTREMIZE && optimization->constraints.view->GetNumberOfActiveConstraints() > 0)
    throw Exception("Optimality Condition 'extremize' is not implemented for constraints");

  LOG_DBG(ocm) << "OptimalityCondition of type " << type.ToString(type_) << " active constraints=" << optimization->constraints.view->GetNumberOfActiveConstraints();

  vault_.Resize(optimization->GetDesign()->data.GetSize());
  evaluate_tmp_.Resize(optimization->GetDesign()->data.GetSize());
  
  optimizer_timer_->Stop();

  PostInitScale(1.0, true);
}

void OptimalityCondition::DescribeProperties(StdVector<std::pair<std::string, std::string> >& map) const
{
  map.Push_back(std::make_pair("move_limit", std::to_string(move_limit_)));
  map.Push_back(std::make_pair("damping", std::to_string(oc_damping_)));
  map.Push_back(std::make_pair("biscect", type.ToString(type_)));
}

void OptimalityCondition::SolveProblem()
{
  // we measure the optimizer in the loops only
  optimizer_timer_->Stop();

  // solve the state problem first
  Optimization::context->pde->GetAssemble()->SetAllReassemble(); // tell assemble that the Design has changed
  optimization->SolveStateProblem();
  optimization->CalcObjective();

  // start with iteration 0 which is the initial design
  bool stalled_err = false; // framed only
  assert(optimization->GetCurrentIteration() == 0);

  LOG_DBG(ocm) << "SP: #active constraints=" << optimization->constraints.view->GetNumberOfActiveConstraints();

  // run at least once such that we have observes and constraints for the stopping criteria evaluated
  do
  {
    // calc gradients to store the results in data[element]...
    // only the gradients are needed for the calculation of the next iteration

    optimization->SolveAdjointProblems();

    eval_grad_obj_timer_->Start();
    optimization->CalcObjectiveGradient(nullptr);
    eval_grad_obj_timer_->Stop();
    
    // reset values of the constraint gradients
    optimization->GetDesign()->Reset(DesignElement::CONSTRAINT_GRADIENT, DesignElement::DEFAULT);

    ConditionContainer& cc = optimization->constraints;

    // if we have more than one constraints
    // we agglomerate all but the chosen "real" constraint in a smooth minimum and combine it with the objective
    // (1-weight_) * objective - weight_ * smoothmin(constraints)

    // here we get all "fake" constraint values -> needed for derivative of smooth minimum
    StdVector<double> constraint_values;
    if(cc.view->GetNumberOfActiveConstraints() > 1)
    {
      for(int i = 0; i < cc.view->GetNumberOfActiveConstraints(); i++)
      {
        if(i != constraint_idx_)
        {
          Condition* c = cc.view->Get(i);
          optimization->CalcConstraint(c);
          constraint_values.Push_back(c->IsLocal() ? ((LocalCondition*)c)->GetValue(): c->GetValue());
        }
      }
      cc.view->Done();

      LOG_DBG2(ocm) << " constraint_values=" << constraint_values.ToString(TS_PYTHON);
      LOG_DBG2(ocm) << " min(con_vals)=" << SmoothMin(constraint_values, beta_, true);
    }

    // calculate the derivative of the "real" constraint and the smooth minimum of the "fake" constraints
    if(cc.view->GetNumberOfActiveConstraints() > 0) {
      StdVector<double> grad(optimization->GetDesign()->data.GetSize());
      grad.window.Set(0, grad.GetSize());
      agglomerated_grad_.Resize(optimization->GetDesign()->data.GetSize());
      agglomerated_grad_.Init(0);

      //SolveAdjointProblemsIfNeeded(n, x, cfs_scale);

      eval_grad_const_timer_->Start();

      int fake_constraint_idx = 0;
      for(int i = 0; i < cc.view->GetNumberOfActiveConstraints(); ++i)
      {
        if(i == constraint_idx_)
          // "real" constraint
          optimization->CalcConstraintGradient(cc.view->Get(i));
        else
        {
          // "fake" constraints
          optimization->CalcConstraintGradient(cc.view->Get(i), &grad);

          LOG_DBG2(ocm) << "e=" << i << " x=" << optimization->GetDesign()->data[i].GetValue(DesignElement::DESIGN, DesignElement::SMART) << " f=" << optimization->CalcConstraint(cc.view->Get(i)) << " df=" << grad.ToString(2);
          // agglomerate all gradients of "fake" constraints, chain rule applies
          // d/drho smoothmin(f_1(rho), ..., f_n(rho)) = sum_i d/df_i smoothmin(f_1(rho), ..., f_n(rho)) * d/drho f_i(rho)
          for(unsigned int j = 0; j < grad.GetSize(); ++j)
          {
            agglomerated_grad_[j] += DerivSmoothMin(constraint_values, beta_, fake_constraint_idx) * grad[j];
            LOG_DBG3(ocm) << " dmin * df_SMART = " << DerivSmoothMin(constraint_values, beta_, fake_constraint_idx) << " * " << grad[j]
                         << " = " << DerivSmoothMin(constraint_values, beta_, fake_constraint_idx) * grad[j];
          }
          ++fake_constraint_idx;
        }
      }
      cc.view->Done();
      eval_grad_const_timer_->Stop();

      LOG_DBG3(ocm) << " agglo_grad=" << agglomerated_grad_.ToString(TS_PYTHON);

      DesignSpace* design = optimization->GetDesign();
      int res_idx = design->GetSpecialResultIndex(DesignElement::GENERIC_ELEM, "ocm_agglomerated_grad_");
      if(res_idx != -1)
      {
        for(unsigned int j = 0; j < grad.GetSize(); ++j)
          design->data[j].specialResult[res_idx] = agglomerated_grad_[j];
      }
    }

    // store iteration 0
    if(optimization->GetCurrentIteration() == 0)
    {
      eval_obj_timer_->Start();
      design_.value = optimization->CalcObjective();   // for output
      eval_obj_timer_->Stop();

      CommitIteration(); // don't assert we are running
      continue; // redo gradients and start optimization
    }

    optimizer_timer_->Start();
    nan_fixes.Push_back(0);
    LOG_DBG(ocm) << "SO: CalcNext(" << type.ToString(type_) << ")Iteration() iter=" << optimization->GetCurrentIteration() << " ml=" << move_limit_;
    // do a SIMP Optimality Condition step -> calc new design vector
    switch(type_)
    {
    case FRAMED:     stalled_err = CalcNextFramedIteration(stalled_err);
                     break;

    case FUMBLE:     CalcNextFumbleIteration();
                     break;

    case TRAJECTORY: CalcNextTrajectoryIteration();
                     break;

    case EXTREMIZE:  CalcNextExtremizeIteration();
                     break;

    default: assert(false);
    }
    optimizer_timer_->Stop();

    // solve the state problem for the new design vector
    // we have to set reassemblence for all pdes
    for(unsigned int i = 0; i < Optimization::manager.context.GetSize(); i++) {
      Optimization::manager.SwitchContext(i);
      Optimization::context->pde->GetAssemble()->SetAllReassemble();
    }
    optimization->SolveStateProblem();

    // calc the objective and "fake" constraint values for the logging in CommitIteration(),
    // for the optimality condition it is not required.

    eval_obj_timer_->Start();
    design_.value = optimization->CalcObjective();
    eval_obj_timer_->Stop();

    // every state problem is an iteration
    // The gradients "point" to this design vector.
    CommitIteration(); // don't assert we are running
  }
  while(!optimization->DoStopOptimization() && optimization->GetCurrentIteration() <= optimization->GetMaxIterations());
  
  PtrParamNode in = optimization->optInfoNode->Get(ParamNode::SUMMARY)->Get("break");
  
  if(optimization->GetCurrentIteration() >= optimization->GetMaxIterations()-1)
  {
    optimization->DoStopOptimizationHelper(false, "Maximum iterations exceeded");
    std::cout << "++ max iterations reached" << std::endl;
  }
  assert(in->GetChildren().GetSize() > 0);
} 


bool OptimalityCondition::CalcNextFramedIteration(bool last_was_stalled_err)
{
  // framed is the standard

  // we store the current densities in the temp variable. Otherwise we cannot
  // find the proper lambda
  optimization->GetDesign()->WriteDesignToExtern(vault_.GetPointer());

  // set the frame borders
  if(lower_ == upper_)
  {
    // this is only reached in the first iteration
    lower_ = start_lower_;
    upper_ = start_upper_;
  }
  else if(last_was_stalled_err)
  {
    // when the move limit was wrong, lambda can be arbitrary nonsense
    // But we also may not start from start_*_ as then we repeat the situation
    // this heuristic is not really cool - better start feasible or adjust the move_limit
    lower_ = lambda_ * 1e-4;
    upper_ = std::min(2*start_upper_, lambda_ * 1e4);
  }
  else
  {
    //FIXME
    assert(always_enlarge_);
    lower_ = lambda_ * (always_enlarge_ ? enlarge_lower_ : 1.0);
    upper_ = lambda_ * (always_enlarge_ ? enlarge_upper_ : 1.0);
  }

  // we count lambda iterations to handle the problem of coming too close to a boundary
  lambda_iters_ = 0;
  // stalled error means that we cannot find a better lamda. Reason can be move limit or lambda bounds
  bool stalled_err = false; // when | err - last_err | is zero plus other conditions

  // we enable move limit enlarging via overcome_deadlock_
  bool no_further_lamda_search = false; // when we cannot enlarge ml any more because of org bound or not enabled

  // this is the move_limit before we eventually enlarge it. To detect enlarging and reset before next iteration
  double iter_org_move_limit = move_limit_;

  double err = -1; // volume bound error

  LOG_DBG2(ocm) << "CNFI: w_s=" << last_was_stalled_err << " l=" << lower_ << " u=" << upper_ << " ml=" << move_limit_ << " se=" << stalled_err << " fix=" << overcome_deadlock_;

  while(!no_further_lamda_search)
  {
    int ml_iter = 0; // own iteration counter within move limits to check stalled_err not too early
    double last_err = -2;   // when move limit is a bound: stop, when err is not changing
    err = -1;
    do
    {
      // calc next lambda
      lambda_ = 0.5 * (upper_ + lower_);

      // evaluate with new lambda
      // restore original density from temp so we always start the calculation
      // on the same base but with different lambda
      // do not write density file
      optimization->GetDesign()->ReadDesignFromExtern(vault_.GetPointer(), false);

      last_err = err;
      err = Evaluate(lambda_);

      // we want to detect, if we have too small move limits. We do this, when for several lambda, the error stays the same
      // for initial iterations it may take time to fine a proper lambda range first
      // offshot lambda can have the same bounded error - especially with small move limit
      int lambda_err_tol = std::min(2 * 1/move_limit_, 20.);

      // when delta err does not change, but err exits
      // and also not to early, as initially err might stay constant in a valid setting (-> Technical/mech_2d has twice err=.2)
      stalled_err = check_stalled_err_ && ml_iter > lambda_err_tol && abs(err) > feasibility_ && abs(err - last_err) < .5 * feasibility_;
      // move frames according to new lambda
      if(err > 0) upper_ = lambda_; // = center
      else lower_ = lambda_;

      ml_iter++;
      lambda_iters_++;

      LOG_DBG2(ocm) << "CNFI: ml=" << move_limit_ << " mli=" << ml_iter << " li=" << lambda_iters_ << " l=" << lambda_ << " err=" << err << " lo=" << lower_ << " up=" << upper_
          << " delta_err=" << abs(err - last_err) << " lt=" << lambda_err_tol << " se=" << stalled_err;
    }
    while(abs(err) > feasibility_  && lambda_iters_ < max_lambda_iters_ && !stalled_err && abs(lambda_) > abs(lambda_min_));

    double max_ml = std::max(initial_move_limit_, 0.5); // could be improved with real bounds
    if(stalled_err && overcome_deadlock_ && 2*move_limit_ <= max_ml)
    {
      std::cout << "OCM cannot further reduce error: lambda=" << lambda_ << " err=" << err << " move_limit=" << move_limit_ << std::endl;
      move_limit_ *= 2;
      lower_ = start_lower_;
      upper_ = start_upper_;
    }
    else
    {
      no_further_lamda_search = true;
    }
  }
  if(move_limit_ != iter_org_move_limit)
  {
    std::stringstream ss;
    ss << "OCM move_limit temporarily increased " << iter_org_move_limit << " to " << move_limit_
       << " in iteration " << optimization->GetCurrentIteration() << " lambda=" << lambda_;
    info_->Get(ParamNode::PROCESS)->SetWarning(ss.str());
    move_limit_ = iter_org_move_limit; // reset such that the next iteration is plain
  }
  assert(IsWithinBounds(evaluate_tmp_));

  // do not reset move_limit_ here but let it be logged by CommitIteration
  // write design to density file
  optimization->GetDesign()->ReadDesignFromExtern(evaluate_tmp_.GetPointer(), true);

  if(abs(lambda_) < abs(lambda_min_))
   std::cout << "Iteration fails due to too small lambda, check feasibility, move_limit or other ocm@type than 'framed'" << std::endl;

  if(stalled_err) {
   std::cout << "Iteration fails due to too stalled error: enlarge move_limit (" << move_limit_ << ") or disable ocm/framed@checkErr" << std::endl;
  }
  if(lambda_iters_ >= max_lambda_iters_)
   std::cout << "Iteration fails to find valid Lagrangian: " << lambda_ << " err: " << err << " check move_limit or bounds in 'ocm/framed'\n";

  return stalled_err;
}

void OptimalityCondition::CalcNextFumbleIteration()
{
  assert(step_ > 0);  
  
  // we store the current densities in the temp variable. Otherwise we cannot
  // find the proper lambda
  optimization->GetDesign()->WriteDesignToExtern(vault_.GetPointer());
  
  // we have a lambda_ and a step_ which is by definition > 0.
  // the we evaluate at lambda_k + 0.5 * step_, lambda_k + 2.0 * step_,  lambda_k - 0.5 * step_, lambda_k - 2.0 * step_,
  // the least error becomes the new lambda_k+1 and step_ becomes |lambda_k+1 - lambda_k|

  lambda_iters_ = 0;
  double min_err;
  
  do
  {
    // if the initial lambda is far away all evaluations are limited and give 
    // the same error. Therefore we assume the initial lambda to be too large
    // and start with the lambda - expand * step check.
    
    // restore original density from temp so we always start the calculation 
    // on the same base but with different lambda
    // start with check lambda_ - expand_ * step_
    optimization->GetDesign()->ReadDesignFromExtern(vault_.GetPointer());    
    min_err = Evaluate(lambda_ - expand_ * step_);
    double fumble = -1.0 * expand_;

    // check with lambda_ - contract_ * step_
    optimization->GetDesign()->ReadDesignFromExtern(vault_.GetPointer());    
    double t = Evaluate(lambda_ - contract_ * step_);
    if(t < min_err) {// abs(t) < abs(min_err) ??
      fumble = -1.0 * contract_;
      min_err = t;
    }

    // check with lambda_ + contract_ * step_
    optimization->GetDesign()->ReadDesignFromExtern(vault_.GetPointer());
    t = Evaluate(lambda_ + contract_ * step_);
    if(t < min_err) {    
      fumble = contract_;
      min_err = t;
    }

    // check lambda_ + expand_ * step_
    optimization->GetDesign()->ReadDesignFromExtern(vault_.GetPointer());    
    t = Evaluate(lambda_ + expand_ * step_);
    if(t < min_err) {
      fumble = expand_;
      min_err = t;
    }
    
    // new lambda:
    lambda_ += fumble * step_;
    step_ = abs(fumble * step_);
    
    lambda_iters_++;

    LOG_DBG2(ocm) << "lambda_iter/lambda/err/step/fumble/-ext/-cont/+con/+ext = " <<  lambda_iters_ << "\t"
                 << lambda_ << "\t" << min_err << "\t" << (fumble < 0 ? -1.0 * step_ : step_) << "\t" << fumble << "\t"
                 << Evaluate(lambda_ - expand_ * step_) << "\t" << Evaluate(lambda_ - contract_ * step_) << "\t"
                 << Evaluate(lambda_ + contract_ * step_) << "\t" <<  Evaluate(lambda_ + expand_ * step_);
  }
  while(abs(min_err) > feasibility_ && lambda_iters_ < max_lambda_iters_);

  if(lambda_iters_ >= max_lambda_iters_)
    std::cout << "Iteration fails to find valid Lagrangian: " << lambda_ << " err: " << min_err << std::endl;
}



void OptimalityCondition::CalcNextTrajectoryIteration()
{
  // we store the current densities in the temp variable. Otherwise we cannot
  // find the proper lambda
  optimization->GetDesign()->WriteDesignToExtern(vault_.GetPointer());
  
  double err = Evaluate(lambda_);

  LOG_DBG(ocm) << "CalcNextIteration: lambda= " << lambda_ << " -> err=" << err;
  
  int    last_dir = 0;
  double factor;
  lambda_iters_ = 0;

  while(abs(err) > feasibility_ && lambda_iters_ < max_lambda_iters_)
  {
    // we have to support to step over zero!
    if(abs(lambda_) < lambda_min_)
    {
      std::cout << "switch lamba from " << lambda_ << " to " << (lambda_ > 0 ? -1.0 : 1.0) << std::endl; 
      lambda_ = lambda_ > 0 ? -1.0 : 1.0;
      LOG_DBG(ocm) << "swap lambda sign to " << lambda_ << " at iter/err: " << lambda_iters_ << "\t" << "\t" << err;
    }
    else    
    {
      // normal operation, not too close to zero
      if(err > 0)
      {
        factor = last_dir > 0 ? 0.5 : 0.75;
        lambda_ *= factor;
        last_dir = 1;
      } 
      else
      {
        factor = last_dir < 0 ? 2.0 : 1.5;
        lambda_ *= factor;
        last_dir = -1;                   
      }
      LOG_DBG2(ocm) << "lambda/err/l_iter = " <<  lambda_<< "\t" << err << "\t" << lambda_iters_
                   << " factor = " << factor << " last_dir = " << last_dir;
    }
    // restore original density from temp so we always start the calculation 
    // on the same base but with different lambda
    optimization->GetDesign()->ReadDesignFromExtern(vault_.GetPointer());

    err = Evaluate(lambda_);
    
    lambda_iters_++;
  }
  
  if(lambda_iters_ >= max_lambda_iters_)
    std::cout << "Iteration fails to find valid Lagrangian: " << lambda_ << " err: " << err << std::endl;
}

void OptimalityCondition::CalcNextExtremizeIteration()
{
  StdVector<DesignElement>& data = optimization->GetDesign()->data;

  // we cannot set the design directly otherwise the filter does not 
  // work (it becomes unsymmetrically for symmetric problems) as all 
  // elements but the first have old and new elements in their filter
  // stencil. Hence we store in evaluate_tmp_
#pragma omp parallel for num_threads(CFS_NUM_THREADS)
  for(int i = 0; i < (int) data.GetSize(); i++)
  {
    DesignElement* de = &data[i];
    // rho_e is the old rho
    double rho_e = de->GetDesign(DesignElement::PLAIN);   
    double obj_grad = de->GetValue(DesignElement::COST_GRADIENT, DesignElement::SMART);
    
    // next is density times b_e which is compared with box constraints and move limit
    double next = rho_e + (obj_grad > 0 ? move_limit_ : -1.0 * move_limit_);        
                    
    double lower = std::max(de->GetLowerBound(), rho_e - move_limit_);
    double upper = std::min(de->GetUpperBound(), rho_e + move_limit_);            

    // we cannot set the design directly - otherwise the filter stencil get violated 
    evaluate_tmp_[i] = next;
    if(next <= lower)
      evaluate_tmp_[i] =lower;
    if(upper <= next)
      evaluate_tmp_[i] =upper;
    
    LOG_DBG3(ocm) << "CalcNextExtremizeIteration:" << de->elem->elemNum << " obj_grad=" << obj_grad
                 << "(" << de->GetValue(DesignElement::COST_GRADIENT, DesignElement::PLAIN) << ")" << " old= " << rho_e
                 << " next=" << next << " lower=" << lower << " upper=" << upper << " new=" << evaluate_tmp_[i];
  }
  
  // store the new values in the design variables
  optimization->GetDesign()->ReadDesignFromExtern(evaluate_tmp_.GetPointer());
}


double OptimalityCondition::Evaluate(double lambda)
{
  // we assume DensityElement.objective_gradient to be set
  // we assume DensityElement.constraint_gradient to be set or disabled
  LOG_DBG2(ocm) << "E: enter lambda=" << lambda << " ml=" << move_limit_;

  if(abs(lambda) < abs(lambda_min_))
  {
   double org_lambda = lambda;
   lambda = lambda < 0 ? -1.0 * lambda_min_ : lambda_min_;
   std::cout << "Optimality Condition evaluates with too small lambda "
             << org_lambda << " adjust to " << lambda << std::endl;
   LOG_DBG(ocm) << "Evaluate: adjust " << org_lambda << " to " << lambda;
  }

  ConditionContainer& cc = optimization->constraints;
  assert(cc.view->GetNumberOfActiveConstraints() > 0);
  Condition* g = cc.view->Get(constraint_idx_);
  StdVector<DesignElement>& data = optimization->GetDesign()->data;

  // we cannot set the design directly otherwise the filter does not
  // work (it becomes unsymmetrically for symmetric problems) as all
  // elements but the first have old and new elements in their filter
  // stencil. Hence we store in evaluate_tmp_
  int nan_cnt = 0;
#pragma omp parallel for num_threads(CFS_NUM_THREADS) reduction(+ : nan_cnt)
  for(int i = 0; i < (int) data.GetSize(); i++)
  {
    DesignElement* de = &data[i];
    // rho_e is the old rho
    double rho_e = de->GetDesign(DesignElement::PLAIN);

    // if filter is enabled we use the filtered value otherwise the plain one
    double smart_obj_grad = de->GetValue(DesignElement::COST_GRADIENT, DesignElement::SMART);

    // "real" objective = (1-weight_) * objective - weight_ * smoothmin(constraints)
    if(cc.view->GetNumberOfActiveConstraints() > 1 && constr_grad_)
      smart_obj_grad = (1-weight_) * smart_obj_grad - weight_ * agglomerated_grad_[de->GetIndex()];

    smart_obj_grad *= optimization->objectives.DoMaximize() ? -1.0 : 1.0;

    double b_e = -1.0 * smart_obj_grad;

    // for compliant mechanism and eigenvalue optimization the gradient can be positive, this is cut
    // -> Bendsoe/Sigmund. p 74, p 97
    // for piezo we might become negative lambdas -> cut the positive!
    b_e = lambda >= 0.0 ? std::max(0.0, b_e) : std::min(0.0, b_e);

    double smart_con_grad = constr_grad_ ? de->GetValue(DesignElement::CONSTRAINT_GRADIENT, DesignElement::SMART, g) : 1.0;
    b_e /= (lambda * smart_con_grad);

    // for Heaviside/tanh projection smart_con_grad can be close to zero
    double b_e_save = b_e;
    if(std::isnan(b_e_save))
    {
      b_e_save = 1e16; // something really big
      nan_cnt++;
    }
    // next is density times b_e which is compared with box constraints and move limit
    double next = rho_e * std::pow(b_e_save, oc_damping_);
    LOG_DBG3(ocm) << "Evaluate:" << de->elem->elemNum
                 << " smart_obj_grad=" << smart_obj_grad << " smart_con_grad=" << smart_con_grad
                 << "(" << Function::type.ToString(g->GetType()) << ")"
                 << " lambda=" << lambda << " b_e=" << b_e << " b_e_save=" << b_e_save << " damping=" << oc_damping_
                 << " rho_e=" << rho_e << " -> next=" << next;

    double lower = std::max(de->GetLowerBound(), rho_e - move_limit_);
    double upper = std::min(de->GetUpperBound(), rho_e + move_limit_);

    // we cannot set the design directly - otherwise the filter stencil get violated
    evaluate_tmp_[i] = std::max(lower, std::min(upper, next));

    LOG_DBG3(ocm) << "Evaluate:" << de->elem->elemNum<< " old=" << rho_e
                 << " next=" << next << " lower=" << lower << " upper=" << upper
                 << " -> new=" << evaluate_tmp_[i];
  }
  assert(nan_fixes.GetSize() > 0);
  nan_fixes[nan_fixes.GetSize()-1] = nan_cnt;
  optimizer_timer_->Stop();

  eval_const_timer_->Start();
  // store the new values in the design variables and evaluate the constraint
  optimization->GetDesign()->ReadDesignFromExtern(evaluate_tmp_.GetPointer(), false);
  double vol = optimization->CalcConstraint(g);
  eval_const_timer_->Stop();

  optimizer_timer_->Start();

  double err = g->GetBoundValue() - vol;

  LOG_DBG2(ocm) << "E: lambda=" << lambda << " v=" << vol << " e=" << err;

  return err;
}

bool OptimalityCondition::IsWithinBounds(const Vector<double>& data)
{
  const StdVector<DesignElement>& des = optimization->GetDesign()->data;

  assert(data.GetSize() == des.GetSize());
  for(unsigned int i = 0; i < data.GetSize(); i++)
    if(data[i] < des[i].GetLowerBound() || data[i] > des[i].GetUpperBound())
      return false;
  return true;
}

void OptimalityCondition::ToInfo(PtrParamNode info)
{
  PtrParamNode in = info_->Get(ParamNode::HEADER);
  in->Get("type")->SetValue(type.ToString(type_));
  in->Get("move_limit")->SetValue(tuned ? "tuned" : std::to_string(move_limit_));
    if(tuned)
      tuned->ToInfo(in->Get("tuned_move_limit"));
  in->Get("damping")->SetValue(oc_damping_);
  in->Get("max_lambda_iters")->SetValue(max_lambda_iters_);
  in->Get("lambda_min")->SetValue(lambda_min_);
  in->Get("constr_grad_")->SetValue(constr_grad_);

  if(type_ == FRAMED)
  {
    PtrParamNode t = in->Get(type.ToString(type_));
    t->Get("alwaysEnlarge")->SetValue(always_enlarge_);
    t->Get("lower")->SetValue(start_lower_);
    t->Get("upper")->SetValue(start_upper_);
    t->Get("enlargeLower")->SetValue(enlarge_lower_);
    t->Get("enlargeUpper")->SetValue(enlarge_upper_);
    t->Get("checkErr")->SetValue(check_stalled_err_);
  }

  if(nan_fixes.GetSize() > 0)
  {
    in = info_->Get(ParamNode::SUMMARY)->Get("nan_fixes");
    in->Get("total")->SetValue(nan_fixes.Sum());
    in->Get("last_iter")->SetValue(nan_fixes.Last());
  }
}

void OptimalityCondition::LogFileHeader(Optimization::Log& log)
{
  BaseOptimizer::LogFileHeader(log);

  log.AddToHeader("lambda");
  log.AddToHeader("lambda_iters");

  if(type_ == FUMBLE)
    log.AddToHeader("step");
}


void OptimalityCondition::LogFileLine(std::ofstream* out, PtrParamNode iteration)
{
  BaseOptimizer::LogFileLine(out, iteration);

  if(out)
    *out << " \t" << lambda_ << " \t" << lambda_iters_;

  iteration->Get("lambda")->SetValue(lambda_);
  iteration->Get("lambda_iters")->SetValue(lambda_iters_);

  if(type_ == FUMBLE)
  { 
    if(out)
      *out << " \t" << step_;
    iteration->Get("step")->SetValue(step_);
  }
}
