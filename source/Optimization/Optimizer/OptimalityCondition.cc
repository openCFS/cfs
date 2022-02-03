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

DEFINE_LOG(oc, "optimalityCondition")

OptimalityCondition::OptimalityCondition(Optimization* optimization, PtrParamNode pn)
 : BaseOptimizer(optimization, pn, Optimization::OPTIMALITY_CONDITION)
{
  type.SetName("OptimalityCondition::Type");
  type.Add(FRAMED, "framed");
  type.Add(FUMBLE, "fumble");
  type.Add(TRAJECTORY, "trajectory");
  type.Add(EXTREMIZE, "extremize");
  
  this->lambda_ = 0; // just to set a value

  // the following values are standard in mech SIMP -> see e.g. the 99 lines paper
  this->move_limit_ = 0.2;
  this->oc_damping_ = 0.5;
  this->lambda_min_ = 1e-20;
  this->lambda_iters_ = 0;
  this->max_lambda_iters_ = 70;
  this->feasibility_    = 1e-6;
  this->type_       = optimization->objectives.Has(Objective::COMPLIANCE) ? FRAMED : FUMBLE;

  // framed
  this->upper_ = 0.0;
  this->lower_ = 0.0;
  this->start_lower_ = 0;
  this->start_upper_ = 100000;
  this->enlarge_lower_ = 0.5;
  this->enlarge_upper_ = 2.0;
  this->always_enlarge_ = true;
  
  // fumble
  this->step_ = 10;
  this->contract_ = 0.49;
  this->expand_ = 1.99;

  // reduce to our actual ParamNode
  pn = pn->Get(Optimization::optimizer.ToString(Optimization::OPTIMALITY_CONDITION), ParamNode::PASS);
  
  // read the xml values
  if(pn != NULL)
  {
    type_       = type.Parse(pn->Get("type")->As<std::string>());
    move_limit_ = pn->Get("move_limit")->As<double>();
    oc_damping_ = pn->Get("damping")->As<double>();
    lambda_min_ = pn->Get("lambda_min")->As<double>();
    feasibility_    = pn->Get("feasibility")->As<double>();
    max_lambda_iters_ = pn->Get("max_lambda_iters")->As<int>();

    // it doesn't harm to read the parameters for all types!
    if(pn->Has(type.ToString(FRAMED)))
    {
      PtrParamNode t = pn->Get(type.ToString(FRAMED));
      // there are defaults in XML
      always_enlarge_   = t->Get("alwaysEnlarge")->As<bool>();
      start_lower_      = t->Get("lower")->As<double>();
      start_upper_      = t->Get("upper")->As<double>();
      enlarge_lower_    = t->Get("enlargeLower")->As<double>();
      enlarge_upper_    = t->Get("enlargeUpper")->As<double>();
      check_stalled_err_= t->Get("checkErr")->As<bool>();

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
  }

  // some plausibility about optimality condition
  if(type_ != EXTREMIZE && optimization->constraints.view->GetNumberOfActiveConstraints() != 1)
    throw Exception("optimality condition is only possible with exactly one constraint");
  if(type_ == EXTREMIZE && optimization->constraints.view->GetNumberOfActiveConstraints() > 0)
    throw Exception("Optimality Condition 'extremize' is not implemented for constraints");

  LOG_TRACE(oc) << "OptimalityCondition of type " << type.ToString(type_);

  vault_.Resize(optimization->GetDesign()->data.GetSize());
  evaluate_tmp_.Resize(optimization->GetDesign()->data.GetSize());
  
  optimizer_timer_->Stop();

  PostInitScale(1.0, true);
}

void OptimalityCondition::SolveProblem()
{
  // we measure the optimizer in the loops only
  optimizer_timer_->Stop();

  // solve the state problem first
  Optimization::context->pde->GetAssemble()->SetAllReassemble(); // tell assemble that the Design has changed
  optimization->SolveStateProblem();

  // start with iteration 0 which is the initial design
  int iter = 0;
  int max_iter = optimization->GetMaxIterations();
  bool stalled_err = false; // framed only

  while(!optimization->DoStopOptimization() && iter <= max_iter)
  {
    // calc gradients to store the results in data[element]...
    // only the gradients are needed for the calculation of the next iteration

    optimization->SolveAdjointProblems();

    eval_grad_obj_timer_->Start();
    optimization->CalcObjectiveGradient(NULL);
    eval_grad_obj_timer_->Stop();
    
    // reset values of the constraint gradients
    optimization->GetDesign()->Reset(DesignElement::CONSTRAINT_GRADIENT, DesignElement::DEFAULT);

    if(optimization->constraints.view->GetNumberOfActiveConstraints() > 0) {
      eval_grad_const_timer_->Start();
      optimization->CalcConstraintGradient(NULL);
      eval_grad_const_timer_->Stop();
    }

    // store iteration 0
    if(iter == 0)
    {
      eval_obj_timer_->Start();
      optimization->CalcObjective();   // for output
      eval_obj_timer_->Stop();
      // the gradients are (here only!) pointing to the next design vector,
      // hence the gradients for iteration "0" and 1 are identical
      optimization->CommitIteration(); // don't assert we are running
      iter++;
      continue; // redo gradients and start optimization
    }

    optimizer_timer_->Start();
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

    // calc the objective for the logging in CommitIteration(),
    // for the optimality condition it is not required.

    eval_obj_timer_->Start();
    optimization->CalcObjective();
    eval_obj_timer_->Stop();

    // every state problem is an iteration
    // The gradients "point" to this design vector.
    optimization->CommitIteration(); // don't assert we are running
    iter++;
  }
  
  PtrParamNode in = optimization->optInfoNode->Get(ParamNode::SUMMARY)->Get("break");
  
  if(iter >= max_iter-1) 
  {
    in->Get("converged")->SetValue("no");
    in->Get("reason/msg")->SetValue("Maximum iterations exceeded");
    std::cout << "++ max iterations reached" << std::endl;
  }
  assert(in->GetChildren().GetSize() > 0);
} 


bool OptimalityCondition::CalcNextFramedIteration(bool last_was_stalled_err)
{
  // we store the current densities in the temp variable. Otherwise we cannot
  // find the proper lambda
  optimization->GetDesign()->WriteDesignToExtern(vault_.GetPointer());

  // set the frame borders
  if(lower_ == upper_) {
    // this is only reached in the first iteration
    lower_ = start_lower_;
    upper_ = start_upper_;
  } else if(last_was_stalled_err) {
    // when the move limit was wrong, lambda can be arbitrary nonsense
    // But we also may not start from start_*_ as then we repeat the situation
    // this heuristic is not really cool - better start feasible or adjust the move_limit
    lower_ = lambda_ * 1e-4;
    upper_ = std::min(2*start_upper_, lambda_ * 1e4);
  } else {
    //FIXME
    assert(always_enlarge_);
    lower_ = lambda_ * (always_enlarge_ ? enlarge_lower_ : 1.0);
    upper_ = lambda_ * (always_enlarge_ ? enlarge_upper_ : 1.0);
  }
  
  LOG_DBG2(oc) << "CNFI: w_s=" << last_was_stalled_err << " l=" << lower_ << " u=" << upper_;

  // we count lambda iterations to handle the problem of coming too close to a boundary
  lambda_iters_ = 0;
  double err = -1;
  double last_err = -2;   // when move limit is a bound, stop, when err is not changing
  bool stalled_err = false; // when | err - last_err | is zero plus other conditions

  do
  {
    // calc next lambda
    lambda_ = 0.5 * (upper_ + lower_);

    // evaluate with new lambda 
    // restore original density from temp so we always start the calculation 
    // on the same base but with different lambda
    optimization->GetDesign()->ReadDesignFromExtern(vault_.GetPointer());

    last_err = err;
    err = Evaluate(lambda_);

    // when delta err does not change, but err exits
    // and also not to early, as initially err might stay constant in a valid setting (-> Technical/mech_2d has twice err=.2)
    stalled_err = check_stalled_err_ && lambda_iters_ > 5 && abs(err) > feasibility_ && abs(err - last_err) < .5 * feasibility_;
    // move frames according to new lambda
    if(err > 0) upper_ = lambda_; // = center
           else lower_ = lambda_;
    
    lambda_iters_++;

    LOG_DBG2(oc) << "CNFI: li=" << lambda_iters_ << " l=" << lambda_ << " e=" << err << " lo=" << lower_ << " up=" << upper_ << " de=" << abs(err - last_err);
   }
   while(abs(err) > feasibility_  && lambda_iters_ < max_lambda_iters_ && !stalled_err && abs(lambda_) > abs(lambda_min_));

   if(abs(lambda_) < abs(lambda_min_))
     std::cout << "Iteration fails due to too small lambda, check feasibility, move_limit or other optimalityCondition@type than 'framed'" << std::endl;

   if(stalled_err) {
     std::cout << "Iteration fails due to too stalled error: enlarge move_limit (" << move_limit_ << ") or disable optimalityCondition/framed@checkErr" << std::endl;
   }
   if(lambda_iters_ >= max_lambda_iters_)
     std::cout << "Iteration fails to find valid Lagrangian: " << lambda_ << " err: " << err << " check move_limit or bounds in 'optimalityCondition/framed'\n";

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
    if(t < min_err) {
//    if(abs(t) < abs(min_err)) {
      fumble = -1.0 * contract_;
      min_err = t;
    }

    // check with lambda_ + contract_ * step_
    optimization->GetDesign()->ReadDesignFromExtern(vault_.GetPointer());
    t = Evaluate(lambda_ + contract_ * step_);
    if(t < min_err) {    
    //if(abs(t) < abs(min_err)) {
      fumble = contract_;
      min_err = t;
    }

    // check lambda_ + expand_ * step_
    optimization->GetDesign()->ReadDesignFromExtern(vault_.GetPointer());    
    t = Evaluate(lambda_ + expand_ * step_);
    if(t < min_err) {
    //if(abs(t) < abs(min_err)) {
      fumble = expand_;
      min_err = t;
    }
    
    // new lambda:
    lambda_ += fumble * step_;
    step_ = abs(fumble * step_);
    
    lambda_iters_++;

    LOG_DBG2(oc) << "lambda_iter/lambda/err/step/fumble/-ext/-cont/+con/+ext = " <<  lambda_iters_ << "\t" 
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

  LOG_DBG(oc) << "CalcNextIteration: lambda= " << lambda_ << " -> err=" << err;
  
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
      LOG_DBG(oc) << "swap lambda sign to " << lambda_ << " at iter/err: " << lambda_iters_ << "\t" << "\t" << err;
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
      LOG_DBG2(oc) << "lambda/err/l_iter = " <<  lambda_<< "\t" << err << "\t" << lambda_iters_ 
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
  for(Integer i = 0; i < (Integer) data.GetSize(); i++)    
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
    if(next <= lower) evaluate_tmp_[i] =lower;
    if(upper <= next) evaluate_tmp_[i] =upper;
    
    LOG_DBG3(oc) << "CalcNextExtremizeIteration:" << de->elem->elemNum << " obj_grad=" << obj_grad
                 << "(" << de->GetValue(DesignElement::COST_GRADIENT, DesignElement::PLAIN) << ")" << " old= " << rho_e
                 << " next=" << next << " lower=" << lower << " upper=" << upper << " new=" << evaluate_tmp_[i];
  }
  
  // store the new values in the design variables
  optimization->GetDesign()->ReadDesignFromExtern(evaluate_tmp_.GetPointer());
}


double OptimalityCondition::Evaluate(double lambda)
{
   // we assume DensityElement.objective_gradient to be set
   // we assume DensityElement.constraint_gradient to be set
  
   if(abs(lambda) < abs(lambda_min_))
   {
     double org_lambda = lambda;
     lambda = lambda < 0 ? -1.0 * lambda_min_ : lambda_min_;
     std::cout << "Optimality Condition evaluates with too small lambda " 
               << org_lambda << " adjust to " << lambda << std::endl;
     LOG_DBG(oc) << "Evaluate: adjust " << org_lambda << " to " << lambda;
   }

   ConditionContainer& cc = optimization->constraints;
   assert(cc.view->GetNumberOfActiveConstraints() == 1);
   Condition* g = cc.view->Get(0);
   StdVector<DesignElement>& data = optimization->GetDesign()->data;

   // we cannot set the design directly otherwise the filter does not 
   // work (it becomes unsymmetrically for symmetric problems) as all 
   // elements but the first have old and new elements in their filter
   // stencil. Hence we store in evaluate_tmp_

#pragma omp parallel for num_threads(CFS_NUM_THREADS)
   for(Integer i = 0; i < (Integer) data.GetSize(); i++)    
   {
     DesignElement* de = &data[i];
     // rho_e is the old rho
     double rho_e = de->GetDesign(DesignElement::PLAIN);   
    
     // if filter is enabled we use the filtered value otherwise the plain one
     double smart_obj_grad = de->GetValue(DesignElement::COST_GRADIENT, DesignElement::SMART);
     smart_obj_grad *= optimization->objectives.DoMaximize() ? -1.0 : 1.0;
     double b_e = -1.0 * smart_obj_grad;

     // ill posed problems have a problem here!  
     if((boost::math::isnan)(b_e)) EXCEPTION("b_e is nan");

     // for compliant mechanism the gradient can be positive, this is cut
     // -> Bendsoe/Sigmund. p 97
     // for piezo we might become negative lambdas -> cut the positive!
     b_e = lambda >= 0.0 ? std::max(0.0, b_e) : std::min(0.0, b_e);
     
     b_e /= (lambda * de->GetValue(DesignElement::CONSTRAINT_GRADIENT, DesignElement::SMART, g));
     
     // next is density times b_e which is compared with box constraints and move limit
     double next = rho_e * std::pow(b_e, oc_damping_);

     double lower = std::max(de->GetLowerBound(), rho_e - move_limit_);
     double upper = std::min(de->GetUpperBound(), rho_e + move_limit_);            

     // we cannot set the design directly - otherwise the filter stencil get violated 
     evaluate_tmp_[i] = next;
     if(next <= lower) evaluate_tmp_[i] =lower;
     if(upper <= next) evaluate_tmp_[i] =upper;
     
     LOG_DBG3(oc) << "Evaluate:" << de->elem->elemNum << " obj_grad=" << smart_obj_grad
                  << "(" << de->GetValue(DesignElement::COST_GRADIENT, DesignElement::PLAIN) << ")"
                  << " const_grad=" << de->GetValue(DesignElement::CONSTRAINT_GRADIENT, DesignElement::SMART, g)
                  << " old= " << rho_e << " next=" << next << " lower=" << lower
                  << " upper=" << upper << " new=" << evaluate_tmp_[i];
   }
   optimizer_timer_->Stop();

   eval_const_timer_->Start();

   // store the new values in the design variables
   optimization->GetDesign()->ReadDesignFromExtern(evaluate_tmp_.GetPointer());
   
   double vol = optimization->CalcConstraint(g);
   eval_const_timer_->Stop();

   optimizer_timer_->Start();

   double err = g->GetBoundValue() - vol;

   LOG_DBG2(oc) << "E: lambda=" << lambda << " v=" << vol << " e=" << err;

   return err;
}

void OptimalityCondition::LogFileHeader(Optimization::Log& log)
{
  log.AddToHeader("lambda");
  log.AddToHeader("lambda_iters");

  if(type_ == FRAMED) {
    log.AddToHeader("lower");
    log.AddToHeader("upper");
  }

  if(type_ == FUMBLE)
    log.AddToHeader("step");

}


void OptimalityCondition::LogFileLine(std::ofstream* out, PtrParamNode iteration)
{
  if(out)
    *out << " \t" << lambda_ << " \t" << lambda_iters_;

  iteration->Get("lambda")->SetValue(lambda_);
  iteration->Get("lambda_iters")->SetValue(lambda_iters_);
  
  if(type_ == FRAMED)
  { 
    if(out)
      *out << " \t" << lower_ << " \t" << upper_;
    iteration->Get("lower")->SetValue(lower_);
    iteration->Get("upper")->SetValue(upper_);
  }
  if(type_ == FUMBLE)
  { 
    if(out)
      *out << " \t" << step_;
    iteration->Get("step")->SetValue(step_);
  }
}
