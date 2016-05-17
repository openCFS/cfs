#include <assert.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/Logging/log.hpp"
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

DECLARE_LOG(oc)
DEFINE_LOG(oc, "optimalityCondition")

OptimalityCondition::OptimalityCondition(Optimization* optimization, PtrParamNode pn)
 : BaseOptimizer(optimization, pn, Optimization::OPTIMALITY_CONDITION)
{
  type.SetName("OptimalityCondition::Type");
  type.Add(FRAMED, "framed");
  type.Add(FUMBLE, "fumble");
  type.Add(TRAJECTORY, "trajectory");
  type.Add(EXTREMIZE, "extremize");
  
  this->lambda_ = 1000; // just to set a value

  // the following values are standard in mech SIMP -> see e.g. the 99 lines paper
  this->move_limit_ = 0.2;
  this->oc_damping_ = 0.5;
  this->lambda_min_ = 1e-30;
  this->lambda_iters_ = 0;
  this->max_lambda_iters_ = 70;
  this->err_eps_    = 1e-3;
  this->type_       = optimization->objectives.Has(Objective::COMPLIANCE) ? FRAMED : FUMBLE;

  // framed
  this->upper_ = 0.0;
  this->lower_ = 0.0;
  this->start_lower_ = 0;
  this->start_upper_ = 1000;
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
    move_limit_ = pn->Get("move_limit")->As<Double>();
    oc_damping_ = pn->Get("damping")->As<Double>();
    lambda_min_ = pn->Get("lambda_min")->As<Double>();
    err_eps_    = pn->Get("err_eps")->As<Double>();
    max_lambda_iters_ = pn->Get("max_lambda_iters")->As<int>();

    // it doesn't harm to read the parameters for all types!
    if(pn->Has(type.ToString(FRAMED)))
    {
      PtrParamNode t = pn->Get(type.ToString(FRAMED));
      // there are defaults in XML
      always_enlarge_ = t->Get("alwaysEnlarge")->As<bool>();
      start_lower_    = t->Get("lower")->As<Double>();
      start_upper_    = t->Get("upper")->As<Double>();
      enlarge_lower_  = t->Get("enlargeLower")->As<Double>();
      enlarge_upper_  = t->Get("enlargeUpper")->As<Double>();

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
      
      step_     = t->Get("step")->As<Double>();
      contract_ = t->Get("contract")->As<Double>();
      expand_   = t->Get("expand")->As<Double>();
      
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
  
  PostInitScale(1.0, true);
}

void OptimalityCondition::SolveProblem()
{
  // solve the state problem first
  Optimization::context->pde->GetAssemble()->SetAllReassemble(); // tell assemble that the Design has changed
  optimization->SolveStateProblem();

  // start with iteration 0 which is the initial design
  int iter = 0;
  int max_iter = optimization->GetMaxIterations();
  
  while(!optimization->DoStopOptimization() && iter <= max_iter)
  {
    // calc gradients to store the results in data[element]...
    // the gradients are based for the calculation of the next iteration
    optimization->SolveAdjointProblems();
    optimization->CalcObjectiveGradient(NULL);
    
    // reset values of the constraint gradients
    optimization->GetDesign()->Reset(DesignElement::CONSTRAINT_GRADIENT, DesignElement::DEFAULT);
    
    if(optimization->constraints.view->GetNumberOfActiveConstraints() > 0)
      optimization->CalcConstraintGradient(NULL);
    
    // store iteration 0
    if(iter == 0)
    {
      optimization->CalcObjective();   // for output
      // the gradients are (here only! )pointing to the next design vector, 
      // hence the gradients for iteration "0" and 1 are identical
      optimization->CommitIteration(); 
      iter++;
      continue; // redo gradients and start optimization
    }
    
    
    // do a SIMP Optimality Condition step -> calc new design vector
    switch(type_)
    {
    case FRAMED:     CalcNextFramedIteration();
                     break;

    case FUMBLE:     CalcNextFumbleIteration();
                     break;

    case TRAJECTORY: CalcNextTrajectoryIteration();
                     break;

    case EXTREMIZE:  CalcNextExtremizeIteration();
                     break;
                     
    default: assert(false); 
    }
    
    // solve the state problem for the new design vector
    Optimization::context->pde->GetAssemble()->SetAllReassemble();
    optimization->SolveStateProblem();

    // calc the objective for the logging in CommitIteration(),
    // for the optimality condition it is not required.
    optimization->CalcObjective();
    
    // every state problem is an iteration 
    // The gradients "point" to this design vector. 
    optimization->CommitIteration();
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


void OptimalityCondition::CalcNextFramedIteration()
{
  // we store the current densities in the temp variable. Otherwise we cannot
  // find the proper lambda
  optimization->GetDesign()->WriteDesignToExtern(vault_.GetPointer());

  // set the frame borders, first iteration when lower_ == upper_
  lower_ = lower_ == upper_ ? start_lower_ : lambda_ * (always_enlarge_ ? enlarge_lower_ : 1.0);
  upper_ = lower_ == upper_ ? start_upper_ : lambda_ * (always_enlarge_ ? enlarge_upper_ : 1.0);
  
  // we count lambda iterations to handle the problem of coming too close to a boundary
  lambda_iters_ = 0;
  int count = 0;
  double err;
   
  do
  {
    // check if we have to enlarge
    if(++count > max_lambda_iters_)
    {
      // enlarge only where we need to 
      if(abs(upper_ - lambda_) < abs(lambda_ - lower_))
        upper_ = lower_ == upper_ ? start_upper_ : lambda_ * enlarge_upper_;
      else
        lower_ = lower_ == upper_ ? start_lower_ : lambda_ * enlarge_lower_;
      
      LOG_DBG(oc) << "enlarge after " << count << " iterations: lower=" << lower_ << " upper=" << upper_;
      count = 0;
    }

    // calc next lambda
    lambda_ = 0.5 * (upper_ + lower_);

    // evaluate with new lambda 
    // restore original density from temp so we always start the calculation 
    // on the same base but with different lambda
    optimization->GetDesign()->ReadDesignFromExtern(vault_.GetPointer());

    err = Evaluate(lambda_);
    // move frames according to new lambda
    if(err > 0) upper_ = lambda_; // = center
           else lower_ = lambda_;
    
    lambda_iters_++;

    LOG_DBG2(oc) << "lambda_iter/lambda/err/lower/upper = " <<  lambda_iters_ << "\t" 
                 << lambda_ << "\t" << err << "\t" << lower_ << "\t" << upper_; 
   }
   while(abs(err) > err_eps_  && lambda_iters_ < max_lambda_iters_);
  
   if(lambda_iters_ >= max_lambda_iters_)
     std::cout << "Iteration fails to find valid lagrangian: " << lambda_ << " err: " << err << std::endl;
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
  while(abs(min_err) > err_eps_ && lambda_iters_ < max_lambda_iters_);

  if(lambda_iters_ >= max_lambda_iters_)
    std::cout << "Iteration fails to find valid lagrangian: " << lambda_ << " err: " << min_err << std::endl;
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

  while(abs(err) > err_eps_ && lambda_iters_ < max_lambda_iters_)
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
    std::cout << "Iteration fails to find valid lagrangian: " << lambda_ << " err: " << err << std::endl;
}

void OptimalityCondition::CalcNextExtremizeIteration()
{
  StdVector<DesignElement>& data = optimization->GetDesign()->data;

  // we cannot set the design directly otherwise the filter does not 
  // work (it becomes unsymmetrically for symmetric problems) as all 
  // elements but the first have old and new elements in their filter
  // stencil. Hence we store in evaluate_tmp_
  
  for(unsigned int i = 0; i < data.GetSize(); i++)    
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
   
   for(unsigned int i = 0; i < data.GetSize(); i++)    
   {
     DesignElement* de = &data[i];
     // rho_e is the old rho
     double rho_e = de->GetDesign(DesignElement::PLAIN);   
    
     // if filter is enabled we use the filtered value otherwise the plain one
     double smart_obj_grad = de->GetValue(DesignElement::COST_GRADIENT, DesignElement::SMART);
     double b_e = -1.0 * smart_obj_grad;

     // ill posed problems have a problem here!  
     if(std::isnan(b_e)) EXCEPTION("b_e is nan");

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
   
   // store the new values in the design variables
   optimization->GetDesign()->ReadDesignFromExtern(evaluate_tmp_.GetPointer());
   
   double vol = optimization->CalcConstraint(g);
   double err = g->GetBoundValue() - vol;
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
