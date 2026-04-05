#include <assert.h>
#include <stddef.h>
#include <iostream>
#include <limits>
#include <string>

#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Domain/Domain.hh"
#include "Domain/ElemMapping/Elem.hh"
#include "General/Enum.hh"
#include "General/defs.hh"
#include "General/Environment.hh"
#include "General/Exception.hh"
#include "MatVec/Vector.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Optimization.hh"
#include "Optimization/Optimizer/GradientCheck.hh"
#include "Optimization/Excitation.hh"
#include "PDE/BasePDE.hh"
#include "Utils/StdVector.hh"
#include "Utils/ToolsFull.hh"

using namespace CoupledField;

EXTERN_LOG(optimizer)

GradientCheck::GradientCheck(Optimization* optimization, PtrParamNode pn) :
  BaseOptimizer(optimization, pn, Optimization::GRADIENT_CHECK)
{
  Context* ctxt = Optimization::context;
  assert(ctxt->GetExcitation() != NULL);
  ctxt->GetExcitation()->reassemble = true;

  opt_timer->Stop(); // we don't spend time here

  // reduce to our actual ParamNode
  pn = pn->Get(Optimization::optimizer.ToString(Optimization::GRADIENT_CHECK),
      ParamNode::PASS);

  h = pn != NULL ? pn->Get("h")->As<Double>() : 1e-6;
  if (h <= 0.0)
    throw Exception("gradient check 'h' needs to be non-negative");
  order = pn != NULL ? (pn->Get("order")->As<std::string>() == "first" ? 1 : 2) : 1;
  design = optimization->GetDesign();

  DesignElement& de = design->data[0];

  // here we modify the variables */
  design_values.Resize(design->GetNumberOfVariables()); // shall include also slack variables, ...
  design->WriteDesignToExtern(design_values.GetPointer());

  finite_diff_result_index_ = design->GetSpecialResultIndex(de.DEFAULT, de.COST_GRADIENT, de.FINITE_DIFF_COST_GRADIENT, de.PLAIN);
  error_result_index_ = design->GetSpecialResultIndex(de.DEFAULT, de.COST_GRADIENT, de.ERROR_COST_GRADIENT, de.PLAIN);

  PostInitScale(1.0, true);

  info_->Get(ParamNode::HEADER)->Get("type")->SetValue(Optimization::optimizer.ToString(Optimization::GRADIENT_CHECK));
  info_->Get(ParamNode::HEADER)->Get("h")->SetValue(h);
  info_->Get(ParamNode::HEADER)->Get("order")->SetValue(order == 1 ? "first" : "second");
}

void GradientCheck::SolveProblem()
{
  std::cout << "Check gradients for initial guess via finite differences ("
      << (order == 1 ? "first" : "second") << " order)" << std::endl;

  if (finite_diff_result_index_ == -1 || error_result_index_ == -1)
    info_->SetWarning("Not all results defined for finite difference (value='costGradient' detail='...')");

  // solve the original problem once!!
  opt_timer->Stop();
  optimization->SolveStateProblem();
  double curr_obj = optimization->CalcObjective();
  optimization->SolveAdjointProblems();
  eval_grad_obj_timer_->Start();
  optimization->CalcObjectiveGradient(NULL);
  eval_grad_obj_timer_->Stop();
  opt_timer->Start();

  // store here the finite difference value
  Vector<double> finite;
  finite.Resize(design->data.GetSize());

  // save here the analytical cost gradient
  Vector<double> analytical;
  analytical.Resize(design->data.GetSize());

  // relative error
  Vector<double> error;
  error.Resize(design->data.GetSize());

  PtrParamNode in = info_->Get(ParamNode::PROCESS)->Get("gradientCheck");
  in->Get("total")->SetValue(design->data.GetSize());
  PtrParamNode cg = in->Get("costGradient", ParamNode::APPEND);

  Context* ctxt = Optimization::context;
  AnalysisID& id = ctxt->GetDriver()->GetAnalysisId();

  // perform a lot of calculations
  for (unsigned int i = 0; i < design->data.GetSize(); i++)
  {
    id.step = i+1;

    DesignElement* de = &design->data[i];

    // expensive!!
    // store only the latest so progress can get checked
    finite[i] = PerformFiniteDifferenceEval(de, curr_obj, cg);
    analytical[i] = de->GetValue(DesignElement::COST_GRADIENT, DesignElement::SMART);
    error[i] = finite[i] != 0.0 ? ((finite[i] - analytical[i]) / finite[i]) : std::numeric_limits<double>::max();

    cg->Get("id")->SetValue(i+1);
    cg->Get("grad")->SetValue(analytical[i]);
    cg->Get("fd_grad")->SetValue(finite[i]);
    cg->Get("diff")->SetValue(finite[i] - analytical[i]);
    if (finite[i] != 0.0) {
      cg->Get("rel_error")->SetValue((finite[i] - analytical[i]) / finite[i]);
      cg->Get("factor")->SetValue(analytical[i] / finite[i]);
    } else {
      cg->Get("rel_error")->SetValue("finite difference is zero");
      cg->Get("factor")->SetValue("finite difference is zero");
    }
  }

  // write the stuff -> The PDE solution is slighly wrong!
  CommitIteration();

  // finish comparison
  double l2 = error.NormL2() / error.GetSize();
  double max = error.NormMax();

  info_->Get(ParamNode::SUMMARY)->Get("max_relative_error")->SetValue(max);
  info_->Get(ParamNode::SUMMARY)->Get("relative_L2_error")->SetValue(l2);
  std::cout << "relative error -> max=" << max << " L2=" << l2 << std::endl;
}

double GradientCheck::PerformFiniteDifferenceEval(DesignElement* de, double f_x, PtrParamNode in)
{
  assert(Optimization::context->GetExcitation()->reassemble == true);

  double f_x1 = 0.0, f_x2 = 0.0, fd = 0.0;
  // in first order we always perform a forward fifference quotient
  // if we are to close to the max value we step back.
  // Two things have to be considered:
  // - make sure f_x can be reused
  // - take care with forward and backward differences for this case.

  double x_org = de->GetDesign(DesignElement::PLAIN);
  assert(x_org == design_values[de->GetIndex()]);

  // consider first order
  double x_eval_1 = x_org + h <= de->GetUpperBound() ? x_org + h : x_org - h;
  int x_dir_1 = x_org + h <= de->GetUpperBound() ? +1 : -1;

  // x_dir_1 = 1  : x_org .. x_eval_1 .. x_max
  // x_dir_1 = -1 : x_eval_1 .. x_org .. x_max

  // in the second order we also have to check against the lower bound
  double x_eval_2 = x_eval_1 - (x_dir_1 == 1 ? 2 : 1) * h;
  int x_dir_2 = +1;
  if (order == 2 && x_eval_2 < de->GetLowerBound())
  {
    // too far left
    assert(x_dir_1 == 1); // then there was no problem to the right side
    x_eval_2 = x_org + h;
    x_eval_1 = x_org + 2 * h;
    x_dir_2 = -1;
  }

  // x_dir_2 = 1  : x_min .. x_eval_2 .. x_org .. x_eval_1
  // x_dir_2 = -1 : x_min .. x_org .. x_eval_2 .. x_eval_1    

  // forward  difference quotient: D_+(x_0) = h^-1 ( f(x_0+h)-f(x_0) )
  // backward difference quotient: D_-(x_0) = h^-1 ( f(x_0)-f(x_0-h) )
  // central difference quotient: D(x_0)   = 2 h^-1 ( f(x_0+h) - f(x_0-h) )

  // note! in second order case we can do only one evaluation for x_dir_1/2 = -1!
  if (order == 2)
  {
    if (x_dir_1 == -1)
      x_eval_1 = std::numeric_limits<double>::max(); // we may not calculate this!
    if (x_dir_2 == -1)
      x_eval_2 = std::numeric_limits<double>::min(); // we may not calculate this!
  }

  // calc forward value not with degenerated second order case
  if (order == 1 || x_dir_1 == 1)
  {
    design_values[de->GetIndex()] = x_eval_1;
    design->ReadDesignFromExtern(design_values); // communicate the change

    LOG_DBG(optimizer) << "PFDE: " << design->ToString(1);

    opt_timer->Stop();
    optimization->SolveStateProblem();
    optimization->SolveAdjointProblems();
    f_x1 = optimization->CalcObjective();
    opt_timer->Start();
  }

  // do not calc degenerated second order case
  assert(!(x_dir_1 == -1 && x_dir_2 == -1));
  if (order == 2 && x_dir_2 == 1)
  {
    design_values[de->GetIndex()] = x_eval_2;
    design->ReadDesignFromExtern(design_values); // communicate the change

    opt_timer->Stop();
    optimization->SolveStateProblem(); // expensive
    optimization->SolveAdjointProblems();
    f_x2 = optimization->CalcObjective();
    opt_timer->Start();
  }

  // reset design
  design_values[de->GetIndex()] = x_org;
  design->ReadDesignFromExtern(design_values); // communicate the change

  if (order == 1)
  {
    fd = x_dir_1 == 1 ? (f_x1 - f_x) / h : (f_x - f_x1) / h;
  }
  if (order == 2)
  {
    if (x_dir_1 == 1 && x_dir_2 == 1)
    {
      assert(close(x_eval_1 - x_eval_2, 2 * h));
      fd = (f_x1 - f_x2) / (2 * h);
    }
    if (x_dir_1 == -1)
    {
      assert(close(x_org - x_eval_2, 2 * h));
      fd = (f_x - f_x2) / (2 * h);
    }
    if (x_dir_2 == -1)
    {
      assert(close(x_eval_1 - x_org, 2 * h));
      fd = (f_x1 - f_x) / (2 * h);
    }
  }

  in->Get("actual_order")->SetValue(order == 1 || x_dir_1 == -1 || x_dir_2 == -1 ? "first" : "second");

  LOG_DBG(optimizer)
    <<  "PFDE: elem=" << de->elem->elemNum
    << " x_org=" << x_org << " order=" << order << " f_x=" << f_x
    << " x_eval_1=" << x_eval_1 << " x_dir_1=" << x_dir_1 << " f_x1=" << f_x1
    << " x_eval_2=" << x_eval_2 << " x_dir_2=" << x_dir_2 << " f_x2=" << f_x2
    << " fd=" << fd << " f_grad=" << de->GetValue(DesignElement::COST_GRADIENT, DesignElement::PLAIN);

  // eventually store value
  if(finite_diff_result_index_ != -1)
    de->specialResult[finite_diff_result_index_] = fd;
  if(error_result_index_ != -1)
    de->specialResult[error_result_index_] = (fd - de->GetValue(DesignElement::COST_GRADIENT, DesignElement::PLAIN)) / fd;

  return fd;
}
