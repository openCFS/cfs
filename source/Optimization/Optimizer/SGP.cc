#include "Optimization/Optimizer/SGP.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "MatVec/Matrix.hh"
#include "Utils/tools.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/Logging/log.hpp"
#include <math.h>
#include "Optimization/TransferFunction.hh"
#include "PDE/StdPDE.hh"
#include "Optimization/OptimizationMaterial.hh"
#include "Domain/CoefFunction/CoefFunctionOpt.hh"
#include "Driver/Assemble.hh"
#include "Optimization/Optimization.hh"
#include <limits>

DECLARE_LOG(sgp)
DEFINE_LOG(sgp, "sgp")

using namespace CoupledField;
using std::pow;
using std::max;
using std::min;
using std::abs;
using std::string;

SGP::SGP(Optimization* opt, PtrParamNode pn) : BaseOptimizer(opt, pn, Optimization::SGP_SOLVER)
{
  obj = NULL;
  n = 0;
  m = 0;

  pmin_vol   = this_opt_pn_->Get("pmin_vol")->As<double>();
  pmin_filt   = this_opt_pn_->Get("pmin_filt")->As<double>();
  pmax_vol   = this_opt_pn_->Get("pmax_vol")->As<double>();
  pmax_filt   = this_opt_pn_->Get("pmax_filt")->As<double>();
  lbound = this_opt_pn_->Get("lower_bound")->As<double>();
  ubound = this_opt_pn_->Get("upper_bound")->As<double>();
  tolerance = this_opt_pn_->Get("tolerance")->As<double>();
  tau = this_opt_pn_->Get("tau")->As<double>();
  derivative_check = this_opt_pn_->Get("derivative_check")->As<bool>();


  ppen_vol = 0.;
  ppen_filt = 0.;
  ppeni = 0.;
  pmaxi = 0.;
  pmini = 0.;
  bisect = this_opt_pn_->Get("bisect_iter")->As<int>();

  merit = 0;

  //pmin = 0.05;
  //pmax = 0.05;

  tf = optimization->GetDesign()->GetTransferFunction(DesignElement::DENSITY, App::MECH);


  //asymptotes.SetName("SGP::Asypmtotes");
  //asymptotes.Add(FIXED, "fixed");
  //asymptotes.Add(MMA, "mma");

  // setup design variables
  n_elem = optimization->GetDesign()->GetNumberOfElements();
  rho_outer.Resize(n_elem);
  theta_outer.Resize(n_elem);
  s1_outer.Resize(n_elem);
  s2_outer.Resize(n_elem);
  E_outer.Resize(n_elem);
  E_inner.Resize(n_elem);

  PostInitScale(1.0);
}

SGP::~SGP()
{
}

void SGP::PostInit()
{
  assert(obj == NULL); // call once
  assert(optimization->objectives.data.GetSize() == 1); // trivial case only

  n = optimization->GetDesign()->GetNumberOfVariables();

  x_outer.Resize(n, 0.0);

  // setup functions
  //Function* f = optimization->objectives.data[0];
  obj = new SGPApproximation(this);
  obj->PostInit();

  /** Setup material tensor E_0 */
  E_0.Resize(3,3);
  SinglePDE * pde = Optimization::context->pde;
  BiLinFormContext* c = pde->GetAssemble()->GetBiLinForm("LinElastInt", optimization->GetDesign()->GetRegionId(), pde, pde, false);
  shared_ptr<CoefFunctionOpt> coef = Optimization::context->mat->GetMatCoef("LinElastInt", c);
  LocPointMapped lpm;
  coef->orgMat->GetTensor(E_0,lpm);

  /** Setup inner variable material Tensor E_inner */
  DesignSpace* space = optimization->GetDesign();
  unsigned int mech11 = space->FindDesign(DesignElement::MECH_11);
  unsigned int mech12 = space->FindDesign(DesignElement::MECH_12);
  unsigned int mech13 = space->FindDesign(DesignElement::MECH_13);
  unsigned int mech22 = space->FindDesign(DesignElement::MECH_22);
  unsigned int mech23 = space->FindDesign(DesignElement::MECH_23);
  unsigned int mech33 = space->FindDesign(DesignElement::MECH_33);
  for (unsigned int i = 0; i < n_elem; i++) {
    E_inner[i].Resize(3,3);
    E_inner[i][0][0] = space->GetDesignElement(mech11*n_elem+i)->GetDesign(DesignElement::PLAIN);
    E_inner[i][0][1] = space->GetDesignElement(mech12*n_elem+i)->GetDesign(DesignElement::PLAIN);
    E_inner[i][0][2] = space->GetDesignElement(mech13*n_elem+i)->GetDesign(DesignElement::PLAIN);
    E_inner[i][1][1] = space->GetDesignElement(mech22*n_elem+i)->GetDesign(DesignElement::PLAIN);
    E_inner[i][1][2] = space->GetDesignElement(mech23*n_elem+i)->GetDesign(DesignElement::PLAIN);
    E_inner[i][2][2] = space->GetDesignElement(mech33*n_elem+i)->GetDesign(DesignElement::PLAIN);
    E_inner[i][1][0] = E_inner[i][0][1];
    E_inner[i][2][0] = E_inner[i][0][2];
    E_inner[i][2][1] = E_inner[i][1][2];
  }

  // setup lower and upper bounds, they might be from design bounds. After this we must not use DesignElement::GetLower/UpperBound() !!
  /**lower_bound.Resize((int) DesignElement::ALL_DESIGNS, 0.0);
  upper_bound.Resize((int) DesignElement::ALL_DESIGNS, 0.0);
  for(unsigned int i = 0; i < optimization->GetDesign()->design.GetSize(); i++)
  {
    DesignElement::Type dt = optimization->GetDesign()->design[i].design;
    assert(dt >= 0);
    DesignElement& de = optimization->GetDesign()->data[optimization->GetDesign()->elements * i]; // fallback
    assert(de.GetType() == dt);
    Condition* g = optimization->constraints.Get(Function::DESIGN_BOUND, dt, Condition::LOWER_BOUND, false);
    assert(!(g != NULL && g->GetBoundValue() < de.GetLowerBound())); // the design bound shall be -inf then!
    lower_bound[dt] = g != NULL ? g->GetBoundValue() : de.GetLowerBound();
    g = optimization->constraints.Get(Function::DESIGN_BOUND, dt, Condition::UPPER_BOUND, false);
    assert(!(g != NULL && g->GetBoundValue() > de.GetUpperBound()));
    upper_bound[dt] = g != NULL ? g->GetBoundValue() : de.GetUpperBound();
    LOG_DBG3(sgp) << "SGP:PI dt=" << dt << "=" << DesignElement::type.ToString(dt) << " lb=" << lower_bound[dt] << " ub=" << upper_bound[dt];
  } */


  // setup lower and upper design bound, steps for inner design variables
  ParamNodeList ivs =  this_opt_pn_->Get("subsolver")->GetList("var");
  inner_variables.Resize(ivs.GetSize());
  for (unsigned int i = 0; i < ivs.GetSize(); i++) {
    inner_variables[i].Read(ivs[i], optimization->GetDesign());
  }

  // set design variable configuration, e.g. DESIGN_ROTANGLE or STIFF1_STIFF2 or FOMO or FMO
  SetConfiguration();
  if (configuration == STIFF1_STIFF2) {
    helper_dm = new DesignMaterial(this_opt_pn_->Get("subsolver/designMaterial"), OptimizationMaterial::system.Parse(optimization->GetDesign()->designMaterial->GetErsatzMaterial()->pn->Get("material")->As<std::string>()), optimization->GetDesign()->design, optimization->GetDesign()->designMaterial->GetErsatzMaterial());
    helper_dm->SetType(DesignMaterial::HOM_RECT_C1);
  }

  // setup SGP asymptotes
  L.Resize(n_elem);
  U.Resize(n_elem);
  for(unsigned int i = 0; i < n_elem; i++)
  {
    L[i].Resize(3,3);
    U[i].Resize(3,3);
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        if (k==j) {
          L[i][j][k] = 1.;
          U[i][j][k] = 1.;
        } else {
          L[i][j][k] = 0.;
          U[i][j][k] = 0.;
        }
      }
    }
    L[i] *= lbound;
    U[i] *= ubound;
  }

  // setup constraints
  ConditionContainer& cc = optimization->constraints;
  m = cc.view->GetNumberOfActiveConstraints();
  constr.Resize(m);
  for(unsigned int i = 0; i < m; i++)
  {
    // this is the only place where we are allowed not to use MMAApproximation::GetCondition()
    Condition* g = cc.view->Get(i);
    //MMAApproximation* func = new MMAApproximation(this, i, false);
    constr[i] = g;
    if (constr[i]->GetType() == Condition::FILTERING_GAP)
      // setup size of the filter constraint derivative
      filter_outer_grad.Resize(obj->outer_grad.GetSize());
  }
  cc.view->Done();

  ToInfo();
}

void SGP::InnerVariable::Read(PtrParamNode pn, DesignSpace* space) {
  type = DesignElement::type.Parse(pn->Get("design")->As<string>());
  int index = space->FindDesign(type,false);
  if (index < 0) {
    EXCEPTION("SGP inner variable '" << pn->Get("design")->As<string>() << "' not defined as design variable.");
  }
  BaseDesignElement* de = space->GetDesignElement(index*space->GetNumberOfElements());
  assert(de->GetType() == type);
  lower_bound = de->GetLowerBound();
  upper_bound = de->GetUpperBound();
  steps = pn->Get("steps")->As<int>();
  inc = (upper_bound - lower_bound)/steps;
}

bool SGP::HasInnerVar(BaseDesignElement::Type tp) const{
  for (unsigned int i = 0; i < inner_variables.GetSize();i++) {
    if (inner_variables[i].type == tp) {
      return true;
    }
  }
  return false;
}

void SGP::SolveProblem()
{
  int max_iter = optimization->GetMaxIterations();
  PtrParamNode summary = optimization->optInfoNode->Get(ParamNode::SUMMARY);

  bool converged = false;
  bool penal_filt = false;
  int iter = 1;

  // initialize merit functions and constraints
  double merit_old = -1.;

  // set penalty parameter
  ppen_vol = 0.5*(pmin_vol+pmax_vol);
  ppen_filt = 0.5*(pmin_filt+pmax_filt);

  // sum of all filtering gaps
  double filtering_gaps = 0;

  // initialize outer gradient
  StdVector<Matrix<double> > df(n_elem);
  for (unsigned int i = 0; i < n_elem;i++) {
    df[i].Resize(3,3);
  }

  // writes design to outer variables rho_outer
  DesignToOuter();
  // updates outer variables rho_outer, writes outer_variables to design, and evaluates objective and constraint functions, gradients
  UpdateToCurrentStep();

  // check derivative if option is turned on by finite differences
  if (derivative_check) {
    StdVector<double> deriv_check(6*n_elem), deriv(6*n_elem);
    double max_error = 0;
    deriv_check = GradientCheck(max_error);
    std::cout<<"Derivative check:"<<std::endl;
    std::cout<<"approx derivative = "<<deriv_check.ToString()<<std::endl;
    GetOuterDerivativeVector(deriv,obj->outer_grad);
    std::cout<<"derivative = "<<deriv<<std::endl;
    std::cout<<"max derivative error = "<<max_error<<std::endl;
  }
  optimization->CommitIteration();

  // outer optimization loop
  while(!optimization->DoStopOptimization() && iter <= max_iter && !converged)
  {
    // saves old merit function value
    merit_old = merit;
    merit = obj->outer_val;

    // calculate objective function including constraint penalty term
    for(unsigned int i = 0; i < m; i++)
    {
      // Add non-physical constraint to merit functions
      if (constr[i]->GetType() == Condition::VOLUME || constr[i]->GetType() == Condition::GLOBAL_TWO_SCALE_VOL) {
        merit += ppen_vol * constr[i]->GetValue()*n_elem;
      } else if (constr[i]->GetType() == Condition::FILTERING_GAP) {
        merit += ppen_filt * constr[i]->GetValue();
      }
    }
    LOG_DBG2(sgp) << "SP: it=" << iter << " merit = " << merit;

    // test if outer iteration progresses
//    if (merit <= merit_old) {
//      worseCounter = 0;
//    } else {
//      worseCounter++;
//      if (worseCounter > 5) {
//        std::cout << "\nSGP stopped due to no progress: " << summary->Get("reason/msg")->As<string>() << std::endl;
//        break;
//      }
//    }
//    std::cout << "\n merit = " << merit << " merit_old = " << merit_old << " worse_counter = " << worseCounter << std::endl;


    // Get outer gradient for subproblem
    GetOuterDerivative(df, obj->outer_grad);
    if (abs(pmin_vol-pmax_vol) < tolerance) {
      // solves sub-model and updates rho_outer and theta_outer
      if (configuration == DENSITY_ROTANGLE) {
        obj->SubSolve_Density_Rotangle(SGPApproximation::FUNC,df,ppen_vol,rho_outer,theta_outer);
      } else if (configuration == STIFF1_STIFF2) {
        obj->SubSolve(SGPApproximation::FUNC,df,ppen_vol,s1_outer,s2_outer,theta_outer);
      } else if (configuration == FOMO) {
        obj->SubSolve_FOMO(SGPApproximation::FUNC,df,ppen_vol,rho_outer,theta_outer);
      }
    } else {
      if (iter < 1000) {
        pmaxi = pmax_vol;
        pmini = pmin_vol;
      } else if (iter < 2000) {
        pmaxi = min(pmax_vol,1.1*ppen_vol);
        pmini = max(pmin_vol,0.9*ppen_vol);
      } else {
        pmaxi = min(pmax_vol, 1.02*ppen_vol);
        pmini = max(pmin_vol, 0.98*ppen_vol);
      }
      ppeni = ppen_vol;
      int ki = 0;
      bool penal_vol = true;
      double cond = 0;
      std::string output_str = "";
      while (ki < bisect && penal_vol) {
        output_str = "s-iteration: "+std::to_string(ki)+" compliance: "+std::to_string(obj->outer_val)+ " merit: "+std::to_string(merit)+" pmaxi: " + std::to_string(pmaxi)+" pmini: "+std::to_string(pmini);
        if (configuration == DENSITY_ROTANGLE) {
          obj->SubSolve_Density_Rotangle(SGPApproximation::FUNC,df,ppeni,rho_outer,theta_outer);
        } else if (configuration == STIFF1_STIFF2) {
          obj->SubSolve(SGPApproximation::FUNC,df,ppeni,s1_outer,s2_outer,theta_outer);
        } else if (configuration == FOMO) {
          obj->SubSolve_FOMO(SGPApproximation::FUNC,df,ppeni,rho_outer,theta_outer);
        }
        // Writes new outer variables to Design
        //OuterToDesign();
        // evaluate all functions such that we have the function values for the current design
        // the subproblem is based on the approximated values only
        // writes design to outer variables rho_outer,... and evaluates objective and constraint functions, gradients
        UpdateToCurrentStep(true);

        // calculate objective function including constraint penalty term
        merit = obj->outer_val;
        filtering_gaps = 0;
        for(unsigned int i = 0; i < m; i++)
        {
          // Add non-physical constraint to merit functions
          if (constr[i]->GetType() == Condition::VOLUME || constr[i]->GetType() == Condition::GLOBAL_TWO_SCALE_VOL) {
            merit += ppeni * constr[i]->GetValue()*n_elem;
            output_str += " volume: "+std::to_string(constr[i]->GetValue()) + " ppeni: "+ std::to_string(ppeni);
          } else if (constr[i]->GetType() == Condition::FILTERING_GAP) {
            merit += ppen_filt * constr[i]->GetValue();
            filtering_gaps += constr[i]->GetValue();
          }
        }
        if (filtering_gaps > 0)
          output_str += " filtering_gaps: "+std::to_string(filtering_gaps) + " ppen_filt: "+ std::to_string(ppen_filt);
        std::cout<<output_str<<std::endl;
        cond = 0;
        if(penal_vol) {
          for(unsigned int i = 0; i < m; i++)
          {
            if (constr[i]->GetType() == Condition::VOLUME || constr[i]->GetType() == Condition::GLOBAL_TWO_SCALE_VOL) {
              cond = constr[i]->GetValue() - constr[i]->GetBoundValue();
            }
          }
          if (abs(cond) < 0.001) {
            penal_vol = false;
          } else {
            if (cond > 0) {
              pmini = ppeni;
              ppeni = 0.5 * (pmaxi + ppeni);
            } else {
              pmaxi = ppeni;
              ppeni = 0.5 * (pmini + ppeni);
            }
          }
        }
        ki++;
      }
      ppen_vol = ppeni;
     }

     // adapt filter gap penalty parameter
     if (abs(filtering_gaps) < 0.001) {
       penal_filt = true;
     } else {
         if (filtering_gaps > 0) {
           pmin_filt = ppen_filt;
           ppen_filt = 0.5 * (pmax_filt + ppen_filt);
         } else {
           pmax_filt = ppen_filt;
           ppen_filt = 0.5 * (pmin_filt + ppen_filt);
         }
         penal_filt = false;
    }

    LOG_DBG2(sgp) << "SP: it=" << iter << " c_grad = [" << obj->outer_grad.ToString() << "]";
    LOG_DBG2(sgp) << "SP: it=" << iter << " x_old = [" << x_outer.ToString() << "]";

    //convergence criterion
    if (abs(merit-merit_old) < tolerance && penal_filt) {
      converged = true;
    }

    // Writes new outer variables to Design
    //OuterToDesign();
    // evaluate all functions such that we have the function values for the current design
    // the subproblem is based on the approximated values only
    // writes design to outer variables rho_outer,... and evaluates objective and constraint functions, gradients
    UpdateToCurrentStep();

    PtrParamNode pn_iter = optimization->CommitIteration();
    ToInfoIter(pn_iter);

    LOG_DBG2(sgp) << "SP: it=" << iter << " x_new_curr = [" << x_outer.ToString() << "]";
    iter = optimization->GetCurrentIteration();
  }

  summary->Get("break/converged")->SetValue(converged ? "yes" : "no");
  if(iter >= max_iter)
    summary->Get("reason/msg")->SetValue("maximum iterations reached");

  if(converged)
    std::cout << "\nSGP converged after " << (iter-1) << " iterations: " << summary->Get("reason/msg")->As<string>() << std::endl;
  else
    std::cout << "\nSGP did not converge: " << summary->Get("reason/msg")->As<string>() << std::endl;
}

StdVector<double> SGP::GradientCheck(double & max_grad_error) {
  DesignSpace* space = optimization->GetDesign();
  space->WriteDesignToExtern(x_outer);
  optimization->GetDesign()->Reset(DesignElement::COST_GRADIENT, DesignElement::DEFAULT);

  // calculate objective and gradient
  StdVector<double> grad(obj->outer_grad.GetSize());
  StdVector<double> constr_grad(m*obj->outer_grad.GetSize());
  EvalGradObjective(n, x_outer.GetPointer(), true, grad);
  EvalGradConstraints(n,x_outer.GetPointer(),m,m * obj->outer_grad.GetSize(),true,true,constr_grad);
  for (unsigned int i = 0;i < obj->outer_grad.GetSize(); i++) {
    for (unsigned int j=0; j < m;j++) {
      if (constr[j]->GetType() == Condition::FILTERING_GAP)
      {
        grad[i] += ppen_filt * constr_grad[j*obj->outer_grad.GetSize() + i];
      }
    }
  }

  // calculate finite difference gradient from objective function with central difference quotient
  StdVector<unsigned int> grad_index(6);
  grad_index[0] = space->FindDesign(DesignElement::MECH_11);
  grad_index[1] = space->FindDesign(DesignElement::MECH_12);
  grad_index[2] = space->FindDesign(DesignElement::MECH_13);
  grad_index[3] = space->FindDesign(DesignElement::MECH_22);
  grad_index[4] = space->FindDesign(DesignElement::MECH_23);
  grad_index[5] = space->FindDesign(DesignElement::MECH_33);

  double eps = 1e-6;
  StdVector<double> diff_grad(grad_index.GetSize() * n_elem), x_outer_right(x_outer.GetSize());
  // Check gradient with central difference quotient
  double obj_right, obj_left, grad_error,count=0;
  StdVector<double> constr_val(m);
  for (unsigned int i = 0; i < n_elem;i++) {
    for (unsigned int j=0; j< grad_index.GetSize();j++) {
      // calculate central difference quotient
        // Eval function value f(x+h) (right)
      unsigned int index = grad_index[j]*n_elem+i;
      assert(index < x_outer.GetSize());
      x_outer[index] += eps;
      obj_right = EvalObjective(n, x_outer.GetPointer(), true);
      EvalConstraints(n,x_outer.GetPointer(),m,true,constr_val.GetPointer(),true);
      for (unsigned int k=0; k < m;k++) {
        if (constr[k]->GetType() == Condition::FILTERING_GAP)
        {
          //int ind = k*obj->outer_grad.GetSize() + i;
          obj_right += ppen_filt * constr_val[k];
        }
      }
      assert(index < x_outer.GetSize());
      x_outer[index] -= 2*eps;

      // Eval function value f(x-h) (left)
      obj_left = EvalObjective(n, x_outer.GetPointer(), true);
      EvalConstraints(n,x_outer.GetPointer(),m,true,constr_val.GetPointer(),true);
      for (unsigned int k=0; k < m;k++) {
        if (constr[k]->GetType() == Condition::FILTERING_GAP)
        {
          obj_left += ppen_filt * constr_val[k];
        }
      }
      //delete constr_val;
      assert(index < x_outer.GetSize());
      x_outer[index] += eps;
      diff_grad[count] = (obj_right-obj_left)/ (2.*eps);
      grad_error = abs(diff_grad[count]-grad[index]);
      count++;
      if (grad_error > max_grad_error)
        max_grad_error = grad_error;
    }
  }
  return diff_grad;
}

void SGP::UpdateToCurrentStep(bool inner)
{
  DesignSpace* space = optimization->GetDesign();
  DesignToOuter(inner,true);
  OuterToDesign();
  space->WriteDesignToExtern(x_outer);

  //LOG_DBG(sgp) << "UTCP x=" << x_outer.ToString();

  if (!inner) {
    // prepare all functions for the present design
    obj->outer_val = EvalObjective(n, x_outer.GetPointer(), true);

    optimization->GetDesign()->Reset(DesignElement::COST_GRADIENT, DesignElement::DEFAULT);

    EvalGradObjective(n, x_outer.GetPointer(), true, obj->outer_grad);
    LOG_DBG3(sgp) << "FP:UTCP obj=" << obj->ToString() << " outer_val=" << obj->outer_val << " grad=" << obj->outer_grad.ToString();
  }

  // reset values of the constraint gradients before the loop
  // as it also contains a loop over all the design elements
  optimization->GetDesign()->Reset(DesignElement::CONSTRAINT_GRADIENT, DesignElement::DEFAULT);
  for(unsigned int i = 0; i < m; i++)
  {
    Condition* g = constr[i];
    constr[i]->SetValue(optimization->CalcConstraint(g));
    if (constr[i]->GetType() == Condition::FILTERING_GAP) {
      EvalGradConstraint(g, 0, true, true, filter_outer_grad);
      for (unsigned int j=0; j < obj->outer_grad.GetSize(); j++)
        obj->outer_grad[j] += ppen_filt * filter_outer_grad[j];
    }

    LOG_DBG3(sgp) << "FP:UTCP g[" << i << "]=" << constr[i]->ToString() << " outer_val=" << constr[i]->GetValue(); //<< " grad=" << constr[i]->outer_grad.ToString();
  }
  optimization->constraints.view->Done(); // reset slope constraint to global mode
}

void SGP::OuterToDesign(bool filter) {
  DesignSpace* space = optimization->GetDesign();
  int e_num;
  for (unsigned int i = 0; i < n; i++) {
      e_num = space->data[i].elem->elemNum - 1.;
      switch (space->data[i].GetType()) {
        case (DesignElement::DENSITY):
          if (!filter)
            space->data[i].SetDesign(rho_outer[e_num]);
          break;
        case DesignElement::ROTANGLE:
          if (!filter)
            space->data[i].SetDesign(theta_outer[e_num]);
          break;
        case DesignElement::STIFF1:
          if (!filter)
            space->data[i].SetDesign(s1_outer[e_num]);
          break;
        case DesignElement::STIFF2:
          if (!filter)
            space->data[i].SetDesign(s2_outer[e_num]);
          break;
        case DesignElement::MECH_11:
          space->data[i].SetDesign(E_outer[e_num][0][0]);
          break;
        case DesignElement::MECH_12:
          space->data[i].SetDesign(E_outer[e_num][0][1]);
          break;
        case DesignElement::MECH_13:
          space->data[i].SetDesign(E_outer[e_num][0][2]);
          break;
        case DesignElement::MECH_22:
          space->data[i].SetDesign(E_outer[e_num][1][1]);
          break;
        case DesignElement::MECH_23:
          space->data[i].SetDesign(E_outer[e_num][1][2]);
          break;
        case DesignElement::MECH_33:
          space->data[i].SetDesign(E_outer[e_num][2][2]);
          break;
        default:
          break;
      }
  }

  // enforce design update and solution of linear system
  if (!filter)
    space->IncrementDesignId();
}

void SGP::GetOuterDerivativeVector(StdVector<double> & out, StdVector<Double> obj_grad) {
  DesignSpace* space = optimization->GetDesign();
  unsigned int mech11 = space->FindDesign(DesignElement::MECH_11);
  unsigned int mech12 = space->FindDesign(DesignElement::MECH_12);
  unsigned int mech13 = space->FindDesign(DesignElement::MECH_13);
  unsigned int mech22 = space->FindDesign(DesignElement::MECH_22);
  unsigned int mech23 = space->FindDesign(DesignElement::MECH_23);
  unsigned int mech33 = space->FindDesign(DesignElement::MECH_33);

  for (unsigned int i = 0; i < 6 * n_elem; i += 6) {
    out[i] = obj_grad[mech11*n_elem+i/6];
    out[i+1] = obj_grad[mech12*n_elem+i/6];
    out[i+2] = obj_grad[mech13*n_elem+i/6];
    out[i+3] = obj_grad[mech22*n_elem+i/6];
    out[i+4] = obj_grad[mech23*n_elem+i/6];
    out[i+5] = obj_grad[mech33*n_elem+i/6];
  }
}
void SGP::GetOuterDerivative(StdVector<Matrix<double> > & out, StdVector<Double> obj_grad) {
  DesignSpace* space = optimization->GetDesign();
  unsigned int mech11 = space->FindDesign(DesignElement::MECH_11);
  unsigned int mech12 = space->FindDesign(DesignElement::MECH_12);
  unsigned int mech13 = space->FindDesign(DesignElement::MECH_13);
  unsigned int mech22 = space->FindDesign(DesignElement::MECH_22);
  unsigned int mech23 = space->FindDesign(DesignElement::MECH_23);
  unsigned int mech33 = space->FindDesign(DesignElement::MECH_33);

  for (unsigned int i = 0; i < n_elem; i++) {
    out[i][0][0] = obj_grad[mech11*n_elem+i];
    out[i][0][1] = obj_grad[mech12*n_elem+i];
    out[i][0][2] = obj_grad[mech13*n_elem+i];
    out[i][1][1] = obj_grad[mech22*n_elem+i];
    out[i][1][2] = obj_grad[mech23*n_elem+i];
    out[i][2][2] = obj_grad[mech33*n_elem+i];
    out[i][1][0] = out[i][0][1];
    out[i][2][0] = out[i][0][2];
    out[i][2][1] = out[i][1][2];
  }
}

void SGP::TestRotation() {
  Matrix<double> unit(3,3);
  unit[0][0] = 1;
  unit[0][1] = 0;
  unit[0][2] = 0;
  unit[1][0] = 0;
  unit[1][1] = 0;
  unit[1][2] = 0;
  unit[2][0] = 0;
  unit[2][1] = 0;
  unit[2][2] = 0;

  Matrix<double> unit2 = unit;
  Matrix<double> unit3 = unit;

  Matrix<double> rot(2,2);

  GetRotationMatrix(0,rot);
  unit.PerformRotation(rot,unit);
  std::cout<<"RotbyZero = "<<unit.ToString(2)<<std::endl;

  GetRotationMatrix(0.78539816339,rot);
  unit.PerformRotation(rot,unit);
  std::cout<<"Rotby0.5 = "<<unit.ToString(2)<<std::endl;
  std::cout<<"Rot = "<<rot.ToString(2)<<std::endl;

  //DesignMaterial::RotateHMStiffnessTensor(unit2,PLANE,DesignElement::NO_DERIVATIVE,0.78539816339,DesignMaterial::HILL_MANDEL);
  //std::cout<<"Rotby0.5 HM = "<<unit2.ToString(2)<<std::endl;

  //DesignMaterial::RotateHMStiffnessTensor(unit3,PLANE,DesignElement::NO_DERIVATIVE,0.78539816339,DesignMaterial::VOIGT);
  //std::cout<<"Rotby0.5 Voigt = "<<unit3.ToString(2)<<std::endl;

  optimization->GetDesign()->designMaterial->RotateTensor(unit2,DesignElement::NO_DERIVATIVE, DesignMaterial::VOIGT, DesignMaterial::CCW, true, 0.78539816339);
  std::cout<<"Rotby0.5 Voigt = "<<unit2.ToString(2)<<std::endl;

}

void SGP::GetRotationMatrix(double alpha, Matrix<Double> & rot) {
  rot[0][0] = cos(alpha);
  rot[0][1] = sin(alpha);
  rot[1][0] = -sin(alpha);
  rot[1][1] = cos(alpha);
}

void SGP::DesignToOuter(bool inner,bool only_update_outer) {
  DesignSpace* space = optimization->GetDesign();

  // 2D rotation matrix
  Matrix<double> Rot;
  Rot.Resize(2,2);
  unsigned int dens = 0., rot = 0., s1 = 0., s2 = 0.;
  Vector<double> p(2);
  //designMaterial->GetTensor(material, dtype, stt, de->elem, de->GetType(), f->GetNotation());
  if (!only_update_outer) {
    if (configuration == DENSITY_ROTANGLE || configuration == FOMO || configuration == FMO) {
      dens = space->FindDesign(DesignElement::DENSITY);
      rot = space->FindDesign(DesignElement::ROTANGLE);
    } else if (configuration == STIFF1_STIFF2) {
      s1 = space->FindDesign(DesignElement::STIFF1);
      s2 = space->FindDesign(DesignElement::STIFF2);
      rot = space->FindDesign(DesignElement::ROTANGLE);

    }
  }
  unsigned int mech11 = space->FindDesign(DesignElement::MECH_11);
  unsigned int mech12 = space->FindDesign(DesignElement::MECH_12);
  unsigned int mech13 = space->FindDesign(DesignElement::MECH_13);
  unsigned int mech22 = space->FindDesign(DesignElement::MECH_22);
  unsigned int mech23 = space->FindDesign(DesignElement::MECH_23);
  unsigned int mech33 = space->FindDesign(DesignElement::MECH_33);


  for (unsigned int i = 0; i < n_elem; i++) {
    if (configuration == DENSITY_ROTANGLE || configuration == FOMO || configuration == FMO) {
      if (!only_update_outer) {
        rho_outer[i] = space->GetDesignElement(dens*n_elem+i)->GetDesign(DesignElement::SMART);
        theta_outer[i] = space->GetDesignElement(rot*n_elem+i)->GetDesign(DesignElement::SMART);
      }
      if (!inner) {
        E_outer[i].Resize(3,3);
        Matrix<double> tmp(E_0);
        if (configuration == FOMO || configuration == FMO) {
          tmp = E_inner[i];
        }
        space->designMaterial->RotateTensor(tmp,DesignElement::NO_DERIVATIVE, DesignMaterial::HILL_MANDEL, DesignMaterial::CCW, true, theta_outer[i]);
        //LOG_DBG3(sgp) << "A:before rho E_outer["<< i << "]= "<< E_outer[i].ToString();
        E_outer[i] = tmp;
        E_outer[i] *= tf->Transform(rho_outer[i]);
      }
    } else if (configuration == STIFF1_STIFF2) {
      if (!only_update_outer) {
        s1_outer[i] = space->GetDesignElement(s1*n_elem+i)->GetDesign(DesignElement::SMART);
        s2_outer[i] = space->GetDesignElement(s2*n_elem+i)->GetDesign(DesignElement::SMART);
        theta_outer[i] = space->GetDesignElement(rot*n_elem+i)->GetDesign(DesignElement::SMART);
      }
      p[0] = s1_outer[i];
      p[1] = s2_outer[i];
      if (!inner) {
        helper_dm->ApplyHomRectC1Tensor(E_outer[i],p,DesignElement::NO_DERIVATIVE,PLANE);
        space->designMaterial->RotateTensor(E_outer[i],DesignElement::NO_DERIVATIVE, DesignMaterial::HILL_MANDEL, DesignMaterial::CW, true, theta_outer[i]);
      }
    }
//    x_outer[mech11*n_elem+i] = E_outer[i][0][0];
//    x_outer[mech12*n_elem+i] = E_outer[i][0][1];
//    x_outer[mech13*n_elem+i] = E_outer[i][0][2];
//    x_outer[mech22*n_elem+i] = E_outer[i][1][1];
//    x_outer[mech23*n_elem+i] = E_outer[i][1][2];
//    x_outer[mech33*n_elem+i] = E_outer[i][2][2];
    LOG_DBG3(sgp) << "A:E_outer["<< i << "]= "<< E_outer[i].ToString();
  }
  // write tensor to design in order to apply filter
  if (!inner)
    OuterToDesign(true);
  for (unsigned int i = 0; i < n_elem; i++) {
    if (!inner) {
      E_outer[i][0][0] = space->GetDesignElement(mech11*n_elem+i)->GetDesign(DesignElement::SMART);
      E_outer[i][0][1] = space->GetDesignElement(mech12*n_elem+i)->GetDesign(DesignElement::SMART);
      E_outer[i][0][2] = space->GetDesignElement(mech13*n_elem+i)->GetDesign(DesignElement::SMART);
      E_outer[i][1][1] = space->GetDesignElement(mech22*n_elem+i)->GetDesign(DesignElement::SMART);
      E_outer[i][1][2] = space->GetDesignElement(mech23*n_elem+i)->GetDesign(DesignElement::SMART);
      E_outer[i][2][2] = space->GetDesignElement(mech33*n_elem+i)->GetDesign(DesignElement::SMART);
      E_outer[i][1][0] = E_outer[i][0][1];
      E_outer[i][2][0] = E_outer[i][0][2];
      E_outer[i][2][1] = E_outer[i][1][2];
    }
    x_outer[mech11*n_elem+i] = E_outer[i][0][0];
    x_outer[mech12*n_elem+i] = E_outer[i][0][1];
    x_outer[mech13*n_elem+i] = E_outer[i][0][2];
    x_outer[mech22*n_elem+i] = E_outer[i][1][1];
    x_outer[mech23*n_elem+i] = E_outer[i][1][2];
    x_outer[mech33*n_elem+i] = E_outer[i][2][2];
  }
}

void SGP::ToInfoIter(PtrParamNode pn_iter)
{
  pn_iter->Get("ppen_vol")->SetValue(ppen_vol);
  pn_iter->Get("ppeni")->SetValue(ppeni);
  pn_iter->Get("merit")->SetValue(merit);
}

void SGP::ToInfo()
{
  PtrParamNode sgp = optimization->optInfoNode->Get(ParamNode::HEADER)->Get("sgp");

  for (unsigned int i=0; i < inner_variables.GetSize();i++) {
    PtrParamNode pn = sgp->Get("innerVariables")->Get("var",ParamNode::APPEND);
    InnerVariable& iv = inner_variables[i];
    pn->Get("type")->SetValue(DesignElement::type.ToString(iv.type));
    pn->Get("upper")->SetValue(iv.upper_bound);
    pn->Get("lower")->SetValue(iv.lower_bound);
    pn->Get("steps")->SetValue(iv.steps);

  }
}

void SGP::SetConfiguration(){
  assert(inner_variables.GetSize() > 0);
  configuration = NO_CONF;
  if (HasInnerVar(DesignElement::DENSITY) && HasInnerVar(DesignElement::ROTANGLE) && HasInnerVar(DesignElement::MECH_11) && HasInnerVar(DesignElement::MECH_12) && HasInnerVar(DesignElement::MECH_22) && HasInnerVar(DesignElement::MECH_33)) {
    configuration = FMO;
  }
  else if (HasInnerVar(DesignElement::DENSITY) && HasInnerVar(DesignElement::ROTANGLE) && HasInnerVar(DesignElement::MECH_11) && HasInnerVar(DesignElement::MECH_12) && HasInnerVar(DesignElement::MECH_22)) {
    configuration = FOMO;
  }
  else if (HasInnerVar(DesignElement::DENSITY) && HasInnerVar(DesignElement::ROTANGLE)) {
    configuration = DENSITY_ROTANGLE;
  }
  else if (HasInnerVar(DesignElement::STIFF1) && HasInnerVar(DesignElement::STIFF2) && HasInnerVar(DesignElement::ROTANGLE)) {
    configuration = STIFF1_STIFF2;
  }
  if (configuration == NO_CONF) {
    throw Exception("Unknown configuration of inner variables.");
  }
}

SGP::InnerVariable& SGP::GetInnerVar(DesignElement::Type tp) {
  for (unsigned int i = 0; i < inner_variables.GetSize();i++) {
    if (inner_variables[i].type == tp) {
      return inner_variables[i];
    }
  }
  throw Exception("Undefined inner variable type.");
}


SGPApproximation::SGPApproximation(SGP* sgp)
{
  this->common = sgp;
  this->outer_val = -1.0;
}

void SGPApproximation::PostInit()
{
  Function*  f = GetFunction();
  jac_pattern = f->GetSparsityPattern();
  outer_grad.Resize(jac_pattern.GetSize());
}

double SGPApproximation::SubSolve(Eval eval, StdVector<Matrix<double> > df, double ppen, StdVector<double> & s1_outer, StdVector<double> & s2_outer, StdVector<double> & theta_outer) {
  double obj = 0.0;

  Matrix<double> dL,BB,E_tmptmp(3,3), E_tmp(3,3);
  double obj_min;

  SGP::InnerVariable& s1_iv = common->GetInnerVar(DesignElement::STIFF1);
  SGP::InnerVariable& s2_iv = common->GetInnerVar(DesignElement::STIFF2);
  SGP::InnerVariable& theta_iv = common->GetInnerVar(DesignElement::ROTANGLE);

  // Loop over all elements
  for (unsigned int i = 0; i < common->n_elem; i++) {
    dL = common->E_outer[i] - common->L[i];
    for (int ii = 0;ii < 3; ii++) {
      for (int jj = 0;jj < 3;jj++) {
        if (ii != jj) {
          df[i][ii][jj] /= 2.;
        }
      }
    }
    BB = -dL * df[i] * dL;
    Matrix<double> tmp(BB);
    LOG_DBG3(sgp) << "Subsolve: dL =" << dL.ToString();
    LOG_DBG3(sgp) << "Subsolve: BB =" << BB.ToString();
    obj_min = std::numeric_limits<double>::infinity();
    //brute force optimization
    for (double theta = theta_iv.lower_bound; theta <= theta_iv.upper_bound; theta += theta_iv.inc) {
      for (double s1 = s1_iv.lower_bound; s1 <= s1_iv.upper_bound; s1 += s1_iv.inc) {
        for (double s2 = s2_iv.lower_bound; s2 <= s2_iv.upper_bound; s2 += s2_iv.inc) {
          CalcE_inner(E_tmptmp,s1,s2,theta,tmp);
          LOG_DBG3(sgp) << "Subsolve BBphi = ["<<tmp.ToString()<<"]";
          E_tmp = E_tmptmp;
          obj = EvalApproximation(s1+s2-s1*s2,eval, BB, E_tmp, ppen,i);
          LOG_DBG3(sgp) << "Subsolve: s1 =" << s1 << " s2 = " << s2 << " theta= " << theta <<" E_tmp = [" << E_tmp.ToString() << "]";
          LOG_DBG3(sgp) << "Subsolve: obj =" << obj;
          if (obj < obj_min) {
            s1_outer[i] = s1;
            s2_outer[i] = s2;
            theta_outer[i] = theta;
            obj_min = obj;
          }
        }
      }
    }
    LOG_DBG3(sgp) << "Subsolve: s1_min =" << s1_outer[i] << " s2_min = " << s2_outer[i] << " theta_min= " << theta_outer[i] << "]";
    LOG_DBG3(sgp) << "Subsolve: obj_min =" << obj_min;
    // int dummy = 1;
  }
  return obj;
}

double SGPApproximation::SubSolve_Density_Rotangle(Eval eval, StdVector<Matrix<double> > df, double ppen, StdVector<double> & rho_outer, StdVector<double> & theta_outer) {
  double obj = 0.0;

  Matrix<double> dL,BB,E_tmptmp(3,3), E_tmp(3,3);
  double obj_min;

  SGP::InnerVariable& theta_iv = common->GetInnerVar(DesignElement::ROTANGLE);
  SGP::InnerVariable& rho_iv = common->GetInnerVar(DesignElement::DENSITY);

  // Loop over all elements
  for (unsigned int i = 0; i < common->n_elem; i++) {
    dL = common->E_outer[i] - common->L[i];
    /*for (int ii = 0;ii < 3; ii++) {
      for (int jj = 0;jj < 3;jj++) {
        if (ii != jj) {
          df[i][ii][jj] /= 2.;
        }
      }
    }*/
    BB = -dL * df[i] * dL;
    LOG_DBG3(sgp) << "Subsolve: dL =" << dL.ToString();
    LOG_DBG3(sgp) << "Subsolve: BB =" << BB.ToString();
    obj_min = std::numeric_limits<double>::infinity();
    //brute force optimization
    for (double theta = theta_iv.lower_bound; theta <= theta_iv.upper_bound; theta += theta_iv.inc) {
      CalcE_inner_Density_Rotangle(E_tmptmp, common->E_0, theta);
      for (double rho = rho_iv.lower_bound; rho <= rho_iv.upper_bound; rho += rho_iv.inc) {
        E_tmp = E_tmptmp;
        E_tmp *= common->tf->Transform(rho);
        obj = EvalApproximation(rho,eval, BB, E_tmp, ppen, i);
        LOG_DBG3(sgp) << "Subsolve: rho =" << rho << " theta = " << theta << " E_tmp = [" << E_tmp.ToString() << "]";
        LOG_DBG3(sgp) << "Subsolve: obj =" << obj;
        if (obj < obj_min) {
          rho_outer[i] = rho;
          theta_outer[i] = theta;
          obj_min = obj;
        }
      }
    }
  }
  return obj;
}

double SGPApproximation::SubSolve_FOMO(Eval eval, StdVector<Matrix<double> > df, double ppen, StdVector<double> & rho_outer, StdVector<double> & theta_outer) {
  double obj = 0.0;

  Matrix<double> dL,BB,E_tmptmp(3,3);
  double obj_min;

  SGP::InnerVariable& theta_iv = common->GetInnerVar(DesignElement::ROTANGLE);

  double rho1 = 0, rho2 = 0, rho,rho1_min,rho2_min;
  Vector<double> ev(2);
  Matrix<double> ev_vector(2,2),ev_vectorT(2,2), ev_vector_min(2,2),ev_vectorT_min(2,2);
  Matrix<double> atmp(2,2),help(2,2);
  double l1_min,l2_min;
  // Loop over all elements
  for (unsigned int i = 0; i < common->n_elem; i++) {
    dL = common->E_outer[i] - common->L[i];

    // TODO: Why is this transformation necessary?
    for (int ii = 0;ii < 3; ii++) {
      for (int jj = 0;jj < 3;jj++) {
        if (ii != jj) {
          df[i][ii][jj] /= 2.;
        }
      }
    }
    BB = -dL * df[i] * dL;

    LOG_DBG3(sgp) << "Subsolve: dL =" << dL.ToString();
    LOG_DBG3(sgp) << "Subsolve: df =" << df[i].ToString();
    LOG_DBG3(sgp) << "Subsolve: BB =" << BB.ToString();
    obj_min = std::numeric_limits<double>::infinity();
    //brute force optimization
    for (double theta = theta_iv.lower_bound; theta <= theta_iv.upper_bound; theta += theta_iv.inc) {
      obj = CalcAnalyticSol_FOMO(rho1, rho2, rho, ev,  ev_vector,  eval, BB, theta, ppen, i);
      ev_vector.Transpose(ev_vectorT);
      LOG_DBG3(sgp) << "Subsolve: theta =" << theta << " obj = " << obj <<" rho = " <<rho<< " rho1 = " << rho1 << " rho2 = " << rho2<< " ppen = " << ppen;

      if (obj < obj_min) {
        atmp[0][0] = rho1;
        atmp[1][1] = rho2;
        atmp[0][1] = 0.;
        atmp[1][0] = 0.;

        //atmp = ev_vector*tmp2*ev_vectorT
        atmp.Mult(ev_vectorT,help);
        ev_vector.Mult(help,atmp);

        rho_outer[i] = rho;
        theta_outer[i] = theta;
        ev_vector_min = ev_vector;
        ev_vectorT_min = ev_vectorT;
        l1_min = ev[0];
        l2_min = ev[1];
        rho1_min = rho1;
        rho2_min = rho2;
        common->E_inner[i][0][0] = atmp[0][0];
        common->E_inner[i][0][1] = atmp[0][1];
        common->E_inner[i][0][2] = 0.;
        common->E_inner[i][1][0] = atmp[1][0];
        common->E_inner[i][1][1] = atmp[1][1];
        common->E_inner[i][1][2] = 0.;
        common->E_inner[i][2][2] = 1.-atmp[0][0]-atmp[1][1];
        common->E_inner[i][2][0] = 0.;
        common->E_inner[i][2][1] = 0.;
        obj_min = obj;
      }
    }
    LOG_DBG3(sgp)<< "Subsolve: theta_min = " << theta_outer[i] << ", obj_min = "<<obj_min<< ", rho_min= "<<rho_outer[i];
    LOG_DBG3(sgp) << "Subsolve: E_inner =" << common->E_inner[i].ToString();
    LOG_DBG3(sgp) << "Subsolve: ev_vector_min =" << ev_vector_min.ToString();
    LOG_DBG3(sgp) << "Subsolve: l1_min =" <<l1_min <<" l2_min = "<<l2_min;
    LOG_DBG3(sgp) << "Subsolve: rho1_min =" <<rho1_min <<" rho2_min = "<<rho2_min;
    LOG_DBG3(sgp) << "Subsolve: atmp =" << atmp.ToString();

  }

  return obj;
}

double SGPApproximation::EvalApproximation(double vol_inner_vars, Eval eval, Matrix<double> BB, Matrix<double> E_tmptmp, double ppen, int index) {

  Matrix<double> E_tmp(3,3), E_tmpinv(3,3),LLL(3,3);
  double result = 0;
  switch(eval)
  {
  case FUNC:
  {
    E_tmp = E_tmptmp-common->L[index];
    E_tmp.Invert(E_tmpinv);
    // use E_tmp again for proximal point term, default: tau = 0
    E_tmp = E_tmptmp - common->E_outer[index];
    E_tmp *= common->tau;
    BB += E_tmp;
    E_tmpinv.Mult(BB,LLL);
    result = LLL[0][0] + LLL[1][1] + LLL[2][2] + ppen * vol_inner_vars;
    LOG_DBG3(sgp) << "A:E f=" << ToString() << " func r_= " << result;
    break;
  }
  case GRAD:
  {
    // diff(fp,xi); diff(fn,xi)
    //double c = grad >= 0 ? 2.0*tau*(xi-xo)/(U-xi) + tau*(xi-xo)*(xi-xo)/((U-xi)*(U-xi))
    //                     : 2.0*tau*(xi-xo)/(xi-L) - tau*(xi-xo)*(xi-xo)/((xi-L)*(xi-L));
    //(*out)[e] = + p/((U-xi)*(U-xi)) - q/((xi-L)*(xi-L)) + c;
    //LOG_DBG3(sgp) << "A:E f=" << ToString() << " c=" << c << " grad out[" << e << "]=" << (*out)[e];
    break;
  }
  default:
    break;
  }

  return result;
}

double SGPApproximation::CalcAnalyticSol_FOMO(double &rho1, double &rho2, double & rho, Vector<double> & ev,  Matrix<double> & ev_vector,  Eval eval, Matrix<double> BB, double theta_inner, double ppen, int index) {

  Matrix<double> E_tmp(2,2), E_tmpinv(3,3),LLL(3,3);
  double result = 0;
  double psimp = common->tf->GetParam();
  switch(eval)
  {
  case FUNC:
  {
    Matrix<double> tmp(BB);
    ev.Resize(2);
    ev_vector.Resize(2,2);
    common->optimization->GetDesign()->designMaterial->RotateTensor(tmp,DesignElement::NO_DERIVATIVE, DesignMaterial::HILL_MANDEL, DesignMaterial::CW, true, theta_inner);
    E_tmp[0][0] = tmp[0][0];
    E_tmp[0][1] = tmp[0][1];
    E_tmp[1][0] = tmp[1][0];
    E_tmp[1][1] = tmp[1][1];
    E_tmp.eigenvaluesWithLapack(ev,&ev_vector);
    double sqrB1 = sqrt(ev[0]);
    double sqrB2 = sqrt(ev[1]);
    double sqrB3 = sqrt(tmp[2][2]);

    rho1 = std::min(std::max(.01, .99*sqrB1/(sqrB1+sqrB2+sqrB3)),.98);
    rho2 = std::min(std::max(.01, .99*sqrB2/(sqrB1+sqrB2+sqrB3)),.98);

    double rtmp = ev[0]/rho1+ev[1]/rho2+tmp[2][2]/(1.0-rho1-rho2);

    rho = pow(psimp * rtmp/ppen,1./(psimp+1.));
    rho = std::min(1.,std::max(.01,rho));
    result =  pow(rho,(-psimp))*rtmp + ppen * rho;
    break;
  }
  default:
    break;
  }

  return result;
}

void SGPApproximation::CalcE_inner(Matrix<double> & E_inner, double s1, double s2, double theta_inner, Matrix<double> & tmp) {
  E_inner.Resize(3,3);
  Vector<double> p(2);
  p[0] = s1;
  p[1] = s2;
  common->helper_dm->ApplyHomRectC1Tensor(E_inner,p,DesignElement::NO_DERIVATIVE,PLANE);
  // Rotate 2-scale tensor
  common->optimization->GetDesign()->designMaterial->RotateTensor(E_inner,DesignElement::NO_DERIVATIVE, DesignMaterial::HILL_MANDEL, DesignMaterial::CW, true, theta_inner);
  // only for debuggin purpose, delete afterwards
  common->optimization->GetDesign()->designMaterial->RotateTensor(tmp,DesignElement::NO_DERIVATIVE, DesignMaterial::HILL_MANDEL, DesignMaterial::CW, true, theta_inner);

}

void SGPApproximation::CalcE_inner_Density_Rotangle(Matrix<double> & E_inner, Matrix<double> E_0, double theta_inner) {
  E_inner.Resize(3,3);
  Matrix<double> tmp(E_0);
  common->optimization->GetDesign()->designMaterial->RotateTensor(tmp,DesignElement::NO_DERIVATIVE, DesignMaterial::VOIGT, DesignMaterial::CCW, true, theta_inner);
  E_inner = tmp;
}

std::string SGPApproximation::ToString()
{
  return Function::type.ToString(GetFunction()->GetType());
}

/*Condition* SGPApproximation::GetCondition(bool determinant)
{
  Function* f = GetFunction(determinant);
  assert(!f->IsObjective());
  return static_cast<Condition*>(f);
}*/

Function* SGPApproximation::GetFunction()
{
    return common->optimization->objectives.data[0];
}








