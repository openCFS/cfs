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

  pmin   = this_opt_pn_->Get("pmin")->As<double>();
  pmax   = this_opt_pn_->Get("pmax")->As<double>();
  lbound = this_opt_pn_->Get("lower_bound")->As<double>();
  ubound = this_opt_pn_->Get("upper_bound")->As<double>();
  tolerance = this_opt_pn_->Get("tolerance")->As<double>();

  ppen = 0;

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
  E_outer.Resize(n_elem);

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

  // setup lower and upper bounds, they might be from design bounds. After this we must not use DesignElement::GetLower/UpperBound() !!
  /** lower_bound.Resize((int) DesignElement::ALL_DESIGNS, 0.0);
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
  }*/

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
  }
  cc.view->Done();
  optimization->constraints.view->Done();
}

void SGP::SolveProblem()
{
  int max_iter = optimization->GetMaxIterations();
  PtrParamNode summary = optimization->optInfoNode->Get(ParamNode::SUMMARY);

  bool converged = false;
  int iter = 1;

  // initialize merit functions and constraints
  double merit_old = -1.;

  // set penalty parameter
  ppen = 0.5*(pmin+pmax);

  // initialize outer gradient
  StdVector<Matrix<double> > df(n_elem);
  for (unsigned int i = 0; i < n_elem;i++) {
    df[i].Resize(3,3);
  }

  // writes design to outer variables rho_outer,... and evaluates objective and constraint functions, gradients
  UpdateToCurrentStep();
  optimization->CommitIteration();

  // outer optimization loop
  while(!optimization->DoStopOptimization() && iter <= max_iter && !converged)
  {
    // solve sub-problem
    PtrParamNode in = info_->Get(ParamNode::PROCESS)->Get("sgp")->Get("subsolver", ParamNode::APPEND);
    in->Get("iter")->SetValue(iter);

    // saves old merit function value
    merit_old = merit;
    merit = obj->outer_val;

    // calculate objective function including constraint penalty term
    for(unsigned int i = 0; i < m; i++)
    {
      // Add non-physical constraint to merit functions
      if (!constr[i]->IsPhysical()) {
        merit += ppen * constr[i]->GetValue();
      }
    }
    LOG_DBG2(sgp) << "SP: it=" << iter << " merit = " << merit;

    // Get outer gradient for subproblem
    GetOuterDerivative(df, obj->outer_grad);

    // solves sub-model and updates rho_outer and theta_outer
    obj->SubSolve(SGPApproximation::FUNC,df,ppen,rho_outer,theta_outer);

    LOG_DBG2(sgp) << "SP: it=" << iter << " c_grad = [" << obj->outer_grad.ToString() << "]";
    LOG_DBG2(sgp) << "SP: it=" << iter << " x_old = [" << x_outer.ToString() << "]";

    bool penal_vol = true;
    double cond = 0;
    int count = 0;
    if(penal_vol) {
      for(unsigned int i = 0; i < m; i++)
      {
        if (constr[i]->IsPhysical()) {
          cond = constr[i]->GetValue() - constr[i]->GetBound();
          count++;
        }
      }
      // only one physical constraint should be used (currently: physical volume)
      assert(count <= 1);
      if (iter > 10 && abs(cond) < 0.02) {
        penal_vol = false;
      }
      if (cond > 0) {
        pmin = ppen;
        ppen = 0.5 * (pmin + ppen);
      } else {
        pmax = ppen;
        ppen = 0.5 * (pmin + ppen);
      }
    }



    //convergence criterion
    if (abs(merit-merit_old) < tolerance) {
      converged = true;
    }

    // Writes new outer variables to Design
    OuterToDesign();
    // evaluate all functions such that we have the function values for the current design
    // the subproblem is based on the approximated values only
    // writes design to outer variables rho_outer,... and evaluates objective and constraint functions, gradients
    UpdateToCurrentStep();
    ToInfo();

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

void SGP::UpdateToCurrentStep()
{
  DesignSpace* space = optimization->GetDesign();
  space->WriteDesignToExtern(x_outer);
  DesignToOuter();

  LOG_DBG(sgp) << "UTCP x=" << x_outer.ToString();

  // prepare all functions for the present design
  obj->outer_val = EvalObjective(n, x_outer.GetPointer(), true);

  optimization->GetDesign()->Reset(DesignElement::COST_GRADIENT, DesignElement::DEFAULT);
  EvalGradObjective(n, x_outer.GetPointer(), true, obj->outer_grad);
  LOG_DBG3(sgp) << "FP:UTCP obj=" << obj->ToString() << " outer_val=" << obj->outer_val << " grad=" << obj->outer_grad.ToString();

  // we have to consider transformation between benson vanderbei and determinant constraints. Therefore me must not call EvalConstraints()
  // reset values of the constraint gradients before the loop
  // as it also contains a loop over all the design elements
  //optimization->GetDesign()->Reset(DesignElement::CONSTRAINT_GRADIENT, DesignElement::DEFAULT);
  for(unsigned int i = 0; i < m; i++)
  {
    Condition* g = constr[i];
    constr[i]->SetValue(optimization->CalcConstraint(g));
    //EvalGradConstraint(g, 0, true, true, constr[i]->outer_grad); // scale and normalize

    LOG_DBG3(sgp) << "FP:UTCP g[" << i << "]=" << constr[i]->ToString() << " outer_val=" << constr[i]->GetValue(); //<< " grad=" << constr[i]->outer_grad.ToString();
  }
  optimization->constraints.view->Done(); // reset slope constraint to global mode
}

void SGP::OuterToDesign() {
  DesignSpace* space = optimization->GetDesign();
  int e_num;
  for (unsigned int i = 0; i < n; i++) {
      e_num = space->data[i].elem->elemNum - 1.;
      switch (space->data[i].GetType()) {
        case (DesignElement::DENSITY):
          space->data[i].SetDesign(rho_outer[e_num]);
          break;
        case DesignElement::ROTANGLE:
          space->data[i].SetDesign(theta_outer[e_num]);
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
  space->IncrementDesignId();
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

void SGP::DesignToOuter() {
  DesignSpace* space = optimization->GetDesign();

  // 2D rotation matrix
  Matrix<double> Rot;
  Rot.Resize(2,2);

  //Function* f = GetFunction(true);

  //designMaterial->GetTensor(material, dtype, stt, de->elem, de->GetType(), f->GetNotation());
  unsigned int dens = space->FindDesign(DesignElement::DENSITY);
  unsigned int rot = space->FindDesign(DesignElement::ROTANGLE);
  unsigned int mech11 = space->FindDesign(DesignElement::MECH_11);
  unsigned int mech12 = space->FindDesign(DesignElement::MECH_12);
  unsigned int mech13 = space->FindDesign(DesignElement::MECH_13);
  unsigned int mech22 = space->FindDesign(DesignElement::MECH_22);
  unsigned int mech23 = space->FindDesign(DesignElement::MECH_23);

  unsigned int mech33 = space->FindDesign(DesignElement::MECH_33);


  for (unsigned int i = 0; i < n_elem; i++) {
    rho_outer[i] = space->GetDesignElement(dens*n_elem+i)->GetDesign(DesignElement::SMART);
    theta_outer[i] = space->GetDesignElement(rot*n_elem+i)->GetDesign(DesignElement::SMART);
    E_outer[i].Resize(3,3);
    Matrix<double> tmp(E_0);
    space->designMaterial->RotateTensor(tmp,DesignElement::NO_DERIVATIVE, DesignMaterial::VOIGT, DesignMaterial::CCW, true, theta_outer[i]);
    //LOG_DBG3(sgp) << "A:before rho E_outer["<< i << "]= "<< E_outer[i].ToString();
    E_outer[i] = tmp;
    E_outer[i] *= tf->Transform(rho_outer[i]);
    LOG_DBG3(sgp) << "A:E_outer["<< i << "]= "<< E_outer[i].ToString();
    x_outer[mech11*n_elem+i] = E_outer[i][0][0];
    x_outer[mech12*n_elem+i] = E_outer[i][0][1];
    x_outer[mech13*n_elem+i] = E_outer[i][0][2];
    x_outer[mech22*n_elem+i] = E_outer[i][1][1];
    x_outer[mech23*n_elem+i] = E_outer[i][1][2];
    x_outer[mech33*n_elem+i] = E_outer[i][2][2];
    /** E_outer[i][0][0] = space->GetDesignElement(mech11+i);
    E_outer[i][0][1] = space->GetDesignElement(mech12+i);
    E_outer[i][0][2] = space->GetDesignElement(mech13+i);
    E_outer[i][1][1] = space->GetDesignElement(mech22+i);
    E_outer[i][1][2] = space->GetDesignElement(mech23+i);
    E_outer[i][2][2] = space->GetDesignElement(mech33+i);
    E_outer[i][1][0] = E_outer[i][0][1];
    E_outer[i][2][0] =  E_outer[i][0][2];
    E_outer[i][2][1] = E_outer[i][1][2]; */
  }
}

void SGP::ToInfo()
{
  PtrParamNode iter = optimization->CommitIteration();
  iter->Get("ppen")->SetValue(ppen);
  iter->Get("merit")->SetValue(merit);
}

SGPApproximation::SGPApproximation(SGP* sgp)
{
  this->common = sgp;
  this->outer_val = -1.0;
}

void SGPApproximation::PostInit()
{
  Function*  f = GetFunction(true);
  jac_pattern = f->GetSparsityPattern();
  outer_grad.Resize(jac_pattern.GetSize());
}

double SGPApproximation::SubSolve(Eval eval, StdVector<Matrix<double> > df, double ppen, StdVector<double> & rho_outer, StdVector<double> & theta_outer) {
  double obj = 0.0;

  Matrix<double> dL,BB,E_tmptmp(3,3), E_tmp(3,3);
  double obj_min;
  // Loop over all elements
  for (unsigned int i = 0; i < common->n_elem; i++) {
    dL = common->E_outer[i] - common->L[i];
    BB = -dL * df[i] * dL;
    LOG_DBG3(sgp) << "Subsolve: dL =" << dL.ToString();
    LOG_DBG3(sgp) << "Subsolve: BB =" << BB.ToString();
    obj_min = 10000000.;
    //brute force optimization
    for (double theta = M_PI/4.; theta <= M_PI; theta += M_PI/4.) {
      CalcE_inner(E_tmptmp, common->E_0, theta);
      for (double rho = 0.1; rho <= 1.; rho += 0.1) {
        E_tmp = E_tmptmp;
        E_tmp *= common->tf->Transform(rho);
        obj = EvalApproximation(rho, theta,eval, common->E_0, BB, E_tmp, common->L[i], ppen);
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

double SGPApproximation::EvalApproximation(double rho_inner, double theta_inner, Eval eval, Matrix<double> E_0, Matrix<double> BB, Matrix<double> E_tmptmp, Matrix<double> L0, double ppen) {

  Matrix<double> E_tmp(3,3), E_tmpinv(3,3),LLL(3,3);
  double result = 0;
  switch(eval)
  {
  case FUNC:
  {
    E_tmp = E_tmptmp-L0;
    E_tmp.Invert(E_tmpinv);
    E_tmpinv.Mult(BB,LLL);
    result = LLL[0][0] + LLL[1][1] + LLL[2][2] + ppen * rho_inner;
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

void SGPApproximation::CalcE_inner(Matrix<double> & E_inner, Matrix<double> E_0, double theta_inner) {
  Matrix<double> Rot;
  Rot.Resize(2,2);

  E_inner.Resize(3,3);
  Matrix<double> tmp(E_0);
  common->optimization->GetDesign()->designMaterial->RotateTensor(tmp,DesignElement::NO_DERIVATIVE, DesignMaterial::VOIGT, DesignMaterial::CCW, true, theta_inner);
  E_inner = tmp;
}

std::string SGPApproximation::ToString(bool determinant)
{
  return Function::type.ToString(GetFunction(determinant)->GetType());
}

Condition* SGPApproximation::GetCondition(bool determinant)
{
  Function* f = GetFunction(determinant);
  assert(!f->IsObjective());
  return static_cast<Condition*>(f);
}

Function* SGPApproximation::GetFunction(bool determinant)
{
  assert(constraint_idx >= -1);
  assert(constraint_idx > 0 || common->optimization->objectives.data.GetSize() == 1);

  if(constraint_idx == -1)
    return common->optimization->objectives.data[0];
  else
    return common->optimization->constraints.view->Get(constraint_idx);
}








