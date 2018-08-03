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
using boost::lexical_cast;

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
  volume_tolerance = this_opt_pn_->Get("volume_tolerance")->As<double>();
  upper_tau = this_opt_pn_->Get("upper_tau")->As<double>();
  derivative_check = this_opt_pn_->Get("derivative_check")->As<bool>();


  ppen_vol = 0.;
  ppen_filt = 0.;
  ppeni = 0.;
  pmaxi = 0.;
  pmini = 0.;
  bisect = this_opt_pn_->Get("bisect_iter")->As<int>();

  merit = 0;
  compliance = -1;
  volume = 0.;
  volume_observe = 0.;
  volume_bound = 0.;
  volume_cfs = 0;
  filtering_gaps = 0.;
  filtering_gaps_observe = 0.;
  filtering_gaps_bound = 0.;

  tf = optimization->GetDesign()->GetTransferFunction(DesignElement::DENSITY, App::MECH);
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

  // set design variable configuration, e.g. DESIGN_ROTANGLE or STIFF1_STIFF2 or FOMO or FOMO_Top or FMO
  SetConfiguration();
  ErsatzMaterial* em = dynamic_cast<ErsatzMaterial*>(optimization);
  assert(em != NULL);
  if (configuration == STIFF1_STIFF2) {
    helper_dm = new DesignMaterial(this_opt_pn_->Get("subsolver/designMaterial"),
            OptimizationMaterial::system.Parse(em->pn->Get("material")->As<std::string>()),
            optimization->GetDesign()->design, optimization->GetDesign());
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
  //m = cc.view->GetNumberOfActiveConstraints();
  m = cc.view->GetNumberOfTotalConstraints();
  constr.Resize(m);

  for(unsigned int i = 0; i < m; i++)
  {
    // this is the only place where we are allowed not to use MMAApproximation::GetCondition()
    Condition* g = cc.view->Get(i);
    //MMAApproximation* func = new MMAApproximation(this, i, false);
    constr[i] = g;
    if (constr[i]->GetType() == Condition::FILTERING_GAP) {
      // setup size of the filter constraint derivative
      filtering_gap_grad.Resize(obj->outer_grad.GetSize());
      filtering_gaps_bound += constr[i]->GetBoundValue();
    }
    if (constr[i]->GetType() == Condition::VOLUME || constr[i]->GetType() == Condition::GLOBAL_TWO_SCALE_VOL || constr[i]->GetType() == Condition::GLOBAL_TENSOR_TRACE)
      volume_bound = constr[i]->GetBoundValue();
      volume_grad.Resize(obj->outer_grad.GetSize());
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
  int iter = 1;
  double Vloc =1.;

  // initialize merit functions and constraints
  double merit_old = -1.;

  // set penalty parameter
  ppen_vol = 0.5*(pmin_vol+pmax_vol);
  ppen_filt = 0.5*(pmin_filt+pmax_filt);

  // initialize outer gradient
  StdVector<Matrix<double> > df(n_elem);
  for (unsigned int i = 0; i < n_elem;i++) {
    df[i].Resize(3,3);
  }

  // writes design to outer variables, e.g. rho_outer
  DesignToOuter(false,true);
  // updates outer variables e.g. rho_outer, writes outer_variables to design, and evaluates objective and constraint functions and gradients
  UpdateToCurrentStep();

  // checks derivatives, if option is turned, on by finite differences
  if (derivative_check) {
    StdVector<double> deriv_check(6*n_elem), deriv(6*n_elem);
    double max_error = 0.;
    deriv_check = GradientCheck(max_error);
    std::cout<<"Derivative check:"<<std::endl;
    // Print finite difference derivative to screen
    std::cout<<"approx derivative = "<<deriv_check.ToString()<<std::endl;
    GetOuterDerivativeVector(deriv,obj->outer_grad);
    std::cout<<"derivative = "<<deriv<<std::endl;
    std::cout<<"max derivative error = "<<max_error<<std::endl;
  }
  optimization->CommitIteration();

  // outer optimization loop
  while(!optimization->DoStopOptimization() && iter <= max_iter && !converged)
  {
    LOG_DBG2(sgp) << "SP: it=" << iter << " merit = " << merit;

    // test if outer iteration progresses
    if (merit <= merit_old) {
      worseCounter = 0;
    } else {
      worseCounter++;
      if (worseCounter > 5) {
        std::cout << "\nSGP stopped due to no progress: " << summary->Get("reason/msg")->As<string>() << std::endl;
        break;
      }
    }
    std::cout << "\n merit = " << merit << " merit_old = " << merit_old << " worse_counter = " << worseCounter << std::endl;


    // Get outer gradient for subproblem
    GetOuterDerivative(df, obj->outer_grad);

    double sub_min = 0;
    int tau = 0;
    bool tau_cond = false;
    Vector<double> l_min(n_elem);

    // Reset L asymptote to zero matrix
    Reset_L();
    // globalization loop, varies asymptote L
    while ( tau < upper_tau && !tau_cond) {
      if (abs(pmin_vol-pmax_vol) < tolerance) {
        // solves sub-model and updates rho_outer and theta_outer
        if (configuration == DENSITY_ROTANGLE) {
          sub_min = obj->SubSolve_Density_Rotangle(SGPApproximation::FUNC,df,ppen_vol,rho_outer,theta_outer);
        } else if (configuration == STIFF1_STIFF2) {
          sub_min = obj->SubSolve(SGPApproximation::FUNC,df,ppen_vol,s1_outer,s2_outer,theta_outer);
        } else if (configuration == FOMO_TOP) {
          sub_min = obj->SubSolve_FOMO_Top(SGPApproximation::FUNC,df,ppen_vol,rho_outer,theta_outer,Vloc);
        } else if (configuration == FOMO) {
          sub_min = obj->SubSolve_FOMO(SGPApproximation::FUNC,df,ppen_vol,theta_outer,Vloc);
        } else if (configuration == FMO) {
          sub_min = obj->SubSolve_FMO(SGPApproximation::FUNC,df,ppen_vol,Vloc);
        } else {
          throw "SGP configuration not known!";
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
        std::string output_str = "";

        //bisection for volume constraint
        while (ki < bisect && penal_vol) {
          output_str = "s-iteration: " + lexical_cast<string>(ki) + " compliance: " + lexical_cast<string>(compliance) + " merit: " + lexical_cast<string>(merit) + " pmaxi: " + lexical_cast<string>(pmaxi) + " pmini: " + lexical_cast<string>(pmini);
          if (configuration == DENSITY_ROTANGLE) {
            sub_min = obj->SubSolve_Density_Rotangle(SGPApproximation::FUNC,df,ppeni,rho_outer,theta_outer);
          } else if (configuration == STIFF1_STIFF2) {
            sub_min = obj->SubSolve(SGPApproximation::FUNC,df,ppeni,s1_outer,s2_outer,theta_outer);
          } else if (configuration == FOMO_TOP) {
            sub_min = obj->SubSolve_FOMO_Top(SGPApproximation::FUNC,df,ppeni,rho_outer,theta_outer,Vloc);
          } else if (configuration == FOMO) {
            sub_min = obj->SubSolve_FOMO(SGPApproximation::FUNC,df,ppeni,theta_outer,Vloc);
          } else if (configuration == FMO) {
            sub_min = obj->SubSolve_FMO(SGPApproximation::FUNC,df,ppeni,Vloc);
          } else {
            throw "SGP configuration not known!";
          }

          // evaluates volume constraint for the current design (FOMO + Top, Density_Rotangle, 2sc)
          // merit function is updated
          // writes design to outer variables rho_outer, but does not write E_outer for bisection steps
          UpdateToCurrentStep(true,ppeni);

          //dirty solution volume fmo + fomo
          if (configuration == FMO || configuration == FOMO) {
            // subtracts old volume from merit (dirty solution)
            merit -= ppeni * volume;
            volume = 0;
            for (unsigned int elem = 0; elem < n_elem; elem++) {
              volume += E_inner[elem][0][0] + E_inner[elem][1][1] + E_inner[elem][2][2];
            }
            volume *= 1./float(n_elem);
            merit += ppeni * volume;
          }

          // Create output for bisection iterations
          if (volume > 0) {
            output_str += " volume: " + lexical_cast<string>(volume/Vloc) + " ppeni: " + lexical_cast<string>(ppeni);//+ " volume_cfs: " + lexical_cast<string>(volume_cfs/scale);
          } else if (volume_observe > 0) {
            output_str += " volume: " + lexical_cast<string>(volume_observe/Vloc);
          }
          if (filtering_gaps > 0 && filtering_gaps_observe == 0) {
            output_str += " filtering_gaps: " + lexical_cast<string>(filtering_gaps) + " pmax_filt: " + lexical_cast<string>(pmax_filt) + " pmin_filt: " + lexical_cast<string>(pmin_filt) +  " ppen_filt: " + lexical_cast<string>(ppen_filt);
          } else if (filtering_gaps > 0 && filtering_gaps_observe > 0) {
            output_str += " filtering_gaps: " + lexical_cast<string>(filtering_gaps) + " pmax_filt: " + lexical_cast<string>(pmax_filt) + " pmin_filt: " + lexical_cast<string>(pmin_filt) +  " ppen_filt: "+ lexical_cast<string>(ppen_filt) + " filtering_gaps_observe = " + lexical_cast<string>(filtering_gaps_observe);
          } else if (filtering_gaps_observe > 0) {
            output_str += " filtering_gaps_observe: "+lexical_cast<string>(filtering_gaps_observe);
          }

          output_str += " subsolve min: " + lexical_cast<string>(sub_min + merit);
          std::cout<<output_str<<std::endl;

          // Update strategy for penalty term of volume constraint
          if(penal_vol && volume > 0) {
            if (abs(volume/(Vloc)-volume_bound) < volume_tolerance) {
              penal_vol = false;
            } else {
              if (volume/(Vloc)-volume_bound > 0) {
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

        // Calculate min eigenvalues per element and update of asymptotes L
        obj->CalcMinEigenvalue(E_outer,l_min);
        Update_L(l_min,tau);

        // globalization criterion
        if (sub_min + compliance + ppen_filt * filtering_gaps > EvalCostFunction()) {
          tau_cond = true;
        } else {
          tau++;
          std::cout << "Widening Step!!!!!"<<std::endl;

        }
     }
     if (filtering_gaps > 0) {
       // Update strategy for filtering gap penalty parameter
       if (abs(filtering_gaps) < filtering_gaps_bound) {
         //penal_filt = true;
       } else {
           if ((filtering_gaps) > filtering_gaps_bound) {
             pmin_filt = ppen_filt;
             ppen_filt = 0.5 * (pmax_filt + ppen_filt);
           } else {
             pmax_filt = ppen_filt;
             ppen_filt = 0.5 * (pmin_filt + ppen_filt);
           }

         //penal_filt = false;
      }
    }
    // saves old merit function value
    merit_old = merit;

    // Writes new outer variables to Design
    //OuterToDesign();
    // evaluate all functions such that we have the function values for the current design
    // the subproblem is based on the approximated values only
    // writes design to outer variables rho_outer,... and evaluates objective and constraint functions, gradients and merit function
    UpdateToCurrentStep();

    PtrParamNode pn_iter = optimization->CommitIteration();
    ToInfoIter(pn_iter);

    LOG_DBG2(sgp) << "SP: it=" << iter << " x_new_curr = [" << x_outer.ToString() << "]";
    iter = optimization->GetCurrentIteration();

    //convergence criterion
    if (abs(merit-merit_old) < tolerance) {
      converged = true;
    }
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

  // calculate finite difference gradient from objective function with central difference quotient
  StdVector<unsigned int> grad_index(6);
  grad_index[0] = space->FindDesign(DesignElement::MECH_11);
  grad_index[1] = space->FindDesign(DesignElement::MECH_12);
  grad_index[2] = space->FindDesign(DesignElement::MECH_13);
  grad_index[3] = space->FindDesign(DesignElement::MECH_22);
  grad_index[4] = space->FindDesign(DesignElement::MECH_23);
  grad_index[5] = space->FindDesign(DesignElement::MECH_33);

  // tolerance for finite difference gradient
  double eps = 1e-6;
  StdVector<double> diff_grad(grad_index.GetSize() * n_elem);
  // Check gradient with central difference quotient
  double obj_right, obj_left, grad_error,count=0;
  StdVector<double> constr_val(m);
  for (unsigned int i = 0; i < n_elem;i++) {
    for (unsigned int j = 0; j < grad_index.GetSize();j++) {
      // calculate central difference quotient
      // Eval function value f(x+h) (right)
      unsigned int index = grad_index[j]*n_elem+i;
      assert(index < x_outer.GetSize());
      // Set to x + eps
      x_outer[index] += eps;
      obj_right = EvalObjective(x_outer, true);
      EvalConstraints(n,x_outer.GetPointer(),m,true,constr_val.GetPointer(),true);
      for (unsigned int k=0; k < m;k++) {
        if (constr[k]->GetType() == Condition::FILTERING_GAP)
        {
          //int ind = k*obj->outer_grad.GetSize() + i;
          obj_right += ppen_filt * constr_val[k];
        } else if ((configuration == FMO || configuration == FOMO) && constr[k]->GetType() == Condition::VOLUME) {
          obj_right += ppen_vol * constr_val[k];
        }
      }
      // Set to x - eps
      x_outer[index] -= 2*eps;

      // Eval function value f(x-h) (left)
      obj_left = EvalObjective(x_outer, true);
      EvalConstraints(n,x_outer.GetPointer(),m,true,constr_val.GetPointer(),true);
      for (unsigned int k=0; k < m;k++) {
        if (constr[k]->GetType() == Condition::FILTERING_GAP)
        {
          obj_left += ppen_filt * constr_val[k];
        } else if ((configuration == FMO || configuration == FOMO) && constr[k]->GetType() == Condition::VOLUME) {
          obj_left += ppen_vol * constr_val[k];
        }
      }
      x_outer[index] += eps;
      diff_grad[count] = (obj_right-obj_left)/ (2.*eps);
      grad_error = abs(diff_grad[count]-obj->outer_grad[index]);
      count++;
      if (grad_error > max_grad_error)
        max_grad_error = grad_error;
    }
  }
  return diff_grad;
}

double SGP::EvalCostFunction(void) {
  // write current solution E_inner to temporary outer variable x
  Vector<double> x(x_outer.GetSize());
  Matrix<double> E_tmp_outer(E_inner[0]);
  Vector<double> p(2);
  DesignSpace* space = optimization->GetDesign();
  unsigned int mech11 = space->FindDesign(DesignElement::MECH_11);
  unsigned int mech12 = space->FindDesign(DesignElement::MECH_12);
  unsigned int mech13 = space->FindDesign(DesignElement::MECH_13);
  unsigned int mech22 = space->FindDesign(DesignElement::MECH_22);
  unsigned int mech23 = space->FindDesign(DesignElement::MECH_23);
  unsigned int mech33 = space->FindDesign(DesignElement::MECH_33);
  for (unsigned int i = 0; i < n_elem; i++) {

    // calculate temporary E_outer
    if (configuration == FOMO || configuration == FOMO_TOP || configuration == DENSITY_ROTANGLE) {
      E_tmp_outer = E_inner[i];
      space->designMaterial->RotateTensor(E_tmp_outer,DesignElement::NO_DERIVATIVE, DesignMaterial::HILL_MANDEL, DesignMaterial::CCW, true, theta_outer[i]);
      if (configuration == DENSITY_ROTANGLE || configuration == FOMO_TOP)
        E_tmp_outer *= tf->Transform(rho_outer[i]);
    } else if ( configuration == FMO) {
      E_tmp_outer = E_inner[i];
    } else if (configuration == STIFF1_STIFF2) {
      p[0] = s1_outer[i];
      p[1] = s2_outer[i];
      helper_dm->ApplyHomRectC1Tensor(E_tmp_outer,p,DesignElement::NO_DERIVATIVE,PLANE);
      space->designMaterial->RotateTensor(E_tmp_outer,DesignElement::NO_DERIVATIVE, DesignMaterial::HILL_MANDEL, DesignMaterial::CW, true, theta_outer[i]);
    } else {
      throw("Configuration not known!");
    }

    // write temporary E_outer to outer variables
    x[mech11*n_elem+i] = E_tmp_outer[0][0];
    x[mech12*n_elem+i] = E_tmp_outer[0][1];
    x[mech13*n_elem+i] = E_tmp_outer[0][2];
    x[mech22*n_elem+i] = E_tmp_outer[1][1];
    x[mech23*n_elem+i] = E_tmp_outer[1][2];
    x[mech33*n_elem+i] = E_tmp_outer[2][2];
  }

  // evaluate objective and constraints for temporary outer variable x
  double obj = EvalObjective(x, true);
  StdVector<double> constr_val(m);
  EvalConstraints(n,x.GetPointer(),m,true,constr_val.GetPointer(),true);
  for (unsigned int k=0; k < m;k++) {
    if (constr[k]->GetType() == Condition::FILTERING_GAP) {
        obj += ppen_filt * constr_val[k];
    } else if (constr[k]->GetType() == Condition::VOLUME) {
        obj += ppen_vol * constr_val[k];
    }
  }
  return obj;
}


void SGP::UpdateToCurrentStep(bool inner, double ppeni, bool widening)
{
  // Update Design for new step
  DesignSpace* space = optimization->GetDesign();
  DesignToOuter(inner,false);
  OuterToDesign();
  space->WriteDesignToExtern(x_outer);

  // Update cost function and cost gradient only for non-inner (non-bisection) iterations
  if (!inner) {
    // prepare all functions for the present design
    compliance = EvalObjective(x_outer, true);

    optimization->GetDesign()->Reset(DesignElement::COST_GRADIENT, DesignElement::DEFAULT);
    EvalGradObjective(x_outer, true, obj->outer_grad);
    LOG_DBG3(sgp) << "FP:UTCP obj=" << obj->ToString() << " compliance=" << compliance << " compliance_grad=" << obj->outer_grad.ToString();
  }

  // Update constraints and constraint gradients
  // reset values of the constraint gradients before the loop
  // as it also contains a loop over all the design elements
  optimization->GetDesign()->Reset(DesignElement::CONSTRAINT_GRADIENT, DesignElement::DEFAULT);
  for(unsigned int i = 0; i < m; i++)
  {
    Condition* g = constr[i];
    if (!inner && constr[i] ->GetType() == Condition::FILTERING_GAP ) {
      constr[i]->SetValue(optimization->CalcConstraint(g));
      EvalGradConstraint(g, 0, true, true, filtering_gap_grad);
      LOG_DBG3(sgp) << "SGP:UTCP outer_grad = "<<obj->outer_grad.ToString();
      if (!constr[i]->IsObservation()) {
        //std::cout<<"n_elem "<<n_elem<<" "<<1./n_elem;
        for (unsigned int j=0; j < obj->outer_grad.GetSize(); j++) {
          obj->outer_grad[j] += (1./n_elem) * ppen_filt * filtering_gap_grad[j];
        }
      }
      LOG_DBG3(sgp) << "SGP:UTCP g[" << i << "]=" << constr[i]->ToString() << "d_type = "<<  constr[i]->GetDesignType() << " constr_value=" << constr[i]->GetValue()<< " grad=" << filtering_gap_grad.ToString() << " observe= "<<constr[i]->IsObservation();

    } else {
      if(constr[i] ->GetType() == Condition::VOLUME || constr[i] ->GetType() == Condition::GLOBAL_TWO_SCALE_VOL || constr[i] ->GetType() == Condition::GLOBAL_TENSOR_TRACE) {
        constr[i]->SetValue(optimization->CalcConstraint(g));
        /**if (configuration == FMO || configuration == FOMO) {
          EvalGradConstraint(g, 0, true, true, volume_grad);
          if (!constr[i]->IsObservation()) {
            for (unsigned int j=0; j < obj->outer_grad.GetSize(); j++) {
              if (ppeni > 0)
                obj->outer_grad[j] += ppeni * volume_grad[j];// * n_elem;
              else
                obj->outer_grad[j] += ppen_vol * volume_grad[j];
            }
            LOG_DBG3(sgp) << "SGP: volume_grad = "<<volume_grad.ToString();
            LOG_DBG3(sgp) << "SGP: ppen_vol = "<<ppen_vol;
          }
        } */
      } else {
        if (!inner)
          throw Exception("Constraint type not handled in SGP Update function!");
      }
    }
    //LOG_DBG3(sgp) << "SGP:UTCP g[" << i << "]=" << constr[i]->ToString() << " constr_value=" << constr[i]->GetValue(); //<< " grad=" << constr[i]->outer_grad.ToString();
  }


  // Reset Merit function
  merit = 0;
  // Reset filtering_gaps for non inner iterations
  if (!inner) {
    filtering_gaps = 0;
    filtering_gaps_observe = 0;
  }
  for(unsigned int i = 0; i < m; i++)
  {
    // Add non-physical constraint to merit functions
    if ((constr[i]->GetType() == Condition::VOLUME || constr[i]->GetType() == Condition::GLOBAL_TWO_SCALE_VOL || constr[i] ->GetType() == Condition::GLOBAL_TENSOR_TRACE)) {
      if (!constr[i]->IsObservation()) {
        volume = constr[i]->GetValue(); //*n_elem;
      } else {
        volume_observe = constr[i]->GetValue();
        volume = 0.;
      }
    // Add filtering constraint value if iteration is not a inner iteration (bisection)
    } else if (!inner && constr[i]->GetType() == Condition::FILTERING_GAP) {
      if (!constr[i]->IsObservation()) {
        //std::cout<<"n_elem "<<n_elem;
        filtering_gaps += (1./n_elem)*  constr[i]->GetValue();
      } else {
        filtering_gaps_observe += (1./n_elem) * constr[i]->GetValue();
      }
    } else {
      if (!inner)
        throw Exception("Constraint type not handled in SGP Update function!");
    }
  }

  // calculate merit function: use cost function value including constraint multiplied by penalty term
  if (ppeni > 0)
    merit = compliance + ppen_filt * filtering_gaps + ppeni * volume;
  else
    merit = compliance + ppen_filt * filtering_gaps + ppen_vol * volume;

  LOG_DBG3(sgp)<<"SGP: UTCP: merit = "<<merit <<" outer_grad = "<<obj->outer_grad.ToString() <<" inner = "<< inner;
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
  int count = 0;
  for (unsigned int i = 0; i < n_elem; i++) {
    out[count] = obj_grad[mech11*n_elem+i];
    out[count+1] = obj_grad[mech12*n_elem+i];
    out[count+2] = obj_grad[mech13*n_elem+i];
    out[count+3] = obj_grad[mech22*n_elem+i];
    out[count+4] = obj_grad[mech23*n_elem+i];
    out[count+5] = obj_grad[mech33*n_elem+i];
    count += 6;
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

void SGP::Reset_L() {
  for(unsigned int i = 0; i < n_elem; i++) {
      L[i].Resize(3,3);
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          if (k==j) {
            L[i][j][k] = 0.;
          }
        }
      }
  }
}

void SGP::Update_L(Vector<double> l_min,int tau) {
  double res = 2;
  for (int i = 2; i < upper_tau-tau; i++) {
    res *= 2;
  }
  if (tau + 1 == upper_tau) {
    res = 1.;
  }
  for(unsigned int i = 0; i < n_elem; i++) {
    L[i].Resize(3,3);
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        if (k==j) {
          L[i][j][k] = l_min[i]/res;
        }
      }
    }
  }
}


void SGP::DesignToOuter(bool inner,bool initial) {
  DesignSpace* space = optimization->GetDesign();
  unsigned int dens = 0., rot = 0., s1 = 0., s2 = 0.;
  Vector<double> p(2);

  // only used for initialization of variables before optimization
  if (initial) {
    if (configuration == DENSITY_ROTANGLE || configuration == FOMO_TOP || configuration == FOMO) {
      rot = space->FindDesign(DesignElement::ROTANGLE);
      if (configuration == DENSITY_ROTANGLE || configuration == FOMO_TOP)
        dens = space->FindDesign(DesignElement::DENSITY);
    } else if (configuration == STIFF1_STIFF2) {
      s1 = space->FindDesign(DesignElement::STIFF1);
      s2 = space->FindDesign(DesignElement::STIFF2);
    }
  }
  unsigned int mech11 = space->FindDesign(DesignElement::MECH_11);
  unsigned int mech12 = space->FindDesign(DesignElement::MECH_12);
  unsigned int mech13 = space->FindDesign(DesignElement::MECH_13);
  unsigned int mech22 = space->FindDesign(DesignElement::MECH_22);
  unsigned int mech23 = space->FindDesign(DesignElement::MECH_23);
  unsigned int mech33 = space->FindDesign(DesignElement::MECH_33);

  // Map inner to outer variables for all outer iterations
  for (unsigned int i = 0; i < n_elem; i++) {
    if (configuration == DENSITY_ROTANGLE || configuration == FOMO_TOP || configuration == FOMO) {
      if (initial) {
        if (configuration == DENSITY_ROTANGLE || configuration == FOMO_TOP) {
          rho_outer[i] = space->GetDesignElement(dens*n_elem+i)->GetDesign(DesignElement::SMART);
        }
        theta_outer[i] = space->GetDesignElement(rot*n_elem+i)->GetDesign(DesignElement::SMART);
      }
      if (!inner) {
        E_outer[i].Resize(3,3);
        // only for initialization
        Matrix<double> tmp(E_0);

        tmp = E_inner[i];
        space->designMaterial->RotateTensor(tmp,DesignElement::NO_DERIVATIVE, DesignMaterial::HILL_MANDEL, DesignMaterial::CCW, true, theta_outer[i]);
        E_outer[i] = tmp;
        if (configuration == DENSITY_ROTANGLE || configuration == FOMO_TOP)
          E_outer[i] *= tf->Transform(rho_outer[i]);
      }
    } else if (configuration == STIFF1_STIFF2) {
      if (initial) {
        s1_outer[i] = space->GetDesignElement(s1*n_elem+i)->GetDesign(DesignElement::SMART);
        s2_outer[i] = space->GetDesignElement(s2*n_elem+i)->GetDesign(DesignElement::SMART);
        theta_outer[i] = space->GetDesignElement(rot*n_elem+i)->GetDesign(DesignElement::SMART);
      }
      if (!inner) {
        p[0] = s1_outer[i];
        p[1] = s2_outer[i];
        helper_dm->ApplyHomRectC1Tensor(E_outer[i],p,DesignElement::NO_DERIVATIVE,PLANE);
        space->designMaterial->RotateTensor(E_outer[i],DesignElement::NO_DERIVATIVE, DesignMaterial::HILL_MANDEL, DesignMaterial::CW, true, theta_outer[i]);
      }
    } else if (configuration == FMO) {
      if (!inner)
        E_outer[i] = E_inner[i];
    }
    LOG_DBG3(sgp) << "A:E_outer["<< i << "]= "<< E_outer[i].ToString();
  }

  // write tensor to design in order to apply filter for outer iterations
  // for inner iterations: write non-changed E_outer to x_outer
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
  if (HasInnerVar(DesignElement::DENSITY) && HasInnerVar(DesignElement::ROTANGLE) && HasInnerVar(DesignElement::MECH_11) && HasInnerVar(DesignElement::MECH_12) && HasInnerVar(DesignElement::MECH_22)) {
    configuration = FOMO_TOP;
  }
  else if (HasInnerVar(DesignElement::DENSITY) && HasInnerVar(DesignElement::ROTANGLE)) {
    configuration = DENSITY_ROTANGLE;
  }
  else if (HasInnerVar(DesignElement::STIFF1) && HasInnerVar(DesignElement::STIFF2) && HasInnerVar(DesignElement::ROTANGLE)) {
    configuration = STIFF1_STIFF2;
  } else if (HasInnerVar(DesignElement::MECH_11) && HasInnerVar(DesignElement::MECH_12) && HasInnerVar(DesignElement::MECH_13) && HasInnerVar(DesignElement::MECH_22) && HasInnerVar(DesignElement::MECH_23) && HasInnerVar(DesignElement::MECH_33)) {
    configuration = FMO;
  } else if (HasInnerVar(DesignElement::ROTANGLE) && HasInnerVar(DesignElement::MECH_11) && HasInnerVar(DesignElement::MECH_12) && HasInnerVar(DesignElement::MECH_22)) {
    configuration = FOMO;
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
}

void SGPApproximation::PostInit()
{
  Function*  f = GetFunction();
  jac_pattern = f->GetSparsityPattern();
  outer_grad.Resize(jac_pattern.GetSize());
}

double SGPApproximation::SubSolve(Eval eval, StdVector<Matrix<double> > df, double ppen, StdVector<double> & s1_outer, StdVector<double> & s2_outer, StdVector<double> & theta_outer) {
  double obj = 0.0;

  Matrix<double> dL,BB,dfdL,E_tmptmp(3,3), E_tmp(3,3);
  double obj_min;

  SGP::InnerVariable& s1_iv = common->GetInnerVar(DesignElement::STIFF1);
  SGP::InnerVariable& s2_iv = common->GetInnerVar(DesignElement::STIFF2);
  SGP::InnerVariable& theta_iv = common->GetInnerVar(DesignElement::ROTANGLE);

  double cost = 0;

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
    dfdL = df[i] * dL;
    BB = -dL * dfdL;
    Matrix<double> tmp(BB);
    LOG_DBG3(sgp) << "Subsolve: dL =" << dL.ToString();
    LOG_DBG3(sgp) << "Subsolve: BB =" << BB.ToString();
    obj_min = std::numeric_limits<double>::infinity();
    //brute force multilevel optimization Level 1
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

    // Level 2
    double s1_iv_lb = ((s1_outer[i] - s1_iv.inc/2.) >= s1_iv.lower_bound) ? s1_outer[i] - s1_iv.inc/2. : s1_iv.lower_bound;
    double s1_iv_ub = ((s1_outer[i] + s1_iv.inc/2.) <= s1_iv.upper_bound) ? s1_outer[i] + s1_iv.inc/2. : s1_iv.upper_bound;
    double s1_inc = (s1_iv_ub - s1_iv_lb)/s1_iv.steps;
    double s2_iv_lb = ((s2_outer[i] - s2_iv.inc/2.) >= s2_iv.lower_bound) ? s2_outer[i] - s2_iv.inc/2. : s2_iv.lower_bound;
    double s2_iv_ub = ((s2_outer[i] + s2_iv.inc/2.) <= s2_iv.upper_bound) ? s2_outer[i] + s2_iv.inc/2. : s2_iv.upper_bound;
    double s2_inc = (s2_iv_ub - s2_iv_lb)/s2_iv.steps;

    for (double theta = theta_iv.lower_bound; theta <= theta_iv.upper_bound; theta += theta_iv.inc/10) {
      for (double s1 = s1_iv_lb; s1 <= s1_iv_ub; s1 += s1_inc) {
        for (double s2 = s2_iv_lb; s2 <= s2_iv_ub; s2 += s2_inc) {
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

    // calculate full model function for globalization
    cost += obj_min + dfdL.Trace();
  }
  return cost;
}

double SGPApproximation::SubSolve_Density_Rotangle(Eval eval, StdVector<Matrix<double> > df, double ppen, StdVector<double> & rho_outer, StdVector<double> & theta_outer) {
  double obj = 0.0;

  Matrix<double> dL,BB,dfdL,E_tmptmp(3,3), E_tmp(3,3);
  double obj_min;

  SGP::InnerVariable& theta_iv = common->GetInnerVar(DesignElement::ROTANGLE);
  SGP::InnerVariable& rho_iv = common->GetInnerVar(DesignElement::DENSITY);

  double cost = 0;
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
    dfdL = df[i] * dL;
    BB = -dL * dfdL;
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
    // calculate full model function for globalization
    cost += obj_min + dfdL.Trace();
  }
  return cost;
}

double SGPApproximation::SubSolve_FOMO_Top(Eval eval, StdVector<Matrix<double> > df, double ppen, StdVector<double> & rho_outer, StdVector<double> & theta_outer, double & Vloc) {
  double obj = 0.0;

  Matrix<double> dL,BB,dfdL,E_tmptmp(3,3);
  double obj_min;

  SGP::InnerVariable& theta_iv = common->GetInnerVar(DesignElement::ROTANGLE);

  double rho1 = 0, rho2 = 0, rho = 0,rho1_min,rho2_min;
  Vector<double> ev(2);
  Matrix<double> ev_vector(2,2),ev_vectorT(2,2), ev_vector_min(2,2),ev_vectorT_min(2,2);
  Matrix<double> atmp(2,2),help(2,2);
  double l1_min,l2_min;
  double cost = 0;
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
    dfdL = df[i] * dL;
    BB = -dL * dfdL;

    LOG_DBG3(sgp) << "Subsolve: dL =" << dL.ToString();
    LOG_DBG3(sgp) << "Subsolve: df =" << df[i].ToString();
    LOG_DBG3(sgp) << "Subsolve: BB =" << BB.ToString();
    obj_min = std::numeric_limits<double>::infinity();
    //brute force optimization
    double Vloc = 0;
    for (double theta = theta_iv.lower_bound; theta <= theta_iv.upper_bound; theta += theta_iv.inc) {
      obj = CalcAnalyticSol_FOMO_Top(rho1, rho2, rho, ev,  ev_vector,  eval, BB, theta, ppen, i,Vloc);
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
        common->E_inner[i][0][0] = atmp[0][0] + common->L[i][0][0];
        common->E_inner[i][0][1] = atmp[0][1] + common->L[i][0][1];
        common->E_inner[i][0][2] = 0. + common->L[i][0][2];
        common->E_inner[i][1][0] = atmp[1][0] + common->L[i][1][0];
        common->E_inner[i][1][1] = atmp[1][1] + common->L[i][1][1];
        common->E_inner[i][1][2] = 0. + common->L[i][1][2];
        common->E_inner[i][2][2] = Vloc-atmp[0][0]-atmp[1][1] + common->L[i][2][2];
        common->E_inner[i][2][0] = 0. + common->L[i][2][0];
        common->E_inner[i][2][1] = 0. + common->L[i][2][1];
        obj_min = obj;
      }
    }
    LOG_DBG3(sgp)<< "Subsolve: theta_min = " << theta_outer[i] << ", obj_min = "<<obj_min<< ", rho_min= "<<rho_outer[i];
    LOG_DBG3(sgp) << "Subsolve: E_inner =" << common->E_inner[i].ToString();
    LOG_DBG3(sgp) << "Subsolve: ev_vector_min =" << ev_vector_min.ToString();
    LOG_DBG3(sgp) << "Subsolve: l1_min =" <<l1_min <<" l2_min = "<<l2_min;
    LOG_DBG3(sgp) << "Subsolve: rho1_min =" <<rho1_min <<" rho2_min = "<<rho2_min;
    LOG_DBG3(sgp) << "Subsolve: atmp =" << atmp.ToString();

    // calculate full model function for globalization
    cost += obj_min + dfdL.Trace();
  }

  return cost;
}

double SGPApproximation::SubSolve_FOMO(Eval eval, StdVector<Matrix<double> > df, double ppen, StdVector<double> & theta_outer, double & Vloc) {
  double obj = 0.0;

  Matrix<double> dL,BB,dfdL;
  double obj_min;

  SGP::InnerVariable& theta_iv = common->GetInnerVar(DesignElement::ROTANGLE);

  double rho1 = 0, rho2 = 0, rho3 = 0, rho1_min,rho2_min, rho3_min;
  Vector<double> ev(3);
  Matrix<double> ev_vector(3,3),ev_vectorT(3,3), ev_vector_min(3,3),ev_vectorT_min(3,3);
  Matrix<double> atmp(3,3),help(3,3);
  double l1_min,l2_min, l3_min;
  double cost = 0;
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
    dfdL = df[i] * dL;
    BB = -dL * dfdL;

    LOG_DBG3(sgp) << "Subsolve: dL =" << dL.ToString();
    LOG_DBG3(sgp) << "Subsolve: df =" << df[i].ToString();
    LOG_DBG3(sgp) << "Subsolve: BB =" << BB.ToString();
    obj_min = std::numeric_limits<double>::infinity();
    //brute force optimization
    for (double theta = theta_iv.lower_bound; theta <= theta_iv.upper_bound; theta += theta_iv.inc) {
      obj = CalcAnalyticSol_FOMO(rho1, rho2, rho3, ev,  ev_vector,  eval, BB,common->L[i], theta, ppen, i, Vloc);
      ev_vector.Transpose(ev_vectorT);
      LOG_DBG3(sgp) << "Subsolve: theta =" << theta << " obj = " << obj << " rho1 = " << rho1 << " rho2 = " << rho2<< " ppen = " << ppen;

      if (obj < obj_min) {
        // E_opt = V*diag(rho_i)*V' + L
        atmp[0][0] = rho1;
        atmp[1][1] = rho2;
        atmp[2][2] = rho3;
        atmp[0][2] = 0.;
        atmp[2][0] = 0.;
        atmp[0][1] = 0.;
        atmp[1][0] = 0.;
        atmp[1][2] = 0.;
        atmp[2][1] = 0.;

        //atmp = ev_vector*tmp2*ev_vectorT
        atmp.Mult(ev_vectorT,help);
        ev_vector.Mult(help,atmp);

        theta_outer[i] = theta;
        ev_vector_min = ev_vector;
        ev_vectorT_min = ev_vectorT;
        l1_min = ev[0];
        l2_min = ev[1];
        l3_min = ev[2];
        rho1_min = rho1;
        rho2_min = rho2;
        rho3_min = rho3;
        common->E_inner[i][0][0] = atmp[0][0] + common->L[i][0][0];;
        common->E_inner[i][0][1] = atmp[0][1] + common->L[i][0][1];;
        common->E_inner[i][0][2] = 0. + common->L[i][0][2];;
        common->E_inner[i][1][0] = atmp[1][0] + common->L[i][1][0];;
        common->E_inner[i][1][1] = atmp[1][1] + common->L[i][1][1];;
        common->E_inner[i][1][2] = 0. + common->L[i][1][2];;
        common->E_inner[i][2][2] = atmp[2][2] + common->L[i][2][2];
        common->E_inner[i][2][0] = 0. + common->L[i][2][0];;
        common->E_inner[i][2][1] = 0. + common->L[i][2][1];;
        obj_min = obj;
      }
    }
    LOG_DBG3(sgp)<< "Subsolve: theta_min = " << theta_outer[i] << ", obj_min = "<<obj_min;
    LOG_DBG3(sgp) << "Subsolve: E_inner =" << common->E_inner[i].ToString();
    LOG_DBG3(sgp) << "Subsolve: ev_vector_min =" << ev_vector_min.ToString();
    LOG_DBG3(sgp) << "Subsolve: l1_min =" <<l1_min <<" l2_min = "<<l2_min <<" l3_min = "<<l3_min;
    LOG_DBG3(sgp) << "Subsolve: rho1_min =" <<rho1_min <<" rho2_min = "<<rho2_min <<" rho3_min = "<<rho3_min;
    LOG_DBG3(sgp) << "Subsolve: atmp =" << atmp.ToString();

    // calculate full model function for globalization
    cost += obj_min + dfdL.Trace();
  }

  return cost;
}

double SGPApproximation::SubSolve_FMO(Eval eval, StdVector<Matrix<double> > df, double ppen, double & Vloc) {
  double obj = 0.0;

  Matrix<double> dL,BB,dfdL;
  double obj_min;

  double rho1 = 0, rho2 = 0,rho3 = 0,rho1_min,rho2_min,rho3_min;
  Vector<double> ev(3);
  Matrix<double> ev_vector(3,3),ev_vectorT(3,3), ev_vector_min(3,3),ev_vectorT_min(3,3);
  Matrix<double> atmp(3,3),help(3,3);
  double l1_min,l2_min,l3_min;
  double cost = 0;
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
    dfdL = df[i] * dL;
    BB = -dL * dfdL;

    LOG_DBG3(sgp) << "Subsolve: dL =" << dL.ToString();
    LOG_DBG3(sgp) << "Subsolve: df =" << df[i].ToString();
    LOG_DBG3(sgp) << "Subsolve: BB =" << BB.ToString();
    obj_min = std::numeric_limits<double>::infinity();
    obj = CalcAnalyticSol_FMO(rho1, rho2, rho3, ev,  ev_vector,  eval, BB, common->L[i],ppen,i,Vloc);
    LOG_DBG3(sgp) << "Subsolve: obj = " << obj << " rho1 = " << rho1 << " rho2 = " <<rho2<< " rho3 = "<<rho3;

    if (obj < obj_min) {
      // E_opt = V*diag(rho_i)*V' + L
      atmp[0][0] = rho1;
      atmp[1][1] = rho2;
      atmp[2][2] = rho3;
      atmp[0][2] = 0.;
      atmp[2][0] = 0.;
      atmp[0][1] = 0.;
      atmp[1][0] = 0.;
      atmp[1][2] = 0.;
      atmp[2][1] = 0.;
      ev_vector.Transpose(ev_vectorT);
      atmp.Mult(ev_vectorT,help);
      ev_vector.Mult(help,atmp);
      ev_vector_min = ev_vector;
      ev_vectorT_min = ev_vectorT;
      l1_min = ev[0];
      l2_min = ev[1];
      l3_min = ev[2];
      rho1_min = rho1;
      rho2_min = rho2;
      rho3_min = rho3;
      common->E_inner[i][0][0] = atmp[0][0] + common->L[i][0][0];
      common->E_inner[i][0][1] = atmp[0][1] + common->L[i][0][1];
      common->E_inner[i][0][2] = atmp[0][2] + common->L[i][0][2];
      common->E_inner[i][1][0] = atmp[1][0] + common->L[i][1][0];
      common->E_inner[i][1][1] = atmp[1][1] + common->L[i][1][1];
      common->E_inner[i][1][2] = atmp[1][2] + common->L[i][1][2];
      common->E_inner[i][2][2] = atmp[2][2] + common->L[i][2][2];
      common->E_inner[i][2][0] = atmp[2][0] + common->L[i][2][0];
      common->E_inner[i][2][1] = atmp[2][1] + common->L[i][2][1];
      obj_min = obj;
    }
    LOG_DBG3(sgp)<< "Subsolve: obj_min = "<<obj_min;
    LOG_DBG3(sgp) << "Subsolve: E_inner =" << common->E_inner[i].ToString();
    LOG_DBG3(sgp) << "Subsolve: ev_vector_min =" << ev_vector_min.ToString();
    LOG_DBG3(sgp) << "Subsolve: l1_min =" <<l1_min <<" l2_min = "<<l2_min<<" l3_min = "<<l3_min;
    LOG_DBG3(sgp) << "Subsolve: rho1_min =" <<rho1_min <<" rho2_min = "<<rho2_min<<" rho3_min = "<<rho3_min;
    LOG_DBG3(sgp) << "Subsolve: atmp =" << atmp.ToString();

    // calculate full model function for globalization
    cost += obj_min + dfdL.Trace();
  }
  return cost;
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
    //E_tmp = E_tmptmp - common->E_outer[index];
    //E_tmp *= common->tau;
    //BB += E_tmp;
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

double SGPApproximation::CalcAnalyticSol_FOMO_Top(double &rho1, double &rho2, double & rho, Vector<double> & ev,  Matrix<double> & ev_vector,  Eval eval, Matrix<double> BB, double theta_inner, double ppen, int index, double & Vloc) {

  Matrix<double> E_tmp(2,2);
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

    DesignSpace* space = common->optimization->GetDesign();
    unsigned int mech11 = space->FindDesign(DesignElement::MECH_11);
    unsigned int mech22 = space->FindDesign(DesignElement::MECH_22);
    //unsigned int mech33 = space->FindDesign(DesignElement::MECH_33);
    unsigned int dens = space->FindDesign(DesignElement::DENSITY);

    double lower_mech11 = space->GetDesignElement(mech11*common->n_elem+index)->GetLowerBound();
    double upper_mech11 = space->GetDesignElement(mech11*common->n_elem+index)->GetUpperBound();
    double lower_mech22 = space->GetDesignElement(mech22*common->n_elem+index)->GetLowerBound();
    double upper_mech22 = space->GetDesignElement(mech22*common->n_elem+index)->GetUpperBound();
    //double lower_mech33 = space->GetDesignElement(mech33*common->n_elem+index)->GetLowerBound();
    //double upper_mech33 = space->GetDesignElement(mech33*common->n_elem+index)->GetUpperBound();
    double lower_dens = space->GetDesignElement(dens*common->n_elem+index)->GetLowerBound();
    double upper_dens = space->GetDesignElement(dens*common->n_elem+index)->GetUpperBound();

    // Calculate Vloc
    //Vloc = upper_mech11 + upper_mech22 + upper_mech33 - lower_mech11 - lower_mech22 - lower_mech33;
    Vloc = 1.;
    double scale = 1. - common->L[index][0][0] - common->L[index][1][1] - common->L[index][2][2];

    // projection to bounds
    rho1 = std::min(std::max(lower_mech11, .99 * scale * sqrB1/(sqrB1+sqrB2+sqrB3)),upper_mech11-lower_mech11);
    rho2 = std::min(std::max(lower_mech22, .99 * scale * sqrB2/(sqrB1+sqrB2+sqrB3)),upper_mech22-lower_mech22);

    double rtmp = ev[0]/rho1+ev[1]/rho2+tmp[2][2]/(Vloc-rho1-rho2);

    rho = pow(psimp * rtmp/ppen,1./(psimp+1.));
    rho = std::min(upper_dens,std::max(lower_dens,rho));
    result =  pow(rho,(-psimp))*rtmp + ppen * rho;
    break;
  }
  default:
    break;
  }

  return result;
}

double SGPApproximation::CalcAnalyticSol_FOMO(double &rho1, double &rho2, double &rho3, Vector<double> & ev,  Matrix<double> & ev_vector,  Eval eval, Matrix<double> BB, Matrix<double> L,double theta_inner, double ppen, int index, double & Vloc) {

  Matrix<double> E_tmp(3,3);
  double result = 0;
    switch(eval)
  {
  case FUNC:
  {
    Matrix<double> tmp(BB);
    ev.Resize(3);
    ev_vector.Resize(3,3);
    common->optimization->GetDesign()->designMaterial->RotateTensor(tmp,DesignElement::NO_DERIVATIVE, DesignMaterial::HILL_MANDEL, DesignMaterial::CW, true, theta_inner);
    E_tmp[0][0] = tmp[0][0];
    E_tmp[0][1] = tmp[0][1];
    E_tmp[0][2] = tmp[0][2];
    E_tmp[1][0] = tmp[1][0];
    E_tmp[1][1] = tmp[1][1];
    E_tmp[1][2] = tmp[1][2];
    E_tmp[2][0] = tmp[2][0];
    E_tmp[2][1] = tmp[2][1];
    E_tmp[2][2] = tmp[2][2];
    E_tmp.eigenvaluesWithLapack(ev,&ev_vector);
    double sqrB1 = sqrt(ev[0]);
    double sqrB2 = sqrt(ev[1]);
    double sqrB3 = sqrt(ev[2]);
    double sqrPpen = sqrt(ppen);

    DesignSpace* space = common->optimization->GetDesign();
    unsigned int mech11 = space->FindDesign(DesignElement::MECH_11);
    unsigned int mech22 = space->FindDesign(DesignElement::MECH_22);
    unsigned int mech33 = space->FindDesign(DesignElement::MECH_33);

    double lower_mech11 = space->GetDesignElement(mech11*common->n_elem+index)->GetLowerBound();
    double upper_mech11 = space->GetDesignElement(mech11*common->n_elem+index)->GetUpperBound();
    double lower_mech22 = space->GetDesignElement(mech22*common->n_elem+index)->GetLowerBound();
    double upper_mech22 = space->GetDesignElement(mech22*common->n_elem+index)->GetUpperBound();
    double lower_mech33 = space->GetDesignElement(mech33*common->n_elem+index)->GetLowerBound();
    double upper_mech33 = space->GetDesignElement(mech33*common->n_elem+index)->GetUpperBound();

    // Calculate Vloc
    Vloc = upper_mech11 + upper_mech22 + upper_mech33 - lower_mech11 - lower_mech22 - lower_mech33;

    rho1 = std::min(std::max(lower_mech11, .99*sqrB1/sqrPpen),upper_mech11-lower_mech11);
    rho2 = std::min(std::max(lower_mech22, .99*sqrB2/sqrPpen),upper_mech22-lower_mech22);
    rho3 = std::min(std::max(lower_mech33, .99*sqrB3/sqrPpen),upper_mech33-lower_mech33);

    result = (ev[0]/rho1+ev[1]/rho2+ev[2]/rho3) + ppen * (rho1 + rho2 + rho3) + ppen * (L[0][0] + L[1][1] + L[2][2]);
    break;
  }
  default:
    break;
  }

  return result;
}

void SGPApproximation::CalcMinEigenvalue(StdVector<Matrix<double> > & E_in, Vector<double> & l_min) {
  Vector<double> ev;
  Matrix<double> ev_vector;
  ev.Resize(3);
  ev_vector.Resize(3,3);
  for (unsigned int i=0; i < common->n_elem; i++) {
    E_in[i].eigenvaluesWithLapack(ev,&ev_vector);
    l_min[i] = std::min(std::min(ev[0],ev[1]),ev[2]);
  }
}

double SGPApproximation::CalcAnalyticSol_FMO(double &rho1, double &rho2, double & rho3, Vector<double> & ev,  Matrix<double> & ev_vector,  Eval eval, Matrix<double> BB, Matrix<double> L, double ppen, int index, double & Vloc) {

  Matrix<double> E_tmp(3,3);
  double result = 0;
  switch(eval)
  {
  case FUNC:
  {
    Matrix<double> tmp(BB);
    ev.Resize(3);
    ev_vector.Resize(3,3);
    E_tmp[0][0] = BB[0][0];
    E_tmp[0][1] = BB[0][1];
    E_tmp[0][2] = BB[0][2];
    E_tmp[1][0] = BB[1][0];
    E_tmp[1][1] = BB[1][1];
    E_tmp[1][2] = BB[1][2];
    E_tmp[2][0] = BB[2][0];
    E_tmp[2][1] = BB[2][1];
    E_tmp[2][2] = BB[2][2];
    E_tmp.eigenvaluesWithLapack(ev,&ev_vector);
    double sqrB1 = sqrt(ev[0]);
    double sqrB2 = sqrt(ev[1]);
    double sqrB3 = sqrt(ev[2]);
    double sqrPpen = sqrt(ppen);

    DesignSpace* space = common->optimization->GetDesign();
    unsigned int mech11 = space->FindDesign(DesignElement::MECH_11);
    unsigned int mech22 = space->FindDesign(DesignElement::MECH_22);
    unsigned int mech33 = space->FindDesign(DesignElement::MECH_33);

    double lower_mech11 = space->GetDesignElement(mech11*common->n_elem+index)->GetLowerBound();
    double upper_mech11 = space->GetDesignElement(mech11*common->n_elem+index)->GetUpperBound();
    double lower_mech22 = space->GetDesignElement(mech22*common->n_elem+index)->GetLowerBound();
    double upper_mech22 = space->GetDesignElement(mech22*common->n_elem+index)->GetUpperBound();
    double lower_mech33 = space->GetDesignElement(mech33*common->n_elem+index)->GetLowerBound();
    double upper_mech33 = space->GetDesignElement(mech33*common->n_elem+index)->GetUpperBound();

    // Calculate Vloc
    Vloc = upper_mech11 + upper_mech22 + upper_mech33 - lower_mech11 - lower_mech22 - lower_mech33;

    rho1 = std::min(std::max(lower_mech11, .99*sqrB1/sqrPpen),upper_mech11 - lower_mech11);
    rho2 = std::min(std::max(lower_mech22, .99*sqrB2/sqrPpen),upper_mech22 - lower_mech22);
    rho3 = std::min(std::max(lower_mech33, .99*sqrB3/sqrPpen),upper_mech33 - lower_mech33);

    //result = 1./(Vloc - L[0][0] - L[1][1] - L[2][2])*(ev[0]/rho1+ev[1]/rho2+ev[2]/rho3) + ppen * Vloc * (rho1+rho2+rho3);
    result = (ev[0]/rho1+ev[1]/rho2+ev[2]/rho3) + ppen * (rho1 + rho2 + rho3) + ppen * (L[0][0] + L[1][1] + L[2][2]);
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
  //common->optimization->GetDesign()->designMaterial->RotateTensor(tmp,DesignElement::NO_DERIVATIVE, DesignMaterial::HILL_MANDEL, DesignMaterial::CW, true, theta_inner);

}

void SGPApproximation::CalcE_inner_Density_Rotangle(Matrix<double> & E_inner, Matrix<double> E_0, double theta_inner) {
  E_inner.Resize(3,3);
  Matrix<double> tmp(E_0);
  //don't know why VOIGT is used here
  common->optimization->GetDesign()->designMaterial->RotateTensor(tmp,DesignElement::NO_DERIVATIVE, DesignMaterial::VOIGT, DesignMaterial::CCW, true, theta_inner);
  E_inner = tmp;
}

std::string SGPApproximation::ToString()
{
  return Function::type.ToString(GetFunction()->GetType());
}

Function* SGPApproximation::GetFunction()
{
    return common->optimization->objectives.data[0];
}
