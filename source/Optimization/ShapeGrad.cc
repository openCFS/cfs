#include "Optimization/ShapeGrad.hh"
#include "Optimization/DesignSpace.hh"
#include "Optimization/DesignElement.hh"
#include "Optimization/OptimizationMaterial.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "Domain/domain.hh"
#include "PDE/mechPDE.hh"
#include "Forms/mechStressStrain.hh"
#include "Utils/nodestoresol.hh"
#include <string>

using namespace CoupledField;


DECLARE_LOG(shapeGrad)
DEFINE_LOG(shapeGrad, "shapeGrad")

DECLARE_LOG(topGrad)
DEFINE_LOG(topGrad, "topGrad")

ShapeGrad::ShapeGrad() : ErsatzMaterial()
{
  ParamNode* sg_pn = pn->Get("shapeGrad");

  topgrad = sg_pn != NULL && sg_pn->Has("topGrad") && sg_pn->Get("topGrad")->Get("enabled")->AsBool();
  
  mech_mat_ = NULL; // to be set in PostInit()
}


void ShapeGrad::PostInit()
{
  ErsatzMaterial::PostInit();
  
  mech_mat_ = dynamic_cast<OptMechMat*>(material); // just created in PostInit()
  assert(material != NULL);
}

void ShapeGrad::CalcObjectiveGradient(double* grad_out)
{
  CalcTopGrad();
}


double ShapeGrad::CalcTopGrad(unsigned int iter, unsigned int elems_removed_per_iteration,
    std::list<TopGradElementValues>& tgvalues, std::list<TopGradElementValues>& remtgvalues)
{
  unsigned int removed_counter = 0;
  unsigned int elements = design->GetNumberOfElements();
  DesignElement* de;

//  if(iter == 1)
//  {
//    SolveStateProblem();
//    CalcObjectiveGradient(NULL);
//
//    for(unsigned int i = 0; i < elements; i++)
//    {
//      // remember first computation of objective gradient for calculation of top_grad
//      de = &design->data[i];
//      de->SetInitialCostGradient(de->GetObjectiveGradient(DesignElement::PLAIN));
//    }
//  }

  // SolveAdjointProblem also computes the forward solution
  assert(excitations.GetSize() == 1); // no multiple load cases!
  int idx = 0;
  SolveAdjointProblem(excitations[idx]);

  // GetElementMatrix(GetForm(mech, mech, "linElastInt"), mechStiffness);

  // now compute the value of the topological gradient and the cost function

  double val = MultiplyForwardAndAdjointSolution(forward.data[idx], adjoint.data[idx], iter, elems_removed_per_iteration);
  cost->SetValue(val);

  // here, the topology gradient should be correctly calculated
  // and saved in the design-elements variable top_gradient

  // if elements < elems_removed_per_iteration*iter we want to remove more elements than available
  assert(elements > elems_removed_per_iteration * iter);

  std::list<TopGradElementValues>::iterator it;

  if(iter == 1)
  {
    TopGradElementValues tgev;
    for(unsigned int e = 0; e < elements; e++)
    {
      de = &design->data[e];
      tgev.elem_num = e;
      tgev.UpdateTopGrad(de->lsn->GetTopGradient());
      tgvalues.push_back(tgev);
    }
  }
  else
  {
    for (it = tgvalues.begin(); it != tgvalues.end(); it++)
    {
      de = &design->data[(*it).elem_num];
      (*it).UpdateTopGrad(de->lsn->GetTopGradient());
    }
  }

  // sort list, low topgradvalues are in front
  tgvalues.sort();

  int negative = 0;
  // TopGradElementValues tgelemtmp;

  for(unsigned int c = 0; c < elems_removed_per_iteration; c++)
  {
    // get element
    it = tgvalues.begin();

    if((*it).topgradvalue < 0.0) negative++;

    de = &design->data[(*it).elem_num];
    // set design to lower bound
    de->SetDesign(de->GetLowerBound());

    // remove element from list of all elements
    // remember which elements have been removed
    // remtgvalues.push_back( (*it) );
    tgvalues.pop_front();

    removed_counter++;
  }

  // volume constraint???

  double percent = (double) (elements - tgvalues.size()) / (double) elements * 100.0;

  std::cout << "iteration ";
  std::printf("%2i", iter);
  std::cout << ": elements removed: ";
  std::printf("%4i", removed_counter);
  std::cout << " (total: ";
  std::printf("%5.2f", percent);
  std::cout << "% of volume) | cost = ";
  std::printf("%14.2f", val);
  std::cout << " | neg: ";
  std::printf("%3i", negative);
  std::cout << " | rem: ";
  std::printf("%2i", (int) remtgvalues.size());
  std::cout << std::endl;

  return percent;
}

double ShapeGrad::MultiplyForwardAndAdjointSolution(Solution* sol_forward,
    Solution* sol_adjoint, unsigned int iter,
    unsigned int elems_removed_per_iteration)
{
  // for now, the application we consider is only mech
  Application app1 = MECH;
  Application app2 = MECH;

  unsigned int elements = design->GetNumberOfElements();
  // for cost function
  double sum = 0.0;

  StdVector<SingleVector*>& sol_forward_elem = sol_forward->elem[app1];
  StdVector<SingleVector*>& sol_adjoint_elem = sol_adjoint->elem[app2];
  Vector<double> tmp;

  assert(sol_forward_elem.GetSize() != 0);
  assert(sol_adjoint_elem.GetSize() != 0);

  // create an element list to gain the iterator in the loop
  ElemList elemList(domain->GetGrid());

  // for calculation of strains = partial derivatives of displacement u
  Vector<double> intPoint;
  Vector<double> forward_strain;
  Vector<double> adjoint_strain;
  Matrix<double> elem_sol_forward_matrix;
  Matrix<double> elem_sol_adjoint_matrix;

  // we need a NodeStoreSol to get a element*matrix*. We use
  // the current one and explicitly set the raw solution vector from the forward problem
  NodeStoreSol<double>* node_store_sol = dynamic_cast<NodeStoreSol<double>* >(ErsatzMaterial::mech->getPDESolution());

  // transform the subType of the pde (FULL, PLANE_STRAIN etc. -> defined in General/environment.hh)
  SubTensorType type;
  String2Enum(mech->GetSubType(),type);

  BaseMaterial* material = mech->getPDEMaterialData()[regionIds[0]]; // TODO: extend for multi-region-optimization if necessary
  double poisson;
  material->GetScalar(poisson, MECH_POISSON, Global::REAL);
  double e_mod;
  material->GetScalar(e_mod, MECH_EMODULUS, Global::REAL);

  // calculate lame parameters (-> Kaltenbacher)
  // from poisson ratio and elasticity modulus
  double lambda = (poisson * e_mod) / ((1 + poisson) * (1 - 2 * poisson));
  double mu     = e_mod / (2 * (1 + poisson));

  // will contain the calculated strains = partial derivatives of displacement u
  MechStressStrain<double>* strain = new MechStressStrain<double>(material,type);

  // calculate topGrad for each element
  for(unsigned int e = 0; e < elements; e++)
  {
    // get solution for this element (here only for evaluation of cost function)
    Vector<double>& sol_forward_vec = dynamic_cast<Vector<double>& >(*sol_forward_elem[e]);
    Vector<double>& sol_adjoint_vec = dynamic_cast<Vector<double>& >(*sol_adjoint_elem[e]);

    // get current design element
    DesignElement* de= &design->data[e];

    // Calculate element strain results for forward and adjoint problem
    elemList.SetElement(de->elem);
    EntityIterator it = elemList.GetIterator();

    it.GetElem()->ptElem->GetCoordMidPoint(intPoint);

    // set element solution to matrix
    node_store_sol->GetElemSolutionAsMatrix(elem_sol_forward_matrix, it, 0, sol_forward->raw[MECH]);
    node_store_sol->GetElemSolutionAsMatrix(elem_sol_adjoint_matrix, it, 0, sol_adjoint->raw[MECH]);

    strain->SetIntPoint(intPoint);

    // calculate strain for forward problem
    strain->SetActElemSol(elem_sol_forward_matrix);
    strain->CalcStrainVec(forward_strain, 1, it);

    // adjoint problem
    strain->SetActElemSol(elem_sol_adjoint_matrix);
    strain->CalcStrainVec(adjoint_strain, 1, it);

    LOG_DBG3(topGrad) << "elem " << de->elem->elemNum << " forward strain: " << forward_strain.ToString();
    LOG_DBG3(topGrad) << "elem " << de->elem->elemNum << " adjoint strain: " << adjoint_strain.ToString();

    double topG = 0.0, div_u = 0.0, div_v = 0.0, tr_sig_u = 0.0, tr_eps_v = 0.0, temp_calc = 0.0;
    const double pi = std::acos( -1.0 );

    switch(type)
    {
    case PLANE_STRAIN: // top-grad 2d *neumann* garreau
      assert(forward_strain.GetSize() == 3);

      div_u = forward_strain[0] + forward_strain[1];
      div_v = adjoint_strain[0] + adjoint_strain[1];

      tr_sig_u = 2 * (lambda + mu) * div_u;
      tr_eps_v = div_v;

      temp_calc = (lambda - mu) * tr_sig_u * tr_eps_v;

      topG  = lambda * div_u * div_v;
      topG += 2 * mu * (forward_strain[0] * adjoint_strain[0] + forward_strain[1] * adjoint_strain[1]);
      topG += 4 * mu * (forward_strain[2] * adjoint_strain[2]);
      topG *= 4 * mu;

      topG += temp_calc;

      topG *= pi * (lambda + 2 * mu) / (2 * (lambda + mu) * mu);
      break;

    case FULL: // top-grad 3d *neumann*
      assert(forward_strain.GetSize() == 6);

      // we have that u = v in case of compliance
      div_u = forward_strain[0] + forward_strain[1] + forward_strain[2];
      div_v = div_u;
      // div_v = adjoint_strain[0] + adjoint_strain[1] + adjoint_strain[2];

      tr_sig_u = (3 * lambda + 2 * mu) * div_u;
      tr_eps_v = div_v;

      temp_calc = (3 * lambda - 2 * mu) * tr_sig_u * tr_eps_v;

      topG  = lambda * div_u * div_v;
      topG += 2 * mu * (forward_strain[0] * adjoint_strain[0] +
          forward_strain[1] * adjoint_strain[1] + forward_strain[2] * adjoint_strain[2]);
      topG += 4 * mu * (forward_strain[3] * adjoint_strain[3] +
          forward_strain[4] * adjoint_strain[4] + forward_strain[5] * adjoint_strain[5]);
      topG *= 20 * mu;
      topG += temp_calc;

      topG *= pi * (lambda + 2 * mu) / ((9 * lambda + 14 * mu) * mu);
      break;

    case AXI:
      EXCEPTION("SubType AXI not implemented!");
    case PLANE_STRESS:
      EXCEPTION("SubType PLANE_STRESS not implemented!");
    default:
      EXCEPTION("SubType not implemented!");
    }

    /* store the value in the top_gradient of the current design-element */
    de->lsn->SetTopGradient(topG);

    // compliance: -u_e^T K_0 u_e
    // mechanism:(-)l_e^T K_0 u_e
    tmp = mech_mat_->MechStiffness(de->elem) * sol_adjoint_vec;
    double sp = sol_forward_vec * tmp;

    // when doing complex Jensen 22.07.07 shows that we always have
    // 2 * Re(lamda * grad S * u)
    // the factor gives the negative sign
    double value = ((std::complex<double>) sp).real();
    sum += value;
  }

  delete strain;

  // sum = value of cost function
  return sum;
}
