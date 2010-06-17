#include "Optimization/PiezoSIMP.hh"
#include "Optimization/SIMP.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/OptimizationMaterial.hh"
#include "Driver/basedriver.hh"
#include "Driver/singleDriver.hh"
#include "Forms/baseForm.hh"
#include "Forms/linearForm.hh"
#include "Forms/linNeumannInt.hh"
#include "Forms/linPressureInt.hh"
#include "Domain/domain.hh"
#include "PDE/mechPDE.hh"
#include "PDE/elecPDE.hh"
#include "Driver/stdSolveStep.hh"
#include "General/exception.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "MatVec/basematrix.hh"
#include "MatVec/stdmatrix.hh"
#include "Utils/nodestoresol.hh"
#include "OLAS/algsys/basesystem.hh"
#include "boost/lexical_cast.hpp"

using namespace CoupledField;

using std::complex;

DECLARE_LOG(simp)


PiezoSIMP::PiezoSIMP()
{
  elec = dynamic_cast<ElecPDE*>(domain->GetSinglePDE("electrostatic"));
  pdes[ELEC] = elec;

  for(unsigned int r = 0; r < regionIds.GetSize(); r++)
  {
    GetForm(regionIds[r], pde, elec, "linPiezoCoupling")->SetSolDependent(true);
    GetForm(regionIds[r], elec, elec, "linElecInt")->SetSolDependent(true);
  }
  // The linear forms (pressure, charge density) are set in SoluctionRef::Init()

  // validate the transfer functions
  if(design->design.Find(DesignElement::DENSITY) >= 0)
  {
    if(design->GetTransferFunction(DesignElement::DENSITY, MECH, false) == NULL)
      throw Exception("miss transfer function for densitiy and mechanic");

    if(design->GetTransferFunction(DesignElement::DENSITY, ELEC, false) == NULL &&
        design->GetTransferFunction(DesignElement::POLARIZATION, ELEC, false) == NULL)
      throw Exception("miss transfer function for densitiy/polarization and electrostatic");

    if(design->GetTransferFunction(DesignElement::DENSITY, PIEZO_COUPLING, false) == NULL)
      throw Exception("miss transfer function for densitiy and coupling");
  }
  if(design->design.Find(DesignElement::POLARIZATION) >= 0)
  {
    if(design->GetTransferFunction(DesignElement::POLARIZATION, PIEZO_COUPLING, false) == NULL)
      throw Exception("miss transfer function for polarization and coupling");
  }

  piezo_mat_ = NULL; // to be set in PostInit()
}

PiezoSIMP::~PiezoSIMP()
{
  if(elecRHS.vec != NULL) { delete elecRHS.vec;   elecRHS.vec = NULL; }
}

void PiezoSIMP::PostInit()
{
  if(harmonic) elecRHS.Init<complex<double> >(design, CHARGE_DENSITY); // mechRHS in SIMP!
          else elecRHS.Init<double>(design, CHARGE_DENSITY);

  SIMP::PostInit();

  // just created in PostInit
  piezo_mat_ = dynamic_cast<OptPiezoMat*>(material);
  assert(piezo_mat_ != NULL);
}


template <class T>
double PiezoSIMP::CalcElecEnergy(Excitation& excite)
{
  // calculate the element sum of p^T K_pp p or p^T K_pp p^*
  // here we do <K_pp p, p> which is equivalent as K_pp is self adjoined
  // (it is real and K_pp = K_pp^T). As the complex scalar product is
  // <x, y> = \sum x_i y_i^* we have 
  // <K_pp p, p> =  p^T K_pp p in the real case and p^T K_pp p^* in complex 
  
  assert(design->design.GetSize() == 1);
  
  DesignElement::Type dt = design->design[0];
  TransferFunction* tf = design->GetTransferFunction(dt, ELEC);
  
  // our solution vectors
  StdVector<SingleVector*>& all_p = forward.Get(excite)->elem[ELEC];
  Matrix<T> mat(all_p[0]->GetSize(), all_p[0]->GetSize());
  
  Vector<T> mat_vec(all_p[0]->GetSize()); // for the temporary K_pp * p result
  
  // traverse over our elements
  // create an element list to gain the iterator in the loop
  ElemList elemList(grid);
  
  double sum = 0.0;
  
  // for ParamMat we need the derivative w.r.t. every design variable, else the base loop is only run once
  for(int i = 0, n = design->data.GetSize(); i < n; i++) 
  {
    Vector<T>& p_vec = dynamic_cast<Vector<T>& >(*all_p[i]);
    DesignElement* de = &design->data[i];

    LOG_DBG3(simp) << "CalcElecEnergy: e=" << i << " 'density'=" << tf->Transform(de) << " p=" << p_vec.ToString();
    
    Assign(mat, piezo_mat_->ElecStiffness(de->elem, 1), tf->Transform(de));
    
    LOG_DBG3(simp) << "CalcElecEnergy: mat: " << mat.ToString();
    
    // see comments above why we calculate <K_pp p, p>
    mat_vec = mat * p_vec;
    T result;
    mat_vec.Inner(p_vec, result); // in the complex case with the conjugate complex p_vec
    LOG_DBG3(simp) << "CalcElecEnergy: " << sum << " + " << result;
    assert( abs(((complex<double>) result).imag()) < 1e-6 * abs(((complex<double>) result).real())); // check imag = real
    
    sum += ((complex<double>) result).real();
  }
  // the objective is actually 0.5 p^T K_pp p
  sum *= 0.5;
  
  LOG_DBG(simp) << "CalcElecEnergy: returns " << sum << " * "
                << excite.GetFactor(objectives.Get(Objective::ELEC_ENERGY)) << " = "
                << sum * excite.GetFactor(objectives.Get(Objective::ELEC_ENERGY));
  
  sum *= excite.GetFactor(objectives.Get(Objective::ELEC_ENERGY));
  return sum;
}

/** Sets -K_pp p or -K_pp p^* for ELEC_ENERGY */
template <class T>
void PiezoSIMP::ConstructAdjointRHS(Excitation& excite, Objective* cost)
{
  assert(objectives.Has(Objective::ELEC_ENERGY));
  assert(design->design.GetSize() == 1);

  DesignElement::Type dt = design->design[0];
  TransferFunction* tf = design->GetTransferFunction(dt, ELEC);

  // our solution vectors
  StdVector<SingleVector*>& all_p = forward.Get(excite)->elem[ELEC];
  Matrix<T> mat(all_p[0]->GetSize(), all_p[0]->GetSize()); // store K_pp

  // define our new adjont RHS
  Vector<T> rhs(forward.Get(excite)->GetVector(Solution::RHS_VECTOR)->GetSize()); // is initialized
  // TODO: up to now the rhs is not PDE specific but for the whole system!!

  Vector<T> mat_vec(all_p[0]->GetSize()); // for the temporary K_pp * p 

  // traverse over our elements
  // create an element list to gain the iterator in the loop
  ElemList elemList(grid);

  shared_ptr<EqnMap> eqn_map = elec->GetEqnMap();

  for(int e = 0, n = design->data.GetSize(); e < n; e++) 
  {
    Vector<T>& p_vec = dynamic_cast<Vector<T>& >(*all_p[e]); // element solution size
    DesignElement* de = &design->data[e];

    // gain +K_pp(rho), the plus because K_pp is only in the piezo -K_pp
    Assign(mat, piezo_mat_->ElecStiffness(de->elem, 1), tf->Transform(de));

    // in the complex case with the conjugate complex
    mat.MultInner(p_vec, mat_vec);
    
    LOG_DBG3(simp) << "CARHS: rhs before setting: " << rhs.ToString();

    // go over all nodes of the element to sum it up on the rhs vector!
    for(unsigned int n = 0, nde =de->elem->connect.GetSize(); n < nde; n++)
    {
      unsigned int nn = de->elem->connect[n]; // nodenumber

      // get equation number within global rhs
      int eqnr = eqn_map->GetNodeEqn(nn, 1); // Equations are 1-based! dof = 1 for phi
      LOG_DBG3(simp) << "CARHS: Node: " << nn << " eqnr (pure): " << eqnr;

      if(eqnr == 0) continue; // this means a homogeneous case and RHS has not to be set

      // negative eqnr is mapped to its positive equivalent eqnr because of constraints
      if(eqnr < 0) eqnr = eqnr * (-1);

      //LOG_DBG3(simp) << "CARHS: Node: " << nn << " assigned to rhs[" << eqnr << "] + mat_vec[n] -> " << (rhs[eqnr-1] - mat_vec[n]);

      rhs[eqnr-1] += -1.0 * mat_vec[n]; // OLAS is 0-based, the equation number is 1 based, the rhs is "-dJ/du")
    }
  }

  // RHS has to be applied
  LOG_DBG2(simp) << "CARHS: final rhs before setting: " << rhs.ToString();
  assemble_->GetAlgSys()->InitRHS(rhs);
  LOG_DBG2(simp) << "CARHS: rhs after setting: " << rhs.ToString();
  
  return;
}


void PiezoSIMP::CalcObjectiveGradient(Excitation& excite, Objective* cost)
{
  double factor = excite.GetWeightedFactor(cost);
  
  switch(cost->GetType())
  {
  case Objective::ELEC_ENERGY:
    factor *= 0.5; // no break! -> J = 0.5 p^T K_pp p
  case Objective::OUTPUT:
  case Objective::DYNAMIC_OUTPUT:
  case Objective::ENERGY_FLUX:
  case Objective::GLOBAL_DYNAMIC_COMPLIANCE:
  {
    LOG_DBG2(simp) << "PS::CalcObjectiveGradient(idx=" << excite.index << ") factor = "
                   << excite.normalized_weight << " * " << excite.GetFactor(cost) << " -> " << factor;

    // we see at mechRHS and elecRHS if we have such an excitation
    SurfaceRef* mech_rhs = mechRHS.valid ? &mechRHS : NULL;
    SurfaceRef* elec_rhs = elecRHS.valid ? &elecRHS : NULL;

    LOG_TRACE2(simp) << "PS:CalcOutputGradient: output -> sensitive rhs's: "
                     << mechRHS.ToString(1) << ", " << elecRHS.ToString(1);

    // normally the gradient is sol^T A' sol^* + (2 Re) lamda^T K' sol.
    // 2 Re{...} only in the harmonic case.
    // lambda = (lambda_u lambda_p)^T
    // sol    = (u p)^T
    // K      = first row: S_uu K_up, second row: K_pu K_pp
    // A      = matrix from objective. Often L or for ElecEnergy K_pp
    
    // sol = displacement and potential and K = K_uu, K_up, K_up^T, K_pp
    // we calcualate the 4 individual vec Mat vec via CalcU1KU2() and sum it up
    for(unsigned int i = 0; i < design->design.GetSize(); i++)
    {
      DesignElement::Type dt = design->design[i];
      // we allow NULL for the transfer functions,
      // then the gradient is 0 what is done via Reset!
      TransferFunction*   tf = NULL;
      int                 res_idx = -1; // we can write the details as special resul by the detail

      // the standard adjoint solution (u1 in the gradient)
      Solution* adj = adjoint.Get(excite);
      // the solution (u2 in the gradient)
      Solution* sol = forward.Get(excite);

      
      // sol^T A' sol^* only for elec energy = p^T K_pp' p (real) or p^T K_pp' p^* (harmonic)
      // we solve <K_pp p, p> as K_pp is real and symmetric. Hence the inner product makes p^* in the harmonic case!
      if(cost->GetType() == Objective::ELEC_ENERGY)
      {
        // p^T K_pp' p^*
        tf = design->GetTransferFunction(dt, ELEC, true);
        res_idx = GetSpecialResultIndex(ELEC, ELEC, CONJ_QUAD);
        CalcU1KU2(tf, sol->elem[ELEC], ELEC, sol->elem[ELEC], elec_rhs, factor, CONJ_QUAD, cost, NULL, res_idx);
      }
      
      // lambda_u * K_uu' * u
      tf = design->GetTransferFunction(dt, MECH, false); // we allow NULL
      res_idx = GetSpecialResultIndex(MECH, MECH);
      if(tf != NULL)
        CalcU1KU2(tf, adj->elem[MECH], MECH, sol->elem[MECH], mech_rhs, factor, STANDARD, cost, NULL, res_idx);

      // lambda_u * K_up' * p
      tf = design->GetTransferFunction(dt, PIEZO_COUPLING, false);
      res_idx = GetSpecialResultIndex(MECH, ELEC);
      if(tf != NULL)
        CalcU1KU2(tf, adj->elem[MECH], PIEZO_COUPLING, sol->elem[ELEC], mech_rhs, factor, STANDARD, cost, NULL, res_idx);

      // lambda_p * (K_up^T)' * u
      tf = design->GetTransferFunction(dt, PIEZO_COUPLING, false);
      res_idx = GetSpecialResultIndex(ELEC, MECH);
      if(tf != NULL)
        CalcU1KU2(tf, adj->elem[ELEC], PIEZO_COUPLING, sol->elem[MECH], elec_rhs, factor, STANDARD, cost, NULL, res_idx);

      // lambda_p * K_pp' * p
      tf = design->GetTransferFunction(dt, ELEC, false);
      res_idx = GetSpecialResultIndex(ELEC, ELEC);
      if(tf != NULL)
        CalcU1KU2(tf, adj->elem[ELEC], ELEC, sol->elem[ELEC], elec_rhs, factor, STANDARD, cost, NULL, res_idx);
    }
  }
  break;

  default:
    SIMP::CalcObjectiveGradient(excite, cost);
  }
}

template <class T>
void PiezoSIMP::SetElementK(DesignElement* de, Application app, DenseMatrix* mat_out, CalcMode calcMode, bool derivative)
{
  TransferFunction* tf = design->GetTransferFunction(de->GetType(), app);
  double factor = derivative ? tf->Derivative(de, DesignElement::SMART) : tf->Transform(de, DesignElement::SMART);

  Matrix<T>& out = dynamic_cast<Matrix<T>& >(*mat_out);

  switch(app)
  {
  case ELEC:
    if(calcMode == CONJ_QUAD)
      Assign(out, piezo_mat_->ElecStiffness(de->elem, 1), factor); // we need the K_pp matrix
    else // STANDARD and GRAD_N
      Assign(out, piezo_mat_->ElecStiffness(de->elem, -1), factor); // we need the -K_pp matrix
    break;

  case PIEZO_COUPLING:
  {
    const Matrix<double>& coupledStiffness = piezo_mat_->CoupledStiffness(de->elem);
    // we see on the size of in if we cave to be transposed!
    if(out.GetNumCols() == coupledStiffness.GetNumCols())
      Assign(out, coupledStiffness, factor);
    else
      Assign(out, piezo_mat_->CoupledStiffnessTransposed(de->elem), factor);
    break;
  }

  default:
    // mech and surface normal matrix are handled in SIMP
    SIMP::SetElementK(de, app, mat_out, calcMode, derivative);
    return; // all calculation done there (or assert!)
  }

  LOG_DBG2(simp) << "PiezoSIMP::SetElementK elem: " << de->elem->elemNum << " app: " << application.ToString(app) << " factor: " << factor;
}




