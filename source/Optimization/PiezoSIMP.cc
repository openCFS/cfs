#include <assert.h>
#include <stddef.h>
#include <cmath>
#include <map>
#include <ostream>
#include <string>

#include "DataInOut/Logging/LogConfigurator.hh"
#include "Domain/Domain.hh"
#include "Domain/ElemMapping/Elem.hh"
#include "Domain/ElemMapping/EntityLists.hh"
#include "Driver/Assemble.hh"
#include "General/Enum.hh"
#include "General/defs.hh"
#include "General/Environment.hh"
#include "General/Exception.hh"
#include "MatVec/SingleVector.hh"
#include "MatVec/exprt/xpr2.hh"
#include "MatVec/Matrix.hh"
#include "MatVec/Vector.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Excitation.hh"
#include "Optimization/OptimizationMaterial.hh"
#include "Optimization/PiezoSIMP.hh"
#include "Optimization/SIMP.hh"
#include "Optimization/TransferFunction.hh"
#include "PDE/SinglePDE.hh"
#include "PDE/ElecPDE.hh"
#include "Utils/StdVector.hh"

namespace CoupledField {
class DenseMatrix;
}  // namespace CoupledField

using namespace CoupledField;

EXTERN_LOG(simp)

using std::complex;

PiezoSIMP::PiezoSIMP()
{
  for(unsigned int r = 0; r < design->GetRegionIds().GetSize(); r++)
  {
    assert(false);
    // GetForm(design->GetRegionIds()[r], pde, elec, "linPiezoCoupling")->SetSolDependent(true);
    // GetForm(design->GetRegionIds()[r], elec, elec, "linGradBDBInt")->SetSolDependent(true);
  }
  // The linear forms (pressure, charge density) are set in SoluctionRef::Init()

  // validate the transfer functions
  if(design->FindDesign(DesignElement::DENSITY, false) >= 0)
  {
    if(design->GetTransferFunction(DesignElement::DENSITY, App::MECH, false) == NULL)
      throw Exception("miss transfer function for densitiy and mechanic");

    if(design->GetTransferFunction(DesignElement::DENSITY, App::ELEC, false) == NULL &&
        design->GetTransferFunction(DesignElement::POLARIZATION, App::ELEC, false) == NULL)
      throw Exception("miss transfer function for densitiy/polarization and electrostatic");

    if(design->GetTransferFunction(DesignElement::DENSITY, App::PIEZO_COUPLING, false) == NULL)
      throw Exception("miss transfer function for densitiy and coupling");
  }
  if(design->FindDesign(DesignElement::POLARIZATION, false) >= 0)
  {
    if(design->GetTransferFunction(DesignElement::POLARIZATION, App::PIEZO_COUPLING, false) == NULL)
      throw Exception("miss transfer function for polarization and coupling");
  }
}

PiezoSIMP::~PiezoSIMP()
{
  if(elecRHS.vec != NULL) { delete elecRHS.vec;   elecRHS.vec = NULL; }
}

void PiezoSIMP::PostInit()
{
  // ignores the SetPDE() framework :(
  if(context->IsComplex())
    elecRHS.Init<complex<double> >(design, App::CHARGE_DENSITY); // mechRHS in SIMP!
  else elecRHS.Init<double>(design, App::CHARGE_DENSITY);

  SIMP::PostInit();
}


template <class T>
double PiezoSIMP::CalcElecEnergy(Excitation& excite)
{
  // calculate the element sum of p^T K_pp p or p^T K_pp p^*
  // here we do <K_pp p, p> which is equivalent as K_pp is self adjoint
  // (it is real and K_pp = K_pp^T). As the complex scalar product is
  // <x, y> = \sum x_i y_i^* we have 
  // <K_pp p, p> =  p^T K_pp p in the real case and p^T K_pp p^* in complex 
  assert(design->design.GetSize() == 1);
  PiezoElecMat* pem = dynamic_cast<PiezoElecMat*>(context->mat); // don't cache!
  
  DesignElement::Type dt = design->design[0].design;
  TransferFunction* tf = design->GetTransferFunction(dt, App::ELEC);
  
  // our solution vectors
  StdVector<SingleVector*>& all_p = forward.Get(excite)->elem[App::ELEC];
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

    LOG_DBG3(simp) << "CalcElecEnergy: e=" << i << " 'density'=" << tf->Transform(de, DesignElement::SMART) << " p=" << p_vec.ToString();
    
    Assign(mat, pem->ElecStiffnessPos(de), tf->Transform(de, DesignElement::SMART));
    
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
void PiezoSIMP::ConstructAdjointRHS(Excitation& excite, Function* f)
{
  assert(objectives.Has(Objective::ELEC_ENERGY));
  assert(design->design.GetSize() == 1);

  DesignElement::Type dt = design->design[0].design;
  TransferFunction* tf = design->GetTransferFunction(dt, App::ELEC);
  PiezoElecMat* pem = dynamic_cast<PiezoElecMat*>(context->mat); // don't cache!

  // our solution vectors
  StdVector<SingleVector*>& all_p = forward.Get(excite)->elem[App::ELEC];
  Matrix<T> mat(all_p[0]->GetSize(), all_p[0]->GetSize()); // store K_pp

  // define our new adjont RHS
  Vector<T> rhs(forward.Get(excite)->GetVector(StateSolution::RHS_VECTOR)->GetSize()); // is initialized
  // TODO: up to now the rhs is not PDE specific but for the whole system!!

  Vector<T> mat_vec(all_p[0]->GetSize()); // for the temporary K_pp * p 

  // traverse over our elements
  // create an element list to gain the iterator in the loop
  ElemList elemList(grid);

  assert(false);
  shared_ptr<EqnMap> eqn_map; // FIXME = elec->GetEqnMap();

  for(int e = 0, n = design->data.GetSize(); e < n; e++) 
  {
    Vector<T>& p_vec = dynamic_cast<Vector<T>& >(*all_p[e]); // element solution size
    DesignElement* de = &design->data[e];

    // gain +K_pp(rho), the plus because K_pp is only in the piezo -K_pp
    Assign(mat, pem->ElecStiffnessPos(de), tf->Transform(de, DesignElement::SMART));

    // in the complex case with the conjugate complex
    mat.MultInner(p_vec, mat_vec);
    
    LOG_DBG3(simp) << "CARHS: rhs before setting: " << rhs.ToString();

    // go over all nodes of the element to sum it up on the rhs vector!
    for(unsigned int n = 0, nde =de->elem->connect.GetSize(); n < nde; n++)
    {
      unsigned int nn = de->elem->connect[n]; // nodenumber

      // get equation number within global rhs
      assert(false);
      int eqnr = -5; // FIXME = eqn_map->GetNodeEqn(nn, 1); // Equations are 1-based! dof = 1 for phi
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
  assert(false);
  // FIXME assemble_->GetAlgSys()->InitRHS(rhs);
  LOG_DBG2(simp) << "CARHS: rhs after setting: " << rhs.ToString();
  
  return;
}


double PiezoSIMP::CalcFunction(Excitation& excite, Function* f, bool derivative)
{
  double factor = excite.GetWeightedFactor(f);

  switch(f->GetType())
  {
  case Function::ELEC_ENERGY:
    if(!derivative)
    {
      return context->IsComplex() ? CalcElecEnergy<std::complex<double> >(excite) : CalcElecEnergy<double>(excite);
    }
    factor *= 0.5; // no break! -> J = 0.5 p^T K_pp p

  // TODO make it smarter to evaluate the functions as piezo functions and not mech functions!
  case Function::OUTPUT:
  case Function::SQUARED_OUTPUT:
  case Objective::DYNAMIC_OUTPUT:
  case Objective::ENERGY_FLUX:
  case Objective::GLOBAL_DYNAMIC_COMPLIANCE:
  case Objective::ABS_OUTPUT:
  {
    if(!derivative)
      return SIMP::CalcFunction(excite, f, derivative);

    LOG_DBG2(simp) << "PS::CF(idx=" << excite.index << ") factor = "
        << excite.normalized_weight << " * " << excite.GetFactor(f) << " -> " << factor;

    // we see at mechRHS and elecRHS if we have such an excitation
    DesignDependentRHS* mech_rhs = mechRHS.valid ? &mechRHS : NULL;
    DesignDependentRHS* elec_rhs = elecRHS.valid ? &elecRHS : NULL;

    LOG_TRACE2(simp) << "PS:CalcOutputGradient: output -> sensitive rhs's: "
        << mechRHS.ToString(1) << ", " << elecRHS.ToString(1);

    // normally the gradient is sol^T A' sol^* + (2 Re) lamda^T K' sol.
    // 2 Re{...} only in the harmonic case.
    // lambda = (lambda_u lambda_p)^T
    // sol    = (u p)^T
    // K      = first row: S_uu K_up, second row: K_pu K_pp
    // A      = matrix from objective. Often L or for ElecEnergy K_pp

    // sol = displacement and potential and K = K_uu, K_up, K_up^T, K_pp
    // we calculate the 4 individual vec Mat vec via CalcU1KU2() and sum it up
    //
    // in the FMO case we run only via one design, the transfer function is identity
    unsigned int size = context->dm != NULL ? 1 : design->design.GetSize();
    for(unsigned int i = 0; i < size; i++)
    {
      DesignElement::Type dt = design->design[i].design;
      // we allow NULL for the transfer functions,
      // then the gradient is 0 what is done via Reset!
      TransferFunction*   tf = NULL;
      int                 res_idx = -1; // we can write the details as special resul by the detail

      // the standard adjoint solution (u1 in the gradient)
      StateSolution* adj = adjoint.Get(excite, f);
      // the solution (u2 in the gradient)
      StateSolution* sol = forward.Get(excite, NULL);


      // sol^T A' sol^* only for elec energy = p^T K_pp' p (real) or p^T K_pp' p^* (harmonic)
      // we solve <K_pp p, p> as K_pp is real and symmetric. Hence the inner product makes p^* in the harmonic case!
      if(f->GetType() == Objective::ELEC_ENERGY)
      {
        // p^T K_pp' p^*
        tf = design->GetTransferFunction(dt, App::ELEC, true);
        res_idx = GetSpecialResultIndex(App::ELEC, App::ELEC, CONJ_QUAD);
        CalcU1KU2(tf, sol->elem[App::ELEC], App::ELEC, sol->elem[App::ELEC], elec_rhs, factor, CONJ_QUAD, f, res_idx);
      }

      // squared output gradient 2 * <u,l> * <u,l>'
      if(f->GetType() == Function::SQUARED_OUTPUT)
      {
         f->SetType(Function::OUTPUT);
         factor *= 2. * SIMP::CalcFunction(excite, f, false);
         f->SetType(Function::SQUARED_OUTPUT);
      }
      // lambda_u * K_uu' * u
      tf = design->GetTransferFunction(dt, App::MECH, false); // we allow NULL
      res_idx = GetSpecialResultIndex(App::MECH, App::MECH);
      if(tf != NULL)
        CalcU1KU2(tf, adj->elem[App::MECH], App::MECH, sol->elem[App::MECH], mech_rhs, factor, STANDARD, f, res_idx);

      // lambda_u * K_up' * p
      tf = design->GetTransferFunction(dt, App::PIEZO_COUPLING, false);
      res_idx = GetSpecialResultIndex(App::MECH, App::ELEC);
      if(tf != NULL)
        CalcU1KU2(tf, adj->elem[App::MECH], App::PIEZO_COUPLING, sol->elem[App::ELEC], mech_rhs, factor, STANDARD, f, res_idx);

      // lambda_p * (K_up^T)' * u
      tf = design->GetTransferFunction(dt, App::PIEZO_COUPLING, false);
      res_idx = GetSpecialResultIndex(App::ELEC, App::MECH);
      if(tf != NULL)
        CalcU1KU2(tf, adj->elem[App::ELEC], App::PIEZO_COUPLING, sol->elem[App::MECH], elec_rhs, factor, STANDARD, f, res_idx);

      // lambda_p * K_pp' * p
      tf = design->GetTransferFunction(dt, App::ELEC, false);
      res_idx = GetSpecialResultIndex(App::ELEC, App::ELEC);
      if(tf != NULL)
        CalcU1KU2(tf, adj->elem[App::ELEC], App::ELEC, sol->elem[App::ELEC], elec_rhs, factor, STANDARD, f, res_idx);
    }
  }
  return 0.0; // we calculated only derivatives

  default: // return below as we don't implement
    break;
  }
  return SIMP::CalcFunction(excite, f, derivative);
}




template <class T1, class T2>
void PiezoSIMP::SetElementK(Function* f, DesignElement* de, const TransferFunction* tf, App::Type app, DenseMatrix* mat_out, bool derivative, CalcMode calcMode, double ev)
{
  PiezoElecMat* pem = dynamic_cast<PiezoElecMat*>(f->ctxt->mat); // don't cache!

  double factor = derivative ? tf->Derivative(de, DesignElement::SMART) : tf->Transform(de, DesignElement::SMART);

  Matrix<T1>& out = dynamic_cast<Matrix<T1>& >(*mat_out);

  switch(app)
  {
  case App::ELEC:
    if(calcMode == CONJ_QUAD)
      Assign(out, pem->ElecStiffnessPos(de), factor); // we need the K_pp matrix
    else // STANDARD and GRAD_N
      Assign(out, pem->ElecStiffnessNeg(de), factor); // we need the -K_pp matrix
    break;

  case App::PIEZO_COUPLING:
  {
    assert(out.GetNumCols() != out.GetNumRows());
    if(out.GetNumRows() > out.GetNumCols())
      Assign(out, pem->CoupledStiffness(de), factor);
    else
      Assign(out, pem->CoupledStiffnessTransposed(de), factor);
    break;
  }

  default:
    // mech and surface normal matrix are handled in SIMP
    SIMP::SetElementK(f, de, tf, app, mat_out, derivative, calcMode, ev);
    return; // all calculation done there (or assert!)
  }

  LOG_DBG2(simp) << "PiezoSIMP::SetElementK elem: " << de->elem->elemNum << " app: " << application.ToString(app) << " factor: " << factor;
}

#ifndef _MSC_VER
  template void PiezoSIMP::ConstructAdjointRHS<double>(Excitation& excite, Function* cost);
  template void PiezoSIMP::ConstructAdjointRHS<Complex>(Excitation& excite, Function* cost);
#else
  template void PiezoSIMP::ConstructAdjointRHS<std::complex<double>>(Excitation& excite, Function* f);
  template void PiezoSIMP::ConstructAdjointRHS<double>(Excitation& excite, Function* f);
#endif



