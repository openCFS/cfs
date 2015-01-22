#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <map>
#include <ostream>
#include <string>

#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/Logging/log.hpp"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Domain/Domain.hh"
#include "Domain/ElemMapping/Elem.hh"
#include "Domain/ElemMapping/EntityLists.hh"
#include "Domain/ElemMapping/SurfElem.hh"
#include "Driver/Assemble.hh"
#include "Driver/SolveSteps/BaseSolveStep.hh"
#include "Driver/FormsContexts.hh"
#include "Forms/LinForms/LinearForm.hh"
#include "General/Enum.hh"
#include "General/defs.hh"
#include "General/Environment.hh"
#include "General/Exception.hh"
#include "MatVec/SingleVector.hh"
#include "MatVec/exprt/xpr1.hh"
#include "MatVec/Matrix.hh"
#include "MatVec/Vector.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Design/DesignStructure.hh"
#include "Optimization/Excitation.hh"
#include "Optimization/Function.hh"
#include "Optimization/OptimizationMaterial.hh"
#include "Optimization/SIMP.hh"
#include "Optimization/StressConstraint.hh"
#include "Optimization/TransferFunction.hh"
#include "PDE/SinglePDE.hh"
#include "PDE/MechPDE.hh"
#include "Utils/StdVector.hh"
#include "Utils/mathParser/mathParser.hh"
#include "Utils/tools.hh"

namespace CoupledField {
class DenseMatrix;
}  // namespace CoupledField

using namespace CoupledField;

using std::complex;


DECLARE_LOG(conditions)

DECLARE_LOG(simp)
DEFINE_LOG(simp, "simp")


SIMP::SIMP() : ErsatzMaterial()
{
}

SIMP::~SIMP()
{
  if(mechRHS.vec != NULL)       { delete mechRHS.vec;  mechRHS.vec = NULL;  }
}

void SIMP::PostInit()
{
  PtrParamNode simp_pn = pn->Get("SIMP", ParamNode::PASS);

  // There might be a filter regularization based on the design element.
  if(simp_pn)
  {
    if(simp_pn->HasByVal("regularization", "type", "filter"))
    {
      ParamNodeList list = simp_pn->Get("regularization")->GetList("filter");
      // this is save for design=polarization
      for(unsigned int i = 0; i < list.GetSize(); i++)
      {
        if(structure_ == NULL)
          structure_ = new DesignStructure(this);
        structure_->SetFilters(list[i], this->optInfoNode);
      }
    }
    else
    {
      if(simp_pn->Has("regularization"))
        throw Exception("regularization not implemented");
    }
  }
  
  if(harmonic) mechRHS.Init<complex<double> >(design, PRESSURE); // in many cases NULL;
          else mechRHS.Init<double>(design, PRESSURE);

  ErsatzMaterial::PostInit();
}

void SIMP::SetElementK(DesignElement* de, const TransferFunction* tf, Application app, DenseMatrix* out, CalcMode calcMode, bool derivative)
{
  if(harmonic){
    if(pde->HasComplexMatData(de->elem->regionId))
      SetElementK<Complex, Complex >(de, tf, app, out, calcMode, derivative);
    else
      SetElementK<Complex, double >(de, tf, app, out, calcMode, derivative);
  }
  else SetElementK<double,double>(de, tf, app, out, calcMode, derivative);
}

template <class T1, class T2>
void SIMP::SetElementK(DesignElement* de, const TransferFunction* tf, Application app, DenseMatrix* mat_out, CalcMode calcMode, bool derivative)
{
  Matrix<T1>& out = dynamic_cast<Matrix<T1>& >(*mat_out);


  switch(app)
  {
  case MECH:
  case ACOUSTIC:
  {
    int mm = de->multimaterial != NULL ? de->multimaterial->index : -1;

    const Matrix<T2>& stiffness = dynamic_cast<Matrix<T2>& >(material->Stiffness(de->elem, false, mm)); // no bimaterial

    // Find the transfer function for K (e.g. DENSITY, MECH)
    T1 k_factor = derivative ? tf->Derivative(de, DesignElement::SMART) : tf->Transform(de, DesignElement::SMART);

    // copy from real mechStiffness to potential complex out and factor the derivative
    Assign(out, stiffness, k_factor);
    // This log is very expensive, it blows up inv_tensor in the debug mode
    LOG_DBG3(simp) << "SetElementK: el=" << de->elem->elemNum << " di=" << de->GetIndex() << " mm=" << mm << " K_org=" <<  stiffness.ToString() << " k_factor " << k_factor << " -> " << out.ToString();

    if(design->GetRegion(de->elem->regionId)->HasBiMaterial())
    {
      const Matrix<T2>& bimat = dynamic_cast<Matrix<T2>& >(material->Stiffness(de->elem, true)); // yes, bimaterial
      // rho^3 * E1 + (1-rho^3) * E2, in the derivative case 3*rho^2 * E1 - 3*rho^2 * E2
      k_factor = !derivative ? 1.0 - k_factor : -1.0 *  k_factor;
      Add(out, k_factor, bimat);
      // LOG_DBG3(simp) << "SetElementK: K_bi_org=" <<  bimat.ToString() << " k_factor " << k_factor << " -> " << out.ToString();
    }

    if(harmonic)
    {
      tf = design->GetTransferFunction(de->GetType(), MASS);
      AddMassToStiffness(tf, de, dynamic_cast<Matrix<complex<double> >& >(out), derivative, false); // no bimaterial

      // LOG_DBG3(simp) << "SetElementK: m_factor " << m_factor << " -> " << out.ToString();

      if(design->GetRegion(de->elem->regionId)->HasBiMaterial())
      {
        // rho^3 * E1 + (1-rho^3) * E2, in the derivative case 3*rho^2 * E1 - 3*rho^2 * E2
        AddMassToStiffness(tf, de, dynamic_cast<Matrix<complex<double> >& >(out), derivative, true); // bimaterial

        // LOG_DBG3(simp) << "SetElementK: m_bi_factor " << m_factor << " -> " << out.ToString();
      }
    }
    break;
  }

  case ELEC:
  {
    Matrix<std::complex<double> >& stiffness = dynamic_cast<ElecMat *>(material)->ElecStiffness(de->elem, false); // no bimaterial

    // Find the transfer function for K (e.g. DENSITY, MECH)
    T1 k_factor = derivative ? tf->Derivative(de, DesignElement::SMART) : tf->Transform(de, DesignElement::SMART);

    // copy from ElecStiffness to out and factor the derivative
    if (harmonic)
      Assign(out, dynamic_cast<Matrix<T1>& >(stiffness), k_factor);
    else
      Assign(out, stiffness.GetPart(Global::REAL), k_factor);

    // This log is very expensive, it blows up inv_tensor in the debug mode
    // LOG_DBG3(simp) << "SetElementK: K_org=" <<  stiffness.ToString() << " k_factor " << k_factor << " -> " << out.ToString();

    if(design->GetRegion(de->elem->regionId)->HasBiMaterial())
    {
      Matrix<std::complex<double> >& bimat = dynamic_cast<ElecMat *>(material)->ElecStiffness(de->elem, true); // yes, bimaterial
      // rho^3 * E1 + (1-rho^3) * E2, in the derivative case 3*rho^2 * E1 - 3*rho^2 * E2
      k_factor = derivative ? tf->Derivative(de, DesignElement::SMART, true) : tf->Transform(de, DesignElement::SMART, -13.456, true);
      if (harmonic)
        Add(out, k_factor, dynamic_cast<Matrix<T1>& >(bimat));
      else
        Add(out, k_factor, bimat.GetPart(Global::REAL));
      // LOG_DBG3(simp) << "SetElementK: K_bi_org=" <<  bimat.ToString() << " k_factor " << k_factor << " -> " << out.ToString();
    }
    break;
  }

  default:
    assert(false); // other cases should be handled in PiezoSIMP
  } // end switch
}


void SIMP::AddMassToStiffness(const TransferFunction* mtf, DesignElement* de, Matrix<complex<double> >& K_in_S_out, bool derivative, bool bimaterial)
{
  // The result matrix is
  // S = K + i*omega*C - omega^2*M
  // with purely imaginary C = alpha_k*K+alpha*M
  // S = K + i*omega*alpha_k*K + i*omega*alpha_m*M - omega^2*M
  // with m_factor
  // S = K + i*omega*alpha_k*K + i*omega*alpha_m*m_factor*M - omega^2*m_factor*M
  // With S = K in the beginning this is
  // S += i*alpha_k*S + (i*alpha_m*m_factor-omega^2*m_factor)*M
  //
  // in case we have pamping (e.g. Sigmund; Morhology; 2007) there is to add
  // j*omega*pamping*rho*(1-rho)*M and for the derivative case
  // j*omega*pamping*rho'*M - j*2*omega*pamping*rho*rho'*M = j*omega*pamping*rho'(1-2*rho)

  double mtv =  mtf->Transform(de, DesignElement::SMART);
  double mdv =  mtf->Derivative(de, DesignElement::SMART);

  double m_factor = derivative ? mdv : mtv;
  if(bimaterial)  // rho^3 * E1 + (1-rho^3) * E2, in the derivative case 3*rho^2 * E1 - 3*rho^2 * E2
    m_factor = !derivative ? 1.0 - m_factor : -1.0 *  m_factor;

  // change name only
  Matrix<complex<double> >& S = K_in_S_out;

  // multimaterial stuff
  int index = de->multimaterial != NULL ? de->multimaterial->index : -1;

  const Matrix<double>& M = material->Mass(de->elem, bimaterial, index);
  assert(S.GetNumRows() == M.GetNumRows() && S.GetNumCols() == M.GetNumCols());

  // find alpha, beta and omega
  double omega = 2.0 * M_PI * pde->GetSolveStep()->GetActFreq() ;  // todo: check with multiple excitation frequencies!
  double alpha_k = 0.0;
  double alpha_m  = 0.0;
  double pamping_m = 0.0; // add on without omega

  assert(false);
  /* FIXME
  // do we have damping (C = alpha*M+beta*K) -> this is pure imaginary!
  RegionIdType regionId = de->elem->regionId;
  
  if(pde->GetDamping(regionId) == RAYLEIGH)
  {
    // the alpha and beta might be calculated and adjusted, get them
    // from the integrators in the form as they are used for the state problem!
    alpha_k = assemble_->GetBiLinForm(regionId, pde, pde, "LinElastInt")->EvalSecMatFac();

    // now alpha_m
    alpha_m = assemble_->GetBiLinForm(regionId, pde, pde, "MassInt")->EvalSecMatFac();

    assert(omega > 0);

    // pamping stuff without omega
    double pamping = design->GetPampingValue(); // 0 if not applicable
    if(!derivative)
      pamping_m = pamping * mtv * (1.0 - mtv);
    else // pamping*rho'(1-2*rho)
      pamping_m = pamping * mdv * (1.0 - 2.0 * mtv);
  }

  */

	const unsigned int srows(S.GetNumRows());
	const unsigned int scols(S.GetNumCols());
  // we first add the K part of C (= pure imaginary)
  for(unsigned int r = 0; r < srows; r++)
    for(unsigned int c = 0; c < scols; c++)
      S[r][c] = complex<double>(S[r][c].real(), omega * alpha_k * S[r][c].real());

  // we the add the M part of C and the real mass part
  complex<double> damp_mass = complex<double>(-1.0 *omega*omega*m_factor, omega*(alpha_m*m_factor  + pamping_m));
  for(unsigned int r = 0; r < srows; r++)
    for(unsigned int c = 0; c < scols; c++)
      S[r][c] += damp_mass * M[r][c];


  LOG_DBG2(simp) << "AddMassToStiffness: d=" << de->elem->elemNum << " der=" << derivative << " bm=" << bimaterial
                 << " m_factor:" << m_factor << " alpha_k: " << alpha_k << " alpha_m: " << alpha_m << " pamping_m:" << pamping_m
                 << " omega: " << omega << " K_img: " << (omega * alpha_k) << " damp_mass: " << damp_mass << " M=" << M.ToString();
}


double SIMP::CalcFunction(Excitation& excite, Function* f, bool derivative)
{
  // this app is for the PDE
  Application app = ToApp(pde);

  if(!derivative)
    return ErsatzMaterial::CalcFunction(excite, f, derivative);

  // only special derivatives, the rest also EM
  switch(f->GetType())
  {
  case Function::STRESS:
  case Function::STRESS_DENSITY:
  {
    TransferFunction* tf = design->GetTransferFunction(DesignElement::Default(pde), TransferFunction::Default(pde), true);
    CalcVonMisesStressGradient(excite, f, tf);
    break;
  }

  case Function::GLOBAL_DYNAMIC_COMPLIANCE:
  case Function::OUTPUT:
  case Function::DYNAMIC_OUTPUT:
  case Function::CONJUGATE_COMPLIANCE:
  case Function::ABS_OUTPUT:
  {
    // synthesis of compliant mechanism: As our adjoint PDE
    // c' = l K' u
    TransferFunction* tf = design->GetTransferFunction(DesignElement::Default(pde), TransferFunction::Default(pde), true, true); // excpetion and use_single
    double weight = excite.GetWeightedFactor(f);
    LOG_DBG(simp) << "CalcFunction(idx=" << excite.index << ") norm_weight= " <<  excite.normalized_weight  << " factor=" << excite.GetFactor(f) << " weight=" << weight;
    CalcU1KU2(tf, adjoint.Get(excite, f)->elem[app], app, forward.Get(excite)->elem[app], NULL, weight, STANDARD, f);
    break;
  }

  default:
    return ErsatzMaterial::CalcFunction(excite, f, derivative);
  }

  return 0.0; // only derivatives evaluated
}


void SIMP::CalcVonMisesStressGradient(Excitation& excite, Function* f, TransferFunction* tf)
{
	// see comment in ErsatzMaterial::CalcVonMisesStressVector()! it's tricky stuff :(
  // For the function we pack the stuff in Function::Local, for the gradient we do it here as the computation are too far
  // away from the other local gradient computations.
  //
	// the gradient is lambda^T * ( K' * u - f')  + alpha * 2 * stress^T * M * (rho^p)' * E_0 * B * u
  //
  // that means, if the stress constraint region is not a design region, we don't add something. But take
  // care, there might also be several stress constraints for a set of design regions.
  //
  // we do NOT weight!

  for(unsigned int i = 0; i < design->data.GetSize(); i++)
    LOG_DBG2(simp) << "CVMSG: f=" << f->ToString(this->me) << " de=" << design->data[i].elem->elemNum << " org=" << design->data[i].GetPlainGradient(f);

  // alpha is from the globalization which is in the form sum max(0, g_i-c)^p and alpha is p*max(0, g_i-c)^(p-1) where g_i is the vonMisesStress
  Vector<double> alpha;
  // 2 * stress^T * M * (rho^p)' * E_0 * B * u
  Vector<double> appendix;

  if(harmonic)
  {
    StressConstraint<complex<double> > sc(&excite, f, this, &forward);
    sc.CalcGlobalizationFactor(alpha);
    sc.CalcGradStresses(appendix);
  }
  else
  {
    StressConstraint<double> sc(&excite, f, this, &forward);
    sc.CalcGlobalizationFactor(alpha);
    sc.CalcGradStresses(appendix);
  }
  assert(appendix.GetSize() == alpha.GetSize());

  DesignDependentRHS rhs;
  rhs.Init<double>(Optimization::STRESS, excite.label);
  // calc lambda^T *  K' * u -> this already stores the results by AddGradient()!
  CalcU1KU2(tf, adjoint.Get(excite, f)->elem[MECH], MECH, forward.Get(excite)->elem[MECH], &rhs, 1.0, STANDARD, f);

  // add the appendix stuff
  for(unsigned int i = 0; i < design->data.GetSize(); i++)
  {
    DesignElement& de = design->data[i];
    // Three cases:
    // a) stress is defined on whole design domain (one or more regions)
    // b) stress region is defined on one of two or more design regions
    // c) stress region is not in any of the design regions.

    int idx = -1; // case c) an sometimes b) never a)

    // case a)
    if(f->region == ALL_REGIONS || design->regions[0].GetSize() == 1)
    {
      assert(de.elem->elemNum == f->elements[i]->elem->elemNum);
      idx = i;
    }
    // case b)
    if(idx != -1 && design->Contains(f->region))
    {
      // we are at a design element and have to find it within stress
      // TODO make it faster
      for(unsigned int e = 0; idx != -1 && e < f->elements.GetSize(); e++)
        if(de.elem->elemNum == f->elements[e]->elem->elemNum)
          idx = e;
    }
    // case c): idx is already -1
    if(idx != -1)
      de.AddGradient(f, alpha[idx] * appendix[idx]);
    LOG_DBG2(simp) << "CVMSG: f=" << f->ToString(this->me) << " de=" << de.elem->elemNum << " idx=" << idx << " alpha="
                   << (idx != -1 ? alpha[i] : -1.0)  << "* app=" << (idx != -1 ? appendix[i] : -1.0) << " -> " << de.GetPlainGradient(f);
	}
}

DesignDependentRHS::DesignDependentRHS()
{
  valid       = false;
  app         = Optimization::NO_APP;
  vec         = NULL;
  elem        = NULL;
  test_strain = MechPDE::NOT_SET;
}

DesignDependentRHS::~DesignDependentRHS()
{
  valid = false;
  if(vec != NULL) { delete vec; vec = NULL; }
}

template <class T>
bool DesignDependentRHS::Init(DesignSpace* design, Optimization::Application app)
{
  assert(app == Optimization::CHARGE_DENSITY || app == Optimization::PRESSURE);
  std::string name = app == Optimization::CHARGE_DENSITY ? "LinNeumannInt" : "PressureLinForm";


  // check if we have a form with the application name
  LinearForm* form = NULL;
  LinearFormContext* actContext = NULL;

  SinglePDE* mech = domain->GetSinglePDE("mechanic", false);
  if(mech == NULL) return false; // wrong pde -> extend if you need it!

  StdVector<LinearFormContext*>* forms = &(mech->GetAssemble()->GetLinForms());

  for(StdVector<LinearFormContext*>::iterator it = forms->Begin(); it != forms->End(); it++)
  {
    // get integrator
    actContext = *it;
    if(actContext->GetIntegrator()->GetName() == name)
    {
      assert(false);
      if(form != NULL) EXCEPTION("linear surface form '" << name << "' not unique");
      form = actContext->GetIntegrator(); // form = dynamic_cast<LinearSurfForm*>(actContext->GetIntegrator());
    }
  }

  LOG_DBG(simp) << "DesignDependentRHS::Init(app = " << Optimization::application.ToString(app) << ") -> form = "
               << (form != NULL ? form->GetName() : "NULL");

  // form is not necessary defined in the xml file!
  if(form == NULL) return false; // no form, no RHS!

  // the context knows the surface elements!
  this->valid = true;
  this->app   = app;
  EntityIterator eit = actContext->GetEntities()->GetIterator();
  elem = eit.GetSurfElem();

  // FIXME
  assert(false);
  // calculate the rhs for the reference element, first store and then extract all but one node
  design->DisableTransferFunctions();
  // FIXME form->SetSurfElem(const_cast<SurfElem*>(elem)); // set the internal actElem_ of the form
  Vector<T> full;
  // FIXME form->CalcElemVector(full,const_cast<EntityIterator&>(eit));
  // enable again our transfer functions
  design->EnableTransferFunctions();


  // check the number of dofs, copy the the first nodes dofs, check that all other are the same
  assert(full.GetSize() >= elem->connect.GetSize());
  int dof = full.GetSize() / elem->connect.GetSize();
  if(dof < 0 || dof > 3)
    EXCEPTION("Surface element of " << elem->connect.GetSize() << " nodes has RHS with "
              << full.GetSize() << " entries");
  assert(vec == NULL);
  // copy the first nodes dofs
  Vector<T>* vt = new Vector<T>(dof);
  // vec is the base SingleVector which has no abstract operators overloaded
  vec = vt;
  for(int i = 0; i < dof; i++) (*vt)[i] = full[i];

  // check
  for(unsigned int n = 0; n < full.GetSize()/dof; n += dof)
    for(int d = 0; d < dof; d++)
      if(!close(full[n+d],(*vt)[d]))
        EXCEPTION("RHS values are not the same for each node: " << full.ToString());

  // store all node numbers in the sorted set
  // do at the end such iterator is valid for CalcElemVector()
  eit.Begin();
  while(!eit.IsEnd())
  {
    StdVector<unsigned int> elem_nodes = eit.GetElem()->connect;
    for(unsigned int n = 0; n < elem_nodes.GetSize(); n++)
      nodes.insert(elem_nodes[n]);
    eit++;
  }

  LOG_DBG(simp) << "DesignDependentRHS::Init -> " << ToString(1);
  LOG_DBG2(simp) << "DesignDependentRHS::Init -> " << ToString(0);

  return true;
}


template <class T>
bool DesignDependentRHS::Init(Optimization::Application app, std::string excite_label)
{
  assert(app == Optimization::STRESS);
  this->app = app;
  this->test_strain = MechPDE::testStrain.IsValid(excite_label) ? MechPDE::testStrain.Parse(excite_label) : MechPDE::NOT_SET;
  return true;
}


std::string DesignDependentRHS::ToString(int level)
{
  std::ostringstream os;
  os << "valid=" << valid;
  if(!valid) return os.str();

  os << " vec=" << vec->ToString();
  os << " elem=" << elem->ToString();

  os << " nodes=";
  if(level == 0)
  {
    for(std::set<unsigned int>::iterator i = nodes.begin(); i != nodes.end(); i++)
      os << *i << ", ";
  }
  else
  {
    os << "#" << nodes.size();
  }
  return os.str();
}

  // Explicit template instantiation
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
template bool DesignDependentRHS::Init<double>(DesignSpace* design, Optimization::Application app);
template bool DesignDependentRHS::Init<complex<double> >(DesignSpace* design, Optimization::Application app);
#endif
