#include "Optimization/SIMP.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/DesignStructure.hh"
#include "Optimization/OptimizationMaterial.hh"
#include "Domain/domain.hh"
#include "Domain/surfElem.hh"
#include "General/exception.hh"
#include "Utils/basenodestoresol.hh"
#include "PDE/SinglePDE.hh"
#include "PDE/mechPDE.hh"
#include "Forms/baseForm.hh"
#include "Forms/linSurfForm.hh"
#include "CoupledPDE/DirectCoupledPDE.hh"
#include "Utils/StdVector.hh"
#include "MatVec/vector.hh"
#include "Utils/result.hh"
#include "Utils/mathParser/mathParser.hh"
#include "Driver/assemble.hh"
#include "Driver/baseSolveStep.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "OLAS/algsys/baseentrymanipulator.hh"
#include "OLAS/algsys/basesystem.hh"
#include "MatVec/basematrix.hh"
#include "MatVec/stdmatrix.hh"

#include <string>

using namespace CoupledField;

using std::complex;

DECLARE_LOG(conditions)

DECLARE_LOG(simp)
DEFINE_LOG(simp, "simp")


SIMP::SIMP() : ErsatzMaterial()
{
  mech_mat_ = NULL; // set in PostInit()
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

    // check for bimaterial and read tensor if available
    if(simp_pn->Has("bimaterial"))
    {
      Matrix<double> t(3, 3);
      bool ok = Function::ReadTensor(simp_pn->Get("bimaterial"), t);
      if(!ok) EXCEPTION("bimaterial specified but no tensor given or incorrect format");

      design->SetBiMatTensor(t);
      LOG_DBG3(simp) << "bimaterial tensor = " << std::endl << design->GetBiMatTensor().ToString(1);
    }
  }
  
  if(harmonic) mechRHS.Init<complex<double> >(design, PRESSURE); // in many cases NULL;
          else mechRHS.Init<double>(design, PRESSURE);

  ErsatzMaterial::PostInit();
  
  // only after ErsatzMaterial::PostInit();
  mech_mat_ = dynamic_cast<OptMechMat*>(material);
  assert(mech_mat_ != NULL);
}

void SIMP::SetElementK(DesignElement* de, Application app, DenseMatrix* out, CalcMode calcMode, bool derivative)
{
  if(harmonic) SetElementK<std::complex<double> >(de, app, out, calcMode, derivative);
  else SetElementK<double>(de, app, out, calcMode, derivative);
}

template <class T>
void SIMP::SetElementK(DesignElement* de, Application app, DenseMatrix* mat_out, CalcMode calcMode, bool derivative)
{
  Matrix<T>& out = dynamic_cast<Matrix<T>& >(*mat_out);

  switch(app)
  {
  case MECH:
  {
    const Matrix<double> &mechStiffness = mech_mat_->MechStiffness(de->elem);
    
    // Find the transfer function for K (e.g. DENSITY, MECH)
    TransferFunction* tf = design->GetTransferFunction(de->GetType(), app);
    double k_factor = derivative ? tf->Derivative(de, DesignElement::SMART) : tf->Transform(de, DesignElement::SMART);

    // copy from real mechStiffness to potential complex out and factor the derivative
    Assign(out, mechStiffness, k_factor);
    // This log is very expensive, it blows up inv_tensor in the debug mode
    //LOG_DBG3(simp) << "SetElementK: org mech " << out.ToString(0);

    if(harmonic)
    {
      tf = design->GetTransferFunction(de->GetType(), MASS);
      AddMassToStiffness(tf->Derivative(de), de, dynamic_cast<Matrix<complex<double> >& >(out));
    }
    break;
  }
  default:
    assert(false); // other cases should be handled in PiezoSIMP
  } // end switch
}


void SIMP::AddMassToStiffness(double m_factor, DesignElement* de, Matrix<complex<double> >& K_in_S_out)
{
  // The result matrix is
  // S = K + i*omega*C - omega^2*M
  // with purly imaginary C = alpha_k*K+alpha*M
  // S = K + i*omega*alpha_k*K + i*omega*alpha_m*M - omega^2*M
  // with m_factor
  // S = K + i*omega*alpha_k*K + i*omega*alpha_m*m_factor*M - omega^2*m_factor*M
  // With S = K in the beginning this is
  // S += i*alpha_k*S + (i*alpha_m*m_factor-omega^2*m_factor)*M

  // change name only
  Matrix<complex<double> >& S = K_in_S_out;
  const Matrix<double>& M = mech_mat_->MechMass(de->elem);
  assert(S.GetNumRows() == M.GetNumRows() && S.GetNumCols() == M.GetNumCols());

  // find alpha, beta and omega
  double omega = 2.0 * M_PI * pde->GetSolveStep()->GetActFreq() ;  // todo: check with multiple excitation frequencies!
  double alpha_k = 0.0;
  double alpha_m  = 0.0;

  // do we have damping (C = alpha*M+beta*K) -> this is pure imaginary!
  RegionIdType regionId = de->elem->regionId;
  
  if(pde->GetDamping(regionId) == RAYLEIGH)
  {
    // we need a math parser (is an unsigned int)
    unsigned int handle = domain->GetMathParser()->GetNewHandle();

    // the alpha and beta might be calculated and adjusted, get them
    // from the integrators in the form as they are used for the state problem!
    std::string txt = assemble_->GetBiLinForm(regionId, pde, pde, "linElastInt")->GetSecMatFac();
    domain->GetMathParser()->SetExpr(handle, txt);
    alpha_k = domain->GetMathParser()->Eval(handle);

    // now alpha_m
    txt = assemble_->GetBiLinForm(regionId, pde, pde, "MassInt")->GetSecMatFac();
    domain->GetMathParser()->SetExpr(handle, txt);
    alpha_m = domain->GetMathParser()->Eval(handle);

    domain->GetMathParser()->ReleaseHandle(handle);
    assert(omega > 0 && alpha_k > 0 && alpha_m > 0);
  }

	const unsigned int srows(S.GetNumRows());
	const unsigned int scols(S.GetNumCols());
  // we first add the K part of C (= pure imaginary)
  for(unsigned int r = 0; r < srows; r++)
    for(unsigned int c = 0; c < scols; c++)
      S[r][c] = complex<double>(S[r][c].real(), omega * alpha_k * S[r][c].real());

  // we the add the M part of C and the real mass part
  complex<double> damp_mass = complex<double>(-1.0 *omega*omega*m_factor, omega*alpha_m*m_factor);
  for(unsigned int r = 0; r < srows; r++)
    for(unsigned int c = 0; c < scols; c++)
      S[r][c] += damp_mass * M[r][c];

  LOG_DBG2(simp) << "AddMassToStiffness: m_factor:" << m_factor << " alpha_k: " << alpha_k << " alpha_m: " << alpha_m
                 << " omega: " << omega << " K_img: " << (omega * alpha_k) << " damp_mass: " << damp_mass;
}


double SIMP::CalcFunction(Excitation& excite, Function* f, bool derivative)
{
  // this implements only the gradients of some functions
  if(!derivative)
    return ErsatzMaterial::CalcFunction(excite, f, derivative);

  TransferFunction* tf = design->GetTransferFunction(DesignElement::DENSITY, MECH, true);
  double weight = excite.GetWeightedFactor(f);
  LOG_DBG(simp) << "CalcFunction(idx=" << excite.index << ") norm_weight= " <<  excite.normalized_weight
                << " factor=" << excite.GetFactor(f) << " weight=" << weight;

  // only special derivatives, the rest also EM
  switch(f->GetType())
  {
  case Function::STRESS:
  {
    CalcVonMisesStressGradient(excite, f, tf);
    break;
  }


  case Function::GLOBAL_DYNAMIC_COMPLIANCE:
  case Function::OUTPUT:
  case Function::DYNAMIC_OUTPUT:
  case Function::CONJUGATE_COMPLIANCE:
  case Function::ABS_DYN_OUTPUT_SQUARED:
    // synthesis of compliant mechanism: As our adjoint PDE
    // c' = l K' u
    CalcU1KU2(tf, adjoint.Get(excite, f)->elem[MECH], MECH, forward.Get(excite)->elem[MECH], NULL, weight, STANDARD, f);
    break;

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
  // we do NOT weight!

  for(unsigned int i = 0; i < design->data.GetSize(); i++)
    LOG_DBG2(simp) << "CVMSG: f=" << f->ToString(this->me) << " de=" << design->data[i].elem->elemNum << " org=" << design->data[i].GetPlainGradient(f);

  // alpha is from the globalization which is in the form sum max(0, g_i-c)^p and alpha is p*max(0, g_i-c)^(p-1) where g_i is the vonMisesStress
  Vector<double> alpha = CalcVonMisesStressGlobalizationFactor(excite, f);
  assert(alpha.GetSize() == design->data.GetSize());

  // 2 * stress^T * M * (rho^p)' * E_0 * B * u can be obtained with a special attributes
  Vector<double> appendix = CalcVonMisesStressVector(excite, f, false, true);
  assert(appendix.GetSize() == alpha.GetSize());

  DesignDependentRHS rhs;
  rhs.Init<double>(Optimization::STRESS, MechPDE::testStrain.Parse(excite.label));
  // calc lambda^T *  K' * u -> this already stores the results by AddGradient()!
  CalcU1KU2(tf, adjoint.Get(excite, f)->elem[MECH], MECH, forward.Get(excite)->elem[MECH], &rhs, 1.0, STANDARD, f);

  // add the appendix stuff
  for(unsigned int i = 0; i < design->data.GetSize(); i++)
	{
    DesignElement& de = design->data[i];
		de.AddGradient(f, alpha[i] * appendix[i]);
		LOG_DBG2(simp) << "CVMSG: f=" << f->ToString(this->me) << " de=" << de.elem->elemNum << " alpha=" << alpha[i] << "* app=" << appendix[i] << " ="
		               << alpha[i] * appendix[i] << " -> " << de.GetPlainGradient(f);
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
  LinearSurfForm* form = NULL;
  LinearFormContext* actContext = NULL;
  StdVector<LinearFormContext*>* forms =
    &(domain->GetSinglePDE("mechanic")->getPDE_assemble()->GetLinForms());

  for(StdVector<LinearFormContext*>::iterator it = forms->Begin(); it != forms->End(); it++)
  {
    // get integrator
    actContext = *it;
    if(actContext->GetIntegrator()->GetName() == name)
    {
      if(form != NULL) EXCEPTION("linear surface form '" << name << "' not unique");
      form = dynamic_cast<LinearSurfForm*>(actContext->GetIntegrator());
    }
  }

  LOG_DBG(simp) << "DesignDependentRHS::Init(app = " << Optimization::application.ToString(app) << ") -> form = "
                  << (form != NULL ? form->GetName() : "NULL");

  // form is not necessary defined int the xml file!
  if(form == NULL) return false; // no form, no RHS!

  // the context knows the surface elements!
  this->valid = true;
  this->app   = app;
  EntityIterator eit = actContext->GetEntities()->GetIterator();
  elem = eit.GetSurfElem();

  // calculate the rhs for the reference element, first store and then extract all but one node
  design->DisableTransferFunctions();
  form->SetSurfElem(const_cast<SurfElem*>(elem)); // set the internal actElem_ of the form
  Vector<T> full;
  form->CalcElemVector(full,const_cast<EntityIterator&>(eit));
  // enable again our tranfer functions
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
  // do at the end such iterater is valid for CalcElemVector()
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
bool DesignDependentRHS::Init(Optimization::Application app, MechPDE::TestStrain test_strain)
{
  assert(app == Optimization::STRESS);
  this->app = app;
  this->test_strain = test_strain;
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
