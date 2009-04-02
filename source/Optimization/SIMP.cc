#include "Optimization/SIMP.hh"
#include "Optimization/DesignSpace.hh"
#include "Optimization/DesignElement.hh"
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
  ParamNode* simp_pn = pn->Get("SIMP", false);

  // There might be a filter regularization based on the design element.
  if(simp_pn != NULL && simp_pn->Has("regularization", "type", "filter"))
  {
    StdVector<ParamNode*> list = simp_pn->Get("regularization")->GetList("filter");
    // this is save for design=polarization
    for(unsigned int i = 0; i < list.GetSize(); i++)
      SIMPElement::InitFilter(design->data, list[i]);
  }
  else
  {
    if(simp_pn != NULL && simp_pn->Has("regularization"))
      throw Exception("regularization not implemented");
  }

  mech_mat_ = NULL; // set in PostInit()
}

SIMP::~SIMP()
{
  if(mechRHS.vec != NULL)       { delete mechRHS.vec;  mechRHS.vec = NULL;  }
}

void SIMP::PostInit()
{
  if(harmonic) mechRHS.Init<complex<double> >(design, PRESSURE); // in many cases NULL;
          else mechRHS.Init<double>(design, PRESSURE);

  ErsatzMaterial::PostInit();
  
  // only after ErsatzMaterial::PostInit();
  mech_mat_ = dynamic_cast<OptMechMat*>(material);
  assert(mech_mat_ != NULL);
}

void SIMP::SetElementK(DesignElement* de, Application app, DenseMatrix* out, CalcMode calcMode)
{
  if(harmonic) SetElementK<std::complex<double> >(de, app, out, calcMode);
  else SetElementK<double>(de, app, out, calcMode);
}

template <class T>
void SIMP::SetElementK(DesignElement* de, Application app, DenseMatrix* mat_out, CalcMode calcMode)
{
  Matrix<T>& out = dynamic_cast<Matrix<T>& >(*mat_out);

  switch(app)
  {
  case MECH:
  {
    const Matrix<double> mechStiffness = mech_mat_->MechStiffness(de->elem);
    out.Resize(mechStiffness.GetNumRows(), mechStiffness.GetNumCols());

    // Find the transferfunction for K (e.g. DENSITY, MECH)
    TransferFunction* tf = design->GetTransferFunction(de->GetType(), app);
    double k_factor = tf->Derivative(de);

    // copy from real mechStiffness to potential complex out and factor the derivative
    Assign(out, mechStiffness, k_factor);
    LOG_DBG3(simp) << "SetElementK: org mech " << out.ToString(0);

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
  double omega = 2.0 * M_PI * mech->GetSolveStep()->GetActFreq() ;  // todo: check with multiple excitation frequencies!
  double alpha_k = 0.0;
  double alpha_m  = 0.0;

  // do we have damping (C = alpha*M+beta*K) -> this is pure imaginary!
  RegionIdType regionId = de->elem->regionId;
  if(mech->GetDampingList().find(regionId) != mech->GetDampingList().end()
      && mech->GetDampingList().find(regionId)->second == RAYLEIGH)
  {
    // we need a math parser (is an unsigned int)
    unsigned int handle = domain->GetMathParser()->GetNewHandle();

    // the alpha and beta might be calculated and adjusted, get them
    // from the integrators in the form as they are used for the state problem!
    std::string txt = assemble_->GetBiLinForm(regionId, mech, mech, "linElastInt")->GetSecMatFac();
    domain->GetMathParser()->SetExpr(handle, txt);
    alpha_k = domain->GetMathParser()->Eval(handle);

    // now alpha_m
    txt = assemble_->GetBiLinForm(regionId, mech, mech, "MassInt")->GetSecMatFac();
    domain->GetMathParser()->SetExpr(handle, txt);
    alpha_m = domain->GetMathParser()->Eval(handle);

    domain->GetMathParser()->ReleaseHandle(handle);
    assert(omega > 0 && alpha_k > 0 && alpha_m > 0);
  }

  // we first add the K part of C (= pure imaginary)
  for(unsigned int r = 0; r < S.GetNumRows(); r++)
    for(unsigned int c = 0; c < S.GetNumCols(); c++)
      S[r][c] = complex<double>(S[r][c].real(), omega * alpha_k * S[r][c].real());

  // we the add the M part of C and the real mass part
  complex<double> damp_mass = complex<double>(-1.0 *omega*omega*m_factor, omega*alpha_m*m_factor);
  for(unsigned int r = 0; r < S.GetNumRows(); r++)
    for(unsigned int c = 0; c < S.GetNumCols(); c++)
      S[r][c] += damp_mass * M[r][c];

  LOG_DBG2(simp) << "AddMassToStiffness: m_factor:" << m_factor << " alpha_k: " << alpha_k << " alpha_m: " << alpha_m
                 << " omega: " << omega << " K_img: " << (omega * alpha_k) << " damp_mass: " << damp_mass;
}

void SIMP::CalcObjectiveGradient(Excitation& excite)
{
  TransferFunction* tf = design->GetTransferFunction(DesignElement::DENSITY, MECH, true);

  int idx = excite.index;
  double weight = excite.GetWeightedFactor(cost);
  LOG_DBG(simp) << "CalcObjectiveGradient(idx=" << excite.index << ") norm_weight= " <<  excite.normalized_weight
                << " factor=" << excite.GetFactor(cost) << " weight=" << weight;

  switch(cost->type)
  {
  case GLOBAL_DYNAMIC_COMPLIANCE:
    // synthesis of compliant mechanism: As our adjoint PDE
    // c' = l K' u
    CalcU1KU2(tf, adjoint.data[idx]->elem[MECH], MECH, forward.data[idx]->elem[MECH], NULL, weight);
    break;

  case OUTPUT:
  case DYNAMIC_OUTPUT:
  case CONJUGATE_COMPLIANCE:
    // synthesis of compliant mechanism: As our adjoint PDE
    // c' = l K' u
    CalcU1KU2(tf, adjoint.data[idx]->elem[MECH], MECH, forward.data[idx]->elem[MECH], NULL, weight);
    break;

  default:
    ErsatzMaterial::CalcObjectiveGradient(excite);
  }
}

SurfaceRef::SurfaceRef()
{
  valid = false;
  app   = Optimization::NO_APP;
  vec = NULL;
}

SurfaceRef::~SurfaceRef()
{
  valid = false;
  if(vec != NULL) { delete vec; vec = NULL; }
}

template <class T>
bool SurfaceRef::Init(DesignSpace* design, Optimization::Application app)
{
  assert(app == Optimization::CHARGE_DENSITY || app == Optimization::PRESSURE);
  std::string name = app == Optimization::CHARGE_DENSITY ? "LinNeumannInt" : "PressureLinForm";

  // check if we have a form with the application name
  LinearSurfForm* form = NULL;
  LinearFormContext* actContext = NULL;
  std::set<LinearFormContext*>* forms =
    domain->GetSinglePDE("mechanic")->getPDE_assemble()->GetLinForms();

  for(std::set<LinearFormContext*>::iterator it = forms->begin(); it != forms->end(); it++)
  {
    // get integrator
    actContext = *it;
    if(actContext->GetIntegrator()->GetName() == name)
    {
      if(form != NULL) EXCEPTION("linear surface form '" << name << "' not unique");
      form = dynamic_cast<LinearSurfForm*>(actContext->GetIntegrator());
    }
  }

  LOG_DBG(simp) << "SurfaceRef::Init(app = " << Optimization::application.ToString(app) << ") -> form = "
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
  // vec is the base SingleVector which has no abstact operators overloaded
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

  LOG_DBG(simp) << "SurfaceRef::Init -> " << ToString(1);
  LOG_DBG2(simp) << "SurfaceRef::Init -> " << ToString(0);

  return true;
}


std::string SurfaceRef::ToString(int level)
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

