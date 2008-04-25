#include "Optimization/SIMP.hh"
#include "Optimization/DesignSpace.hh"
#include "Optimization/DesignElement.hh"
#include "Domain/domain.hh"
#include "Domain/surfElem.hh"
#include "General/exception.hh"
#include "Utils/basenodestoresol.hh"
#include "PDE/SinglePDE.hh"
#include "PDE/mechPDE.hh"
#include "Forms/baseForm.hh"
#include "Forms/linSurfForm.hh"
#include "Forms/SurfaceNormalInt.hh"
#include "CoupledPDE/DirectCoupledPDE.hh"
#include "Utils/StdVector.hh"
#include "Utils/vector.hh"
#include "Utils/result.hh"
#include "Utils/mathParser/mathParser.hh"
#include "Driver/assemble.hh"
#include "Driver/baseSolveStep.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "OLAS/olas.hh"
#include "OLAS/algsys/baseentrymanipulator.hh"
#include "OLAS/matvec/basematrix.hh"
#include "OLAS/matvec/stdmatrix.hh"

#include <string>

using namespace CoupledField;

using OLAS::StdMatrix;
using std::complex;

DECLARE_LOG(conditions)

DECLARE_LOG(simp)
DEFINE_LOG(simp, "simp")

Enum<SIMP::System> SIMP::system;

SIMP::SIMP() : ErsatzMaterial()
{
  ParamNode* simp_pn = pn->Get("SIMP");  
  system_ = system.Parse(simp_pn->Get("system"));

  // There might be a filter regularization based on the design element.
  if(simp_pn->Has("regularization", "type", "filter"))
  {
    StdVector<ParamNode*> list = simp_pn->Get("regularization")->GetList("filter");
    // this is save for design=polarization
    for(unsigned int i = 0; i < list.GetSize(); i++)
      SIMPElement::InitFilter(design->data, list[i]);
  }
  else 
  {
    if(simp_pn->Has("regularization"))
      throw Exception("regularization not implemented"); 
  }

  // has to be done after we have mech and we have design
  GetElementMatrix(GetForm(mech, mech, "linElastInt"), mechStiffness);
  if(harmonic)
    GetElementMatrix(GetForm(mech, mech, "MassInt"), mechMass);
  // handle surfaceNormal in PostInit()
}

SIMP::~SIMP()
{
  if(mechRHS.vec != NULL)       { delete mechRHS.vec;  mechRHS.vec = NULL;  }
}

void SIMP::PostInit()
{
  if(harmonic) mechRHS.Init<complex<double> >(design, PRESSURE); // in many cases NULL;
          else mechRHS.Init<double>(design, PRESSURE);
  
  // check
  if(cost->type == RADIATION && !harmonic)
    throw Exception("objective 'radiation' is only defined for harmonic");

  // all other actions in ErsatzMaterial
  if(cost->type == RADIATION)
  {
    // note, that the surface is not the optimization region with is all material
    ParamNode* pn = param->Get("optimization")->Get("costFunction");
    std::string region = pn->Get("radiation")->Get("surfRegion")->Get("name")->AsString();
    RegionIdType sreg = domain->GetGrid()->RegionNameToId(region);
    BaseForm* form = assemble_->GetBiLinForm(sreg, mech, mech, "SurfaceNormalInt")->GetIntegrator();
    // the surface elements for the matrix are not design elements! 
    // Therefore provide an explicit element
    StdVector<Elem*> elems;
    domain->GetGrid()->GetElems(elems, sreg);
    assert(elems.GetSize() != 0);
    GetElementMatrix(form, surfaceNormal, elems[0]);
  }

  ErsatzMaterial::PostInit();  
}

BiLinFormContext* SIMP::CreateSurfaceNormalMatrix(SinglePDE* pde, BaseMaterial* baseMat, shared_ptr<ResultInfo> result)
{
  // killme! this is called in mechPDE, it would be cooler if one could really add this
  // bilinear form in the optimization PostInit() ;(
  
  // check if we have this setting, otherwise return NULL
  if(!param->Has("optimization")) return NULL;
  // using the Enum stuff would be cooler but yet there is not prior SetEnums() calling
  if(param->Get("optimization")->Get("costFunction")->Get("type")->AsString() != "radiation") return NULL;

  // read what to optimize
  Grid* grid = domain->GetGrid();

  // killme works only for one region
  std::string region = param->Get("optimization")->Get("SIMP")->Get("radiation")->Get("surfRegion")->Get("name")->AsString();
  shared_ptr<EntityList> list = 
    grid->GetEntityList(EntityList::SURF_ELEM_LIST, region, EntityList::REGION);
  
  // will be automatically deleted in the BiLinFormContext destructor by assemble
  SurfaceNormalInt* sni = new SurfaceNormalInt(baseMat); 

  // this bilinear form is to be assembled in an extra matrix. it is
  // only used to calculate the rhs for radiation optimization
  BiLinFormContext * bilifc = new BiLinFormContext(sni, AUXILIARY);
  bilifc->SetPtPdes(pde, pde); // actually we don't need the pde
  
  bilifc->SetResults(result, result, list, list);

  return bilifc;
}

template <class T>
void SIMP::SetElementK(DesignElement* de, Application app, CFSMatrix* mat_out)
{
  Matrix<T>& out = dynamic_cast<Matrix<T>& >(*mat_out);
  
  switch(app)
  {
  case SURFACE_NORMAL:
    // this does not depend on a transfer fuction!
    // also surfaceNormal is always real and might become complex with imag = 0
    Assign(out, surfaceNormal, 1.0);
    break;
    
  case MECH:
  {
    out.Resize(mechStiffness.GetSizeRow(), mechStiffness.GetSizeCol());
    
    // Find the transferfunction for K (e.g. DENSITY, MECH)
    TransferFunction* tf = design->GetTransferFunction(de->GetType(), app);
    double k_factor = tf->Derivative(de);
    
    // copy from real mechStiffness to potential complex out and factor the derivative
    Assign(out, mechStiffness, k_factor);
    LOG_DBG3(simp) << "SetElementK: org mech " << out.ToString(0);
    
    if(harmonic)
    {
      tf = design->GetTransferFunction(de->GetType(), MASS);
      AddMassToStiffness(tf->Derivative(de), dynamic_cast<Matrix<complex<double> >& >(out));
    }
    break;
  }
  default:
    assert(false); // other cases should be handled in PiezoSIMP
  } // end switch 
}


void SIMP::AddMassToStiffness(double m_factor, Matrix<complex<double> >& K_in_S_out)
{
  // The result matrix is 
  // S = K + i*omega*C - omega^2*M
  // with purly imaginary C = alpha*M+beta*K
  // S = K + i*alpha_k*K + i*alpha_,*M - omega^2*M
  // with m_factor
  // S = K + i*alpha_k*K + i*alpha_m*m_factor*M - omega^2*m_factor*M
  // With S = K in the beginning this is
  // S += i*alpha_k* + S + (i*alpha_m*m_factor-omega^2*m_factor)*M 
  
  // change name only
  Matrix<complex<double> >& S = K_in_S_out;
  Matrix<double>& M = mechMass;
  assert(S.GetSizeRow() == M.GetSizeRow() && S.GetSizeCol() == M.GetSizeCol());
    
  // find alpha, beta and omega
  double omega = mech->GetSolveStep()->GetActFreq();
  double alpha_k = 0.0;
  double alpha_m  = 0.0;

  // do we have damping (C = alpha*M+beta*K) -> this is pure imaginary!
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
  for(unsigned int r = 0; r < S.GetSizeRow(); r++)
    for(unsigned int c = 0; c < S.GetSizeCol(); c++)
      S[r][c] = complex<double>(S[r][c].real(), omega * alpha_k * S[r][c].real());
  
  // we the add the M part of C and the real mass part 
  complex<double> damp_mass = complex<double>(-1.0 *omega*omega*m_factor, omega*alpha_m*m_factor);
  for(unsigned int r = 0; r < S.GetSizeRow(); r++)
    for(unsigned int c = 0; c < S.GetSizeCol(); c++)
      S[r][c] += damp_mass * mechMass[r][c];

  LOG_DBG2(simp) << "AddMassToStiffness: m_factor:" << m_factor << " alpha_k: " << alpha_k << " alpha_m: " << alpha_m 
                 << " omega: " << omega << " K_img: " << (omega * alpha_k) << " damp_mass: " << damp_mass;
}


double SIMP::CalcRadiation()
{
  design->WriteDesignToExtern(last_evaluation.GetPointer());

  // get surface normal matrix
  StdMatrix* snm = assemble_->GetAlgSys()->GetSysMat(AUXILIARY);
  assert(snm != NULL);
  // assert(out.G)assert(snm->GetMaxDiag() != 0.0);

  // Make an OLAS Vector and feed it with the solution vector
  Vector<complex<double> >& sol   = dynamic_cast<Vector<complex<double> >& >(*forward_->raw[MECH]);
  OLAS::Vector<complex<double> > olas_sol;
  // OLAS is one based!!
  olas_sol.Replace(sol.GetSize(), sol.GetPointer()-1, false);

  // The result is also an olas vector
  OLAS::Vector<complex<double> > olas_prod;
  olas_prod.Resize(olas_sol.GetSize());

  // multiply (S_n U)
  snm->Mult(olas_sol, olas_prod);
  LOG_DBG(simp) << "radiation objective: ||S_n*U||=" << olas_prod.NormEuclid();

  // scalar product U (S_n U)
  complex<double>  tsp; 
  olas_sol.Inner(olas_prod, tsp);
  LOG_DBG(simp) << "radiation objective: U*S_n*U=" << tsp;
  // check that the scalar product is (numerical) real
  assert(((complex<double>) tsp).imag() < 1e-10 * ((complex<double>) tsp).real());
  double sp = ((complex<double>) tsp).real();

  // the objective (Du, Olhoff 2007) is pi=0.5 * gamma * c * omega^2* U*S_n*U
  double omega = mech->GetSolveStep()->GetActFreq();
  double c     = pn->Get("radiation")->Get("soundSpeed")->AsDouble();
  double gamma = pn->Get("radiation")->Get("specificMass")->AsDouble();

  double result = 0.5 * gamma * c * omega * omega * sp;
  LOG_DBG(simp) << " radiation objective: 0.5 * " << gamma << " * " << c  
  << " *  " << omega << "^2 * " << sp << " -> " << result;
  cost->SetValue(result);

  return cost->GetValue();
}

template <class T>
void SIMP::CalcSurfaceNormalTimesSolution(OLAS::Vector<T>& olas_prod)
{
  // get surface normal matrix
  StdMatrix* snm = assemble_->GetAlgSys()->GetSysMat(AUXILIARY);
  assert(snm != NULL);
  assert(snm->GetMaxDiag() != 0.0);

  // Make an OLAS Vector and feed it with the solution vector
  Vector<T>& sol   = dynamic_cast<Vector<T>& >(*forward_->raw[MECH]);
  OLAS::Vector<T> olas_sol;
  // OLAS is one based!!
  olas_sol.Replace(sol.GetSize(), sol.GetPointer()-1, false);
  LOG_DBG3(simp) << "Solution 'backcasted' to OLAS::Vector -> " << olas_sol.ToString(); 

  // The result is also an olas vector
  olas_prod.Resize(olas_sol.GetSize());
  
  // multiply (S_n U)
  snm->Mult(olas_sol, olas_prod);
  LOG_DBG(simp) << "||S_n*U||=" << olas_prod.NormEuclid();
  LOG_DBG3(simp) << "S_n * U -> " << olas_prod.ToString();
}

void SIMP::CalcObjectiveGradient(double* grad_out)
{
  TransferFunction* tf = design->GetTransferFunction(DesignElement::DENSITY, MECH, true);
  
  switch(cost->type)
  {
  case COMPLIANCE:
    // calculate the compliance which is according to 
    // "A 99 line topology optimization code written in Matlab"; O.Sigmund, 2001
    // -> dc/dx_e = -p * x_e ^(p-1) u_e^T k_0 u_e
    // Here we use a constant element stiffness matrix
    CalcU1KU2(tf, forward_->elem[MECH], MECH, forward_->elem[MECH]);
    break;

  case GLOBAL_DYNAMIC_COMPLIANCE:
    CalcU1KU2(tf, forward_->elem[MECH], MECH, forward_->elem[MECH]);
    break;

  case OUTPUT:
    // synthesis of compliant mechanism: As our adjoint PDE
    // c' = l K' u
    CalcU1KU2(tf, adjoint_->elem[MECH], MECH, forward_->elem[MECH]);
    break;
    
  case RADIATION:
  {
    // copy & paste from CalcObjective
    double omega = mech->GetSolveStep()->GetActFreq();
    double c     = pn->Get("radiation")->Get("soundSpeed")->AsDouble();
    double gamma = pn->Get("radiation")->Get("specificMass")->AsDouble();
    
    // for constant rhs: grad = gamma * c * omega^2 * -1 * U_s^T * S' * U
    // U_s is adjoint solution, S' is derivative of S = stiffness, masss (, damping)
    double factor = -1.0 * gamma * c * omega * omega; 
    
    CalcU1KU2(tf, adjoint_->elem[MECH], MECH, forward_->elem[MECH], NULL, false, factor);
    break;
  }
  default: throw Exception("objective gradient no handled");
  }

  if(grad_out != NULL) design->WriteGradientToExtern(grad_out, DesignElement::DENSITY, 
                               DesignElement::COST_GRADIENT, DesignElement::SMART);
}

void SIMP::CalcRadiationAdjointRHS()
{
  // the rhs for radiation is S_n * U

  // We have to use an OLAS Vector here
  OLAS::Vector<std::complex<double> > rhs;
  CalcSurfaceNormalTimesSolution(rhs); // does S_n * U as we also need it in the objective

  // rhs is OLAS 1-based but handles this internally
  assemble_->GetAlgSys()->InitRHS(const_cast<const OLAS::Vector<std::complex<double> >* >(&rhs));
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
  // vec is the base CFSVector which has no abstact operators overloaded
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

