#include <assert.h>
#include <cmath>
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
#include "Driver/FormsContexts.hh"
#include "Forms/LinForms/LinearForm.hh"
#include "Forms/BiLinForms/BiLinearForm.hh"
#include "Forms/BiLinForms/BDBInt.hh"
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
#include "PDE/LatticeBoltzmannPDE.hh"
#include "Utils/StdVector.hh"
#include "Utils/mathParser/mathParser.hh"
#include "Utils/tools.hh"

namespace CoupledField {
class DenseMatrix ;
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
  ErsatzMaterial::PostInit();
  
  // FIXME
  if(context->IsComplex()) mechRHS.Init<complex<double> >(design, App::PRESSURE); // in many cases NULL;
                      else mechRHS.Init<double>(design, App::PRESSURE);
}

void SIMP::SetElementK( Function* f, DesignElement* de, const TransferFunction* tf, App::Type app, DenseMatrix* out, bool derivative, CalcMode calcMode, double ev)
{
  if(f->ctxt->IsComplex())
  {
    if(f->ctxt->mat->ComplexElementMatrix(de->elem->regionId)) // handles also bloch which real material but complex BOp
      SetElementK<Complex, Complex >( f, de, tf, app, out, derivative, calcMode, ev);
    else
      SetElementK<Complex, double >( f, de, tf, app, out, derivative, calcMode, ev);
  }
  else
    SetElementK<double,double>( f, de, tf, app, out, derivative, calcMode, ev);
}

template <class T1, class T2>
void SIMP::SetElementK(Function* f, DesignElement* de, const TransferFunction* tf, App::Type app, DenseMatrix* mat_out, bool derivative, CalcMode calcMode, double ev)
{
  assert(f->ctxt->mat != NULL);
  OptimizationMaterial* mat = f->ctxt->mat;
  Matrix<T1>& out = dynamic_cast<Matrix<T1>& >(*mat_out);
  //std::cout << "out= " << out.ToString() << std::endl;

  //assert(app != App::MAG); // shall be in MagSIMP.cc

  switch(app)
  {
  case App::MECH:
  case App::ACOUSTIC:
  case App::HEAT:
  {
    int mm = de->multimaterial != NULL ? de->multimaterial->index : -1;

    // element matrix with org material, might be cached -> local_element_cache
    const Matrix<T2>& stiffness = dynamic_cast<const Matrix<T2>& >(mat->Stiffness(de->elem, false, mm)); // no bimaterial

    // Find the transfer function for K (e.g. DENSITY, App::MECH)
    T1 k_factor = derivative ? tf->Derivative(de, DesignElement::SMART, false) : tf->Transform(de, DesignElement::SMART);// not the bimat case

    // copy from real mechStiffness to potential complex out and factor the derivative
    Assign(out, stiffness, k_factor); // out = k_factor * stiffness
    // This log is very expensive, it blows up inv_tensor in the debug mode
    // LOG_DBG3(simp) << "SetElementK: el=" << de->elem->elemNum << " di=" << de->GetIndex() << " mm=" << mm << " K_org=" <<  stiffness.ToString() << " k_factor " << k_factor << " -> " << out.ToString();

    if(design->GetRegion(de->elem->regionId)->HasBiMaterial())
    {
      const Matrix<T2>& bimat = dynamic_cast<const Matrix<T2>& >(mat->Stiffness(de->elem, true)); // yes, bimaterial
      // rho^3 * E1 + (1-rho)^3 * E2, in the derivative case 3*rho^2 * E1 - 3*(1-rho)^2 * E2
      k_factor = !derivative ? 1.0 - k_factor : -1.0 *  k_factor;
      Add(out, k_factor, bimat);
      // LOG_DBG3(simp) << "SetElementK: bimat k_factor " << k_factor << " bimat=" << bimat.ToString() << " -> " << out.ToString();

      // LOG_DBG3(simp) << "SetElementK: K_bi_org=" <<  bimat.ToString() << " k_factor " << k_factor << " -> " << out.ToString();
    }

    if(f->ctxt->IsComplex())
    {
      tf = design->GetTransferFunction(de->GetType(), App::MASS);
      AddMassToStiffness(f->ctxt, tf, de, dynamic_cast<Matrix<complex<double> >& >(out), derivative, false, calcMode, ev); // no bimaterial

      // LOG_DBG3(simp) << "SetElementK: m_factor " << m_factor << " -> " << out.ToString();

      if(design->GetRegion(de->elem->regionId)->HasBiMaterial())
      {
        // rho^3 * E1 + (1-rho^3) * E2, in the derivative case 3*rho^2 * E1 - 3*rho^2 * E2
        AddMassToStiffness(f->ctxt, tf, de, dynamic_cast<Matrix<complex<double> >& >(out), derivative, true, calcMode, ev); // bimaterial

        // LOG_DBG3(simp) << "SetElementK: m_bi_factor " << m_factor << " -> " << out.ToString();
      }
    }
    break;
  }

  case App::ELEC:
  {
    const Matrix<std::complex<double> >& stiffness = dynamic_cast<ElecMat *>(mat)->ElecStiffness(de->elem, false); // no bimaterial

    // Find the transfer function for K (e.g. DENSITY, App::MECH)
    T1 k_factor = derivative ? tf->Derivative(de, DesignElement::SMART) : tf->Transform(de, DesignElement::SMART);

    // copy from ElecStiffness to out and factor the derivative
    if(f->ctxt->IsComplex())
      Assign(out, dynamic_cast<const Matrix<T1>& >(stiffness), k_factor);
    else
      Assign(out, stiffness.GetPart(Global::REAL), k_factor);

    // This log is very expensive, it blows up inv_tensor in the debug mode
    // LOG_DBG3(simp) << "SetElementK: K_org=" <<  stiffness.ToString() << " k_factor " << k_factor << " -> " << out.ToString();

    if(design->GetRegion(de->elem->regionId)->HasBiMaterial())
    {
      const Matrix<std::complex<double> >& bimat = dynamic_cast<ElecMat *>(mat)->ElecStiffness(de->elem, true); // yes, bimaterial
      // rho^3 * E1 + (1-rho^3) * E2, in the derivative case 3*rho^2 * E1 - 3*rho^2 * E2
      k_factor = derivative ? tf->Derivative(de, DesignElement::SMART, true) : tf->Transform(de, DesignElement::SMART, true);
      if(f->ctxt->IsComplex())
        Add(out, k_factor, dynamic_cast<const Matrix<T1>& >(bimat));
      else
        Add(out, k_factor, bimat.GetPart(Global::REAL));
      // LOG_DBG3(simp) << "SetElementK: K_bi_org=" <<  bimat.ToString() << " k_factor " << k_factor << " -> " << out.ToString();
    }
    break;
  }

  // u1^T (K' u2 - f') -> find "K'"
  case App::MAG:
  {
    MagMat* mag = dynamic_cast<MagMat*>(mat);
    assert(mag != NULL);

    assert(derivative);
    double d_rho = tf->Derivative(de, DesignElement::SMART, false);

    // element matrix with org material, might be cached -> local_element_cache
    const Matrix<T2>& stiffness = dynamic_cast<const Matrix<T2>& >(mag->Stiffness(de->elem));
    LOG_DBG3(simp) << "e=" << de->elem->elemNum << " K_0=" << stiffness.ToString(2);

    double nu_r = mag->nu_r[de->elem->regionId];
    double nu_0 = mag->nu_0;

    //assert(nu_r > 0);
    //assert(mag->nu_0 > 0);

    // simulation: BDB with D=(d 0; 0 d) with d = nu_0*nu_r
    // optimization: d = nu_0 * nu_r * f(rho) - nu_0 * f(rho) + nu_0
    // derivative: d = nu_0 * nu_r * f'(rho) - nu_0 * f'(rho)

    // in difference to elasticity where K' = my(rho)'*K_0
    // we have K' = alpha * K_0 where alpha = (nu_0*nu_r - nu_0) / (nu_0 * nu_r)

    double alpha = (nu_0*nu_r*d_rho - nu_0) / (nu_0 * nu_r);
    LOG_DBG3(simp) << "e=" << de->elem->elemNum << " "
        "nu_0=" << nu_0 << " nu_r=" << nu_r << " rho=" << de->GetDesign(DesignElement::SMART) << " d_rho=" << d_rho << " alpha=" << alpha;

    Assign(out, stiffness, alpha); // out = alpha * stiffness

    LOG_DBG3(simp) << "out=" << out.ToString(2);
    break;
  }

  default:
    assert(false); // other cases should be handled in PiezoSIMP
  } // end switch
}


double SIMP::CalcFunction(Excitation& excite, Function* f, bool derivative)
{
  // this app is for the PDE
  App::Type app = f->ctxt->ToApp();

  if(!derivative)
    return ErsatzMaterial::CalcFunction(excite, f, derivative);

  // only special derivatives, the rest also EM
  switch(f->GetType())
  {
  case Function::STRESS:
  case Function::STRESS_DENSITY:
  {
    TransferFunction* tf = design->GetTransferFunction(DesignElement::Default(f->ctxt), TransferFunction::Default(f->ctxt), true);
    CalcVonMisesStressGradient(excite, f, tf);
    break;
  }
  case Function::SQR_MAG_FLUX_DENS_X:
  case Function::SQR_MAG_FLUX_DENS_Y:
  {
    TransferFunction* tf = design->GetTransferFunction(DesignElement::Default(f->ctxt), TransferFunction::Default(f->ctxt), true, true);
    CalcMagFluxDensGradient(excite, f, tf);
    break;
  }

  case Function::GLOBAL_DYNAMIC_COMPLIANCE:
  case Function::OUTPUT:
  case Function::SQUARED_OUTPUT:
  case Function::DYNAMIC_OUTPUT:
  case Function::CONJUGATE_COMPLIANCE:
  case Function::ABS_OUTPUT:
  {
    // synthesis of compliant mechanism: As our adjoint PDE
    // c' = l K' u
    TransferFunction* tf = design->GetTransferFunction(DesignElement::Default(f->ctxt), TransferFunction::Default(f->ctxt), true, true); // exception and use_single
    double weight = excite.GetWeightedFactor(f);

    // squared output gradient 2 * <u,l> * <u,l>'
    if (f->GetType() == Function::SQUARED_OUTPUT)
    {
      f->SetType(Function::OUTPUT);
      weight *= 2. * SIMP::CalcFunction(excite, f, false);
      f->SetType(Function::SQUARED_OUTPUT);
    }
    LOG_DBG(simp) << "CalcFunction(idx=" << excite.index << ") norm_weight= " <<  excite.normalized_weight  << " factor=" << excite.GetFactor(f) << " weight=" << weight;
    CalcU1KU2(tf, adjoint.Get(excite, f)->elem[app], app, forward.Get(excite)->elem[app], NULL, weight, STANDARD, f);
    break;
  }

  case Function::PRESSURE_DROP:
  {
    LatticeBoltzmannPDE* lbmPde = f->ctxt->GetLatticeBoltzmannPDE();
    assert(lbmPde != NULL);
    lbmPde->SensitivityAnalysis(design->GetTransferFunction(f->elements[0]), f, design);
    break;
  }

  default:
    return ErsatzMaterial::CalcFunction(excite, f, derivative);
  }

  return 0.0; // only derivatives evaluated
}


void SIMP::CalcVonMisesStressGradient(Excitation& excite, Function* f, TransferFunction* tf)
{
  assert(excite.sequence == f->ctxt->sequence);
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
    LOG_DBG2(simp) << "CVMSG: f=" << f->ToString() << " de=" << design->data[i].elem->elemNum << " org=" << design->data[i].GetPlainGradient(f);

  // alpha is from the globalization which is in the form sum max(0, g_i-c)^p and alpha is p*max(0, g_i-c)^(p-1) where g_i is the vonMisesStress
  Vector<double> alpha;
  // 2 * stress^T * M * (rho^p)' * E_0 * B * u
  Vector<double> appendix;

  if(f->ctxt->IsComplex())
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
  rhs.Init<double>(App::STRESS, excite.label);
  // calc lambda^T *  K' * u -> this already stores the results by AddGradient()!
  CalcU1KU2(tf, adjoint.Get(excite, f)->elem[App::MECH], App::MECH, forward.Get(excite)->elem[App::MECH], &rhs, 1.0, STANDARD, f);

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
    LOG_DBG2(simp) << "CVMSG: f=" << f->ToString() << " de=" << de.elem->elemNum << " idx=" << idx << " alpha="
                   << (idx != -1 ? alpha[i] : -1.0)  << "* app=" << (idx != -1 ? appendix[i] : -1.0) << " -> " << de.GetPlainGradient(f);
	}
}

void SIMP::CalcMagFluxDensGradient(Excitation& excite, Function* f, TransferFunction* tf)
{

  // the context->GetExcitation() is now the last one as we solve and store all excitations first before calculating the gradients
  //Transform* trans = f != NULL && f->GetExcitation() != NULL ? f->GetExcitation()->transform : NULL; // even ->transform might be NULL

  // the gradient is < lambda^T, K' * A >
  assert(excite.sequence == f->ctxt->sequence);

  DesignDependentRHS rhs;
  // calc lambda^T *  K' * A -> this already stores the results by AddGradient()!
  CalcU1KU2(tf, adjoint.Get(excite, f)->elem[App::MAG], App::MAG, forward.Get(excite)->elem[App::MAG], &rhs, 1.0, STANDARD, f);

}

DesignDependentRHS::DesignDependentRHS()
{
  valid       = false;
  app         = App::NO_APP;
  vec         = NULL;
  elem        = NULL;
  test_strain = MechPDE::NOT_SET;
  isInterfaceDriven_ = false;
}

DesignDependentRHS::~DesignDependentRHS()
{
  valid = false;
  if(vec != NULL) { delete vec; vec = NULL; }
}

template <class T>
bool DesignDependentRHS::Init(DesignSpace* design, App::Type app)
{
  assert(app == App::CHARGE_DENSITY || app == App::PRESSURE || app == App::HEAT);

  if (app == App::HEAT) {
    valid = true;
    isInterfaceDriven_ = true;
    return true;
  }

  std::string name = app == App::CHARGE_DENSITY ? "LinNeumannInt" : "PressureLinForm";


  // check if we have a form with the application name
  LinearForm* form = NULL;
  LinearFormContext* actContext = NULL;

  SinglePDE* mech = Optimization::context->ToPDE(App::MECH, false);
  if(mech == NULL) return false; // wrong pde -> extend if you need it!

  StdVector<LinearFormContext*>& forms = mech->GetAssemble()->GetLinForms();

  for(StdVector<LinearFormContext*>::iterator it = forms.Begin(); it != forms.End(); it++)
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
bool DesignDependentRHS::Init(App::Type app, std::string excite_label)
{
  assert(app == App::STRESS);
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
template bool DesignDependentRHS::Init<double>(DesignSpace* design, App::Type app);
template bool DesignDependentRHS::Init<complex<double> >(DesignSpace* design, App::Type app);
#endif
