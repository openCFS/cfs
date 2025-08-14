#include <assert.h>
#include <cmath>
#include <stddef.h>
#include <map>
#include <ostream>
#include <string>

#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Domain/Domain.hh"
#include "Domain/ElemMapping/Elem.hh"
#include "Domain/ElemMapping/EntityLists.hh"
#include "Domain/ElemMapping/SurfElem.hh"
#include "Driver/Assemble.hh"
#include "Driver/BucklingDriver.hh"
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
  
  // if this hurts you, overwrite InitSecondMaterialCache()
  InitSecondMaterialCache();

  // FIXME
  if(context->IsComplex()) mechRHS.Init<complex<double> >(design, App::PRESSURE); // in many cases NULL;
                      else mechRHS.Init<double>(design, App::PRESSURE);
}

void SIMP::InitSecondMaterialCache()
{
  // we are mechanic, elec or lbm in SIMP. The others are MagSIMP, AcouSIMP, PiezoSIMP 
  // TODO make ElecSIMP once we do electric bi/ground material. 

  AddSecondMaterialCache(MaterialClass::MECHANIC, MaterialType::MECH_STIFFNESS_TENSOR);
  AddSecondMaterialCache(MaterialClass::MECHANIC, MaterialType::DENSITY);
}

void SIMP::AddSecondMaterialCache(MaterialClass mc, MaterialType mt)
{
  for (StdVector<DesignSpace::DesignRegion>& drv: design->regions)
  {
    for (DesignSpace::DesignRegion& dr: drv)
    {
      if (dr.HasScndMaterial())
      {
        // this also does some caching such that subsequent GetScndMaterial() 
        // are possible from threads
        PtrCoefFct coef = dr.GetScndMaterial(mc, mt);
        const Elem* elem = design->data[dr.base].elem; // first element
        shared_ptr<ElemShapeMap> esm = grid->GetElemShapeMap(elem);
       
        assert(elem->extended != nullptr);
        assert(elem->extended->barycenter.data.GetSize() >= 2);
        LocPoint lp(elem->extended->barycenter.data);
        LocPointMapped lpm;
        lpm.Set(lp, esm);
        
        LOG_DBG(simp) << "ASMC r=" << dr.regionId << " d=" << dr.design << " mc=" << MaterialClassEnum.ToString(mc) 
                      << " mt=" << MaterialTypeEnum.ToString(mt) << " dim=" << CoefFunction::coefDimType.ToString(coef->GetDimType())
                      << " cplx=" << coef->IsComplex() << " coef=" << coef->ToString();

        assert(coef->GetDimType() == CoefFunction::SCALAR || coef->GetDimType() == CoefFunction::VECTOR || coef->GetDimType() == CoefFunction::TENSOR);

        if(coef->GetDimType() == CoefFunction::TENSOR)
        {
          // we don't ask the pde if it is complex, as for most harmonic pdes, the material is real
          if(!coef->IsComplex())
          {
            Matrix<double> val;
            coef->GetTensor(val, lpm);
            dr.scnd_material_cached[mc][mt] = val;
          }
          else
          {
            Matrix<complex<double>> val;
            coef->GetTensor(val, lpm);
            dr.scnd_material_cached[mc][mt] = val;
          }
        }
        else
        {
          if(!coef->IsComplex()) 
          {
            // it is so stupid to have no void Get*() functions!
            double val;
            coef->GetScalar(val, lpm);
            dr.scnd_material_cached[mc][mt] = val;
          }
          else
          {
            complex<double> val;
            coef->GetScalar(val, lpm);
            dr.scnd_material_cached[mc][mt] = val;
         }
        }
      }
      else
      {
        LOG_DBG(simp) << "ASMC r=" << dr.regionId << " d=" << dr.design << " has no second material";
      }
    }
  }
}

void SIMP::SetElementK( Function* f, DesignElement* de, const TransferFunction* tf, App::Type app, DenseMatrix* out, bool derivative, CalcMode calcMode, double ev)
{
  if(f->ctxt->IsComplex())
  {
    if(f->ctxt->mat->ComplexElementMatrix(de->elem->regionId)) // handles also bloch with real material but complex BOp
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

  //assert(app != App::MAG); // shall be in MagSIMP.cc

  switch(app)
  {
  case App::MECH:
  case App::ACOUSTIC:
  case App::HEAT:
  case App::BUCKLING:
  {
    int mm = de->multimaterial != NULL ? de->multimaterial->index : -1;

    // element matrix with org material, might be cached -> local_element_cache
    const Matrix<T2>& stiffness = dynamic_cast<const Matrix<T2>& >(mat->Stiffness(de->elem, false, mm)); // no bimaterial
    LOG_DBG3(simp) << "stiffness=" << stiffness.ToString();

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

    if(app == App::BUCKLING)
    {
      assert(f->ctxt->DoBuckling());

      if (f->ctxt->GetBucklingDriver()->IsInverseProblem())
        out *= -ev;

      tf = design->GetTransferFunction(de->GetType(), App::BUCKLING);
      AddGeometricStiffnessToStiffness(f->ctxt, tf, de, dynamic_cast<Matrix<Complex>& >(out), derivative, false, calcMode, ev); // no bimaterial
      // LOG_DBG3(simp) << "SetElementK: GeoStiff out -> " << out.ToString();

      if(design->GetRegion(de->elem->regionId)->HasBiMaterial())
      {
        // rho^3 * E1 + (1-rho^3) * E2, in the derivative case 3*rho^2 * E1 - 3*rho^2 * E2
        AddGeometricStiffnessToStiffness(f->ctxt, tf, de, dynamic_cast<Matrix<Complex>& >(out), derivative, true, calcMode, ev); // bimaterial
        // LOG_DBG3(simp) << "SetElementK: GeoStiff bimat out " -> " << out.ToString();
      }
    }

    if(f->ctxt->IsComplex() && !f->ctxt->DoBuckling())
    {
      tf = design->GetTransferFunction(de->GetType(), App::MASS);
      AddMassToStiffness(f->ctxt, tf, de, dynamic_cast<Matrix<complex<double> >& >(out), derivative, false, calcMode, ev); // no bimaterial
      // LOG_DBG3(simp) << "SetElementK: Mass out -> " << out.ToString();

      if(design->GetRegion(de->elem->regionId)->HasBiMaterial())
      {
        // rho^3 * E1 + (1-rho^3) * E2, in the derivative case 3*rho^2 * E1 - 3*rho^2 * E2
        AddMassToStiffness(f->ctxt, tf, de, dynamic_cast<Matrix<complex<double> >& >(out), derivative, true, calcMode, ev); // bimaterial
        // LOG_DBG3(simp) << "SetElementK: Mass bimat out " -> " << out.ToString();
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

  default:
    assert(false); // other cases should be handled in PiezoSIMP
  } // end switch
}


double SIMP::CalcFunction(Excitation& excite, Function* f, bool derivative)
{
  // this app is for the PDE
  App::Type app = f->ctxt->ToApp();

  // an exceptional case where we also calculate the function value as we share a lot code with the gradient
  if(f->GetType() == Function::GLOBAL_STRESS)
  {
    if(f->ctxt->IsComplex())
      return CalcGlobalVonMisesStressOrLoadFactor<Complex>(excite, f, derivative);
    else
      return CalcGlobalVonMisesStressOrLoadFactor<double>(excite, f, derivative);
  }

  // microscopic load factor is interpolated value divided by local stress
  if(f->GetType() == Function::LOCAL_STRESS || f->GetType() == Function::LOCAL_BUCKLING_LOAD_FACTOR)
    return ErsatzMaterial::CalcLocalVonMisesStressOrLoadFactor(excite, f, derivative);

  if(!derivative)
    return ErsatzMaterial::CalcFunction(excite, f, derivative);

  // only special derivatives, the rest also EM
  switch(f->GetType())
  {
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

  case Function::REFLECTED_WAVE:
  {
    // J = (u - z)^T L (u - z)*
    // J' = 2*Re(lambda K' (u-z))
    // K^T lambda = - L (u-z)
    // Similar to dynamicOutput, but the solution gets shifted by a constant value z. 
    // In acoustics, z is the pressure caused by the excitation, so we only see the reflected wave.
    
    TransferFunction* tf = design->GetTransferFunction(DesignElement::Default(f->ctxt), TransferFunction::Default(f->ctxt), true, true); // exception and use_single
    double weight = excite.GetWeightedFactor(f);

    LOG_DBG(simp) << "CalcFunction(idx=" << excite.index << ") norm_weight= " <<  excite.normalized_weight  << " factor=" << excite.GetFactor(f) << " weight=" << weight;

    // get the z paramter and the region of excitation
    RegionIdType reg = GetExcitationRegion(f);
    Complex z = GetExcitationPressure(f);
    StdVector<SingleVector*>& u2 = forward.Get(excite)->elem[app];
    StdVector<Vector<complex<double>>> umz(u2.GetSize());

    // This code is only for excitation inside the optimization region and not verified
    // get all nodeIds for the ecitation
    StdVector<unsigned int> regNodes;
    grid->GetNodesByRegion(regNodes, reg);
    LOG_DBG(simp) << "CF: reg=" << reg << " -> nodes=" << regNodes.ToString(); 
    assert(umz.GetSize() == design->data.GetSize());
    // loop over all design elements
    for(DesignElement& de : design->data){
      unsigned int idx = de.GetIndex();
      Vector<complex<double>>* vec = dynamic_cast<Vector<complex<double>>*>(u2[idx]);
      umz[idx] = *vec;
      assert(umz[idx].GetSize() == de.elem->connect.GetSize());
      LOG_DBG3(simp) << "CF: Set u-z for: elem idx=" << idx << " test for=" << de.elem->connect.ToString(); 
      // loop over all nodes of the design elements
      for(unsigned int node_idx = 0; node_idx < de.elem->connect.GetSize(); node_idx++){
        // check if one of the nodes is in the excitation region -> substract z
        if(regNodes.Contains(de.elem->connect[node_idx])) {
          throw Exception("Excitation inside optimization region is not verified.");
          umz[idx][node_idx] -= z;
          LOG_DBG3(simp) << "CF: Set u-z for: elem idx=" << idx << " node idx=" << node_idx << " -> " << umz[idx][node_idx]; 
        }
      }
    }
    StdVector<SingleVector*> umzptr(u2.GetSize());
    for(unsigned int i = 0; i < u2.GetSize(); i++)
      umzptr[i] = &umz[i];

    CalcU1KU2(tf, adjoint.Get(excite, f)->elem[app], app, umzptr, NULL, weight, STANDARD, f);
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


DesignDependentRHS::DesignDependentRHS(App::Type my_app)
{
  app         = my_app;
  valid       = app == App::MAG; // MAG needs no init
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
bool DesignDependentRHS::Init(DesignSpace* design, App::Type my_app)
{
  assert(!(app != App::NO_APP && app != my_app && my_app != App::NO_APP));

  if(my_app != App::NO_APP)
    this->app = my_app;

  assert(app == App::CHARGE_DENSITY || app == App::PRESSURE || app == App::HEAT || app == App::MAG);

  if (app == App::HEAT) {
    this->valid = true;
    isInterfaceDriven_ = true;
    return true;
  }

  std::string name = app == App::CHARGE_DENSITY ? "LinNeumannInt" : "PressureLinForm";


  // check if we have a form with the application name
  LinearForm* form = NULL;
  LinearFormContext* actContext = NULL;

  SinglePDE* mech = Optimization::context->ToPDE(App::MECH, false);
  if(mech == NULL)
	return false; // wrong pde -> extend if you need it!

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
  if(form == NULL)
	return false; // no form, no RHS!

  // the context knows the surface elements!
  this->valid = true; // doubled code?!

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
bool DesignDependentRHS::Init(std::string excite_label, App::Type my_app)
{
  assert(!(app != App::NO_APP && app != my_app && my_app != App::NO_APP));
  if(my_app != App::NO_APP)
    app = my_app;

  assert(app == App::STRESS);
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
template bool DesignDependentRHS::Init<double>(DesignSpace* design, App::Type app);
template bool DesignDependentRHS::Init<complex<double> >(DesignSpace* design, App::Type app);
template bool DesignDependentRHS::Init<double>(std::string excite_label, App::Type my_app);
