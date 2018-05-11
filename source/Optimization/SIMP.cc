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

void SIMP::SetElementK(Context* ctxt, DesignElement* de, const TransferFunction* tf, App::Type app, DenseMatrix* out, bool derivative, CalcMode calcMode, double ev)
{
  if(ctxt->IsComplex())
  {
    if(ctxt->mat->ComplexElementMatrix(de->elem->regionId)) // handles also bloch which real material but complex BOp
      SetElementK<Complex, Complex >(ctxt, de, tf, app, out, derivative, calcMode, ev);
    else
      SetElementK<Complex, double >(ctxt, de, tf, app, out, derivative, calcMode, ev);
  }
  else
    SetElementK<double,double>(ctxt, de, tf, app, out, derivative, calcMode, ev);
}

template <class T1, class T2>
void SIMP::SetElementK(Context* ctxt, DesignElement* de, const TransferFunction* tf, App::Type app, DenseMatrix* mat_out, bool derivative, CalcMode calcMode, double ev)
{
  assert(ctxt->mat != NULL);
  OptimizationMaterial* mat = ctxt->mat;
  Matrix<T1>& out = dynamic_cast<Matrix<T1>& >(*mat_out);
  std::cout << "out= " << out.ToString() << std::endl;

  //assert(app != App::MAG); // shall be in MagSIMP.cc

  switch(app)
  {
  case App::MAG:
  case App::MECH:
  case App::ACOUSTIC:
  case App::HEAT:
  {
    int mm = de->multimaterial != NULL ? de->multimaterial->index : -1;
    //std::cout << "mm= " << mm << std::endl;

    const Matrix<T2>& stiffness = dynamic_cast<const Matrix<T2>& >(mat->Stiffness(de->elem, false, mm)); // no bimaterial
    //std::cout << "stifness= " << stiffness.ToString() << std::endl;

    // Find the transfer function for K (e.g. DENSITY, App::MECH)
    T1 k_factor = derivative ? tf->Derivative(de, DesignElement::SMART, false) : tf->Transform(de, DesignElement::SMART);// not the bimat case

    //std::cout << "k_factor= " << k_factor << std::endl;
    // copy from real mechStiffness to potential complex out and factor the derivative
    Assign(out, stiffness, k_factor); // out = k_factor * stiffness
    //std::cout << "stiffness= " << stiffness << std::endl;
    //std::cout << "out= " << out << std::endl;
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

    if(ctxt->IsComplex())
    {
      tf = design->GetTransferFunction(de->GetType(), App::MASS);
      AddMassToStiffness(ctxt, tf, de, dynamic_cast<Matrix<complex<double> >& >(out), derivative, false, calcMode, ev); // no bimaterial

      // LOG_DBG3(simp) << "SetElementK: m_factor " << m_factor << " -> " << out.ToString();

      if(design->GetRegion(de->elem->regionId)->HasBiMaterial())
      {
        // rho^3 * E1 + (1-rho^3) * E2, in the derivative case 3*rho^2 * E1 - 3*rho^2 * E2
        AddMassToStiffness(ctxt, tf, de, dynamic_cast<Matrix<complex<double> >& >(out), derivative, true, calcMode, ev); // bimaterial

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
    if (ctxt->IsComplex())
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
      if(ctxt->IsComplex())
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
  case Function::MAG_FLUX_DENS_X:
  case Function::MAG_FLUX_DENS_Y:
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
  Transform* trans = f != NULL && f->GetExcitation() != NULL ? f->GetExcitation()->transform : NULL; // even ->transform might be NULL

  // the gradient is < lambda^T, K' * A > + < J^T, B A' >

  bool derivative = true;
  Vector<double> appendix;
  assert(appendix.GetSize() == 0); // if not, we can skip resizing
  const Vector<double>& stateSol = forward.Get(excite,NULL)->GetRealVector(StateSolution::RAW_VECTOR);
  appendix.Resize(stateSol.GetSize(),0.0);

  assert(excite.sequence == f->ctxt->sequence);

  // the stored element solution vector
  StdVector<SingleVector*>& sol = forward.Get(excite)->elem[App::MAG];

  DesignDependentRHS rhs;
  // calc lambda^T *  K' * A -> this already stores the results by AddGradient()!
  CalcU1KU2(tf, adjoint.Get(excite, f)->elem[App::MAG], App::MAG, forward.Get(excite)->elem[App::MAG], &rhs, 1.0, STANDARD, f);
  //std::cout << f->ToString() << std::endl;

  // add < J^T, B dA/drho  >
  // loop over all elements in region f->region
  // get B and set it to the nodal positions via eqnMap from element nodes

  // get elements
  //int elems = design->GetNumberOfElements();
  StdVector<Elem*> elems;
  domain->GetGrid()->GetElems(elems, design->GetRegionId());
  assert(!elems.IsEmpty());

  // annoying entity iterator got hold the elem
  ElemList el(elems[0],domain->GetGrid());

  Matrix<double> b_mat; // this holds the curl-operator for the whole element. For 2D rectangular case it shall be 2 rows, 4 columns
  Vector<double> mag_flux_density;
  Vector<double> res; // in relation to the integration point
  Vector<double> resElem; // in relation to the element
  Vector<double> twoDim(2);
  Vector<double> threeDim(3);
  Vector<double> var(elems.GetSize());
  int Lauf = 0;
  double base = 0;

  for(unsigned int e = 0; e < elems.GetSize(); e++)
  {
    Elem* elem = elems[e];
    el.SetElement(elem);
    DesignElement* org = &design->data[e + base];
    DesignElement* de = design->ApplyTransformations(org, org, trans);
    double k_factor = derivative ? tf->Derivative(de, DesignElement::SMART, false) : tf->Transform(de, DesignElement::SMART);// not the bimat case
    //std::cout << "k_factor= " << k_factor << std::endl;
    //std::cout << "sol.GetSize " << sol.GetSize() << std::endl;
    //std::cout << "elem " << elem->elemNum << std::endl;
    assert(sol.GetSize() > e);
    Vector<double>* vec = dynamic_cast<Vector<double>*>(sol[e]);
    //assert(sol.GetSize() > elem->elemNum);
    //Vector<double>* vec = dynamic_cast<Vector<double>*>(sol[elem->elemNum]); // or -1 for 0-based???
    assert(vec != NULL);
    Vector<double>& a = *vec; // a stands for the Vectorpotential in the element
    //std::cout << "a= " <<a.ToString() << std::endl;
    Assign(a,a,k_factor); // a' = k_factor * a
    //std::cout << "aNeu= " <<a.ToString() << std::endl;
    //assert(!(domain->GetGrid()->IsRegionRegular(f->region) && a.GetSize() != 4)); // for regular 2D grid!!!
    //LOG_DBG2(em) << "CMDF: i=" << e << " e=" << elem->elemNum << " -> " << a.ToString();

    // we need a lot of similar stuff as in BDBInt::CalcElementMatrix().
    // get the form first
    BiLinearForm* form = context->pde->GetAssemble()->GetBiLinForm("CurlCurlIntegrator", f->region, context->pde)->GetIntegrator();
    BDBInt<>* bdb = dynamic_cast<BDBInt<>*>(form);
    assert(bdb != NULL);

    // get B-mat
    StdVector<LocPoint> intPoints; // Get integration Points
    LocPointMapped lp;
    StdVector<double> weights;
    IntegOrder order;
    IntScheme::IntegMethod method = IntScheme::UNDEFINED;
    BaseFE* ptFe = form->GetFeSpace1()->GetFe(el.GetIterator(), method, order );
    form->GetIntScheme()->GetIntPoints(Elem::GetShapeType(elem->type), method, order, intPoints, weights );
    //LOG_DBG2(em) << "CMFD i=" << e << " e=" << elem->elemNum << " method=" << method << " order=" << order.ToString() << " iP=" << intPoints;
    assert(method != IntScheme::UNDEFINED);
    assert(!intPoints.IsEmpty());

    // Get shape map from grid
    shared_ptr<ElemShapeMap> esm =domain->GetGrid()->GetElemShapeMap(elem);

    // Loop over all integration points
    for(unsigned int ip = 0; ip < intPoints.GetSize(); ip++)
    {
      // Calculate for each integration point the LocPointMapped
      lp.Set( intPoints[ip], esm, weights[ip] );

      // Call the CalcBMat()-method
      bdb->GetBOp()->CalcOpMat(b_mat, lp, ptFe);
      assert(b_mat.GetNumCols() == a.GetSize());
      assert(b_mat.GetNumRows() == domain->GetGrid()->GetDim());

      // Initialize if the problem is two or three dimensional
      if (b_mat.GetNumRows() == 2 && ip == 0)
      {
        res = twoDim;
        resElem = twoDim;
      }
      else if (b_mat.GetNumRows() == 3)
      {
        res = threeDim;
        resElem = threeDim;
      }

      // Calculation of (B * dA/drho) over the integration points
      // b_mat consists of the derivative of the shape function, a consists of the magnetic vector potential
      b_mat *= weights[ip];
      b_mat.Mult(a, res);
      std::cout << "res= " << res.ToString() << std::endl;
      resElem += res;
      std::cout << "resElem= " << resElem.ToString() << std::endl;
      if (ip+1 == intPoints.GetSize())
      {
        resElem /= intPoints.GetSize();
      }
    } // end loop elems
    if(f->GetType() == Function::MAG_FLUX_DENS_X)
    {
      var[Lauf] = resElem[0];
    }
    else if (f->GetType() == Function::MAG_FLUX_DENS_Y)
    {
      var[Lauf] = resElem[1];
    }
    else
    {
      assert(false);
    }
    Lauf ++;
    std::cout << "var= " << var << std::endl;

    de->AddGradient(f, var[e]);
  }
  /*for(unsigned int i = 0; i < design->data.GetSize(); i++)
        {
          DesignElement& de = design->data[i];
          // Three cases:
          // a) stress is defined on whole design domain (one or more regions)
          // b) stress region is defined on one of two or more design regions
          // c) stress region is not in any of the design regions.

          int idx = -1; // case c) an sometimes b) never a)



          // case b)
          if(design->Contains(f->region)  )
          {
            // we are at a design element and have to find it within stress
            // TODO make it faster
            for(unsigned int e = 0; idx != -1 && e < f->elements.GetSize(); e++)
              if(de.elem->elemNum == f->elements[e]->elem->elemNum)
                idx = e;
          }
          // case c): idx is already -1
          if(idx != -1){

          }
        }*/

  /*
  // ugly copy and paste from CalcMagFluxDensity()
  StdVector<Elem*> elems;
  domain->GetGrid()->GetElems(elems, f->region);
  // annoying entity iterator got hold the elem
  ElemList el(elems[0],domain->GetGrid());
  Matrix<double> b_mat; // this holds the curl-operator for the whole element. for 2D rectangular it shall be 2 rows, 4 columns
  BiLinearForm* form = context->pde->GetAssemble()->GetBiLinForm("CurlCurlIntegrator", f->region, context->pde)->GetIntegrator();

  // this gets the equation numbers for the element
  StdVector<int> eqn;

  // select first or second row of b_mat
  assert(f->GetType() == Function::MAG_FLUX_DENS_X || f->GetType() == Function::MAG_FLUX_DENS_Y);
  bool x_comp = f->GetType() == Function::MAG_FLUX_DENS_X;

  for(unsigned int e = 0; e < elems.GetSize(); e++)
    {
    Elem* elem = elems[e];
    el.SetElement(elem);
    DesignElement& de = design->data[e];

    Vector<double>* vec = dynamic_cast<Vector<double>*>(sol[e]);
    //assert(sol.GetSize() > elem->elemNum);
    //Vector<double>* vec = dynamic_cast<Vector<double>*>(sol[elem->elemNum]); // or -1 for 0-based???
    assert(vec != NULL);
    Vector<double>& a = *vec; // a stands for the Vectorpotential in the element

    // we need a lot of similar stuff as in BDBInt::CalcElementMatrix().
    // get the form first
    BDBInt<>* bdb = dynamic_cast<BDBInt<>*>(form);
    assert(bdb != NULL);
    // get B-mat
    StdVector<LocPoint> intPoints; // Get integration Points
    LocPointMapped lp;
    StdVector<double> weights;
    IntegOrder order;
    IntScheme::IntegMethod method = IntScheme::UNDEFINED;
    BaseFE* ptFe = form->GetFeSpace1()->GetFe(el.GetIterator(), method, order );
    form->GetIntScheme()->GetIntPoints(Elem::GetShapeType(elem->type), method, order, intPoints, weights );
    //LOG_DBG2(em) << "CMFD i=" << e << " e=" << elem->elemNum << " method=" << method << " order=" << order.ToString() << " iP=" << intPoints;
    assert(method != IntScheme::UNDEFINED);
    assert(!intPoints.IsEmpty());

    // Get shape map from grid
    shared_ptr<ElemShapeMap> esm =domain->GetGrid()->GetElemShapeMap(elem);

    // Loop over all integration points
    for(unsigned int ip = 0; ip < intPoints.GetSize(); ip++)
      {
      // Calculate for each integration point the LocPointMapped

      lp.Set( intPoints[ip], esm, weights[ip] );

      // Call the CalcBMat()-method
      bdb->GetBOp()->CalcOpMat(b_mat, lp, ptFe);
      //LOG_DBG3(em) << "CMDF: ip=" << ip << " w=" << weights[ip] << " b_mat=" << b_mat.ToString(0, false);
      //std::cout << "b_mat= " << b_mat << std::endl;

      // traverse nodes
      form->GetFeSpace1()->GetElemEqns(eqn, elem);
      assert(eqn.GetSize() == elem->connect.GetSize());
      assert(b_mat.GetNumCols() == elem->connect.GetSize());
      for(unsigned int n = 0; n < eqn.GetSize(); n++)
      {
        appendix[eqn[n]] += b_mat[x_comp ? 0 : 1][n];
        std::cout << "appendix= " << appendix[eqn[n]] << std::endl;
        de.AddGradient(f, appendix[eqn[n]]); // is it correct???
      }
      }

    } // end loop elems*/


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
