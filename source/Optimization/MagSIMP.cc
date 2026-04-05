#include <assert.h>
#include <stddef.h>
#include <cmath>
#include <map>
#include <ostream>
#include <string>

#include "Optimization/MagSIMP.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/TransferFunction.hh"
#include "Optimization/Excitation.hh"
#include "Domain/CoefFunction/CoefFunctionConst.hh"
#include "Domain/CoefFunction/CoefFunctionOpt.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Driver/Assemble.hh"
#include "Forms/BiLinForms/BDBInt.hh"
#include "FeBasis/H1/H1Elems.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/LinForms/LinearForm.hh"
#include "Utils/ToolsFull.hh"
#include "Domain/Domain.hh"
#include "PDE/MagneticPDE.hh"
#include "PDE/MagEdgePDE.hh"

namespace CoupledField {
class DenseMatrix;
}  // namespace CoupledField

using namespace CoupledField;

DEFINE_LOG(ms, "magSimp")

double MagSIMP::nu_0 = 1.0/(4 * M_PI * 1e-7);


MagSIMP::MagSIMP()
{
  assert(close(nu_0,1.0/(4 * M_PI * 1e-7)));

  sel_x_.Resize(domain->GetDim(), domain->GetDim());
  sel_x_.Init();
  sel_x_[0][0] = 1;

  sel_y_.Resize(sel_x_.GetNumRows(), sel_x_.GetNumCols());
  sel_y_.Init();
  sel_y_[1][1] = 1;

  sel_xy_.Resize(sel_x_.GetNumRows(), sel_x_.GetNumCols());
  sel_xy_.Init();
  sel_xy_[0][0] = 1;
  sel_xy_[1][1] = 1;

  lin_nu_r_.Resize(domain->GetGrid()->GetNumRegions() + 1, -1.0);
}

MagSIMP::~MagSIMP()
{
}

void MagSIMP::PostInit()
{
  SIMP::PostInit();

  assert(manager.context.GetSize() == 1);

  magRHS.Init<double>(design, App::MAG);

  nonlin_.Resize(domain->GetGrid()->GetNumRegions() + 1, false);

  PtrParamNode pn = optInfoNode->Get(ParamNode::HEADER)->Get("magOptRegions");

  for(unsigned int i = 0; i < nonlin_.GetSize(); i++)
  {
    nonlin_[i] = context->pde->GetAssemble()->GetBiLinForm("CurlCurlIntegrator-NL", i, NULL, NULL, true) != NULL;
    bool opt = design->Contains(i, true); // with pseudo-design
    bool pseudo = opt && !design->Contains(i, false);
    if(opt)
    {
      PtrParamNode pnr = pn->Get("region", ParamNode::APPEND);
      pnr->Get("name")->SetValue(domain->GetGrid()->GetRegion().ToString(i));
      pnr->Get("id")->SetValue(i);
      pnr->Get("analysis")->SetValue(nonlin_[i] ? "nonlinear" : "linear");
      pnr->Get("pseudo_design")->SetValue(pseudo);
    }
  }
}

bool MagSIMP::FillRealAdjointRHS(Excitation& excite, Function* f, Vector<double>& rhs)
{
  switch(f->GetType())
  {
  case Function::SQR_MAG_FLUX_DENS_X:
  case Function::SQR_MAG_FLUX_DENS_Y:
  case Function::SQR_MAG_FLUX_DENS_RZ:
    CalcMagFluxAdjRHS(excite, f, rhs);
    return true;
    break;
  case Function::LOSS_MAG_FLUX_RZ:
    CalcMagFluxLossesAdjRHS(excite, f, rhs);
    return true;
    break;
  case Function::MAG_COUPLING:
    CalcCouplingAdjRealRHS(excite, f, rhs);
    return true;
    break;

  default:
    break;
  }
  return false;
}

bool MagSIMP::FillComplexAdjointRHS(Excitation& excite, Function* f, Vector<Complex>& rhs)
{
  switch(f->GetType())
  {
  case Function::SQR_MAG_FLUX_DENS_X:
  case Function::SQR_MAG_FLUX_DENS_Y:
  case Function::SQR_MAG_FLUX_DENS_RZ:
  case Function::LOSS_MAG_FLUX_RZ:
    EXCEPTION("Complex RHS for MagFluxDensity not yet implemented!");
    break;
  case Function::MAG_COUPLING:
    CalcCouplingAdjComplexRHS(excite, f, rhs);
    return true;
    break;

  default:
    break;
  }
  return false;
}

double MagSIMP::CalcRelactivity(const Elem* elem, UInt dim)
{
  assert(manager.context.GetSize() == 1); // otherwise we would need it
  OptimizationMaterial* mat = context->mat;
  CoefFunctionOpt* coef = mat->GetMatCoef(mat->stiff, elem->regionId);

  return ExtractRelactivity(coef->orgMat.get(), nonlin_[elem->regionId] ? elem : NULL, dim);
}

double MagSIMP::ExtractRelactivity(CoefFunction* org_mat, const Elem* elem, UInt dim)
{
  double nu_0_nu_r = -1; // in the end we only need the scalar

  if(elem == NULL)
  {
    // in the linear 2D case we have BDBInt with a diagonal dim x dim tensor with identical coefficients, in 3D it is a BBInt with a scalar value
    CoefFunctionConst<double>* cfc = dynamic_cast<CoefFunctionConst<double>*>(org_mat);
    assert(cfc);
    if(dim == 2)
    {
      nu_0_nu_r = cfc->GetTensor()[0][0]; // copy constructor :(
    }
    else
    {
      nu_0_nu_r = cfc->GetScalar();
    }
  }
  else
  {
    // we are here in the linear BDBInt case (first time) and nonlinear scalar BBInt case
    LocPointMapped lpm;
    shared_ptr<ElemShapeMap> esm = domain->GetGrid()->GetElemShapeMap(elem);
    LocPoint lp = Elem::GetShape( Elem::GetShapeType( elem->type) ).midPointCoord;
    lpm.Set(lp, esm); // element constant we need no weight

    if(org_mat->GetDimType() == CoefFunction::SCALAR)
      org_mat->GetScalar(nu_0_nu_r, lpm);
    else
    {
      assert(org_mat->GetDimType() == CoefFunction::TENSOR);
      Matrix<double> mat;
      org_mat->GetTensor(mat, lpm); // TODO: change this to a reference case
      nu_0_nu_r = mat[0][0];
    }
  }

  return nu_0_nu_r/ nu_0;
}

double MagSIMP::CalcFunction(Excitation& excite, Function* f, bool derivative)
{
  switch(f->GetType())
  {
  case Function::SQR_MAG_FLUX_DENS_X:
  case Function::SQR_MAG_FLUX_DENS_Y:
  case Function::SQR_MAG_FLUX_DENS_RZ:
    if(!derivative)
      return CalcMagFluxDensity(excite, f);
    else
    {
      CalcMagFluxDensGradient(excite, f);
      return 0.0;
    }
  case Function::LOSS_MAG_FLUX_RZ:
    if(!derivative)
      return CalcMagFluxDensityLosses(excite, f);
    else
    {
      CalcMagFluxDensGradientLosses(excite, f);
      return 0.0;
    }

  case Function::MAG_COUPLING:
    if(!derivative)
    {
      if(f->ctxt->IsComplex())
        return CalcMagCouplingComplex(excite, f);
      else
        return CalcMagCouplingReal(excite, f);
    }
    else
    {
      TransferFunction* tf = design->GetTransferFunction(DesignElement::Default(f->ctxt), TransferFunction::Default(f->ctxt), true, true);//TODO
      CalcCouplingGradient(excite, f, tf);
      return 0.0;
    }
  default: // return below as we don't implement
    break;
  }
  return SIMP::CalcFunction(excite, f, derivative);
}

double MagSIMP::CalcMagFluxDensity(Excitation& excite, Function* f)
{
  if(f->region == ALL_REGIONS)
    throw Exception("For function " + f->ToString() + " the attribute 'region' is mandatory");

  // Calculation of the normed squared directional magnetic flux density
  // J = sum_e 1/N B^T s B
  // s = (0 1)^T or (1 0)^T
  // B = curl A = M*A where M = the discrete curl operator curl N (=b_mat)
  // (B^T s B) is evaluated at the integration points j:
  // sum_j jacDet * w(j) * (M_j*A_e)^T * s * (M_j*A_e)

  // we need a lot of similar stuff as in BDBInt::CalcElementMatrix().
  // get the form first
  BaseBDBInt* bdb = dynamic_cast<BaseBDBInt*>(context->pde->GetAssemble()->GetBiLinForm("CurlCurlIntegrator", f->region, context->pde)->GetIntegrator());
  assert(bdb != NULL);
  LOG_DBG2(ms) << "CMFD enter: bdb=" << bdb->GetName();

  // annoying entity iterator got hold the elem
  ElemList el(domain->GetGrid());

  // the stored element solution vector
  StdVector<SingleVector*>& sol = forward.Get(excite)->elem[App::MAG];

  Matrix<double> M; // this holds the curl-operator for the whole element at an integration point. For 2D rectangular case it shall be 2 rows, 4 columns
  Vector<double> flux_dens; // M*a or in magnetic symbols B = curl A = curl N * A = M * A
  const Matrix<double>& S = GetSelectionMatrix(f);
  Vector<double> S_flux_dens; // S*flux_dens = vector of dim

  double result = 0;

  StdVector<LocPoint> intPoints; // Get integration Points
  LocPointMapped lp;
  StdVector<double> weights;
  IntegOrder order;
  IntScheme::IntegMethod method = IntScheme::UNDEFINED;
  // calculate volume of whole optimization domain by adding up element volumes
  double vol; // element volume
  double volume = 0; // accumulated volume whole optimization domain

  assert(!f->elements.IsEmpty());
  assert(f->region != ALL_REGIONS);
  for(unsigned int e = 0; e < f->elements.GetSize(); e++)
  {
    DesignElement* de = f->elements[e];

    // element solution
    LOG_DBG3(ms) << "CMFD ind=" << de->GetElementSolutionIndex() << " sol=" << sol[de->GetElementSolutionIndex()];
    Vector<double>* vec = dynamic_cast<Vector<double>*>(sol[de->GetElementSolutionIndex()]);
    LOG_DBG3(ms) << "CMFD e=" << e << " el=" << de->elem->elemNum << " esi=" << de->GetElementSolutionIndex();
    assert(vec != NULL);
    const Vector<double>& a = *vec; // a = the vector potential in the element
    LOG_DBG3(ms) << "CMFD e=" << e << " el=" << de->elem->elemNum << " esi=" << de->GetElementSolutionIndex() << " a=" << a.ToString() << " -> " << result;

    // prepare to get the curl operator
    el.SetElement(de->elem);

    BaseFE* ptFe = bdb->GetFeSpace1()->GetFe(el.GetIterator(), method, order );

    bdb->GetIntScheme()->GetIntPoints(Elem::GetShapeType(de->elem->type), method, order, intPoints, weights );
    assert(method != IntScheme::UNDEFINED);
    assert(!intPoints.IsEmpty());
    // Get shape map from grid
    shared_ptr<ElemShapeMap> esm = domain->GetGrid()->GetElemShapeMap(de->elem);

    vol = esm->CalcVolume();

    LOG_DBG2(ms) << "CMFD i=" << e << " el=" << de->elem->elemNum << " method=" << method << " order=" << order.ToString() << " iP=" << intPoints.ToString() << " v=" << vol;
    // add element volume to volume of whole domain
    volume += vol;
    LOG_DBG3(ms) << "CMFD accumulated volume =" << volume;

    double el_val = 0.0;

    for(unsigned int ip = 0; ip < intPoints.GetSize(); ip++)
    {
      // Calculate for each integration point the LocPointMapped
      lp.Set(intPoints[ip], esm, weights[ip]);

      // Call the CalcBMat()-method
      assert(bdb->GetBOp());
      bdb->GetBOp()->CalcOpMat(M, lp, ptFe);
      assert(M.GetNumCols() == a.GetSize());
      assert(M.GetNumRows() == domain->GetDim());
      LOG_DBG3(ms) << "CMFD: e= " << e << " ip=" << ip << "/(" << intPoints[ip].coord.ToString() << ") w=" << weights[ip] << " jacDet=" << lp.jacDet << " M_" << ip << "=" << M.ToString();

      // flux_denx = M * a
      flux_dens.Resize(dim);
      M.Mult(a, flux_dens);
      assert(flux_dens.GetSize() == dim);

      // S_flux_dens = S * flux_dens
      S_flux_dens.Resize(dim);
      S.Mult(flux_dens, S_flux_dens);

      if(dim == 2)
      {
        // scalar = flux_dens * S_flux_dens
        el_val += weights[ip] * lp.jacDet * S_flux_dens.Inner(flux_dens);
      }
      else
      {
        el_val += weights[ip] * lp.jacDet * S_flux_dens.Inner(flux_dens) * 2; // * 2 because of edge element?
      }

      LOG_DBG3(ms) << "CMFD: e= " << e << " flux_dens=" << flux_dens.ToString() << " Sfd=" << S_flux_dens << " inner=" << S_flux_dens.Inner(flux_dens) << " el -> " << el_val;
    } // end ip

    result += el_val;
  } // end loop elems
  if(dim == 3)
  {
    result /= intPoints.GetSize();
  }

  // set optimization volume if not already set
  if (opt_vol_ == -1)
  {
    opt_vol_ = volume;
    LOG_DBG2(ms) << "CMFD: calculated volume =" << opt_vol_;
  }

  LOG_DBG(ms) << "CMFD: exit normed -> " << result;
  return result;
}
double MagSIMP::CalcMagFluxDensityLosses(Excitation& excite, Function* f)
{
  double result = CalcMagFluxDensity(excite, f);
  //TODO Scale result by Volume*density
  return result;
}

void MagSIMP::CalcMagFluxDensGradient(Excitation& excite, Function* f)
{
  assert(f->GetExcitation() != NULL);
  assert(f->GetExcitation()->transform == NULL); // don' do the complicated stuff yet


  // the gradient is < lambda^T, K' * A >
  assert(excite.sequence == f->ctxt->sequence);
  assert(f->ctxt->ToApp() == App::MAG);

  // we use this just as indicator for coil optimization.
  DesignDependentRHS rhs(App::MAG);

  // for the two design case, we loop over both variables. In the single design case this also works
  for(unsigned int d = 0; d < design->design.GetSize(); d++)
  {
    DesignElement::Type dt = design->design[d].design;
    TransferFunction* tf = design->GetTransferFunction(dt, App::MAG, true);

    DesignDependentRHS* ptr_rhs = dt == DesignElement::RHS_DENSITY ? &rhs : NULL;

    double factor = 1.0/ f->elements.GetSize(); // factor for norming the gradient; same as in objective function
    // calc lambda^T *  K' * A -> this already stores the results by AddGradient()!
    CalcU1KU2(tf, adjoint.Get(excite, f)->elem[App::MAG], App::MAG, forward.Get(excite)->elem[App::MAG], ptr_rhs, factor, STANDARD, f);
  }
}

void MagSIMP::CalcMagFluxDensGradientLosses(Excitation& excite, Function* f)
{
  assert(f->GetExcitation() != NULL);
  assert(f->GetExcitation()->transform == NULL); // don' do the complicated stuff yet


  // the gradient is < lambda^T, K' * A >
  assert(excite.sequence == f->ctxt->sequence);
  assert(f->ctxt->ToApp() == App::MAG);

  // we use this just as indicator for coil optimization.
  DesignDependentRHS rhs(App::MAG);

  // for the two design case, we loop over both variables. In the single design case this also works
  for(unsigned int d = 0; d < design->design.GetSize(); d++)
  {
    DesignElement::Type dt = design->design[d].design;
    TransferFunction* tf = design->GetTransferFunction(dt, App::MAG, true);

    DesignDependentRHS* ptr_rhs = dt == DesignElement::RHS_DENSITY ? &rhs : NULL;

    double factor = 1.0/ f->elements.GetSize(); // factor for norming the gradient; same as in objective function
    // calc lambda^T *  K' * A -> this already stores the results by AddGradient()!
    CalcU1KU2(tf, adjoint.Get(excite, f)->elem[App::MAG], App::MAG, forward.Get(excite)->elem[App::MAG], ptr_rhs, factor, STANDARD, f);
  }
  //TODO add design dependent term
}

void MagSIMP::CalcN(LinearFormContext* form, Vector<double>& N)
{
  // Calculates the vector N to perform the numerical integration over the region encoded in form by <N, A>
  assert(N.GetSize() > 0); // needs to be set to full state solution size such that we can index the entries by the equation map

  shared_ptr<EntityList> el = form->GetEntities();
  assert(el->GetRegion() != NO_REGION_ID);
  assert(el->GetRegion() >= 0);
  LOG_DBG3(ms) << "CN: ft=" << el->GetType() << " r=" << el->GetRegion() << " fn="  << form->ToString();

  // get bilinear form to gain fe-space from it
  BaseBDBInt* bdb = dynamic_cast<BaseBDBInt*>(context->pde->GetAssemble()->GetBiLinForm("CurlCurlIntegrator", el->GetRegion(), context->pde)->GetIntegrator());
  assert(bdb != NULL);
  LOG_DBG3(ms) << "CN: bdb=" << bdb->GetName();

  StdVector<LocPoint> intPoints; // Get integration Points
  LocPointMapped lpm;
  StdVector<double> weights;
  Vector<double> shapes; // shape function coefficients for specific ip
  IntegOrder order;
  IntScheme::IntegMethod method = IntScheme::UNDEFINED;
  // this gets the equation numbers for the element
  StdVector<int> eqn;

  for(EntityIterator iter = el->GetIterator(); !iter.IsEnd(); iter++)
  {
    // from the element we get the local shape functions
    BaseFE* bfe = bdb->GetFeSpace1()->GetFe(iter, method, order );
    FeH1* h1 = dynamic_cast<FeH1*>(bfe);
    assert(h1 != NULL);

    // the shape map knows about the real element size
    shared_ptr<ElemShapeMap> esm = domain->GetGrid()->GetElemShapeMap(iter.GetElem());

    // the equations assigned to the element
    bdb->GetFeSpace1()->GetElemEqns(eqn, iter.GetElem());

    // set integration points and weights
    bdb->GetIntScheme()->GetIntPoints(Elem::GetShapeType(iter.GetElem()->type), method, order, intPoints, weights );
    assert(method != IntScheme::UNDEFINED);
    assert(!intPoints.IsEmpty());

    LOG_DBG3(ms) << "CN: e=" << iter.GetElem()->elemNum << " w=" << weights.ToString()  << " o=" << order.ToString() << " m=" << method << " ns=" << bfe->GetNumFncs() << " fe=" << bfe->FeType() << " eqn=" << eqn.ToString();

    for(unsigned int ip = 0; ip < intPoints.GetSize(); ip++)
    {
      // Calculate for each integration point the determinant of the Jacobian
      lpm.Set(intPoints[ip], esm, weights[ip]);

      // get the shape functions of the ip, same values for all ips in one element but permuted
      h1->GetShFnc(shapes, lpm.lp, iter.GetElem());
//      LOG_DBG3(ms) << "CN: e=" << iter.GetElem()->elemNum << " ip=" << ip << "=" << intPoints[ip].coord.ToString() << " J=" << lpm.jacDet << " s=" << shapes.ToString();
      assert(shapes.GetSize() == eqn.GetSize());
      for(unsigned int s = 0; s < shapes.GetSize(); s++)
      {
        double v = weights[ip] * lpm.jacDet * shapes[s];
        // the equation number is 1 based with 0 indicating HDBC and constrained nodes for negative indices. The equation index is 0-based!
        int eqn_nbr = eqn[s];
        if(eqn_nbr <= 0)
        {
          LOG_DBG3(ms) << "CN: s=" << s << " eqn_nbr=" << eqn_nbr << " -> skip RHS node";
        }
        else
        {
          unsigned int eqn_idx = eqn_nbr-1;
//          LOG_DBG3(ms) << "CN: s=" << s << " eqn_idx=" << eqn_idx << " N[" << eqn_idx << "] = " << N[eqn_idx] << " + " << v << " -> " << (N[eqn_idx] + v);
          N[eqn_idx] += v;
         }
      } // end shape fcts
    } // end ip
  } // end elem
}

double MagSIMP::CalcMagCouplingReal(Excitation& excite, Function* f)
{
  if(GetMultipleExcitation()->excitations.GetSize() < 2)
    throw Exception("'magCoupling' requires two coils and enabled multiple_excitations");

  // cfs makes two coupling functions for two excitations with a weight of 0.5 each.
  // we return value 0 for the first excitation and 2*coupling for the second excitation
  assert(GetMultipleExcitation()->excitations.GetSize() >= 2);
  if(excite.index != 1)
  {
    LOG_DBG(ms) << "CMC: excitation index is " << excite.index <<" not 1; nothing to do here; returning 0.0";
    return 0.0;
  }
  assert(excite.index == 1);
  Excitation& excite_A = GetMultipleExcitation()->excitations[0];
  Excitation& excite_B = GetMultipleExcitation()->excitations[1];
  // Coupling can only be calculated for exactly two coils
  assert(excite_A.index == 0 && excite_B.index == 1);
  // we need both forms to get the different regions
  StdVector<LinearFormContext*>& forms_A = excite_A.forms;
  StdVector<LinearFormContext*>& forms_B = excite_B.forms;
  assert((forms_A.GetSize() == 1) && (forms_B.GetSize() == 1));
  LinearFormContext* form_A = forms_A[0];
  LinearFormContext* form_B = forms_B[0];
  // get the two As
  const Vector<double>& A_a = forward.Get(excite_A,NULL)->GetRealVector(StateSolution::RAW_VECTOR);
  const Vector<double>& A_b = forward.Get(excite_B,NULL)->GetRealVector(StateSolution::RAW_VECTOR);
  assert(A_a.GetSize() > 0);
  assert(A_b.GetSize() > 0);

  LOG_DBG3(ms) << "CMC: size of A_a = " << A_a.GetSize();
  LOG_DBG3(ms) << "CMC: size of A_b = " << A_b.GetSize();

  Vector<double> N1(A_a.GetSize());
  Vector<double> N2(A_a.GetSize());

  CalcN(form_A, N1);
  CalcN(form_B, N2);

  double sp_N1_Aa = N1.Inner(A_a);
  double sp_N1_Ab = N1.Inner(A_b);
  double sp_N2_Ab = N2.Inner(A_b);

  LOG_DBG3(ms) << "CMC: <N1, A_a>=" << sp_N1_Aa;
  LOG_DBG3(ms) << "CMC: <N1, A_b>=" << sp_N1_Ab;
  LOG_DBG3(ms) << "CMC: <N2, A_b>=" << sp_N2_Ab;

  double k = (pow(sp_N1_Ab, 2))/(sp_N1_Aa*sp_N2_Ab);
  LOG_DBG(ms) << "CMC: Coupling = " << k;
  // return 2-times the function value since we have two excitations and it will be halved afterwards
  return k;
  }

double MagSIMP::CalcMagCouplingComplex(Excitation& excite, Function* f)
{
  /* for the complex coupling we take the amplitude of A
   * k = \frac{(N_1^T (A_R^B .* A_R^B) + N_1^T (A_I^B .* A_I^B))^2} {(N_1^T (A_R^A .* A_R^A) + N_1^T (A_I^A .* A_I^A))*(N_2^T (A_R^B .* A_R^B) + N_2^T (A_I^B .* A_I^B))}
   */
  //  if(GetMultipleExcitation()->excitations.GetSize() != 2)
  //    throw Exception("'magCoupling' requires two coils and enabled multiple_excitations");

  // cfs makes two coupling functions for two excitations with a weight of 0.5 each.
  // we return value 0 for the first excitation and 2*coupling for the second excitation
  //  assert(GetMultipleExcitation()->excitations.GetSize() == 2);
  if(excite.index == 0)
  {
    LOG_DBG2(ms) << "CMC: excitation index is 0; nothing to do here; returning 0.0";
    return 0.0;
  }
  assert(excite.index == 1);
  Excitation& excite_A = GetMultipleExcitation()->excitations[0];
  Excitation& excite_B = GetMultipleExcitation()->excitations[1];
  // Coupling can only be calculated for exactly two coils
  assert(excite_A.index == 0 && excite_B.index == 1);
  // we need both forms to get the different regions
  StdVector<LinearFormContext*>& forms_A = excite_A.forms;
  StdVector<LinearFormContext*>& forms_B = excite_B.forms;
  assert((forms_A.GetSize() == 1) && (forms_B.GetSize() == 1));
  LinearFormContext* form_A = forms_A[0];
  LinearFormContext* form_B = forms_B[0];
  // get the two As
  const Vector<Complex>& A_a = forward.Get(excite_A,NULL)->GetComplexVector(StateSolution::RAW_VECTOR);
  const Vector<Complex>& A_b = forward.Get(excite_B,NULL)->GetComplexVector(StateSolution::RAW_VECTOR);
  assert(A_a.GetSize() > 0);
  assert(A_b.GetSize() > 0);

  LOG_DBG3(ms) << "CMC: size of A_a = " << A_a.GetSize();
  LOG_DBG3(ms) << "CMC: size of A_b = " << A_b.GetSize();

  Vector<double> N1(A_a.GetSize());
  Vector<double> N2(A_a.GetSize());

  CalcN(form_A, N1);
  CalcN(form_B, N2);

  // N1_ABR = <N1, sqrt(real(AB).^2 + imag(AB).^2)>
  double N1_AA = InnerHelper(N1, A_a);

  double N1_AB = InnerHelper(N1, A_b);

  double N2_AB = InnerHelper(N2, A_b);

  LOG_DBG(ms) << "CMC: N1_AA = " << N1_AA;
  LOG_DBG(ms) << "CMC: N1_AB = " << N1_AB;
  LOG_DBG(ms) << "CMC: N2_AB = " << N2_AB;
  double k = pow(N1_AB, 2)/(N1_AA*N2_AB);
  LOG_DBG(ms) << "CMC: Coupling = " << k;
  // return 2-times the function value since we have two excitations and it will be halved afterwards
  return 2*k;
}

void MagSIMP::CalcCouplingGradient(Excitation& excite, Function* f, TransferFunction* tf)
{
  //TODO
  LOG_DBG2(ms) << "CCG: excitation: " << excite.index;
  CalcU1KU2(tf, adjoint.Get(excite, f)->elem[App::MAG], App::MAG, forward.Get(excite)->elem[App::MAG], NULL, 1, STANDARD, f);
}

void MagSIMP::CalcMagFluxAdjRHS(Excitation& excite, Function* f, Vector<double>& out)
{
  // B = 2 rows, 4 columns, A = 4 rows, 1 column, D = 2 row, 2 col,
  // D is the selection vector for x or y component [0 0; 0 1] or [1 0; 0 0]
  // d(<B*A,D*B*A>)/dA = -2*(B^T*D*B)*A = -2*M*A -> 4 rows, 1 col
  // do this for all elements A and add with factor 1/N to rhs
  assert(out.GetSize() == 0);
  const Vector<double>& stateSol = forward.Get(excite,NULL)->GetRealVector(StateSolution::RAW_VECTOR);
  out.Resize(stateSol.GetSize(),0.0);

  MagMat* mag = dynamic_cast<MagMat*>(f->ctxt->mat);
  SinglePDE* pde = f->ctxt->pde;

  SolutionType solt = MAG_POTENTIAL;
  shared_ptr<BaseFeFunction> fe = context->pde->GetFeFunction(solt);

  // We use BDBInt to compute M, but in the nonlinear case we have only BBInt stored in assemble
  // therefore we construct a own one
    // new BDBInt<>(new CurlOperator<FeH1,2,Double>(), curCoef,factor, updatedGeo_); in MagneticPDE.cc
  // take the curl operator from the BDBInt/BBInt for our region
  assert(f->region != ALL_REGIONS);
  BiLinFormContext* c = pde->GetAssemble()->GetBiLinForm(mag->stiff.integrator, f->region);
  BaseBDBInt* c_bdb = dynamic_cast<BaseBDBInt*>(c->GetIntegrator());
  assert(c_bdb);

  // this selects the orientation, similar to the dmat in BDBInt but w/o material
  const Matrix<double>& S =GetSelectionMatrix(f);
  shared_ptr<CoefFunctionConst<double> > coef_S(new CoefFunctionConst<double>());
  coef_S->SetTensor(S);

  BaseBDBInt* bdb = new BDBInt<double>(c_bdb->GetBOp()->Clone(), coef_S, 1.0, c_bdb->IsCoordUpdate());
  bdb->SetFeSpace(c_bdb->GetPtrFeSpace1());

  // the stored element solution vector
  StdVector<SingleVector*>& sol = forward.Get(excite)->elem[App::MAG];

  // this gets the equation numbers for the element
  StdVector<int> eqn;
  Matrix<double> M;
  Vector<double> rhs_el;

  ElemList elemList(domain->GetGrid());
  EntityIterator it = elemList.GetIterator();

  // calculate volume of whole optimization domain by adding up element volumes
  double vol; // element volume
  double volume = 0; // accumulated volume whole optimization domain

  for(unsigned int e = 0; e < f->elements.GetSize(); e++)
  {
    DesignElement* de = f->elements[e];
    elemList.SetElement(de->elem);

    // M = B^T S B as in BDBInt::CalcElementMatrix, just S from above instead of D with material and density
    bdb->CalcElementMatrix(M, it, it);
    LOG_DBG3(ms) << "CMFAR e=" << e << " el=" << de->elem->elemNum<< "coef_S= " << coef_S->GetTensor().ToString() << " M=" << M.ToString();

    // now the part -2 * M * A
    Vector<double>* vec = dynamic_cast<Vector<double>*>(sol[de->GetElementSolutionIndex()]);
    LOG_DBG3(ms) << "CMFAR e=" << e << " esi=" << de->GetElementSolutionIndex() << " nodal values=" << vec->ToString();
    assert(vec != NULL);
    Vector<double>& a = *vec; // a = the vector potential in the element
    if(dim == 3)
    {
      a /= 2; // /2 because of edge element?
    }

    rhs_el.Resize(a.GetSize());
    M.Mult(a, rhs_el);
    rhs_el *= -2;

    LOG_DBG3(ms) << "CMFAR e=" << e << " a=" << a.ToString() << " -> -2 M*a: r=" << rhs_el.ToString();

    c->GetIntegrator()->GetFeSpace1()->GetElemEqns(eqn, de->elem);
    LOG_DBG3(ms) << "CMFAR e=" << e << " eqn=" << eqn.ToString();
    assert(rhs_el.GetSize() == eqn.GetSize());

    // add element volume to volume of whole domain
    // Get shape map from grid
    shared_ptr<ElemShapeMap> esm = domain->GetGrid()->GetElemShapeMap(de->elem);
    vol = esm->CalcVolume();
    volume += vol;
    LOG_DBG3(ms) << "CMFAR accumulated volume =" << volume;
    for(unsigned int n = 0; n < eqn.GetSize(); n++)
    {
      // the equation number is 1 based with 0 indicating HDBC and constrained nodes for negative indices. The equation index is 0-based!
      int eqn_nbr = eqn[n];

      if(eqn_nbr <= 0) {
        LOG_DBG2(ms) << "CMFAR: n=" << n << " eqn_nbr=" << eqn_nbr << " -> skip RHS node";
      }
      else
      {
        unsigned int eqn_idx = eqn_nbr-1;

        out[eqn_idx] += rhs_el[n];
        LOG_DBG2(ms) << "CMFAR: n=" << n << " eqn_idx=" << eqn_idx << " normed_rhs_el[n]= " << rhs_el[n] << " out[eqn_idx]=" << out[eqn_idx];
      }
    }
  } // end loop elements

  // set optimization volume if not already set
  if (opt_vol_ == -1)
  {
    opt_vol_ = volume;
    LOG_DBG2(ms) << "CMFAR: calculated volume =" << opt_vol_;
  }

  delete bdb;
}

void MagSIMP::CalcMagFluxLossesAdjRHS(Excitation& excite, Function* f, Vector<double>& out)
{
  CalcMagFluxAdjRHS(excite, f, out);
  //TODO Scale out
}

void MagSIMP::CalcCouplingAdjRealRHS(Excitation& excite, Function* f, Vector<double>& out)
{
  /* The right hand sides look as follows (real case):
   * case A:
   * (<N_1, A_B>^2)/(<N_1, A_A>^2 <N_2, A_B>) * N_1^T
   * case B:
   * <N_1, A_B>/(<N_1, A_A><N_2, A_B>) * N_1^T - (<N_1, A_B>)^2/(<N_1, A_A><N_2, A_B^2>) * N_2^T
   */
  // same as in CalcCoupling
  if(GetMultipleExcitation()->excitations.GetSize() < 2)
  {
    throw Exception("'magCoupling' requires two coils and enabled multiple_excitations");
  }
  assert(GetMultipleExcitation()->excitations.GetSize() >= 2);
  Excitation& excite_A = GetMultipleExcitation()->excitations[0];
  Excitation& excite_B = GetMultipleExcitation()->excitations[1];
  // Coupling can only be calculated for exactly two coils
  assert(excite_A.index == 0 && excite_B.index == 1);
  // we need both forms to get the different regions
  StdVector<LinearFormContext*>& forms_A = excite_A.forms;
  StdVector<LinearFormContext*>& forms_B = excite_B.forms;
  assert((forms_A.GetSize() == 1) && (forms_B.GetSize() == 1));
  LinearFormContext* form_A = forms_A[0];
  LinearFormContext* form_B = forms_B[0];

  // get the two As
  const Vector<double>& A_a = forward.Get(excite_A,NULL)->GetRealVector(StateSolution::RAW_VECTOR);
  const Vector<double>& A_b = forward.Get(excite_B,NULL)->GetRealVector(StateSolution::RAW_VECTOR);
  assert(A_a.GetSize() > 0);
  assert(A_b.GetSize() > 0);

  LOG_DBG3(ms) << "CCAR: size of A_a = " << A_a.GetSize();
  LOG_DBG3(ms) << "CCAR: size of A_b = " << A_b.GetSize();

  Vector<double> N1(A_a.GetSize());
  Vector<double> N2(A_a.GetSize());

  CalcN(form_A, N1);
  CalcN(form_B, N2);
  LOG_DBG3(ms) << "CMC: N1 = " << N1.ToString();
  LOG_DBG3(ms) << "CMC: N2 = " << N2.ToString();

  double sp_N1_Aa = Inner(N1, A_a);
  double sp_N1_Ab = Inner(N1, A_b);
  double sp_N2_Ab = Inner(N2, A_b);

  LOG_DBG2(ms) << "CCAR: <N1, A_a>=" << sp_N1_Aa;
  LOG_DBG2(ms) << "CCAR: <N1, A_b>=" << sp_N1_Ab;
  LOG_DBG2(ms) << "CCAR: <N2, A_b>=" << sp_N2_Ab;

  out.Resize(A_a.GetSize(),0.0);
  if (excite.index == 0) {
    //calc rhs first case
    const double factor = pow(sp_N1_Ab, 2) / ( pow(sp_N1_Aa, 2) * sp_N2_Ab );
    out.Set(factor, N1); // out = factor * N1
    LOG_DBG2(ms) << "CCAR: first excitation. f=" << factor << " |out|=" << out.NormL2();
  } else if (excite.index == 1) {
    //calc rhs second case
    const double factor_N1 = -(2*sp_N1_Ab)/(sp_N1_Aa*sp_N2_Ab);
    const double factor_N2 = (pow(sp_N1_Ab, 2))/(sp_N1_Aa*pow(sp_N2_Ab, 2));
    out.Add(factor_N1, N1, factor_N2, N2); // out = -factor_N1 * N1 + factor_N2 * N2
    LOG_DBG2(ms) << "CCAR: second excitation fN1=" << factor_N1 << " fN2=" <<factor_N2 << " |out|=" << out.NormL2();
  } else {
    EXCEPTION("There should only be two excitations");
  }
}

void MagSIMP::CalcCouplingAdjComplexRHS(Excitation& excite, Function* f, Vector<Complex>& out)
{
  /* The right hand sides look as follows (complex case):
   * In general:
   * -1/2 * ( \frac{\partial J}{\partial u_R} - j*\frac{\partial J}{\partial u_I} )^T
   * case A:
   * factor =  -(N1_ABR+N1_ABI)^2/((N1_AAR+N1_AAI)^2 * (N2_ABR+N2_ABI))
   * - factor * N_1 .* A_R^A + j * 1/2 * factor * N_1 .* A_I^A
   * case B:
   * factor =
   */
  // same as in CalcCoupling

  LOG_DBG2(ms) << "CCAR: num excitations: " << GetMultipleExcitation()->excitations.GetSize();
  assert(f->ctxt->IsComplex());
  Excitation& excite_A = GetMultipleExcitation()->excitations[0];
  LOG_DBG2(ms) << "CCAR: excA label: " << excite_A.label <<" frequency "<<excite_A.frequency << " index " << excite_A.index;
  Excitation& excite_B = GetMultipleExcitation()->excitations[1];
  // Coupling can only be calculated for exactly two coils
  //assert(excite_A.index == 0 && excite_B.index == 1);
  // we need both forms to get the different regions
  StdVector<LinearFormContext*>& forms_A = excite_A.forms;
  LOG_DBG2(ms) << "CCAR: num forms in A: " << forms_A.GetSize();
  StdVector<LinearFormContext*>& forms_B = excite_B.forms;
  LOG_DBG2(ms) << "CCAR: num forms in B: " << forms_B.GetSize();
  assert((forms_A.GetSize() == 1) && (forms_B.GetSize() == 1));
  LinearFormContext* form_A = forms_A[0];
  LinearFormContext* form_B = forms_B[0];

  // get the two As
  const Vector<Complex>& A_a = forward.Get(excite_A,NULL)->GetComplexVector(StateSolution::RAW_VECTOR);
  const Vector<Complex>& A_b = forward.Get(excite_B,NULL)->GetComplexVector(StateSolution::RAW_VECTOR);

  LOG_DBG3(ms) << "CCAR: size of A_a = " << A_a.GetSize();
  LOG_DBG3(ms) << "CCAR: size of A_b = " << A_b.GetSize();
  assert(A_a.GetSize() > 0);
  assert(A_b.GetSize() > 0);

  Vector<double> N1(A_a.GetSize());
  Vector<double> N2(A_a.GetSize());

  CalcN(form_A, N1);
  CalcN(form_B, N2);
  LOG_DBG3(ms) << "CCAR: N1 = " << N1.ToString();
  LOG_DBG3(ms) << "CCAR: N2 = " << N2.ToString();

  // N1_AB = <N1, sqrt(real(AB).^2 + imag(AB).^2)>
  double N1_AA = InnerHelper(N1, A_a);
  double N1_AB = InnerHelper(N1, A_b);
  double N2_AB = InnerHelper(N2, A_b);

  Complex j = Complex(0,1);
  out.Resize(A_a.GetSize(),0.0);
  if (excite.index == 0) {
    //calc rhs first case
    const Complex factor = -(N1_AB*N1_AB)/((N1_AA*N1_AA) * N2_AB);
    LOG_DBG2(ms) << "CCAR: N1_AB " <<N1_AB;
    LOG_DBG2(ms) << "CCAR: N1_AA " <<N1_AA;
    LOG_DBG2(ms) << "CCAR: pow(N1_AB, 2) " <<pow(N1_AB, 2);
    LOG_DBG2(ms) << "CCAR: pow(N1_AA, 2) " <<pow(N1_AA, 2);
    LOG_DBG2(ms) << "CCAR: factor complex case excitation A =" << factor;
    // out = -0.5factor * N1 *Re(A_a) + j*0.5*factor* N1 * Imag(A_a)
    HadamardHelper(out, -0.5*factor, N1, j*0.5*factor, A_a);
    LOG_DBG2(ms) << "CCAR: first excitation. f=" << factor << " |out|=" << out.NormL2();
  } else if (excite.index == 1) {
    //calc rhs second case
    Complex fac1 = 2*N1_AB/(N1_AA*N2_AB);
    Complex fac2 = -(N1_AB*N1_AB)/(N1_AA*N2_AB*N2_AB);
    LOG_DBG2(ms) << "CCAR: fac N1 complex excitation B =" << fac1;
    LOG_DBG2(ms) << "CCAR: fac N2 complex excitation B =" << fac2;
    // out = -factor_N1 * N1 .* Re(AB) + j*factor_N1 * N1 .* Imag(AB) - 0.5 * factor_N2 * N2 * Re(AB) + 0.5*j*factor_N2 * N2 * Imag(AB)
    HadamardHelper2(out, fac1, N1, fac2, N2, A_b);
    LOG_DBG2(ms) << "CCAR: second excitation fN1=" << fac1 << " fN2=" <<fac2 << " |out|=" << out.NormL2();
  } else {
    EXCEPTION("There should only be two excitations");
  }
}

double MagSIMP::InnerHelper(const Vector<double>& N, const Vector<Complex>& A)
{
  // res = <N, .sqrt(real(A).^2 + imag(A).^2)>
  assert(N.GetSize() == A.GetSize());
  LOG_DBG2(ms) << "IH: N norm:" << N.NormL2();
  LOG_DBG2(ms) << "IH: A norm:" << A.NormL2();

  double sum = 0;

  for(unsigned int i = 0; i < N.GetSize(); i++) {
      sum += N[i] * std::sqrt(std::pow(A[i].real(), 2) +  std::pow(A[i].imag(), 2));
  }
  return sum;
}

void MagSIMP::HadamardHelper(Vector<Complex>& out, Complex factor_N1, const Vector<double>& N1,  Complex factor_N2, const Vector<Complex>& A)
{
  // out = factor_N1 * N1 .* Re(A) + factor_N2 * N1 .* Imag(A)
  assert((N1.GetSize() == A.GetSize()));

  for(unsigned int i = 0; i < N1.GetSize(); i++) {
    double AR = A[i].real();
    double AI = A[i].imag();
    Complex frac = 1.0 / (std::sqrt(AR*AR + AI*AI));
    LOG_DBG3(ms) << "HH: frac: " << frac;
    out[i] = factor_N1 * N1[i] * frac * AR + factor_N2 * N1[i] * frac * AI;
    LOG_DBG3(ms) << "HH: out: " << out[i];
  }
}

void MagSIMP::HadamardHelper2(Vector<Complex>& out, Complex factor_N1, const Vector<double>& N1,  Complex factor_N2, const Vector<double>& N2, const Vector<Complex>& AB)
{
  // out = -factor_N1 * N1 .* Re(AB) + j*factor_N1 * N1 .* Imag(AB) - factor_N2 * N2 * Re(AB) + j*factor_N2 * N2 * Imag(AB)
  Complex j = Complex(0,1);
  assert((N1.GetSize() == N2.GetSize()) && (N1.GetSize() == AB.GetSize()));

  for(unsigned int i = 0; i < N1.GetSize(); i++) {
    double ABR = AB[i].real();
    double ABI = AB[i].imag();
    Complex frac = 1.0 / (std::sqrt(ABR*ABR + ABI*ABI));
    LOG_DBG3(ms) << "HH2: frac: " << frac;
    out[i] = -0.5*factor_N1 * N1[i] * frac * ABR + 0.5*j*factor_N1 *N1[i] * frac * ABI - factor_N2 * 0.5 * N2[i] * frac * ABR + j*factor_N2 * 0.5* N2[i] * frac * ABI;
    LOG_DBG3(ms) << "HH2: out: " << out[i];
  }
}

template <class T1, class T2>
void MagSIMP::SetElementK(Function* f, DesignElement* de, const TransferFunction* tf, App::Type app, DenseMatrix* mat_out, bool derivative, CalcMode calcMode, double ev)
{
  OptimizationMaterial* mat = f->ctxt->mat;
  Matrix<T1>& out = dynamic_cast<Matrix<T1>& >(*mat_out);

  switch(app)
  {
  // u1^T (K' u2 - f') -> find "K'"
  case App::MAG:
  {
    MagMat* mag = dynamic_cast<MagMat*>(mat);
    assert(mag != NULL);

    assert(derivative);
    double d_rho = tf->Derivative(de, DesignElement::SMART, false);

    // element matrix with org material, might be cached -> local_element_cache
    const Matrix<T2>& stiffness = dynamic_cast<const Matrix<T2>& >(mag->Stiffness(de->elem));

    // Overwrite again with the stuff from StateSolution.cc line 353 and 445

    //LOG_DBG3(ms) << "e=" << de->elem->elemNum << " K_0=" << stiffness.ToString();

    double nu_r = GetRelactivity(de->elem, domain->GetDim());
    // simulation: BDB with D=(d 0; 0 d) with d = nu_0*nu_r
    // optimization: d = nu_0 * nu_r * f(rho) - nu_0 * f(rho) + nu_0
    // derivative: d = nu_0 * nu_r * f'(rho) - nu_0 * f'(rho)

    // in difference to elasticity where K' = my(rho)'*K_0
    // we have K' = alpha * K_0 where alpha = (nu_0*nu_r - nu_0) / (nu_0 * nu_r)

    double alpha = (nu_0*nu_r*d_rho - nu_0) / (nu_0 * nu_r);
    //LOG_DBG3(ms) << "e=" << de->elem->elemNum << " "
     //   "nu_0=" << nu_0 << " nu_r=" << nu_r << " rho=" << de->GetDesign(DesignElement::SMART) << " d_rho=" << d_rho << " alpha=" << alpha;

    Assign(out, stiffness, alpha); // out = alpha * stiffness

    //LOG_DBG3(ms) << "out=" << out.ToString();
    break;
  }

  default:
    SIMP::SetElementK(f, de, tf, app, mat_out, derivative, calcMode, ev);
    return; // other cases should be handled in SIMP
  }  // end switch
}


void MagSIMP::SubstractCalcU1KU2RHS(Function* f, TransferFunction* tf, DesignElement* de, DesignDependentRHS* rhs, SingleVector* mat_vec)
{
  assert(f->ctxt->ToApp() == App::MAG);
  assert(rhs != NULL);
  assert(dynamic_cast<MagneticPDE*>(f->ctxt->pde)->DoCoilOptimization());
  assert(f->GetType() == Function::SQR_MAG_FLUX_DENS_X || f->GetType() == Function::SQR_MAG_FLUX_DENS_Y);

  Vector<double>& out = dynamic_cast<Vector<double>& >(*mat_vec);

  const Vector<double>& vec = dynamic_cast<MagMat*>(context->mat)->MagExcitationRHS("CoilIntegrator",de->elem);

  double d_rho = tf->Derivative(de, DesignElement::SMART, false);
  Vector<double> in_out(vec.GetSize());
  in_out.Add(-1 * d_rho, vec);

  out += in_out;

  LOG_DBG3(ms) << "SCUKU: de=" << de->ToString() << " tf=" << tf->ToString() << " d_rho= " << d_rho << " in_out=" << in_out.ToString() << " mat_vec= " << mat_vec->ToString();

}

const Matrix<double>& MagSIMP::GetSelectionMatrix(const Function* f) const
{
  switch(f->GetType())
  {
  case Function::SQR_MAG_FLUX_DENS_X:
    return sel_x_;
  case Function::SQR_MAG_FLUX_DENS_Y:
    return sel_y_;
  case Function::SQR_MAG_FLUX_DENS_RZ:
  case Function::LOSS_MAG_FLUX_RZ:
    return sel_xy_;
  default:
    assert(false);
    return sel_y_;
  }
}
