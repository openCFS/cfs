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
#include "DataInOut/Logging/log.hpp"
#include "Driver/Assemble.hh"
#include "Forms/BiLinForms/BDBInt.hh"
#include "Utils/tools.hh"
#include "Domain/Domain.hh"

namespace CoupledField {
class DenseMatrix;
}  // namespace CoupledField

using namespace CoupledField;

DECLARE_LOG(ms)
DEFINE_LOG(ms, "magSimp")

double MagSIMP::nu_0 = 1.0/(4 * M_PI * 1e-7);

Vector <double> real_nu;

MagSIMP::MagSIMP()
{
  assert(close(nu_0,1.0/(4 * M_PI * 1e-7)));

  sel_x_.Resize(domain->GetGrid()->GetDim(), domain->GetGrid()->GetDim());
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

  nonlin_.Resize(domain->GetGrid()->GetNumRegions() + 1, -1.0);

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
  case Function::MAG_COUPLING:
    CalcCouplingAdjRHS(excite, f, rhs);
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
      TransferFunction* tf = design->GetTransferFunction(DesignElement::Default(f->ctxt), TransferFunction::Default(f->ctxt), true, true);
      CalcMagFluxDensGradient(excite, f, tf);
      return 0.0;
    }
  case Function::MAG_COUPLING:
    if(!derivative)
      return CalcMagCoupling(excite, f);
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
    LOG_DBG3(ms) << "CMFD e=" << e << " el=" << de->elem->elemNum << " esi=" << de->GetElementSolutionIndex() << " a=" << a.ToString(2) << " -> " << result;

    // prepare to get the curl operator
    el.SetElement(de->elem);
    BaseFE* ptFe = bdb->GetFeSpace1()->GetFe(el.GetIterator(), method, order );

    bdb->GetIntScheme()->GetIntPoints(Elem::GetShapeType(de->elem->type), method, order, intPoints, weights );
    assert(method != IntScheme::UNDEFINED);
    assert(!intPoints.IsEmpty());
    // Get shape map from grid
    shared_ptr<ElemShapeMap> esm = domain->GetGrid()->GetElemShapeMap(de->elem);

    vol = esm->CalcVolume();

    LOG_DBG(ms) << "CMFD i=" << e << " el=" << de->elem->elemNum << " method=" << method << " order=" << order.ToString() << " iP=" << intPoints.ToString(2) << " v=" << vol;
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
      assert(M.GetNumRows() == domain->GetGrid()->GetDim());
      LOG_DBG3(ms) << "CMFD: e= " << e << " ip=" << ip << "/(" << intPoints[ip].coord.ToString() << ") w=" << weights[ip] << " jacDet=" << lp.jacDet << " M_" << ip << "=" << M.ToString(2);

      // flux_denx = M * a
      flux_dens.Resize(dim);
      M.Mult(a, flux_dens);
      assert(flux_dens.GetSize() == dim);

      // S_flux_dens = S * flux_dens
      S_flux_dens.Resize(dim);
      S.Mult(flux_dens, S_flux_dens);

      // scalar = flux_dens * S_flux_dens
      el_val += weights[ip] * lp.jacDet * S_flux_dens.Inner(flux_dens);

      LOG_DBG3(ms) << "CMFD: e= " << e << " flux_dens=" << flux_dens.ToString(2) << " Sfd=" << S_flux_dens << " inner=" << S_flux_dens.Inner(flux_dens) << " el -> " << el_val;
    } // end ip

    result += el_val * vol; // Norm with volume of element
  } // end loop elems

  // set optimization volume if not already set
  if (opt_vol_ == -1)
  {
    opt_vol_ = volume;
    LOG_DBG2(ms) << "CMFD: calculated volume =" << opt_vol_;
  }
  // norm by the volume of the optimization domain
  result *= 1.0/opt_vol_;
  LOG_DBG(ms) << "CMFD: exit normed -> " << result;
  return result;
}

void MagSIMP::CalcMagFluxDensGradient(Excitation& excite, Function* f, TransferFunction* tf)
{

  // the context->GetExcitation() is now the last one as we solve and store all excitations first before calculating the gradients
  //Transform* trans = f != NULL && f->GetExcitation() != NULL ? f->GetExcitation()->transform : NULL; // even ->transform might be NULL

  // the gradient is < lambda^T, K' * A >
  assert(excite.sequence == f->ctxt->sequence);

  double factor = 1.0/ opt_vol_; // factor for norming the gradient; same as in objective function
  // calc lambda^T *  K' * A -> this already stores the results by AddGradient()!
  CalcU1KU2(tf, adjoint.Get(excite, f)->elem[App::MAG], App::MAG, forward.Get(excite)->elem[App::MAG], NULL, factor, STANDARD, f);
}

void MagSIMP::CalcN(LinearFormContext* form, Vector<double>& N)
{
  shared_ptr<EntityList> el = form->GetEntities();
  assert(el->GetRegion() != NO_REGION_ID);
  assert(el->GetRegion() >= 0);

  LOG_DBG(ms) << "CN: ft=" << el->GetType() << " r=" << el->GetRegion() << " fn="  << form->ToString();

  /*
  BaseBDBInt* bdb = dynamic_cast<BaseBDBInt*>(form);
  assert(bdb);

  StdVector<LocPoint> intPoints; // Get integration Points
  LocPointMapped lp;
  StdVector<double> weights;
  IntegOrder order;
  IntScheme::IntegMethod method = IntScheme::UNDEFINED;

  form->

  SolutionType solt = MAG_POTENTIAL;
  shared_ptr<BaseFeFunction> fe = form->GetPde()->GetFeFunction(solt);

  form->Get()->

  fe->

  // get the form first
  LinearForm* integrator = form->GetIntegrator();
  //BaseBDBInt* bdb = dynamic_cast<BaseBDBInt*>(context->pde->GetAssemble()->GetBiLinForm("CurlCurlIntegrator", REGION, context->pde)->GetIntegrator());
  assert(integrator != NULL);
  //LOG_DBG2(ms) << "CMFD enter: bdb=" << bdb->GetName();

  shared_ptr<EntityList> ent = form->GetEntities();
  EntityIterator iter = ent->GetIterator();

  // Calculate volume

  double vol_sum = 0.0;
  // loop over elements
  while (!iter.IsEnd())
  {
    // TODO shall be made faster
    const Elem* elem = iter.GetElem();
    shared_ptr<ElemShapeMap> esm = domain->GetGrid()->GetElemShapeMap(elem);
    double vol = esm->CalcVolume();
    // loop over integration points
      //N(elem) = vol* weights[ip] * lp.jacDet * Shapfct(=1?)

    vol_sum += vol;
    iter++;
  }
  return 0.0; //return vol_sum
  */
}

double MagSIMP::CalcMagCoupling(Excitation& excite, Function* f)
{
  if(GetMultipleExcitation()->excitations.GetSize() != 2)
    throw Exception("'magCoupling' requires two coils and enabled multiple_excitations");

  // cfs makes two coupling functions for two excitations with a weight of 0.5 each.
  // we return value 0 for the first excitation and 2*coupling for the second excitation
  assert(GetMultipleExcitation()->excitations.GetSize() == 2);
  if(excite.index == 0)
    return 0.0;

  assert(excite.index == 1);
  Excitation excite_A = GetMultipleExcitation()->excitations[0];
  Excitation excite_B = GetMultipleExcitation()->excitations[1];
  // Coupling can only be calculated for exactly two coils
  assert(excite_A.index == 0 && excite_B.index == 1);
  StdVector<LinearFormContext*>& forms_A = excite_A.forms;
  StdVector<LinearFormContext*>& forms_B = excite_B.forms;
  // get the two As
  // const Vector<double>& A = forward.Get(excite,NULL)->GetRealVector(StateSolution::RAW_VECTOR);
  // StdVector<SingleVector*>& A = forward.Get(excite)->elem[App::MAG]; // TODO is forward unknown here?

  assert((forms_A.GetSize() == 1) && (forms_B.GetSize() == 1));
  LinearFormContext* form_A = forms_A[0];
  LinearFormContext* form_B = forms_B[0];

  Vector<double> N1;
  CalcN(form_A, N1);

  /*Vector<double> N2 = zeros(size(A))
   vol_2 = CalcN(form_B, N2);
  Vector<double> N2 = CalcN(form_B);
  double N1AB = N1.inner(A_B);
  double N1AA = N1.inner(A_A);
  double N2AB = N2.inner(A_B);
  double k = (vol_2/vol_1)*((N1AB*N1AB)/(N1AA*N2AB));//
   */
  //LOG_DBG(ms) << "CCARH: ent A=" << ent_A->GetName();
  //LOG_DBG(ms) << "CCARH: ent B=" << ent_B->GetName();
  //assert(false);
  }

void MagSIMP::CalcCouplingGradient(Excitation& excite, Function* f, TransferFunction* tf)
{
  //TODO
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

  // Write the solution from the algebraic system to the vector real_nu, real_nu will be needed in SetElementK to get the correct nu
  SolutionType solt = MAG_POTENTIAL;
  shared_ptr<BaseFeFunction> fe = context->pde->GetFeFunction(solt);
  real_nu = dynamic_cast<Vector <double>& >(*(fe->GetSingleVector()));

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
    LOG_DBG2(ms) << "CMFAR e=" << e << " el=" << de->elem->elemNum<< "coef_S= " << coef_S->GetTensor().ToString(2) << " M=" << M.ToString(2);

    // now the part -2 * M * A
    Vector<double>* vec = dynamic_cast<Vector<double>*>(sol[de->GetElementSolutionIndex()]);
    LOG_DBG3(ms) << "CMFAR e=" << e << " esi=" << de->GetElementSolutionIndex() << " nodal values=" << vec->ToString(2);
    assert(vec != NULL);
    const Vector<double>& a = *vec; // a = the vector potential in the element

    rhs_el.Resize(a.GetSize());
    M.Mult(a, rhs_el);
    rhs_el *= -2;

    LOG_DBG3(ms) << "CMFAR e=" << e << " a=" << a.ToString(2) << " -> -2 M*a: r=" << rhs_el.ToString(2);

    c->GetIntegrator()->GetFeSpace1()->GetElemEqns(eqn, de->elem);
    LOG_DBG3(ms) << "CMFAR e=" << e << " eqn=" << eqn.ToString(2);
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

        out[eqn_idx] += rhs_el[n] * vol; // norm with volume of element
        LOG_DBG2(ms) << "CMFAR: n=" << n << " eqn_idx=" << eqn_idx << " normed_rhs_el[n]= " << vol * rhs_el[n] << " out[eqn_idx]=" << out[eqn_idx];
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

void MagSIMP::CalcCouplingAdjRHS(Excitation& excite, Function* f, Vector<double>& out)
{
  /* The right hand sides lood as follows:
   * case A:
   * (2/(Vol_A)^3)* <N_1, A_B>/(<N_1, A_A><N_2, A_B>) * N_1^T - (<N_1, A_B>)^2/(<N_1, A_A><N_2, A_B^2>) * N_2^T
   * case B:
   * (<N_1, A_B>^2)/(<N_1, A_A>^2 <N_2, A_B>) * N_1^T
   */
  assert(GetMultipleExcitation()->excitations.GetSize() == 2);
  assert(excite.index == 0 || excite.index == 1);

  // we get the region from the coil name in the forms
//  StdVector<LinearFormContext*>& forms = excite.forms;
//  assert(forms.GetSize() == 1);
//  LinearFormContext* form = forms[0];
//  shared_ptr<EntityList> ent = form->GetEntities();
//  LOG_DBG(ms) << "CCARH: ent=" << ent->GetName();

  assert(out.GetSize() == 0);
  // stupid dummy code :((
  Excitation exc_A = GetMultipleExcitation()->excitations[0];
  Excitation exc_B = GetMultipleExcitation()->excitations[1];
  const Vector<double>& A_A = forward.Get(exc_A,NULL)->GetRealVector(StateSolution::RAW_VECTOR);
  const Vector<double>& A_B = forward.Get(exc_B,NULL)->GetRealVector(StateSolution::RAW_VECTOR);

  StdVector<LinearFormContext*>& forms = exc_A.forms;
  assert(forms.GetSize() == 1);
  LinearFormContext* form_A = forms[0];
  forms = exc_B.forms;
  assert(forms.GetSize() == 1);
  LinearFormContext* form_B = forms[0];
  Vector<double> N1;
  Vector<double> N2;
  CalcN(form_A, N1);
  CalcN(form_B, N2);
  if (excite.index == 0) {
    //calc rhs...
  } else if (excite.index == 1) {
    //calc rhs...
  } else {
    throw Exception("There should only be two excitations");
  }

  out.Resize(A_A.GetSize(),0.0);
  out[0] = 1;
}

template <class T1, class T2>
void MagSIMP::SetElementK(Function* f, DesignElement* de, const TransferFunction* tf, App::Type app, DenseMatrix* mat_out, bool derivative, CalcMode calcMode, double ev)
{
  OptimizationMaterial* mat = f->ctxt->mat;
  Matrix<T1>& out = dynamic_cast<Matrix<T1>& >(*mat_out);

  //assert(app != App::MAG); // shall be in MagSIMP.cc

  // Not 100% sure about this but with it we get the correct nu
  // Define a new Vector buffer_store containing the information from StateSolution
  // Undo the overwriting from StateSolution.cc line 353 and 445 to get the right nu
  Vector <double> buffer_store;
  SolutionType solt = MAG_POTENTIAL;
  shared_ptr<BaseFeFunction> fe = context->pde->GetFeFunction(solt);
  buffer_store = dynamic_cast<Vector <double>& >(*(fe->GetSingleVector()));
  dynamic_cast<Vector<double>& >(*(fe->GetSingleVector())) = real_nu;

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
    dynamic_cast<Vector<double>& >(*(fe->GetSingleVector())) = buffer_store;
    //LOG_DBG3(ms) << "e=" << de->elem->elemNum << " K_0=" << stiffness.ToString(2);

    double nu_r = GetRelactivity(de->elem, domain->GetGrid()->GetDim());
    // simulation: BDB with D=(d 0; 0 d) with d = nu_0*nu_r
    // optimization: d = nu_0 * nu_r * f(rho) - nu_0 * f(rho) + nu_0
    // derivative: d = nu_0 * nu_r * f'(rho) - nu_0 * f'(rho)

    // in difference to elasticity where K' = my(rho)'*K_0
    // we have K' = alpha * K_0 where alpha = (nu_0*nu_r - nu_0) / (nu_0 * nu_r)

    double alpha = (nu_0*nu_r*d_rho - nu_0) / (nu_0 * nu_r);
    //LOG_DBG3(ms) << "e=" << de->elem->elemNum << " "
     //   "nu_0=" << nu_0 << " nu_r=" << nu_r << " rho=" << de->GetDesign(DesignElement::SMART) << " d_rho=" << d_rho << " alpha=" << alpha;

    Assign(out, stiffness, alpha); // out = alpha * stiffness

    //LOG_DBG3(ms) << "out=" << out.ToString(2);
    break;
  }

  default:
    SIMP::SetElementK(f, de, tf, app, mat_out, derivative, calcMode, ev);
    return; // other cases should be handled in SIMP
  }  // end switch
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
    return sel_xy_;
  default:
    assert(false);
    return sel_y_;
  }
}
