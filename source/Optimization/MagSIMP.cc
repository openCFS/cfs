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
#include "Forms/LinForms/BUInt.hh"
#include "Forms/LinForms/LinearForm.hh"
#include "Utils/tools.hh"
#include "Domain/Domain.hh"
#include "PDE/MagneticPDE.hh"
#include "PDE/MagEdgePDE.hh"

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
    if(!derivative)
      return CalcMagFluxDensity(excite, f);
    else
    {
      CalcMagFluxDensGradient(excite, f);
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

  assert(!f->elements.IsEmpty());
  assert(f->region != ALL_REGIONS);
  for(unsigned int e = 0; e < f->elements.GetSize(); e++)
  {
    DesignElement* de = f->elements[e];

    // element solution
    Vector<double>* vec = dynamic_cast<Vector<double>*>(sol[de->GetElementSolutionIndex()]);
    assert(vec != NULL);
    const Vector<double>& a = *vec; // a = the vector potential in the element
    LOG_DBG3(ms) << "CMFD e=" << e << " el=" << de->elem->elemNum << " esi=" << de->GetElementSolutionIndex() << " a=" << a.ToString(2);

    // prepare to get the curl operator
    el.SetElement(de->elem);

    BaseFE* ptFe = bdb->GetFeSpace1()->GetFe(el.GetIterator(), method, order );

    bdb->GetIntScheme()->GetIntPoints(Elem::GetShapeType(de->elem->type), method, order, intPoints, weights );
    LOG_DBG2(ms) << "CMFD i=" << e <<  " method=" << method << " order=" << order.ToString() << " iP=" << intPoints;
    assert(method != IntScheme::UNDEFINED);
    assert(!intPoints.IsEmpty());
    // Get shape map from grid
    shared_ptr<ElemShapeMap> esm = domain->GetGrid()->GetElemShapeMap(de->elem);

    double el_val = 0.0;

    for(unsigned int ip = 0; ip < intPoints.GetSize(); ip++)
    {
      // Calculate for each integration point the LocPointMapped
      lp.Set(intPoints[ip], esm, weights[ip]);

      // Call the CalcBMat()-method
      assert(bdb->GetBOp());
      // std::cout << "bdb->B " << bdb->GetBOp()->GetName() << "\n";
      bdb->GetBOp()->CalcOpMat(M, lp, ptFe);
      assert(M.GetNumCols() == a.GetSize());
      assert(M.GetNumRows() == domain->GetGrid()->GetDim());
      LOG_DBG3(ms) << "CMFD: e= " << e << " ip=" << ip << "/" << intPoints[ip].coord.ToString() << " w=" << weights[ip] << " jacDet=" << lp.jacDet << " M_" << ip << "=" << M.ToString(2);

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


      LOG_DBG3(ms) << "CMFD: e= " << e << " flux_dens=" << flux_dens.ToString(2) << " Sfd=" << S_flux_dens.ToString(2) << " el_val= " << el_val;
    } // end ip

    result += el_val;

  } // end loop elems
  if(dim == 3)
  {
    result /= intPoints.GetSize();
  }

  result /= f->elements.GetSize(); // norm the function
  LOG_DBG(ms) << "CMFD: exit -> " << result;
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

    DesignDependentRHS* ptr_rhs = dt == DesignElement::NONFERRITE_DENSITY ? &rhs : NULL;

    double factor = 1.0/ f->elements.GetSize(); // factor for norming the gradient; same as in objective function
    // calc lambda^T *  K' * A -> this already stores the results by AddGradient()!
    CalcU1KU2(tf, adjoint.Get(excite, f)->elem[App::MAG], App::MAG, forward.Get(excite)->elem[App::MAG], ptr_rhs, factor, STANDARD, f);
  }

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

  for(unsigned int e = 0; e < f->elements.GetSize(); e++)
  {
    DesignElement* de = f->elements[e];
    elemList.SetElement(de->elem);

    // M = B^T S B as in BDBInt::CalcElementMatrix, just S from above instead of D with material and density
    bdb->CalcElementMatrix(M, it, it);
    LOG_DBG3(ms) << "CMFAR e=" << e << " el=" << de->elem->elemNum<< "coef_S= " << coef_S->GetTensor().ToString(2) << " M=" << M.ToString(2);

    // now the part -2 * M * A
    Vector<double>* vec = dynamic_cast<Vector<double>*>(sol[de->GetElementSolutionIndex()]);
    LOG_DBG3(ms) << "CMFAR e=" << e << " esi=" << de->GetElementSolutionIndex() << " nodal values=" << vec->ToString(2);
    assert(vec != NULL);
    Vector<double>& a = *vec; // a = the vector potential in the element
    if(dim == 3)
    {
      a /= 2; // /2 because of edge element?
    }

    rhs_el.Resize(a.GetSize());
    M.Mult(a, rhs_el);
    rhs_el *= -2;

    LOG_DBG3(ms) << "CMFAR e=" << e << " a=" << a.ToString(2) << " -> -2 M*a: r=" << rhs_el.ToString(2);

    c->GetIntegrator()->GetFeSpace1()->GetElemEqns(eqn, de->elem);
    LOG_DBG3(ms) << "CMFAR e=" << e << " eqn=" << eqn.ToString(2);
    assert(rhs_el.GetSize() == eqn.GetSize());

    for(unsigned int n = 0; n < eqn.GetSize(); n++)
    {
      // the equation number is 1 based with 0 indicating HDBC and constrained nodes for negative indices. The equation index is 0-based!
      int eqn_nbr = eqn[n];
      if(eqn_nbr <= 0){
        LOG_DBG2(ms) << "CMFAR: n=" << n << " eqn_nbr=" << eqn_nbr << " -> skip RHS node";
      }
      else
      {
        unsigned int eqn_idx = eqn_nbr-1;

        out[eqn_idx] += rhs_el[n]; // we don't norm here, we moved the norming to CalcMagFluxDensGradient
        LOG_DBG2(ms) << "CMFAR: n=" << n << " eqn_idx=" << eqn_idx << " normed_rhs_el[n]= " << ((1.0/f->elements.GetSize()) * rhs_el[n]) << " out[eqn_idx]=" << out[eqn_idx] << " normed= " << 1.0/f->elements.GetSize();
      }
    }
  } // end loop elements
  delete bdb;
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
  //Vector <double> buffer_store;
  //SolutionType solt = MAG_POTENTIAL;
  //shared_ptr<BaseFeFunction> fe = context->pde->GetFeFunction(solt);
  //buffer_store = dynamic_cast<Vector <double>& >(*(fe->GetSingleVector()));
  //dynamic_cast<Vector<double>& >(*(fe->GetSingleVector())) = real_nu;

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
    //dynamic_cast<Vector<double>& >(*(fe->GetSingleVector())) = buffer_store;
    //LOG_DBG3(ms) << "e=" << de->elem->elemNum << " K_0=" << stiffness.ToString(2);

    double nu_r = GetRelactivity(de->elem, domain->GetGrid()->GetDim());
    // simulation: BDB with D=(d 0; 0 d) with d = nu_0*nu_r
    // optimization: d = nu_0 * nu_r * f(rho) - nu_0 * f(rho) + nu_0
    // derivative: d = nu_0 * nu_r * f'(rho) - nu_0 * f'(rho)

    // in difference to elasticity where K' = my(rho)'*K_0
    // we have K' = alpha * K_0 where alpha = (nu_0*nu_r - nu_0) / (nu_0 * nu_r)

    double alpha = (nu_0*nu_r*d_rho - nu_0) / (nu_0 * nu_r);
    LOG_DBG3(ms) << "e=" << de->elem->elemNum << " "
        "nu_0=" << nu_0 << " nu_r=" << nu_r << " rho=" << de->GetDesign(DesignElement::SMART) << " d_rho=" << d_rho << " alpha=" << alpha;

    Assign(out, stiffness, alpha); // out = alpha * stiffness

    LOG_DBG3(ms) << "out=" << out.ToString(2);
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

  LOG_DBG3(ms) << "SCUKU: de=" << de->ToString() << " tf=" << tf->ToString() << " d_rho= " << d_rho << " in_out=" << in_out.ToString(2) << " mat_vec= " << mat_vec->ToString();

}

const Matrix<double>& MagSIMP::GetSelectionMatrix(const Function* f) const
{
  switch(f->GetType())
  {
  case Function::SQR_MAG_FLUX_DENS_X:
    return sel_x_;
  case Function::SQR_MAG_FLUX_DENS_Y:
    return sel_y_;
  default:
    assert(false);
    return sel_y_;
  }
}



