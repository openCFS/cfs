#include <assert.h>
#include <cstdlib>
#include <map>
#include <ostream>
#include <string>

#include "DataInOut/Logging/LogConfigurator.hh"
#include "Domain/CoefFunction/CoefFunctionOpt.hh"
#include "Domain/Domain.hh"
#include "Domain/ElemMapping/Elem.hh"
#include "Domain/ElemMapping/ElementAccess.hh"
#include "Domain/Mesh/Grid.hh"
#include "Driver/Assemble.hh"
#include "Forms/BiLinForms/BDBInt.hh"
#include "FeBasis/BaseFE.hh"
#include "General/defs.hh"
#include "General/Environment.hh"
#include "General/Exception.hh"
#include "MatVec/SingleVector.hh"
#include "MatVec/opdefs.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Excitation.hh"
#include "Optimization/Function.hh"
#include "Optimization/StressConstraint.hh"
#include "PDE/SinglePDE.hh"
#include "PDE/MechPDE.hh"
#include "Utils/tools.hh"

namespace CoupledField {
class BaseMaterial;
}  // namespace CoupledField

using namespace std;

DEFINE_LOG(sc, "stressConstraint")


namespace CoupledField
{

template<typename T>
StressConstraint<T>::StressConstraint(Excitation* excite, Function* f, ErsatzMaterial* em, StateContainer* forward) :
    elemList(domain->GetGrid())
{
  Function::Type type = f->GetType();
  assert(type == Function::GLOBAL_STRESS || type == Function::LOCAL_STRESS
      || type == Function::LOCAL_BUCKLING_LOAD_FACTOR || type == Function::GLOBAL_BUCKLING_LOAD_FACTOR);

  this->excite = excite;
  this->f = f;
  this->em = em;
  this->forward = forward;

  space = em->GetDesign();

  this->form = NULL;
  this->u1_elem_ptr = NULL;
  this->u2_elem_ptr = NULL;

  // global initializations
  if(type == Function::LOCAL_BUCKLING_LOAD_FACTOR || type == Function::GLOBAL_BUCKLING_LOAD_FACTOR)
    // for the local buckling load factor we need the norm of the stress
    M = dynamic_cast<MechPDE*>(em->context->ToPDE(App::MECH, true))->GetHillMandelMatrix(domain->GetDim());
  else
    M = dynamic_cast<MechPDE*>(em->context->ToPDE(App::MECH, true))->GetVonMisesMatrix(domain->GetDim());

  if(f->region != ALL_REGIONS && !space->Contains(f->region))
    tf = TransferFunction(App::NO_APP, TransferFunction::FULL, 0.0, f->GetDesignType());
  else
    tf = *(space->GetTransferFunction(DesignElement::DENSITY, App::STRESS, true)); // for qp = rho^3/rho^2.8 use SIMP with 0.2

  LOG_DBG2(sc) << "SC: tf=" << tf.ToString() << " M=" << M.ToString();

}

template<typename T>
Vector<double> StressConstraint<T>::CalcStresses()
{
  // output of penalized von mises stresses? Only used in !adjoint_rhs and !grad_contrib
  int res_idx = space->GetSpecialResultIndex(DesignElement::DEFAULT, DesignElement::QUADRATIC_VM_STRESS, DesignElement::NONE, DesignElement::PLAIN, excite->label);

  Vector<double> out;

  CalcStresses(STRESS, res_idx, out);

  return out;
}

template<typename T>
void StressConstraint<T>::CalcGradStresses(Vector<double>& out)
{
  CalcStresses(GRAD_STRESS, -1, out);
}


template<typename T>
void StressConstraint<T>::CalcStresses(Mode mode, int res_idx, Vector<double>& out)
{
  assert(mode == STRESS || mode == GRAD_STRESS);

  out.Resize(f->elements.GetSize(), 0.0);

  Vector<T> M_stress2; // temporary

  double fm = mode == STRESS ? 1.0 : 2.0;

  StdVector<pair<App::Type, App::Type> > apps = GetApplications();

  for(unsigned int a = 0; a < apps.GetSize(); a++)
  {
    pair<App::Type, App::Type>& app = apps[a];

    all_u1_elem = &(forward->Get(excite)->elem[app.first == App::MECH ? App::MECH : App::ELEC]);
    all_u2_elem = &(forward->Get(excite)->elem[app.second == App::MECH ? App::MECH : App::ELEC]); // for the adjoint rhs we need no app2 solution

    ElementAccess ea(em->context->GetBiLinFormContext(f->elements[0]->elem->regionId, app.first, app.second, true));

    // It might be that stress sensitive region is not within the design domain itself
    for(unsigned int e = 0, en = f->elements.GetSize(); e < en; e++)
    {
      DesignElement* de = f->elements[e];
      ea.SetElem(de->elem);
      if (mode == STRESS)
        SetupElement(&ea, de, app.first, mode);
      else {
        assert(mode == GRAD_STRESS);
        SetupElement(&ea, de, app.first, mode, de->GetType());
      }
      double elem_vol = ea.esm->CalcVolume();

      // we integrate over the element by averages summation and then multiplying with the volume
      for(unsigned int ip = 0; ip < ea.intPoints.GetSize(); ip++)
      {
        EvalIP(mode, &ea, ip);

        // normal von Mises stress element results < stress, M*stress>
        M_stress2 = M * stress2;

        T inner = stress1.Inner(M_stress2);

        // do the (de)normalization stuff outside of the inner product
        // ea.lpm.jacDet is elem vol and ref vol, by canceling elem_vol we have a field property
        double factor = fm * ea.weights[ip] * ea.lpm.jacDet / elem_vol;

        out[e] += factor * Real(inner);

        LOG_DBG2(sc) << "CS: de=" << de->ToString() << " ip=" << ip << " inner=" << inner << " w=" << ea.weights[ip] << " f=" << factor
                     << " stress1=" << stress1.ToString() << " M_stress2=" << M_stress2.ToString() << " -> " << (factor * Real(inner));
      }

      // output von mises stress? Note, that this is excitation specific! Only for STRESS
      if(res_idx != -1)
      {
        de->specialResult[res_idx] = out[e];
        // LOG_DBG3(em) << "CS:sr de=" << de->ToString() << " res_idx=" << res_idx << " v=" << out[e];
      }

      LOG_DBG2(sc) << "CS: de=" << de->ToString() << " rho=" << de->GetDesign(DesignElement::SMART) << " evol=" << ea.esm->CalcVolume()
          << " mode=" << mode << " sMs=" << out[e] << " trans=" <<  tf.Transform(de, DesignElement::SMART);
    }
  }
}


template<typename T>
double StressConstraint<T>::CalcElementStress(Mode mode, int res_idx, DesignElement* de)
{
  assert(mode == STRESS || mode == GRAD_STRESS);

  Vector<T> M_stress2; // temporary

  double res = 0.0;

  double fm = mode == STRESS ? 1.0 : 2.0;

  StdVector<pair<App::Type, App::Type> > apps = GetApplications();

  for(unsigned int a = 0; a < apps.GetSize(); a++)
  {
    pair<App::Type, App::Type>& app = apps[a];

    all_u1_elem = &(forward->Get(excite)->elem[app.first == App::MECH ? App::MECH : App::ELEC]);
    all_u2_elem = &(forward->Get(excite)->elem[app.second == App::MECH ? App::MECH : App::ELEC]); // for the adjoint rhs we need no app2 solution


    ElementAccess ea(em->context->GetBiLinFormContext(f->elements[0]->elem->regionId, app.first, app.second, true));

    ea.SetElem(de->elem);
    if (mode == STRESS)
      SetupElement(&ea, de, app.first, mode);
    else {
      assert(mode == GRAD_STRESS);
      SetupElement(&ea, de, app.first, mode, de->GetType());
    }
    double elem_vol = ea.esm->CalcVolume(); // see ::CalcStresses()

    // we integrate over the element by averages summation and then multiplying with the volume
    for(unsigned int ip = 0; ip < ea.intPoints.GetSize(); ip++)
    {
      EvalIP(mode, &ea, ip);

      // normal von Mises stress element results < stress, M*stress>
      M_stress2 = M * stress2;

      T inner = stress1.Inner(M_stress2);

      // do the (de)normalization stuff outside of the inner product
      double factor = fm * ea.weights[ip] * ea.lpm.jacDet / elem_vol;

      res += factor * Real(inner);

      LOG_DBG2(sc) << "CS: de=" << de->ToString() << " ip=" << ip << " inner=" << inner << " w=" << ea.weights[ip] << " f=" << factor
          << " stress1=" << stress1.ToString() << " M_stress2=" << M_stress2.ToString() << " -> " << (factor * Real(inner));
    }

    // output von mises stress? Note, that this is excitation specific! Only for STRESS
    if(res_idx != -1)
    {
      de->specialResult[res_idx] = res;
      // LOG_DBG3(em) << "CS:sr de=" << de->ToString() << " res_idx=" << res_idx << " v=" << out[e];
    }
    LOG_DBG2(sc) << "CES: de=" << de->ToString() << " rho=" << de->GetDesign(DesignElement::SMART) << " evol=" << ea.esm->CalcVolume() << " sMs=" << res << " trans=" <<  tf.Transform(de, DesignElement::SMART);
  }

  return res;
}


template<typename T>
void StressConstraint<T>::CalcAdjointRHS(Vector<T>& out)
{
  // elastic adjoint rhs w.r.p to one element static : - 2* (rho^p*E_0*B*u)^T * M * rho^p E_0 B
  // elastic adjoint rhs w.r.p to one element dynamic: - 1* (rho^p*E_0*B*u^*)^T * M * rho^p E_0 B
  // for the globalization it is multiplied by alpha
  // In the piezoelectric case it is (dynamic)
  // -1*alpha*(E1*B1*u1^*)^T*M*E2*B2 with app1=mech and app2=mech
  // +1*alpha*(E1*B1*u1^*)^T*M*E2*B2 with app1=piezo and app2=mech
  // -1*alpha*(E1*B1*u1^*)^T*M*E2*B2 with app1=piezo and app2=piezo
  // +1*alpha*(E1*B1*u1^*)^T*M*E2*B2 with app1=mech and app2=piezo
  // app2 determines the pde. The elastic case is simply App::MECH, App::MECH

  const Vector<double> stress = CalcStresses();

  Vector<double> local_values = stress;
  if(f->GetType() == Function::GLOBAL_BUCKLING_LOAD_FACTOR)
  {
    // it might be that stress sensitive region is not within the design domain itself
    for(unsigned int e = 0, en = f->elements.GetSize(); e < en; e++)
    {
      double vol = f->elements[e]->GetDesign(DesignElement::SMART);
      double ev_interpolated = em->GetMicroLoadFactor(vol);
      local_values[e] = ev_interpolated / std::sqrt(local_values[e]);
    }
  }

  Vector<double> alpha = CalcGlobalizationFactor(local_values, true);

  SingleVector* sv = em->forward.Get(excite)->GetVector(StateSolution::RAW_VECTOR);
  assert(out.GetSize() == 0); // if not, check the logic
  assert(sv->GetSize() > 0);
  out.Resize(sv->GetSize());
  out.Init(0);

  // It might be that stress sensitive region is not within the design domain itself
  for(unsigned int e = 0, en = f->elements.GetSize(); e < en; e++)
  {
    DesignElement* de = f->elements[e];

    if(f->GetType() == Function::GLOBAL_BUCKLING_LOAD_FACTOR)
    {
      double vol = de->GetDesign(DesignElement::SMART);
      double ev_interpolated = em->GetMicroLoadFactor(vol);
      alpha[e] *= - ev_interpolated * 0.5 / std::sqrt(stress[e]) / stress[e];
      LOG_DBG(sc) << "CARHS: de=" << de->ToString() << " ev=" << ev_interpolated << " stress=" << stress[e] << " -> " << alpha[e];
    }

    CalcElemAdjointRHS(de, alpha[e], out);
  }
}


template<typename T>
void StressConstraint<T>::CalcElemAdjointRHS(DesignElement* de, double alpha, Vector<T>& out_set)
{
  bool harmonic = em->context->IsComplex();

  Matrix<T> stress_transp(1, domain->GetDim() == 2 ? 3 : 6);
  Matrix<T> rhs_transp;

  // any bilinear shall do, hence take any region
  RegionIdType reg = f->region != ALL_REGIONS ? f->region : de->elem->regionId;
  ElementAccess ea(em->context->pde->GetAssemble()->GetBiLinForm(em->context->mat->stiff.integrator, reg));

  // out needs to be set already!
  assert(out_set.GetSize() > 0);

  StdVector<pair<App::Type, App::Type> > apps = GetApplications();

  for(unsigned int a = 0; a < apps.GetSize(); a++)
  {
    pair<App::Type, App::Type>& app = apps[a];

    all_u1_elem = &(forward->Get(*excite)->elem[app.first == App::MECH ? App::MECH : App::ELEC]);
    all_u2_elem = all_u1_elem; // for the adjoint rhs we need no app2 solution

    ea.SetElem(de->elem);
    LOG_DBG3(sc) << "CEAR: " << ea.ToString(2);
    SetupElement(&ea, de, app.first, ADJOINT_RHS);
    double elem_vol = ea.esm->CalcVolume(); // see ::CalcStresses()

    for(unsigned int ip = 0; ip < ea.intPoints.GetSize(); ip++)
    {
      // Calculate for each integration point the LocPointMapped
      EvalIP(ADJOINT_RHS, &ea, ip);

      assert(stress1.GetSize() == stress_transp.GetNumCols());
      // we have to transpose the stress (3*1 or 6*1) manually and do the conjugate stuff
      for(unsigned int si = 0; si < stress1.GetSize(); si++)
        stress_transp[0][si] = conj(stress1[si]); // nothing changes in real case

      // finally :)
      rhs_transp = stress_transp * M_E2_B2;
      LOG_DBG3(sc) << "CEAR: stress=" << stress_transp.ToString() << " MEB=" << M_E2_B2.ToString() << " rhs=" << rhs_transp;

      // include all the factors
      double factor = (harmonic ? -1.0 : -2.0) * ea.weights[ip] * ea.lpm.jacDet / elem_vol;

      LOG_DBG3(sc) << "CEAR: ip=" << ip << " w=" << ea.weights[ip] << " jD=" << ea.lpm.jacDet << " vol=" << ea.esm->CalcVolume() << " -> " << factor;

      // there is a factor from the globalization function, which is the gradient of the glob function(func_val)
      rhs_transp *= factor * alpha;

      LOG_DBG3(sc) << "CEAR: " << ea.ToString(3) << " alpha=" << alpha << " factor=" << factor;
      LOG_DBG3(sc) << "CEAR: idx=" << ea.elem_eqn_idx.ToString() << " rhs=" << rhs_transp.ToString();

      for(unsigned int n = 0; n < ea.elem_eqn_idx.GetSize(); n++)
      {
        // the equation number is 1 based with 0 indicating HDBC and constrained nodes for negative indices. The equation index is 0-based!
        if(ea.elem_eqn_idx[n] >= 0)
        {
          out_set[ea.elem_eqn_idx[n]] += rhs_transp[0][n];
          // LOG_DBG3(sc) << "CEAR: de=" << de->ToString() << " alpha=" << alpha << " factor=" << factor << " n=" << n << " val=" << rhs_transp[0][n] << " -> " << out_set[ea.elem_eqn_idx[n]];
        }
      }
    } // ip
  } // apps
}

template<typename T>
Vector<double> StressConstraint<T>::CalcGlobalizationFactor(const Vector<double>& local_values, bool gradient)
{
  Function::Local* local = f->GetLocal();
  assert(local != NULL);

  Vector<double> out(local_values.GetSize());

  // this is a complicated way to calc element wise power*max(0,stress)^(power-1) but that way
  // we are open for further globalization functions and it is not that expensive

  StdVector<Function::Local::Identifier>& vem = local->virtual_elem_map;
  assert(vem.GetSize() == local_values.GetSize());

  double org_value = f->GetValue(); // we use it for the local evaluation but don't want to spoil info.xml and plot.dat output

  for(unsigned int i = 0; i < vem.GetSize(); i++)
  {
    Function::Local::Identifier& id = vem[i];
    assert(id.neighbor.GetSize() == 0);
    f->SetValue(local_values[i]); // in case of objective LocalCondition::GetValue() is not available
    double gfv = id.EvalFunction(local, gradient);
    out[i] = gfv;
  }
  f->SetValue(org_value);

  LOG_DBG(sc) << "CGF: ex=" << excite->index << " alpha=" << out.ToString();

  return out;
}

template<typename T>
void StressConstraint<T>::SetupElement(ElementAccess* ea, DesignElement* de, App::Type app, Mode mode, DesignElement::Type direction)
{
  assert(app != App::PIEZO_COUPLING); // code deleted. See in pre-FE-Space
  OptimizationMaterial* mat = em->context->mat;
  BiLinFormContext* blfc = em->context->GetBiLinFormContext(de->elem->regionId, app, App::NO_APP, true);
  form = dynamic_cast<BaseBDBInt*>(blfc->GetIntegrator());

  // we need to be careful to use the right index!!
  u1_elem_ptr = dynamic_cast<Vector<T>* >((*all_u1_elem)[de->GetElementSolutionIndex()]);
  u2_elem_ptr = dynamic_cast<Vector<T>* >((*all_u2_elem)[de->GetElementSolutionIndex()]);

  Vector<T>& u1_elem = *u1_elem_ptr;
  Vector<T>& u2_elem = *u2_elem_ptr;

  LOG_DBG3(sc) << "SE: de=" << de->elem->elemNum << " a=" << app << " m=" << mode << " u1=" << u1_elem.ToString() << " u2=" << u2_elem.ToString();

  // set the element matrices
  ea->SetIP(0); // set an arbitrary integration point
  if (em->GetMethod() == ErsatzMaterial::SIMP_METHOD) {
    mat->GetOrgMatCoef(form)->GetTensor(E1, ea->lpm);
    E2 = E1; // TODO on the coupling case this comes from form2
    E1 *= tf.Transform(de, DesignElement::SMART);

    if(mode == GRAD_STRESS)
      E2 *= tf.Derivative(de, DesignElement::SMART);
    else
      E2 *= tf.Transform(de, DesignElement::SMART);
  } else {
    CoefFunctionOpt* coef = mat->GetMatCoef(blfc);

    assert(coef->GetMaterialDerivative() == DesignElement::NO_DERIVATIVE && direction != DesignElement::NO_MULTIMATERIAL);

    coef->GetTensor(E1, ea->lpm);

    E2 = E1;

    if(mode == GRAD_STRESS) {
      coef->SetToMaterialDerivative(direction);
      coef->GetTensor(E2, ea->lpm);
      coef->SetToOptimization();
    }
  }

  LOG_DBG3(sc) << "SE: de=" << de->elem->elemNum << " E1=" << E1.ToString() << " E2=" << E2.ToString();
}


template<typename T>
void StressConstraint<T>::EvalIP(Mode mode, ElementAccess* ea, unsigned int ip)
{
  assert(ea != NULL);

  // sets lpm
  ea->SetIP(ip);

  Vector<T>& u1_elem = *u1_elem_ptr;
  Vector<T>& u2_elem = *u2_elem_ptr;

  // B1 and B2
  form->GetBOp()->CalcOpMat(B1, ea->lpm, ea->CurrBaseFE());
  B2 = B1; // clearly this is from form2 in the coupling case
  LOG_DBG3(sc) << "SE: de=" << ea->CurrElem()->elemNum << " B1=" << B1.ToString();

  // left side stress
  strain1 = B1 * u1_elem;
  stress1 = E1 * strain1;
  LOG_DBG3(sc) << "SE: de=" << ea->CurrElem()->elemNum << " strain1=" << strain1.ToString() << " stress1=" << stress1.ToString();

  // right side stress - ignored for adjoint but clean that way
  if(mode == ADJOINT_RHS)
  {
    // right side
    E2_B2 = E2 * B2;
    M_E2_B2 = M * E2_B2;
  }
  else
  {
    strain2 = B2 * u2_elem;
    stress2 = E2 * strain2; // E2 might be dE2
    LOG_DBG3(sc) << "SE: de=" << ea->CurrElem()->elemNum << " strain2=" << strain2.ToString() << " stress2=" << stress2.ToString();
  }
}

template<typename T>
StdVector<pair<App::Type, App::Type> >  StressConstraint<T>::GetApplications()
{
  StdVector<pair<App::Type, App::Type> > result;

  if(f->GetStressType() == Function::MECH)
   result.Push_back(std::make_pair(App::MECH, App::MECH));
  else
  {
    assert(false);
    /* FIXME
    // one of three piezo case - is the stress constraint defined on a piezo region ?!
    RegionIdType reg = f->elements[0]->elem->regionId;

    BaseForm* form = em->GetForm(reg, App::MECH, App::ELEC, false);
    if(form == NULL)
      throw Exception("piezoelectric stress constraint not defined on a piezoelectric region");
*/
    switch(f->GetStressType())
    {
    case Function::ONLY_COUPLING: // special case only
      result.Push_back(std::make_pair(App::PIEZO_COUPLING, App::PIEZO_COUPLING));
      break;

    case Function::PIEZO: // standard piezo case
      result.Push_back(std::make_pair(App::MECH, App::MECH));
      result.Push_back(std::make_pair(App::PIEZO_COUPLING, App::MECH));
      result.Push_back(std::make_pair(App::PIEZO_COUPLING, App::PIEZO_COUPLING));
      result.Push_back(std::make_pair(App::MECH, App::PIEZO_COUPLING));
      break;

    default:
      assert(false);
      break;
    }
  }

  return result;

}

template<typename T>
typename StressConstraint<T>::DKuCache& StressConstraint<T>::GetdKuCache() {
  return dKuCache;
}


// Explicit template instantiation
template class StressConstraint<double> ;
template class StressConstraint<complex<double> > ;
template<class T> typename StressConstraint<T>::DKuCache StressConstraint<T>::dKuCache;


}

