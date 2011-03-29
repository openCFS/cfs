#include "Optimization/StressConstraint.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Design/DesignElement.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "Domain/domain.hh"
#include "Forms/baseForm.hh"
#include "Forms/linElastInt.hh"
#include "Driver/assemble.hh"
#include "MatVec/opdefs.hh"

using namespace std;

DECLARE_LOG(sc)
DEFINE_LOG(sc, "stressConstraint")


namespace CoupledField
{

template<typename T>
StressConstraint<T>::StressConstraint(Excitation* excite, Function* f, ErsatzMaterial* em, ErsatzMaterial::Solutions* forward) :
    elemList(domain->GetGrid())
{
  assert(f->GetType() == Function::STRESS);

  this->excite = excite;
  this->f = f;
  this->em =em;
  this->forward = forward;

  space = em->GetDesign();

  // global initializations
  M = dynamic_cast<MechPDE*>(em->ToPDE(Optimization::MECH))->GetVonMisesMatrix(domain->GetGrid()->GetDim());

  if(f->region != ALL_REGIONS && !space->Contains(f->region))
    tf = TransferFunction(Optimization::NO_APP, TransferFunction::FULL, 0.0, f->design);
  else
    tf = *(space->GetTransferFunction(DesignElement::DENSITY, Optimization::STRESS)); // for qp = rho^3/rho^2.8 use SIMP with 0.2


  LOG_DBG2(sc) << "SC: tf=" << tf.ToString() << " M=" << M.ToString();

}

template<typename T>
void StressConstraint<T>::CalcStresses(Vector<double>& out)
{
  // output of penalized von mises stresses? Only used in !adjoint_rhs and !grad_contrib
  int res_idx = space->GetSpecialResultIndex(DesignElement::DEFAULT, DesignElement::PENALIZED_STRESS, DesignElement::NONE, DesignElement::PLAIN, excite->label);

  CalcStresses(STRESS, res_idx, out);
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

  double factor = mode == STRESS ? 1.0 : 2.0;


  StdVector<pair<Optimization::Application, Optimization::Application> > apps = GetApplications();

  for(unsigned int a = 0; a < apps.GetSize(); a++)
  {
    pair<Optimization::Application, Optimization::Application>& app = apps[a];

    all_u1_elem = &(forward->Get(*excite)->elem[app.first == Optimization::MECH ? Optimization::MECH : Optimization::ELEC]);
    all_u2_elem = &(forward->Get(*excite)->elem[app.second == Optimization::MECH ? Optimization::MECH : Optimization::ELEC]); // for the adjoint rhs we need no app2 solution

    // It might be that stress sensitive region is not within the design domain itself
    for(unsigned int e = 0, en = f->elements.GetSize(); e < en; e++)
    {
      DesignElement* de = f->elements[e];

      Setup(de, app.first, app.second, mode);

      // normal von Mises stress element results < stress, M*stress>
      M_stress2 = M * stress2;

      out[e] += Real(stress1.Inner(M_stress2)) * factor;

      LOG_DBG2(sc) << "CS de=" << de->ToString() << " stress1" << stress1.ToString() << " M_stress2=" << M_stress2.ToString() << " -> " << stress1.Inner(M_stress2) << "(" << out[e] << ")";

      // output von mises stress? Note, that this is excitation specific! Only for STRESS
      if(res_idx != -1)
      {
        de->specialResult[res_idx] = out[e];
        // LOG_DBG3(em) << "CS:sr de=" << de->ToString() << " res_idx=" << res_idx << " v=" << out[e];
      }
      LOG_DBG2(sc) << "CS de=" << de->ToString() << " rho=" << de->GetDesign(DesignElement::SMART) << " sMs=" << out[e] << " trans=" <<  tf.Transform(de, DesignElement::SMART);
    }
  }
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
  // app2 determines the pde. The elastic case is simply MECH, MECH


  // evaluate in the piezo case 4 times. Better nice code!
  Vector<double> alpha;
  CalcGlobalizationFactor(alpha);

  Matrix<T> stress_transp(1, domain->GetGrid()->GetDim() == 2 ? 3 : 6);
  Matrix<T> rhs_transp;

  bool harmonic = em->IsHarmonic();
  unsigned int dim = domain->GetGrid()->GetDim();

  out.Resize(forward->Get(*excite)->GetVector(ErsatzMaterial::Solution::RAW_VECTOR)->GetSize(), T());

  StdVector<pair<Optimization::Application, Optimization::Application> > apps = GetApplications();

  for(unsigned int a = 0; a < apps.GetSize(); a++)
  {
    pair<Optimization::Application, Optimization::Application>& app = apps[a];

    shared_ptr<EqnMap> eqnMap = em->ToPDE(app.second)->GetEqnMap(); // mech our coupling

    all_u1_elem = &(forward->Get(*excite)->elem[app.second == Optimization::MECH ? Optimization::MECH : Optimization::ELEC]);
    all_u2_elem = all_u1_elem; // for the adjoint rhs we need no app2 solution

    double factor = (app.first != app.second ? -1.0 : 1.0) * (harmonic ? -1.0 : -2.0); // see piezo adjoint rhs

    // It might be that stress sensitive region is not within the design domain itself
    for(unsigned int e = 0, en = f->elements.GetSize(); e < en; e++)
    {
      DesignElement* de = f->elements[e];

      Setup(de, app.first, app.second, ADJOINT_RHS);

      assert(stress1.GetSize() == stress_transp.GetNumCols());
      // we have to transpose the stress (3*1 or 6*1) manually and to the conjugate stuff
      for(unsigned int si = 0; si < stress1.GetSize(); si++)
        stress_transp[0][si] = Conj<T>(stress1[si]); // nothing changes in real case

      // final!
      rhs_transp = stress_transp * M_E2_B2;

      // there is a factor from the globalization function, which is the gradient of the glob function(func_val)
      rhs_transp *= factor * alpha[e];

      // sum it up to the global rhs vector
      assert(de->elem->connect.GetSize() * dim == rhs_transp.GetNumCols());
      for(unsigned int n = 0; n < de->elem->connect.GetSize(); n++)
      {
        unsigned int node =  de->elem->connect[n];

        for(unsigned int dof = 1; dof <= dim; dof++)
        {
          // fuck 1-based in GetNodeEqn() !!
          int eqn_nr = std::abs(eqnMap->GetNodeEqn(node,dof)); // map periodic bc to the master equation
          assert(eqn_nr >= 0);
          int eqn_idx = eqn_nr -1; // fuck 1-based!!
          // don't set the homogeneous dirichlet boundary conditions :)
          if(eqn_idx >= 0)
          {
            out[eqn_idx] += rhs_transp[0][dim * n + (dof-1)];

            LOG_DBG3(sc) << "CAR de=" << de->ToString() << " alpha=" << alpha[e] << " factor=" << factor
                << " n=" << n << " node=" << node
                << " dof=" << dof << " eqn_idx=" << eqn_idx << " idx=" << (dim * n + (dof-1)) << " val="
                << rhs_transp[0][dim * n + (dof-1)] << " -> " << out[eqn_idx];
          }
        } // dof
      } // node
    } // elements
  } // apps
}

template<typename T>
void StressConstraint<T>::CalcGlobalizationFactor(Vector<double>& out)
{
  const Function::Local* local = f->GetLocal();
  assert(local != NULL);

  Vector<double> stress;
  CalcStresses(stress);
  out.Resize(stress.GetSize());

  // this is a complicated way to calc element wise power*max(0,stress)^(power-1) but that way
  // we are open for further globalization functions and it is not that expensive

  const StdVector<Function::Local::Identifier>& vem = local->virtual_elem_map;
  assert(vem.GetSize() == stress.GetSize());

  for(unsigned int i = 0; i < vem.GetSize(); i++)
  {
    const Function::Local::Identifier& id = vem[i];
    assert(id.neighbor.GetSize() == 0);
    double gfv = id.EvalFunction(local, true, stress[i]);
    //int idx = design->Find(id.element);
    // unsigned int idx = id.element->GetIndex();
    out[i] = gfv;
  }

  LOG_DBG2(sc) << "CGF ex=" << excite->index << " alpha=" << out.ToString();
}

template<typename T>
void StressConstraint<T>::Setup(DesignElement* de, Optimization::Application app1, Optimization::Application app2, Mode mode)
{
  BaseForm* form1 = em->GetForm(de->elem->regionId, app1, Optimization::NO_APP, true, true); // global!
  BaseForm* form2 = em->GetForm(de->elem->regionId, app2, Optimization::NO_APP, true, true); // global!

  BaseMaterial* bimat1 = space->GetBiMaterial(de->elem->regionId, app1, false);
  BaseMaterial* bimat2 = space->GetBiMaterial(de->elem->regionId, app2, false);

  static Vector<double> intPoint;
  static Matrix<double> E2_B2; // adjoint case only
  static Matrix<double> tmp;   // for transposing [e] to [e]^T


  // we need to be careful to use the right index!!
  Vector<T>& u1_elem = dynamic_cast<Vector<T>& >(*((*all_u1_elem)[de->GetElementSolutionIndex()]));
  Vector<T>& u2_elem = dynamic_cast<Vector<T>& >(*((*all_u2_elem)[de->GetElementSolutionIndex()]));

  // apply our own physical densities!
  form1->GetScaledMaterial(tf.Transform(de, DesignElement::SMART), false, bimat1, E1);
  if(app1 == Optimization::PIEZO_COUPLING)
  {
    tmp = E1;
    tmp.Transpose(E1);
  }


  if(mode == GRAD_STRESS)
    form2->GetScaledMaterial(tf.Derivative(de, DesignElement::SMART), false, bimat2, E2);
  else
    form2->GetScaledMaterial(tf.Transform(de, DesignElement::SMART), false, bimat2, E2);

  if(app2 == Optimization::PIEZO_COUPLING)
  {
    tmp = E2;
    tmp.Transpose(E2);
  }


  // B1 and B2
  de->elem->ptElem->GetCoordMidPoint(intPoint);
  form1->CalcBMatOnly(B1, intPoint, de->elem);
  form2->CalcBMatOnly(B2, intPoint, de->elem);

  // left side stress
  strain1 = B1 * u1_elem;
  stress1 = E1 * strain1;
  LOG_DBG3(sc) << "SE de=" << de->ToString() << " strain1=" << strain1.ToString() << " stress1=" << stress1.ToString();

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
    LOG_DBG3(sc) << "SE de=" << de->ToString() << " strain2=" << strain2.ToString() << " stress2=" << stress2.ToString();
  }
}

template<typename T>
StdVector<pair<Optimization::Application, Optimization::Application> >  StressConstraint<T>::GetApplications()
{
  StdVector<pair<Optimization::Application, Optimization::Application> > result;

  if(f->GetStressType() == Function::MECH)
   result.Push_back(std::make_pair(Optimization::MECH, Optimization::MECH));
  else
  {
    // one of three piezo case - is the stress constraint defined on a piezo region ?!
    RegionIdType reg = f->elements[0]->elem->regionId;
    // global as we might have no PiezoSIMP (e.g. harvester)
    BaseForm* form = em->GetForm(reg, Optimization::MECH, Optimization::ELEC, false, true);
    if(form == NULL)
      throw Exception("piezoelectric stress constraint not defined on a piezoelectric region");

    switch(f->GetStressType())
    {
    case Function::ONLY_MECH_PIEZO: // special case only
      result.Push_back(std::make_pair(Optimization::MECH, Optimization::MECH));
      break;

    case Function::ONLY_PIEZO_PIEZO: // special case only
      result.Push_back(std::make_pair(Optimization::PIEZO_COUPLING, Optimization::PIEZO_COUPLING));
      break;

    case Function::PIEZO: // standard piezo case
      result.Push_back(std::make_pair(Optimization::MECH, Optimization::MECH));
      result.Push_back(std::make_pair(Optimization::PIEZO_COUPLING, Optimization::MECH));
      result.Push_back(std::make_pair(Optimization::PIEZO_COUPLING, Optimization::PIEZO_COUPLING));
      result.Push_back(std::make_pair(Optimization::MECH, Optimization::PIEZO_COUPLING));
      break;

    default:
      assert(false);
    }
  }

  return result;

}


// Explicit template instantiation
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
template class StressConstraint<double> ;
template class StressConstraint<complex<double> > ;
#endif



}

