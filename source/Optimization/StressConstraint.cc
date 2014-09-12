#include <assert.h>
#include <stdlib.h>
#include <map>
#include <ostream>
#include <string>

#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/Logging/log.hpp"
#include "Domain/Domain.hh"
#include "Domain/ElemMapping/Elem.hh"
#include "Domain/Mesh/Grid.hh"
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

DECLARE_LOG(sc)
DEFINE_LOG(sc, "stressConstraint")


namespace CoupledField
{

template<typename T>
StressConstraint<T>::StressConstraint(Excitation* excite, Function* f, ErsatzMaterial* em, ErsatzMaterial::Solutions* forward) :
    elemList(domain->GetGrid())
{
  assert(f->GetType() == Function::STRESS || f->GetType() == Function::STRESS_DENSITY);

  this->excite = excite;
  this->f = f;
  this->em =em;
  this->forward = forward;

  space = em->GetDesign();

  this->form1 = NULL;
  this->form2 = NULL;
  this->u1_elem_ptr = NULL;
  this->u2_elem_ptr = NULL;


  // global initializations
  assert(true);
  // FIXME M = dynamic_cast<MechPDE*>(em->ToPDE(Optimization::MECH))->GetVonMisesMatrix(domain->GetGrid()->GetDim());

  if(f->region != ALL_REGIONS && !space->Contains(f->region))
    tf = TransferFunction(Optimization::NO_APP, TransferFunction::FULL, 0.0, f->GetDesignType());
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

  double fm = mode == STRESS ? 1.0 : 2.0;

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

      // the element volume is actually only required if we no NOT want the stress density
      double elem_vol = domain->GetGrid()->GetElemShapeMap(de->elem, false)->CalcVolume();

      assert(false);
      Vector<double> weights;
      // FIXME const Vector<double>& weights = de->elem->ptElem->GetIntWeights();

      SetupElement(de, app.first, app.second, mode);

      // we integrate over the element by averages summation and then multiplying with the volume
      for(unsigned int ip = 1, ipn = 4/* FIXME de->elem->ptElem->GetNumIntPoints()*/; ip <= ipn; ip++) // 1-based!! :(
      {
        SetupIntPoint(de->elem, ip, mode);

        // normal von Mises stress element results < stress, M*stress>
        M_stress2 = M * stress2;

        T inner = stress1.Inner(M_stress2);

        // do the (de)normalization stuff outside of the inner product
        assert(false);
        double jac_det = 0.0; // FIXME de->elem->ptElem->CalcJacobianDetAtIp(ip, coords, de->elem);
        double factor = fm * weights[ip-1] * jac_det / (f->GetType() == Function::STRESS ? elem_vol : 1.0); // fuck 1-based!

        out[e] += factor * Real(inner);

        LOG_DBG2(sc) << "CS de=" << de->ToString() << " ip=" << ip << " inner=" << inner << " w=" << weights[ip-1] << " f=" << factor
                     << " stress1=" << stress1.ToString() << " M_stress2=" << M_stress2.ToString() << " -> " << (factor * Real(inner));
      }

      // output von mises stress? Note, that this is excitation specific! Only for STRESS
      if(res_idx != -1)
      {
        de->specialResult[res_idx] = out[e];
        // LOG_DBG3(em) << "CS:sr de=" << de->ToString() << " res_idx=" << res_idx << " v=" << out[e];
      }

      LOG_DBG2(sc) << "CS de=" << de->ToString() << " rho=" << de->GetDesign(DesignElement::SMART) << " ev=" << elem_vol << " sMs=" << out[e] << " trans=" <<  tf.Transform(de, DesignElement::SMART);
    }
  }

}

template<typename T>
void StressConstraint<T>::CalcAdjointRHS(Vector<T>& out)
{
  assert(false);
  /* FIXME

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
  Matrix<double> coords; // for the Jacobi determinant we need the coordinates

  bool harmonic = em->IsHarmonic();

  out.Resize(forward->Get(*excite)->GetVector(ErsatzMaterial::Solution::RAW_VECTOR)->GetSize(), T());

  StdVector<pair<Optimization::Application, Optimization::Application> > apps = GetApplications();

  for(unsigned int a = 0; a < apps.GetSize(); a++)
  {
    pair<Optimization::Application, Optimization::Application>& app = apps[a];

    // the second app determines u or phi.
    unsigned int dof = app.second ==  Optimization::MECH ? domain->GetGrid()->GetDim() : 1;
    shared_ptr<EqnMap> eqnMap = em->ToPDE(app.second == Optimization::PIEZO_COUPLING ? Optimization::ELEC : Optimization::MECH)->GetEqnMap();

    all_u1_elem = &(forward->Get(*excite)->elem[app.first == Optimization::MECH ? Optimization::MECH : Optimization::ELEC]);
    all_u2_elem = all_u1_elem; // for the adjoint rhs we need no app2 solution

    // It might be that stress sensitive region is not within the design domain itself
    for(unsigned int e = 0, en = f->elements.GetSize(); e < en; e++)
    {
      DesignElement* de = f->elements[e];

      // the element volume is actually only required if we no NOT want the stress density
      domain->GetGrid()->GetElemNodesCoord(coords, de->elem->connect, false); // no updated coordinates
      double elem_vol = de->elem->ptElem->CalcVolume(coords, false); // by default no axis symmetry!

      const Vector<double>& weights = de->elem->ptElem->GetIntWeights();

      SetupElement(de, app.first, app.second, ADJOINT_RHS);

      for(unsigned int ip = 1, ipn = de->elem->ptElem->GetNumIntPoints(); ip <= ipn; ip++) // 1-based!! :(
      {
        SetupIntPoint(de->elem, ip, ADJOINT_RHS);

        assert(stress1.GetSize() == stress_transp.GetNumCols());
        // we have to transpose the stress (3*1 or 6*1) manually and do the conjugate stuff
        for(unsigned int si = 0; si < stress1.GetSize(); si++)
          stress_transp[0][si] = Conj<T>(stress1[si]); // nothing changes in real case

        // finally :)
        rhs_transp = stress_transp * M_E2_B2;

        // include all the factors
        double jac_det = de->elem->ptElem->CalcJacobianDetAtIp(ip, coords, de->elem);
        double factor = (harmonic ? -1.0 : -2.0) * weights[ip-1] * jac_det  / (f->GetType() == Function::STRESS ? elem_vol : 1.0);

        // there is a factor from the globalization function, which is the gradient of the glob function(func_val)
        rhs_transp *= factor * alpha[e];

        // sum it up to the global rhs vector
        assert(de->elem->connect.GetSize() * dof == rhs_transp.GetNumCols());
        for(unsigned int n = 0; n < de->elem->connect.GetSize(); n++)
        {
          unsigned int node =  de->elem->connect[n];

          for(unsigned int d = 1; d <= dof; d++)
          {
            // fuck 1-based in GetNodeEqn() !!
            int eqn_nr = std::abs(eqnMap->GetNodeEqn(node,d)); // map periodic bc to the master equation
            assert(eqn_nr >= 0);
            int eqn_idx = eqn_nr -1; // fuck 1-based!!
            // don't set the homogeneous dirichlet boundary conditions :)
            if(eqn_idx >= 0)
            {
              out[eqn_idx] += rhs_transp[0][dof * n + (d-1)];

              LOG_DBG3(sc) << "CAR de=" << de->ToString() << " alpha=" << alpha[e] << " factor=" << factor
                  << " n=" << n << " node=" << node
                  << " d=" << d << " eqn_idx=" << eqn_idx << " idx=" << (dof * n + (d-1)) << " val="
                  << rhs_transp[0][dof * n + (d-1)] << " -> " << out[eqn_idx];
            }
          } // dof
        } // node
      } // ip
    } // elements
  } // apps */
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
void StressConstraint<T>::SetupElement(DesignElement* de, Optimization::Application app1, Optimization::Application app2, Mode mode)
{
  assert(false);
  /* FIXME
  form1 = em->GetForm(de->elem->regionId, app1, Optimization::NO_APP, true);
  form2 = em->GetForm(de->elem->regionId, app2, Optimization::NO_APP, true);

  BaseMaterial* bimat1 = space->GetBiMaterial(de->elem->regionId, app1, false);
  BaseMaterial* bimat2 = space->GetBiMaterial(de->elem->regionId, app2, false);

  static Vector<double> intPoint;
  static Matrix<double> tmp;   // for transposing [e] to [e]^T


  // we need to be careful to use the right index!!
  u1_elem_ptr = dynamic_cast<Vector<T>* >((*all_u1_elem)[de->GetElementSolutionIndex()]);
  u2_elem_ptr = dynamic_cast<Vector<T>* >((*all_u2_elem)[de->GetElementSolutionIndex()]);

  Vector<T>& u1_elem = *u1_elem_ptr;
  Vector<T>& u2_elem = *u2_elem_ptr;

  LOG_DBG3(sc) << "S: de=" << de->elem->elemNum << " a1=" << app1 << " a2=" << app2 << " m=" << mode << " u1=" << u1_elem.ToString() << " u2=" << u2_elem.ToString();

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
  */
}
template<typename T>
void StressConstraint<T>::SetupIntPoint(Elem* elem, unsigned int ip, Mode mode)
{
  Vector<T>& u1_elem = *u1_elem_ptr;
  Vector<T>& u2_elem = *u2_elem_ptr;

  assert(false);
  /* FIXME
  // we need the integration points for B
  Vector<Double>* intPoints = elem->ptElem->GetIntPoints();

  // B1 and B2
  // elem->ptElem->GetCoordMidPoint(intPoint);
  form1->CalcBMatOnly(B1, intPoints[ip-1], elem); // fucking 1-based
  form2->CalcBMatOnly(B2, intPoints[ip-1], elem);
*/
  // left side stress
  strain1 = B1 * u1_elem;
  stress1 = E1 * strain1;
  LOG_DBG3(sc) << "SE de=" << elem->elemNum << " ip=" << ip << " strain1=" << strain1.ToString() << " stress1=" << stress1.ToString();

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
    LOG_DBG3(sc) << "SE de=" << elem->elemNum << " ip=" << ip << " strain2=" << strain2.ToString() << " stress2=" << stress2.ToString();
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
    assert(false);
    /* FIXME
    // one of three piezo case - is the stress constraint defined on a piezo region ?!
    RegionIdType reg = f->elements[0]->elem->regionId;

    BaseForm* form = em->GetForm(reg, Optimization::MECH, Optimization::ELEC, false);
    if(form == NULL)
      throw Exception("piezoelectric stress constraint not defined on a piezoelectric region");
*/
    switch(f->GetStressType())
    {
    case Function::ONLY_COUPLING: // special case only
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
      break;
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

