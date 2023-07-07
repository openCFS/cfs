#ifndef STRESSCONSTRAINT_HH_
#define STRESSCONSTRAINT_HH_

#include <complex>
#include <utility>

#include "Domain/ElemMapping/EntityLists.hh"
#include "MatVec/Matrix.hh"
#include "MatVec/Vector.hh"
#include "Optimization/ErsatzMaterial.hh"
#include "Optimization/Optimization.hh"
#include "Optimization/TransferFunction.hh"
#include "Utils/StdVector.hh"

namespace CoupledField {
class BaseForm;
class DesignElement;
class DesignSpace;
class Excitation;
class Function;
class SingleVector;
struct ElementAccess;
struct Elem;

/** The calculation of the von Mises Stress function vms, the derivative dvms/drho and the adjoint RHS are
 * all rather involved and very similar.
 * To share the code and handle also the complexity of piezo stresses, the implementation is class based.
 * The implementation is not tuned for speed and data reuse - it is complicated enough!
 *
 * TODO: the code worked in pre-FE-Space. Now, piezo is not enabled any more. To to this, change form -> form1, form2.
 * Some code is still prepared (B1, B2 app.first, app.second, but always same)
 */
template<class T>
class StressConstraint
{
public:
  /** Create the object for each excitation */
  StressConstraint(Excitation* excite, Function* f, ErsatzMaterial* em, StateContainer* forward);

  /** The stress values for every f->elements */
  Vector<double> CalcStresses();

  /** element version for local stress */
  double CalcElementStress(DesignElement* de) {
    int res_idx = space->GetSpecialResultIndex(DesignElement::DEFAULT, DesignElement::QUADRATIC_VM_STRESS, DesignElement::NONE, DesignElement::PLAIN, excite->label);
    return CalcElementStress(STRESS, res_idx, de);
  }

  /** The design gradient of the stress values for every f->elements */
  void CalcGradStresses(Vector<double>& out);

  /** element version for local stress */
  double CalcGradElementStress(DesignElement* de) {
    return CalcElementStress(GRAD_STRESS, -1, de);
  }

  /** globalizes rhs by Calling CalcElemAdjointRHS in a loop
   * @param out will be set, resized and filled */
  void CalcAdjointRHS(Vector<T>& out);

  /** helper for the globalized function and for LOCAL_STRESS and LOCAL_BUCKLING_LOAD_FACTOR
   * @param out_set needs to be a set output rhs where we add our stuff ad the dofs of de */
  void CalcElemAdjointRHS(DesignElement* de, double alpha, Vector<T>& out_set);

  /** is actually the max(0, local_values-c)^2 applied per element. if the bound c is too large it is all zero */
  Vector<double> CalcGlobalizationFactor(const Vector<double>& local_values, bool gradient);

  /** this is a cache for K' * u on element level
   *  for the local stress gradients we have to calculate lambda^T * K' * u,
   *  where only lamda changes for each element */
  struct DKuCache {
    int currentIteration = -100;
    StdVector<Vector<T>> dKu;
  };

  /** return the cache for K' * u */
  DKuCache& GetdKuCache();

private:
  /** This are the three modes of operation */
  typedef enum { STRESS, GRAD_STRESS, ADJOINT_RHS } Mode;

  /** Set up element data which is integration point independent (E1) */
  /** Set up the data for the general formula (E1 ... stress2). E2 is mode dependent, but not u1 */
  void SetupElement(ElementAccess* ea, DesignElement* de, App::Type app1, Mode mode, DesignElement::Type direction = DesignElement::NO_DERIVATIVE);

  /** Set up integration point dependent element data after SetupElement is called!
   * calculates: strain1, B1, B2, ...
   * depending on mode: M_E2_B2, strain2, stress2 */
  void EvalIP(Mode mode, ElementAccess* ea, unsigned int ip);

  /** common for CalcStresses and CalcGradStresses() */
  void CalcStresses(Mode mode, int res_idx, Vector<double>& out);

  /** Element extraction from CalcStresses() */
  double CalcElementStress(Mode mode, int res_idx, DesignElement* de);

  /** Determines if we have mech stress or piezoelectric stress.
   * @return mech/mech or mech/mech, piezo/mech, piezo/piezo, mech/piezo */
  StdVector<std::pair<App::Type, App::Type> > GetApplications();

  /** this is a cache for K' * u, which is used in ErsatzMaterial::CalcLocalVonMisesStressOrLoadFactor */
  static DKuCache dKuCache;

  /** This is the general formula (E1*B1*u1)^T*M*(E2*B2*u2).
   * In the adjoint case u1 = u1^* and u2 is omitted, in the grad case E2 = grad_E2 */
  Matrix<double> M;
  Matrix<double> E1;
  Matrix<double> E2; // in the grad_stress case this is dE2
  Matrix<double> B1;
  Matrix<double> B2;
  Vector<double> alpha; // the globalization factors form GRAD_STRESS and ADJOINT_RHS

  StdVector<SingleVector*>* all_u1_elem = NULL;
  StdVector<SingleVector*>* all_u2_elem = NULL;

  /** This are helper data elements */
  Vector<T> strain1; // B*u
  Vector<T> stress1; // E*B*u

  /** Not calculated in the adjoint case */
  Vector<T> strain2;
  Vector<T> stress2;

  /** Only calculated in the adjoint case */
  Matrix<double> M_E2_B2;

  /** These are set up by SetupElement() to be used in SetupIp(). There was a form2 for the coupling case */
  BaseBDBInt* form = NULL;

  // we need to be careful to use the right index!!
  Vector<T>* u1_elem_ptr = NULL;
  Vector<T>* u2_elem_ptr = NULL;

  Matrix<double> E2_B2; // adjoint case only

  /** technical global stuff */
  TransferFunction tf; // either stress from xml file or implicitly FULL for off-design stresses */
  ElemList elemList;

  Excitation* excite = NULL;
  Function* f = NULL;
  ErsatzMaterial* em = NULL;
  DesignSpace* space = NULL; // shortcut
  StateContainer* forward = NULL;
};

} // end of namespace

#endif /* STRESSCONSTRAINT_HH_ */
