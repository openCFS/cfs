#ifndef STRESSCONSTRAINT_HH_
#define STRESSCONSTRAINT_HH_

#include "Optimization/ErsatzMaterial.hh"
#include "Optimization/TransferFunction.hh"

namespace CoupledField
{


/** The calculation of the von Mises Stress function vms, the derivative dvms/drho and the adjoint RHS are
 * all rather involved and very similar.
 * To share the code and handle also the complexity of piezo stresses, the implementation is class based.
 * The implementation is not tuned for speed and data reuse - it is complicated enough!
 */
template<class T>
class StressConstraint
{
public:
  /** Create the object for each excitation */
  StressConstraint(Excitation* excite, Function* f, ErsatzMaterial* em, ErsatzMaterial::Solutions* forward);

  /** The stress values for every f->elements */
  void CalcStresses(Vector<double>& out);

  /** The design gradient of the stress values for every f->elements */
  void CalcGradStresses(Vector<double>& out);

  void CalcAdjointRHS(Vector<T>& out);

  /** is actually the max(0, stress-c)^2 applied per element. if the bound c is too large it is all zero */
  void CalcGlobalizationFactor(Vector<double>& out);

private:
  /** This are the three modes of operation */
  typedef enum { STRESS, GRAD_STRESS, ADJOINT_RHS } Mode;

  /** Set up the data for the general formula (E1 ... stress2). E2 is mode dependent, but not u1 */
  void Setup(DesignElement* de, Optimization::Application app1, Optimization::Application app2, Mode mode);

  /** common for CalcStresses and CalcGradStresses() */
  void CalcStresses(Mode mode, int res_idx, Vector<double>& out);

  /** Determines if we have mech stress or piezoelectric stress.
   * @return mech/mech or mech/mech, piezo/mech, piezo/piezo, mech/piezo */
  StdVector<std::pair<Optimization::Application, Optimization::Application> > GetApplications();

  /** This is the general formula E1*B1*u1)^T*M*(E2*B2*u2).
   * In the adjoint case u1 = u1^* and u2 is omitted, in the grad case E2 = grad_E2 */
  Matrix<double> M;
  Matrix<double> E1;
  Matrix<double> E2; // in the grad_stress case this is dE2
  Matrix<double> B1;
  Matrix<double> B2;
  Vector<double> alpha; // the globalization factors form GRAD_STRESS and ADJOINT_RHS

  StdVector<SingleVector*>* all_u1_elem;
  StdVector<SingleVector*>* all_u2_elem;

  /** This are helper data elements */
  Vector<T> strain1; // B*u
  Vector<T> stress1; // E*B*u

  /** Not calculated in the adjoint case */
  Vector<T> strain2;
  Vector<T> stress2;

  /** Only calculated in the adjoint case */
  Matrix<double> M_E2_B2;


  /** technical global stuff */
  TransferFunction tf; // either stress from xml file or implicitly FULL for off-design stresses */
  ElemList elemList;

  Excitation* excite;
  Function* f;
  ErsatzMaterial* em;
  DesignSpace* space; // shortcut
  ErsatzMaterial::Solutions* forward;
};


}

#endif /* STRESSCONSTRAINT_HH_ */
