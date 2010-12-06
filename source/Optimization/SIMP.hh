#ifndef SIMP_HH_
#define SIMP_HH_

#include "Optimization/ErsatzMaterial.hh"
#include "Domain/bcs.hh"
#include "PDE/mechPDE.hh"
#include "MatVec/SingleVector.hh"
#include "MatVec/vector.hh"
#include "MatVec/matrix.hh"
#include <map>
#include <set>

namespace CoupledField
{
class StdPDE;
class SinglePDE;
class ElecPDE;
class BaseForm;
class BiLinFormContext;
class BaseMaterial;
class Condition;
class Assemble;
class TransferFunction;
class SurfElem;
class OptMechMat;

template <class TYPE> class StdVector;
template <class TYPE> class Vector;
template <class TYPE> class Matrix;


/** In SIMP many sensitivities are with a non-sensitive RHS. This is not
 * necessary true for PiezoSIMP when we excite with a charge density. But
 * even in pure mech SIMP we can excite with a pressure.
 * Then in CalcU1KU2 for the grad case we actually calc \f$<l,K'u-f'>\f$.
 * This helper stores the information we need to calculate this. The tricky stuff
 * is, that in CalcU1KU2 we generally have volume elements as design elements
 * but the charge density and pressure (inhomogeneus Neumann boundary conditions)
 * are defined on surface elements which are one dimension lower.
 * We assue uniform elements and hence each node values for all surface elements
 * are the same. We store only one node-value wich is scalar (potential) or
 * a vector (displacement).
 * We use the design variable from the (volume) element and kind of project it
 * on the rhs which comes from the surface excitation. One has to check all volume
 * nodes if they are part of the surface.
 * Note, that test strain excitation is handled directly by CalcU1KU2() but for homogenization
 * this object is used to indicate the excitation/test_strain. */
class DesignDependentRHS
{
public:
  DesignDependentRHS();
  ~DesignDependentRHS();

  /** This is kind of constructor. The return value/status is reflected in valid.
   * @param app is either PRESSURE or CHARGE_DENSITY
   * @return true if the linear form was found and the variables are init. */
  template <class T>
  bool Init(DesignSpace* design, Optimization::Application app);

  /** In this mode the test strain is kept.
   * @param app needs to be STRESS
   * @param test_strain taken from the excitation by MechPDE::testStrain.Parse(excitation.label) */
  template <class T>
  bool Init(Optimization::Application app, std::string excite_label);

  /** kind of inhom Neumbann. From Init() */
  Optimization::Application app;

  /** bool are we set? */
  bool valid;

  /** This is one nodal result (scalar or displacement vector) (real/complex) */
  SingleVector* vec;

  /** This is out reference element */
  const SurfElem*  elem;

  /** this holds the test_strain when the proper Init() was called, otherwise NOT_SET */
  MechPDE::TestStrain test_strain;

  /** This are all node numbers of all surface elements in question. We have
   * to check with the volume element nodes, if they are part of the surface
   * and hence sensitive - otherwise the (volume) RHS does not depend on the
   * design variable. */
  std::set<unsigned int> nodes;

  /** give debug information#
   * @param level 0 is full detail, >0 less detail */
  std::string ToString(int level=0);
};


/** Holds a SIMP (Solid Isotropic Material with Penaltization) optimization.
 *  Actually holds the elements of region to optimize, its densities and
 *  global parameters. Also the cost function and does the looping.
 *  The Reference is Bendsoe, Sigmund; Topology Optimization; Springer Verlag; 2003.
 *  All page numbers refer to this issue. */
class SIMP : public ErsatzMaterial
{
public:
  /** Up to now w/o parameters */
  SIMP();

  /** e.g. closing the exportDesign */
  virtual ~SIMP();


  /** Adds validation stuff here to keep out of long constructor */
  virtual void PostInit();

protected:

  /** overwrites the ErsatzMaterial version, is overwritten in PiezoSIMP */
  virtual double CalcFunction(Excitation& excite, Function* f, bool derivative);

  /** Calculate the Stress gradient. The weight is always 1 as the stress needs to be per excitation */
  void CalcVonMisesStressGradient(Excitation& excite, Function* f,  TransferFunction* tf);

  /** This is a helper for CalcU1KU2 to determine the "K" which in most cases includes a
   * derivative. It also includes mechanical damping and mass matrix via AddMassToStiffness().
   * The templated stuff is private, as C++ does not allow virtual templates. */
  virtual void SetElementK(DesignElement* de, Application app, DenseMatrix* out, CalcMode calcMode, bool derivative = true);

  /** the mechanical element rhs, complex or real */
  DesignDependentRHS mechRHS;

private:

  /** This private, as no virtual templates are possible with C++ */
  template <class T>
  void SetElementK(DesignElement* de, Application app, DenseMatrix* out, CalcMode, bool derivative = true);

  /** This is a helper for SetElementK() which adds for MECH in the harmonic case damping and mass */
  void AddMassToStiffness(double m_factor, DesignElement* de, Matrix<std::complex<double> >& K_in_S_out);

  /** This is a shortcut to ErsatzMaterial::material */
  OptMechMat* mech_mat_;
};


} // namespace


#endif /*SIMP_HH_*/
