#ifndef SIMP_HH_
#define SIMP_HH_

#include <map>
#include <set>
#include <string>

#include "Optimization/ErsatzMaterial.hh"
#include "Optimization/Optimization.hh"
#include "PDE/MechPDE.hh"

namespace CoupledField {
class DenseMatrix;
class DesignElement;
class DesignSpace;
class Excitation;
class Function;
class SingleVector;
}  // namespace CoupledField

namespace CoupledField
{
class TransferFunction;
struct SurfElem;
template <class TYPE> class Matrix;


/** In SIMP many sensitivities are with a non-sensitive RHS. This is not
 * necessarily true for PiezoSIMP when we excite with a charge density. But
 * even in pure mech SIMP we can excite with a pressure.
 * Then in CalcU1KU2 for the grad case we actually calc \f$<l,K'u-f'>\f$.
 * This helper stores the information we need to calculate this. The tricky stuff
 * is, that in CalcU1KU2 we generally have volume elements as design elements
 * but the charge density and pressure (inhomogeneus Neumann boundary conditions)
 * are defined on surface elements which are one dimension lower.
 * We assume uniform elements and hence each node values for all surface elements
 * are the same. We store only one node-value which is scalar (potential) or
 * a vector (displacement).
 * We use the design variable from the (volume) element and kind of project it
 * on the rhs which comes from the surface excitation. One has to check all volume
 * nodes if they are part of the surface.
 * Note, that test strain excitation is handled directly by CalcU1KU2() but for homogenization
 * this object is used to indicate the excitation/test_strain. */
class DesignDependentRHS
{
public:
  /** Only MAG needs no init, all other are not valid up to init */
  DesignDependentRHS(App::Type app = App::NO_APP);
  ~DesignDependentRHS();

  /** This is kind of constructor. The return value/status is reflected in valid.
   * @param app is either PRESSURE or App::CHARGE_DENSITY or left when set in constructor
   * @return true if the linear form was found and the variables are init. */
  template <class T>
  bool Init(DesignSpace* design, App::Type app = App::NO_APP);

  /** In this mode the test strain is kept.
   * @param app needs to be STRESS or left when set in constructor
   * @param test_strain taken from the excitation by MechPDE::testStrain.Parse(excitation.label) */
  template <class T>
  bool Init(std::string excite_label, App::Type app = App::NO_APP);

  /** kind of inhomogeneous Neumann. From Init() */
  App::Type app;

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

  bool isInterfaceDriven_;

};


/** Holds a SIMP (Solid Isotropic Material with Penaltization) optimization.
 *  Actually holds the elements of region to optimize, its densities and
 *  global parameters.
 *  This is the single pde version, where typically mech is the pde of choice
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

  /** this calls SetupSecondMaterialCache() for specific MC, MT combinations.
   * Overwrite to change */
  virtual void InitSecondMaterialCache();

  /** this fills DesignSpace::DesignRegion::scnd_material_cached 
   * for all design regions. Type is double, complex, Matrix<double>, Matrix<complex<double> > */
  void AddSecondMaterialCache(MaterialClass, MaterialType);

  /** overwrites the ErsatzMaterial version, is overwritten in PiezoSIMP */
  virtual double CalcFunction(Excitation& excite, Function* f, bool derivative);

  /** This is a helper for CalcU1KU2 to determine the "K" which in most cases include a
   * derivative. It also includes mechanical damping and mass matrix via AddMassToStiffness().
   * The templated stuff is private, as C++ does not allow virtual templates.
   * @param tf for heat and acoustic we canot uniquely identify the transfer function by app therefore give it. */
  virtual void SetElementK(Function* f, DesignElement* de, const TransferFunction* tf, App::Type app, DenseMatrix* out, bool derivative = true, CalcMode calcMode = STANDARD, double ev = -1.0);

  /** the mechanical element rhs, complex or real */
  DesignDependentRHS mechRHS;

private:

  /** This private, as no virtual templates are possible with C++
   * T1 is the out matrix type
   * T2 is the element matrix type */
  template <class T1, class T2>
  void SetElementK(Function* f, DesignElement* de, const TransferFunction* tf, App::Type app, DenseMatrix* out, bool derivative = true, CalcMode mode = STANDARD, double ev = -1.0);

};


} // namespace


#endif /*SIMP_HH_*/
