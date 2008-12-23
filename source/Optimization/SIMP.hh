#ifndef SIMP_HH_
#define SIMP_HH_

#include "Optimization/ErsatzMaterial.hh"
#include "Domain/bcs.hh"
#include "Utils/cfsvector.hh"
#include "OLAS/matvec/vector.hh"
#include <map>

namespace CoupledField
{
class StdPDE;
class SinglePDE;
class MechPDE;
class ElecPDE;
class BaseForm;
class BiLinFormContext;
class BaseMaterial;
class Condition;
class Assemble;
class TransferFunction;
class SurfElem;

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
 * nodes if they are part of the surface. */
class SurfaceRef
{
public:
  SurfaceRef();
  ~SurfaceRef();

  /** This is kind of constructor. The return value/status is reflected in valid.
   * @param app is either PRESSURE or CHARGE_DENSITY
   * @return true if the linear form was found and the variables are init. */
  template <class T>
  bool Init(DesignSpace* design, Optimization::Application app);

  /** kind of inhom Neumbann. From Init() */
  Optimization::Application app;

  /** bool are we set? */
  bool valid;

  /** This is one nodal result (scalar or displacement vector) (real/complex) */
  CFSVector* vec;

  /** This is out reference element */
  const SurfElem*  elem;

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

  /** e.g. closing the ersatzMaterialFile */
  virtual ~SIMP();

  /** compute the value of the cost function */
  double CalcObjective()
  {
    return cost->type == RADIATION ? CalcRadiation() : ErsatzMaterial::CalcObjective();
  }

  /** The systems we can SIMP */
  typedef enum { PIEZO, MECHANIC } System;

  /** Here we store the system enum */
  static Enum<System> system;

  System GetSystem() const { return system_; }

  /** Adds validation stuff here to keep out of long constructor */
  virtual void PostInit();

  /** if we do radiation optimization this gives the bilinear form which is
   * added in DefineIntegrators of MechPDE early, before we construct Optimization.
   * It would clearly be nicer if one could add this stuff in SIMP::PostInit()
   * by making OLAS more flexible. Actually this is a low priority todo
   * @return NULL if we do not do such "radiation" optimization stuff */
  static BiLinFormContext* CreateSurfaceNormalMatrix(SinglePDE* pde, BaseMaterial* baseMat, shared_ptr<ResultInfo> result);

protected:

  /** Does mech DENSITY gradients, COMPLIANCE is done in ErsatzMaterial */
  virtual void CalcObjectiveGradient(Excitation& excite);


  /** This is a helper for CalcU1KU2 to determine the "K" which in most cases includes a
   * derivative. It also includes mechanical damping and mass matrix via AddMassToStiffness().
   * The templated stuff is private, as C++ does not allow virtual templates. */
  virtual void SetElementK(DesignElement* de, Application app, CFSMatrix* out);

  /** calculates the radiation, an approximation to emitted sound (Du & Olhoff),
   * only harmomic. Is specific to SIMP, hence not in ErsatzMaterial.
   * Includes the assembly of the SurfaceNormalMatrix */
  double CalcRadiation();

  /** values read from xml in case we optimize for radiation */
  double radiation_c_;
  double radiation_gamma_;

  /** The system we optimize. PIEZO or MECHANIC */
  System system_;

  /** The surface normal matrix following Du/Olhoff 2007. Note that this assumes a single plane surface */
  Matrix<double> surfaceNormal;

  /** the mechanical element rhs, complex or real */
  SurfaceRef mechRHS;

private:

  /** This private, as no virtual templates are possible with C++ */
  template <class T>
  void SetElementK(DesignElement* de, Application app, CFSMatrix* out);

  /** This is a helper for SetElementK() which adds for MECH in the harmonic case damping and mass */
  void AddMassToStiffness(double m_factor, DesignElement* de, Matrix<std::complex<double> >& K_in_S_out);

  /** Specific implementation of for RADIATION, in the other cases calls ErsatzMaterial stuff */
  void AdjustComplexAdjointRHS(Excitation& excite);


  /** Calculates the product of the (system) surface normal matrix with the solution already in OLAS.
   * Note that we have to use 1 based OLAS vectors as the sparse system matrix is from OLAS .
   * This calculation is done for the adjoint rhs and also for calculate the radiation objective.
   * It shall be cheap enough to calc here twice! */
  template <class T>
  void CalcSurfaceNormalTimesSolution(OLAS::Vector<T>& olas_prod, Excitation& excite);

};


} // namespace


#endif /*SIMP_HH_*/
