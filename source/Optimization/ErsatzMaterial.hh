#ifndef ERSATZ_MATERIAL_HH_
#define ERSATZ_MATERIAL_HH_

#include <stddef.h>
#include <complex>
#include <map>
#include <string>
#include <utility>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Domain/BCs.hh"
#include "General/Enum.hh"
#include "General/defs.hh"
#include "General/Environment.hh"
#include "General/Exception.hh"
#include "MatVec/Matrix.hh"
#include "MatVec/Vector.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Function.hh"
#include "Optimization/Optimization.hh"
#include "Optimization/StateSolution.hh"
#include "Utils/StdVector.hh"
#include <tuple>

namespace CoupledField {
class DenseMatrix;
class Excitation;
class Objective;
class SingleVector;
struct Elem;
}  // namespace CoupledField

namespace CoupledField
{
class Assemble;
class Condition;
class DensityFile;
class DesignDependentRHS;
class DesignStructure;
class OptimizationMaterial;
class TransferFunction;
struct SurfElem;


/** Base for optimization where the design variable is correlated to finite elements.
 * The classical case is SIMP, where one optimizes for a pseudo density. The sub-classes
 * extend this idea to more complex stuff.
 * Here the common stuff is kept as solving the adjoint pde and calculating the objective.
 * The implementation of gradients, ... is for the subclasses. */
class ErsatzMaterial: public Optimization
{
public:

  /** initially Solution was a subclass of ErsatzMaterial but then separated */
  friend class StateSolution;

  /** Up to now w/o parameters */
  ErsatzMaterial();

  /** e.g. closing the exportDesign */
  virtual ~ErsatzMaterial();

  /** This solves and stores the forward problem. Eventually the adjoint problem is called
   * implicitly.
   * Processes multiple excitations.
   * @param ev_only_excite EvaluateOnly has the special feature to perform a sweep over multiple
   *        frequencies and calculate the mono-harmonic objective. Only EvaluateOnly shall use this
   *        parameter. */
  virtual void SolveStateProblem(Excitation* ev_only_excite = NULL);

  /** This solves all Adjoint problems */
  void SolveAdjointProblems(Excitation* ev_only_excite = NULL);

  /** Here we also write the density files */
  PtrParamNode CommitIteration();

  /** Adds validation stuff here to keep out of long constructor */
  virtual void PostInit();

  /** Types of ersatz material optimization methods, the strings are read from the xml file */
  typedef enum
  {
    NO_METHOD, SIMP_METHOD, PARAM_MAT, SHAPE_GRAD, SHAPE_OPT, SHAPE_PARAM_MAT, SHAPE_MAP,
    SPAGHETTI, SPAGHETTI_PARAM_MAT, SPLINE_BOX, FEATURE_MAPPING, FEATURE_MAPPING_PARAM_MAT
  } Method;

  /** Actually everything with has material parameterization, e.g. rotated anisotropy as in spaghetti */
  static bool IsParamMat(Method test) { return test == PARAM_MAT || test == SHAPE_PARAM_MAT || test == SPAGHETTI_PARAM_MAT || test == FEATURE_MAPPING_PARAM_MAT; }

  /** True for spaghetti and spaghettiParamMat */
  static bool IsSpaghetti(Method test) { return test ==  SPAGHETTI || test == SPAGHETTI_PARAM_MAT; }

  static Enum<Method> method;

  /** This is the current homogenized tensor.
   * Evaluated by HOM_TRACKING, HOM_TENSOR, HOM_FROBENIUS_PRODUCT (as objective only).
   * MechPDE reads it when "homogenizedTensor" is a region result!
   * It is a vector for the variants for transformation and robust */
  StdVector<Matrix<double> > homogenizedTensor;

  DensityFile* GetDensityFile() { return densityFile; }

  Method GetMethod() { return method_; }

  StateContainer& GetForwardStates() { return forward; }

  inline const DesignStructure* GetDesignStructure(){
    assert(structure_ != NULL);
    return structure_;
  }

  // calculates for given node res = interface * (stateSol - trackVal)^2
  double CalcStateTrackingAtNode(int node);

  // calculates for given node res = load * stateSol
  double CalcTempAtInterface(int node);

  void CalcStressesForBucklingHomogenization(Matrix<Double>& S, const LocPointMapped* lpm);

  /** this is the optimization->ersatzMaterial XML element */
  PtrParamNode pn;

  /** Here we store the solution of the problem. Multiple solutions for multiple loadcases */
  StateContainer forward;

  /** Here we store the solution of the adjoint problem. */
  StateContainer adjoint;

  /** Calculates a local load factor for normalized stress by evaluating
   * an interpolated function (stored in the material catalogue file).
   * Also does a quadratic extrapolation to avoid too large values for vol -> 1. */
  double GetMicroLoadFactor(double vol, bool derivative = false);

protected:
  
  /** overwrites the Optimization version. To be called within Optimization::CommitIteration() */
  virtual void LogFileLine(std::ofstream* out, PtrParamNode iteration);

  /** When "commit" is set, we write "forward"/"adjoint" or "both_cases" */
  virtual void StoreResults(double step_val);

  /** @see Optimization::GetIterationFrequency() */
  std::string GetIterationFrequency();

  /** This are the modes for CalcU1KU2(). */
  enum CalcMode
  {
    STANDARD = 0, /*!< add u1^T (K' u2  - f') or2 * Re{ u1^T (K' u2 - f')} in the harmonic case  */
    CONJ_QUAD,
    EIGENFREQ,    /*!< The derivative is <u, (K' - ev M') u> which implies CONJ_QUAD */
    BUCKLING      /*!< The derivative is <u, (K' - ev G') u> */
  };

  /*!< add <u, K' u> which is in the real case as STANDARD
   and for the harmonic case u^T K' u^* (conj. complex). u1 = u2 = u!! */
  /** Calculate the sum of  \f$ l^T K'u - f'\f$ or \f$ 2 Re{l^T K'u} - f'\f$ or \f$ <K'l,u> - f'\f$.
   * This is controlled via CalcMode. 

   *
   * When adjoint vector u1/l is not calculated with a negative rhs, one
   * can put in the minus sign as an explicit factor.
   *
   * The K is automatically determined by the application and the size of
   * the vector u(2).
   *
   * @param tf the transfer function defines the design variable and multiplier of K
   * @param u1 for derivatives the lagrange vector l
   * @param k the application determines the stiffness matrix
   * @param u2 the solution or u in \f$<l,K'u-f'>\f$
   * @param rhs if one want to do \f$<l,K'u-f'>\f$ this contains the info for \f$-f'\f$.
   * @param factor see above, more complex in radiation case.
   * @param calcMode how to solve the product.
   * @param res_idx store in de->specialResult. use ErsatzMaterial::GetSpecialResultIndex() -1 is no special result
   * @param ev in the eigenvalue case factor for M' which is the eigenvalue. Then calcMode is EIGENVALUE */
  double CalcU1KU2(TransferFunction* tf, StdVector<SingleVector*>& u1,
      App::Type k, StdVector<SingleVector*>& u2, DesignDependentRHS* rhs,
      double factor, CalcMode calcMode, Function* f, int res_idx = -1, double ev = -1.0);

  /** Helper calling CalcU1KU2()
   * If there is a result with value='costGradient' or 'constraintGradient' it is checked for detail='mech_mech',
   * 'elec_elec', 'elec_elec_quad', 'elec_mech', 'mech_elec' */
  int GetSpecialResultIndex(App::Type app1, App::Type app2,
      CalcMode calcMode = STANDARD, Condition* constraint = NULL);

  /** This is a helper for CalcU1KU2 to determine the "K" which in most cases includes a
   * derivative. It also includes mechanical damping and mass matrix via AddMassToStiffness().
   * Also bi-material is nicely considered.
   * The template stuff is private, as C++ does not allow virtual templates. */
  virtual void SetElementK(Function* f, DesignElement* de, const TransferFunction* tf, App::Type app, DenseMatrix* out, bool derivative = true, CalcMode mode = STANDARD, double ev = -1.0)
  {
    throw Exception("not implemented");
  }

  /** This is a helper for CalcU1KU2() and called when CalcU1KU2 got a DesignDependentRHS* rhs value. To be specified form MagSIMP, ...
   * It substracts a from mat_vec */
  virtual void SubstractCalcU1KU2RHS(Function* f, TransferFunction* tf, DesignElement* de, DesignDependentRHS* rhs, SingleVector* mat_vec);

  /** Helper that asks MechanicMaterial. Works only for a single region.
   * works also with multiple regions if grid is regular
   * @param excitation we need to make sure the excitation is the active one for robust
   * @return empty if multiple regions */
  StdVector<std::pair<std::string, double> > GetOrthotropeProperties(const Matrix<double>& tensor, Excitation* ex);

  /** This is an extension to SolveStateProblem() where the forward problem is solved and stored.
   * Depending on the objective function SolveAdjointProblem() is called to additionally solve and store the
   * adjoint problem.
   * It works for both (mechanical) SIMP and PiezoSIMP. 
   * gradient is used to calculate some adjoints only for gradient calculations, some for function evaluations */
  void SolveAdjointProblem(Excitation* excite, Function* f);

  /** Determines the selection vector by a "pseudo loading" for output like objectives.
   * Stores in adjoint.select. Used by SolveAdjointProblem()
   * @param alter_rsh false if you want only selection and the system shall not be changed! */
  void ConstructSelection(Excitation& excite, Function* f, bool alter_rhs);

  /** This is helper SolveAdjointProblem().
   * There is a template method (which cannot be virtual) with distinct implementation.
   * Assumes adjont.select is set (by ConstructSelection())
   * This is for output loads or general real/complex rhs. */
  virtual void ConstructAdjointRHS(Excitation& excite, Function* f);

  /** A simple variant of IntegrateDesignVariable() which works also for non-simp transfer function.
   * Handles also non-regular and physical
   * @param normalized @see CalcVolume() */
  double CalcTrivialVolume(Function* f, bool derivative, bool normalized);

  /** calculates the integral over a design variable (note that volume is a special case of this (with all standard values) @see CalcVolume
   * regularization is using this usually with normalized = true, scale = true, square = true, factor = "the regularization parameter"
   * @param dtype The design variable to integrate over, use TENSOR_TRACE for tensor trace or DEFAULT for single design
   * @param derivative calculate the derivative instead of the integral
   * @param f set if objective otherwise NULL
   * @param g as f but constraint
   * @param normalized see CalcVolume
   * @param scale calculate the integral over the scaled design variable (<design scale="true">), i.e. the variable as in the optimizer
   * @param square calculate the integral over the square of the DesignVariable */
  virtual double IntegrateDesignVariable(Objective* f, Condition* g,
      bool derivative, DesignElement::Type dtype, bool normalized = true,
      bool scale = false, double exponent = 1.0);

  /** Handles the Volume constraint. Has a constraint and constraint derivative mode
   * @param derivative if false the return value is calculated. Otherwise the value in
   *                   the design element is set.
   * @param constraint if set, calculate as given constraint, if null calculate as objective
   * @param grad_out if derivative is set and grad_out is not null it is set.
   * @param normalized if set use normalized "volume" i.e. sum(design * area) / sum(area) for each element, else only sum(design * area)
   * @return invalid in derivative case*/
  virtual double CalcVolume(Objective* f, Condition* g, bool derivative,
      bool normalized = true);

  /** Handles the Compliance constraint/objective. Has a objective, objective derivative,
   * constraint and constraint derivative mode
   * @param excite for multiple excitations, defines which solution and rhs to use
   *               (default the first any only one)
   * @param derivative if false the return value is calculated. Otherwise the value in
   *                   the design element is set.
   * @param constraint if set calculate as given constraint, if null calculate as objective
   * @param grad_out if derivative is set and grad_out is not null it is set.
   * @return invalid in derivative case*/
  virtual double CalcCompliance(Excitation& excite, Objective* f, Condition* g, bool derivative);

  /** Calculates the objective only, no derivative */
  double CalcGlobalDynamicCompliance(Excitation& excite, Function* f);

  /** Calculates heat energy as an equivalence to compliance in lin elasticity */
  virtual double CalcHeatEnergy(Excitation& excite, Objective* f, Condition* g, bool derivative);

  /** call a python function for evaluation. Implementation in PythonFunction.cc */
  double CalcPython(Excitation& excite, Function* f, bool derivative);

  /** Calculates <l,u> or <conj(u) L, u> where l/L is adjoint[idx]->rhs */
  template<class T>
  double CalcOutput(Excitation& excite, Function* f);

  /** Handles the Tracking constraint/objective. Has a objective, objective derivative, 
   * constraint and constraint derivative mode
   * @param excite  the excitation used 
   * @param derivative if false the return value is calculated. Otherwise the value in
   *                   the design element is set.
   * @param constraint if set calculate as given constraint, if null calculate as objective
   * @param solveproblem solve the tracking problem, e.g. shapeopt does solve the same problem already
   * @return invalid in derivative case*/
  virtual double CalcTracking(Excitation& excite, Objective* f, Condition* g, bool derivative);

  /** Handles tracking constraint/objective for given temperature at interfaces between solid and void
   *  @param excite The used excitation
   *  @param f function
   *  @param derivative flag for calculating derivative
   *  @param trackVal the value that we want to track
   *  @return sum over all tracked values at interface nodes
   */
  virtual double CalcStateTrackingAtInterface(Excitation& excite, Function* f, bool derivative, double trackVal);

  /**
   * Calculates and sets adjoint rhs for temperature (or any scalar state) at interfaces between solid an void
   * K*l^T = -2 * F' * (u - u_track)
   */
  virtual void CalcAdjointRHSStateTracking(Excitation& excite, Function* f, double trackVal, Vector<double>& out);

  /** magnetic flux density */
  void CalcMagFluxAdjRHS(Excitation& excite, Function* f, Vector<double>& out);

  /** Calculate the energy flux through a surface region: 1/2*Re{j*u^T Q u^*} where
   * Q is the grad operator in z direction. Only for acoustic but easy to extend!*/
  double CalcEnergyFlux(Excitation& excite, Objective* f);

  /** at the moment meant for param <= alpha +/- slack. To be extended by mathparser when required */
  double CalcExpression(Condition* g, bool derivative);

  /** Standard Eigenfrequency (ev -> ef) problem.
   * This problem is better scaled than the eigenvalue problem and it matches eigenfrequency output
   * @param excite the wave vector for the function when bloch=full or the last for bloch=extremal*/
  double CalcEigenfrequency(Excitation& excite, Function* f, bool derivative);

  void CalcEigenvalueDerivativeBuckling(Excitation& excite, Function* f, StateSolution* sol, Double ev);

  /** calculates the global stresses and the gradients. The weight is always 1 as the stress needs to be per excitation */
  template<class TYPE>
  double CalcGlobalVonMisesStressOrLoadFactor(Excitation& excite, Function* f, bool gradient);

  /** calculates local stress value or local microscopic load factor or gradient for any of those two.
   * We do not go via LocalFunction here! */
  double CalcLocalVonMisesStressOrLoadFactor(Excitation& excite, Function* f, bool gradient);

  /** bandgap is the difference between two eigenfrequency problems in the bloch mode.
   * It would make sense to have a generic gap function between two independent functions */
  double CalcBandGap(Excitation& excite, Function* f, bool derivative);

  /** calculates the variants of slack functions */
  double CalcSlackFunction(Function* f, bool derivative);

  /** This is a helper with the common part for CalcEnergyFlux and the adjoint RHS.
   * Determines the global vector Q*u^* or (Q - Q^T)^T*u^* in the adjoint case.
   * @param f the cost function as we need the ParamNode
   * @param u_glob the complete equation index global solution vector
   * @param adjoint switch action
   * @param q_u_glob output depending on adjoint flag. see u_glob */
  void SetEnergyFluxVector(Function* f, const Vector<std::complex<double> >& u_glob, bool adjoint, Vector<std::complex<double> >& q_u_glob);

  /** Find the node numbers which are common from a surface element and a volume element.
   * This maps from a surface element to the volume element.
   * @param common_nodes e.g. a center surface element node is lost on a 20-hex-volume element */
  void FindCommonNodes(const SurfElem* se, const Elem* vol,  StdVector<unsigned int>& common_nodes) const;

  /** does the substep of solving K z = Proj(u - u0) for z */
  void SolveTrackingProblem(Excitation& excite, bool designelem = true, bool gridelem = false);

  /** converts the teststrain vector in voigt notation to the corresponding matrix
   * @param app application type (PDE)
   * @param matrix output
   * @param vec input */
  void SetTestStrainMatrix(App::Type app, Matrix<double>& matrix, const Vector<double>& vec);

  /** takes the result of the test strain computations and calculates the homogenized 
   *  material tensor (see Bendsoe/Sigmund: Topology Optimization, p. 122ff.
   *  It must be called only for the last excitation when all test strains are known.
   *  @param f to handle transformaiton and robustness. Extracts the excitation from f */
  virtual Matrix<double> CalcHomogenizedTensor(Function* f);

  /** Calculates the gradient of the homogenization tracking. When J = 0.5 * || E^* - E^H ||^2
   * the this calulates -1 (E^* - E^H) * d(E^H)/d(rho_e) using a matrix scalar product
   * @param target E^* what we want
   * @param hom the pre calculated tensor E^H */
  virtual void CalcHomogenizedTrackingGradient(const Matrix<double>& target, const Matrix<double>& hom, Function* f);

  /** Calculates the gradient for the Frobenius inner prodcut. */
  void CalcHomFrobeniusProductGradient(const Matrix<double>& target,  const Matrix<double>& hom, Function* f);

  /** Calculates the gradient if the constraints E^H = E^* where for each interested
   * tensor entry a own HOM_TENSOR constraint is required.
   * @param derivative this sets d(E^H)/d(rho_e) for the current tensor entry
   * @param g the constraint is mandatory. It defines in the coord pair the tensor entry
   * @return the E^H tensor entry if !derivative */
  double CalcHomogenizedTensorConstraint(Condition* g, bool derivative);

  /** Calculates a single homogenized tensor entry or it's derivative.
   * This is a helper for CalcHomogenizedTensorConstraint() but more general.
   * CalcHomogenizedTensor() could multiple times calls this but this would be slower.
   * @param entry the 1-based!!!! entry within the tensor
   * @param derivative this sets d(E^H)/d(rho_e) for the current tensor entry
   * @param out_grad of derivative it is resized and the gradients are set otherwise it is untouched
   * @param meta the meta excitation index (rotations, robust) or 0 for standard case
   * @return the E^H tensor entry if !derivative or 0 */
  double CalcHomogenizedTensorEntry(Function* f, const std::tuple<int, int, double> entry, bool derivative, StdVector<double>& grad_out, unsigned int meta);

  /** Calculates globalized local functions: globalSlope and globalCheckerboard.
   * When g_i is the slope function x_i - x_i+1 -c and g_i+1 = x_1+1 - x_i - c
   * the global slope is sum max(0, g_i)^2, hence we need NEXT_AND_REVERSE locality */
   double CalcGlobalFunction(Excitation& excite, Function* f, bool derivative);

  /** Evaluates objective and constraint functiond and gradient.
   * Overloaded in PiezoSIMP for own objectives.
   * @param grad_out only used in derivative case
   * @return zero for derivative */
  virtual double CalcFunction(Excitation& excite, Function* f, bool derivative);

  /** Store the results from the forward/adjoint problem. Handles multiple excitations
   * @param timestep_mode_local transient time step, bloch mode oder local stress virtual index
   * @param read_sol store solution (maybe one would only like to save rhs)
   * @param read_rhs is only interesting for the forward problem
   * @param save_sol set this in the adjoint problem -> see Solution::Read()
   * @param comment is just to LOG_DBG */
  virtual void StorePDESolution(StateContainer& solutions, Excitation& excite,
      Function* f, int timestep_mode_local, bool read_sol, bool read_rhs,
      bool save_sol, TimeDeriv derivative, const std::string& comment);

  // virtual void TimeStepCalculated(UInt timeStep, AdjointParameters* adjParams);
  // virtual void RhsCalculated(AdjointParameters* adjParams);

  /** Helper that gives the physical material tensor considers bi-material */
  void GetPhysicalMaterial(BiLinForm* form, const DesignElement* de, const TransferFunction* tf, bool derivative, Matrix<double>& out);

  /** This is a helper for SetElementK() which adds for App::MECH in the harmonic case damping and mass
   * @param bimaterial describes only the material, the factor needs to be set as rho^3 or 1-rho^3 already! */
  void AddMassToStiffness(Context* ctxt, const TransferFunction* mtf, DesignElement* de, Matrix<std::complex<double> >& K_in_S_out, bool derivative, bool bimaterial, CalcMode mode = STANDARD, double ev = -1.0);

  /** This is a helper for SetElementK() which adds for App::MECH in the buckling case geometric stiffness
   * @param bimaterial describes only the material, the factor needs to be set as rho^3 or 1-rho^3 already! */
  void AddGeometricStiffnessToStiffness(Context* ctxt, const TransferFunction* tf, DesignElement* de, Matrix<Complex>& K_in_S_out, bool derivative, bool bimaterial, CalcMode mode = STANDARD, double ev = -1.0);

  /** For derived optimization to fill their contribution to ErsatzMaterial::ConstructRealAdjointRHS()
   * @return true if function is handled */
  virtual bool FillRealAdjointRHS(Excitation& excite, Function* f, Vector<double>& rhs) { return false; }

  /** For derived optimization to fill their contribution to ErsatzMaterial::ConstructComplexAdjointRHS()
   * @return true if function is handled */
  virtual bool FillComplexAdjointRHS(Excitation& excite, Function* f, Vector<Complex>& rhs) { return false; }

  /** Obtain the complex pressure amplitude of the excitation, e.g. caused by the normalVelocity b.c.
   * @param f Current function.
   * @return Nodal (constant) pressure of the incident wave. */
  virtual const Complex GetExcitationPressure(Function* f) { throw Exception("not implemented"); }

  /** Obtain the excitation surfRegion region
   * @param f Current function.
   * @param attr The name of the region attribute ("surfRegion" by default, can be "volumeNeighbour").
   * @return The retrieved region ID.*/
  virtual const RegionIdType GetExcitationRegion(Function* f, const std::string& attr = "surfRegion") { throw Exception("not implemented"); }

  /** do we do SIMP or FreeMat or ... */
  Method method_;

  /** cache the 1.0 / complete volume of the domain */
  double volumeFraction_;

  /** The DesignStructure is required by SIMP for filters and by Condition for slope constraints
   * and checkerboard. They share this element. It can only be created by PostInit(), hence every
   * PostInit() who needs the structure needs to check if it was created before. Deleted by ~EM */
  DesignStructure* structure_ = NULL;

  /** This is just a shortcut for the actual dimensions (2 or 3) */
  const unsigned int dim;

  /** Convenience class for writing the pseudo density file*/
  DensityFile* densityFile;

  /** do we perform homogenization for both permittivity and permeability? */
  bool bitensor_;

private:

  /** This is a helper for the calculation of the homogenized tensor or the derivative of it.
   * This is the inner of the sum for the homogenized tensor or the derivative formulation
   * in Bendsoe/Sigmund - Topology Optimization page 124
   * @param u1 the element solution vector
   * @return the product test strain diff * (K or K') * test strain diff  */
  static double CalcHomogenizedElementProduct(ErsatzMaterial* obj, Function* f, DesignElement* de, bool derivative, Vector<double>& u1, Vector<double>& u2, UInt ij, UInt kl);

  static Complex CalcU1KU2(ErsatzMaterial* obj, DesignElement* de, bool derivative, Vector<Complex> u1_vec, Vector<Complex> u2_vec);

  /** See the non-template version for documentation! */
  template<class T>
  double CalcU1KU2(TransferFunction* tf,
      StdVector<SingleVector*>& u1, App::Type k, StdVector<SingleVector*>& u2,
      DesignDependentRHS* ref, double factor, CalcMode calcMode, Function* f,
      int res_idx, double ev);

  /** for design dependent interface driven excitation as for heat .... calculate element based f'
   * for all design elements de:
   *  run over all neighbor nodes of design de
   *  f'=4*d_rho_i/d_rho_j *(1-2*rho_i), where rho_i is the node based density calculated via averaging the densities of neighboring elements
   * @param function f */
  template<class T>
  void CalcAndStoreInterfaceDrivenGrad(Function* f, TransferFunction* tf);

  template<class T>
  void SubstractInterfaceDrivenGradRHS(Function* f, TransferFunction* tf, const DesignElement* de, Vector<T>& in_out);

  template<class T>
  void SubstractCalcU1KU2RHS(Function* f, TransferFunction* tf, DesignElement* de, DesignDependentRHS* rhs, Vector<T>& mat_vec);

  /** Handles sensitive RHS, e.g. when we have sensitive Neuman boundary condition (elect surface charge).
   * SurfaceRef is  given to CalcU1KU2 and this method does from \f$<l,K'u-f'>\f$ the \f$-f'\f$ part.
   * It checks if any nodes of the design element are part of the surface and
   * substracts for all dof of that node only */
  template<class T> void SubtractGradSurfaceRHS(DesignElement* de, TransferFunction* tf, DesignDependentRHS* ref, Vector<T>& in_out);

  /** Called by CalcU1KU2 on IfStrainExcitedSystem()
   * @param rhs might be null but if not takes the test_strain. */
  template<class T> void SubtractGradStrainRHS(DesignElement* de, TransferFunction* tf, DesignDependentRHS* rhs, Vector<T>& in_out);

  /** Calculates a scalar product of two vectors and the derivative of the right hand side newmark update,
   * used for transient optimization derivative calculation
   * @param excite excitation to consider
   * @param forward first vector
   * @param adjoint second vector
   * @param factor factor to multiply the value by (can be excitation weight)
   * @param f objective the result is to be stored with
   * @param g constraint the result is to be stored with */
  void CalcNewmarkDerivative(Excitation& excite, StateContainer& forward, StateContainer& adjoint, double factor, Objective* f, Condition* g);

  /** This solves the adjoint problem problem only and stores all relevant data. Calls SetAndSolveAdjointRHS().
   * For local functions every adjoint is stored at index 0 and overwrites the previous one! */
  template<class T> void SolveAdjointProblem(Excitation* excite, Function* f);

  /** Set the rhs for the adjoint equation, called by assemble */
  // virtual void SetAdjointRhs(AdjointParameters* adjointParams);

  /** Set the rhs for the tracking adjoint */
  // void SetTrackingAdjointRhs(Excitation& excite, int ts);

  /** Takes care for making CFS solving the adjoint PDE. Sets the rhs as  adjoint[excite.index]->rhs[App::MECH] */
  void SetAndSolveAdjointRHS(Excitation& excite,  Function* cost);

  /** return value of PrepareAdjointSystem */
  struct SystemState
  {
    StdVector<LinearFormContext*> forms;
    IdBcList  idbc;
  };

  /** Preparing the adjoint system constructs the rhs and the system (IDBC- > HDBC).
   * It needs afterwards to be solved by assemble->GetAlgSys()->Solve()
   * Then the adjoint state can be read and after this the system can be resored to the state system */
  SystemState PrepareAdjointSystem(Excitation& excite, Function* f);

  /** restores state system after PrepareAdjointSystem(), solve and StateSolution::Read(). To have IDBC->HDBC in the adjoint
   * requires first StateSolution::Read() before we reset proper IDBC. */
  void RestoreStateSystem(SystemState& sys);

  /** Helper for CommitIteration. Appends or replaces a design line */
  void SetCurrentExportDesign();

  /** In ErsatzMaterial: Saves the original loads and sets the output nodes as adjoint pde rhs
   * Has distinct implementations for complex and real part.
   * @see virtual ConstructAdjointRHS()
   * @see FillReadAdjointRHS() ! */
  void ConstructRealAdjointRHS(Excitation& excite, Function* f);
  void ConstructComplexAdjointRHS(Excitation& excite, Function* f);

  /** construct the right hand side for the adjoint equation of the buckling problem
   * The adjoint system is \f$Kv = phi^T \frac{\partial G(u)}{\partial u} phi\f$. */
  void ConstructAdjointRHSBuckling(Function* f, Vector<Complex>& mode, Vector<Complex>& rhs);

  /** Calculates the Grayness. Works as constraint extremely poor, is more for observation or penalty in multi objective.
   * non physical: g(x) = 4 * x * (1 - x) with x = (rho-lb) / (ub-lb)
   * physical: g(x) = s * t(x) * t(1-x) 
   * g(0.5) := 1, s is accordingly found (see gray_scale in .info.xml) */
  double CalcGreyness(Function* f, bool derivative);

  /** Calculates the difference between filtered and non-filtered stiffness tensor.
   * @param derivative if false the return value is calculated. Otherwise the value in
   *        the design element is set.
   * */
  double CalcFilteringGap(Condition* g, bool derivative);

  /** Evaluates virtually blown up local constraints based on the Function::Local neighborhood.
   * E.g. slope and mole. Note, that there are also the globalized variants.
   * @see CalcGlobalFunction() */
  double CalcLocalConstraint(Condition* g, bool derivative);

  /** IntegrateDesignVariables() can do a lot, but no one wants to extend it to hande the derivative
   * case of the gap constraint: volume - penalized volume */
  void CalcRegularGapConstraint(Function* f, DesignElement::Type dt);

  /** The design tracking matches a given pattern of physcial design by a certain extend (e.g. 90 %). T
   * The designTarget attribute is mandatory. If the region attribute is not given the periodic border
   * elements are used. If the region is given it might make sense to use the pattern definition in within the domain. */
  double CalcDesignTracking(Condition* g, bool derivative);

  /** Homogenization objective/ constraint.
   * Is once evaluate only! */
  double CalcPoissonsRatioAndYoungsModulus(Function* f, bool derivative);

  /** for CalcFunction() */
  double CalcHomTensor(Objective* c, Condition* g, bool derivative);

  /** Calculates the product of the (system) surface normal matrix with the solution already in OLAS.
   * Note that we have to use 1 based OLAS vectors as the sparse system matrix is from OLAS .
   * This calculation is done for the adjoint rhs and also for calculate the radiation objective.
   * It shall be cheap enough to calc here twice! */
  template<class T>
  void CalcSurfaceNormalTimesSolution(Vector<T>& olas_prod);

  /** Get the prescribed macro stress from the xml for buckling homogenization. */
  Vector<Double> GetMacroStress();

  /** Calculates a local load factor for a given stress or the derivative */
  double CalcLoadFactor(Excitation& excite, Function* f, DesignElement* de, double stress, bool gradient = false, Vector<double>* alpha = nullptr);

  /** Have we already calculated gradient of interface driven load gradient for each design element?*/
  bool interfaceDrivenGradCalc_ = false;

  Function* trackingFunc_ = NULL;

  bool printProgressBar_ = false;

  shared_ptr<Timer> calc_u1ku2_timer_;

};

} // namespace


#endif /*ERSATZ_MATERIAL_HH_*/
