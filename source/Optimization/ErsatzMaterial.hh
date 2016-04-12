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
#include "boost/tuple/tuple.hpp"

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
class SinglePDE;
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
  PtrParamNode CommitIteration(bool keep_iteration_number = false);

  /** Adds validation stuff here to keep out of long constructor */
  virtual void PostInit();

  /** Have all design elements the same size? -> same local element matrices */
  bool IsDomainStructured()
  {
    return assume_constant_element_matrices_;
  }

  /** Types of ersatz material optimization methods, the strings are read from the xml file */
  typedef enum
  {
    NO_METHOD, SIMP_METHOD, PARAM_MAT, SHAPE_GRAD, SHAPE_OPT, SHAPE_PARAM_MAT
  } Method;

  static Enum<Method> method;

  /** This is the current homogenized tensor.
   * Evaluated by HOM_TRACKING, HOM_TENSOR, HOM_FROBENIUS_PRODUCT (as objective only).
   * MechPDE reads it when "homogenizedTensor" is a region result!
   * It is a vector for the variants for transformation and robust */
  StdVector<Matrix<double> > homogenizedTensor;

  DensityFile* GetDensityFile() { return densityFile; }

  Method GetMethod() { return method_; }

  StateContainer& GetForwardStates() { return forward; }

  /** this is the optimization->ersatzMaterial XML element */
  PtrParamNode pn;

  inline const DesignStructure& GetDesignStructure(){
    assert(structure_ != NULL);
    return *structure_;
  }

protected:
  
  /** When "commit" is set, we write "forward"/"adjoint" or "both_cases" */
  virtual void StoreResults(double step_val);

  /** @see Optimization::GetIterationFrequency() */
  std::string GetIterationFrequency();

  /** This are the modes for CalcU1KU2(). */
  enum CalcMode
  {
    STANDARD = 0, /*!< add u1^T (K' u2  - f') or2 * Re{ u1^T (K' u2 - f')} in the harmonic case  */
    CONJ_QUAD,
    EIGENFREQ    /*!< The derivative is <u, (K' - ev M') u> which implies CONJ_QUAD */
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



  /** for Virgininies stuff */
  double CalcU1KU2_mapping(TransferFunction* tf, StdVector<SingleVector*>& u1,
      App::Type k, StdVector<SingleVector*>& u2, DesignDependentRHS* rhs,
      double factor, CalcMode calcMode, Function* f, int res_idx = -1);

  /** for Virgininies stuff */
  double CalcU1KU2_mapping2(TransferFunction* tf, StdVector<SingleVector*>& u1,
      App::Type k, StdVector<SingleVector*>& u2, DesignDependentRHS* rhs,
      double factor, CalcMode calcMode, Function* f, int res_idx = -1);


  /** Helper calling CalcU1KU2()
   * If there is a result with value='costGradient' or 'constraintGradient' it is checked for detail='mech_mech',
   * 'elec_elec', 'elec_elec_quad', 'elec_mech', 'mech_elec' */
  int GetSpecialResultIndex(App::Type app1, App::Type app2,
      CalcMode calcMode = STANDARD, Condition* constraint = NULL);

  /** This is a helper for CalcU1KU2 to determine the "K" which in most cases includes a
   * derivative. It also includes mechanical damping and mass matrix via AddMassToStiffness().
   * Also bi-material is nicely considered.
   * The template stuff is private, as C++ does not allow virtual templates. */
  virtual void SetElementK(Context* ctxt, DesignElement* de, const TransferFunction* tf, App::Type app, DenseMatrix* out, bool derivative = true, CalcMode mode = STANDARD, double ev = -1.0)
  {
    throw Exception("not implemented");
  }


  /** this is for Virginies stuff */
  virtual void SetElementKMapping(DesignElement* de, BaseDesignElement::Type type, const TransferFunction* tf, App::Type app, DenseMatrix* out, CalcMode calcMode, bool derivative = true)
  {
    throw Exception("not implemented");
  }


  virtual void SetElementRHS(DesignElement* de, const TransferFunction* tf,
      App::Type app, SingleVector* out, CalcMode calcMode, bool derivative =
          true)
  {
    throw Exception("not implemented");
  }

  /** Helper that asks MechanicMaterial. Works only for a single region.
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
  virtual double CalcCompliance(Excitation& excite, Objective* f, Condition* g,
      bool derivative);
  /** Calculates the objective only, no derivative */
  double CalcGlobalDynamicCompliance(Excitation& excite, Function* f);
  /** Calculates heat energy as an equivalence to compliance in lin elasticity */
  virtual double CalcHeatEnergy(Excitation& excite, Objective* f, Condition* g, bool derivative);
  /** Calculates <l,u> or <conj(u) L, u> where l/L is adjoint[idx]->rhs */
  template<class T> double CalcOutput(Excitation& excite, Function* f);
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
   *  @param f objectiv
   *  @param g condition
   *  @param derivative flag for calculating derivative
   *  @param refVal the reference value that we want to track
   *  @param adjoint Are we calculating adjoint rhs?
   *  @param adjointRHS If adjoint RHS should be calculated, this is the output
   *  @return sum over all tracked values at interface nodes
   */
  virtual void CalcTempTrackingAtInterface(Excitation& excite, Objective* f, Condition* g, double refVal, Vector<double>& res, bool adjoint=false);

  /** Calculate the energy flux through a surface region: 1/2*Re{j*u^T Q u^*} where
   * Q is the grad operator in z direction. Only for acoustic but easy to extend!*/
  double CalcEnergyFlux(Excitation& excite, Objective* f);

  /** at the moment meant for param <= alpha +/- slack. To be extended by mathparser when required */
  double CalcExpression(Condition* g, bool derivative);

  /** Standard Eigenfrequency (ev -> ef) problem.
   * This problem is better scaled than the eigenvalue problem and it matches eigenfrequency output
   * @param excite the wave vector for the function when bloch=full or the last for bloch=extremal*/
  double CalcEigenfrequency(Excitation& excite, Function* f, bool derivative);

  /** bandgap is the difference between two eigenfrequency problems in the bloch mode.
   * It would make sense to have a generic gap function between two independent functions */
  double CalcBandGap(Excitation& excite, Function* f, bool derivative);

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
   * @param matrix output
   * @param vec input */
  void SetTestStrainMatrix(Matrix<double>& matrix, const Vector<double>& vec);

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
  double CalcHomogenizedTensorEntry(Context* ctxt, const boost::tuple<int, int, double> entry, bool derivative, StdVector<double>& grad_out, unsigned int meta);

  /** Calculates globalized local functions. globalSlope and globalCheckerboard.
   * When g_i is the slope function x_i - x_i+1 -c and g_i+1 = x_1+1 - x_i - c
   * the global slope is sum max(0, g_i)^2, hence we need NEXT_AND_REVERSE locality
   * @param von_mises_stress set only for f == STRESS for derivative and not derivative */
  double CalcGlobalFunction(Function* f, bool derivative,  const Vector<double>* von_mises_stress = NULL);

  /** Evaluates objective and constraint functiond and gradient.
   * Overloaded in PiezoSIMP for own objectives.
   * @param grad_out only used in derivative case
   * @return zero for derivative */
  virtual double CalcFunction(Excitation& excite, Function* f, bool derivative);
  /** Store the results from the forward/adjoint problem. Handles multiple excitations
   * @param read_sol store solution (maybe one would only like to save rhs)
   * @param read_rhs is only interesting for the forward problem
   * @param save_sol set this in the adjoint problem -> see Solution::Read()
   * @param comment is just to LOG_DBG */
  virtual void StorePDESolution(StateContainer& solutions, Excitation& excite,
      Function* f, int timestep_mode, bool read_sol, bool read_rhs,
      bool save_sol, TimeDeriv derivative, const std::string& comment);

  // virtual void TimeStepCalculated(UInt timeStep, AdjointParameters* adjParams);
  // virtual void RhsCalculated(AdjointParameters* adjParams);

  /** Helper that gives the physical material tensor considers bi-material */
  void GetPhysicalMaterial(BiLinForm* form, const DesignElement* de, const TransferFunction* tf, bool derivative, Matrix<double>& out);

  /** Here we store the solution of the problem. Multiple solutions for multiple loadcases */
  StateContainer forward;

  /** Here we store the solution of the adjoint problem. */
  StateContainer adjoint;

  /** do we do SIMP or FreeMat or ... */
  Method method_;

  /** true, if assuming regular grid, and only optimizing density, not DesignMaterial */
  bool assume_constant_element_matrices_;

  /** cache the 1.0 / complete volume of the domain */
  double volume_fraction_;

  /** This is a helper for SetElementK() which adds for App::MECH in the harmonic case damping and mass
   * @param bimaterial describes only the material, the factor needs to be set as rho^3 or 1-rho^3 already! */
  void AddMassToStiffness(Context* ctxt, const TransferFunction* mtf, DesignElement* de, Matrix<std::complex<double> >& K_in_S_out, bool derivative, bool bimaterial, CalcMode mode = STANDARD, double ev = -1.0);
  /** The DesignStructure is required by SIMP for filters and by Condition for slope constraints
   * and checkerboard. They share this element. It can only be created by PostInit(), hence every
   * PostInit() who needs the structure needs to check if it was created before. Deleted by ~EM */
  DesignStructure* structure_;

  /** This is just a shortcut for the actual dimensions (2 or 3) */
  const unsigned int dim;

  /** Convenience class for writing the pseudo density file*/
  DensityFile* densityFile;

  /** do we perform homogenization for both permittivity and permeability? */
  bool bitensor_;

private:

  /** Checks if the G-Matrix from model reduction shall be written as special result - and does it in case */
  void OutputModRedGTensor();

  /** This is a helper for the calculation of the homogenized tensor or the derivative of it.
   * This is the inner of the sum for the homogenized tensor or the derivative formulation
   * in Bendsoe/Sigmund - Topology Optimization page 124
   * @param u1 the element solution vector
   * @return the product test strain diff * (K or K') * test strain diff  */
  static double CalcHomogenizedElementProduct(ErsatzMaterial* em, Context* ctxt, DesignElement* de, bool derivative, Vector<double>& u1,
      Vector<double>& u2, Matrix<double>& test_strain_matrix_ij, Matrix<double>& test_strain_matrix_kl);

  static Complex CalcU1KU2(ErsatzMaterial* obj, DesignElement* de, bool derivative, Vector<Complex> u1_vec, Vector<Complex> u2_vec);

  /** See the non-template version for documentation! */
  template<class T> double CalcU1KU2(TransferFunction* tf,
      StdVector<SingleVector*>& u1, App::Type k, StdVector<SingleVector*>& u2,
      DesignDependentRHS* ref, double factor, CalcMode calcMode, Function* f,
      int res_idx, double ev);

  /** for design dependent interface driven excitation f=4*rho*(1-rho) as for heat .... calculate element based f'
   * @param de for this element we compute "out"
   * @param out gets size of nodes of de elem and contains d_f/d_de */
  template<class T> void CalcInterfaceDrivenGradRHS(const DesignElement* de, Vector<T>& out);

  template<class T> void SubstractInterfaceDrivenGradRHS(const DesignElement* de, Vector<T>& in_out);

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

  /** This solves the adjoint problem problem only and stores all relevant data. Calls SetAndSolveAdjointRHS() */
  template<class T> void SolveAdjointProblem(Excitation* excite, Function* f);

  /** Set the rhs for the adjoint equation, called by assemble */
  // virtual void SetAdjointRhs(AdjointParameters* adjointParams);

  /** Set the rhs for the tracking adjoint */
  // void SetTrackingAdjointRhs(Excitation& excite, int ts);

  /** Takes care for making CFS solving the adjoint PDE. Sets the rhs as  adjoint[excite.index]->rhs[App::MECH] */
  template<class T> void SetAndSolveAdjointRHS(Excitation& excite,  Function* cost);

  /** Helper for CommitIteration. Appends or replaces a design line */
  void SetCurrentExportDesign();

  /** In ErsatzMaterial: Saves the original loads and sets the output nodes as adjoint pde rhs
   * Has distinct implementations for complex and real part.
   * @see virtual ConstructAdjointRHS() */
  void ConstructRealAdjointRHS(Excitation& excite, Function* f);
  void ConstructComplexAdjointRHS(Excitation& excite, Function* f);

  /** Calculates the Greyness OR gauss-greyness! and the derivative of the (gauss) greyness.
   * @param derivative if false the return value is calculated. Otherwise the value in
   *                   the design element is set. Optionally also grad_out
   * @param grad_out if derivative is set and grad_out is not null it is set.
   * @return invalid in derivative case*/
  double CalcGreyness(Condition* g, bool derivative);

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
};

} // namespace


#endif /*ERSATZ_MATERIAL_HH_*/
