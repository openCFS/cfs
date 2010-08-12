#ifndef ERSATZ_MATERIAL_HH_
#define ERSATZ_MATERIAL_HH_

#include <map>

#include "Optimization/Optimization.hh"
#include "Domain/bcs.hh"
#include "Utils/result.hh"
#include "MatVec/vector.hh"
#include "MatVec/matrix.hh"

namespace CoupledField
{
class StdPDE;
class SinglePDE;
class BaseForm;
class BiLinFormContext;
class BaseMaterial;
class Condition;
class Assemble;
class TransferFunction;
class SurfElem;
class SurfaceRef;
class OptimizationMaterial;
class DesignStructure;
class DensityFile;

template<class TYPE> class StdVector;
template<class TYPE> class Vector;
template<class TYPE> class Matrix;

/** Base for optimization where the design variable is correlated to finite elements.
 * The classical case is SIMP, where one optimizes for a pseudo density. The sub-classes
 * extend this idea to more complex stuff.
 * Here the common stuff is kept as solving the adjoint pde and calculating the objective.
 * The implementation of gradients, ... is for the subclasses. */
class ErsatzMaterial: public Optimization
{

protected:
  // forward declaration
  class Solutions;

public:

  /** Up to now w/o parameters */
  ErsatzMaterial();

  /** e.g. closing the exportDesign */
  virtual ~ErsatzMaterial();

  /** Compute the value of the cost functions.
   *  Do Overwrite the excitation version, not this one!
   * To perform multiple excitation it calls CalcObjective(Objective*, Excitation&) and
   * performs the averaging. Handles objectives with single evaluation */
  double CalcObjective();

  /** Evaluates the gradient of the cost function. Saves the result to data.objective_gradient.
   * The real work is done by CalcObjectiveGradient(Excitation&) which is overloaded.
   * @param grad_out if not NULL copy write the (design space size) results.
   * @see CalcObjectiveGradient(Excitation&)
   * If filtering is enabled, this is automatically filtered. */
  void CalcObjectiveGradient(StdVector<double>* grad_out);

  /** Calculates the constraint(s) for all excitations */
  double CalcConstraint(Condition* constraint = NULL);

  /** The jacobian of the gradient here as a vector with only one constraint! */
  void CalcConstraintGradient(Condition* constraint = NULL, StdVector<double>* grad_out = NULL);

  /** This solves and stores the forward problem. Eventually the adjoint problem is called
   * implicitly.
   * Processes multiple excitations.
   * @param ev_only_excite EvaluateOnly has the special feature to perform a sweep over multiple
   *        frequencies and calculate the mono-harmonic objective. Only EvaluateOnly shall use this
   *        parameter. */
  void SolveStateProblem(Excitation* ev_only_excite = NULL);

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

  /** Helper to convert from natural solution/design to application
   * @param DesignElement::DENSITY -> MECH, DesignElement::POLARIZATION -> ELEC */
  static Application ToApp(DesignElement::Type dt);

  /** Helper that converts from mechPDE to MECH and elecPDE to ELEC, ...
   * @throws if neither mechPDE nor elecPDE
   *  @see ToPDE()
   * @see SetPDEs() */
  Application ToApp(SinglePDE* pde) const;

  /** Find our PDE in SIMP by application from the pdes map
   * @see ToApp()*/
  SinglePDE* ToPDE(Application app, bool throw_exception = true) const;

  /** Helper which extracts the FormContext from assemble using the optimization region
   * @param regionId the corresponding region
   * @param pde1 the first pde (e.g. mech)
   * @param pde2 this is either the same as pde1 or the coupling partner
   * @param integrator there is no nice enum yet :( e.g. linElastInt, MechInt, ... */
  BiLinFormContext* GetFormContext(RegionIdType regionId, StdPDE* pde1,
      StdPDE* pde2, const std::string& integrator);

  /** Helper which extracts the Form from assemble using the optimization region
   * @param regionId the corresponding region
   * @param pde1 the first pde (e.g. mech)
   * @param pde2 this is either the same as pde1 or the coupling partner
   * @param integrator there is no nice enum yet :( e.g. linElastInt, MechInt, ... */
  BaseForm* GetForm(RegionIdType regionId, StdPDE* pde1, StdPDE* pde2,
      const std::string& integrator);

  /** Types of ersatz material optimization methods, the strings are read from the xml file */
  typedef enum
  {
    NO_METHOD, SIMP_METHOD, PARAM_MAT, SHAPE_GRAD, SHAPE_OPT, SHAPE_PARAM_MAT
  } Method;

  static Enum<Method> method;

  /** The order of the pdes is not defined, Therefore we use the map
   * @see ToApp()
   * @see ToPDE() */
  std::map<Application, SinglePDE*> pdes;

  /** This is simple one SinglePDE from pdes. */
  SinglePDE* pde;

  /** Here we have the set of excitations. Only relevant for the multiple excitations
   * case (multiple loads or frequencies) */
  StdVector<Excitation> excitations;

  /** The region to optimize */
  StdVector<RegionIdType> regionIds;

  /** This is the current homogenized tensor.
   * Evaluated by HOMOGENIZATION_TRACKING and HOMOGENIZED_TENSOR (as objective only).
   * MechPDE reads it when "homogenizedTensor" is a region result! */
  Matrix<double> homogenizedTensor;

 protected:

  /** This class holds the solution of the PDE. It is in a class such that it
   * helps to encapsulate real and complex solutions. Note that the Piezo
   * has other solutions! */
  class Solution
  {
  public:

    Solution(ErsatzMaterial* em);

    ~Solution();

    typedef enum
    {
      ELEMENT_VECTORS, RAW_VECTOR, RHS_VECTOR, SEL_VECTOR, GRIDELEM_VECTORS
    } StorageType;

    /** Copies the solution for the pde in our own storage.
     * In the ELEMENT_VECTORS case make sure, that the solution is in the PDE!
     * For manual adjoint stuff you might have to do SaveSolution() first!
     * @param st if we copy the vector as RAW_VECTOR or element wise
     * @param pde will me mech in SIMP and also elec in PiezoSIMP
     * @param app redundant to pde. NO_APP for non ELEMENT_VECTOR as long as OLAS makes no difference
     * @param save_sol when an adjoint system was solved, one has to call
     *        "SaveSolution()" in the pde such that we can extract it element wise.
     *        Only relevant for st = ELEMENT_VECTORS
     * @return NULL if st = ELEMENT_VECTOR, otherwise it is the vector */
    SingleVector* Read(StorageType st, StdPDE* pde, Application app = NO_APP, bool save_sol = false);

    /** Writes the solution (raw vector) back to the pde */
    void Write(StdPDE* pde);

    static void Write(StdPDE* pde, SingleVector* vec);

    /** average the raw solutions by the excitations.
     * @param excitations average or solution by one entry only. Is strictly speaking already known
     * by the em_ parameter but is more explicit this way.  */
    static void Write(StdPDE* pde, ErsatzMaterial::Solutions& sol, Function* f, int time_step, StdVector<Excitation>& excitations);

    /** return an existing nodal vector.
     * As the type is not known we cannot create on the fly.
     * @param st RHS_VECTOR, RAW_VECTOR, SEL_VECTOR
     * asser() if vector exists (debug mode, NULL in release) */
    SingleVector* GetVector(StorageType st);

    /** return (eventually create) a nodal vector.
     *  creates what is desired when the vector does not exist yet and hence always return a vector.
     * @see GetVector() */
    Vector<double>& GetRealVector(StorageType st);

    /** @see GetRealVector() */
    Vector<std::complex<double> >& GetComplexVector(StorageType st);

    /** This is an element wise storage of the solution
     * the Application shall be MECH or ELEC */
    std::map<Application, StdVector<SingleVector*> > elem;

    /** This is an element wise storage of the solution
     * considering all elements from the grid instead of all design elements only
     * needed by shape optimization */
    std::map<Application, StdVector<SingleVector*> > gridelem;

  private:
    template<class T>
    SingleVector* Read(StorageType st, StdPDE* pde, Application app,
        bool save_sol = false);

    template<class T>
    void Write(StdPDE* pde);

    template<class T>
    static void Write(StdPDE* pde, ErsatzMaterial::Solutions& sol, Function* f, int time_step, StdVector<Excitation>& excitations);

    /** common helper for the Get*Vector() stuff */
    template<class T>
    SingleVector* GetVector(StorageType st, bool create);

    /** This is the algsys solution vector. */
    SingleVector* raw;
    //std::map<Application, SingleVector* > raw;

    /** This stores the right hand sides in raw format.
     * In the adjoint case this might be based on select for output like objectives.
     * In the forward case this stores the rhs from the forward excitation to perform multiple
     * load cases. */
    SingleVector* rhs;
    // std::map<Application, SingleVector* > rhs;

    /** For output like objectives this stores the selection l resp. L which is found by
     * an artificial load-case in the adjoint section. Otherwise it is not used.
     * This kind of base for the rhs information which is additionally multiplied by -1, -1 u^*, ...*/
    SingleVector* select;
    //std::map<Application, SingleVector* > select;

    /** Reference to our optimization problem */
    ErsatzMaterial* em_;
  };

  /** As the solutions come for multiple excitations in sets we store the list and the
   * averaged (when mutiple excitations are enabled). W/o is just some overhead with data size = 1 */
  class Solutions
  {
  public:
    friend class Solution;

    Solutions();

    ~Solutions();

    /** init when em is known */
    void Init(ErsatzMaterial* em);

    /** The solution is identified by Function, excitation index (0-based) and time step.
     * @param f the function is NULL for the forward problems, for the adjoints it needs to be given! */
    Solution* Get(Excitation& excitation, Function* f = NULL, unsigned int timestep = 0);

    Solution* Get(int excitation_index,  Function* f = NULL, unsigned int timestep = 0);

    /** Returns the currently stored functions. Empty for forward */
    StdVector<Function*> GetFunctions() const;

  private:

    /** On the fly init when the function has not been used before */
    void Init(Function* f);

    /** Contain the excitations and summarized multiple data for one problem.
     * For almost all cases there is only one problem. */
    struct Unit
    {
      Unit() { }; // only for compiler

      Unit(ErsatzMaterial* em);

      ~Unit();

      /* Contains at least one element, more when doing multiple excitations.*/
      StdVector<Solution*> data;
    };

    /** Stores the excitation by function and by time stamp.
     * Stored are units which contains eventually multiple excitations.
     * @see Unit() */
    std::map<Function*, StdVector<Unit*> > data_;

    ErsatzMaterial* em_;
  };


  /** When "commit" is set, we write "forward"/"adjoint" or "both_cases" */
  virtual void StoreResults(double step_val);

  /** @see Optimization::GetIterationFrequency() */
  std::string GetIterationFrequency();

  /** This are the modes for CalcU1KU2(). */
  enum CalcMode
  {
    STANDARD = 0, /*!< add u1^T (K' u2  - f') or2 * Re{ u1^T (K' u2 - f')} in the harmonic case  */
    CONJ_QUAD
  /*!< add <u, K' u> which is in the real case as STANDARD
   and for the harmonic case u^T K' u^* (conj. complex). u1 = u2 = u!! */
  };

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
   * @param f if set then an objective gradient with the set index is stored
   * @param g as f for constrained gradient
   * @param res_idx store in de->specialResult. use ErsatzMaterial::GetSpecialResultIndex() -1 is no special result*/
  double CalcU1KU2(TransferFunction* tf, StdVector<SingleVector*>& u1, Application k, StdVector<SingleVector*>& u2,
                     SurfaceRef* rhs, double factor, CalcMode calcMode, Function* f, int res_idx = -1);

  /** Helper calling CalcU1KU2()
   * If there is a result with value='costGradient' or 'constraintGradient' it is checked for detail='mech_mech',
   * 'elec_elec', 'elec_elec_quad', 'elec_mech', 'mech_elec' */
  int GetSpecialResultIndex(Application app1, Application app2,
      CalcMode calcMode = STANDARD, Condition* constraint = NULL);

  /** This is a helper for CalcU1KU2 to determine the "K" which in most cases includes a
   * derivative. It also includes mechanical damping and mass matrix via AddMassToStiffness().
   * The templated stuff is private, as C++ does not allow virtual templates. */
  virtual void SetElementK(DesignElement* de, Application app,
      DenseMatrix* out, CalcMode calcMode, bool derivative = true) { throw Exception("not implemented"); }

  /** Get the ErsatzMaterialTensor as the Tensor itself, not the stiffness matrix
   * @param mat holds the tensor
   * @param elem the Element for which the tensor should be returned
   * @param direction if given return derivative in that direction*/
  void GetErsatzMaterialTensor(Matrix<double>& mat, Elem* elem,
      DesignElement::Type direction = DesignElement::NO_DERIVATIVE);

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
  double CalcGlobalDynamicCompliance(Excitation& excite, Objective* f);

  /** Calculates <l,u> or <conj(u) L, u> where l/L is adjoint[idx]->rhs */
  template<class T>
  double CalcOutputObjective(Excitation& excite, Objective* f);

  /** Handles the Tracking constraint/objective. Has a objective, objective derivative, 
   * constraint and constraint derivative mode
   * @param excite  the excitation used 
   * @param derivative if false the return value is calculated. Otherwise the value in
   *                   the design element is set.
   * @param constraint if set calculate as given constraint, if null calculate as objective
   * @param solveproblem solve the tracking problem, e.g. shapeopt does solve the same problem already
   * @return invalid in derivative case*/
  virtual double CalcTracking(Excitation& excite, Objective* f, Condition* g,
      bool derivative);

  /** Calculate the energy flux through a surface region: 1/2*Re{j*u^T Q u^*} where
   * Q is the grad operator in z direction. Only for acoustic but easy to extend!*/
  double CalcEnergyFlux(Excitation& excite, Objective* f);

  /** The stress constraints, following Kucvara and Stingl; 2007, are implemented as global local function
   * (like globalSlope, globalOscillation) such that we can share the globalization techniques. The
   * element values are calculated here by using MechPDE::CalcVonMisesStress() and Function::Local::Idenitifies::CalcStress()
   * will decide if the value is used or set to zero or whatever. Therefore Function::Local::Idenitifies::EvalFunction()
   * gets the result of this method. It is a template to make it easier to work with MechPDE only */
  template <class TYPE>
  StdVector<double> CalcStress(Excitation& excite, Function* f);

  /** This is a helper with the common part for CalcEnergyFlux and the adjoint RHS.
   * Determines the global vector Q*u^* or (Q - Q^T)^T*u^* in the adjoint case.
   * @param f the cost function as we need the ParamNode
   * @param u_glob the complete equation index global solution vector
   * @param adjoint switch action
   * @param q_u_glob output depending on adjoint flag. see u_glob */
  void SetEnergyFluxVector(Function* f,
      const Vector<std::complex<double> >& u_glob, bool adjoint, Vector<
          std::complex<double> >& q_u_glob);

  /** Find the node numbers which are common from a surface element and a volume element.
   * This maps from a surface element to the volume element.
   * @param common_nodes e.g. a center surface element node is lost on a 20-hex-volume element */
  void FindCommonNodes(const SurfElem* se, const Elem* vol, StdVector<
      unsigned int>& common_nodes) const;

  /** does the substep of solving K z = Proj(u - u0) for z */
  void SolveTrackingProblem(Excitation& excite, bool designelem = true,
      bool gridelem = false);

  /** converts the teststrain vector in voigt notation to the corresponding matrix
   * @param matrix output
   * @param vec input */
  void SetTestStrainMatrix(Matrix<double> &matrix, const Vector<double> &vec);

  /** takes the result of the test strain computations and calculates the homogenized 
   *  material tensor (see Bendsoe/Sigmund: Topology Optimization, p. 122ff.
   *  It must be called only for the last excitation when all test strains are known.
   *  Writes the tensor to info.log */
  Matrix<double> CalcHomogenizedTensor();

  /** Calculates the gradient of the homogenization tracking. When J = 0.5 * || E^* - E^H ||^2
   * the this calulates -1 (E^* - E^H) * d(E^H)/d(rho_e) using a matrix scalar product
   * @param target E^* what we want
   * @param hom the pre calculated tensor E^H
   * @param g the HOMOGENIZATION_TRACKING or NULL if for objective function */
  void CalcHomogenizedTrackingGradient(const Matrix<double>& target,
      const Matrix<double>& hom, Objective* f, Condition* g);

  /** Calculates the gradient if the constraints E^H = E^* where for each interested
   * tensor entry a own HOMOGENIZATION_TENSOR constraint is required.
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
   * @return the E^H tensor entry if !derivative or 0 */
  double CalcHomogenizedTensorEntry(const std::pair<int, int> entry,
      bool derivative, StdVector<double>& grad_out);

  /** This is to be overwritten for any case there are other PDEs in ErsatzMaterial::pdes to be set.
   * PiezoSIMP does it simply in the constructor */
  virtual void SetPDEs();

  /** Here we store the solution of the problem. Multiple solutions for multiple loadcases */
  Solutions forward;

  /** Here we store the solution of the adjoint problem. */
  Solutions adjoint;

  /** do we do SIMP or FreeMat or ... */
  Method method_;

  /** this is the optimization->simp XML element */
  PtrParamNode pn;

  /** The assemble class for our PDE */
  Assemble* assemble_;

  /** true, if assuming regular grid, and only optimizing density, not DesignMaterial */
  bool assume_constant_element_matrices_;

  /** cache the 1.0 / complete volume of the domain */
  double volume_fraction_;

  /** This contains our concrete material class */
  OptimizationMaterial* material;

protected:
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
  virtual void StorePDESolution(Solutions& solutions, Excitation &excite, Function* f,
      unsigned int timestep, bool read_sol, bool read_rhs, bool save_sol, const std::string& comment);

  virtual void TimeStepCalculated(UInt timeStep, AdjointParameters* adjParams);

  virtual void RhsCalculated(AdjointParameters* adjParams);

  /** The DesignStructure is required by SIMP for filters and by Condition for slope constraints
   * and checkerboard. They share this element. It can only be created by PostInit(), hence every
   * PostInit() who needs the structure needs to check if it was created before. Deleted by ~EM */
  DesignStructure* structure_;

  /** This is just a shortcut for the actual dimensions (2 or 3) */
  const unsigned int dim;

  /** Convenience class for writing the pseudo density file*/
  DensityFile* densityFile;

private:

  /** This is a helper for the calculation of the homogenized tensor or the derivative of it.
   * This is the inner of the sum for the homogenized tensor or the derivative formulation
   * in Bendsoe/Sigmund - Topology Optimization page 124
   * @param u1 the element solution vector
   * @return the product test strain diff * (K or K') * test strain diff
   */
  static double CalcHomogenizedElementProduct(ErsatzMaterial* em,
      DesignElement* de, bool derivative, Vector<double>& u1,
      Vector<double>& u2, Matrix<double>& test_strain_matrix_ij,
      Matrix<double>& test_strain_matrix_kl);

  /** See the non-template version for documentation! */
  template<class T>
  double
      CalcU1KU2(TransferFunction* tf, StdVector<SingleVector*>& u1,
          Application k, StdVector<SingleVector*>& u2, SurfaceRef* ref,
          double factor, CalcMode calcMode, Function* f, int res_idx);

  /** Calculates a scalar product of two vectors and the derivative of the right hand side newmark update,
   * used for transient optimization derivative calculation
   * @param excite excitation to consider
   * @param forward first vector
   * @param adjoint second vector
   * @param factor factor to multiply the value by (can be excitation weight)
   * @param f objective the result is to be stored with
   * @param g constraint the result is to be stored with */
  void CalcNewmarkDerivative(Excitation& excite, Solutions& forward,
      Solutions& adjoint, double factor, Objective* f, Condition* g);

  /** Handles sensitive RHS, e.g. when we have sensitive Neuman boundary condition (elect surface charge).
   * SurfaceRef is  given to CalcU1KU2 and this method does from \f$<l,K'u-f'>\f$ the \f$-f'\f$ part.
   * It checks if any nodes of the design element are part of the surface and
   * substracts for all dof of that node only */
  template<class T>
  void SubstractGradSurfaceRHS(DesignElement* de, TransferFunction* tf,
      SurfaceRef* ref, Vector<T>& in_out);

  /** This solves the adjoint problem problem only and stores all relevant data. Calls SetAndSolveAdjointRHS() */
  template<class T>
  void SolveAdjointProblem(Excitation* excite, Function* f);

  /** Set the rhs for the adjoint equation, called by assemble */
  virtual void SetAdjointRhs(AdjointParameters* adjointParams);

  /** Set the rhs for the tracking adjoint */
  void SetTrackingAdjointRhs(Excitation& excite, int ts);

  /** Takes care for making CFS solving the adjoint PDE. Sets the rhs as  adjoint[excite.index]->rhs[MECH] */
  template<class T>
  void SetAndSolveAdjointRHS(Excitation& excite, Function* cost);

  /** Helper for CommitIteration. Appends or replaces a design line */
  void SetCurrentExportDesign();

  /** Stores the IDBC and sets the to HDBC for calculating the adjoint PDE. Also applies
   * @return this are the original IDBC values for RestHDBC*/
  StdVector<pair<SinglePDE*, IdBcList> > SetHDBC();

  /** Rests HDBC after adjoint PDE is solved */
  void ResetHDBC(StdVector<pair<SinglePDE*, IdBcList> >& org_idbc);

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

  /** Calculates globalized local functions. globalSlope and globalCheckerboard.
   * When g_i is the slope function x_i - x_i+1 -c and g_i+1 = x_1+1 - x_i - c
   * the global slope is sum max(0, g_i)^2, hence we need NEXT_AND_REVERSE locality
   * @param stress only for stress to be set */
  double CalcGlobalFunction(Function* f, bool derivative, const StdVector<double>* stress = NULL);

  /** IntegrateDesignVariables() can do a lot, but no one wants to extend it to hande the derivative
   * case of the gap constraint: volume - penalized volume */
  void CalcRegularGapConstraint(Objective* f, Condition* g,   DesignElement::Type dt);

  /** Homogenization objective/ constraint.
   * Is once evaluate only! */
  double CalcPoissonsRatioAndYoungsModulus(Objective* cost, Condition* g,  bool derivative);

  /** Calculates the product of the (system) surface normal matrix with the solution already in OLAS.
   * Note that we have to use 1 based OLAS vectors as the sparse system matrix is from OLAS .
   * This calculation is done for the adjoint rhs and also for calculate the radiation objective.
   * It shall be cheap enough to calc here twice! */
  template<class T>
  void CalcSurfaceNormalTimesSolution(Vector<T>& olas_prod);

  /** Handle multiple excitations (loads/frquencies). By defefinition the size is almost 1, even
   * if there is no load (e.g. static piezo with inhomgeneous Dirichlet BC. */
  void PrepareMultipleExcitations();

  /** Helper for PrepareMultipleExcitations(). Excitations are set with hard coded test strains */
  int SetHomogenizationTestStrains();

  /** Helper for PrepareMultipleExcitations(). Excitations are set with hard coded polarization matrix excitations */
  int SetPolarizationMatrixExcitations();

  /** For doing adjust weights when doing multiple excitation with meta objective, this method
   * does the job. It requires the cost entries in excitations to be set. 
   * The \f$w_k^p=const\;\sum w_k = 1\f$ condition is fulfilled here. */
  void NormalizeMultipleExcitations();



  /** When we optimize output we store here the nodes */
  LoadList output_nodes_;

  /** do we perform homogenization induced by any of the objective or constraints? */
  bool homogenization_;
};

} // namespace


#endif /*ERSATZ_MATERIAL_HH_*/
