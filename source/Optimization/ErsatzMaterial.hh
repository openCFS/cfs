#ifndef ERSATZ_MATERIAL_HH_
#define ERSATZ_MATERIAL_HH_

#include <map>

#include "Optimization/Optimization.hh"
#include "Domain/bcs.hh"
#include "MatVec/vector.hh"
#include "MatVec/matrix.hh"

namespace CoupledField
{
class LoadBc;
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
class SurfaceRef;

template <class TYPE> class StdVector;
template <class TYPE> class Vector;
template <class TYPE> class Matrix;


/** Base for optimization where the design variable is correlated to finite elements.
 * The classical case is SIMP, where one optimizes for a pseudo density. The sub-classes
 * extend this idea to more complex stuff.
 * Here the common stuff is kept as solving the adjoint pde and calculating the objective.
 * The implementation of gradients, ... is for the subclasses. */
class ErsatzMaterial : public Optimization
{
public:
  /** Up to now w/o parameters */
  ErsatzMaterial();

  /** e.g. closing the ersatzMaterialFile */
  virtual ~ErsatzMaterial();

  /** compute the value of the cost function */
  virtual double CalcObjective();


  /** Evaluates the gradient of the cost function. Saves the result to data.objective_gradient.
   * The real work is done by CalcObjectiveGradient(Excitation&) which is overloaded.
   * @param grad_out if not NULL copy write the (design space size) results.
   * @see CalcObjectiveGradient(Excitation&)
   * If filtering is enabled, this is automatically filtered. */
  void CalcObjectiveGradient(double* grad_out);


  /** This solves and stores the forward problem. Eventually the adjoint problem is called
   * implicitly.
   * Processes multiple excitations.
   * @param ev_only_excite EvaluateOnly has the special feature to perform a sweep over multiple
   *        frequencies and calculate the mono-harmonic objective. Obly EvaluateOnly shall use this
   *        parameter. */
  void SolveStateProblem(Excitation* ev_only_excite = NULL);

  /** Calculates the constraint(s) */
  double CalcConstraint(Condition* constraint = NULL)
  {
    return CalcConstraint(constraint, false); // no gradient
  }

  /** The jacobian of the gradient here as a vector with only one constraint! */
  void CalcConstraintGradient(Condition* constraint = NULL, double* grad_out = NULL)
  {
    CalcConstraint(constraint, true, grad_out);
  }

  /** Here we also write the density files */
  InfoNode* CommitIteration(bool keep_iteration_number = false);

  /** Adds validation stuff here to keep out of long constructor */
  virtual void PostInit();

  /** Types of ersatz material optimization methods, the strings are read from the xml file */
  typedef enum { SIMP_METHOD, PARAM_MAT, SHAPE_GRAD, NO_METHOD } Method;

  static Enum<Method> method;

  /** the commit mode defines what of the iterations is to be written to gid, ... */
  typedef enum { FORWARD, ADJOINT, BOTH } CommitMode;

  static Enum<CommitMode> commitMode;

  /** Here we have the set of excitations. Only relevant for the mulitple excitations
   * case (multiple loads or frequencies) */
  StdVector<Excitation> excitations;

  typedef LoadBc TrackingBc;

  typedef LoadList TrackingList;

protected:

  /** This class holds the solution of the PDE. It is in a class such that it
   * helps to encapulate real and complex solutions. Note that the Piezo
   * has other solutions! */
  class Solution
  {
  public:
    Solution(ErsatzMaterial* em);

    ~Solution();

    typedef enum { ELEMENT_VECTORS, RAW_VECTOR, RHS_VECTOR } StorageType;

    /** Copies the solution for the pde in our own storage.
     * In the ELEMENT_VECTORS case make sure, that the solution is in the PDE!
     * For manual adjoint stuff you might have to do SaveSolution() first!
     * @param st if we copy the vector as RAW_VECTOR or element wise
     * @param pde will me mech in SIMP and also elec in PiezoSIMP
     * @param app redundant to pde. Either MECH or ELEC
     * @param save_sol when an adjoint system was solved, one has to call
     *        "SaveSolution()" in the pde such that we can extract it elementwise.
     *        Only relevant for st = ELEMENT_VECTORS
     * @return NULL if st = ELEMENT_VECTOR, otherwise it is the vector */
    SingleVector* Read(StorageType st, StdPDE* pde, Application app, bool save_sol = false)
    {
      if(em_->harmonic) return Read<std::complex<double> >(st, pde, app, save_sol);
                   else return Read<double>(st, pde, app, save_sol);
    }

    /** Writes the solution (raw vector) back to the pde */
    void Write(StdPDE* pde, Application app)
    {
      if(em_->harmonic) Write<std::complex<double> >(pde, app);
      else Write<double>(pde, app);
    }

    static void Write(StdPDE* pde, SingleVector* vec);


    /** Helper method */
    SingleVector* GetVector(StorageType st, Application app);
    Vector<std::complex<double> >* GetComplexPointer(StorageType st, Application app);
    Vector<double>*                GetRealPointer(StorageType st, Application app);

    /** This is an element wise storage of the solution
     * the Application shall be MECH or ELEC */
    std::map<Application, StdVector<SingleVector* > > elem;

    /** This is the algsys solution vector */
    std::map<Application, SingleVector* > raw;

    /** This stores the right hand sides in raw format.
     * In the adjoint case this is the constructed rhs form the artifical loads.
     * Note, that this vector contains the loads but the actual rhs is * -1.0
     * and in the harmonic case also multiplied with conj(u). Check by exporting
     * the linear system.
     * In the forward case this stores the rhs from the forward excitation to perform multiple
     * loadcases. */
    std::map<Application, SingleVector* > rhs;

  private:
    template <class T>
    SingleVector* Read(StorageType st, StdPDE* pde, Application app, bool save_sol = false);

    template <class T>
    void Write(StdPDE* pde, Application app);

    /** Reference to our optimization problem */
    ErsatzMaterial* em_;
  };

  /** As the solutions come for multiple excitations in sets we store the list and the
   * averaged (when mutiple excitations are enabled). W/o is just some overhead with data size = 1 */
  struct Solutions
  {
    Solutions();

    ~Solutions();

    /** @param size number of excitations */
    void Init(int size, ErsatzMaterial* em);

    /** Contains at leas one element, more when doing multiple excitations.*/
    StdVector<Solution*> data;

    /** When we want to average this stuff */
    SingleVector* multiple;
  };

  /** When "commit" is set, we write "forward"/"adjoint" or "both_cases" */
  virtual void StoreResults(double step_val);

  /** @see Optimization::GetIterationFrequency() */
  std::string GetIterationFrequency();

  /** This implementes the actual gradient evaluations not handled by
   *  ErsatzMaterial::CalcObjectiveGradient(double*) itself.
   * @see CalcObjectiveGradient(double*) the calling method */
  virtual void CalcObjectiveGradient(Excitation& exite)
  {
    EXCEPTION("not implemented here";)
  }


  /** switches to the proper constraint, also for gradient case.
   * @param design if not gradient ignored
   * @param grad_out only for gradient and even then optional if not for extern optimizer
   * @return not defined in the gradient case */
  double CalcConstraint(Condition* constraint, bool gradient, double* grad_out = NULL);

  /** Helper which extracts the Form from assemble using the optimization region
   * @param regionId the corresponding region
   * @param pde1 the first pde (e.g. mech)
   * @param pde2 this is either the same as pde1 or the coupling partner
   * @param integrator there is no nice enum yet :( e.g. linElastInt, MechInt, ... */
  BaseForm* GetForm(RegionIdType regionId, StdPDE* pde1, StdPDE* pde2, const std::string& integrator);

  /** <p>Get the original element matrix (stiffness, mass, ...)
   * which is constant for all isotripic elements.
   * This method is not only for mechanical SIMP but is also used by PiezoSIMP,
   * therefore it is generic.</p>
   * <p>If no elemen is given, the one from the first design element is used.</p>
   * <p>All transfer functions are disabled during this method. Call only for
   * enabled transfer functions (default)</p>
   * @param form to be extracted via GetForm()
   * @param out here the element stiffness matrix written. e.h. K_uu which is \int B E B
   * @param elem if not given the first design element is used, otherwise the provided one
   * @param factor in piezoelectricity K_pp is -1* BDBInt */
  void GetElementMatrix(BaseForm* form, Matrix<double>& out, Elem* elem = NULL,
                        const DesignElement::Type direction = DesignElement::NO_DERIVATIVE, double factor = 1.0);


  /** Calculate the derivative form \f$<l,K'u-f'>\f$.
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
   * @param rhs if one want to do \f$<l,K'u-f'>\f$ this containts the info for \f$-f'\f$.
   * @param add if true we sum up to the desing elements cost gradient
   * @param factor see above, more complex in radiation case. */
  double CalcU1KU2(TransferFunction* tf, StdVector<SingleVector*>& u1, Application k,
      StdVector<SingleVector*>& u2, SurfaceRef* rhs = NULL,
      bool add = false, double factor = 1.0, Condition* constraint = NULL)
  {
    if(harmonic) return CalcU1KU2<std::complex<double> >(tf, u1, k, u2, rhs, add, factor, constraint);
            else return CalcU1KU2<double>(tf, u1, k, u2, rhs, add, factor, constraint);
  }


  /** This is a helper for CalcU1KU2 to determine the "K" which in most cases includes a
   * derivative. It also includes mechanical damping and mass matrix via AddMassToStiffness().
   * The templated stuff is private, as C++ does not allow virtual templates. */
  virtual void SetElementK(DesignElement* de, Application app, DenseMatrix* out)
  {
    throw Exception("not implemented");
  }
  /** Get the ElementStiffness Matrix for this element, this is the region constant version
   * @param elem the Element for which the Matrix should be returned
   * @param direction if given, calculate derivative of Stiffness Matrix instead
   * @return a pointer to the Element Stiffness Matrix*/
  const Matrix<double>& MechStiffness(Elem* elem, const DesignElement::Type direction = DesignElement::NO_DERIVATIVE);

  /** Get the ElementMass Matrix for this element, this is the region constant version
   * @param elem the Element for which the Matrix should be returned
   * @return a pointer to the Element Mass Matrix*/
  const Matrix<double>& MechMass(Elem* elem);

  /** Get the ErsatzMaterialTensor as the Tensor itself, not the stiffness matrix
   * @param mat holds the tensor
   * @param elem the Element for which the tensor should be returned
   * @param direction if given return derivative in that direction*/
  void GetErsatzMaterialTensor(Matrix<double>& mat, Elem* elem, DesignElement::Type direction = DesignElement::NO_DERIVATIVE);


  /** This sets the objective specific RHS for adjoint problems.
   * For DYNAMIC_OUTPUT this is a post processing of of ConstructOutputRHS.
   * For RADIATION there is a specific implementation in SIMP. */
  virtual void AdjustComplexAdjointRHS(Excitation& excite);

  /** This is an extension to SolveStateProblem() where the forward problem is solved and stored.
   * Depending on the objective function SolveAdjointProblem() is called to additionally solve and store the
   * adjoint problem.
   * It works for both (mechanical) SIMP and PiezoSIMP. */
  void SolveAdjointProblem(Excitation& excite)  {
    if(harmonic) SolveAdjointProblem<std::complex<double> >(excite);
            else SolveAdjointProblem<double>(excite);
  }

  /** This is the tempated instance of the virtual CalcObjective().
   * To perform multiple excitation it calls CalcObjective(Excitation& excite) and
   * performs the averaging. */
  template <class T>
  double CalcObjective();

  double CalcObjective(Excitation& excite) {
    if(harmonic) return CalcObjective<std::complex<double> >(excite);
            else return CalcObjective<double>(excite);
  }

  /** This calculates the objective for the given excitation. The result is also stored
   * in excite.cost. It includes the objective factor (e.g. omega^2) but not the weighting */
  template <class T>
  double CalcObjective(Excitation& excite);


  /** Handles the Volume constraint. Has a constraint and constraint derivative mode
   * @param derivative if false the return value is calculated. Otherwise the value in
   *                   the design element is set.
   * @param constraint if set, calculate as given constraint, if null calculate as objective
   * @param grad_out if derivative is set and grad_out is not null it is set.
   * @return invalid in derivative case*/
  double CalcVolume(bool derivative, Condition* constraint);

  /** Handles the Compliance constraint/objective. Has a objective, objective derivative,
   * constraint and constraint derivative mode
   * @param excite for multiple excitations, defines which solution and rhs to use
   *               (default the first any only one)
   * @param derivative if false the return value is calculated. Otherwise the value in
   *                   the design element is set.
   * @param constraint if set calculate as given constraint, if null calculate as objective
   * @param grad_out if derivative is set and grad_out is not null it is set.
   * @return invalid in derivative case*/
  double CalcCompliance(Excitation& excite, bool derivative, Condition* constraint);

  /** Calculates the objective only, no derivative */
  double CalcGlobalDynamicCompliance(Excitation& excite);

  /** Calculates <l,u> or <conj(u) L, u> where l/L is adjoint[idx]->rhs */
  template <class T>
  double CalcOutputObjective(Excitation& excite);

  /** Handles the Tracking constraint/objective. Has a objective, objective derivative,
   * constraint and constraint derivative mode
   * @param derivative if false the return value is calculated. Otherwise the value in
   *                   the design element is set.
   * @param constraint if set calculate as given constraint, if null calculate as objective
   * @param grad_out if derivative is set and grad_out is not null it is set.
   * @return invalid in derivative case*/
  double CalcTracking(bool derivative, Condition* constraint);

  /** Helper that converts from mechPDE to MECH and elecPDE to ELEC
   * @throws if neither mechPDE nor elecPDE */
  Application ToApp(SinglePDE* pde) const;

  /** Find our pde in SIMP by application */
  SinglePDE* ToPDE(Application app, bool throw_exception = true) const;


  /** We have always a MechPDE for SIMP. This (and elec) is also in pde! */
  MechPDE*           mech;

  /** This vector is an alternative access to the mech and elec pointers.
   * mech is guaranteed to be index 0 and elec would be 1.
   * @see ToApp()
   * @see ToPDE() */
  StdVector<SinglePDE*> pde;

  /** The region to optimize */
  StdVector<RegionIdType> regionIds;

  /** Here we store the solution of the problem. Multiple solutions for multiple loadcases */
  Solutions forward;

  /** Here we store the soluton of the adjoint problem. */
  Solutions adjoint;

  /** do we do SIMP or FreeMat or ... */
  Method method_;

  /** this is the optimization->simp XML element */
  ParamNode* pn;

  /** The assemble class for our PDE */
  Assemble* assemble_;

  /** The mechanical element stiffness matrix is constant */
  std::map<RegionIdType, Matrix<double> > mechStiffness_map;

  /** The mechanical element mass matrix is also constant. Only for harmonic! */
  std::map<RegionIdType, Matrix<double> > mechMass_map;

  /** Here we store the solution of the tracking subproblem. */
  Solution* tracking_;

  /**  Here we store all tracking nodelists. */
  TrackingList tracks_;

  /** whether we can assume a regular grid or not
   * note, that only volume, compliance and tracking objective/constraints are working with non-regular grids */
  bool assume_regular_grid_;

  /** true, if assuming regular grid, and only optimizing density, not DesignMaterial */
  bool assume_constant_element_matrices_;

  /** cache the 1.0 / complete volume of the domain */
  double volume_fraction_;

private:

  /** Creates the ErsatzMaterialFile and writes the header */
  void CreateErsatzMaterialFile(const std::string& filename, StdVector<ParamNode*>& des, StdVector<ParamNode*>& tfs);

  /** See the non-template version for documentation! */
  template <class T>
  double CalcU1KU2(TransferFunction* tf, StdVector<SingleVector*>& u1,
      Application k, StdVector<SingleVector*>& u2, SurfaceRef* ref, bool add, double factor, Condition* constraint);

  /** Handles sensitive RHS, e.g. when we have sensitive Neuman boundary condition (elect surface charge).
   * Another application is  RADIATION. Then SurfaceRef is
   * given to CalcU1KU2 and this method does from \f$<l,K'u-f'>\f$ the \f$-f'\f$ part.
   * It checks if any nodes of the design element are part of the surface and
   * substracts for all dof of that node only */
  template <class T>
  void SubstractGradSurfaceRHS(DesignElement* de, TransferFunction* tf, SurfaceRef* ref, Vector<T>& in_out);

  /** This solves the adjoint problem problem only and stores all relevant data. Calls SetAndSolveAdjointRHS() */
  template <class T>
  void SolveAdjointProblem(Excitation& excite);

  /** Store the results from the forward/adjoint problem. Handles multiple excitations
   * @param read_rhs is only interesting for the forward problem
   * @param save_sol set this in the adjoint problem -> see Solution::Read()
   * @param comment is just to LOG_DBG */
  void StorePDESolution(Excitation &excite, Solutions& solutions, bool read_rhs, bool save_sol, const std::string& comment);

  /** Takes care for making CFS solving the adjoint PDE. Sets the rhs as  adjoint[excite.index]->rhs[MECH] */
  template <class T>
  void SetAndSolveAdjointRHS(Excitation& excite);

  /** Stores the IDBC and sets the to HDBC for calculating the adjoint PDE
   * @return this are the original IDBC values for RestHDBC*/
  StdVector<IdBcList> SetHDBC();

  /** Rests HDBC after adjoint PDE is solved */
  void ResetHDBC(StdVector<IdBcList> org_idbc);

  /** Saves the original loads and sets the output nodes as adjoint pde rhs */
  template <class T>
  LoadList ConstructOutputRHS(Excitation& excite);

  /** Calculates the Greyness OR gauss-greyness! and the derivative of the (gauss) greyness.
   * @param derivative if false the return value is calculated. Otherwise the value in
   *                   the design element is set. Optionally also grad_out
   * @param grad_out if derivative is set and grad_out is not null it is set.
   * @return invalid in derivative case*/
  double CalcGreyness(bool derivative, Condition* constraint);

  /** Calculates the product of the (system) surface normal matrix with the solution already in OLAS.
   * Note that we have to use 1 based OLAS vectors as the sparse system matrix is from OLAS .
   * This calculation is done for the adjoint rhs and also for calculate the radiation objective.
   * It shall be cheap enough to calc here twice! */
  template <class T>
  void CalcSurfaceNormalTimesSolution(Vector<T>& olas_prod);

  /** Handle multiple excitations (loads/frquencies). By defefinition the size is almost 1, even
   * if there is no load (e.g. static piezo with inhomgeneous Dirichlet BC. */
  void PrepareMultipleExcitations();

  /** For doing adjust weights when doing multiple excitation with meta objective, this method
   * does the job. It requires the cost entries in excitations to be set. 
   * The \f$w_k^p=const\;\sum w_k = 1\f$ condition is fulfilled here. */ 
  void NormalizeMultipleExcitations();

  CommitMode commitMode_;

  /** flag indicating if we write the ersatz material file. see destructor */
  std::ofstream*     ersatzMaterialFile;

  /** When we optimize output we store here the nodes */
  LoadList output_nodes_;
};

} // namespace


#endif /*ERSATZ_MATERIAL_HH_*/
