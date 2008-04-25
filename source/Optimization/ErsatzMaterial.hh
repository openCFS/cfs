#ifndef ERSATZ_MATERIAL_HH_
#define ERSATZ_MATERIAL_HH_

#include "Optimization/Optimization.hh"
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
  virtual double CalcObjective()
  {
    return harmonic ? CalcObjective<std::complex<double> >() 
        : CalcObjective<double>();
  }

  /** with the output objective we also have to solve the adjoint problem. */
  virtual void SolveStateProblem(); 

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
  void CommitIteration();

  /** Adds validation stuff here to keep out of long constructor */
  virtual void PostInit();

  /** Types of ersatz material optimization methods, the strings are read from the xml file */
  typedef enum { SIMP_METHOD, FREE_MAT, SHAPE_GRAD, NO_METHOD } Method;
  
  static Enum<Method> method;
  
protected:

  /** This class holds the solution of the PDE. It is in a class such that it 
   * helps to encapulate real and complex solutions. Note that the Piezo 
   * has other solutions! */
  class Solution
  {
  public:
    Solution(ErsatzMaterial* em);

    ~Solution();

    typedef enum { ELEMENT_VECTORS, RAW_VECTOR } StorageType;

    /** Copies the solution for the pde in our own storage.
     * In the ELEMENT_VECTORS case make sure, that the solution is in the PDE!
     * For manual adjoint stuff you might have to do SaveSolution() first!
     * @param st if we copy the vector as RAW_VECTOR or element wise
     * @param pde will me mech in SIMP and also elec in PiezoSIMP
     * @param app redundant to pde. Either MECH or ELEC */
    void ReadSolution(StorageType st, StdPDE* pde, Application app)
    {
      if(em->harmonic) ReadSolution<std::complex<double> >(st, pde, app);
      else ReadSolution<double>(st, pde, app);
    }

    /** Writes the solution (raw vector) back to the pde */
    void WriteSolution(StdPDE* pde, Application app)
    {
      if(em->harmonic) WriteSolution<std::complex<double> >(pde, app);
      else WriteSolution<double>(pde, app);
    }

    /** This is an element wise storage of the solution 
     * the Application shall be MECH or ELEC */
    std::map<Application, StdVector<CFSVector* > > elem;

    /** This is the algsys solution vector */
    std::map<Application, CFSVector* > raw;

  private:
    template <class T>
    void ReadSolution(StorageType st, StdPDE* pde, Application app);

    template <class T>
    void WriteSolution(StdPDE* pde, Application app);

    /** Reference to our optimization problem */
    ErsatzMaterial* em;
  };

  /** When "commit" is set, we write "forward"/"adjoint" or "both_cases" */
  virtual void StoreResults(double step_val);


  /** switches to the proper constraint, also for gradient case.
   * @param design if not gradient ignored
   * @param grad_out only for gradient and even then optional if not for extern optimizer
   * @return not defined in the gradient case */
  double CalcConstraint(Condition* constraint, bool gradient, double* grad_out = NULL);

  /** Helper which extracts the Form from assemble using the optimization region
   * @param pde1 the first pde (e.g. mech)
   * @param pde2 this is either the same as pde1 or the coupling partner
   * @param integrator there is no nice enum yet :( e.g. linElastInt, MechInt, ... */ 
  BaseForm* GetForm(StdPDE* pde1, StdPDE* pde2, const std::string& integrator);

  /** <p>Get the original element matrix (stiffness, mass, ...) 
   * which is constant for all isotripic elements.
   * This method is not only for mechanical SIMP but is also used by PiezoSIMP,
   * therefore it is generic.</p>
   * <p>If no elemen is given, the one from the first design element is used.</p>
   * <p>All transfer functions are disabled during this method. Call only for
   * enabled transfer functions (default)</p> 
   * @param form to be extracted via GetForm()  
   * @param out here the element stiffness matrix written. e.h. K_uu which is \int B E B 
   * @param elem if not given the first design element is used, otherwise the provided one*/
  void GetElementMatrix(BaseForm* form, Matrix<double>& out, Elem* elem = NULL);


  /** Calculate the derivative form \f$<l,K'u-f'>\f$. 
   * 
   * When adjoint vector u1/l is not calculated with a negative rhs, one
   * can put in the minus sign as an explicit factor.
   * 
   * The K is automatically determined by the application and the size of 
   * the vector u(2).
   * 
   * @param tf the transfer function defines the design variable and multiplier of K
   * @param u1 for derivatives the lagrange vector l\f$
   * @param k the application determines the stiffness matrix
   * @param u2 the solution or u in \f$<l,K'u-f'>\f$  
   * @param rhs if one want to do \f$<l,K'u-f'>\f$ this containts the info for \f$-f'\f$.
   * @param add if true we sum up to the desing elements cost gradient
   * @param factor see above, more complex in radiation case. */
  double CalcU1KU2(TransferFunction* tf, StdVector<CFSVector*>& u1, Application k,
      StdVector<CFSVector*>& u2, SurfaceRef* rhs = NULL,
      bool add = false, double factor = -1.0)
  {
    if(harmonic) return CalcU1KU2<std::complex<double> >(tf, u1, k, u2, rhs, add, factor);
            else return CalcU1KU2<double>(tf, u1, k, u2, rhs, add, factor);
  }


  /** This is a helper for CalcU1KU2 to determine the "K" which in most cases includes a 
   * derivative. It also includes mechanical damping and mass matrix via AddMassToStiffness().
   * The templated stuff is private, as C++ does not allow virtual templates. */
  virtual void SetElementK(DesignElement* de, Application app, CFSMatrix* out) 
  {
    throw Exception("not implemented");
  }

  
  
  /** Set the RHS for the adjoint problem (templates!). For RADIATION there
   * is an implementation in SIMP */
  virtual void CalcAdjointRHS()
  {
    if(harmonic) CalcAdjointRHS<std::complex<double> >();
            else CalcAdjointRHS<double>();
  }
  
  
  /** This does a forward/ajoint step for the mechanism optimization. It works for
   * both (mechanical) SIMP and PiezoSIMP. */
  void SolveAdjointProblem()
  {
    if(harmonic)
      SolveAdjointProblem<std::complex<double> >();
    else
      SolveAdjointProblem<double>();
  }

  template <class T>
  double CalcObjective();

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
  RegionIdType       regionId;

  /** Here we store the solution of the problem */
  Solution* forward_;

  /** Here we store the soluton of the adjoint problem. */
  Solution* adjoint_;

  /** do we do SIMP or FreeMat or ... */
  Method method_;

  /** this is the optimization->simp XML element */
  ParamNode* pn;

  /** The assemble class for our PDE */
  Assemble* assemble_;

private:

  /** Creates the ErsatzMaterialFile and writes the header */
  void CreateErsatzMaterialFile(const std::string& filename, StdVector<ParamNode*>& des, StdVector<ParamNode*>& tfs);

  /** Handles the Volume constraint. Has a constraint and constraint derivative mode 
   * @param derivative if false the return value is calculated. Otherwise the value in
   *                   the design element is set. Optionally also grad_out
   * @param grad_out if derivative is set and grad_out is not null it is set.
   * @return invalid in derivative case*/ 
  double CalcVolume(bool derivative, Condition* constraint, double* grad_out);            

  /** See the non-template version for documentation! */
  template <class T>
  double CalcU1KU2(TransferFunction* tf, StdVector<CFSVector*>& u1, 
      Application k, StdVector<CFSVector*>& u2, SurfaceRef* ref, bool add, double factor);

  /** Handles sensitive RHS, e.g. when we have sensitive Neuman boundary condition (elect surface charge). 
   * Another application is  RADIATION. Then SurfaceRef is 
   * given to CalcU1KU2 and this method does from \f$<l,K'u-f'>\f$ the \f$-f'\f$ part.
   * It checks if any nodes of the design element are part of the surface and 
   * substracts for all dof of that node only */
  template <class T> 
  void SubstractGradSurfaceRHS(DesignElement* de, TransferFunction* tf, SurfaceRef* ref, Vector<T>& in_out);
    
  /** This solves the forward and adjoint problem and stores all relevant data. Calls SetAndSolveAdjointRHS() */
  template <class T>         
  void SolveAdjointProblem();

  /** Takes care for making CFS solving the adjoint PDE. Sets output_vector_ */
  template <class T>
  void SetAndSolveAdjointRHS();

  /** Stores the IDBC and sets the to HDBC for calculating the adjoint PDE */
  StdVector<IdBcList> IDBC_2_HDBC();

  /** Rests HDBC after adjoint PDE is solved */
  void ResetHDBC(StdVector<IdBcList> org_idbc);

  /** Saves the original loads and sets the output nodes as adjoint pde rhs */
  template <class T>
  LoadList SetOutputRHS();

  /** Calcs the specific rhs for the adjoint pde by objective. 
   * Specific implementations (e.g. RADIATION) in child classes */
  template <class T>
  void CalcAdjointRHS();

  /** Calculates the Greyness OR gauss-greyness! and the derivative of the (gauss) greyness.
   * @param derivative if false the return value is calculated. Otherwise the value in
   *                   the design element is set. Optionally also grad_out
   * @param grad_out if derivative is set and grad_out is not null it is set.
   * @return invalid in derivative case*/
  double CalcGreyness(bool derivative, Condition* constraint, double* grad_out); 

  /** Calculates the product of the (system) surface normal matrix with the solution already in OLAS.
   * Note that we have to use 1 based OLAS vectors as the sparse system matrix is from OLAS .
   * This calculation is done for the adjoint rhs and also for calculate the radiation objective.
   * It shall be cheap enough to calc here twice! */
  template <class T>
  void CalcSurfaceNormalTimesSolution(OLAS::Vector<T>& olas_prod);         

  /** flag indicating if we write the ersatz material file. see destructor */
  std::ofstream*     ersatzMaterialFile;  

  /** When we optimize output we store here the nodes */
  LoadList output_nodes_;

  /** This is the vector corresponding to the output nodes for <l,u>. To be deleted in destructor! */
  CFSVector* output_vector_;

};

} // namespace


#endif /*ERSATZ_MATERIAL_HH_*/
