#ifndef STATE_SOLUTION_HH_
#define STATE_SOLUTION_HH_

#include <stddef.h>
#include <complex>
#include <map>
#include <string>
// #include <utility>

#include "General/Enum.hh"
#include "General/defs.hh"
#include "General/Environment.hh"
//#include "General/Exception.hh"
#include "MatVec/Vector.hh"
#include "Optimization/Optimization.hh"
#include "Utils/StdVector.hh"

namespace CoupledField
{
class SinglePDE;
class StateSolution;
class SingleVector;

/** As the solutions come for multiple excitations in sets we store the list and the
 * averaged (when multiple excitations are enabled). W/o is just some overhead with data size = 1 */
class StateSolutions
{
public:
  friend class StateSolution;

  StateSolutions();

  ~StateSolutions();

  /** init when em is known */
  void Init(ErsatzMaterial* em);

  /** The solution is identified by Function, excitation index (0-based) and time step.
   * @param f the function is NULL for the forward problems, for the adjoints it needs to be given! */
  StateSolution* Get(Excitation& excitation, Function* f = NULL, unsigned int timestep = 0, const DERIVType derivative = NO_DERIVTYPE);

  StateSolution* Get(int excitation_index,  Function* f = NULL, unsigned int timestep = 0, const DERIVType derivative = NO_DERIVTYPE);

  /** Returns the currently stored functions. Empty for forward */
  StdVector<Function*> GetFunctions() const;

  /** Return whether this is the Solution of the forward problem */
  bool IsForward(){ return(isForward); };

  /** Set whether this Solution is the Solution of the forward problem */
  void SetIsForward(bool forward){ isForward = forward; };

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
    StdVector<StateSolution*> data;
  };

  /** Stores the excitation by function and by derivative and time step.
   * data_[function][derivative][timestep]
   * Stored are units which contains eventually multiple excitations.
   * @see Unit() */
  std::map<Function*, StdVector<StdVector<Unit*> > > data_;

  // if this Solutions is forward, it does not use the value in function in Get
  bool isForward;

  // Pointer to data[NULL] to speed up things
  StdVector<StdVector<Unit*> >* forward_data_;

  ErsatzMaterial* em_;
};

/** This class holds the solution of the PDE. It is in a class such that it
 * helps to encapsulate real and complex solutions. Note that the Piezo
 * has other solutions! */
class StateSolution
{
public:

  StateSolution(ErsatzMaterial* em);

  ~StateSolution();

  typedef enum
  {
    ELEMENT_VECTORS = 0, RAW_VECTOR, RHS_VECTOR, SEL_VECTOR, GRIDELEM_VECTORS
  } StorageType;

  static SolutionType GetSolutionType(SinglePDE* pde, Optimization::Application app = Optimization::NO_APP);

  /** Copies the solution for the pde in our own storage.
   * In the ELEMENT_VECTORS case make sure, that the solution is in the PDE!
   * For manual adjoint stuff you might have to do SaveSolution() first!
   * In ELEMENT_VECTORS also the registered pseudo elements are considered.
   * @param st if we copy the vector as RAW_VECTOR or element wise
   * @param pde will me mech in SIMP and also elec in PiezoSIMP
   * @param app redundant to pde. NO_APP for non ELEMENT_VECTOR as long as OLAS makes no difference
   * @param save_sol when an adjoint system was solved, one has to call
   *        "SaveSolution()" in the pde such that we can extract it element wise.
   *        Only relevant for st = ELEMENT_VECTORS ---------- FIXE! still required!
   * @return NULL if st = ELEMENT_VECTOR, otherwise it is the vector */
  SingleVector* Read(StorageType st, SinglePDE* pde, Optimization::Application app = Optimization::NO_APP, bool save_sol = false, DERIVType derivative = NO_DERIVTYPE);

  template<class T>
  SingleVector* Read(StorageType st, SinglePDE* pde, Optimization::Application app,  bool save_sol = false, DERIVType derivative = NO_DERIVTYPE);

  /** Writes the solution (raw vector) back to the pde */
  void Write(SinglePDE* pde);

  static void Write(SinglePDE* pde, SingleVector* vec);

  /** average the raw solutions by the excitations.
   * @param excitations average or solution by one entry only. Is strictly speaking already known
   * by the em_ parameter but is more explicit this way.  */
  static void Write(SinglePDE* pde, StateSolutions& sol, Function* f, int time_step, StdVector<Excitation>& excitations);

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

  /** For debug output */
  std::string ToString();

  /** This is an element wise storage of the solution
   * the Application shall be MECH or ELEC */
  std::map<Optimization::Application, StdVector<SingleVector*> > elem;

  /** This is an element wise storage of the solution
   * considering all elements from the grid instead of all design elements only
   * needed by shape optimization */
  std::map<Optimization::Application, StdVector<SingleVector*> > gridelem;

private:

  template<class T>
  void Write(SinglePDE* pde);

  template<class T>
  static void Write(SinglePDE* pde, StateSolutions& sol, Function* f, int time_step, StdVector<Excitation>& excitations);

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

} // namespace


#endif /*STATE_SOLUTION_HH_*/
