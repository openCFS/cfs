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
class SingleVector;


/** This class holds the solution of the PDE. It is in a class such that it
 * helps to encapsulate real and complex solutions. Note that the Piezo
 * has other solutions! */
class StateSolution
{
public:

  StateSolution();

  ~StateSolution();

  typedef enum
  {
    ELEMENT_VECTORS = 0, RAW_VECTOR = 1, RHS_VECTOR = 2, SEL_VECTOR = 3, GRIDELEM_VECTORS = 4
  } StorageType;

  static SolutionType GetSolutionType(SinglePDE* pde, App::Type app = App::NO_APP);

  /** Copies the solution for the pde in our own storage.
   * In the ELEMENT_VECTORS case make sure, that the solution is in the PDE!
   * For manual adjoint stuff you might have to do SaveSolution() first!
   * In ELEMENT_VECTORS also the registered pseudo elements are considered.
   * @param st if we copy the vector as RAW_VECTOR or element wise
   * @param pde will me mech in SIMP and also elec in PiezoSIMP
   * @param app redundant to pde. App::NO_APP for non ELEMENT_VECTOR as long as OLAS makes no difference
   * @param save_sol when an adjoint system was solved, one has to call
   *        "SaveSolution()" in the pde such that we can extract it element wise.
   *        Only relevant for st = ELEMENT_VECTORS ---------- FIXE! still required!
   * @return NULL if st = ELEMENT_VECTOR, otherwise it is the vector */
  SingleVector* Read(StorageType st, SinglePDE* pde, App::Type app = App::NO_APP, bool save_sol = false, TimeDeriv derivative = NO_DERIVTYPE);

  template<class T>
  SingleVector* Read(StorageType st, SinglePDE* pde, App::Type app,  bool save_sol = false, TimeDeriv derivative = NO_DERIVTYPE);

  /** do we contain the state solution? Already read? */
  bool ContainsState() const { return set_; }

  /** Writes the solution (raw vector) back to the pde */
  void Write(SinglePDE* pde);

  static void Write(SinglePDE* pde, SingleVector* vec);

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
   * the App::Type shall be App::MECH or App::ELEC */
  std::map<App::Type, StdVector<SingleVector*> > elem;

  /** This is an element wise storage of the solution
   * considering all elements from the grid instead of all design elements only
   * needed by shape optimization */
  std::map<App::Type, StdVector<SingleVector*> > gridelem;

  /** For identification: for adjoint problems the corresponding function. NULL indicates forward state problem */
  const Function* func;

  /** For identification: the derivative time is only for transient problems */
  TimeDeriv derivative;

  /** For identification: the timestep is only for transient problems, the mode is for eigenmodes. We use this as a union.
   * For eigenvalue problems (hence also Bloch) one solution has as many states as
   * eigenmodes and eigenvalues are calculated. With Bloch mode all modes share one excitation */
  int timestep_mode;

  /** For identification: This is the active excitation when storing the state. Contains also the sequence */
  const Excitation* excitation;

  /** for eigenvalue problems, this is the frequency (not the eigenvalue)
   * @see StateSolutions::CollectEigenfrequencies() */
  double eigenfreq;

private:

  template<class T>
  void Write(SinglePDE* pde);

  /** common helper for the Get*Vector() stuff */
  template<class T>
  SingleVector* GetVector(StorageType st, bool create);

  /** This is the algsys solution vector. */
  SingleVector* raw;
  //std::map<App::Type, SingleVector* > raw;

  /** This stores the right hand sides in raw format.
   * In the adjoint case this might be based on select for output like objectives.
   * In the forward case this stores the rhs from the forward excitation to perform multiple
   * load cases. */
  SingleVector* rhs;
  // std::map<App::Type, SingleVector* > rhs;

  /** For output like objectives this stores the selection l resp. L which is found by
   * an artificial load-case in the adjoint section. Otherwise it is not used.
   * This kind of base for the rhs information which is additionally multiplied by -1, -1 u^*, ...*/
  SingleVector* select;
  //std::map<App::Type, SingleVector* > select;

  /** read already performed? */
  bool set_;
};


/** As the solutions come for multiple excitations in sets we store the list and the
 * averaged (when multiple excitations are enabled). W/o is just some overhead with data size = 1.
 *
 * For eigenvalue problems the modes are seen as timestep_mode in StateSolutions. The excitations are here only is we
 * do Bloch mode optimization where an excitation is a wave_vector.  */
class StateContainer
{
public:
  StateContainer();

  /** The solution is identified by Function, excitation index (0-based) and time step.
   * Similar to ParamNode the requested item is created if it does not exist yet!.
   * If data_ was not initialized (capacity set, does it on the fly)
   * @param f the function is NULL for the forward problems, for the adjoints it needs to be given!
   * @param timestep_mode only for transient or eigenvalue problems */
  StateSolution* Get(const Excitation& excitation, const Function* f = NULL, int timestep_mode = -1, TimeDeriv derivative = NO_DERIVTYPE);
  StateSolution* Get(const Excitation* excitation, const Function* f = NULL, int timestep_mode = -1, TimeDeriv derivative = NO_DERIVTYPE);

  /** Averages the specified states over all excitations.
   * @param sequence 1-based context
   * @param func NULL for forward state solution */
  void WriteAverage(SinglePDE* pde, int sequence, const Function* func = NULL);

  /** Returns the currently stored functions (uniquely). Empty for forward */
  StdVector<const Function*> GetFunctions() const;

  /** Helper which collects the eigenfrequencies for a wave_vector. If no Bloch the wave_vector is 0.
   * @param sequence 1-based*/
  StdVector<double> CollectEigenfrequencies(Excitation& ex);

private:

  /** returns all matching items.
   * @param excitation might be NULL for WriteAverage()
   * @param sequence the context as encoded in StateSolution::excitation->sequence. -1 == unspecific. For WriteAverage() where no specific excitation can be given */
  StdVector<StateSolution*> Search(const Excitation* excitation, int sequence = -1, const Function* f = NULL, int timestep_mode = -1, TimeDeriv derivative = NO_DERIVTYPE);

  /** @see WriteAverage() */
  template<class T>
  void WriteAverage(SinglePDE* pde, int sequence, const Function* func);


  /** the data is actually conditionally multidimensional; by adjoint function, by excitation, by
   * eigenmode, by timestep in transient optimization, by time drivative in transient optimization.
   * It proved too complex to map this for direct access, therefore we simply make an array and perform
   * sequential search. This shall be fast enough as searching the state is negligible compared to
   * calculating the state :). As we return pointers to the states they must not be re-mapped! Therefore
   * it is essential to reserve a sufficiently large array, done in Get() */
  StdVector<StateSolution> data_;
};



} // namespace


#endif /*STATE_SOLUTION_HH_*/
