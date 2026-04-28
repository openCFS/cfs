#ifndef OPTIMIZATION_CONTEXT_HH_
#define OPTIMIZATION_CONTEXT_HH_

//include <cstddef>
#include "Utils/StdVector.hh"
#include "PDE/BasePDE.hh"

namespace CoupledField
{
class ContextManager;
class DesignMaterial;
class Excitation;
class SingleDriver;
class SinglePDE;
class BucklingDriver;
class EigenFrequencyDriver;
class EigenValueDriver;
class HarmonicDriver;
class Exctiation;
class MultipleExcitation;
class Function;
class OptimizationMaterial;
class LatticeBoltzmannPDE;
class BiLinFormContext;

struct App
{
  /** The App::Type type identifies the PDE to use.
   *  A subset of the values are PDE identifiers for Context::ToPDE() and Context::ToApp().*/
  typedef enum { MECH, ELEC, PIEZO_COUPLING, PRESSURE, CHARGE_DENSITY, MASS, HEAT, ACOUSTIC, LAPLACE, STRESS, LBM, MAG, BUCKLING, NO_APP} Type;
};


/** The context describes where we are within the current optimization process.
 *
 * A basic information is the currently applied Excitation or for a multisequence
 * problem the current sequence interest.
 *
 * @see static Optimization::context
 * @see ContextManager
 *
 * By purpose we have minimal includes to avoid circular includes */
class Context
{
  public:

  Context();

  ~Context();

  /** the real optimizer to be called once we have domain initialized.
   * In the multisequence case the drivers are only initialized on request */
  void Setup(ContextManager* manager, BasePDE::AnalysisType analyis, PtrParamNode node, int sequence_step);

  /** We initialize after each context switch to get the current driver and pdes */
  void Update();

  /** this holds the current excitation. */
  Excitation* GetExcitation() { return excitation_; }

  /** this might trigger SwitchSequence() */
  void SetExcitation(Excitation* ex);

  /** Gets the excitation based on the meta level. This allows to traverse the meta labels easily
   * @param base e.g. for homogenization the number of the teststrain, typically 0
   * @param meta e.g. the number of the */
  Excitation* GetExcitation(unsigned int base, unsigned int meta);

  /** Gets the excitation based on the meta level. This allows to traverse the meta labels easily
   * @param base e.g. for homogenization the number of the teststrain, typically 0
   * @param meta needs to be a number */
  Excitation* GetExcitation(unsigned int base, const std::string& meta);

  /** The excitation is not that easy to find if we have loads/homogenization/frequencies and concurrently robustness and transformations.
   * The functions have excitations for the later but not necessarily for the first
   * @param base the "normal" index of test strains, ...
   * @param f checks for transformation and robustness in the excitation of the function.
   * @see GetExcitation(unsigned int, Transform*) */
   Excitation* GetExcitation(unsigned int base, Function* f);

  /** might change in the multi sequence case so don't cache! */
  SingleDriver* GetDriver() { return driver; }

  /** shortcut for dynamic_cast to eigenfrequency driver. Does no checking, so check for return NULL
   * TODO make this stuff inline!
   * @return NULL if other driver. */
  EigenFrequencyDriver* GetEigenFrequencyDriver();

  BucklingDriver* GetBucklingDriver();

  EigenValueDriver* GetEigenValueDriver();

  HarmonicDriver* GetHarmonicDriver();

  LatticeBoltzmannPDE* GetLatticeBoltzmannPDE();

  /** Do we have a harmonic problem? Then we are complex. Even if not, we might be eigenvalue and also complex*/
  bool IsHarmonic() const { return harmonic_; }

  /** we are complex in the harmonic or eigenvalue case */
  bool IsComplex() const { return complex_; }

  /** do we solve eigenvalue problems? Then we are complex! */
  bool IsEigenvalue() const { return eigenvalue_; }

  /** Do we Bloch Mode analysis? */
  bool DoBloch() const { assert((bloch_ && num_bloch_wave_vectors > 0) || (!bloch_ && num_bloch_wave_vectors == 0)); return bloch_; }

  bool DoLBM() const {return (ToApp() == App::LBM);}

  bool DoBuckling() const { return buckling_; }

  /** the driver steps: 1 for static, numFreq for harmonic and wave numbers for bloch */
  unsigned int GetDriverSteps() const { assert(driver_steps_ > 0); return driver_steps_; }

  /** are we within an multi sequence optimization */
  bool DoMultiSequence() const;

  /** Is the current system test strain excited? True for special test case and homogenization */
  bool IsStrainExcitedSystem() const;

  /** Set the multiple excitations specific information. The excitations array is set independently. */
  void SetMultiExcitations(MultipleExcitation* me, unsigned int basic_excitaions);

  /** Helper that converts from mechPDE to App::MECH and elecPDE to App::ELEC, ...
   * @param from heat and acoustic the application for the transfer function is laplace, this is indicated by the flag if
   *        we do not want a marker for the pde but the transfer function. Sorry, very messy !! :((
   * @throws if neither mechPDE nor elecPDE
   * @see ToPDE()
   * @see SetPDEs() */
  static App::Type ToApp(const SinglePDE* pde);

  App::Type ToApp() const { assert(pde != NULL); return ToApp(pde); }

  /** Find our PDE in SIMP by application from the pdes map
   * @see ToApp()*/
  SinglePDE* ToPDE(App::Type app, bool throw_exception = true);

  /** Service function. A linear form would have a similar implementation */
  BiLinFormContext* GetBiLinFormContext(const RegionIdType reg, App::Type app1, App::Type app2, bool throw_exception);

  /** the corresponding 1-based multi sequence step or -1 of not set yes */
  int sequence;

  /** is this the currently active context? */
  bool active;

  /** the 0-based index within ContextManager::context. Is sequence -1 if set. -1 if not set yet */
  int context_idx;

  /** Is only set when active, otherwise it is zero, don't cache! In the multiple sequence case the driver is always newly created! */
  SingleDriver* driver;

  /** the pde dependend optimization material. Due to multiple sequence needs to wait for Up */
  OptimizationMaterial* mat;

  /** pde dependend material for parametric optimization, owned by DesignSpace **/
  DesignMaterial* dm = NULL;

  /** our analysis type */
  BasePDE::AnalysisType analysis;

  /** Do we do homogenization? */
  bool homogenization;

  /** number of frequencies in harmonic driver case */
  unsigned int num_harm_freq;

  /** The number of loads, either by bcsAndLoads or multipleExciatation*/
  unsigned int num_loads;

  /** number of bloch wave vectors in bloch eigenfrequency case */
  unsigned int num_bloch_wave_vectors;

  /** number of eigenmodes - also in the non-bloch case! */
  unsigned int num_eigenmodes;

  /** the SubTensorType applies only for MechPDE */
  SubTensorType stt;

  /** context specific excitations (at least one).
   * Is a reference to portions of Optimization::MultipleExcitation::excitation and set in MultipleExcitation::FinalizeMultipleExcitation() */
  StdVector<Excitation*> excitations;

  /** The pdes from the current sequence state, don't cache!
   * Note that for multiple sequence optimization the pdes are always newly created and we must not store the pointer!
   * The order of the pdes is not defined, Therefore we use the map. Set via SetPDEs()
   * @see ToApp()
   * @see ToPDE() */
  std::map<App::Type, SinglePDE*> pdes;

  /** This is simple one SinglePDE from pdes, don't cache! */
  SinglePDE* pde;

  PtrParamNode infoNode;

private:

  /** make shortcuts for the currently available PDEs in pdes.
   * @return fals if it was too early and the pdes are not set in the system yet */
  bool SetPDEs();

  /** are we harmonic/EV or static/transient? */
  bool complex_;

  /** only for the driver, not for complex_! */
  bool harmonic_;

  /** do we solve an eigenvalue problem. Includes bloch mode problems */
  bool eigenvalue_;

  /** bloch mode analysis is also eigenvalue but special due to the wave vectors encapsulated in excitations */
  bool bloch_;

  /** buckling analysis is also eigenvalue but special due to the system matrix of the second excitation depending on the first*/
  bool buckling_;

  /** we read the driver steps even without driver object to allow PrepareMultipleExcitation() */
  unsigned int driver_steps_;

  /** the currently active excitation */
  Excitation* excitation_;

  ContextManager* manager_;

  /** this is the number of basic excitations without meta (robust/ transformation). Hence 3/6 for homogenization.
   * It cannot be larger than excitations.GetSize() and equals excitations.GetSize() without robust/ transformation */
  unsigned int basic_excitations_;

  /** Link to MultipleExcitations. Actually we could own our own instance ... */
  MultipleExcitation* me_;

};

/** we have only one static instance of the context manager in Optimization::contextManager.
 * We have the optimization case and the DesignSpace w/o optimization case.
 * the MultiSequenceDriver creates the sequence (driver/ pdes) on demand and destroys them afterwards. */
class ContextManager
{
public:
  /** called via static initialization before there is a domain with driver */
  ContextManager();

  ~ContextManager();

  /** initialize after cfs has initialized the drivers */
  void Init();

  bool IsInitialized() const { return initialized_; }

  /** sets and initializes the context as active.
   * Resets the information from the other context as the single driver and the pde(s) will be newly created!
   * @param index 0-based */
  void SwitchContext(int index);

  /** just to prevent the 0-based index issue. Also sets the excitation
   * @see SwitchContext(int) */
  void SwitchContext(Excitation* excitation);

  /** gives the context corresponding to the function.
   * Simply used the 1-based sequence attribute of the function */
  Context& GetContext(const Excitation* ex);

  /** returns the homogenization context or NULL if there is none */
  Context* GetHomogenization();

  StdVector<Context> context;

  /** this is a max norm of all context */
  struct AnyContext
  {
    AnyContext();

    bool bloch;
    bool eigenvalue;
    bool harmonic;
    bool complex;
    bool homogenization;
  };

  /** freshly generated summary of all contexts. Fast */
  const AnyContext any() const;

private:

  // we need a default Excitation in the case of load ersatzmaterial only for the default filter index
  Excitation* default_excitation_;

  bool initialized_;
};

} // end of namespace


#endif /* OPTIMIZATION_CONTEXT_HH_ */
