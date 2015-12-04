#ifndef OPTIMIZATION_CONTEXT_HH_
#define OPTIMIZATION_CONTEXT_HH_

//include <cstddef>
#include "Utils/StdVector.hh"
#include "PDE/BasePDE.hh"

namespace CoupledField
{
class ContextManager;
class Excitation;
class SingleDriver;
class SinglePDE;
class EigenFrequencyDriver;
class HarmonicDriver;
class Exciation;
class Function;
class OptimizationMaterial;

struct App
{
  /** The App::Type type identifies the PDE to use.
   *  A subset of the values are PDE identifiers for Context::ToPDE() and Context::ToApp().
   * The heat and acoustic transfer functions are Laplace! */
  typedef enum { MECH, ELEC, PIEZO_COUPLING, PRESSURE, CHARGE_DENSITY, MASS, HEAT, ACOUSTIC, LAPLACE, STRESS, LBM, NO_APP} Type;
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

  /** might change in the multi sequence case so don't cache! */
  SingleDriver* GetDriver() { return driver; }

  /** shortcut for dynamic_cast to eigenfrequency driver. Does no checking, so check for return NULL
   * @return NULL if other driver. */
  EigenFrequencyDriver* GetEigenFrequencyDriver();

  HarmonicDriver* GetHarmonicDriver();

  /** Do we have a harmonic problem? Then we are complex. Even if not, we might be eigenvalue and also complex*/
  bool IsHarmonic() const { return harmonic_; }

  /** we are complex in the harmonic or eigenvalue case */
  bool IsComplex() const { return complex_; }

  /** do we solve eigenvalue problems? Then we are complex! */
  bool IsEigenvalue() const { return eigenvalue_; }

  /** Dow we Bloch Mode analysis? */
  bool DoBloch() const { return bloch_; }

  /** the driver steps: 1 for static, numFreq for harmonic and wave numbers for bloch */
  unsigned int GetDriverSteps() const { assert(driver_steps_ > 0); return driver_steps_; }

  /** are we within an multi sequence optimization */
  bool DoMultiSequence() const;

  /** Helper that converts from mechPDE to App::MECH and elecPDE to App::ELEC, ...
   * @param from heat and acoustic the application for the transfer function is laplace, this is indicated by the flag if
   *        we do not want a marker for the pde but the transfer function. Sorry, very messy !! :((
   * @throws if neither mechPDE nor elecPDE
   * @see ToPDE()
   * @see SetPDEs() */
  static App::Type ToApp(const SinglePDE* pde);

  App::Type ToApp() const { return ToApp(pde); }

  /** Find our PDE in SIMP by application from the pdes map
   * @see ToApp()*/
  SinglePDE* ToPDE(App::Type app, bool throw_exception = true);

  /** the corresponding 1-based multi sequence step or -1 of not set yes */
  int sequence;

  /** the 0-based index within ContextManager::context. Is sequence -1 if set. -1 if not set yet */
  unsigned int context_idx;

  /** Is only set when active, otherwise it is zero, don't cache! In the multiple sequence case the driver is always newly created! */
  SingleDriver* driver;

  /** the pde dependend optimization material. Due to multiple sequence needs to wait for Up */
  OptimizationMaterial* mat;

  /** our analysis type */
  BasePDE::AnalysisType analysis;

  /** number of frequencies in harmonic driver case */
  unsigned int num_harm_freq;

  /** number of bloch wave vectors in bloch eigenfrequency case */
  unsigned int num_bloch_wave_vectors;

  /** context specific excitations (at least one).
   * Is a reference to portions of Optimization::MultipleExcitation::excitation and set in MultipleExcitation::FinalizeMultipleExcitation() */
  StdVector<Excitation*> excitation;

  /** The pdes from the current sequence state, don't cache!
   * Note that for multiple sequence optimization the pdes are always newly created and we must not store the pointer!
   * The order of the pdes is not defined, Therefore we use the map. Set via SetPDEs()
   * @see ToApp()
   * @see ToPDE() */
  std::map<App::Type, SinglePDE*> pdes;

  /** This is simple one SinglePDE from pdes, don't cache! */
  SinglePDE* pde;

private:

  /** make shortcuts for the currently available PDEs in pdes */
  void SetPDEs();

  /** are we harmonic/EV or static/transient? */
  bool complex_;

  /** only for the driver, not for complex_! */
  bool harmonic_;

  /** do we solve an eigenvalue problem. Includes block mode problems */
  bool eigenvalue_;

  /** bloch mode analyis is also eigenvalue but special due to the wave vectors encapsualted in excitations */
  bool bloch_;

  /** we read the driver steps even without driver object to allow PrepareMultipleExcitation() */
  unsigned int driver_steps_;

  Excitation* excitation_;

  ContextManager* manager_;
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
  const Context& GetContext(const Function* f) const;

  StdVector<Context> context;

  /** this is a max norm of all context */
  struct AnyContext
  {
    bool bloch;
    bool eigenvalue;
    bool harmonic;
    bool complex;
  };

  /** this read ony struct holds a max norm like summary of all context */
  const AnyContext& any() const { return any_; };

private:

  // we need a default Excitation in the case of load ersatzmaterial only for the default filter index
  Excitation* default_excitation_;

  bool initialized_;

  AnyContext any_;
};

} // end of namespace


#endif /* OPTIMIZATION_CONTEXT_HH_ */
