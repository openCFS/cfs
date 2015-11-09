#ifndef OPTIMIZATION_CONTEXT_HH_
#define OPTIMIZATION_CONTEXT_HH_

#include <cstddef>

namespace CoupledField
{

class Excitation;
class SingleDriver;
class EigenFreuqencyDriver;

/** The context describes where we are within the current optimization process.
 *
 * A basic information is the currently applied Excitation or for a multisequence
 * problem the current sequence interest.
 * @see static Optimization::context
 *
 * By purpose we have minimal includes to avoid circular includes */
class Context
{
  public:

  Context();

  ~Context();

  /** the real optimizer to be called once we have domain initialized. */
  void Init();

  bool IsInititialized() const { return driver_ != NULL; }

  /** this holds the current excitation. */
  Excitation* GetExcitation() { return excitation_; }

  /** this might trigger SwitchSequence() */
  void SetExcitation(Excitation* ex);

  /** might change in the multi sequence case so don't cache! */
  SingleDriver* GetDriver() { return driver_; }

  /** shortcut for dynamic_cast to eigenfrequency driver. Does no checking, so check for return NULL
   * @return NULL if other driver. */
  EigenFrequencyDriver* GetEigenFrequencyDriver() { return dynamic_cast<EigenFrequencyDriver*>(driver_); }

  /** Do we have a harmonic problem? Then we are complex. Even if not, we might be eigenvalue and also complex*/
  bool IsHarmonic() const { return harmonic_; }

  /** we are complex in the harmonic or eigenvalue case */
  bool IsComplex() const { return complex_; }

  /** do we solve eigenvalue problems? Then we are complex! */
  bool IsEigenvalue() const { return eigenvalue_; }

  /** Dow we Bloch Mode analysis? */
  bool DoBloch() const { return bloch_; }

  /** are we within an multi sequence optimization */
  bool MultiSequence() const;

private:

  /** triggered by SetExcitiation */
  void SwitchSequence();

  /** sets the current context stuff based on current driver */
  void ParseSequence();

  /** are we harmonic/EV or static/transient? */
  bool complex_;

  /** only for the driver, not for complex_! */
  bool harmonic_;

  /** do we solve an eigenvalue problem. Includes block mode problems */
  bool eigenvalue_;

  /** bloch mode analyis is also eigenvalue but special due to the wave vectors encapsualted in excitations */
  bool bloch_;

  Excitation* excitation_;

  /** to have not memory responsibility issues */
  Excitation* default_excitation_;

  SingleDriver* driver_;

};

} // end of namespace


#endif /* OPTIMIZATION_CONTEXT_HH_ */
