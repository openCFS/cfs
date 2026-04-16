#ifndef EXCITATION_HH_
#define EXCITATION_HH_

#include <string>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Domain/BCs.hh"
#include "Driver/HarmonicDriver.hh"
#include "General/Enum.hh"
#include "MatVec/Vector.hh"
#include "PDE/MechPDE.hh"
#include "Utils/StdVector.hh"
#include "Optimization/Context.hh"


namespace CoupledField
{
class Function;
class LinearFormContext;
class ObjectiveContainer;
class SinglePDE;
class Transform;

/** For multiple loads (compliance or multiple frequency optimization) we use
 * the summarized term multiple excitations. This object encapsulates a excitation inclusive weight. */
class Excitation
{
public:

  /** default constructor for StdVector() */
  Excitation();

  ~Excitation();

  /** This method makes the current load active.
   * For multiple frequencies it updates the mathparser, e.g. when we use freq dependent coefFunctions
   * @param switch_context perform a context switch if another sequence is active.
   *    This triggers the deletion of the current pdes and single driver and the creation of new ones. So know what you are doing!
   *    A good reason not to switch the context is when we apply the current robust filter for function evaluation
   * @return true if a context switch happened (what is only possible with switch_context set */
  bool Apply(bool switch_context);

  /** Find the fixed factor, does ignore weighting and does not apply it. */
  double GetFactor(Function* f) const;

  /** Returns GetFactor() * normalized_weight */
  double GetWeightedFactor(Function* f) const;

  /** Gets the current omege =  2 * pi * f */
  double GetOmega() const;

  /** read the tracking node list from XML */
  void ReadTrackings(PtrParamNode ts);

  /** read the loads from XML */
  void ReadLoads(Context* ctxt, PtrParamNode ls);

  void ReadTestStrain(Context* ctxt, MechPDE::TestStrain ts);

  /** does not label the frequency, load or test strain but the rotation or robust case if present
   * @return "" if not present */
  std::string GetMetaLabel() const;

  /** the full label is meta label plus base label or it is simply label */
  std::string GetFullLabel() const;

  /** Is this a Bloch Excitatztion? */
  bool DoBloch() const;

  /** the 0-based wave number is for the first sequence the index */
  int GetWaveNumber() const;

  /** In the case of buckling, we have two excitations, each with one pde,
   * first linear elasticity (real), then buckling (complex eigenvalue problem).
   * In optimization we take the solution of the first pde, convert it to
   * complex and inject it in the second one to get the stresses which we need
   * to assemble the buckling pde.
   *  */
  std::map<RegionIdType, PtrCoefFct> GetStressCoefFctFromExcitation(UInt index);
  void SetStressCoefFct(std::map<RegionIdType, PtrCoefFct> stress_map);

  /** the index of this excitation in the excitations array. If -1 something went wrong */
  int index;

  /** the meta index */
  int meta_index;

  /** the corresponding sequence step. 1-base. @see ContextManager */
  int sequence;

  /** For several loads, we need to store the form context with the entities but also the dof and the value!
   * When no excitations are given in the optimization part the bcsAndLoads are used.
   * For the bcsAndLoads case there is one form per excitation and the ownership (destruction) is taken
   * from ~Assemble() to ~Excitation() via Assemble::lin_forms_given_.
   * When excitations are used the original Assemble::linForms_ are deleted and multiple forms are allowed per
   * excitation. see Excitation::ReadLoads() */
  StdVector<LinearFormContext*> forms;

  /** This is a link to the Frequency description from the harmonic driver.
   * It is used for calling the HarmonicDriver to solve the problems */
  HarmonicDriver::Frequency* f_link;

  /** the frequency index, should correspond to the freq_step of the harmonic driver */
  unsigned int freq_idx;

  /** For multiharmonic excitation. -1.0 by default */
  double frequency;

  /** For Bloch mode analysis. The excitation index is the wave vector index.  */
  Vector<double> wave_vector;

  /** this is the weight from the xml file.
   * @see normalized_weight */
  double weight;

  /** do we need to reassemble the system? True for frequency, robust and transformation */
  bool reassemble;

  /** this is the normalized weight (sum of all weights of all excitations is 1).
   * Note that for functions with a single excitation (e.g. stress, volume) it shall be 1 */
  double normalized_weight;

  /** A label denoting the excitation, depends on kind */
  std::string label;

  /** Here we store the calculated objective value, including costFunction/factor
   * to enable metaObjective. In the Optimization::cost struct we store a copy/average */
  double cost;

  /** Here we might store the "original" (@see objective) gradient for analysis/output */
  StdVector<double> cost_gradient;

  /** for the calculation of a homogenized material tensor, we use the test
   * strains defined in ErsatzMaterial::SetHomogenizationTestStrains() */
  Vector<double> test_strain;

  /** For transformation, this is the applied transformation. NULL means we don't do transformation*/
  Transform* transform;

  /** this is the robust index if we do index. Otherwise it is zero which is the standard filter index for non-robust */
  unsigned int robust_filter_idx;

  /** When we do robust, the meta_index is only the robust index when we do no concurrent transformation */
  bool robust;
};

/** This struct stores the multiple excitation Information. It contains the
 * same defaults as the xml schema as the whole element is optional. */
class MultipleExcitation
{
public:
  /** @param enabled comes from the attribute.
   * @param pn if NULL only the defaults are set */
  MultipleExcitation(bool enabled, PtrParamNode pn);

  void ToInfo(PtrParamNode in) const;

  typedef enum { NO_TYPE, FIXED_WEIGHT, META_OBJECTIVE, HOMOGENIZATION_TEST_STRAINS} Type;

  static Enum<Type> type;

  /** Do we do multiple excitation? Implicitly true for bloch for the matching sequence
   * @param sequence -1 at all? otherwise for the specified sequence */
  bool IsEnabled(int sequence = -1) const;

  /** the sequence the multiple excitation is defined for. Extend for ME for more than one sequence.
   * For bloch me is defined and enabled implicitly */
  int GetSequence() const { return sequence_; }

  /** To be called prior to PrepareMultipleExcitations() */
  void InitializeMultipleExcitations(Optimization* opt, ContextManager* manager);

  /** Handle multiple excitations (loads/frequencies). By definition the size is almost 1, even
   * if there is no load (e.g. static piezo with inhomgeneous Dirichlet BC.
   * @param ctxt an own version of Context to setup a multi sequence system */
  void PrepareMultipleExcitations(Optimization* opt, Context* ctxt);

  /** To be called after PrepareMultipleExcitations() */
  void FinalizeMultipleExcitations(Optimization* opt, ContextManager* manager, bool eval_inital_design);

  /** Do we do adjust weights */
  bool DoAdjustWeights() const { return type_ == META_OBJECTIVE; }

  bool DoHomogenization() const { return type_ == HOMOGENIZATION_TEST_STRAINS; }

  /** The number of homogenization test strains. Important when we do also transform */
  unsigned int GetNumberHomogenization(App::Type app) const;

  /** apply excitation specific transformation (rotation) */
  bool DoTransform() const { return num_trans_ > 0; }

  /** Do we do robust */
  bool DoRobust(const Context* ctxt) const;

  /** Do we do real meta excitation?
   * @param ctxt might be null the we check for any */
  bool DoMetaExcitation(const Context* ctxt) const;
  bool DoMetaExcitation(int sequence) const;

  /** handles transform and robust.
   * @param ctxt if NULL search max
   * @param minimum_one if false take care as num_robust can be 0 or 1 w/o robust */
  unsigned int GetNumberMeta(const Context* ctxt, bool minimum_one = false) const;

  /** The number of transformations. Important when we do homogenization */
  unsigned int GetNumberTransform(bool mininum_one = false) const { return mininum_one ? std::max(num_trans_, 1) : num_trans_; }

  /** @param ctxt if NULL search for max (is either none or the same number */
  unsigned int GetNumberRobust(const Context* ctxt, bool mininum_one = false) const;

  /** count the number of unique frequencies for the context. Static is 0 */
  unsigned int GetUniqueFrequencies() const;

  /** Search for the excitation label.
   * @param quiet if true NULL is returned when the label is not found instead of an exception */
  Excitation* GetExcitation(const std::string& label, bool quiet = false);

  /** For doing adjust weights when doing multiple excitation with meta objective, this method
   * does the job. It requires the cost entries in excitations to be set.
   * The \f$w_k^p=const\;\sum w_k = 1\f$ condition is fulfilled here. */
  void NormalizeMultipleExcitations(ObjectiveContainer* objectives);

  /** Here we have the set of excitations. Only relevant for the multiple excitations
   * case (multiple loads or frequencies). Set e.g. by ErsatzMaterial */
  StdVector<Excitation> excitations;

  /** The stride for adjust weights: 1 = every iteration, 2 = every second ... */
  int stride;

  /** the maximal span between the largest (1) and smallest weight as factor */
  double max_gain;

  /** The exponent d in w_k^p J_k = const */
  double damping;

private:

  /** Helper for PrepareMultipleExcitations(). Excitations are set with hard coded test strains
   * @param context_base the size of excitations intially for the current context (0 in single or first sequence case) */
  int SetHomogenizationTestStrains(unsigned int context_base, Context* ctxt);

  /** Helper which sets up the robust filters based on any exciting excitations (e.g. test strains), which are wrapped and multiplied */
  void ApplyRobust(const Context* ctxt);

  /** Helper which sets up the transformation based on any exciting excitations (e.g. test strains) including robust!!!, which are wrapped and multiplied */
  void ApplyTransformations(const Context* ctxt, DesignSpace* space);

  /** Helper to find the corresponding form to id out of forms */
  LinearFormContext* SearchFormByCoilId(StdVector<LinearFormContext*>& forms, const string& id);

  /** Set load for coil case */
  void SetCoils(unsigned int base, Assemble* ass, const ParamNodeList& pn_ex, int num_loads, MathParser* parser, unsigned int handle);

  /** Reads loads from the boundary conditions or from the optimization part.
   * Handles the case that we have multiple loads (e.g. magnetic coils) for a single frequency (num_freq = 1) */
  void SetLoadCases(Context* ctxt, unsigned int context_base, const ParamNodeList& pn_ex, int num_loads, Optimization* opt);

  /** call this only for the last context */
  void WriteInInfo(PtrParamNode in);

  /** Sets the "loads" for multifrequency harmonic problems
   * Single frequency is handled in SetLoadCases */
  void SetHarmonic(Context* ctxt, unsigned int context_base, const ParamNodeList& pn_ex, int num_freq);

  /** sweet little helper for SetHarmonic() */
  void SetHarmonicExcitation(Context* ctxt, Excitation& ex, int freq_idx);

  /** @see SetHomogenizationTestStrains() */
  void SetBlochWaves(Context* ctxt, unsigned int context_base, int num_wave);

  int ValidateTransformation(Optimization* opt);

  /** count the current excitations by sequence */
  unsigned int CountExcitations(const Context* ctxt) const;

  /** return the indices of the matching excitations */
  StdVector<unsigned int> GetExcitations(const Context* ctxt) const;

  /** do we do multiple excitation at all? */
  bool multiple_excitation_;

  Type type_;

  /** the optional sequence attribute for multi sequence optimization. Note we assume only one multiExcitation xml element for a single sequence
   * in the multi sequence case. E.G. bloch mode as sequence one and homogenization as sequence two */
  int sequence_; // 1-based!!

  /** the principle number of excitations (loads, test strains, frequencies) to be multiplied by transformations and robustness */
  unsigned int principle_;

  /** number of transformations in DesignSpace::transform. This is a meta level*/
  int num_trans_;

  /** for every sequence we have an entry. Maps to multipleExcitation/robust in xml */
  struct Robust
  {
    int num_robust = 0; // 0 means disabled
    int alt_filter = -1; // the filter to use when we are not robust. -1 for no
  };

  StdVector<Robust> robust_;
};


} // namespace



#endif /* EXCITATION_HH_ */
