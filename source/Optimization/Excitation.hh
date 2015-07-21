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


namespace CoupledField
{
class Assemble;
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
   * For multiple frequencies it does nothing. The actual frequency is chosen by default.  */
  void Apply();

  /** Find the fixed factor, does ignore weighting and does not apply it. */
  double GetFactor(Function* f) const;

  /** Returns GetFactor() * normalized_weight */
  double GetWeightedFactor(Function* f) const;

  /** Gets the current omege =  2 * pi * f */
  double GetOmega() const;

  /** read the tracking node list from XML */
  void ReadTrackings(PtrParamNode ts);

  /** read the loads from XML */
  void ReadLoads(PtrParamNode ls);

  void ReadTestStrain(MechPDE::TestStrain ts);

  void ReadTestCharges(const Vector<double>& vec);

  /** does not label the frequency, load or test strain but the rotation or robust case if present
   * @return "" if not present */
  std::string GetMetaLabel() const;

  /** the full label is meta label plus base label or it is simply label */
  std::string GetFullLabel() const;

  /** the index of this excitation in the excitations array. If -1 something went wrong */
  int index;

  /** the meta index */
  int meta_index;

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

  Assemble* assemble;
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
  /** Do we do multiple excitation at all? */
  bool IsEnabled() const { return multiple_excitation_; }

  /** Handle multiple excitations (loads/frquencies). By definition the size is almost 1, even
   * if there is no load (e.g. static piezo with inhomgeneous Dirichlet BC. */
  void PrepareMultipleExcitations(Optimization* opt, bool eval_inital_design);

  /** Do we do adjust weights */
  bool DoAdjustWeights() const { return type_ == META_OBJECTIVE; }

  bool DoHomogenization() const { return type_ == HOMOGENIZATION_TEST_STRAINS; }

  /** The number of homogenization test strains. Important when we do also transform */
  unsigned int GetNumberHomogenization() const { return DoHomogenization() ? total_base_ : 0; }

  bool DoBloch() const { return domain->GetDriver()->DoBlochModeEigenfrequency(); }

  /** apply excitation specific transformation (rotation) */
  bool DoTransform() const { return num_trans_ > 0; }

  /** handles transform and robust */
  unsigned int GetNumberMeta(bool minimum_one = false) const { return GetNumberTransform(minimum_one); }

  /** The number of transformations. Important when we do homogenization */
  unsigned int GetNumberTransform(bool mininum_one = false) const { return mininum_one ? std::max(num_trans_, 1) : num_trans_; }

  /** Search for the excitation label.
   * @param quiet if true NULL is returned when the label is not found instead of an exception */
  Excitation* GetExcitation(const std::string& label, bool quiet = false);

  /** Gets the excitation and handles transform
   * @param base e.g. for homogenization the number of the teststrain*/
  Excitation* GetExcitation(unsigned int base, Transform* trans);

  /** Gets the excitation based on the meta level. This allows to traverse the meta labels easily
   * @param base e.g. for homogenization the number of the teststrain, typically 0
   * @param meta e.g. the number of the */
  Excitation* GetExcitation(unsigned int base, unsigned int meta);

  /** Gets the excitation based on the meta level. This allows to traverse the meta labels easily
   * @param base e.g. for homogenization the number of the teststrain, typically 0
   * @param meta needs to be a number */
  Excitation* GetExcitation(unsigned int base, const std::string& meta);

  /** The excitation index is not that easy if we have loads/homogenization/frequencies and concurrently robustness and transformations.
   * The functions have excitations for the later but not necessarily for the first
   * @param base the "normal" index of test strains, ...
   * @param f checks for transformation and robustness in the excitation of the function.
   * @see GetExcitation(unsigned int, Transform*) */
   unsigned int GetExcitationIndex(unsigned int base, Function* f);

  /** The meta excitation index considers only the meta level (transformation, robustness) not the base level (frequency, wave, test strain)
   * @return 0 if we have no meta stuff
   * You may also aks Excitation::meta_index*/
  // unsigned int GetMetaExcitationIndex(Function* f);

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

  /** Helper for PrepareMultipleExcitations(). Excitations are set with hard coded test strains */
  int SetHomogenizationTestStrains();

  /** Helper which sets up the transformation based on any exciting excitations (e.g. test strains), which are wrapped and multiplied */
  void ApplyTransformations(DesignSpace* space);

  void SetLoadCases(const ParamNodeList& pn_ex, int num_loads, Optimization* opt);

  void WriteInInfo(int num_freq, bool eval_inital_design, double weight_sum,  Optimization* opt);

  void SetHarmonic(int num_freq);

  void SetBlochWaves(int num_wave);

  int ValidateTransformation(Optimization* opt);

  /** do we do multiple excitation at all? */
  bool multiple_excitation_;

  Type type_;

  /** the base number of excitations (loads, test strains, frequencies) to be multiplied by transformations and robustness */
  unsigned int total_base_;

  /** number of transformations in DesignSpace::transform. This is a meta level*/
  int num_trans_;
};


} // namespace



#endif /* EXCITATION_HH_ */
