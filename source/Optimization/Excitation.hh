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

class Function;
class LinearFormContext;
class ObjectiveContainer;
class SinglePDE;

/** For multiple loads (compliance or multiple frequency optimization) we use
 * the summarized term multiple excitations. This object encapsulates a excitation inclusive weight. */
class Excitation
{
public:

  /** default constructor for StdVector() */
  Excitation();

  /** This method makes the current load active.
   * For multiple frequencies it does nothing. The actual frequency is chosen by default. */
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

  void AddLinFormsFromAssemble();

  /** return pointer to linForms, used by Shape-Optimization */
  StdVector<LinearFormContext*>& GetLinForms() { return *linForms; }

  /** the index of this excitation in the excitations array. If -1 something went wront */
  int index;

  /** For static/single harmonic optimization with different load-cases. */
  StdVector<shared_ptr<EntityList> > loads;

  /** For static optimization with different pressure or regionLoads */
  StdVector<LinearFormContext*>* linForms;

  /** This is a link to the Frequency description from the harmonic driver.
   * It is used for calling the HarmonicDriver to solve the problems */
  HarmonicDriver::Frequency* f_link;

  /** For multiharmonic excitation. -1.0 by default */
  double frequency;

  /** this is the weight from the xml file.
   * @see normalized_weight */
  double weight;

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
  void PrepareMultipleExcitations(SinglePDE* pde, PtrParamNode optInfoNode, bool harmonic, bool eval_inital_design);

  /** Do we do adjust weights */
  bool DoAdjustWeights() const { return type_ == META_OBJECTIVE; }

  bool DoHomogenization() const { return type_ == HOMOGENIZATION_TEST_STRAINS; }

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

  /** Helper for PrepareMultipleExcitations(). Excitations are set with hard coded test strains */
  int SetHomogenizationTestStrains();

  /** do we do multiple excitation at all? */
  bool multiple_excitation_;
  Type type_;
};


} // namespace



#endif /* EXCITATION_HH_ */
