#ifndef OPTIMIZATION_TUNE_HH_
#define OPTIMIZATION_TUNE_HH_

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Utils/StdVector.hh"
#include "General/Enum.hh"

namespace CoupledField
{
class Optimization;
struct GlobalFilter;
class StoppingRule;

/** allows parameter tuning of iterations, e.g. beta for three field projection or transferFunction penalization.
 * See Dunning and Wein, Achieving near binary topology optimization solutions via automatic projection parameter increase, 2024
 */
class Tune
{
public:
  typedef enum {
    NO_METHOD,
    OBJ,          /** grow times relative objective change, see paper */
    MULT,         /** multiply grow (2 is old fashioned doubling at stride 50) */
    ADD } Method; /** add up grow */

  typedef enum {
    NO_USAGE,
    BETA,            /** beta for density projection */
    PENALTY,         /** param for transfer function */
    FUNC_SCALE } Usage; /** scale for function in multio objective case (0 tolerant!) */

  /** empty constructor to allow instances. We have with GlobalFilter and TransferFunction
   * objects which are copied. Use Init() to init and Register() to activate. */
  Tune() {}

  /** the actual constructor
   * @param pn 'tune' element to be parsed. If not set, an descriptive error is thrown */
  void Init(PtrParamNode pn, Usage use);

  bool IsSet() const { return usage_ != NO_USAGE; }

  /** this is a PostInit() to avoid the GlobalFilter copy constructor issue.
   * @param value pointer to the value to be tuned. Is immediately set to start value
   * @param opt here we register ourself to get the iteration update and we set grayness */
  void Register(double* value, Optimization* opt, GlobalFilter* gf = nullptr);

  /** we check with Optimization::tunes to not mix up with a copied objection which was registered */
  bool IsRegistered() const;

  /** e.g. other robust filters connect to this Tune. Sets the value content */
  void Append(double* value, GlobalFilter* f = nullptr);

  /** opposite of Append() */
  void Remove(double* value, GlobalFilter* f = nullptr);

  void ToInfo(PtrParamNode in) const;

  /** regular update to be called from Optimization::CommitIteration() to set value */
  void Update(unsigned int iter);

  Usage GetUsage() { return usage_; }

  /** the value we drive */
  double GetValue() const { return *value; }

  /** service to BaseOptimizer::Tuned -> based on inverse tanh projection, see paper.
   * Assumes always eta=.5 as the approximation via inverse is poor for eta=0 and 1. */
  double CalcTransistionZone(double kappa) const;

  static Enum<Method> method;

  static Enum<Usage> usage;

  /** we need the start value to properly scale for multiplicative grow with start < 1 */
  double start = -1;

  // <tune method="obj/mult/add" start="1" end="256" grow="1e-4" obj_max_grow="0.2" stride="1" stopping_greyness="true" />
private:

  /** evaluates the grayness rule, false also if not set */
  bool SufficientlyGray();

  /** helper for constructor */
  void FindGraynessStoppingRule();

  /** here we store the actual values.  */
  double* value = nullptr;

  /** here we store additional values (other filters in the robust case)
   * @see Append() and Remove() */
  StdVector<double*> external;

  double end = -1;

  /** prevent too early stopping by graynesss */
  static constexpr double OFF = -4711;
  double minimal = OFF;

  /** how often do we update? 1 is every time, 50 is every 50 iterations */
  unsigned int stride = 1;

  /** for all methods. Default depending on method */
  double grow = -1;

  /** only for OBJ */
  double max_grow_rate = 0.2;

  /** OBJ and ADD are multiplicative, this fails for values < 1. Therefore we conditionally scale the value in a range [1, end] */
  bool one_scale = false;

  /** when we want to stop for grayness, this is the corresponding rule */
  StoppingRule* grayness = nullptr;

  /** here we store the xml attribute for Register() */
  bool stopping_greyness_ = false;

  /** when we once stopped, we don't continue growing and do not further check the stopping critera */
  bool once_stopped_ = false;

  Optimization* opt = nullptr;

  /** when we are linked to a GlobalFilter we call SetNonLinCorrection() */
  StdVector<GlobalFilter*> gf;

  Method method_ = NO_METHOD;

  Usage usage_ = NO_USAGE;
};

}



#endif /* OPTIMIZATION_TUNE_HH_ */
