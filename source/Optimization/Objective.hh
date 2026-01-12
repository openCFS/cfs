#ifndef OBJECTIVE_HH_
#define OBJECTIVE_HH_

#include <string>
#include <tuple>
#include "boost/math/special_functions/sign.hpp"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/Enum.hh"
#include "MatVec/Vector.hh"
#include "Optimization/Function.hh"
#include "Optimization/Tune.hh"
#include "Utils/StdVector.hh"

namespace CoupledField {
class DesignSpace;
class DesignStructure;
class ConditionContainer;

}  // namespace CoupledField

namespace CoupledField
{
class MultipleExcitation;
class ObjectiveContainer;

/** gathered by some of the costFunction attributes in XML, the defaults are in the XML-Schema */
class StoppingRule
{
public:

  StoppingRule() {};

  /** @param pn might be NULL */
  void Init(PtrParamNode pn);

  typedef enum { DESIGN_CHANGE, REL_COST_CHANGE, ABOVE_FUNCTION, BELOW_FUNCTION, MAX_HOURS, OSCILLATIONS } Type;

  /** a single sufficient rule will break. Without sufficient, all rules need to be necessary the same time to break. */
  typedef enum { SUFFICIENT, NECESSARY } Condition;

  static Enum<Type> type;

  static Enum<Condition> condition;

  /** test for some criteria if we shall stop
   * @return if empty string, no. Otherwise the string is the msg */
  string DoStop(ObjectiveContainer* oc, ConditionContainer* cc, StdVector<double>& time);

  Type GetType() const { return type_; }

  Condition GetCondition() const { return condition_; }

  /** stopping rules value */
  double value = -1;

  /** stopping rule queue length */
  int queue = 1;

  /** function name for ABOVE_FUNCTION/ABOVE_FUNCTION to compare with objective/constraint info.xml name */
  std::string function;
private:

  Type type_ = DESIGN_CHANGE;

  Condition condition_ = SUFFICIENT;
};


/** We combine the cost function in a set to handle multiple of it.
 * It contains static const elements (and  working stuff).
 * MultipleExciation is in the XML file part of the objective but in class part of Optimization.
 * Note that ObjectiveContainer is a friend and has access to private data!*/
class Objective : public Function
{
  friend class ObjectiveContainer;

  public:

    /** Reads and processes the costFunction element but not the type information.
     * @param pn the costFunction pointer
     * @param pn_type either the costFunction pointer (for type attribute)
     *                or a multiObjective/objective element
     * @param index our position within the objectives for the design element costGradient */
    Objective(PtrParamNode pn, PtrParamNode pn_type, unsigned int index);
    
    /** Construct dummy function such that we can calculate the volume */
    Objective(Type type, double parameter = 0.0, Function::Access acc = FILTERED);

    virtual ~Objective() {};

    /** overwrites Function::IsObjective() */
    bool IsObjective() const { return true; }

    /** The name eventually enriched by the coord information for HOM_TENSOR */
    std::string GetName() const;

    /** Adds the value but don't touch scale */
    void AddValue(double add) { value_ += add; }

    /** Resets the Value to 0.0 */ // TODO move!
    void ResetValue() { value_ = 0.0; }

    /** is a homogenization tensor coord set */
    bool HasHomogenizationEntry() const { return std::get<0>(coord) != -1; }

    /** overloads Function::ToInfo() */
    void ToInfo(PtrParamNode info);

    /** for a given function value calculate the panalty method value max(0, value - parameter)^2 */
    double CalcPenalty(double value) const
    {
      double tmp = std::max(0.0, value - GetParameter());
      return tmp * tmp;
    }

    /** ln1p is ln damping wich is good for values < 0 and negative values */
    double CalcLn1p(double value) const
    {
       return boost::math::sign(value) * std::log(1.0 + std::abs(value));
    }

    /** penalty method = scale * max(x-parameter)^2
     * ln1p method = scale * sign(x) * ln(1 + |x|) good for x < 1 (even negative)
     * power method = scale * x^parameter where parameter=0.5 or 1/3 standard damping */
    typedef enum { LINEAR, PENALTY, LN1P, POWER } Term;

    static Enum<Term> term;

    Term GetTerm() const { return term_; }

    /** This defines the optional coord pair for HOM_TRACKING, HOM_FROBENIUS_PRODUCT.
     *  e.g. (1,1) for tensor entry (0,0). For Condition this is a list! The double shall be by default 1.0 */
    std::tuple<int, int, double> coord;

    /** Here we store our ParamNode such we can more easily access it in ErsatzMaterial */
    PtrParamNode pn;

    /** optionally we can tune (in the multi objective case) the scaling via tune */
    Tune tune;

    /** we add scale * term(function) */
    double scale = 1.0;

  private:  

    Term term_ = LINEAR;
};


/** This is a service container for the Objectives. We are friend of Objective and therefore have access to
 * private data! */
class ObjectiveContainer
{
public:
  ObjectiveContainer();

  ~ObjectiveContainer();

  /** Calls Function::PostProc() */
  void PostProc(DesignSpace* space, DesignStructure* structure, MultipleExcitation* me);

  /** actual constructor */
  void Read(PtrParamNode obj_node);

  /** Check for given objective */
  bool Has(Objective::Type type) const;

  Objective* Get(Objective::Type type, bool throw_exception = true);
  Objective* Get(const std::string& name, bool throw_exception = true);

  /** Sums up the history results of multiple objectives.
   * @param if true the scale value is included
   * @param history the positive or negative index within history. Negative in Python style (-1 last, ...)*/
  double GetHistoryValue(bool scale = true, int history = -1);

  unsigned int GetHistorySize();

  /** current values go to history */
  void PushBackHistory();

  /** calculate the distance of the current design to the last one and keep the current design */
  void PushBackDesign(const DesignSpace* space);

  /** Calls Objective::ToInfo() */
  void ToInfo(PtrParamNode in);

  bool DoMinimize() const { return minimize_; }
  bool DoMaximize() const { return !minimize_; }

  /** the history is calculated, the design change is simply stored! */
  StdVector<double> design_change;

  /** This are our stopping rule parameters */
  StdVector<StoppingRule> stop;

  StdVector<Objective*> data;

private:

  /** here we keep the last design, exported from DesignSpace - must stay uninitialized! */
  Vector<double> last_iteration_;

  /** this is the last design id -> for validation! */
  int last_design_;

  /** The task is the direction of a cost function (MINIMIZE, MAXIMIZE) */
  bool minimize_;
};

} // namespace


#endif /* OBJECTIVE_HH_*/
