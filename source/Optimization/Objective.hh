#ifndef OBJECTIVE_HH_
#define OBJECTIVE_HH_

#include "Optimization/Function.hh"
#include "Utils/StdVector.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField
{
class ObjectiveContainer;
class MultipleExcitation;

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
     * @param index our position withn the objectives for the design element costGradient */
    Objective(PtrParamNode pn, PtrParamNode pn_type, unsigned int index);
    
    /** Construct dummy function such that we can calculate the volume */
    Objective(Type type, double parameter = 0.0, bool physical = false);

    virtual ~Objective() {};

    /** overwrites Function::IsObjective() */
    bool IsObjective() const { return true; }

    /** The name eventually enriched by the coord information for HOMOGENIZATION_TENSOR */
    std::string GetName() const;

    /** Adds the value but don't touch penalty */
    void AddValue(double add) { value_ += add; }

    /** Resets the Value to 0.0 */ // TODO move!
    void ResetValue() { value_ = 0.0; }

    /** is a homogenization tensor coord set */
    bool HasHomogenizationEntry() const { return get<0>(coord) != -1; }

    double GetPenalty() const { return penalty_; }

    /** overloads Function::ToInfo() */
    void ToInfo(PtrParamNode info);

    /** This defines the optional coord pair for HOMOGENIZATION_TRACKING.
     *  e.g. (1,1) for tensor entry (0,0). For Condition this is a list! The double shall be by default 1.0 */
    tuple<int, int, double> coord;

    /** Here we store our ParamNode such we can more easily access it in ErsatzMaterial */
    PtrParamNode pn;

  private:

    /** This vector stores the cost functions of the iterations. Written in GetObjective() */
    StdVector<double> history_;

    /** by default 1.0 if not multiObjective */
    double penalty_;
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

  /** Sums up the history results of multiple objectives.
   * @param costs all objecties
   * @param if true the penalty value is included
   * @param history the index within history or -1 for the last value */
  double GetHistoryValue(bool penalty = true, int history = -1);

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

  /** gathered by some of the costFunction attributes in XML, the defaults are in the XML-Schema */
  class StoppingRule
  {
  public:

    StoppingRule();

    /** @param pn might be NULL */
    void Init(PtrParamNode pn);

    typedef enum { DESIGN_CHANGE, REL_COST_CHANGE } Type;

    static Enum<Type> type;

    Type GetType() const { return type_; }

    /** stopping rules value */
    double value;

    /** stopping rule queue length */
    unsigned int queue;

  private:

    Type type_;
  };

  /** This are our stopping rule parameters */
  StoppingRule stop;

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
