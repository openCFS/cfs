#ifndef OBJECTIVE_HH_
#define OBJECTIVE_HH_

#include "Optimization/Function.hh"
#include "Utils/StdVector.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField
{
class ObjectiveContainer;

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
    
    virtual ~Objective() {};

    /** The name eventually enriched by the coord information for HOMOGENIZATION_TENSOR */
    std::string GetName() const;

    /** Adds the value but don't touch penalty */
    void AddValue(double add) { value_ += add; }

    /** Resets the Value to 0.0 */ // TODO move!
    void ResetValue() { value_ = 0.0; }

    /** is a homogenization tensor coord set */
    bool HasHomogenizationEntry() const { return coord.first != -1; }

    double GetPenalty() const { return penalty_; }

    /** Shall harmonic optimization multiply with omega^2.
    * This makes "u L conj(u)" to actually calc "v L conj(v)" with v = du/dt. -> approximatates sound intensity */
    bool FactorOmegaOmega() const { return omega_omega_; }

    /** index within all objectives for design element gradient */
    unsigned int GetIndex() const { return index_; }

    /** gathered by some of the costFunction attributes in XML, the defaults are in the XML-Schema */
    class StoppingRule
    {
    public:
      /** stopping rules value */
      double value;

      /** stopping rule queue length */
      unsigned int queue;
    };

    /** This are our stopping rule parameters */
    StoppingRule stop;

    /** This defines the optional coord pair for HOMOGENIZATION_TRACKING.
     *  e.g. (1,1) for tensor entry (0,0). For Condition this is a list! */
    std::pair<int, int> coord;

    /** Here we store our ParamNode such we can more easily access it in ErsatzMaterial */
    PtrParamNode pn;

  private:

    /** This vector stores the cost functions of the iterations. Written in GetObjective() */
    StdVector<double> history_;

    /** @see FactorOmegaOmega() */
    bool omega_omega_;

    bool harmonic_;

    /** by default 1.0 if not multiObjective */
    double penalty_;

    /** our index within the objectives for gradient storage in DesignElement */
    unsigned int index_;
};


/** This is a service container for the Objectives. We are friend of Objective and therefore have access to
 * private data! */
class ObjectiveContainer
{
public:
  ObjectiveContainer();

  ~ObjectiveContainer();

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

  void ToInfo(PtrParamNode in) const;

  bool DoMinimize() const { return minimize_; }
  bool DoMaximize() const { return !minimize_; }

  StdVector<Objective*> data;

private:

  /** The task is the direction of a cost function (MINIMIZE, MAXIMIZE) */
  bool minimize_;
};

} // namespace


#endif /* OBJECTIVE_HH_*/
