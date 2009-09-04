#ifndef OBJECTIVE_HH_
#define OBJECTIVE_HH_

#include "Utils/StdVector.hh"
#include "MatVec/matrix.hh"

namespace CoupledField
{

class ParamNode;
class InfoNode;
class ObjectiveContainer;

/** We combine the cost function in a set to handle multiple of it.
 * It contains static const elements (and  working stuff).
 * MultipleExciation is in the XML file part of the objective but in class part of Optimization.
 * Note that ObjectiveContainer is a friend and has access to private data!*/
class Objective
{
  friend class ObjectiveContainer;

  public:
    /** Reads and processes the costFunction element but not the type information.
     * @param pn the costFunction pointer
     * @param pn_type either the costFunction pointer (for type attribute)
     *                or a multiObjective/objective element
     * @param index our position withn the objectives for the design element costGradient */
    Objective(ParamNode* pn, ParamNode* pn_type, unsigned int index);

    /** Known types of cost functions */
    typedef enum {
      MULTI_OBJECTIVE,           /*!< Special type, not to be evaluated but trigger only */
      COMPLIANCE,                /*!< (u,f) the opposite of stiffness */
      OUTPUT,                    /*!< Re(u,l) maximize solution where vector l is not 0 */
      DYNAMIC_OUTPUT,            /*!< (u, L conj(u)) as OUTPUT but complex */
      ABS_DYN_OUTPUT_SQUARED,    /*!< |<u,l>|^2 = | sum u_l |^2 = < sum u_l, sum u_l> harmonic */
      CONJUGATE_COMPLIANCE,      /*!< (u, F conj(u)) as DYNAMIC_OUTPUT with trace of L is f */
      GLOBAL_DYNAMIC_COMPLIANCE, /*!< (u, I conj(u)) as DYNAMIC_OUTPUT with L is I (everywhere) */
      ELEC_ENERGY,               /*!< p^T K_pp p or p^T K_pp p^* */
      VOLUME,
      TRACKING,
      HOMOGENIZATION_TENSOR,     /*!< match a (partly) given tensor element wise  */
      HOMOGENIZATION_TRACKING,   /*!< match a given tensor by L2 norm  */
      TYCHONOFF,                 /*!< int(|| design ||^2) is a regularization form material opt. */
      TEMPERATURE                /*!< for optimization of poisson and heat conduction pde */
    } Type;

    /** to convert string/enum for this type */
    static Enum<Type> type;

    Type GetType() const { return type_; }

    /** The the current objective value */
    double GetValue() const;

    /** @param val set an unscaled value here -> it is scaled in the return */
    void SetValue(double val);

    /** Adds the value but don't touch penalty */
    void AddValue(double add) { value_ += add; }

    /** Resets the Value to 0.0 */ // TODO move!
    void ResetValue() { value_ = 0.0; }


    double GetPenalty() const { return penalty_; }

    /** Shall harmonic optimization multiply with omega^2.
    * This makes "u L conj(u)" to actually calc "v L conj(v)" with v = du/dt. -> approximatates sound intensity */
    bool FactorOmegaOmega() const { return omega_omega_; }

    /** index within all objectives for design element gradient */
    unsigned int GetIndex() const { return index_; }

    /** Shall/must we evaluate this objective only of the last excitation? */
    bool DoEvaluateOnce() const { return evaluateOnce_; }

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

    /** Here we store our ParamNode such we can more easily access it in ErsatzMaterial */
    ParamNode* pn;

    /** for HOMOGENIZATION_TRACKING this is the target tensor. */
    Matrix<double> tensor;

  private:

    /** The actual kind of cost function. */
    Type type_;

    /** The current value */
    double value_;

    /** This vector stores the cost functions of the iterations. Written in GetObjective() */
    StdVector<double> history_;

    /** @see FactorOmegaOmega() */
    bool omega_omega_;

    /** Is this type only possible/necessary for the last excitation?
     * Then it is only in that case evaluated and the excitation weight is ignored */
    bool evaluateOnce_;

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
  void Read(ParamNode* obj_node);

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

  void ToInfo(InfoNode* in) const;

  bool DoMinimize() const { return minimize_; }
  bool DoMaximize() const { return !minimize_; }

  StdVector<Objective*> data;

private:

  /** The task is the direction of a cost function (MINIMIZE, MAXIMIZE) */
  bool minimize_;
};

} // namespace


#endif /* OBJECTIVE_HH_*/
