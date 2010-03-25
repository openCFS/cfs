#ifndef FUNCTION_HH_
#define FUNCTION_HH_

#include "General/Enum.hh"
#include "MatVec/matrix.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

using std::pair;


namespace CoupledField
{
class Condition;
class Objective;

/** A Function is the (abstact) base class of Objective and Condition (which is a constraint but the name was
 * already used)
 */
class Function
{
    public:

    /** Dummy constructor for StdVector */
    Function() {};
    
    /** virtual dtor because base class */
    virtual ~Function() {};

    /** A Function is too stupid to do any useful - it is just a common base to avoid code dupliciy
     * @param pn our own element */
    Function(PtrParamNode pn);

    /** Different function types - some only objective, some only constraint some both */
    typedef enum {
      // This are exclusive objective functions
      MULTI_OBJECTIVE,           /*!< Special type, not to be evaluated but trigger only */
      OUTPUT,                    /*!< Re(u,l) maximize solution where vector l is not 0 */
      DYNAMIC_OUTPUT,            /*!< (u, L conj(u)) as OUTPUT but complex */
      ABS_DYN_OUTPUT_SQUARED,    /*!< |<u,l>|^2 = | sum u_l |^2 = < sum u_l, sum u_l> harmonic */
      CONJUGATE_COMPLIANCE,      /*!< (u, F conj(u)) as DYNAMIC_OUTPUT with trace of L is f */
      GLOBAL_DYNAMIC_COMPLIANCE, /*!< (u, I conj(u)) as DYNAMIC_OUTPUT with L is I (everywhere) */
      ELEC_ENERGY,               /*!< p^T K_pp p or p^T K_pp p^* */
      ACOU_NEAR_FIELD,           /*!< -Re{j*w*psi^T L grad_n psi^*} */

      // This is objective and constraint together
      COMPLIANCE,                /*!< (u,f) the opposite of stiffness */
      VOLUME,                    /*!< normalized sum of original design elements */
      PENALIZED_VOLUME,          /*!< normalized sum of design elements penalized by parameter */
      GAP,                       /*!< The stingl constraint: volume - penalized volume */
      TRACKING,
      HOMOGENIZATION_TENSOR,     /*!< optimize for the given coordinate if coord is set or console print of tensor.
                                      For a constraint it might blow up to several HOMOGENIZATION_TENSOR if a tensor is given */

      HOMOGENIZATION_TRACKING,   /*!< match a given tensor by L2 norm  */
      POISSONS_RATIO,            /*!< Poisson's Ration (\nu) within homogenization */
      YOUNGS_MODULUS,            /*!< Young's Modulus (E) within homogenization */
      TYCHONOFF,                 /*!< int(|| design ||^2) is a regularization form material opt. */
      TEMPERATURE,               /*!< for optimization of poisson and heat conduction pde */

      // This is constraint only!
      GREYNESS,                  /*!< inaccurate - best for observation only */
      REALVOLUME,
      ISOTROPY,                  /*!< blow up to several HOMOGENITATION_TENSOR constraints with different coords */
      SLOPE,                     /*!< Implementation of a grad rho constraint */
      CHECKERBOARD               /*!< Meassure for the checkerboard, up to now only observe! */
    } Type;

    /** to convert string/enum for this type */
    static Enum<Type> type;

    /** See ToString() for string conversion! */
    Type GetType() const { return type_; }

    /** The real label might be an extended type string. E.g. by "physical_".
     * Check if better use this than type.ToString(GetType()).
     * Is overloaded in Condition */
    virtual std::string ToString() const;

    /** for historical reasons there are Condition and Objective pointers used concurrently. This is a
     * little helper. asserts that only of function is set. */
    static Function* GetFunction(Objective* f, Condition* g);

    /** Get the parameter, if it was set */
    double GetParameter() const
    {
      assert(parameter_ != -1.0);
      return parameter_;
    }


    /** The evaluates function values. -1.0 if not set. */
    virtual double GetValue() const
    {
      return value_;
    }

    /** overloaded in SlopeCondition */
    virtual void SetValue(double val)
    {
      value_ = val;
    }

    /** Some functions can have a physical counterpart. Which means e.g. for volume or greyness
     * the design variable with applied transfer function - hence as the FEM/physics sees the design.
     * One could call this penalized but physical is more exact and includes also density filtering.
     * The label becomes the appendix physical and the calculations are by interpolated values. */
    bool IsPhysical() const
    {
      return physical_;
    }


    /** Shall/must we evaluate this objective only of the last excitation? */
    bool DoEvaluateOnce() const { return evaluateOnce_; }

    /** Requires an objective homogenization */
    bool IsHomogenization() const;

    /** the tensor exists only in the homogenization constraint case */
    Matrix<double>& GetTensor() { return tensor_; }
    
    /** index within all objectives for design element gradient */
    int GetIndex() const { return index_; }

    /** Here we store our ParamNode such we can more easily access it in ErsatzMaterial */
    PtrParamNode pn;

    /** Read the tensor if it is given, otherwise sets to 1.1
     * @param pn might contain a "tensor" child
     * @param matrix where to store the data
     * @return true if the tensor was read */
    static bool ReadTensor(PtrParamNode pn, Matrix<double>& matrix);

  protected:
    
    /** extract the "coord" element and parse it to coord */
    static void ParseCoord(PtrParamNode pn, std::pair<int, int>& coord);

    /** The actual kind of cost function. */
    Type type_;

    /** for HOMOGENIZATION_TRACKING this is the target tensor. */
    Matrix<double> tensor_;

    /** The current function value */
    double value_;

    /** Some special functions use a parameter: slope constraint and penalized volume */
    double parameter_;

    /** Is this type only possible/necessary for the last excitation?
     * Then it is only in that case evaluated and the excitation weight is ignored */
    bool evaluateOnce_;

    /** @see IsPhysical() */
    bool physical_;

    /** this index is the position in the Optimization list and is used to
     * identify the constraint gradient in DesignElement. Only relevant for type = active */
    int index_;
};

} // namespace


#endif /* FUNCTION_HH_ */
