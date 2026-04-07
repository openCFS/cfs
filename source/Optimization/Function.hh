#ifndef FUNCTION_HH_
#define FUNCTION_HH_

#include <assert.h>
#include <string>
#include <utility>
#include <tuple>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Design/DesignElement.hh"
#include "Design/DesignMaterial.hh"
#include "Design/DesignStructure.hh"
#include "Driver/FormsContexts.hh"
#include "General/Enum.hh"
#include "General/Environment.hh"
#include "MatVec/Matrix.hh"
#include "MatVec/Vector.hh"
#include "Utils/ApproxData.hh"
#include "Utils/StdVector.hh"
#include "Optimization/EigenInfo.hh"

namespace CoupledField {
class ErsatzMaterial;
}  // namespace CoupledField


using std::pair;
//using std::get;


namespace CoupledField
{
class Condition;
class DesignSpace;
class Excitation;
class MultipleExcitation;
class Objective;
class ShapeDesign;
class Context;

/** A Function is the (abstract) base class of Objective and Condition (which is a constraint but the name was
 * already used)
 */
class Function
{
    public:

    /** Dummy constructor for StdVector */
    Function();
    
    /** A Function is too stupid to do any useful - it is just a common base to avoid code dupliciy
     * @param pn our own element */
    Function(PtrParamNode pn);

    /** virtual dtor because base class */
    virtual ~Function();

    /** once we won't have this difference any more */
    static Function* Cast(Objective* c, Condition* g);

    /** PostProc called be the containers */
    virtual void PostProc(DesignSpace* space, DesignStructure* structure, ErsatzMaterial* = NULL);

    /** Different function types - some only objective, some only constraint some both */
    typedef enum {
      NO_TYPE = -1,
      // This are exclusive objective functions
      MULTI_OBJECTIVE,           /*!< Special type, not to be evaluated but trigger only */
      SLACK,                     /*!< for min max problems like min alpha s.th. compliance smaller alpha. Not really a function but triggers AuxDesign instead of DesignSpace. */
      SLACK_FNCT,                /*!< Indicates a formula of alpha and slack given in the attribute "function" and type of SlackFunciton */
      BANDGAP,                   /*!< bloch mode eigenfrequency band gap maximization. Requires gap element with the two eigenmode-numbers*/

      // This is objective and constraint together
      OUTPUT,                    /*!< Re(u,l) maximize solution where vector l is not 0 */
      SQUARED_OUTPUT,            /*!< Re(u,l)^2 maximize solution where vector l is not 0 */
      DYNAMIC_OUTPUT,            /*!< (u, L conj(u)) as OUTPUT but complex */
      REFLECTED_WAVE,            /*!< ((u-z), L conj(u-z)) as DYNAMIC_OUTPUT, but a complex value z is substracted. 
                                      E.g. when u is the total pressure and z is the pressure resulting from a normalVelocity b.c., u-z is the pressure of the reflected wave */
      ABS_OUTPUT,                /*!< |<u,l>| harmonic is implemented, real valued easy to add */
      CONJUGATE_COMPLIANCE,      /*!< (u, F conj(u)) as DYNAMIC_OUTPUT with trace of L is f */
      GLOBAL_DYNAMIC_COMPLIANCE, /*!< (u, I conj(u)) as DYNAMIC_OUTPUT with L is I (everywhere) */
      ELEC_ENERGY,               /*!< p^T K_pp p or p^T K_pp p^* */
      ENERGY_FLUX,               /*!< Re{j*w*u^T L grad_n u^*} */
      COMPLIANCE,                /*!< (u,f) the opposite of stiffness */
      VOLUME,                    /*!< normalized sum of original (but filtered) design elements */
      PENALIZED_VOLUME,          /*!< normalized sum of design elements penalized by parameter */
      GAP,                       /*!< The stingl constraint: volume - penalized volume */
      TRACKING,
      HOM_TENSOR,                /*!< optimize for the given coordinate if coord is set or console print of tensor.
                                      For a constraint it might blow up to several HOM_TENSOR if a tensor is given */
      HOM_TRACKING,              /*!< match a given tensor by L2 norm  */
      HOM_FROBENIUS_PRODUCT,     /*!< The Frobenius inner product for a given tensor sum_ij E^H_ij*D_ij. From Michael. The idea is
                                      maximize the inner mech energy <S, E^H S> for strains from a macro-problem. D = S*S^T */
      POISSONS_RATIO,            /*!< Poisson's Ration (\nu) within homogenization */
      YOUNGS_MODULUS,            /*!< Young's Modulus (E) within homogenization */
      YOUNGS_MODULUS_E1,         /*!< Young's Modulus (E1) within orthotrope homogenization */
      YOUNGS_MODULUS_E2,         /*!< Young's Modulus (E2) within orthotrope homogenization */
      TYCHONOFF,                 /*!< int(|| design ||^2) is a regularization form material opt. */
      TEMPERATURE,               /*!< for optimization of Poisson and heat conduction pde */
      HEAT_ENERGY,               /*!< for optimization in heat conduction pde, equivalent to compliance in linear elasticity */
      SQR_MAG_FLUX_DENS_X,       /*!< for optimization in squared magnetics Bx component */
      SQR_MAG_FLUX_DENS_Y,       /*!< for optimization in squared magnetics By component */
      SQR_MAG_FLUX_DENS_RZ,      /*!< for optimization in squared magnetics Br and Bz component */
      LOSS_MAG_FLUX_RZ,          /*!< for optimization in squared magnetics Br and Bz component, scaled by density*volume */
      MAG_COUPLING,              /*!< for optimization of inductive components */
      TEMP_TRACKING_AT_INTERFACE,/*!< tracking temperature at interfaces between solid and void elements */
      GLOBAL_SLOPE,              /*!< different implementation from local slopes */
      GLOBAL_MOLE,               /*!< see mole */
      GLOBAL_OSCILLATION,        /*!< see oscillation */
      GLOBAL_JUMP,
      GLOBAL_CURVATURE,
      GLOBAL_DESIGN,
      PERIMETER,                 /*!< perimeter constraint is a globalization of the (not meaningful local perimeter) */
      GLOBAL_STRESS,             /*!< global stress constraint: Kocvara and Stingl; 2007. Has adjoint! */
      EIGENFREQUENCY,            /*!< with the attribute ev for the number of the eigenfrequency/ eigenvalue */
      BUCKLING_LOAD_FACTOR,      /*!< with the attribute ev for the number of the load factor/ eigenvalue */
      LOCAL_BUCKLING_LOAD_FACTOR,/*!< microscopic load factor/ eigenvalue for two scale optimization*/
      GLOBAL_BUCKLING_LOAD_FACTOR,/*!< globalized microscopic load factor/ eigenvalue for two scale optimization*/
      ARC_OVERLAP,               /*!< prevents overlapping arc segments for spaghetti optimization */
      PYTHON_FUNCTION,           /*!< python global function */
      LOCAL_PYTHON_FUNCTION,     /*!< python local function */

      // External Solvers
      PRESSURE_DROP,             /*!< LBM Pressure Drop */

      // This is constraint only!
      GREYNESS,                  /*!< inaccurate - best for observation only */
      FILTERING_GAP,             /* !< sum of the difference between filtered and non-filtered design */
      REALVOLUME,
      ISOTROPY,                  /*!< blow up to several HOM_TENSOR constraints with different coords */
      ISO_ORTHOTROPY,            /*!< relaxed form of isotropy without fixing shear moduli */
      ORTHOTROPY,                /*!< just some 0 constraints */
      SLOPE,                     /*!< Implementation of a grad rho constraint */
      LOCAL_STRESS,              /*!< local constraint <sigma, M sigma> with stress transfer function */
      MOLE,                      /*!< Feature size control from T. Poulsen */
      OSCILLATION,               /*!< Feature size control by Fabian W. :) */
      JUMP,                      /*!< Weak greyness control by Fabian W. :) */
      BUMP,                      /*!< Prevent intermediate change of slope ('hobbala'). Multiplies slope with prev and with next */
      CURVATURE,                 /*!< Second derivative (prev, this, next) timing h=1 */
      PERIODIC,                  /*!< local constraint right minus left, meant for shape mapping */
      OVERHANG_VERT,             /*!< Overhang constraint for vertical (dof=x) shape mapping structures for additive manufacturing.  */
      OVERHANG_HOR,              /*!< Overhang constraint for horizontal (dof=y) shape mapping structures for additive manufacturing */
      DISTANCE,                  /*!< distance between start and end nodes of a noode (spaghetti design) */
      BENDING,                   /*!< noodle curvature divided by distance between start and end nodes (spaghetti design) */
      CONES,                     /*!< Cone constraints for injectivity of spline box in feature mapping */
      DESIGN_TRACKING,           /*!< Tracking against physical densities in designTarget. Either for region or periodic (constraint nodes) elements */
      SUM_MODULI,                /*!< the sum of the elasticity and shear moduli in parametrized elasticity tensor formulations */
      GLOBAL_SUM_MODULI,         /*!< global resource constraint, see sum_moduli */
      ORTHOTROPIC_TENSOR_TRACE,  /*!< tensor trace in the (DENSITY_TIMES_)ORTHOTROPIC parametrizations */
      GLOBAL_ORTHOTROPIC_TENSOR_TRACE, /*!< global constraint on the tensor trace in (DENSITY_TIMES_)ORTHOTROPIC parametrizations */
      TENSOR_TRACE,              /*!< local constraint on the tensor trace for FMO, laminates, hom_rect. Elasticity or dielec */
      TENSOR_NORM,               /*!< local squared L2 norm of the tensor coefficients (sum of the squared coefficients). For piezo-coupling in piezo FMP */
      TWO_SCALE_VOL,             /*!< Volume constraint / cost function for laminates and hom_rect parametrization */
      GLOBAL_TWO_SCALE_VOL,      /*!< global volume constraint / cost function for laminates and hom_rect parametrization */
      GLOBAL_TENSOR_TRACE,       /*!< global constraint on the tensor trace for fmo or laminates */
      PARAM_PS_POS_DEF,          /*!< constraint to ensure positive definiteness in parametrized elasticity tensor formulation (plane stress). Choose > 0*/
      POS_DEF_DET_MINOR_1,       /*!< 1st minor constraint for FMO positive definiteness by positive determinants */
      POS_DEF_DET_MINOR_2,       /*!< 2nd minor constraint for FMO positive definiteness by positive determinants */
      POS_DEF_DET_MINOR_3,       /*!< 3rd minor constraint for FMO positive definiteness by positive determinants */
      BENSON_VANDERBEI_1,        /*!< 1st minor constraint for numerical problemantic FMO pos def constraint */
      BENSON_VANDERBEI_2,        /*!< 2st minor constraint for numerical problemantic FMO pos def constraint */
      BENSON_VANDERBEI_3,        /*!< 3st minor constraint for numerical problemantic FMO pos def constraint */
      DESIGN,                    /*!< local design bound */
      MULTIMATERIAL_SUM,         /*!< local sum of multimaterial designs */
      SHAPE_INF,                 /*!< In Shape Optimization, there might be restrictions (not only box constraints) for shape parameters, this is the inf-norm version which splits nicely */
      EXPRESSION                 /*!< e.g. value smaller alpha+/-slack to be extended via mathparser when needed */
    } Type; // in ConditionContainer::VirtualView::Refresh() we assume a maximal value for the type. Check!!

    /** to convert string/enum for this type */
    static Enum<Type> type;

    /** See ToString() for string conversion! */
    Type GetType() const { return type_; }

    void SetType(Type type) {type_ = type;}

    typedef enum {
      NO_FUNCTION,               /*!< indicates we have not SLACK_FNCT function Type */
      ALPHA_SLACK_QUOTIENT,      /*!< quotient of the two slack variables alpha and slack: a/s */
      REL_BANDGAP,               /*!< relative band gap formulation 2*s/(a-s) based on 'alpha+/-slack' eigenfrequency bounds */
      NORM_BANDGAP,              /*!< normalized band gap (2*s/a) */
      ALPHA_MINUS_SLACK          /*!< slack variable distance a-s */
    } SlackFnct;

    static Enum<SlackFnct> slackFnct;

    SlackFnct GetSlackFnct() const { return slackFnct_; }

    typedef enum {
      WEIGHTED_SUM,             /*!< weighted sum of all objectives */
      SMOOTH_MIN,               /*!< smooth minimum of all objectives */
      SMOOTH_MAX                /*!< smooth maximum of all objectives */
    } MultiObjType;

    static Enum<MultiObjType> multiObjType;

    /** The real label might be an extended type string. E.g. by "access_".
     * Check if better use this than type.ToString(GetType()).
     * Is overloaded in Condition
     * @param me is for Condition */
    virtual std::string ToString() const;

    /** for historical reasons there are Condition and Objective pointers used concurrently. This is a
     * little helper. asserts that only of function is set. */
    static Function* GetFunction(Objective* f, Condition* g);

    /** are we objective or condition/constraint */
    virtual bool IsObjective() const = 0;

    /** Get the parameter, if it was set */
    double GetParameter() const { return parameter_;  }

    /** The evaluated function to process value. In the multiple objective case this contains weight and term handling! -1.0 if not set.
     * @see GetOrgValue() */
    virtual double GetValue() const { return value_; }

    /** overloaded in LocalCondition */
    virtual void SetValue(double val) { value_ = val; }

    /** the real function value without scaling and term handling for the multiple objective case (e.g. penalty max(0, v - parameter)^2. */
    double GetOrgValue() const { return org_value_; }
    void   SetOrgValue(double v) { org_value_ = v; }


    /** Access of the design variable in (local) functions. Both sensitivity and density.
     * There is no filtered sensitivity value, only gradient.
     * PLAIN: design variable, FILTERED: filtered design variable, PHYSICAL: filtered and penalized design variable
     * DEFAULT: set in DefaultAccess() and shall not be default at runtime.
     * FIXME: Usage is only for very few functions checked for consistent usage. Eeasily the gradient does handling but not the function */
    typedef enum { NO_ACCESS = -1, PLAIN, FILTERED, PHYSICAL, DEFAULT } Access;

    /** to convert string/enum for this type */
    static Enum<Access> access;

    /** This shall actually not be default but replaced by DefaultAccess() */
    Access GetAccess() const { return access_; }

    /** Is the access the programmatic default access (Function::DefaultAccess()) or manually changed in xml?
     * good hint for labeling output. (volume -> plain_volume, grayness -> physical_graness, ...) */
    bool IsDefaultAccess() const { return access_ == DefaultAccess(type_); }

    /** Some functions can have a physical counterpart. Which means e.g. for volume or greyness
     * the design variable with applied transfer function - hence as the FEM/physics sees the design.
     * One could call this penalized but physical is more exact and includes also density filtering.
     * The label becomes the appendix physical and the calculations are by interpolated values. */
    bool IsPhysical() const { return (access_ == PHYSICAL); }

    /** Is filtered but not necessary a filter is activated. Check filter_ for sensitivity/density if interested */
    bool isFiltered() const { return access_ == FILTERED || access_ == PHYSICAL; }

    /** Shall harmonic optimization multiply with omega^2.
    * This makes "u L conj(u)" to actually calc "v L conj(v)" with v = du/dt. -> approximatates sound intensity */
    bool FactorOmegaOmega() const { return omega_omega_; }

    /** The number of the eigenvalue (mode), one based!
     * @return 0 if it does not apply (1-based!) */
    unsigned int GetEigenValueID() { return eigenvalue_id_; }

    /** Shall/must we evaluate this objective at this excitation?
     * Sets the attribute excite_
     * Stress constraints in homogenization are triggered for a single constraint only.
     * @param excite_index -2 is uninitialized/auto, -1 is always */
    void SetExcitation(MultipleExcitation* me, int excite_index = -2);

    /** Get at least one excitation which applies to this function. For excite_ == -1 this might be one sample */
    Excitation* GetExcitation() { return sample_excitation_; }
    const Excitation* GetExcitation() const { return sample_excitation_; }

    /** Evaluate at this excitation? */
    bool DoEvaluate(const Excitation* excite) const;

    /** Evaluate for all excitations if there are multiple?
     * If we would so (excite_ == -1) we do it only for the sequence.
     * Never true for different sequence*/
    bool DoEvaluateAlways(int sequence) const;

    /** Is this function state dependent? */
    bool IsStateDependent() const;

    /** Requires this function an adjoint solution for the gradient? */
    bool IsAdjointBased() const;

    /** Are we generally excitation sensitive? E.g. stress */
    bool IsExcitationSensitive() const { return excite_sensitive_; }

    /** is this a slack type function ? */
    bool IsSlackFunction() const;

    /** checks the history - which contains only 0 for local functions */
    int CountOscillations() const;

    /** to be overwritten in Condition */
    virtual bool HasDenseJacobian() const { return true; }

    /** Is this a local function type */
    static bool IsLocal(Type t);

    bool IsLocal() const { return IsLocal(this->type_); }

    /** The relative position within this local constraints
     * @return -1 if it does not apply or >= 0 for local functions */
    virtual int GetCurrentRelativePosition() const { return -1; }

    /** Gives the sparsity pattern of the Jacobian. It gives the sorted, 0-based indices which have
     * values. For the dens case this is 0, 1, ... m.
     * This works only after ConditionContainer::PostProc() is called as otherwise the design is not known yet.
     * Is overwritten for local constraint which actually has spare patterns. */
    virtual StdVector<unsigned int>& GetSparsityPattern();

    /** LocalCondition::GetSparsityPattern() is always evaluated, this allows to a speed up */
    virtual unsigned int GetSparsityPatternSize() const { return jac_sparsity_.GetSize(); }

    /** Gives the sparsity pattern of the Hessian. Only for special nonlinear local functions */
    virtual Matrix<unsigned int>& GetHessianSparsityPattern();

    /** FeasPP can use the original functions for its strictly feasible MMA approximation. Makes only sense for
     * some special non-linear local constraints.
     * @param out needs to be the size of rows of GetHessianSparsityPattern()
     * @param factor -1 for normalizing lower bound constraints to c <= 0 */
    virtual void CalcHessian(StdVector<double>& out, double factor);

    /** Requires the function evaluation an selection vector associated to the adjoint RHS?.
     * Is an important for the solution of the state problem, if partial stuff from the adjoint setup is required. */
    bool NeedsSelectionVector() const;

    /** Requires an objective homogenization */
    bool IsHomogenization() const;

    /** Is this a linear function? E.g. SnOpt can handle them more efficiently */
    bool IsLinear() const { return linear_; }

    /** Snopt allows to set lower and upper bounds for functions. We make use of it for linear sparse abs() functions like SLOPE or CURVATURE.
     * Then we also need a NEXT or PREV_NEXT locality. When we do not make use of the double bounds (as for scpip or evaluate) we need a
     * NEXT_REVERSE or PREV_NEXT_REVERSE locality. Here we indicate if the function is meant for double bound. This is not the case for
     * vertical OVERHANG constraints (but for horizontal). Technically this is not restricted to local sparse functions but we do not consider it
     * for general functions.
     * True means, that the function is meant for double bounds AND there is no reverse locality. false if not local, other function, ... */
    bool IsDoubleBounded() const;

    /** is the local function intended for double bounded when the optimizer (snopt) supports it. Use with Local::IsReverse() */
    static bool CouldDoubleBounded(Type type);

    /** the tensor exists only in the homogenization constraint case */
    Matrix<double>& GetTensor() { return tensor_; }
    
    /** index within all objectives for design element gradient */
    int GetIndex() const { return index_; }

    /** Grayness scaling for a normalized variables is 4 in the non-physical case as 4 * 0.5 * (1-0.5) = 1 */
    double CalcGraynessScaling() const;

    /** Read the tensor if it is given, otherwise sets to 1.1
     * @param f_ctxt we call this during the constructor an therefore cannot use Function::ctxt
     * @param pn might contain a "tensor" child
     * @param matrix where to store the data
     * @return true if the tensor was read */
    static bool ReadTensor(Context* f_ctxt, PtrParamNode pn, Matrix<double>& matrix);

    /** @see StressConstraint::GetApplications */
    typedef enum { MECH, PIEZO, ONLY_COUPLING } StressType;

    static Enum<StressType> stressType;

    StressType GetStressType() { return stressType_; }

    /** for volume to check the notation in the FMO case with tensor_trace design. */
    MaterialTensorNotation GetNotation() const { return notation_; }

    /** for the bandgap function. Could clearly be a general gap between two functions. This could then handle
     * the old gap function from Michael (volume - penalized volume) */
    struct BandGap
    {
      BandGap() {lower_ev = -1; upper_ev = -1; }
      int lower_ev;
      int upper_ev;

      EigenInfo lower;
      EigenInfo upper;
    };

    BandGap bandgap;


    /** A function can be a local function when it is calculated by the local neighborhood state.
     * This does NOT mean, that the function may not be a global function, e.g. when a the L2 norm
     * of the local information is used!
     * Due to the (current) separation of Objective and Condition the local function property is not
     * a derived "is local property" but it "has the local property". Very similar to the DesignElement
     * structure. */
    class Local
    {
    public:

      /** Initialized the neighborhood */
      Local(Function* func, DesignSpace* space);

      ~Local();

      /** second constructor step. Required as locality initialization calls ShapeMapDesign() where we ask for periodic
       * and the constructor needs to be finished then */
      void PostInit();

      /** Number of identifiers per design element. Usually dim or dim *2, ... */
      int GetElementDimension() const { return element_dimension_; }

      /** The local type, essentially important for slopes. There should be no need to set
       * it as user. */
      typedef enum {
        DEFAULT,                 /*!< Function::PostProc() finds proper value */
        NEXT,                    /*!< x_i and x_i+1 */              /* x, x+, y+ and xy +*/
        NEXT_AND_REVERSE,        /*!< x_i and x_i+1 plus x_i+1 PLUS the x_i for classical slope */
        NEXT_DIAG,
        PREV_NEXT,
        PREV_NEXT_AND_REVERSE,   /*!< x_i-1 and x_i+1 with different sign for small oscillation */
        DEG_45_STAR,             /*!< Different notation. prev_next but also diagonals */
        DEG_45_STAR_AND_REVERSE, /*!< The doubled variant of DEG_45_STAR for oscillation */
        BOUNDARY,                /*!< For a neighbor definition the first and last element (JUMP) */
        CYCLIC,                  /*!< The periodic element for the the periodic constraint */
        ELEMENT,                 /*!< For stress there is no neighborhood, only the element itself */
        MULT_DESIGNS_ELEMENT,    /*!< ELEMENT for multiple different designs - only parametrized PLANE_STRESS for now */
        SHAPE,                   /*!< SHAPE, the sparsity pattern is read from file */
        MULT_DESIGNS_NEXT,                    /*!< x_i and x_i+1 */
        MULT_DESIGNS_NEXT_AND_REVERSE,        /*!< x_i and x_i+1 plus x_i+1 PLUS the x_i for classical slope */
        MULT_DESIGNS_PREV_NEXT,
        MULT_DESIGNS_PREV_NEXT_AND_REVERSE,  /*!< ELEMENT for multiple different designs - only parametrized PLANE_STRESS for now */
        FUNCTION_SPECIFIC,        /*!< tuned for single complex functions: distance, bending */
        FUNCTION_SPECIFIC_TWO_SIGNS,  /*!< reverse case */
        EXTERNALLY_DEFINED        /* for local python function we get the full neighborbood from there */
      } Locality;

      static Enum<Locality> locality;

      Locality GetLocality() const { return locality_; }

      static bool IsReverse(Locality loc);

      static bool RequiresEps(Function::Type type);
      static bool RequiresBeta(Function::Type type);

      /** The phase for oscillation constraint only to define two constraints with different
       * feature sizes for material and void */
      typedef enum {
        BOTH = -1000,   // sync the values with the NO_SIGN, VOID_SIGN and VOID_MATERIAL constants
        VOID_MAT = -1,
        MATERIAL = 1
      } Phase;

      static Enum<Phase> phase;

      /** are we periodic? Only if enabled in local/periodic AND with periodic pde */
      bool periodic;

      /** Data structure for the interpolation coefficients for latticeVol3D*/
      Matrix<double> vol_coeff_;
      Vector<double> vol_a_;
      Vector<double> vol_b_;
      Vector<double> vol_c_;

      ApproxData* volumeInterpolator_;

      /** total volume for CalcLaminatesVol in the unregular grid case*/
      double total_vol_;

      Phase GetPhase() const { return phase_; }

      /** The beta value for smoothing min/max, checks if its set. */
      double GetBeta() const { assert(beta_ != -3.14); return beta_; }

      /** The eps value form smoothing abs, checks if its set. */
      double GetEps() const { assert(eps_ >= 0); return eps_; }

      /** The power for the globalized local sum( max(0, g_i(x) - c)^p) */
      double GetPower() const { return power_; }

      /** is globalied local? */
      bool IsGlobalized() const { return globalized_; }

      /** normalize globalSlope, globalCheckerboard, ... by number of local constraints */
      bool DoNormalizeGlobal() const { return normalize_; }

      /** This is to be called from Local::Local() as Function::ToInfo() is called too early */
      void ToInfo(PtrParamNode info);

      /** the mapping from a relative slope constraint number (0-based) to the actual
       * constraint. This allows to remove constraints for elements which have no (full)
       * neighborhood */
      struct Identifier
      {
      public:
        const static int NO_SIGN;
        /** alias for sign == -1 to VOID Phase for oscillation only */
        const static int VOID_SIGN;
        /** alias for sign == -1 to MATERIAL Phase for oscillation only */
        const static int MATERIAL_SIGN;

        /** default constructor for StdVector() */
        Identifier() {};

        /** @param prev if NONE neighbor is size 1 otherwise size two */
        Identifier(BaseDesignElement* elem, BaseDesignElement* prev, BaseDesignElement* next, int si = NO_SIGN);

        /** Identifier when we have a neighborhood defined by a radius - eg mole */
        Identifier(BaseDesignElement* elem, StdVector<BaseDesignElement*> buddies, int si = NO_SIGN);

        /** spline box specific variant (cones) */
        Identifier(BaseDesignElement* elem, StdVector<BaseDesignElement*> buddies, StdVector<BaseDesignElement*> sb_buddies, int si = NO_SIGN);

        /** Returns the element
         * @param idx == -1 for elem, otherwise from neighbors */
        BaseDesignElement* GetElement(int idx) {
          return idx == -1 ? element : neighbor[idx];
        }

        const BaseDesignElement* GetElement(int idx) const {
          return const_cast<const BaseDesignElement*>(idx == -1 ? element : neighbor[idx]);
        }

        /** identifies the element by the design type. Works only for special neighborhoods! */
        const BaseDesignElement* GetElementByType(BaseDesignElement::Type type) const;

        /** returns design value by the design type.
         * @param get_parameter == true works for ParamMat parameters. Works only for special neighborhoods! */
        double GetDesign(BaseDesignElement::Type type, const Local* local, const DesignElement::Access access = DesignElement::SMART, bool get_parameter = true) const;

        /** Calculates the local function based on function->type.
         * In the *_STRESS case uses GetValue(), Calls SetValue().
         * Is very fast for grad_glob and power == 1
         * @param grad_glob only active when globalized. If grad_glob is active EvalFunction is called as in EvalGrad in order to calculated
         * the gradient. */
        double EvalFunction(const Local* local, bool grad_glob = false);

        /** Service function. Calculates all gradients for this and the neighbors. Only for real local function!.
         * Note, that the von Mises Stress gradient is NOT calculated here but in SIMP::CalcVonMisesStressGradient()!
         * It does the proper constraint_gradient reset first! */
        void EvalGradient(const Local* local);

        /** calculates the slope identified by this neighbor. When sign is not set assumes sign=1.
         * "Petersson, Sigmund; Slope Constrained Topology Optimization; 1998" */
        double CalcSlope() const;
        /** calculate the slope gradient for a given element
         * @param neigh_idx for -1 for the own element, otherwise the neighbor */
        double CalcSlopeGradient(int neigh_idx) const;

        /** calculate the overhang constraint for shape mapping variables for use in additive manufacturing */
        double CalcOverhang(Function::Type ft, double eps) const;
        double CalcOverhangGradient(int neigh_idx, Function::Type ft, double eps) const;

        /** distance of start and end nodes of a noodle for spaghetti optimization */
        double CalcDistance(int neigh_idx, bool grad) const;

        /** noodle curvature divided by distance of start and end node for spaghetti optimization */
        double CalcBending(int neigh_idx, bool grad) const;

        /** implemented in PythonFunction.cc */
        double CalcLocalPythonFunc(const Local* local) const;
        /** we return the "full" local python function gradient at once. Implemented in PythonFunction.cc */
        void CalcLocalPythonGrad(Vector<double>& grad, const Local* local) const;

        /** calculate the cone constraint for spline box variables */
        double CalcCones(const Local* local) const;
        double CalcConesGradient(int neigh_idx, const Local* local) const;

        /** calculates the design bound as constraint. */
        double CalcDesignBound(Function* f, const Local* l, bool derivative) const;


        /** the perimeter is similar to the slope constraint but always globalized (sum) */
        double CalcPerimeter(double eps, double l_k) const;
        double CalcPerimeterGradient(int neigh_idx, double eps, double l_k) const;

        /** periodic means right end minus left end. Meant for shape mapping. A relaxed equal needs both bounds */
        double CalcPeriodic() const;
        double CalcPeriodicGradient(int neigh_idx) const;

        /** calculates the checkerboard value. The sign determines if the smaller or larger value is evaluated
         * @param beta < 0 is real max, otherwise it is a Kreiselmeier Steinhauser approximation */
        double CalcCheckerboard(double beta) const;
        double CalcCheckerboardGradient(int neigh_idx, double beta);

        /** T. Poulsen's feature size control */
        double CalcMole(double eps) const;
        double CalcMoleGradient(int neigh_idx, double eps);

        /** Oscillation is a feature size control which is by variable radius a generalization of a checkerboard constraint */
        double CalcOscillation(double beta) const;
        double CalcOscillationGradient(int neigh_idx, double beta);

        /** weak formulation of a greyness control */
        double CalcJump() const;
        double CalcJumpGradient(int neigh_id) const;

        /** no change of slope sign. Positive if the prev and next slope have the same sign (getting larger or smaller) */
        double CalcBump() const;
        double CalcBumpGradient(int neigh_idx) const;

        /** Curvature Calculation (simplest second derivative) assuming h=1 */
        double CalcCurvature() const;
        double CalcCurvatureGradient(int neigh_idx) const;

        /** sum of elasticity and shear moduli in parametrized elasticity tensor formulations */
        double CalcSumModuli(const Local* local, DesignElement::Access access = DesignElement::PLAIN, int neigh_idx = -1, bool derivative = false) const;

        /** tensor trace of the material tensor in (DENSITY_TIMES_)ORTHOTROPIC parametrizations */
        double CalcOrthotropicTensorTrace(const Local* local, DesignElement::Access access = DesignElement::PLAIN, int neigh_idx = -1, bool derivative = false) const;

        /** volume of material (strong phase for plane strain) in laminate homogenization and two_scale formulas */
        double CalcTwoScaleVolume(const Local* local, DesignElement::Access access = DesignElement::PLAIN, int neigh_idx = -1, bool derivative = false) const;

        /** volume of material from homogenized lattice structure in 3D */
        double CalcLatticeVolume3D(const Local* local, DesignElement::Access access, int neigh_idx=-1, bool derivative = false) const;

        /** to ensure positive definiteness of the material tensor E3-E1*nu31^2 > 0 has to hold */
        double CalcParamPSPosDef(const Local* local, DesignElement::Access access, int neigh_idx, bool derivative) const;

        /** local tensor trace for FMO */
        double CalcTensorTrace(int neigh_idx, const Local* local, bool derivative) const;

        /** squared L2 norm of all tensor entries. Meant for the piezoelectric coupling tensor */
        double CalcTensorNorm(int neigh_idx, const Local* local, bool derivative) const;

        /** sum of all multimaterial designs */
        double CalcMultiMaterialSum(int neigh_idx, const Local* local, bool derivative) const;

        /** local FMO positive definiteness of (E-val*I) >= param via determinants
         * @param type the type we want to evaluate. Might be different from local->func->type_ in Approximation::TransformMultiplyer() */
        double CalcPosDefDeterminant(int neigh_idx, const Local* local, bool derivative, Type type) const;

        /** local FMO positive definiteness of (E-val*I) >= param via Benson Vanderbei constraints */
        double CalcBensonVanderbei(int neigh_idx, const Local* local, bool derivative, Type type) const;

        /* @param type the type we want to evaluate. Might be different from local->func->type_ in Approximation::TransformMultiplyer() */
        //double CalcPosDefDeterminant(int neigh_idx, const Local* local, bool derivative, Type type) const;

        /** CalcStress() and the gradient are actually done in EM/SIMP */
        
        /** Shape Parameter Constraints */
        double CalcShape(Function* f, const Local* l) const;
        double CalcShapeGradient(Function* f, const Local* l, int neigh_idx) const;

        /** debug output */
        std::string ToString() const;

        BaseDesignElement* element = nullptr; // this represents DesignSpace::data[element_idx]

        /** this are all design elements for the local function. For slopes a spatial neighborhood, for
         * globalLaminatesVolumes the variables stiff1, stiff2 for the same FE-Element, ...
         * @see GetElement() */
        StdVector<BaseDesignElement*> neighbor;

        /** a spline box neighbor - specific for cones */
        StdVector<BaseDesignElement*> sb_neighbor;

        /** sign is only needed if we treat slope constraints as two separate constraints.
         *  in case we do not do this, sign will be -1000, else -1 for X_N, 1 for X_P.
         *  for spline box design with cones constraints, sign will be between -12 and 12
         *  and indicate which equation of linear cone to use */
        int sign = NO_SIGN;

        StdVector<int> signs;

        /** specific to bending for spaghetti with variable Normal: 0 N 0, 0 N N, N N 0, N N N */
        typedef enum { NO_BENDING, ZNZ, ZNN, NNZ, NNN } BendingType;

        BendingType bending = NO_BENDING;

      private:
        /** to be reused */
        static Vector<double> tmp1;
        static Vector<double> tmp2;
      }; // end of struct Identifier


      /** Elements with no full neighborhood are not stored. If they would be stored
       * we could easily calculate the virtual element number.
       * This vector maps from the relative virtual constraint number (0 based)
       * to the relative element index (also 0-based). If we have all periodic b.c.
       * then this is a 1:1 mapping, otherwise this list is smaller than 2*dim*n by
       * 2*dim<not full neighborhood> */
      StdVector<Identifier> virtual_elem_map;

      /** We need the design space to access the values */
      DesignSpace* space = nullptr;

      /** Store the local values. At least used for LOCAL_STRESS. See LocalCondition::Get/SetValue() */
      StdVector<double> local_values;

    private:

      /** Service method for the constructor
       * @param phase see SetupStarLocalityElementMap() */
      void SetupVirtualElementMap(Phase phase = BOTH);

      void SetupVirtualStarLocalElementMap(const Function*); //only used in NEXT_DIAG

      /** Special implementation for DEG_45_STAR[_AND_REVERSE] locality.
       * @param phase for oscillation we can separate void and material which is the sign convention */
      void SetupStarLocalityElementMap(Phase phase = BOTH);

      /** for a defined neighborhood only the most prev and next element, not this element */
      void SetupBoundaryElementMap();

      /** trival case form ELEMENT (stress) -> on the element itself */
      void SetupSingularElementMap();

      /** multiple designs on one element for paramMat
       * @param for FMO_POS_DEF_* we need to know which minor */
      void SetupMultDesignsElementMap(const Function* f = NULL);
      
      /** Shape Element Maps */
      void SetupShapeElementMap(const Function* func, ShapeDesign* design);

      /** Multiple designs on several elements for paramMat*/
      void SetupMultDesignsVirtualElementMap(const Function* f = NULL);

      /** small helper to determine the number of neighbors in each (diagonal)
       * direction if we use a neighborhood. Parses the whole stuff */
      struct NeighborhoodStructure
      {
      public:
        NeighborhoodStructure(Local* local, PtrParamNode local_node);

        /** number of element in the X_P, Y_P or Z_P direction. Total length is (2* X_P + 1) */
        StdVector<unsigned int> orthogonal;

        /** the diagonal size in elements. For 2D only one diagonal */
        StdVector<unsigned int> diagonal;

        /** writes the stuff. E.g. in a sub-element */
        void ToInfo(PtrParamNode info);

      private:
        double radius = -1.0;
        double value  = -1.0;
        Filter::FilterSpace fs = Filter::NO_FILTER;
      }; // end of struct NeighborhoodStructure


      NeighborhoodStructure* structure_ = nullptr;

      Function* func_ = nullptr;

      Locality locality_ = DEFAULT;

      /** Functions based on a relaxed max formulation have beta for the Kreiselmeier/Steinhauser
       * continuation. This is (global) checkerboard. -1 is real max = infinity */
      double beta_ = -1.0;

      /** relaxation parameter to smooth abs by A(x) = sqrt(x^2 + eps^2) - eps. */
      double eps_ = -1.0;

      /** power for globalization */
      double power_ = -1.0;

      /** For oscillation we can define if we want the constraint for void, material or both.
       * Such different feature sized can be defined */
      Phase phase_ = BOTH;

      /** normalize global function */
      bool normalize_ = false;

      /** are we a local condition or globalized (condition or objective) */
      bool globalized_ = false;

      /** @see GetElementDimension() */
      int element_dimension_ = 0;
    };


    /** Give the local information. Check for NULL */
    Local* GetLocal() { return local; }

    /** The design type is by default DEFAULT :) */
    BaseDesignElement::Type GetDesignType() const {return design_; }

    /** This are the elements the Function is defined on. Either references to the
     * elements within the design space or to dummy elements if the region is not within the design (stress)
     * @param region as long as only the Condition has this stuff it is an parameter*/
    void SetElements(DesignSpace* space, RegionIdType region);

    /** We also store here the info ptr. When overload, call also this. */
    virtual void ToInfo(PtrParamNode info);

    /** for the python function get_opt_function_properties() */
    virtual void DescribeProperties(StdVector<std::pair<string,string> >& map) const;

    /** to export the (global) python functions to ErsatzMaterial.
     * The local functions are defined in Function::Local::Identifier::CalcLocalPython*()
     * @param eval if false the grad is called */
    PyObject* CallPythonFunction(bool eval);

    /** Here we store our ParamNode such we can more easily access it in ErsatzMaterial */
    PtrParamNode pn;

    /** If condition supports restriction to one region. Currently ALL_REGIONS for objectives */
    RegionIdType region;

    /** real or pseudo design elements defined by the region.
     * if region is ALL_REGIONS this points to the standard design space.
     * Otherwise it is a sub set pointing to the design space or if region is not within the design (stress constraint)
     * it is filled by DesignElements with dummy values.
     * Created on request */
    StdVector<DesignElement*> elements;

    /** When we optimize output we store here the rhs loads */
    StdVector<LinearFormContext*> output_forms;

    /** the multiple sequence step we belong to.
     * @see ContextManager */
    Context* ctxt;

    /** This vector stores the function values by iteration. Written in *Container::PushBackHistory() */
    Vector<double> history;

  protected:

    /** common constructor stuff. To be called from special Objective constructor, too */
    void Init();

    /** Is reentrant save. Initialize the local variable
     * @return either a new Local or the old one */
    Local* InitLocal(DesignSpace* space);
    
    /** extract the "coord" element and parse it to coord */
    static void ParseCoord(PtrParamNode pn, std::tuple<int, int, double>& coord);

    /** By the size of DesignSpace::GetNumberOfVariables() which might include slack - to be handled in AuxDesign.
     * the sparse patterns are determined on the fly by LocalCondition::GetSparsityPattern() */
    void SetDenseSparsityPattern(DesignSpace* space);

    /** Implementation in PythonFunction.cc. Needs to be called late, when we have DesignSpace which might contain the script
     * @param pn of full function which python as child (to be checked) */
    void InitPythonFunction(PtrParamNode pn, DesignSpace* space);

    /** obtain from python the sparsity pattern */
    void SetLocalPythonVirtualElementMap(StdVector<Function::Local::Identifier>& virtual_element_map, DesignSpace* space);

    /** matrices for polynomial coefficients and discretization steps of the interpolation for volume calculation in 3D with cross shaped base cells*/

    /** This is DEFAULT (= applies always) if not defined */
    BaseDesignElement::Type design_ = BaseDesignElement::NO_TYPE;

    /** The actual kind of cost function. */
    Type type_ = NO_TYPE;

    /** The slack function type */
    SlackFnct slackFnct_ = NO_FUNCTION;

    /** for HOM_TRACKING this is the target tensor. For HOM_FROBENIUS_PRODUCT this is the parameter */
    Matrix<double> tensor_;

    /** The current function value */
    double value_ = 0.0;

    /** for the multiple objecitve case where we set to value_ something complex in Optimization::CalcObjective() */
    double org_value_ = -1.0;

    /** Some special functions use a parameter: slope constraint and penalized volume */
    double parameter_;

    /** considered by some function: wheather plain value of the design variable, the filtered or the physical value is used */
    Access access_ = NO_ACCESS;

    /** this index is the position in the Optimization list and is used to
     * identify the constraint gradient in DesignElement. Only relevant for type = active */
    int index_;

    /** (sample) excitation. For excite_ -1 this is only an exemplaric excitation */
    Excitation* sample_excitation_;

    /** Is this function excitation sensitive? */
    int excite_sensitive_;

    /** @see FactorOmegaOmega() */
    bool omega_omega_;

    /** the "ev" parameter for the eigenvalue function. 1-based!
     * @see Condition::bloch_extremal_ */
    int eigenvalue_id_ = -1;

    /** Conditions mark themselves as (non) linear -> no power in the design variable, ...*/
    bool linear_;

    /** Do we have local information? E.G. (global) slopes */
    Local* local;

    /** Here we store our info node. Set only by ToInfo() *after* PostProc. Use preInfo_ instead. Is null when not set, yet. */
    PtrParamNode info_;

    /** In case info_ is not set, yet. ToInfo applies this information. Is initialized in the constructor */
    PtrParamNode preInfo_;

    StressType stressType_;

    /** the sparsity pattern of the Jacobian to be set by ConditionContainer::PostProc() via SetSparsity() */
    StdVector<unsigned int> jac_sparsity_;

    /** the sparsity pattern of the Hessian for special local functions to be used by FeasPP.
     * @see jac_sparsity_
     * @see Approximation::hess_pattern */
    Matrix<unsigned int> hess_sparsity_;

    /** only for tensor trace and volume */
    MaterialTensorNotation notation_;

    /** All python functions have a single argument an dict with all xml values as flat list */
    /** the given options enriched by all xml function values */
    PyObject* py_arg_ = NULL;

    /** this is the name the python function gives itself. Default is python (schema) */
    std::string py_name_;

    /** the python function to evaluate. Having it set indicates, we are a python function */
    PyObject* py_eval_ = NULL;
    PyObject* py_grad_ = NULL;
    PyObject* py_sparsity_ = NULL;
    /** here we store the gradient vector for local python grad evaluations */
    CfsTLS<Vector<double> > py_local_grad_;

  private:

    /** destructor for python part, implemented in PythonFunction.cc */
    void DeletePythonFunction();

    /** to replace default access in xml with a real value */
    Access DefaultAccess(Type ft) const;

    /** special value for excite_ value.
    * -1 for all excitations within this sequence!!. Most interesting for stress constraints.
    * -2 is for unset!
    * -3 for excitation "0_1" */
    typedef enum { ALL_EX = -1, UNSET_EX = -2, COMBINED_0_1_EX = -3 } ExciteIndex;

    /** Excitation index for evaluation.
     * Note that the index is unique over all sequences!
     * >= 0 for the actual excitation, ExciteIndex for other cases */
    int excite_ = -100;

    /** We need this as argument for DefaultAccess() when called by IsDefaultAccess()
     * as we usually have no space and don't want to store it (whyever). */
    Filter::Type filter_ = Filter::NO_FILTERING;
};


} // namespace


#endif /* FUNCTION_HH_ */
