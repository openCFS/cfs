#ifndef FUNCTION_HH_
#define FUNCTION_HH_

#include "General/Enum.hh"
#include "MatVec/matrix.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/DesignStructure.hh"


using std::pair;


namespace CoupledField
{
class Condition;
class Objective;
class DesignSpace;

/** A Function is the (abstract) base class of Objective and Condition (which is a constraint but the name was
 * already used)
 */
class Function
{
    public:

    /** Dummy constructor for StdVector */
    Function() {};
    
    /** virtual dtor because base class */
    virtual ~Function();

    /** A Function is too stupid to do any useful - it is just a common base to avoid code dupliciy
     * @param pn our own element */
    Function(PtrParamNode pn);

    /** PostProc called be the containers */
    virtual void PostProc(DesignSpace* space, DesignStructure* structure);

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
      ENERGY_FLUX,               /*!< Re{j*w*u^T L grad_n u^*} */

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
      GLOBAL_SLOPE,              /*!< different implementation from local slopes */
      GLOBAL_MOLE,               /*!< see mole */
      GLOBAL_OSCILLATION,        /*!< see oscillation */
      GLOBAL_JUMP,

      // This is constraint only!
      GREYNESS,                  /*!< inaccurate - best for observation only */
      REALVOLUME,
      ISOTROPY,                  /*!< blow up to several HOMOGENITATION_TENSOR constraints with different coords */
      SLOPE,                     /*!< Implementation of a grad rho constraint */
      MOLE,                      /*!< Feature size control from T. Poulsen */
      OSCILLATION,               /*!< Feature size control by Fabian W. :) */
      JUMP                       /*!< Weak greyness control by Fabian W. :) */
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

    /** are we objective of condition/constraint */
    virtual bool IsObjective() const = 0;

    /** Get the parameter, if it was set */
    double GetParameter() const { return parameter_;  }

    /** The evaluates function values. -1.0 if not set. */
    virtual double GetValue() const { return value_; }

    /** overloaded in SlopeCondition */
    virtual void SetValue(double val) { value_ = val; }

    /** Some functions can have a physical counterpart. Which means e.g. for volume or greyness
     * the design variable with applied transfer function - hence as the FEM/physics sees the design.
     * One could call this penalized but physical is more exact and includes also density filtering.
     * The label becomes the appendix physical and the calculations are by interpolated values. */
    bool IsPhysical() const { return physical_; }

    /** Shall harmonic optimization multiply with omega^2.
    * This makes "u L conj(u)" to actually calc "v L conj(v)" with v = du/dt. -> approximatates sound intensity */
    bool FactorOmegaOmega() const { return omega_omega_; }

    /** Shall/must we evaluate this objective only of the last excitation? */
    bool DoEvaluateOnce() const { return evaluateOnce_; }

    /** Requires an objective homogenization */
    bool IsHomogenization() const;

    /** Is this function sensitive for density filtering?
     * This implies the gradient with post differentiation.
     * The function does *NOT* know if design filtering is enabled!
     * Volume is always sensitive - physical is for penalization! */
    bool ForDensityFiltering() const;

    /** If we would do sensitivity filtering, do it also for this function?
     * @see ForDensityFiltering() */
    bool ForSensitivityFiltering() const;

    /** the tensor exists only in the homogenization constraint case */
    Matrix<double>& GetTensor() { return tensor_; }
    
    /** index within all objectives for design element gradient */
    int GetIndex() const { return index_; }

    /** Read the tensor if it is given, otherwise sets to 1.1
     * @param pn might contain a "tensor" child
     * @param matrix where to store the data
     * @return true if the tensor was read */
    static bool ReadTensor(PtrParamNode pn, Matrix<double>& matrix);


    /** A function can be be a local function when it is calculated by the local neighborhood state.
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

      /** Number of identifiers per design element. Usually dim or dim *2, ... */
      int GetElememtDimension() const { return element_dimension_; }

      /** The local type, essentially important for slopes. There should be no need to set
       * it as user. */
      typedef enum {
        DEFAULT,                 /*!< Function::PostProc() finds proper value */
        NEXT,                    /*!< x_i and x_i+1 */
        NEXT_AND_REVERSE,        /*!< x_i and x_i+1 plus x_i+1 PLUS the x_i for classical slope */
        PREV_NEXT,
        PREV_NEXT_AND_REVERSE,   /*!< x_i-1 and x_i+1 with different sign for small oscillation */
        DEG_45_STAR,             /*!< Different notation. prev_next but also diagonals */
        DEG_45_STAR_AND_REVERSE, /*!< The doubled variant of DEG_45_STAR for oscillation */
        BOUNDARY                 /*!< For a neighbor definition the first and last element (JUMP) */
      } Locality;

      static Enum<Locality> locality;

      Locality GetLocality() const { return locality_; }

      /** The pase for oscillation constraint only to define two constraints with different
       * feature sizes for material and void */
      typedef enum {
        BOTH = -1000, // syn the values with the NO_SIGN, VOID_SIGN and VOID_MATERIAL constants
        VOID = -1,
        MATERIAL = 1
      } Phase;

      static Enum<Phase> phase;

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
        Identifier() : sign(NO_SIGN) {}

        /** @param prev if NONE neighbor is size 1 otherwise size two */
        Identifier(DesignElement* elem, DesignElement* prev, DesignElement* next, int si = NO_SIGN);

        /** Identifier when we have a neighborgood defined by a radius - eg mole */
        Identifier(DesignElement* elem, StdVector<DesignElement*> buddies, int si = NO_SIGN);

        /** Returns the element
         * @param idx == -1 for elem, otherwise form neighbors */
        DesignElement* GetElement(int idx) {
          return idx == -1 ? element : neighbor[idx];
        }

        const DesignElement* GetElement(int idx) const {
          return const_cast<const DesignElement*>(idx == -1 ? element : neighbor[idx]);
        }

        /** Service function. Calculates the actual objective, based on function->type */
        double EvalFunction(const Local* local) const;

        /** Service function. Calculates all gradients for this and the neighbors. Only for real local function!.
         * It does the proper constraint_gradient reset first! */
        void EvalGradient(const Local* local);

        /** calculates the slope identified by this neighbor. When sign is not set assumes sign=1.
         * "Petersson, Sigmund; Slope Constrained Topology Optimization; 1998" */
        double CalcSlope() const;

        /** calculate the slope gradient for a given element
         * @param neigh_idx for -1 for the own element, otherwise the neighbor */
        double CalcSlopeGradient(int neigh_idx) const;

        /** calculates the checkerboard value. The sign determines if the smaller or larger value is evaluated
         * @param beta < 0 is real max, otherwise it is a Kreiselmeier Steinhauser approximation */
        double CalcCheckerboard(double beta) const;

        /** calculates the gradient for the checkerbord
         * @see CalcSlopeGradient() */
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

        DesignElement* element; // this represents DesignSpace::data[element_idx]
        StdVector<DesignElement*> neighbor;

        /** sign is only needed if we treat slope constraints as two separate constraints
         *  in case we do not do this, sign will be -1000, else -1 for X_N, 1 for X_P */
        int sign;
      private:
        /** to be reused */
        static StdVector<double> tmp1;
        static StdVector<double> tmp2;
      };

      /** Elements with no full neighborhood are not stored. If they would be stored
       * we could easily calculate the virtual element number.
       * This vector maps from the relative virtual constraint number (0 based)
       * to the relative element index (also 0-based). If we have all periodic b.c.
       * then this is a 1:1 mapping, otherwise this list is smaller than 2*dim*n by
       * 2*dim<not full neighborhood> */
      StdVector<Identifier> virtual_elem_map;

      /** We need the design space to access the values */
      DesignSpace* space;

      /** Store the local values. */
      StdVector<double> values;

    private:
      /** Service method for the constructor
       * @param phase see SetupStarLocalityElementMap() */
      void SetupVirtualElementMap(Phase phase = BOTH);

      /** Special implementation for DEG_45_STAR[_AND_REVERSE] locality.
       * @param phase for oscillation we can separate void and material which is the sign convention */
      void SetupStarLocalityElementMap(Phase phase = BOTH);

      /** for a defined neighborhood only the most prev and next element, not this element */
      void SetupBoundaryElementMap();

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
        double radius;
        double value;
        DesignStructure::FilterSpace fs;
      };

      NeighborhoodStructure* structure_;

      Function* func_;

      Locality locality_;

      /** Functions based on a relaxed max formulation have beta for the Kreiselmeier/Steinhauser
       * continuation. This is (global) checkerboard. -1 is real max = infinity */
      double beta_;

      /** relaxation parameter to smooth abs by A(x) = sqrt(x^2 + eps^2) - eps. For (global) mole only */
      double eps_;

      /** power for globalization */
      double power_;

      /** For oscillation we can define if we want the constraint for void, material or both.
       * Such different feature sized can be defined */
      Phase phase_;

      /** normalize global function */
      bool normalize_;

      /** are we a local condition or globalized (condition or objective) */
      bool globalized_;

      /** @see GetElementDimension() */
      int element_dimension_;
    };

    /** Give the local information. Check for NULL */
    Local* GetLocal() { return local; }

    /** We also store here the info ptr. When overload, call also this. */
    virtual void ToInfo(PtrParamNode info);

    /** Here we store our ParamNode such we can more easily access it in ErsatzMaterial */
    PtrParamNode pn;

  protected:

    /** Is reentrant save. Initialize the local variable
     * @return either a new Local or the old one */
    Local* InitLocal(DesignSpace* space);
    
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

    /** @see FactorOmegaOmega() */
    bool omega_omega_;

    bool harmonic_;

    /** Do we have local information? E.G. (global) slopes */
    Local* local;

    /** Here we store our info node */
    PtrParamNode info_;

};


} // namespace


#endif /* FUNCTION_HH_ */
