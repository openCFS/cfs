#ifndef FUNCTION_HH_
#define FUNCTION_HH_

#include <assert.h>
#include <stddef.h>
#include <string>
#include <utility>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/Enum.hh"
#include "General/defs.hh"
#include "General/Environment.hh"
#include "MatVec/Matrix.hh"
#include "MatVec/Vector.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/DesignStructure.hh"
#include "Optimization/Design/DesignMaterial.hh"
#include "Utils/StdVector.hh"
#include "boost/tuple/tuple.hpp"

namespace CoupledField {
class ErsatzMaterial;
}  // namespace CoupledField


using boost::tuple;
using std::pair;
using boost::get;


namespace CoupledField
{
class Condition;
class DesignSpace;
class Excitation;
class MultipleExcitation;
class Objective;
class ShapeDesign;

/** A Function is the (abstract) base class of Objective and Condition (which is a constraint but the name was
 * already used)
 */
class Function
{
    public:

    /** Dummy constructor for StdVector */
    Function() {};
    
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
      // This are exclusive objective functions
      MULTI_OBJECTIVE,           /*!< Special type, not to be evaluated but trigger only */
      OUTPUT,                    /*!< Re(u,l) maximize solution where vector l is not 0 */
      DYNAMIC_OUTPUT,            /*!< (u, L conj(u)) as OUTPUT but complex */
      ABS_OUTPUT,                /*!< |<u,l>| harmonic is implemented, real valued easy to add */
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
      TEMPERATURE,               /*!< for optimization of poisson and heat conduction pde */
      GLOBAL_SLOPE,              /*!< different implementation from local slopes */
      GLOBAL_MOLE,               /*!< see mole */
      GLOBAL_OSCILLATION,        /*!< see oscillation */
      GLOBAL_JUMP,
      STRESS,                    /*!< global stress constraint: Kocvara and Stingl; 2007. Has adjoint! */
      STRESS_DENSITY,            /*!< global stress divided by volume */
      PROJECTION,                /*!< Michael's idea: sum_i || nu(rho_i) - H_eta_beta(rho_i) ||^2 <= eps */

      // External Solvers
      PRESSURE_DROP,             /*!< LBM Pressure Drop */

      // This is constraint only!
      GREYNESS,                  /*!< inaccurate - best for observation only */
      REALVOLUME,
      ISOTROPY,                  /*!< blow up to several HOM_TENSOR constraints with different coords */
      ISO_ORTHOTROPY,            /*!< relaxed form of isotropy without fixing shear moduli */
      ORTHOTROPY,                /*!< just some 0 constraints */
      SLOPE,                     /*!< Implementation of a grad rho constraint */
      MOLE,                      /*!< Feature size control from T. Poulsen */
      OSCILLATION,               /*!< Feature size control by Fabian W. :) */
      JUMP,                      /*!< Weak greyness control by Fabian W. :) */
      BUMP,                      /*!< Prevent intermediate change of slope ('hobbala') by Fabian W. */
      DESIGN_TRACKING,           /*!< Tracking against physical densities in designTarget. Either for region or periodic (constraint nodes) elements */
      SUM_MODULI,                /*!< the sum of the elasticity and shear moduli in parametrized elasticity tensor formulations */
      GLOBAL_SUM_MODULI,         /*!< global resource constraint, see sum_moduli */
      ORTHOTROPIC_TENSOR_TRACE,  /*!< tensor trace in the (DENSITY_TIMES_)ORTHOTROPIC parametrizations */
      GLOBAL_ORTHOTROPIC_TENSOR_TRACE, /*!< global constraint on the tensor trace in (DENSITY_TIMES_)ORTHOTROPIC parametrizations */
      TENSOR_TRACE,              /*!< local constraint on the tensor trace for FMO, laminates, hom_rect. Elasticity or dielec */
      TENSOR_NORM,               /*!< local squared L2 norm of the tensor coefficients (sum of the squared coefficients). For piezo-coupling in piezo FMP */
      LAMINATES_VOL,             /*!< Volume constraint / cost function for laminates parametrization */
      GLOBAL_LAMINATES_VOL,      /*!< global volume constraint / cost function for laminates parametrization */
      GLOBAL_TENSOR_TRACE,       /*!< global constraint on the tensor trace for fmo or laminates */
      PARAM_PS_POS_DEF,          /*!< constraint to ensure positive definiteness in parametrized elasticity tensor formulation (plane stress). Choose > 0*/
      POS_DEF_DET_MINOR_1,       /*!< 1st minor constraint for FMO positive definiteness by positive determinants */
      POS_DEF_DET_MINOR_2,       /*!< 2nd minor constraint for FMO positive definiteness by positive determinants */
      POS_DEF_DET_MINOR_3,       /*!< 3rd minor constraint for FMO positive definiteness by positive determinants */
      BENSON_VANDERBEI_1,        /*!< 1st minor constraint for numerical problemantic FMO pos def constraint */
      BENSON_VANDERBEI_2,        /*!< 2st minor constraint for numerical problemantic FMO pos def constraint */
      BENSON_VANDERBEI_3,        /*!< 3st minor constraint for numerical problemantic FMO pos def constraint */
      DESIGN_BOUND,              /*!< local design bound */
      MULTIMATERIAL_SUM,         /*!< local sum of multimaterial designs */
      SLACK,                     /*!< for min max problems like min alpha s.th. compliance smaller alpha. Not really a function
                                      but triggers AuxDesign instead of DesignSpace. */
      SHAPE_INF                  /*!< In Shape Optimization, there might be restrictions (not only box constraints) for shape parameters, this is the inf-norm version which splits nicely */
    } Type; // in ConditionContainer::VirtualView::Refresh() we assume a maximal value for the type. Check!!

    /** to convert string/enum for this type */
    static Enum<Type> type;

    /** See ToString() for string conversion! */
    Type GetType() const { return type_; }


    /** The real label might be an extended type string. E.g. by "physical_".
     * Check if better use this than type.ToString(GetType()).
     * Is overloaded in Condition
     * @param me is for Condition */
    virtual std::string ToString(MultipleExcitation* me = NULL) const;

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

    /** Shall/must we evaluate this objective at this excitation?
     * Stress constraints in homogenization are triggered for a single constraint only.
     * @param excite_index -2 is uninitialized/auto, -1 is always */
    void SetExcitation(MultipleExcitation* me, int excite_index = -2);

    /** Evaluate at this excitation? */
    bool DoEvaluate(const Excitation* excite) const;

    /** Evaluate for all excitations if there are multiple? */
    bool DoEvaluateAlways() const;

    /** Are we generally excitation sensitive? E.g. stress */
    bool IsExcitationSensitive() const;

    /** to be overwritten in Condition */
    virtual bool HasDenseJacobian() const { return true; }

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

    /** Requires this function an adjoint solution for the gradient? */
    bool IsAdjointBased() const;

    /** Requires the function evaluation an selection vector associated to the adjoint RHS?.
     * Is an important for the solution of the state problem, if partial stuff from the adjoint setup is required. */
    bool NeedsSelectionVector() const;

    /** Requires an objective homogenization */
    bool IsHomogenization() const;

    /** Is this a linear function? E.g. SnOpt can handle them more efficiently */
    bool IsLinear() const { return linear_; }

    /** Is this a local function type */
    static bool IsLocal(Type t);

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

    /** @see StressConstraint::GetApplications */
    typedef enum { MECH, PIEZO, ONLY_COUPLING } StressType;

    static Enum<StressType> stressType;

    StressType GetStressType() { return stressType_; }

    /** for volume to check the notation in the FMO case with tensor_trace design. */
    DesignMaterial::Notation GetNotation() const { return notation_; }

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
      int GetElementDimension() const { return element_dimension_; }

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
        BOUNDARY,                /*!< For a neighbor definition the first and last element (JUMP) */
        ELEMENT,                 /*!< For stress there is no neighborhood, only the element itself */
        MULT_DESIGNS_ELEMENT,    /*!< ELEMENT for multiple different designs - only parametrized PLANE_STRESS for now */
        SHAPE                    /*!< SHAPE, the sparsity pattern is read from file */
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

      /** Data structure for the interpolation coefficients for latticeVol3D*/
      Matrix<double> vol_coeff_;
      Matrix<double> vol_a_;
      Matrix<double> vol_b_;
      Matrix<double> vol_c_;

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
        Identifier() : sign(NO_SIGN) {}

        /** @param prev if NONE neighbor is size 1 otherwise size two */
        Identifier(BaseDesignElement* elem, BaseDesignElement* prev, BaseDesignElement* next, int si = NO_SIGN);

        /** Identifier when we have a neighborgood defined by a radius - eg mole */
        Identifier(BaseDesignElement* elem, StdVector<BaseDesignElement*> buddies, int si = NO_SIGN);

        /** Returns the element
         * @param idx == -1 for elem, otherwise form neighbors */
        BaseDesignElement* GetElement(int idx) {
          return idx == -1 ? element : neighbor[idx];
        }

        const BaseDesignElement* GetElement(int idx) const {
          return const_cast<const BaseDesignElement*>(idx == -1 ? element : neighbor[idx]);
        }

        /** identifies the element by the design type. Works only for special neighborhoods! */
        BaseDesignElement* GetElement(DesignElement::Type type);

        /** Service function. Calculates the actual objective, based on function->type.
         * Is very fast for grad_glob and power == 1
         * @param grad_glob only active when globalized. If grad_glob is active EvalFunction is called as in EvalGrad in order to calculated
         * the gradient.
         * @param von_mises_stresss only used when the function is stress -> determined by ErsatzMaterial::CalcVonMisesStressGlobalizationFactor() */
        double EvalFunction(const Local* local, bool grad_glob = false, double von_mises_stresss = -1.0) const;

        /** Service function. Calculates all gradients for this and the neighbors. Only for real local function!.
         * Note, that the von Mises Stress gradient is NOT calculated here but in SIMP::CalcVonMisesStressGradient()!
         * It does the proper constraint_gradient reset first! */
        void EvalGradient(const Local* local);

        /** calculates the slope identified by this neighbor. When sign is not set assumes sign=1.
         * "Petersson, Sigmund; Slope Constrained Topology Optimization; 1998" */
        double CalcSlope() const;

        /** calculates the design bound as constraint. */
        double CalcDesignBound(bool derivative) const;

        /** calculate the slope gradient for a given element
         * @param neigh_idx for -1 for the own element, otherwise the neighbor */
        double CalcSlopeGradient(int neigh_idx) const;

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

        /** no change of slope sign */
        double CalcBump() const;
        double CalcBumpGradient(int neigh_idx) const;

        /** sum of elasticity and shear moduli in parametrized elasticity tensor formulations */
        double CalcSumModuli(DesignElement::Access access = DesignElement::PLAIN, int neigh_idx = -1, bool derivative = false) const;

        /** tensor trace of the material tensor in (DENSITY_TIMES_)ORTHOTROPIC parametrizations */
        double CalcOrthotropicTensorTrace(const Local* local, DesignElement::Access access = DesignElement::PLAIN, int neigh_idx = -1, bool derivative = false) const;

        /** volume of material of the homogenized cross shaped structure in 3D including derivatives */
        //double Calc3DCrossVolume(double stiff1, double stiff2, double stiff3, bool derivative, double der) const;

        /** Function returns/interpolates the volume in 3D for cross shaped base cell*/
        double Interpolate_Volume3D(Vector<double>& p, const Matrix<double>& vol_a, const Matrix<double>& vol_b, const Matrix<double>& vol_c, const Matrix<double>& vol_coeff,
            double direction) const;

        /** Function evaluates the interpolation polynomial used for volume calculation in 3D for cross shaped base cell*/
        double EvaluateC1Interpolation_3D( Vector<double>& p, const Matrix<double>& vol_a, const Matrix<double>& vol_b, const Matrix<double>& vol_c,const Matrix<double> & vol_coeff, double & da, double & db,
            double & dc, int & j, int & k, int & l, int & m, int & n, int &o) const;

        /** Function calculates the derivative of the interpolation polynomial with respect to stiffness number, specified by variable direction*/
        double EvaluateC1Interpolation_Deriv_3D(Vector<double>& p, const Matrix<double> & vol_a, const Matrix<double> & vol_b, const Matrix<double>& vol_c, const Matrix<double> & vol_coeff, double & da, double & db,
            double & dc, int & j, int & k, int & l, int & m, int & n, int & o,
            double direction) const;

        /** volume of material (strong phase for plane strain) in laminate homogenization formulas */
        double CalcLaminatesVolume(const Local* local, DesignElement::Access access = DesignElement::PLAIN, int neigh_idx = -1, bool derivative = false) const;

        /** volume of material from homogenized lattice structure in 3D */
        double CalcLatticeVolume3D(const Local* local, DesignElement::Access access, int neigh_idx=-1, bool derivative = false) const;


        /** to ensure positive definiteness of the material tensor E3-E1*nu31^2 > 0 has to hold */
        double CalcParamPSPosDef(int neigh_idx, bool derivative) const;

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

        /** CalcStress() and the gradient are actually done in EM/SIMP */
        
        /** Shape Parameter Constraints */
        double CalcShape(Function* f, const Local* l) const;
        double CalcShapeGradient(Function* f, const Local* l, int neigh_idx) const;

        BaseDesignElement* element; // this represents DesignSpace::data[element_idx]

        /** this are all design elements for the local function. For slopes a spatial neighborhood, for
         * globalLaminatesVolumes the variables stiff1, stiff2 for the same FE-Element, ...
         * @see GetElement() */
        StdVector<BaseDesignElement*> neighbor;

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
      Vector<double> values;

      /** Here ErsatzMaterial::CalcGlobalFunction() stores the number of infeasible element functions and prints the
       * value in Optimization::LogFileLine() -> just a service */
      int infeasible;

    private:
      /** Service method for the constructor
       * @param phase see SetupStarLocalityElementMap() */
      void SetupVirtualElementMap(Phase phase = BOTH);

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

    /** The design type is by default DEFAULT :) */
    DesignElement::Type GetDesignType() const {return design_; }

    /** Give the projection data */
    StdVector<DesignElement>& GetProjectionDesignClone();

    /** This are the elements the Function is defined on. Either references to the
     * elements within the design space to to dummy elements if the region is not within the design (stress)
     * @param region as long as only the Condition has this stuff it is an parameter*/
    void SetElements(DesignSpace* space, RegionIdType region);

    /** We also store here the info ptr. When overload, call also this. */
    virtual void ToInfo(PtrParamNode info);

    /** Here we store our ParamNode such we can more easily access it in ErsatzMaterial */
    PtrParamNode pn;

    /** If condition supports constriction to one region. Currently ALL_REGIONS for objectives */
    RegionIdType region;

    /** real or pseudo design elements defined by the region.
     * if region is ALL_REGIONS this points to the standard design space.
     * Otherwise it is a sub set pointing to the design space or if region is not within the design (stress constraint)
     * it is filled by DesignElements with dummy values.
     * Created on request */
    StdVector<DesignElement*> elements;

    /** When we optimize output we store here the nodes */
    // FIXME LoadList output_nodes;

  protected:

    /** common constructor stuff. To be called from special Objective constructor, too */
    void Init();

    /** Is reentrant save. Initialize the local variable
     * @return either a new Local or the old one */
    Local* InitLocal(DesignSpace* space);
    
    /** extract the "coord" element and parse it to coord */
    static void ParseCoord(PtrParamNode pn, tuple<int, int, double>& coord);

    /** By the size of DesignSpace::GetNumberOfVariables() which might include slack - to be handled in AuxDesign.
     * the sparse patterns are determined on the fly by LocalCondition::GetSparsityPattern() */
    void SetDenseSparsityPattern(DesignSpace* space);

    /** matrices for polynomial coefficients and discretization steps of the interpolation for volume calculation in 3D with cross shaped base cells*/


    /** This is DEFAULT (= applies always) if not defined */
    DesignElement::Type design_;

    /** The actual kind of cost function. */
    Type type_;

    /** for HOM_TRACKING this is the target tensor. For HOM_FROBENIUS_PRODUCT this is the parameter */
    Matrix<double> tensor_;

    /** The current function value */
    double value_;


    /** Some special functions use a parameter: slope constraint and penalized volume */
    double parameter_;


    /** @see IsPhysical() */
    bool physical_;

    /** this index is the position in the Optimization list and is used to
     * identify the constraint gradient in DesignElement. Only relevant for type = active */
    int index_;

    /** Excitation index for evaluation. -1 for all excitations. Most interesting for stress constraints.
     * -2 is for unset! */
    int excite_;

    /** Is this function excitation sensitive? */
    int excite_sensitive_;

    /** @see FactorOmegaOmega() */
    bool omega_omega_;

    bool harmonic_;

    /** Conditions mark themselves as (non) linear -> no power in the design variable, ...*/
    bool linear_;

    /** Do we have local information? E.G. (global) slopes */
    Local* local;

    /** Do we have a filter element? Only for the projection function we store an alternative design set which we can filter */
    DesignSpace* projectionDesign_;

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
    DesignMaterial::Notation notation_;

};


} // namespace


#endif /* FUNCTION_HH_ */
