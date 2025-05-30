#ifndef CONDITION_HH_
#define CONDITION_HH_

#include <stddef.h>
#include <string>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/Enum.hh"
#include "General/Environment.hh"
#include "MatVec/Vector.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Function.hh"
#include "Utils/StdVector.hh"

namespace CoupledField {
class DesignStructure;
class ErsatzMaterial;
}  // namespace CoupledField

namespace CoupledField
{
   class DesignSpace;
   class MultipleExcitation;

   /** our constraint criteria. Can be filled directly from XML */
   class Condition : public Function
   {
     friend class ConditionContainer;

     public:
       /** Empty constructor for StdVector */
       Condition() { };
       /** virtual dtor because base class */
       virtual ~Condition() { };

       /** Helper constructor for AddCondition */
       Condition(PtrParamNode pn);

       /** overwrites Function::IsObjective() */
       bool IsObjective() const { return false; }

       /** to be overwritten in LocalCondition. */
       virtual bool IsLocalCondition() const { return false; }

       /** Overwrites and calls Function::PostProc() */
       void PostProc(DesignSpace* space, DesignStructure* structure, ErsatzMaterial* em);

       /** Call this method to append a Condition. This calls the actual (private) constructor.
        * Index is set with position of the relevant list.
        * If it is a homogenization constraint there might be a blow up resulting in several
        * constraints if either multiple entries are given or the entry position is 'all'.
        * @param pn determines active mode
        * @param constraints stuff is added here if the mode is constraint
        * @param observation stuff is added here in observation mode */
       static void AddCondition(PtrParamNode pn, StdVector<Condition*>& constraints,int i = -1, std::string entName = "");

       /** usually for constraints plus globalSlope for objective */
       typedef enum { EQUAL, LOWER_BOUND, UPPER_BOUND } Bound;

       static Enum<Bound> bound;

       /** Be sure not to mix up with Name! */             
       Bound GetBound() const { return bound_; }

       /** The bound value for inhomogeneous constraints.
        * In slack bound case it returns 0 as the constraints g <= slack is transformed to g - slack <= 0
        * @see IsSlackBound() */
       double GetBoundValue() const { return HasGeneralSlackBound() ? 0.0 : boundValue_; } // all slack constraints g <= slack need to be g - slack <= 0

       /** allows to compare with a special bound value as GetBoundValue() would return 0
        * @param compare SLACK_VALUE, ALPHA_VALUE ALPHA_PLUS_SLACK, ALPHA_MINUS_SLACK */
       double IsGeneralSlackBound(double compare) const { return boundValue_ == compare; }

       /** is the bound value one of the four special slack cases? includes alpha only */
       bool HasGeneralSlackBound() const { return boundValue_ == SLACK_VALUE || boundValue_ == ALPHA_VALUE || boundValue_ == ALPHA_MINUS_SLACK_VALUE || boundValue_ == ALPHA_PLUS_SLACK_VALUE; }

       /** Little helper to check if the bounds are violated (up to an eps) */
       bool IsFeasible() const;

       /** some own norm as we have no lagrange multipliers */
       double CalcFeasibility() const;

       /** Is this a feasibility constraint for FeasSCP */
       bool IsFeasibilityConstraint() const;

       /** Is this observation or active */
       bool IsObservation() const { return observation_; }

       /** active not in a active set optimization sense but !observation */
       bool IsActive() const { return !IsObservation(); }
       
       /** Only the local constraint in local mode are virtual */
       virtual bool IsVirtual() const { return false; }

       /** Check whether condition should be calculated for given region */
       bool IsForRegion(RegionIdType regionId);

       /** This is a nice statement for output which adds delta_logging and details for result output.
        * Contains the virtual element for slope */
       virtual std::string ToString() const;

       /** log to info.xml. Overloads Function::ToInfo() */
       void ToInfo(PtrParamNode in);
       
       /** Shall the scaling be linked to the objective scaling */
       bool DoObjectiveScaling() const { return objective_scaling_; }

       /** Is the gradient dense or sparse. Almost all local conditions and slack obj are sparse.
        * The only dense local function is localStress which is state dependent
        * @see Function::HasDenseJacobian() */
       bool HasDenseJacobian() const { return !IsLocalCondition() || type_ == LOCAL_STRESS  || type_ == LOCAL_BUCKLING_LOAD_FACTOR; } // todo: what is with slack?
       
       /** Is it a constraint on the imaginary part? */
       bool IsImag() const { return imag_; }

       /** Is it a constraint on the permeability? */
       bool IsBiisotropy() const { return biisotropy_; }

       /** When we do bloch, do we full bloch for all wave vectors? */
       bool DoFullBloch() const { return !bloch_extremal_; }

       // int GetFMOPosDefMinor() const { return fmo_pos_def_minor_; }

       /** creates an xml attribute name compatible string representation for coords */
       static std::string ToString(const StdVector<boost::tuple<int, int, double> >&);

       void DescribeProperties(StdVector<std::pair<string,string> >& map) const;

       /** The scaling is evaluated for external optimizers, not in OC!
        * This is the manual set scaling value - in objective_scaling_ case this value is ignored! */
       double manual_scaling_value = 0.0;

       /**The penalty formulation allows to add constraints via this penalty term to the objective.
        * Actually a penelty method finds iteratively the right value, in practice it is a given
        * parameter. Currently this is *only* implemented for the *Level-Set* method! */
       double penalty = 0.0;
       
       /** Shall delta constraints be printed? Is only true if a value is given! */
       bool delta_logging = false;

       /** Used for caching 1.0 / complete_volume per region */
       double volume_fraction = -1.0;

       /** For the homogenization tensor constraint this gives the actual position within the matrix_.
        * The first entry is for homogenization always set.
        * In the case of a "smart" isotropy constraint also E11-E22 = 0 and
        * E11-E12-2E33 = E11-E12-E33-E33 = 0 = (E11,1) + (E12,-1) + (E33,-2) are generated. Then coord is 2 or 3 entries.
        * Note, that the entries are 1-based!!!
        * the factor for ErsatzMaterial::CalcHomogenizedTensorEntry() */
       StdVector<boost::tuple<int, int, double> > coords;

       /** For design tracking, this are the elements we have to track. Function::elements is resized accordingly!
        * The vector is empty when we do not do design tracking */
       StdVector<double> pattern;

       /** for bloch eigenvalues which are extremal (searched and not full) */
       EigenInfo bloch;

       static double SLACK_VALUE;
       static double ALPHA_VALUE;
       static double ALPHA_MINUS_SLACK_VALUE;
       static double ALPHA_PLUS_SLACK_VALUE;

       // number of displacement constraints realized by multiple output constraints
       UInt output_multiple_nodes = 0;

    protected:
      /** Reads the coord attribute and sets the coord pair if value is not 'all'
       * @return false if 'all' and the coord pair is not set */
      bool ReadCoord(PtrParamNode pn);

      /** Add a subcondition with only index and value set (to zero) */
      Condition* AppendSubCondition(StdVector<Condition*>& list, bool biisotropy = false, bool imag = false);

      /** Create a new homogenization constraint with the given tensor position
       * @param base the base of cloning. Needs to contain a tensor!
       * @param list where to append the child. Index is set
       * @return the appendend child */
      Condition* AppendSubCondition(StdVector<Condition>& list, PtrParamNode entry_pos);

      /** @see other AppendSubCondition() */
      Condition* AppendSubCondition(StdVector<Condition*>& list, int pos_x, int pos_y);

      /** Bound stuff for condition and globalSlope also for objective */
      Bound bound_ = EQUAL;

      /** the bound value, the value_ attribute contains the function value */
      double boundValue_ = 0.0;

      bool delta_logging_ignored_ = false;

      /** shall we scale the condition with the objective scaling? */
      bool objective_scaling_ = false;

      /** Is this an observation constraint only. */
      bool observation_ = false;

      /** Some special constraints are automatically blown up - like isotropy. But
       * even then the first of the entries is NOT blown up!
       * Set by AppendSubCondition() */
      bool blown_up_ = false;

      /** This is only needed for biisotropy constraints, meaning that both permittitvity and permeability shall be isotropic
       *  biisotropy == true indicates isotropy constraints on the permeability  */
      bool biisotropy_ = false;

      /** Information if the constraint is set for the imaginary part
       *  default=false  */
      bool imag_ = false;

      /** this is the virtual base index of this condition w.r.t. all conditions.
       * For normal condition this is simple the virtual index, for local conditions this is the base*/
      int virtual_base_index_ = -1;

    private:

      /** Read the pattern for design tracking. pattern has in the end the same size as Function::elements.
       * Needs to be called after SetElements() */
      void ReadDesignTrackingPattern(DesignSpace* space, DesignStructure* structure);

      /** Helper for AddCondition().
       * Adds the conditions for isotropy or iso-orthotropy which is isotropy without fixing the
       * shear moduli */
      static void AddXtropyConstraints(PtrParamNode pn, StdVector<Condition*>& list, Condition* g);

      /** Helper for AddCondition().
       * Adds the number i to name of output nodes, necessary for displacement constraints */
      static void AddOutputConstraints(PtrParamNode pn, StdVector<Condition*>& list, Condition* g, int i, std::string entName);

      /** Helper for AddCondition() */
      static void AddHomogenizationTensorConstraints(PtrParamNode pn, StdVector<Condition*>& list, Condition* g);

      /** if in list a stress constraint is found, it is enlarged by the excitations and each is associated
       * to an own excitation */
      static void AddExcitationStressConstraints(StdVector<Condition*>& list, MultipleExcitation* me);

      /** for bloch mode each constraint is multiplied by wave vector which corresponds to excitation */
      static void AddBlochEigenConstraints(StdVector<Condition*>& list, MultipleExcitation* me);

      /** In the bloch case shall we have a constraint for every wave vector or search for the extremal */
      bool bloch_extremal_ = false;

   };

   /** This handles local constraints which exist only virtually - hence the optimizer sees them but
    * within optimization (e.g. log output) it is save to traverse all (real) elements.
    * The key element is the ConditionContainer() and the virtual local mode.
    * Examples are slope conditions and mole.
    * It is based on the local neighborhood information within Function::Local. Note, that there
    * are also global local functions like globalSlope using the Function::Local information. Only
    * when the neighborhood is used as local constraint with potentially many thousands elements,
    * this specialization of Condition is used. */
   class LocalCondition : public Condition
   {
   public:
     /** Helper constructor for AddCondition */
     LocalCondition(PtrParamNode pn);
     
     virtual ~LocalCondition() {};

     /** overwrites Condition::IsLocalCondition() */
     bool IsLocalCondition() const { return true;}

     /** PostInit when we have the design space */
     void PostInit(DesignSpace* space);

     /** The active index within the ConditionContainer::VirtualView iterator blowing up the slope constraints.
      * Requires (base) index_ and there are dim constraint s per design element.
      * @param view_index ranging from base_index plus dim * design.size.
      *        Set back to -1 after traversing! */
     void SetCurrentViewIndex(int view_index) {
       current_view_index_ = view_index;
     }

     /** The relative position within this local constraints
      * @return starts from 0 if valid, -1 if not valid
      * @see Function::GetCurrentRelativePosition() */
     virtual int GetCurrentRelativePosition() const {
       return std::max(current_view_index_ - virtual_base_index_, -1);
     }

     /** The absolute position of the local constraint
      * @return starts from 0 */
     int GetCurrentPosition() const {
       return current_view_index_;
     }

     int GetVirtualBaseIndex() const { return virtual_base_index_; }

     /** The number of slope constraints. */
     unsigned int GetConstraintSize() const {
       return local->local_values.GetSize();
     }

     /** The local mode has current_view_index_ set. For -1 we are in global mode.
      * Do no confuse with Function::IsLocalCondition() */
     bool IsLocal() const {
       return current_view_index_ != -1;
     }

     /** overwrites the base method */
     bool IsVirtual() const {
       return IsLocal();
     }

     /** overloaded version which gives in the local case only the 0 (no full neighborhood)
      * or 2 indices */
     StdVector<unsigned int>& GetSparsityPattern();

     /** @see Function::GetSparsityPatternSize() */
     unsigned int GetSparsityPatternSize() const;

     /** overloaded version which implements the functionality for the determinant functions */
     Matrix<unsigned int>& GetHessianSparsityPattern();

     /** @see Function::CalcHessian() */
     void CalcHessian(StdVector<double>& out, double factor);

     /** This is the local context currently requested by the optimizer */
     Function::Local::Identifier& GetCurrentVirtualContext();
     const Function::Local::Identifier& GetCurrentVirtualContext() const;

     /** Calculates the mean |value| based on local_values */
     double CalcMeanAbsValue() const;

     /** Calculates the max |value| if upper bound or equal and min if lower bounded */
     double CalcMinMaxAbsValue() const;

     /** count the infeasible elements where the value is out of bound. This can be problematic for linear functions! */
     int CountInfeasibles() const;

     /** Overloads the base method. If in special mode element value is returned. Otherwise
      * the max norm is returned (calculated on the fly */
     double GetValue() const;

     /** Overloaded to set the local blown up values if in local mode.  */
     void SetValue(double val);

     /** overloads ToString() to add local information if in local mode. For debug logging */
     std::string ToString(MultipleExcitation* me = NULL) const;

   private:

     /** To be set via SetCurrentViewIndex(). Negative if not initialized. */
     int current_view_index_ = -1;
   };

   /** This is a container for the conditions within Optimization.
    * It holds the active and inactive (observation) conditions.
    * Constraints are called conditions as the is an older type "Constraint" in CFS
    */
   class ConditionContainer
   {
   public:

     ConditionContainer();

     ~ConditionContainer();

     /** The slope constraint is a super constraint which represent from the external
      * optimizers point of view up to several thousand constraints. Therefore
      * the optimizers shall only access it by this class. */
     class VirtualView
     {
     public:

       VirtualView(ConditionContainer* constraints);

       /** The constraints have been changed (slope constraints initialized) */
       void Refresh();

       /** call this after traversing via Get() to switch a potential slope constraint back to global mode */
       void Done();

       /** handles slope constraints. Note that observe is always after active!
        * @param index from 0 to NumberOfTotalConstraints() which might be several thousands if there
        *        is a slope constraint
        * @return in the slope constraint case a "tunded" SlopeCondition object */
       Condition* Get(int index);

       /** When there is no slope constraint this is active otherwise several thousands */
       int GetNumberOfActiveConstraints() const {
         return virtual_active_size_;
       }

       /** Adds observe to active - the slope can also be an observe!
        * Is never smaller NumberOfActiveConstraints(). */
       int GetNumberOfTotalConstraints() const {
         return virtual_total_size_;
       }


     private:
       /** This are the real condition indices of local conditions. Sorted. Used to calculate and navigate in the virtual
        * condition space */
       StdVector<unsigned int> local_cond_index_;

       ConditionContainer* container_;

       int virtual_active_size_;
       int virtual_total_size_;
     };

     /** Process the xml parameters. To be called only once. PostProc() can be called
      * more often within the Optimization subclasses.
      * The slope constraint cannot be processed w/o DesignSpace. It needs to be initialized
      * by PostProc() later
      * @param pn_cond the list of "condition" from the xml file
      * @see PostProc() */
     void Read(ParamNodeList pn_cond);

     /** Set up the all vector, does a (re)-indexing of all constraints and sets up/ refreshes the
      * virtual view. Call any time when the number of constraints has been changed. Save for first-time
      * call, is then a setup */
     void Refresh();

     /** The slope constraints can only be initialized when the design exists.
      * Requires ToInfo() to be called prior such that we can do the info output to the stored info
      * @structure is from ErsatzMaterial */
     void PostProc(DesignSpace* space, DesignStructure* structure, MultipleExcitation* me, ErsatzMaterial* em);

     /** Log the head information. The InfoNode is stored such that PostProc can do the info output
      * if already set. */
     void ToInfo(PtrParamNode in);

     /** Searches in active constraints or all constraints !
      *  @param design NO_TYPE ignores this criteria. DEFAULT would be problematic for
      *                this purpose as it is a valid value
      * @return check active flag! Not NULL! */
     Condition* Get(Condition::Type type = Condition::VOLUME, DesignElement::Type design = DesignElement::NO_TYPE, bool only_active = true) {
       return Get(type, design, only_active, true);
     }
     Condition* Get(Condition::Type type, DesignElement::Type design, Condition::Bound bound, bool throw_exception = true);

     /** Search for the condition by name */
     Condition* Get(const std::string name, bool throw_exception = true);

     StdVector<Condition*> GetList(Condition::Type type, DesignElement::Type design = DesignElement::NO_TYPE, bool only_active = true, Function::Access access = Function::NO_ACCESS);

     /** query before Get() throws an exception */
     bool Has(Condition::Type type = Condition::VOLUME, DesignElement::Type design = DesignElement::NO_TYPE, bool only_active = true);

     /** if so we need to add the bound to the name */
     bool RequiresBoundForUniqueness(const Condition* g);

     /** is at least one constraint state sensitive? */
     bool IsAllStateDependent() const;

     /** are the active constraints feasible */
     bool IsFeasible() const;

     void PushBackHistory();

     /** All external optimizers should only work with this view.
      * It make the special handling for the slope constraints */
     VirtualView* view;

     /** This is a virtual container combining standard + observe */
     StdVector<Condition*> all;

     /** The real constraints which are evaluated by the optimizer. Has nothing to do with
      * "active" in the sense of active sets but means not observe.
      * Be sure that you should not use the VirtualView!! (external optimizers!) */
     StdVector<Condition*> active;

     /** The "inactive" constraints with are only evaluated for printing */
     StdVector<Condition*> observe;

   private:

     Condition* Get(Condition::Type type, DesignElement::Type design, bool only_active, bool throw_exception);

     /* has unique bounds */
     bool HasUniqueBounds(const StdVector<Condition*>&);

     /** save for maxSlope output */
     DesignSpace* space_;
   };

} // namespace


#endif /*CONDITION_HH_*/
