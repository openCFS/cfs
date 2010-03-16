#ifndef OPTIMIZATION_HH_
#define OPTIMIZATION_HH_

#include "General/Enum.hh"
#include "Utils/StdVector.hh"
#include "MatVec/vector.hh"
#include "Driver/harmonicDriver.hh"
#include "Domain/bcs.hh"
#include "Optimization/DesignElement.hh"
#include "Optimization/Condition.hh"
#include "Optimization/Objective.hh"

namespace CoupledField
{
   class Elem;
   class ParamNode;
   class ParamNode;
   class DesignSpace;
   class OptimalityCondition;
   class BaseOptimizer;
   class SinglePDE;
   class Excitation;
   class LinearFormContext;

  /** This is a general optimization object. The optimiziation loop is around
   *  domain->GetDriver()->SolveProblem() and as such general. Note convention,
   * that Optimization solves an optimization problem and an Optimizer is the
   * solver to be used.*/
  class Optimization
  {
      public:

        /** it is very important to have the base destructor virtual, otherwise the
         * derived destructor might not be executed */
        virtual ~Optimization();


         /** This is an implementation of the factory pattern. The method
          * reads the xml file.
          * @return either the described optimization problem or null if nothing in the xml file */
         static Optimization* CreateInstance();

         /** Solves the problem by looping over driver->SolveProblem() up to
          * minimum or max iterations is reached. Overwrite on demand*/
         void SolveProblem();

         /** Where we apply the transformation.
          * A subset of the values are PDE identifiers for ToPDE() and ToApp(). */
         typedef enum { MECH, ELEC, PIEZO_COUPLING, PRESSURE, CHARGE_DENSITY, MASS, HEAT, NO_APP} Application;

         /** Not the optimization problem but the solver! */
         typedef enum { OPTIMALITY_CONDITION, IPOPT_SOLVER, SCPIP_SOLVER, SNOPT_SOLVER, SHAPE_SOLVER, 
												EVALUATE_INITIAL_DESIGN, GRADIENT_CHECK } Optimizer;

         /** to convert string/enum for this type */
         static Enum<Optimizer> optimizer;

         static Enum<Application> application;

         /** the commit mode defines what of the iterations is to be written to gid, ... */
         typedef enum { FORWARD, ADJOINT, BOTH } CommitMode;

         static Enum<CommitMode> commitMode;
         

         /** The DesignSpace element is a container for the complex DesignElements.
          * There are service methods in the container to exchange with external
          * design space arrays as for external optimizers */
        DesignSpace* GetDesign() { return design; }

        /** A cost function is (here) always a scalar type. The value is recalculated when needed!
         * This also stores the value in history of the CostFunction */
        virtual double CalcObjective() = 0;

        /** Evaluates the cost-function gradient w.r.t. the design space. Apply SetDesignSpace() first!
         * Writes to DesingElement.objective_gradient.
         * Does a multiplication with Excitation::GetWeightedFactor(). Note, that 
         * the Calc* methods for the objective do only Excitation::GetFactor().
         * @param grad_out size is GetDesignSpaceSize(). If null only DesingElement.objective_gradient */
        virtual void CalcObjectiveGradient(double* grad_out) = 0;

        /** Determines the constraint.
         * @param which constraint to calc? Default is the only one! */
        virtual double CalcConstraint(Condition* constraint = NULL) = 0;

        /** evaluates the gradient of the cost function by the desing element
         * Writes to DesignElement.contraint_gradient
         * @see CalcObjectiveGradient() for parameter description */
        virtual void CalcConstraintGradient(Condition* constraint = NULL, double* grad_out = NULL) = 0;

        /** This is brute force debug method which calculates the symmetry of a sqared
         * model with horizontal symmetry axis with lexicographic order (at least works for gid).
         * If there is a special result index for vs, access the relative element symmetry errors
         * are written there.
         * @return the average relative symmetry error. */
        double CalcSymmetry(DesignElement::Type de, DesignElement::ValueSpecifier vs, DesignElement::Access access);

        /** This method checks if special results requiring special evaluation are there.
         * This works for now only for CalcSymmetry() */
        void EvaluateSpecialResults();

        /** Evaluates the state problem, does not increment the iteration counter.
         * @param excite provide for multiharmonic steps */
        virtual void SolveStateProblem(Excitation* excite = NULL);

        /** The maximal number of iterations */
        int GetMaxIterations() { return maxIterations; }

        /** The current iteration */
        int GetCurrentIteration() { return currentIteration; }

        /** External optimizers will solve more state problems than doing
         * iterations -> e.g. for doing a line search. With this it is signalled,
         * that the last SolveStateProblem() was actually an iteration.
         * <p>Normally, this will increment the iteration counter and trigger any output
         * stuff (CFS-output-writing, console log, log-file, density xml file)</p>
         * <p>Also makes a push back of cost.value to cost.history</p>
         * For the multiharmonic case one can commit each frequency but keep the
         * iteration number.
         * @see FinalizeStoreResults() for calling after the last CommitIteration()
         * @param keep_iteration_number will keep currentIteration and not rest problemsWithinIteration
         * @return the iteration element we added */
        virtual PtrParamNode CommitIteration(bool keep_iteraton_number = false);

        /** the break condition for the optimization loop.
         * Checks the stopping rule from the XML file an searches for an HALTOPT file.  */
        virtual bool DoStopOptimization();

        /** are we in the harmonic case? */
        bool IsHarmonic() const { return harmonic; }

        /** This struct stores the multiple excitation Information. It contains the
         * same defaults as the xml schema as the whole element is optional. */
        class MultipleExcitation
        {
        public:
          /** @param enabled comes from the attribute.
           * @param pn if NULL only the defaults are set */
          MultipleExcitation(bool enabled, PtrParamNode pn);

          void ToInfo(PtrParamNode in) const;

          typedef enum { NO_TYPE, FIXED_WEIGHT, META_OBJECTIVE, HOMOGENIZATION_TEST_STRAINS } Type;

          static Enum<Type> type;
          /** Do we do multiple excitation at all? */
          bool IsEnabled() const { return multiple_excitation_; }

          /** Dow we do adjust weigts */
          bool DoAdjustWeights() const { return type_ == META_OBJECTIVE; }

          bool DoHomogenization() const { return type_ == HOMOGENIZATION_TEST_STRAINS; }

          /** The stride for adjust weights: 1 = every iteration, 2 = every second ... */
          int stride;

          /** the maximal span between the largest (1) and smallest weight as factor */
          double max_gain;

          /** The exponent d in w_k^p J_k = const */
          double damping;

        private:
          /** do we do multiple excitation at all? */
          bool multiple_excitation_;
          Type type_;
        };

        /** The current multiple excitation state -> check with IsEnabled() */
        MultipleExcitation* GetMultipleExcitation() const { return multiple_excitation; }

        /** The constraints we have */
        StdVector<Condition> constraints;

        /** The "inactive" constraints with output_only mode in xml */
        StdVector<Condition> outputs;

        /** Searches in active and output only constraints!
         * TODO: make name default for only one constraint
         *  @param design NO_TYPE ignores this criteria. DEFAULT would be problematic for
         *                this purpose as it is a valid value
         * @return check active flag!  */
        Condition& GetConstraint(Condition::Name name, DesignElement::Type design = DesignElement::NO_TYPE);

        /** set the (static) enums - if they are used outside optimization, make this method public */
        static void SetEnums();

        /** Our base ParamNode pointer, pointing to a plain <optimization> */
        PtrParamNode optInfoNode;

        /** This is the list of concurrent objective functions.
         * It is guaranteed to have at least one entry.
         * The objectives are almost identical (share most attributes) */
        ObjectiveContainer objectives;

      protected:
        /** Set up the optimization system e.g. prepare the domain for optimization. called
         * exclusively by CreateInstance() -> don't forget to call PostInit() afterwards! */
        Optimization();

        /** Appends to the current analysis_id of the driver a child and
         * sets analysis_id and excite accordingly */
        PtrParamNode CreateAdjointAnalysisIdNode();
        
        /** This tells the driver to store the last solved problem (gid, ...). Called in
         * CommitIteration(). For PiezoSIMP we can save more often and there this method
         * is overwritten and might do nothing.
         * @param step_val the "label" of the "transient" step. -1 is the integer counter */
        virtual void StoreResults(double step_val = -1.0);

        /** called by CommitIteration(), to be overwritten if additional data should be
         * written, then logFileHeader should also be set. Don't add a new-line here!!.
         * @param iteration a duplicate of the log file output to the info xml file */
        virtual void LogFileLine(std::ofstream* out, PtrParamNode iteration);

        /** PostInit is to be called after the constructor. 
         * PostInit  does not really solve something. */
        virtual void PostInit();

        /** This is the second phase of post initialization. It creates the Optimizer tools */
        virtual void PostInitSecond();
        
        
        /** Gives back the current frequency for printing. This is not the current frequency
         * in multifrequency case. Not a fast method!
         * @return an empty string in the non-harmonic case */
        virtual std::string GetIterationFrequency() { return "not implemented"; }

        /** The way the date is stored. Forward/ adjoint/ both and stride. 
         * Set in the <commit> element */
        CommitMode commitMode_;

         /* Set in the <commit> element.
          * stride = 0 corresponds to infinity (first and last)
         * stride = 1 is every iteration
         * stride = 2 is the first, the third, the fifth and the last */ 
        int commitStride; 

        /** The current iteration, 0 is the first run. Note that the state problem might be
         * executed more often (-> line search).
         * @see problemSolvedCounter. */
        int currentIteration;

        /** This checks how often the state problem is solved. This is not necessary equal
         * to the iterations, e.g. for line search of an external optimizer. Incremented by
         * SolveStateProblem(). */
        int problemSolvedCounter;

        /** The number of problems solved in this iteration */
        int problemWithinIteration;

        /** the maximum number of iterations, have a default */
        int maxIterations;

        /** Our MultipleExcitation objecte - by default disabled */
        MultipleExcitation* multiple_excitation;

        /** The actual kind of optimizer.  */
        Optimizer optimizer_;

        /** Here we contain our design space. The domain gets a reference to it to perform
         * the ersatz material ansatz */
        DesignSpace* design;

        /** Here we keep the last iterations design space */
        Vector<double>  last_iteration;

        /** Here we keep the last evaluation design space */
        Vector<double>  last_evaluation;

        /** are we harmonic or static? */
        bool harmonic;

      private:
        /** CommitIteration() does not necessary store the results when we have a stride
         * set in the <commit> element. In case this method makes a StoreResults (not commit)
         * if necessary. Should always be called when finished (in the good and bad sense) */
        void FinalizeStoreResults();

        /** Read the objective functions. Fills the objectives list with almost identical entries.
         * Handles multiObjectives */
        void ReadObjectives(PtrParamNode obj_node);

        /** Encapsulates Logging information */
        class Log
        {
        public:
          /** Sets to meaningul defaults (don not much :) ) */
          Log();

          /** Closes the file */
          ~Log();

          /** @param log_name is interpreted. If allows a file, the logFile is created.
           * @param pn_log pointer to the 'log' element. Might be NULL */
          void Init(const std::string& log_name, PtrParamNode pn_log);

          /** The header of the logFile_, to be overwritten if LogFileLine() is overwritten. CommitIteration()
           * writes this string to logFile_ a the first execution */
           std::string fileHeader;

           /** if set write the design to the logfile */
           bool design;

           /** if set write the gradient of the design to logfile */
           bool designGradient;

           /** optional log the iterations and cost value to a file to gnuplot it */
           std::ofstream* file;
        };

        /** Keeps all logging relevant stuff */
        Log log;

        /** When did we store the last result via CommitIteration() due to stride */
        int lastStoredResult_;

        /** This holds our optimer instance. */
        BaseOptimizer* baseOptimizer_;
  };

  typedef LoadBc TrackingBc;

  typedef LoadList TrackingList;

  /** For multiple loads (compliance or multiple frequency optimization) we use
   * the summarized term multiple excitations. This object encapsulates the such a excitation.
   * This excitations are weighted.  */
  class Excitation
  {
  public:

    /** default constructor for StdVector() */
    Excitation();

    /** This method makes the current load active.
     * For multiple frequencies it does nothing. The actual frequency is choosen by default. */
    void Apply();

    /** Find the fixed factor, does ignore weighting and does not apply it. */
    double GetFactor(Objective* cost);

    /** Returns GetFactor() * normalized_weight */
    double GetWeightedFactor(Objective* cost);
    
    /** Gets the current omege =  2 * pi * f */
    double GetOmega();

    /** @return omega^2 */
    double GetOmegaOmega();
    
    /** read the tracking node list from XML */
    void ReadTrackings(PtrParamNode ts);
    
    /** read the loads from XML */
    void ReadLoads(PtrParamNode ls);
    
    void ReadTestStrain(const Vector<double>& vec);

    /** return pointer to linForms, used by Shape-Optimization */
    std::set<LinearFormContext*>* GetLinForms() { return &linForms; }

    /** the index of this excitation in the excitations array. If -1 something went wront */
    int index;

    /** For static/monoharmonic optimization with different load-cases. Now allowing also multiple loads in one case.
     * empty and apply_linForms=false if not applicable */
    LoadList loads;
    
    /** For static optimization with different pressure or regionLoads */
    std::set<LinearFormContext*> linForms;
    
    /** if linForms are to be applied 
     * set true in multiple load/pressure/regionLoad per load-case or tracking */
    bool apply_linForms;
    
    /** Different possible trackings for tracking objective */
    TrackingList trackings;

    /** This is a link to the Frequency description from the harmonic driver.
     * It is used for calling the HarmonicDriver to solve the problems */
    HarmonicDriver::Frequency* f_link;

    /** For multiharmonic excitation. -1.0 by default */
    double frequency;

    /** this is the weight from the xml file */
    double weight;

    /** this is the normalized weight (sum of all weights of all excitations is 1) */
    double normalized_weight;

    /** Here we store the calculated objective value, including costFunction/factor
     * to enable metaObjective. In the Optimization::cost struct we store a copy/average */
    double cost;

    /** Here we might store the "original" (@see objective) gradient for analysis/output */
    StdVector<double> cost_gradient;
    
    /** for the calculation of a homogenized material tensor, one can specify test strains in
     *  xml file. */
    Vector<double> test_strain;
  };

} // namespace


#endif /*OPTIMIZATION_HH_*/

