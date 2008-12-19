#ifndef OPTIMIZATION_HH_
#define OPTIMIZATION_HH_

#include "General/Enum.hh"
#include "Utils/StdVector.hh"
#include "Driver/harmonicDriver.hh"
#include "Domain/bcs.hh"
#include "Optimization/DesignElement.hh"
#include "Optimization/Condition.hh"


namespace CoupledField
{
   class Elem;
   class ParamNode;
   class InfoNode;
   class DesignSpace;
   class OptimalityCondition;
   class BaseOptimizer;
   class SinglePDE;
   class Excitation;

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

         /** Where we apply the transformation */
         typedef enum { MECH, ELEC, PIEZO_COUPLING, PRESSURE, CHARGE_DENSITY, MASS, SURFACE_NORMAL, NO_APP} Application;

         /** The taks for a cost function either to minimize or maximize */
         typedef enum { MINIMIZE, MAXIMIZE } ObjectiveTask;

         /** Known types of cost functions */
         typedef enum { COMPLIANCE,                /*!< (u,f) the opposite of stiffness */
                        OUTPUT,                    /*!< (u,l) maximize solution where vector l is not 0 */
                        DYNAMIC_OUTPUT,            /*!< (u, L conj(u)) as OUTPUT but complex */
                        CONJUGATE_COMPLIANCE,      /*!< (u, F conj(u)) as DYNAMIC_OUTPUT with trace of L is f */
                        GLOBAL_DYNAMIC_COMPLIANCE, /*!< (u, I conj(u)) as DYNAMIC_OUTPUT with L is I (everywhere) */
                        RADIATION,                 /*!< maximized v on surface elements (Du and Olhoff) */
                        VOLUME, TRACKING } ObjectiveType;

         /** Not the optimization problem but the solver! */
         typedef enum { OPTIMALITY_CONDITION, IPOPT_SOLVER, SCPIP_SOLVER, LEVEL_SET, EVALUATE_INITIAL_DESIGN } Optimizer;

         /** to convert string/enum for this type */
         static Enum<ObjectiveType> objectiveType;

         /** to convert string/enum for this type */
         static Enum<Optimizer> optimizer;

         static Enum<Application> application;

         /** We combine the cost function in a set to handle multiple of it.
          * It contains static const elements (and  working stuff).
          * MultipleExciation is in the XML file part of the objective but in class part of Optimization*/
         class Objective
         {
           public:
             /** @param harmonic are we in a harmonic simulation? */
             Objective(ParamNode* pn, bool harmonic);

             /** The actual kind of cost function. */
             ObjectiveType type;

             /** The task is the direction of a cost function (MINIMIZE, MAXIMIZE) */
             ObjectiveTask task;

             /** This vector stores the cost functions of the iterations. Written in GetObjective() */
             StdVector<double> history;

             /** The objective value */
             double GetValue();

             /** @param val set an unscaled value here -> it is scaled in the return */
             void SetValue(double val);

             void ToInfo(InfoNode* in) const;

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

             /** Here we store out ParamNode such we can more easily access it in ErsatzMaterial */
             ParamNode* pn;

             /** Shall harmonic optimization multiply with omega^2.
              * This makes "u L conj(u)" to actually calc "v L conj(v)" with v = du/dt. -> approximatates sound intensity */
             bool FactorOmegaOmega() {
               return omega_omega_;
             }

           private:

             /** The current value */
             double value_;

             /** @see FactorOmegaOmega() */
             bool omega_omega_;

             bool harmonic_;
         };

         /** The cost function */
         Objective* GetObjective() { return cost; }

         /** The DesignSpace element is a container for the complex DesignElements.
          * There are service methods in the container to exchange with external
          * design space arrays as for external optimizers */
        DesignSpace* GetDesign() { return design; }

        /** A cost function is (here) always a scalar type. The value is recalculated when needed!
         * This also stores the value in history of the CostFunction */
        virtual double CalcObjective() = 0;

        /** Evaluates the cost-function gradient w.r.t. the design space. Apply SetDesignSpace() first!
         * Writes to DesingElement.objective_gradient.
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
         * iteration number
         * @param keep_iteration_number will keep currentIteration and not rest problemsWithinIteration
         * @return the iteration element we added */
        virtual InfoNode* CommitIteration(bool keep_iteraton_number = false);

        /** the break condition for the optimization loop */
        virtual bool IsMinimumReached();

        /** are we in the harmonic case? */
        bool IsHarmonic() const { return harmonic; }

        /** This struct stores the multiple excitation Information. It contains the
         * same defaults as the xml schema as the whole element is optional. */
        class MultipleExcitation
        {
        public:
          /** @param enabled comes from the attribute.
           * @param pn if NULL only the defaults are set */
          MultipleExcitation(bool enabled, ParamNode* pn);

          void ToInfo(InfoNode* in) const;

          /** Do we do multiple excitation at all? */
          bool IsEnabled() const { return multiple_excitation_; }

          /** Do we java a meta_objective */
          bool DoMetaObjective() const { return meta_objective_; }

          /** Dow we do adjust weigts */
          bool DoAdjustWeights() const { return meta_objective_; }

          /** The stride for adjust weights: 1 = every iteration, 2 = every second ... */
          int stride;

          /** the maximal span between the largest (1) and smalles weight as factor */
          double max_gain;

          /** The exponent d in w_k^p J_k = const */
          double damping;
          
        private:
          /** do we do multuple excitation at all? */
          bool multiple_excitation_;
          bool meta_objective_;
        };

        /** The current multiple excitation state -> check with IsEnabled() */
        MultipleExcitation* GetMultipleExcitation() const { return multiple_excitation; }

        /** The constraits we have */
        StdVector<Condition> constraints;

        /** The "inactive" constraits with output_only mode in xml */
        StdVector<Condition> outputs;

        /** Searches in active and output only constrints!
         * TODO: make name default for only one constraint
         *  @param design NO_TYPE ignores this criteria. DEFAULT would be problematic for
         *                this purpose as it is a valid value
         * @return check active flag!  */
        Condition& GetConstraint(Condition::Name name, DesignElement::Type design = DesignElement::NO_TYPE);

        /** set the (static) enums - if they are used outside optimization, make this method public */
        static void SetEnums();

        /** Our base InfoNode pointer, pointing to a plain <optimization> */
        InfoNode* optInfoNode;

      protected:
        /** Set up the optimization system e.g. prepare the domain for optimization. called
         * excusively by CreateInstance() -> don't forget to call PostInit() afterwards! */
        Optimization();

        /** This is the comment for Driver::SolveProblem() which becomes part of the
         * filenames when expoiting the linear system.
         * @param excite if given extracts frequency or load step */
        std::string GetSolveComment(Excitation* excite = NULL);

        /** This tells the driver to store the last solved problem (gid, ...). Called in
         * CommitIteration(). For PiezoSIMP we can save more often and there this method
         * is overwritten and might do nothing.
         * @param step_val the "label" of the "transient" step. -1 is the integer counter */
        virtual void StoreResults(double step_val = -1.0);

        /** called by CommitIteration(), to be overwritten if additional data should be
         * written, then logFileHeader should also be set. Don't add a new-line here!!.
         * @param iteration a duplicate of the log file output to the info xml file */
        virtual void LogFileLine(std::ofstream* out, InfoNode* iteration);

        /** PostInit is to be called after the constructor. This is required for the
         * optional IPOPT which expects the Optimization during construction. PostInit
         * does not really solve something. */
        virtual void PostInit();

        /** Gives back the current frequency for printing. This is not the current frequency
         * in multifrequency case. Not a fast method!
         * @return an empty string in the non-harmonic case */
        virtual std::string GetIterationFrequency() { return "not implemented"; }

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

         /** Up to now we have only one cost function */
         Objective* cost;

         /** Here we contain our design space. The domain gets a reference to it to perform
          * the ersatz material ansatz */
         DesignSpace* design;

         /** The header of the logFile_, to be overwritten if LogFileLine() is overwritten. CommitIteration()
          * writes this string to logFile_ a the first execution */
          std::string logFileHeader;

         /** Here we keep the last iterations design space */
         Vector<double>  last_iteration;

         /** Here we keep the last evaluation design space */
         Vector<double>  last_evaluation;

         /** are we harmonic or static? */
         bool harmonic;

     private:

        /** optional log the iterations and cost value to a file to gnuplot it */
        std::ofstream*     logFile_;

        /** This holds our optimer instance. */
        BaseOptimizer* baseOptimizer_;
  };

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

    /** Gets the current omege =  2 * pi * f */
    double GetOmega();

    /** @return omega^2 */
    double GetOmegaOmega();

    /** the index of this excitation in the excitations array. If -1 something went wront */
    int index;

    /** For static/monoharmonic optimization with different load-cases.
     * NULL if not applicable. */
    shared_ptr<LoadBc> load;

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
  };



} // namespace


#endif /*OPTIMIZATION_HH_*/

