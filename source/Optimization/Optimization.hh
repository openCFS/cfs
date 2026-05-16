#ifndef OPTIMIZATION_HH_
#define OPTIMIZATION_HH_

#include <stddef.h>
#include <iosfwd>
#include <string>
#include <tuple>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/Enum.hh"
#include "General/defs.hh"
#include "MatVec/Vector.hh"
#include "Optimization/Condition.hh"
#include "Optimization/Context.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Objective.hh"
#include "Optimization/OptimizationMaterial.hh"
#include "Utils/StdVector.hh"
#include "Utils/Timer.hh"

namespace CoupledField
{
   class Assemble;
   class BaseOptimizer;
   class BiLinFormContext;
   class DesignSpace;
   class Excitation;
   class Function;
   class Grid;
   class LinearFormContext;
   class MultipleExcitation;
   class StdPDE;
   class Tune;

   // FIXME: this is originally from timestepping.hh and has to be replaced
   typedef enum {NO_DERIVTYPE = 0, FIRST_DERIV = 1, SECOND_DERIV = 2} TimeDeriv;

   /** This is a simple class used as a parameter to SetAdjointRhs called by assemble
    * only for tracking and transient
   class AdjointParameters {
   public:
     AdjointParameters(Function* f, Excitation* e) {
       function = f;
       excite = e;
     }
     Function* GetFunction() {
       return function;
     }
     Excitation* GetExcitation() {
       return excite;
     }
   private:
     Function* function;
     Excitation* excite;
   };*/

  /** This is a general optimization object. The optimiziation loop is around
   *  Driver::SolveProblem() and as such general. Note convention,
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

         /** just queries optimization/ersatzMaterial/material */
         static OptimizationMaterial::System ParseSystem();

         /** PostInit is to be called after the constructor. */
         virtual void PostInit();

         /** This is the second phase of post initialization. It creates the Optimizer tools */
         virtual void PostInitSecond();

         /** Solves the problem by looping over driver->SolveProblem() up to
          * minimum or max iterations is reached. Overwrite on demand*/
         void SolveProblem();

         /** Not the optimization problem but the solver! */
         typedef enum { OCM_SOLVER, IPOPT_SOLVER, SCPIP_SOLVER, SNOPT_SOLVER, PYTHON_SOLVER, DUMAS_MMA, DUMAS_GCMMA,
                        FEAS_PP_SOLVER, MMA_SOLVER, SGP_SOLVER, SHAPE_SOLVER, EVALUATE_INITIAL_DESIGN, GRADIENT_CHECK  } Optimizer;

         /** to convert string/enum for this type */
         static Enum<Optimizer> optimizer;

         static Enum<App::Type> application;

         /** the commit mode defines what of the iterations is to be written to gid, ... */
         typedef enum { FORWARD, ADJOINT, BOTH, EACH_FORWARD, EACH_ADJOINT } CommitMode;

         static Enum<CommitMode> commitMode;

         /** The DesignSpace element is a container for the complex DesignElements.
          * There are service methods in the container to exchange with external
          * design space arrays as for external optimizers */
        DesignSpace* GetDesign() { return design; }

        /** A cost function is (here) always a scalar type. The value is recalculated when needed!
         * This also stores the value in history of the CostFunction */
        double CalcObjective(Excitation* ev_only_excite = NULL);

        /** Evaluates the cost-function gradient w.r.t. the design space. Apply SetDesignSpace() first!
         * Writes to DesingElement.objective_gradient.
         * Does a multiplication with Excitation::GetWeightedFactor(). Note, that
         * the Calc* methods for the objective do only Excitation::GetFactor().
         * @param grad_out size is GetDesignSpaceSize(). If null only DesingElement.objective_gradient */
        virtual void CalcObjectiveGradient(StdVector<double>* grad_out, Excitation* ev_only_excite = NULL);

        /** Determines the constraint.
         * The function value is stored in value_
         * @param which constraint to calc? Default is the only one! */
        virtual double CalcConstraint(Condition* constraint, Excitation* ev_only_excite = NULL);

        /** evaluates the gradient of the cost function by the desing element
         * Writes to DesignElement.contraint_gradient
         * @see CalcObjectiveGradient() for parameter description */
        virtual void CalcConstraintGradient(Condition* constraint = NULL, StdVector<double>* grad_out = NULL, Excitation* ev_only_excite = NULL);

        /** Evaluates objective and constraint function and gradient.
         * Implemented in ErsatzMaterial
         * @param grad_out only used in derivative case
         * @return zero for derivative */
        virtual double CalcFunction(Excitation& excite, Function* f, bool derivative) = 0;

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
         * This is actually a mere helper for the overwritten version in ErsatzMaterial!
         * @param excite provide for multiharmonic steps */
        virtual void SolveStateProblem(Excitation* excite = NULL);

        /** Traverses all objective and active constraint functions (non local) and calls SolveAdjointProblem()
         * if the functions needs it.
         * The reason for separating state and adjoint problem is that we need the adjoint for gradient evaluation
         * and in the line search case we need only states and not the adjoint
         * @see Function::IsAdjointBased() */
        virtual void SolveAdjointProblems(Excitation* excite = NULL);

        /** Solves the Adjoint problem, for given excite and objective
         * This does the real work
         * Only Bastian had an implementation for transient and tracking in Optimization.cc
         * @param excite multi-excitation
         * @param cost multi-objective */
        virtual void SolveAdjointProblem(Excitation* excite, Function* f) { assert(false); }

        /** Sets the rhs for the adjoint, called by assemle */
        // only for transient and tracking
        // virtual void SetAdjointRhs(AdjointParameters* adjointParams) = 0;

        /** The maximal number of iterations */
        int GetMaxIterations() { return maxIterations; }

        /** The current iteration */
        int GetCurrentIteration() { return currentIteration; }

        /** External optimizers will solve more state problems than doing
         * iterations -> e.g. for doing a line search. With this it is signaled,
         * that the last SolveStateProblem() was actually an iteration.
         * <p>Normally, this will increment the iteration counter and trigger any output
         * stuff (CFS-output-writing, console log, log-file, density xml file)</p>
         * <p>Also makes a push back of cost.value to cost.history</p>
         * For the multiharmonic case one can commit each frequency but keep the
         * iteration number.
         * @see FinalizeStoreResults() for calling after the last CommitIteration()
         * @return the iteration element we added */
        virtual PtrParamNode CommitIteration();

        /** the break condition for the optimization loop.
         * Checks the stopping rule from the XML file an searches for an HALTOPT file.
         * Shall be called after CommitIteration() !
         * user_stop_reason_ contains the reason if true is returned */
        bool DoStopOptimization();

        /** for python get_opt_stopping_rules() */
        StdVector<std::pair<string,string> > GetStoppingRules() const;

        /** allow Python to cancel the current optimization.
         * After the call, the next DoStopOptimization() will return true.
         * Uses user_break_message and user_break_converged.
         * @param args tuple with boolean if converged and message  */
        void PythonStopOptimization(PyObject* args);

        /** set the break status
         * @param reason if reason is not set, also converged is ignored */
        PtrParamNode DoStopOptimizationHelper(bool converged = false, const string& reason = "");

        /** allow Python to query the current optimizer and some of its properties
         * @return string/string dict with 'optimizer' and general/specific properties. Only few can be set
         * @see BaseOptimizer::PythonSetProperty() */
        PyObject* PythonGetOptimizerProperties() const;

        /** are we in transient optimization?
         * FIXE -> Context */
        static bool IsTransient();

        /** in transient, first step can be static, so that start displacement can depend on material parameters */
        bool IsFirstTransientStepStatic() const {return firstStepStatic; };

        /** in transient, first step can be static and have a different weight for the objective, this is the weight
         * converted so that it just has to be multiplied */
        double GetStepWeight(unsigned int ts) const;

        /** The current multiple excitation state -> check with IsEnabled() */
        MultipleExcitation* GetMultipleExcitation() const { return me; }

        /** Get the list of excitations given in multipleExcitations in the costFunction */
        ParamNodeList GetMultipleExcitionsNodes();

        /** set the (static) enums - if they are used outside optimization, make this method public */
        static void SetEnums();

        /** implement the PythonKernelFunction::get_opt_function_values()
         * Returns a tuple of string/string dicts with function name (as in .info.xml) and value for objectives and constraints */
        PyObject* PythonFunctionValues() const;

        /** implement PythonKernelFunction::get_opt_function_properties()
         * @param args name as in .info.xml for objective/constraint/observe */
        PyObject* PythonFunctionProperties(PyObject* args);

        /** Returns all functions. Does not blow up local constraints. Combines objective and constraints.
         * Always creates the list, so use only rarely.
         * @param only_active use only active constraints */
        StdVector<Function*> GetFunctions(bool only_active) const;

        /** return a function (objective and constraint/observe) by name */
        Function* GetFunction(const std::string& name, bool throw_exception = true);

        Tune* SearchTune(Tune::Usage usage, bool silent = false);

        /** Our base ParamNode pointer, pointing to input <optimization> */
        PtrParamNode optParamNode;

        /** Our base ParamNode pointer, pointing to a plain <optimization> */
        PtrParamNode optInfoNode;

        /** This is the list of concurrent objective functions.
         * It is guaranteed to have at least one entry.
         * The objectives are almost identical (share most attributes) */
        ObjectiveContainer objectives;

        /** This contains the constraints in their three forms: standard, hidden and observe
         * and several virtual views on that */
        ConditionContainer constraints;

        /** The current context with our optimization. Very important for meta excitation as it holds the excitation.
         * A context corresponds to a multi sequence step, e.g. bloch mode optimization with concurrent stiffness control.
         * @see ContextManager */
        static Context* context;

        /** Manages the context. Enables to deal with multi sequence optimization, e.g. for different pdes */
        static ContextManager manager;

        /** this is the list of Tune objects, they are updated after CommitIteration() */
        StdVector<Tune*> tunes;

        /** is called from transientDriver after each time step is finished, to store the solution */
//        virtual void TimeStepCalculated(UInt timeStep, AdjointParameters* adjParams) = 0;

        /** is called from assemble, after the calculation of the right-hand side, to get the rhs without Update from Newmark */
  //      virtual void RhsCalculated(AdjointParameters* adjParams) = 0;

        /** Helper to convert from natural solution/design to application
         * @param DesignElement::DENSITY -> App::MECH, DesignElement::POLARIZATION -> App::ELEC */
//        static App::Type ToApp(DesignElement::Type dt);

        /** Default standard design type (not mass) by PDE */
        DesignElement::Type ToDesign(const SinglePDE* pde) const;

        /** optimizer type */
        Optimizer GetOptimizerType() const { return optimizer_; }

        BaseOptimizer* GetOptimizerInstance() { return baseOptimizer_; }

        /** Get the combination of functions in multi objective case and write beta */
        Objective::MultiObjType GetMOType(double& beta) const { beta = this->multiObjectiveBeta_; return multiObjectiveType_; }

        bool IsMultiObjective() const { return isMultiObjective_; }

        bool CalcObjectiveCalled() { return calcObjIteration_ == this->GetCurrentIteration(); }

        /** This tells the driver to store the last solved problem (gid, ...). Called in
         * CommitIteration(). For PiezoSIMP we can save more often and there this method
         * is overwritten and might do nothing.
         * @param step_val the "label" of the "transient" step. -1 is the integer counter */
        virtual void StoreResults(double step_val = -1.0);

        /** call not later than PostInit2(). Add the property to file and info.xml iteration log.
         * python name: opt_register_log_property() */
        void RegisterAuxLogValue(const std::string& name, const std::string& initial);
        void RegisterAuxLogValue(const std::string& name, double initial) { RegisterAuxLogValue(name, std::to_string(initial)); }

        /** update the value to be used by next CommitIteration(). The name should be already set by RegisterAuxLogValue()
         * python name: opt_set_log_property() */
        void SetAuxLogValue(const std::string& name, const std::string& value);
        void SetAuxLogValue(const std::string& name, double value) { SetAuxLogValue(name, std::to_string(value)); }


        /** Encapsulates file Logging information - independent to info.xml logging */
        class Log
        {
        public:
          /** Sets to meaningful defaults (don not much :) ) */
          Log();

          /** Closes the file */
          ~Log();

          /** @param log_name is interpreted. If allows a file, the logFile is created.
           * @param pn_log pointer to the 'log' element. Might be NULL */
          void Init(Optimization* opt, const std::string& log_name, PtrParamNode pn_log);

          /** append an item to the fileHeader and adds the index to the label */
          void AddToHeader(const std::string& label);

          /** The header of the logFile_, to be overwritten if LogFileLine() is overwritten. CommitIteration()
           * writes this string to logFile_ a the first execution.
           * Don't manipulate manually but use AddToHeader() */
           std::string fileHeader;

           /** if set write the design to the logfile */
           bool design;

           /** shall the ev constraints be written to plot.dat? Not in full bloch case when not in detail mode */
           bool plot_ev;

           /** here ErsatzMaterial::CommitIteration() stores bloch information for plot.dat writing.
            * First entry is bandgap then the ev_min or ev_max info for each ev constraint. First label then the value */
           StdVector<std::tuple<std::string, double> > bloch_info;

           /** if set write the gradient of the design to logfile */
           bool designGradient;

           /** if set write the constraint gradient of the design to the logfile */
           bool designConstraintGradients;

           /** write the gradient norms */
           bool gradNorm;

           /** optional log the iterations and cost value to a file to gnuplot it */
           std::ofstream* file = NULL;
        private:

           /** counter for AddToHeader() */
           int columns_;
        };

        /** Keeps all file logging relevant stuff */
        Log log;

        /** Our MultipleExcitation object - by default disabled. Even if we have potentially more than one
         * "multipleExcitation" element in the xml problem file in the case of multi sequence optimization we have
         * only one MultipleExcitaiton object. However some of the information is stored in the corresponding context
         * @see Optimization::contextManager */
        MultipleExcitation* me = NULL;

        /** the reason we did a user break in DoStopIteration(). Set e.g. by PythonStopOptimization(). If set, then we stop. */
        std::string user_break_reason;
        /** did user break converge? E.g. for PythonStopOptimization() */
        bool user_break_converged = false;

        /** Here we contain our design space. The domain gets a reference to it to perform
         * the ersatz material ansatz */
        DesignSpace* design = NULL;

      protected:
        /** Set up the optimization system e.g. prepare the domain for optimization. called
         * exclusively by CreateInstance() -> don't forget to call PostInit() afterwards! */
        Optimization();



        /** called by CommitIteration(), to be overwritten if additional data should be
         * written, then logFileHeader should also be set. Don't add a new-line here!!.
         * @param iteration a duplicate of the log file output to the info xml file */
        virtual void LogFileLine(std::ofstream* out, PtrParamNode iteration);


        /** Gives back the current frequency for printing. This is not the current frequency
         * in multifrequency case. Not a fast method!
         * @return an empty string in the non-harmonic case */
        virtual std::string GetIterationFrequency() { return "not implemented"; }

        /** Determines if the adjoint problem shall be solved directly after the state problem.
         * Note that we have one adjoint for every adjoint sensitive function.
         * Usually they are separated as the adjoint is only necessary for the gradient and in the line search case we don't need it.
         * However it is necessary to solve the adjoint directly after the state in case of multiple harmonic excitations
         * (each frequency assembles a new system matrix) and in the case of multiple sequences as for each sequence step a
         * new system matrix is assembled (here we do not check if there is actually an adjoint computed for the sequence).
         * In the harmonic magnetic coupling we have multiple harmonic extiations of the same frequencies, however we need all forward states first, hence return false */
        bool DoSolveAdjointWithState() const;

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
        unsigned int currentIteration;

        /** This checks how often the state problem is solved. This is not necessary equal
         * to the iterations, e.g. for line search of an external optimizer. Incremented by
         * SolveStateProblem(). */
        int problemSolvedCounter;

        /** The number of problems solved in this iteration */
        int problemWithinIteration;

        /** the maximum number of iterations, have a default */
        int maxIterations;

        /** The actual kind of optimizer.  */
        Optimizer optimizer_;

        /** Here we keep the last iterations design space */
        Vector<double>  last_iteration;

        /** is the first step static */
        bool firstStepStatic;

        /** if the first step is static, this weight specifies how much the other steps add to the functional */
        double otherStepWeight;

        /** shortcut to domain->GetGrid() */
        Grid* grid = NULL;

        /** This holds our optimizer instance. */
        BaseOptimizer* baseOptimizer_ = NULL;

      private:
        /** CommitIteration() does not necessary store the results when we have a stride
         * set in the <commit> element. In case this method makes a StoreResults (not commit)
         * if necessary. Should always be called when finished (in the good and bad sense) */
        void FinalizeStoreResults();

        /** counts the written steps, which can be higher than currentIteration if adjoints or multiples are written */
        int writeCounter_;

        /** When did we store the last result via CommitIteration() due to stride */
        int lastStoredResult_;

        /** checks the time of the iterations to be written and used as stopping criteria.
         * This is a shortcut the the main cfs timer, taken from infoNode */
        shared_ptr<Timer> cfs_timer_;

        /** The time in seconds when. Set in CommitIteration() */
        StdVector<double> time_;

        /** If we have a multiObjective, its type */
        Objective::MultiObjType multiObjectiveType_ = Function::WEIGHTED_SUM;

        bool isMultiObjective_;

        /** If we have a multiObjective and smoothMin or smoothMax, the corresponding beta (default = 1) */
        double multiObjectiveBeta_;

        /** In DesignElement::GetPlainCostGradient we need to know, if the objective values have been set */
        int calcObjIteration_ = -1;

        /** Auxiliary logging information, e.g. from Python for file logging and info.xml
         * first is name, second value. The current value is used when writing the current log line
         * @see RegisterAuxLogValue()
         * @see SetAuxLogValue() */
        StdVector<std::pair<std::string, std::string> > aux_log;
  };

} // namespace


#endif /*OPTIMIZATION_HH_*/
