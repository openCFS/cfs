#ifndef OPTIMIZATION_HH_
#define OPTIMIZATION_HH_

#include <stddef.h>
#include <iosfwd>
#include <string>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/Enum.hh"
#include "General/defs.hh"
#include "MatVec/Vector.hh"
#include "Optimization/Condition.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Objective.hh"
#include "Optimization/OptimizationMaterial.hh"
#include "Utils/StdVector.hh"

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

   // FIXME: this is originally from timestepping.hh and has to be replaced
   typedef enum {NO_DERIVTYPE = 0, FIRST_DERIV = 1, SECOND_DERIV = 2} DERIVType;

   /** This is a simple class used as a parameter to SetAdjointRhs called by assemble */
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
   };
   
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

         /** just queries optimization/ersatzMaterial/material */
         static OptimizationMaterial::System ParseSystem();

         /** PostInit is to be called after the constructor. */
         virtual void PostInit();

         /** This is the second phase of post initialization. It creates the Optimizer tools */
         virtual void PostInitSecond();


         /** Solves the problem by looping over driver->SolveProblem() up to
          * minimum or max iterations is reached. Overwrite on demand*/
         void SolveProblem();

         /** Where we apply the transformation.
          * A subset of the values are PDE identifiers for ToPDE() and ToApp().
          * The heat and acoustic transfer functions are Laplace! */
         typedef enum { MECH, ELEC, PIEZO_COUPLING, PRESSURE, CHARGE_DENSITY, MASS, HEAT, ACOUSTIC, LAPLACE, STRESS, LBM, NO_APP} Application;

         /** Not the optimization problem but the solver! */
         typedef enum { OPTIMALITY_CONDITION, IPOPT_SOLVER, SCPIP_SOLVER, SNOPT_SOLVER, KNITRO_SOLVER,
                        FEAS_PP_SOLVER, SHAPE_SOLVER, EVALUATE_INITIAL_DESIGN, GRADIENT_CHECK  } Optimizer;

         /** to convert string/enum for this type */
         static Enum<Optimizer> optimizer;

         static Enum<Application> application;

         /** the commit mode defines what of the iterations is to be written to gid, ... */
         typedef enum { FORWARD, ADJOINT, BOTH, EACH_FORWARD, EACH_ADJOINT } CommitMode;

         static Enum<CommitMode> commitMode;
         

         /** The DesignSpace element is a container for the complex DesignElements.
          * There are service methods in the container to exchange with external
          * design space arrays as for external optimizers */
        DesignSpace* GetDesign() { return design; }

        /** A cost function is (here) always a scalar type. The value is recalculated when needed!
         * This also stores the value in history of the CostFunction */
        virtual double CalcObjective();

        /** Evaluates the cost-function gradient w.r.t. the design space. Apply SetDesignSpace() first!
         * Writes to DesingElement.objective_gradient.
         * Does a multiplication with Excitation::GetWeightedFactor(). Note, that 
         * the Calc* methods for the objective do only Excitation::GetFactor().
         * @param grad_out size is GetDesignSpaceSize(). If null only DesingElement.objective_gradient */
        virtual void CalcObjectiveGradient(StdVector<double>* grad_out);

        /** Determines the constraint.
         * The function value is stored in value_
         * @param which constraint to calc? Default is the only one! */
        virtual double CalcConstraint(Condition* constraint);

        /** evaluates the gradient of the cost function by the desing element
         * Writes to DesignElement.contraint_gradient
         * @see CalcObjectiveGradient() for parameter description */
        virtual void CalcConstraintGradient(Condition* constraint = NULL, StdVector<double>* grad_out = NULL);

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
        
        /** Solves the Adjoint problem, for given excite and objective
         * This does the real work
         * @param excite multi-excitation
         * @param cost multi-objective */
        virtual void SolveAdjointProblem(Excitation* excite, Function* f);
        
        /** Traverses all objective and active constraint functions (non local) and calls SolveAdjointProblem()
         * if the functions needs it
         * @see Function::IsAdjointBased() */
        virtual void SolveAdjointProblems(Excitation* excite = NULL);
        
        /** Sets the rhs for the adjoint, called by assemle */
        virtual void SetAdjointRhs(AdjointParameters* adjointParams) = 0;

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
         * Checks the stopping rule from the XML file an searches for an HALTOPT file.
         * Shall be called after CommitIteration() ! */
        virtual bool DoStopOptimization();

        /** Do we have a harmonic problem? Then we are complex. Even if not, we might be eigenvalue and also complex*/
        bool IsHarmonic() const { return harmonic_; }
        
        /** we are complex in the harmonic or eigenvalue case */
        bool IsComplex() const { return complex_; }

        /** do we solve eigenvalue problems? Then we are complex! */
        bool IsEigenvalue() const { return eigenvalue_; }

        /** are we in transient optimization? */
        static bool IsTransient();
        
        /** in transient, first step can be static, so that start displacement can depend on material parameters */
        bool IsFirstTransientStepStatic() const {return firstStepStatic; };
        
        /** in transient, first step can be static and have a different weight for the objective, this is the weight
         * converted so that it just has to be multiplied */
        double GetStepWeight(unsigned int ts) const;
        
        /** The current multiple excitation state -> check with IsEnabled() */
        MultipleExcitation* GetMultipleExcitation() const { return me; }

        /** set the (static) enums - if they are used outside optimization, make this method public */
        static void SetEnums();

        /** Returns all active functions. Does not blow up local constraints. Combines objective and constraints.active.
         * Always creates the list, so use only rarely. */
        StdVector<Function*> GetActiveFunctions() const;

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

        /** The applied excitation */
        Excitation* applied_excitation;
        
        /** is called from transientDriver after each time step is finished, to store the solution */
        virtual void TimeStepCalculated(UInt timeStep, AdjointParameters* adjParams) = 0;
        
        /** is called from assemble, after the calculation of the right-hand side, to get the rhs without Update from Newmark */
        virtual void RhsCalculated(AdjointParameters* adjParams) = 0;
        
        Assemble* GetAssemble() { return assemble_; }

        /** Helper to convert from natural solution/design to application
         * @param DesignElement::DENSITY -> MECH, DesignElement::POLARIZATION -> ELEC */
        static Application ToApp(DesignElement::Type dt);

        /** Default standard design type (not mass) by PDE */
        DesignElement::Type ToDesign(const SinglePDE* pde) const;

        /** Helper that converts from mechPDE to MECH and elecPDE to ELEC, ...
         * @param from heat and acoustic the application for the transfer function is laplace, this is indicated by the flag if
         *        we do not want a marker for the pde but the transfer function. Sorry, very messy !! :((
         * @throws if neither mechPDE nor elecPDE
         * @see ToPDE()
         * @see SetPDEs() */
        Application ToApp(const SinglePDE* pde) const;

        /** Find our PDE in SIMP by application from the pdes map
         * @see ToApp()*/
        SinglePDE* ToPDE(Application app, bool throw_exception = true) const;

        /** This is to be overwritten for any case there are other PDEs in ErsatzMaterial::pdes to be set.
         * PiezoSIMP does it simply in the constructor */
        virtual void SetPDEs(OptimizationMaterial::System sys);

        /** Get the standard integrators */
        BiLinFormContext* GetBiLinForm(const RegionIdType reg, Application app1, Application app2 = NO_APP, bool throw_exception = true);

        /** optimizer type */
        Optimizer GetOptimizerType() const { return optimizer_; }

        /** Encapsulates Logging information */
        class Log
        {
        public:
          /** Sets to meaningful defaults (don not much :) ) */
          Log();

          /** Closes the file */
          ~Log();

          /** @param log_name is interpreted. If allows a file, the logFile is created.
           * @param pn_log pointer to the 'log' element. Might be NULL */
          void Init(const std::string& log_name, PtrParamNode pn_log);

          /** append an item to the fileHeader and adds the index to the label */
          void AddToHeader(std::string label);

          /** The header of the logFile_, to be overwritten if LogFileLine() is overwritten. CommitIteration()
           * writes this string to logFile_ a the first execution.
           * Don't manipulate manually but use AddToHeader() */
           std::string fileHeader;

           /** if set write the design to the logfile */
           bool design;

           /** if set write the gradient of the design to logfile */
           bool designGradient;

           /** if set write the constraint gradient of the design to the logfile */
           bool designConstraintGradients;

           /** write mean and max for local constraints */
           bool localDetail;

           /** write the gradient norms */
           bool gradNorm;

           /** optional log the iterations and cost value to a file to gnuplot it */
           std::ofstream* file;
        private:

           /** counter for AddToHeader() */
           int columns_;
        };

        /** Keeps all logging relevant stuff */
        Log log;

        /** The order of the pdes is not defined, Therefore we use the map
         * @see ToApp()
         * @see ToPDE() */
        std::map<Application, SinglePDE*> pdes;

        /** This is simple one SinglePDE from pdes. */
        SinglePDE* pde;

        /** Our MultipleExcitation objecte - by default disabled */
        MultipleExcitation* me;


      protected:
        /** Set up the optimization system e.g. prepare the domain for optimization. called
         * exclusively by CreateInstance() -> don't forget to call PostInit() afterwards! */
        Optimization();

        /** Evaluates objective and constraint function and gradient.
         * Implemented in ErsatzMaterial
         * @param grad_out only used in derivative case
         * @return zero for derivative */
        virtual double CalcFunction(Excitation& excite, Function* f, bool derivative) = 0;

        /** This tells the driver to store the last solved problem (gid, ...). Called in
         * CommitIteration(). For PiezoSIMP we can save more often and there this method
         * is overwritten and might do nothing.
         * @param step_val the "label" of the "transient" step. -1 is the integer counter */
        virtual void StoreResults(double step_val = -1.0);

        /** called by CommitIteration(), to be overwritten if additional data should be
         * written, then logFileHeader should also be set. Don't add a new-line here!!.
         * @param iteration a duplicate of the log file output to the info xml file */
        virtual void LogFileLine(std::ofstream* out, PtrParamNode iteration);


        /** Gives back the current frequency for printing. This is not the current frequency
         * in multifrequency case. Not a fast method!
         * @return an empty string in the non-harmonic case */
        virtual std::string GetIterationFrequency() { return "not implemented"; }

        /** The assemble class for our PDE */
        Assemble* assemble_;

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

        /** The actual kind of optimizer.  */
        Optimizer optimizer_;

        /** Here we contain our design space. The domain gets a reference to it to perform
         * the ersatz material ansatz */
        DesignSpace* design;

        /** Here we keep the last iterations design space */
        Vector<double>  last_iteration;

        /** are we harmonic/EV or static/transient? */
        bool complex_;
        
        /** only for the driver, not for complex_! */
        bool harmonic_;

        /** do we solve an eigenvalue problem. Includes block mode problems */
        bool eigenvalue_;

        /** bloch mode analyis is also eigenvalue but special due to the wave vectors encapsualted in excitations */
        bool bloch_;

        /** is the first step static */
        bool firstStepStatic;
        
        /** if the first step is static, this weight specifies how much the other steps add to the functional */
        double otherStepWeight;

        /** shortcut to domain->GetGrid() */
        Grid* grid;

      private:
        /** CommitIteration() does not necessary store the results when we have a stride
         * set in the <commit> element. In case this method makes a StoreResults (not commit)
         * if necessary. Should always be called when finished (in the good and bad sense) */
        void FinalizeStoreResults();

        /** counts the written steps, which can be higher than currentIteration if adjoints or multiples are written */
        int writeCounter_;

        /** When did we store the last result via CommitIteration() due to stride */
        int lastStoredResult_;

        /** This holds our optimer instance. */
        BaseOptimizer* baseOptimizer_;
  };

} // namespace


#endif /*OPTIMIZATION_HH_*/

