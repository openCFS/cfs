#ifndef OPTIMALITYCONDITION_HH_
#define OPTIMALITYCONDITION_HH_

#include <iosfwd>
#include <string>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/Enum.hh"
#include "MatVec/Vector.hh"
#include "Optimization/Optimizer/BaseOptimizer.hh"

namespace CoupledField {
class Optimization;
}  // namespace CoupledField

namespace CoupledField
{
  /** This is an implementation of Optimality Condition for SIMP by Bendsoe and Sigmund.
   * It is such far generalized that it handles also negative gradients by setting them to 0.
   * We ignore any scaling! */
  class OptimalityCondition : public BaseOptimizer
  {
    public:
      /** @param optimization the problem we optimize
       * @param pn here we can have options - might be NULL! */
      OptimalityCondition(Optimization* optimization, PtrParamNode pn);
      
      /** This destructor does nothing but Optimizer has a virtual destructor */
      ~OptimalityCondition() { } 
      
      virtual void LogFileHeader(Optimization::Log& log);

      virtual void LogFileLine(std::ofstream* out, PtrParamNode iteration);

    protected:

      /** Solves the problem. All stuff, including evaluations of the state problem is done
       * within this method. 
       * @throws exception when not ok! */ 
      void SolveProblem();
      
    private:
      /** This are our both bisection types. Framed is what Sigmund did, Trajectory is own breed but
       * most likely not unique */
      typedef enum { FRAMED, FUMBLE, TRAJECTORY, EXTREMIZE } Type;
      
      /** Determines the next rho by finding the proper lambda with the framed bisection.
       * @return if delta_err was too small (and checked) - indicates too small move limit and corrects bound in that case  */
      bool CalcNextFramedIteration(bool last_was_stalled_err);

      /** This shall be a quite stable but expensive bisection variant. It is based on a lambda
       * and a step. We further have a contract factor (e.g. 0.49) and an expand factor (e.g. 1.99).
       * The we check the error for the four variants of lambda +/- [contract/expand] * step.
       * The least error determines the new step variable. This shall help for optimizations 
       * where lambda might move around the 0. The odd values for the contract and expand factors
       * shall prevent lambda to become 0 as we have to divide by lambda. This is also checked in
       * Evaluate() */
      void CalcNextFumbleIteration();
      
      /** Do a bisection iteration with a free moving lambda */
      void CalcNextTrajectoryIteration();
      
      
      /** For simple plate problems with mixed signed gradiens w/o constraints. Follows the
       * gradient by move_limit w/o any sorting! */
      void CalcNextExtremizeIteration();
      
      /** Calculates ans sets the temporary densities according to (2) in the
       * "99-lines paper".
       * The goal is to find the lambda s.th. the (volume) constraint holds. This
       * method does the calculation for a given lambda.
       * !!! assumes the objective_gradient values to be set!!!
       * Checks for lambda smaller lambda_min! 
       * @param lambda this parameter is to be determined
       * @return the error assuming a single volume constraint */
      double Evaluate(double lambda);

      /** What type of bisection do we do? It becomes more subtle when not doing compliance! */
      Type type_;
      
      /** reads the type from xml */
      Enum<Type> type;
      
      /** this is the damping factor of the objective gradient for the optimality condtion.
       * The power is taken - default is 0.5 which is the square root */
      double oc_damping_;
      
      /** minimal densitiy variation (see book (1.12)) = move limit for optimality condition */
      double move_limit_;       
      
      /** This lambda is optimized in CalcNextDensity() and serves as a starting value */
      double lambda_; 

      /** The maximum number of lambda iterations */
      int max_lambda_iters_;
      
      /** The number of lambda_iterations for this iteration only */
      int lambda_iters_;
      
      /** lambda_min from xml. Might not become smaller than this value. in trajectory we
       * switch sign then. */
      double lambda_min_;
      
      /** err_eps for bisection iteration is fine. This is actually the feasibility tolerance */
      double feasibility_;
      
      /** for the framed mode the starting lower border */
      double start_lower_;
      double start_upper_;
      
      /** for the framed mode the enlargement when we start iteration or exceeded max_iters. 
       * This are factors (lambda_ * enlarge_x_)*/
      double enlarge_lower_;
      double enlarge_upper_;

      /** For the framed type this means if we enlarge the borders in every iteration */
      bool always_enlarge_;

      /** this are the borders for the framed bisection */
      double lower_;
      double upper_;

      /** triggers framed check for stalled errors */
      bool check_stalled_err_ = true;
      
      /** For fumble this is the current step */
      double step_;
      
      /** @see CalcNextFumbleIteration() */
      double contract_;
      
      /** @see CalcNextFumbleIteration() */
      double expand_;
      
      /** Here we save a design space such that we can play with different lambdas  */
      Vector<double> vault_;
      
      /** This is a temporay design space such that it can be used in evalue */
      Vector<double> evaluate_tmp_;
  };



} // end of namespace

#endif /*OPTIMALITYCONDITION_HH_*/
