#ifndef SIMP_HH_
#define SIMP_HH_

#include "Optimization/Optimization.hh"

namespace CoupledField
{
  class MechPDE;
  class StdPDE;
  class Condition;
  template <class TYPE> class StdVector;
  template <class TYPE> class Vector;  
  template <class TYPE> class Matrix;
  
  /** Holds a SIMP (Solid Isotropic Material with Penaltization) optimization.
   *  Actually holds the elements of region to optimize, its densities and
   *  global parameters. Also the cost function and does the looping.
   *  The Reference is Bendsoe, Sigmund; Topology Optimization; Springer Verlag; 2003.
   *  All page numbers refer to this issue. */
  class SIMP : public Optimization
  {
      public:
         /** Up to now w/o parameters */
         SIMP();
        
         /** e.g. closing the ersatzMaterialFile */
         ~SIMP();

         /** compute the value of the cost function */
         virtual double CalcObjective();
         
         /** Evaluates the gradient of the cost function. Saves the result to data.objective_gradient  
          * @param grad_out if not NULL copy write the (design space size) results. 
          * If filtering is enabled, this is automatically filtered. */
         virtual void CalcObjectiveGradient(double* grad_out);

         /** Calculates the constraint(s) */
         double CalcConstraint(Condition* constraint = NULL)
         {
            return CalcConstraint(constraint, false); // no gradient
         }
         
         /** The jacobian of the gradient here as a vector with only one constraint! */
         void CalcConstraintGradient(Condition* constraint = NULL, double* grad_out = NULL)
         {
            CalcConstraint(constraint, true, grad_out);
         }

         /** The systems we can SIMP */
         typedef enum { PIEZO, MECHANIC } System;

         /** Here we also write the density files */ 
         void CommitIteration();

         /** Here we store the system enum */
         static Enum system;
         
         System GetSystem() const { return system_; }

      protected:
         /** switches to the proper constraint, also for gradient case.
          * @param design if not gradient ignored
          * @param grad_out only for gradient and even then optional if not for extern optimizer
          * @return not defined in the gradient case */
         double CalcConstraint(Condition* constraint, bool gradient, double* grad_out = NULL);


         /** <p>Set the element stiffness matrix which is constant for all elements.
          * This method is not only for mechanical SIMP but is also used by PiezoSIMP,
          * therefore it is generic.</p>
          * <p>The code calls for the first element of the reagion the CalcElementMatrix()
          * method of the integrator (form). Thefore temporarily the penalties set to 1</p> 
          * @param pde1 either mech in the PiezoSIMP case elec
          * @param pde2 in mech SIMP same as pde1! only for koupled PiezoSIMP elec 
          *        in the coupled case the order must match what CFS is doing internaly
          * @param out here the element stiffness matrix written. e.h. K_uu which is \int B E B */
         void SetElementStiffness(StdPDE* pde1, StdPDE* pde2, Matrix<double>& out);

         /** called from Optimizer() in each loop after the state problem is called */
         void ExecuteOptimizationStep();
   
         /** implementation of this mechanical cost function */
         double CalcCompliance();

         /** In the implementation based on CalcUKU() as CalcCompliance() */
         void CalcComplianceGradient(double* grad_out);

         /** In the case of Optimality Condition we log the number of lambda iterations */
         void LogFileLine(std::ofstream* out);

         /** The system we optimize. PIEZO or MECHANIC */
         System system_;

         /** We have always a MechPDE for SIMP */
         MechPDE*           mech;

         /** The mechanical element stiffness matrix is constant */
         Matrix<double> mechStiffness;

         /** The region to optimize */
         RegionIdType       regionId;
         
      private:
      
         /** Creates the ErsatzMaterialFile and writes the header */
         void CreateErsatzMaterialFile(const std::string& filename, StdVector<ParamNode*>& des, StdVector<ParamNode*>& tfs);
      
         /** Gets the elasticity tensor from the mechPDE - as long as
          * we are isotrop and linear there is no element id. */
         void GetElasticityTensor(Matrix<double>& matrix_out);
      
         /** Determins the next rho by finding the proper lambda */
         void CalcNextDensity();

         /** The compliance is the sum over $\rpo^p * u^T * K * u$, the gradient
          * is $p * \rpo^{p-1} * u^T * K * u$ over all elements. This is calculated
          * together in this common method. The method gets its own transfer function (MECH, DENSITY)
          * @param derivite decides which method to use from the transfer function
          *        if derivative the the objective derivative is set in the design 
          *        elements */
         double CalcUKU(bool derivative); 
         
         /** Handles the Volume constraint. Has a constraint and constraint derivative mode 
          * @param derivative if false the return value is calculated. Otherwise the value in
          *                   the design element is set. Optionally also grad_out
          * @param grad_out if derivative is set and grad_out is not null it is set.
          * @return invalid in derivative case*/ 
         double CalcVolume(bool derivative, Condition* constraint, double* grad_out);            

         /** Calculates the Greyness OR gauss-greyness! and the derivative of the (gauss) greyness.
          * @param derivative if false the return value is calculated. Otherwise the value in
          *                   the design element is set. Optionally also grad_out
          * @param grad_out if derivative is set and grad_out is not null it is set.
          * @return invalid in derivative case*/
         double CalcGreyness(bool derivative, Condition* constraint, double* grad_out); 

         /** Calculates ans sets the temporary densities according to (2) in the
          *  "99-lines paper".
          *  The goal is to find the lambda s.th. the (volume) constraint holds. This
          *  method does the calculation for a given lambda.
          *  !!! assumes the objective_gradient values to be set!!! 
          *  @param lambda this parameter is to be determined */
         void OptimalityConditionStep(double lambda);

         /** flag indicating if we write the ersatz material file. see destructor */
         std::ofstream*     ersatzMaterialFile;  

         /** This lambda is optimized in CalcNextDensity() and serves as a starting value */
         double             lambda_; 

         /** The number of iterations required to calculate lambda to be used in LogFileLine() */
         int                lambda_iterations_;
         
         /** Here we store a temporary design space für OptimalityCondition AND move calculation! */
         Vector<double>  design_tmp_;
         
  };

} // namespace


#endif /*SIMP_HH_*/
