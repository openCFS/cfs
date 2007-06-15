#ifndef PIEZOSIMP_HH_
#define PIEZOSIMP_HH_

#include "Optimization/SIMP.hh"

namespace CoupledField
{
  class ElecPDE;
    
  /** <p>Extends the (mechanical, compliance minimization) SIMP to the Piezo SIMP as described
   *  by [1]: Kögl M., Silva E.; Topology optmization of smart structures: design of 
   *  piezoelectric plate and shell actuators; 2005 Smart Mater. Struct. 14 387-399</p>
   * <p>The reciprocal theroem is described in [2]: Silva E., Nishiwaki S., Kikuchi N.; 
   * Topology Optimization Design of Fextensional Actuators; 2000 IEEE Transactions on Ultrasonics,
   * Ferroelectrics, And Frequency Control, Vol. 47, No. 3</p>
   * <p>The basic idea is, that there is for each iteration a direct and indirect piezo
   * simulation to compute the transduction and another iteration for the mean compliance.
   * Therefore the activate tag in the bcs and loads is used.</p> */
  class PiezoSIMP : public SIMP
  {
    public:
      PiezoSIMP();
        
      /** We solve in each iteration three state problems and store the results. */
      void SolveStateProblem(); 
        
      /** compute the value of the cost function */
      double CalcObjective();
     
      /** Evaluates the gradient of the cost function. Saves the result to data.objective_gradient  
       * @param grad_out if not NULL copy write the (design space size) results. 
       * If filtering is enabled, this is automatically filtered. */
      void CalcObjectiveGradient(double* grad_out);
     
      /** In the case of transduction we write the transduction in more detail */
      void LogFileLine(std::ofstream* out);
     
      /** storage settings for piezo 
       * COMMIT: then Optimization::CommitIteration() does all which is after commit and
       * not after each problem the case2 (mech_case) is stored.
       * All other cases do a storage in each problem solving case. This is for SCPIP
       * mostly equal to the iteration, in IPOPT there can be much more line searches = 
       * problems within one iteration. The iteration counter is increased after
       * each storage and 
       * ELCE_*/
      typedef enum { COMMIT, ELEC_CASE, MECH_CASE, BOTH_CASES } Storage;

      /** Storage translation */
      static Enum storage;
     
     Storage GetStorage() const { return storage_; };
     
    protected:
      /** This is called every Optimization::CommitIteration(). When we have special storing
       *  rules (-> gid, gmv, ...) this does nothing as we store in SolveStateProblem then */ 
      void StoreResults(double step_val = -1.0); 
     
    private:
     /** Solves a subproblem of the transduction, only then we can call CalcTransduction().
      * Has to be call first for ELEC, then for MECH!*/
     void SolveTransductionSubProblem(Application subProblem);
        
     /** The transduction (to be maximized) gives the coupling between the electical and mechanical forces based
      * on the reciprocal theorem. Two state problems are to be calculated as shown in [1], A good introduction
      * to the reciprocal theorem is in [2]. */
     double CalcTransduction();
        
     /** Implementation according to Kögl and Silva in [1]
      *  @see CalcObjectiveGradient() */   
     void CalcTransductionGradient(double* grad_out);  
        
     /** We have for the direct coupled piezo the MechPDE which is in SIMP.hh and 
     * here the corresponding ElecPDE */
     ElecPDE* elec;

     /** For the evaluation of the mean transduction cost function two state problems have to
      * be solved in each iterations. The results are copied to the <mech|elec>_sol_<1|2> vectors */            
     StdVector<Vector<double> > mech_sol_1;           
     StdVector<Vector<double> > mech_sol_2;
     StdVector<Vector<double> > elec_sol_1; 
     StdVector<Vector<double> > elec_sol_2;
     
     /** The elec stiffness matrix $K_{\phi \phi}$. $K_{uu}$ is already as mechStiffness in SIMP.hh */
     Matrix<double> elecStiffness;

     /** The coupling stiffness matrix $K_{u \phi}$ */
     Matrix<double> coupledStiffness;
     
     /** The transposed coupling stiffness matrix $K_{u \phi}^T$ */         
     Matrix<double> coupledStiffnessTransposed;

     /** This is a simple class for setting the transduction values */
     class Transduction
     {
       public:
         Transduction();
       
         void Set(ParamNode* pn);
      
         double value;
         /** either load or pressure or surface density */
         bool isLoad;
         /** internal use to validate if the property is given in the PDE description */
         bool found; 
     };
 
     /** This is the load or integrator for the electrical part of the transduction */
     Transduction elec_;
     
     /** This is the load or integrator for the mechanical part of the transduction */
     Transduction mech_;
     
     /** This are logging variables for LogFileLine */
     double log_coupled_;
     double log_coupled_simp_;
     double log_elec_;
     double log_elec_simp_;
     
     /** the actual storage setting */
     Storage storage_;
  };

} // end of namespace


#endif /*PIEZOSIMP_HH_*/
