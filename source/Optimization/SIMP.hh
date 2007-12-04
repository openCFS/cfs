#ifndef SIMP_HH_
#define SIMP_HH_

#include "Optimization/Optimization.hh"
#include "Domain/bcs.hh"
#include "Utils/cfsvector.hh"
#include "OLAS/matvec/vector.hh"
#include <map>

namespace CoupledField
{
  class StdPDE;
  class SinglePDE;
  class MechPDE;
  class BaseForm;
  class BiLinFormContext;
  class BaseMaterial;
  class Condition;
  class Assemble;
  class TransferFunction;
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
         virtual ~SIMP();

         /** compute the value of the cost function */
         virtual double CalcObjective()
         {
           return harmonic ? CalcObjective<std::complex<double> >() 
                           : CalcObjective<double>();
         }
         
         /** with the output objective we also have to solve the adjoint problem. */
         virtual void SolveStateProblem(); 
         
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
         static Enum<System> system;
         
         System GetSystem() const { return system_; }

         /** Adds validation stuff here to keep out of long constructor */
         void PostInit();
         
         /** if we do radiation optimization this gives the bilinear form which is 
          * added in DefineIntegrators of MechPDE early, before we construct Optimization.
          * It would clearly be nicer if one could add this stuff in SIMP::PostInit()
          * by making OLAS more flexible. Actually this is a low priority todo
          * @return NULL if we do not do such "radiation" optimization stuff */ 
         static BiLinFormContext* CreateSurfaceNormalMatrix(SinglePDE* pde, BaseMaterial* baseMat, shared_ptr<ResultInfo> result);
         
      protected:
         /** switches to the proper constraint, also for gradient case.
          * @param design if not gradient ignored
          * @param grad_out only for gradient and even then optional if not for extern optimizer
          * @return not defined in the gradient case */
         double CalcConstraint(Condition* constraint, bool gradient, double* grad_out = NULL);

         /** Helper which extracts the Form from assemble using the optimization region
          * @param pde1 the first pde (e.g. mech)
          * @param pde2 this is either the same as pde1 or the coupling partner
          * @param integrator there is no nice enum yet :( e.g. linElastInt, MechInt, ... */ 
         BaseForm* GetForm(StdPDE* pde1, StdPDE* pde2, const std::string& integrator);
         
         /** <p>Get the original element matrix (stiffness, mass, ...) 
          * which is constant for all isotripic elements.
          * This method is not only for mechanical SIMP but is also used by PiezoSIMP,
          * therefore it is generic.</p>
          * <p>If no elemen is given, the one from the first design element is used.</p>
          * <p>All transfer functions are disabled during this method. Call only for
          * enabled transfer functions (default)</p> 
          * @param form to be extracted via GetForm()  
          * @param out here the element stiffness matrix written. e.h. K_uu which is \int B E B 
          * @param elem if not given the first design element is used, otherwise the provided one*/
         void GetElementMatrix(BaseForm* form, Matrix<double>& out, Elem* elem = NULL);

         /** This class holds the solution of the PDE. It is in a class such that it 
          * helps to encapulate real and complex solutions. Note that the Piezo 
          * has other solutions! */
         class Solution
         {
           public:
             Solution(SIMP* simp);
             
             ~Solution();
             
             typedef enum { ELEMENT_VECTORS, RAW_VECTOR } StorageType;
             
             /** Copies the solution for the pde in our own storage.
              * In the ELEMENT_VECTORS case make sure, that the solution is in the PDE!
              * For manual adjoint stuff you might have to do SaveSolution() first!
              * @param st if we copy the vector as RAW_VECTOR or element wise
              * @param pde will me mech in SIMP and also elec in PiezoSIMP
              * @param app redundant to pde. Either MECH or ELEC */
             void ReadSolution(StorageType st, StdPDE* pde, Application app)
             {
               if(simp->harmonic)
                 ReadGenericSolution<std::complex<double> >(st, pde, app);
               else
                 ReadGenericSolution<double>(st, pde, app);
             }
             
             /** Writes the solution (raw vector) back to the pde */
             void WriteSolution(StdPDE* pde, Application app)
             {
               if(simp->harmonic)
                 WriteGenericSolution<std::complex<double> >(pde, app);
               else
                 WriteGenericSolution<double>(pde, app);
             }
             
             /** This is an element wise storage of the solution 
              * the Application shall be MECH or ELEC */
             std::map<Application, StdVector<CFSVector* > > elem;

             /** This is the algsys solution vector */
             std::map<Application, CFSVector* > raw;
             
           private:
             template <class T>
             void ReadGenericSolution(StorageType st, StdPDE* pde, Application app);

             template <class T>
             void WriteGenericSolution(StdPDE* pde, Application app);

             /** Reference to our optimization problem */
             SIMP* simp;
         };
         
         /** This method calculates the compliance and compliance gradient
          * plus the mechanism optimization (output) gradiend for SIMP and PiezoSIMP.<br>
          * The General form is SUM over factor * u1 * ny(x,k) * u2. Here ny(x,k) is
          * the transfer function over the design variable and k. This can be the derivative.
          * In compliance u1 = u2. In mechanism optimization u1 is the lagrange state
          * @param tf the transfer function, e.g. for DENSITY and MECH
          * @param app1 which element vector part of u1 (MECH or ELEC)
          * @param u1 the solution or the lagrange multiplier.
          * @param k a stiffness matrix
          * @param derivite decides which method to use from the transfer function
          *        if derivative the the objective derivative is set/added in the design 
          *        elements
          * @param add if true the derivative is summed up in the design elements and not set
          * @param factor modify the vec mat vec product. */
         double CalcU1KU2(TransferFunction* tf, StdVector<CFSVector*>& u1, 
                                Application k, StdVector<CFSVector*>& u2, 
                                bool derivative, bool add = false, double factor = -1.0)  
         {
           if(harmonic) return CalcU1KU2<std::complex<double> >(tf, u1, k, u2, derivative, add, factor);
                   else return CalcU1KU2<double>(tf, u1, k, u2, derivative, add, factor);
         }
         
         
         
         /** <p>We move the K*u2 part out of CalcU1KU2 as this might be quite complex in the damping case
          * 'K' = (K+i*alpha_k*K+i*alpha_m*M-omega^2*M) and K and M have the transfer function applied to it.
          * It also serves for the Surface Normal Matrix (radiation).</p>
          * <p>We have to wrap the templated version here because an templated method cannot
          * be virtual :(/</p>
          * @param app here 'only' MECH and SURFACE_NORMAL_MATRIX is implemented. Piezo stuff -> PiezoSIMP */ 
         virtual void CalcElementKU2(const CFSVector* in, Application app, DesignElement* de, bool derivative, CFSVector* out) {
           if(harmonic) CalcElementKU2<std::complex<double> >(in, app, de, derivative, out);
                   else CalcElementKU2<double>(in, app, de, derivative, out);
         }

         /** This does a forward/ajoint step for the mechanism optimization. It works for
          * both (mechanical) SIMP and PiezoSIMP.
          * @param elec is NULL in SIMP case */
         void SolveAdjointProblem(StdPDE* mech, StdPDE* elec = NULL)
         {
           if(harmonic)
             SolveAdjointProblem<std::complex<double> >(mech, elec);
           else
             SolveAdjointProblem<double>(mech, elec);
         }
         
         template <class T>
         double CalcObjective();
         
         /** The system we optimize. PIEZO or MECHANIC */
         System system_;

         /** We have always a MechPDE for SIMP */
         MechPDE*           mech;

         /** The mechanical element stiffness matrix is constant */
         Matrix<double> mechStiffness;

         /** The mechanical element mass matrix is also constant. Only for harmonic! */
         Matrix<double> mechMass;
         
         /** The surface normal matrix following Du/Olhoff 2007. Note that this assumes a single plane surface */
         Matrix<double> surfaceNormal;
         
         /** The region to optimize */
         RegionIdType       regionId;

         /** Here we store the solution of the problem */
         Solution* forward_;
         
         /** Here we store the soluton of the adjoint problem. */
         Solution* adjoint_;
         
      private:

        /** Creates the ErsatzMaterialFile and writes the header */
         void CreateErsatzMaterialFile(const std::string& filename, StdVector<ParamNode*>& des, StdVector<ParamNode*>& tfs);
         
         /** Handles the Volume constraint. Has a constraint and constraint derivative mode 
          * @param derivative if false the return value is calculated. Otherwise the value in
          *                   the design element is set. Optionally also grad_out
          * @param grad_out if derivative is set and grad_out is not null it is set.
          * @return invalid in derivative case*/ 
         double CalcVolume(bool derivative, Condition* constraint, double* grad_out);            

         /** See the non-template version for documentation! */
         template <class T>
         double CalcU1KU2(TransferFunction* tf, StdVector<CFSVector*>& u1, 
             Application k, StdVector<CFSVector*>& u2, 
             bool derivative, bool add, double factor = -1.0);

         /** see the non-template version! */
         template <class T>
         void CalcElementKU2(const CFSVector* in, Application app, DesignElement* de, bool derivative, CFSVector* out);
         
         
         /** This solves the forward and adjoint problem and stores all relevant data. Calls SetAndSolveAdjointRHS() */
         template <class T>         
         void SolveAdjointProblem(StdPDE* mech, StdPDE* elec = NULL);
         
         /** Takes care for making CFS solving the adjoint PDE. Sets output_vector_ */
         template <class T>
         void SetAndSolveAdjointRHS();
         
         /** Calculates the Greyness OR gauss-greyness! and the derivative of the (gauss) greyness.
          * @param derivative if false the return value is calculated. Otherwise the value in
          *                   the design element is set. Optionally also grad_out
          * @param grad_out if derivative is set and grad_out is not null it is set.
          * @return invalid in derivative case*/
         double CalcGreyness(bool derivative, Condition* constraint, double* grad_out); 

         /** Calculates the product of the (system) surface normal matrix with the solution already in OLAS.
          * Note that we have to use 1 based OLAS vectors as the sparse system matrix is from OLAS .
          * This calculation is done for the adjoint rhs and also for calculate the radiation objective.
          * It shall be cheap enough to calc here twice! */
         template <class T>
         void CalcSurfaceNormalTimesSolution(OLAS::Vector<T>& olas_prod);         
         
         /** flag indicating if we write the ersatz material file. see destructor */
         std::ofstream*     ersatzMaterialFile;  
         
         /** Here we store a temporary design space for move calculation! */
         Vector<double>  design_tmp_;
         
         /** When we optimize output we store here the nodes */
         LoadList output_nodes_;
         
         /** This is the vector corresponding to the output nodes for <l,u>. To be deleted in destructor! */
         CFSVector* output_vector_;
         
         /** this is the optimization->simp XML element */
         ParamNode* pn;
         
         /** The assemble class for our PDE */
         Assemble* assemble_;
         
         /** This is a math parser handle for CalcElementKU2() */
         unsigned int mathParserHandle_;
  };

} // namespace


#endif /*SIMP_HH_*/
