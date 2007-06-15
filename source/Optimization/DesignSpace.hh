#ifndef DESIGN_SPACE_HH_
#define DESIGN_SPACE_HH_

#include "Optimization/TransferFunction.hh"
#include "Optimization/DesignElement.hh"
//#include "Optimization/Optimization.hh"
#include "Forms/baseForm.hh"
#include "General/Enum.hh"
#include "boost/shared_ptr.hpp"


namespace CoupledField
{
  template <class TYPE> class StdVector;
  class SinglePDE;
  class ParamNode;
  class Elem;
  class BaseResult;
  class ResultInfo;
  class RegionTypeId;
  
  /** This is the container of DesingElements which also holds the transferFunctions.
   * It can be initialized by Optimization of can contain the ersatz material stuff. */
  class DesignSpace
  {
    public:
     /** Constructor for SIMP type Optimization - there we lay on a region which containts also n# elements
      * @param result the result description list  */
     DesignSpace(int regionId, StdVector<ParamNode*>& design, StdVector<ParamNode*>& transfer, StdVector<ParamNode*>& result);
    
     /** PostInit as usual when not all can be stuffed into the constructor
      * @param constraints the number of constraints to initialize constraintGradients */
     void PostInit(int constraints); 
    
     /** This gives the ersatz material factor for an element. 
      * This fulfills the trick, that there might be more transfer function for
      * a single element -> as in the coupling term of Piezo SIMP. */
     double GetErsatzMaterialFactor(const Elem* elem, const BaseForm* form);

     /** This is just a variant.
      * @see GetErsatzMaterialFactor(const Elem*, const BaseForm*) for documentation. 
      * @param design_index must be not bigger than elements. */
     double GetErsatzMaterialFactor(unsigned int design_index, Optimization::Application applic);
     
     /** This gets back a uniquely defined transfer function.
      * @param throw_exception if false NULL is returned when nothing is found! */
     TransferFunction* GetTransferFunction(DesignElement::Type design, Optimization::Application application, bool throw_exception = true); 

     /**<p>check the optResult_1/2/3 from the optimization/simp/result elementes against
      * element results in the pde and conditionally add it as store results to the pde.</p>
      * @param pde in the current implementation a MechPDE */ 
     void AppendOptimizationResults(SinglePDE* pde);   
     
     /** <p>Copies the relevant data from the design element to the result such that it can be written
      * to the output (e.g. gid). This is called from the CalcResults() methods of the relevant pdes
      * when the result is defined in the xml file.</p>
      * Note, that the result might come from the 'predefined' MECH_PSEUDO_DENSITY or ELEC_PSEUDO_POLARIZATION
      * or from the OPT_RESULT_1/2/3 where the actual result is defined in a xml element under 
      * optimization/SIMP</p> */
     void ExtractResults(shared_ptr<BaseResult> result);  
     
     /** <p>Take the design space from an external optimizer and apply it on the problem.
      * Do not solve the problem, but call this method before CalcCostFunction() :)</p>
      * <p>The background is, that external optimizers might have a own design space array
      * (only doubles) where we have a StdVector of the complex DesignElement.</p> 
      * @param space_in the design space (in variable). Size is GetDesignSpaceSize() */
     void ReadDesignFromExtern(const double* space_in);        
           
     /** gives the initial guess (for the design space) 
      * @param space_out to this array of GetDesignSpaceSize() the initial guess is wrtitten to.
      * @see SetDesignSpace() */
     void WriteDesignToExtern(double* space_out);      
     
     /** This disables the transfer functions -> sets them to IDENTITY. This is used
      * in SIMP to calculate the original stiffnes matrices.
      * The setting from the XML file is stored -> to be undone with EnableTranferFunctions() */
     void DisableTransferFunctions();
     
     /** Enables the transfer functions -> sets again to the xml settings after
      * temporarily disabled via DisableTransferFunctions() */
     void EnableTransferFunctions();

     /** Service method to find our index in the design vector
      * @return -1 if not throw_exception and not found */
     int FindDesign(DesignElement::Type dt, bool throw_exception = true);

     /** Service method to find a specific design element by element number and design type */
     DesignElement* Find(unsigned int elemNum, DesignElement::Type dt); 

     /** When we have more design types this is a divisor of data.GetSieze() */
     unsigned int GetNumberOfElements() { return elements_; }
     
     /** This is our real design data, a set of DesignElements */
     StdVector<DesignElement> data;    
     
     /** Our transfer functions */
     StdVector<TransferFunction> transfer;
     
     /** Here we store the design we have. Check with Find() for the element */
     StdVector<DesignElement::Type> design;
     
     /** Might be nonsens if our constructor is no simp or ersatz material one! */
     int GetRegionId() { return regionId_; }
     
     /** We can fine define results in the optimization part: <pre>
      * <result id="optResult_2" design="density" access="plain" value="costGradient" detail="pKu" />            
      * <result id="optResult_3" design="density" access="plain" value="objective" detail="pKp" /></pre>                  
      * In the above example the design parameter for objective makes no sense. One shall choose the first.
      * This is for PiezoSIMP to allow to make the transduction (objective, detail=none) and optionally 
      * the intermediate steps (objective, detail=uKp, pKp). Same details are for the transduction
      * 'costGradient, with the details uKu, uKp, pKu, pKp.
      * @return 0, 1, 2 is the index (optResult_x+1) for DesignElement::specialResult_[]. -1 if no
      * result is specified in XML */
     int GetSpecialResultIndex(DesignElement::Type design, DesignElement::ValueSpecifier value, DesignElement::Detail detail);    
   private:

     /** Sets a ResultInfo object by a result type optResult_1/2/3 where the detailed
      * description is in the optimization/SIMP elemenet in the xml file
      * @param solutionType switch to typedef when environment.hh is gone
      * @return null if this solution type was not given in xml or a filled version.
      *         The deletion of the object is the duty of the caller (or shared pointer) */
     ResultInfo* GetResultInfo(ResultDescription& rd);   
   
     /** finds the index of the desing element in desing.data for the element.
      * to be optimized later if needed. Searches only in the first parameters set! */
     unsigned int Find(unsigned int elemNum);
   
     /** for SIMP type constructor we have the region of this DesignSpace
      * -1 is for other constructor */
     int regionId_;
     
     /** for SIMP type constructor we have a number of elments, 
      * data size = num of design * num region elements
      * -1 is for other constructor */
     unsigned int elements_;    

     unsigned int last_find_index_;
   
     /** This transforms FromName to Application -> so the enum is Application */
     Enum applicationForm;
     
     /** Here we store result descriptions as defined in DesignElement.hh */
     StdVector<ResultDescription> resultDescriptions_; 
     
  };  

} // end of namespace

#endif /*DESIGN_SPACE_HH_*/

