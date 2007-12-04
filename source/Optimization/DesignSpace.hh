#ifndef DESIGN_SPACE_HH_
#define DESIGN_SPACE_HH_

#include "Optimization/TransferFunction.hh"
#include "Optimization/DesignElement.hh"
#include "Forms/baseForm.hh"
#include "General/Enum.hh"

// we need it for the template implementation 
#include "Domain/resultInfo.hh"
#include "Utils/result.hh"

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
     void ExtractResults(shared_ptr<BaseResult> result, bool isComplex)
     {
       if(isComplex) ExtractResults<std::complex<double> >(result);
                else ExtractResults<double>(result);
     }
     
     /** <p>Take the design space from an external optimizer and apply it on the problem.
      * Do not solve the problem, but call this method before CalcCostFunction() :)</p>
      * <p>The background is, that external optimizers might have a own design space array
      * (only doubles) where we have a StdVector of the complex DesignElement.</p> 
      * @param space_in the design space (in variable). Size is GetDesignSpaceSize() 
      * @return the design_id which is the old one if space_in did not change the design. */
     int ReadDesignFromExtern(const double* space_in);        
           
     /** gives the initial guess (for the design space) 
      * @param space_out to this array of GetDesignSpaceSize() the initial guess is wrtitten to.
      * @return the internal design_id as calculated by ReadDesignFromExtern()
      * @see SetDesignSpace() */
     int WriteDesignToExtern(double* space_out) const;      

     /** Similar but more general as WriteDesignToExtern()
      * @param de if NO_TYPE or DEFAULT we take the value and do not check!*/ 
     void WriteGradientToExtern(double* out, DesignElement::Type de, 
         DesignElement::ValueSpecifier vs, DesignElement::Access access) const;
     
     /** Sets the value of the described design element to 0 
      * @param design with design elements to set
      * @param vs what values to set. Not all make sense -> exception */
     void Reset(DesignElement::Type design, DesignElement::ValueSpecifier vs);
     
     /** This disables the transfer functions -> sets them to NO_TYPE. This is used
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

     /** Here we store result descriptions as defined in DesignElement.hh */
     StdVector<ResultDescription> resultDescriptions; 
     
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
     int GetSpecialResultIndex(DesignElement::Type design, DesignElement::ValueSpecifier value, 
                               DesignElement::Detail detail, DesignElement::Access access = DesignElement::SMART);    
   private:
     template <class T>
     void ExtractResults(shared_ptr<BaseResult> result);     
     
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

     /** This number indentifies the design space. It is always incremented if ReadDesignFromExtern() reads
      * a different design */
     int design_id;
     
     /** This transforms FormName (Integrator) to Application -> so the enum is actually Application 
      * but here int because of circular includes :( */
     Enum<int> applicationForm;
  };  

  
  /** somehow this has to be in the header file, otherwise it doesn't link :( */
  template <class T>
  void DesignSpace::ExtractResults(shared_ptr<BaseResult> base_result)
  {
    // our results are up to now scalar!
    Result<T>& result = dynamic_cast<Result<T> &>(*base_result);
    Vector<T>& result_data = result.GetVector();
    
    // set the result as we need it
    result_data.Resize(elements_);
    
    // the description of the result
    shared_ptr<ResultInfo> ri = result.GetResultInfo();
      
    // Work with a result description. This is either a result description from the
    // xml file or when using the "predefined" *_PSEUDO_* we set it here.
    ResultDescription def;
    // set the defaults to be maybe replaced by a resultDescription
    def.solutionType = ri->resultType;
    // this is clearly nonsense if the result/solution tyoe is OPT_RESULT_1/2/3
    def.design = ri->resultType == MECH_PSEUDO_DENSITY ? DesignElement::DENSITY : DesignElement::POLARIZATION;
    def.access = DesignElement::PLAIN;
    def.value  = DesignElement::DESIGN;
    
    ResultDescription& descr = def;
    // ignore defaults if there is a result description for the OPT_RESULT_1/2/3 case  
    for(unsigned int i = 0; i < resultDescriptions.GetSize(); i++)
      if(resultDescriptions[i].solutionType == ri->resultType) 
        descr = resultDescriptions[i];
    
    // search where in data we are
    int base = FindDesign(descr.design);

    // loop over elements from result. We have to do it this way as the the connection
    // of design element and result element is the element(->elemeNum) but we cannot
    // search in the result for an element.
    EntityIterator it = result.GetEntityList()->GetIterator();
    for ( it.Begin(); !it.IsEnd(); it++ ) 
    {
      // for elements not in the design region we set to one
      double result_value = 1.0;
       
      if(it.GetElem()->regionId == regionId_)
      {
        // note that the index is from the first design set!
        unsigned int base_index = Find(it.GetElem()->elemNum);
        // base=0 is first!
        unsigned int data_index = (base * elements_) + base_index; 
        DesignElement& de = data[data_index];
        result_value = de.GetValue(&descr);      

        #ifdef CHECK_INDEX
          if(de.elem->elemNum != it.GetElem()->elemNum) 
            EXCEPTION("mixed up indices:" << de.elem->elemNum << "!=" << it.GetElem()->elemNum
                << " base_index=" << base_index << " data_index=" << data_index << " it.Pos()=" << it.GetPos());
        #endif      
      }
      result_data[it.GetPos()] = result_value;
    }
  }
  
  
} // end of namespace

#endif /*DESIGN_SPACE_HH_*/

