#ifndef DESIGN_SPACE_HH_
#define DESIGN_SPACE_HH_

#include "Optimization/TransferFunction.hh"
#include "Optimization/DesignElement.hh"
#include "Optimization/ErsatzMaterial.hh"
#include "Optimization/DesignMaterial.hh"
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
  class InfoNode;
  class Elem;
  class BaseResult;
  class ResultInfo;
  class RegionTypeId;
  class BaseOptimizer;

  /** This is the container of DesingElements which also holds the transferFunctions.
   * It can be initialized by Optimization of can contain the ersatz material stuff. */
  class DesignSpace
  {
    public:
     /** Constructor for SIMP type Optimization - there we lay on a region which containts also n# elements
      * @param result the result description list  */
     DesignSpace(StdVector<RegionIdType>& regionIds, StdVector<ParamNode*>& design, StdVector<ParamNode*>& transfer, StdVector<ParamNode*>& result,
         ErsatzMaterial::Method method = ErsatzMaterial::NO_METHOD);

     ~DesignSpace();
    
     /** creates the corresponding DesignSpace object depending on the method */
     static DesignSpace* CreateInstance(StdVector<RegionIdType> regionIds, StdVector<ParamNode*>& design, StdVector<ParamNode*>& transfer, StdVector<ParamNode*>& result,
              ErsatzMaterial::Method method = ErsatzMaterial::NO_METHOD);

     /** PostInit as usual when not all can be stuffed into the constructor
      * @param constraints the number of constraints to initialize constraintGradients */
     void PostInit(int constraints);

     /** Set the DesignMaterial this is only used in parametric material optimization and therefore not in constructor
      * @param dm ParamNode in XML */
     void SetDesignMaterial(ParamNode* dm);

     /** Set the optimizer, required for level set give the level set values as nodal values.
      * Otherwise not required to be called */
     void SetOptimizer(BaseOptimizer* bo) { optimizer_ = bo; }

     /** This gives the ersatz material factor for an element.
      *  This fulfills the trick, that there might be more transfer function for
      *  a single element -> as in the coupling term of Piezo SIMP.
      * @param design_index use Find(Elem*, bool) to find your index -> is complicated, check it!
      * @param applic finds the real transfer function, see  GetErsatzMaterialFactor(unsigned int, const BaseForm*)
      * @return a good factor or an exception is thrown */
     double GetErsatzMaterialFactor(unsigned int design_index, Optimization::Application applic);

     /** Convenience version
      * @see  GetErsatzMaterialFactor(unsinged int, Optimization::Application) */
     double GetErsatzMaterialFactor(unsigned int design_index, const BaseForm* form)
     {
       return GetErsatzMaterialFactor(design_index, (Optimization::Application) applicationForm.Parse(form->GetName()));
     }

     /** Returns true if optimization does provide a complete tensor, not just a density */
     bool HasErsatzMaterialTensor(){
       return designMaterial != NULL;
     }

     /** Calculates the corresponding ErsatzMaterialTensor for the given element
      * @param t holds the resulting MaterialTensor
      * @param subTensor classifies the kind of Tensor needed
      * @param elem Element
      * @param design_index index of designElement
      * @param direction if !=DEFAULT calculate derivative of Tensor instead of Tensor 
      * @returns whether the given element is subject to optimization and the tensor therefore could be retrieved */
     bool GetErsatzMaterialTensor(Matrix<double>& t, SubTensorType subTensor, const Elem* elem, DesignElement::Type direction);

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
     virtual int ReadDesignFromExtern(const double* space_in);
     int ReadDesignFromExtern(const StdVector<double>& space);
     
           
     /** gives the initial guess (for the design space) 
      * @param space_out to this array of GetDesignSpaceSize() the initial guess is wrtitten to.
      * @param scaling false to return the unscaled design variables (for logging), 
      * true to return the variables as scaled for the optimizer 
      * @return the internal design_id as calculated by ReadDesignFromExtern()
      * @see SetDesignSpace() */
     virtual int WriteDesignToExtern(double* space_out, bool scaling = true) const;
     int WriteDesignToExtern(StdVector<double>& space_out, bool scaling = true) const;

     /** Similar but more general as WriteDesignToExtern() */
     virtual void WriteGradientToExtern(double* out, DesignElement::ValueSpecifier vs,
                                DesignElement::Access access, Condition* constraint = NULL, bool scaling = true) const;

     /** provide the upper and lower bounds on the design variables to the optimizer */
     virtual void WriteBoundsToExtern(double* x_l, double* x_u) const;
     
     /** Sets the value of the described design element to 0
      * @param vs what values to set. Not all make sense -> exception
      * @param design with design elements to set, DEFAULT applies dor all design types */
     virtual void Reset(DesignElement::ValueSpecifier vs, DesignElement::Type design = DesignElement::DEFAULT);
     
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
     DesignElement* Find(unsigned int elemNum, DesignElement::Type dt, bool throw_exception = true);

     /** Searches for the elment idx.
      * @param elem checks region of elem or region of vol elemens if elem is a SurfElem
      * @param throw_region_exception if no region matches returns -1 or throws exception
      * @return either the element index for GetErsatzMaterialFactor() or -1 if no region matches
      * @exception throw_region_exception suppresses only exceptions on non-matching regions! */
     int Find(const Elem* elem, bool throw_region_exception);

     /** When we have more design types this is a divisor of data.GetSieze() */
     unsigned int GetNumberOfElements() { return elements_; }

     /** The number of optimization variables */
     virtual unsigned int GetNumberOfVariables() const;
     
     /** This is our real design data, a set of DesignElements */
     StdVector<DesignElement> data;

     /** Our transfer functions */
     StdVector<TransferFunction> transfer;

     /** Here we store the design we have. Check with Find() for the element */
     StdVector<DesignElement::Type> design;

     /** Reference to DesignMaterial */
     DesignMaterial* designMaterial;

     /** Here we store result descriptions as defined in DesignElement.hh */
     StdVector<ResultDescription> resultDescriptions;

     /** Might be nonsens if our constructor is no simp or ersatz material one! */
     int GetRegionId()
     {
       assert(regions_.GetSize() == 1);
       return regions_[0].regionId;
     }

     /** We can fine define results in the optimization part: <pre>
      * <result id="optResult_2" design="density" access="plain" value="costGradient"  />
      * <result id="optResult_3" design="density" access="plain" value="objective" /></pre>
      * @return 0, 1, 2 is the index (optResult_x+1) for DesignElement::specialResult_[]. -1 if no
      * result is specified in XML */
     int GetSpecialResultIndex(DesignElement::Type design, DesignElement::ValueSpecifier value,
                               DesignElement::Detail detail = DesignElement::NONE, DesignElement::Access access = DesignElement::SMART);

     /** Dumps the design space */
     std::string ToString();

     /** Writes summary information about design variables and transfer functions into the node */
     void ToInfo(InfoNode* in);
     
     struct DesignRegion{
       RegionIdType regionId;
       unsigned int base;
       unsigned int elements;
       bool constant;
     };
     
     StdVector<DesignRegion> regions_;

     /** save parameters for scaling the design to [0..1] in the optimizer: 
      * our design = scaling * optimizer_design + translation 
      * given for every design, for every region */
     StdVector<StdVector<double> > scale_design;
     StdVector<StdVector<double> > translate_design;

     /** for SIMP type constructor we have a number of elments,
      * data size = num of design * num region elements
      * -1 is for other constructor */
     unsigned int elements_;

   protected:

     /** This number indentifies the design space. It is always incremented if ReadDesignFromExtern() reads
      * a different design */
     int design_id;

   private:
     /** Extracts a nodal values */
     double GetNodalValue(unsigned int nodeNumber, DesignElement::ValueSpecifier vs);


     template <class T>
     void ExtractResults(shared_ptr<BaseResult> result);

     /** Helper method for ExtraxtResults() */
     template <class T>
     void FillElementResults(Result<T>& result, ResultDescription& descr);

     template <class T>
     void FillNodeResults(Result<T>& result, ResultDescription& descr);

     /** Sets a ResultInfo object by a result type optResult_1/2/3 where the detailed
      * description is in the optimization/SIMP elemenet in the xml file
      * @param solutionType switch to typedef when environment.hh is gone
      * @return null if this solution type was not given in xml or a filled version.
      *         The deletion of the object is the duty of the caller (or shared pointer) */
     ResultInfo* GetResultInfo(ResultDescription& rd);

     /** finds the index of the desing element in desing.data for the element.
      * to be optimized later if needed. Searches only in the first parameters set! */
     int Find(unsigned int elemNum, bool throw_exception = true);

     /* as regionIds_ does not exist as StdVector anymore, a simple replacement for regionIds_.Find() */
     int FindRegion(RegionIdType regionId);

     unsigned int last_find_index_;

     /** This transforms FormName (Integrator) to Application -> so the enum is actually Application 
      * but here int because of circular includes :( */
     Enum<int> applicationForm;

     /** We have to know the level set method to map to nodal values */
     BaseOptimizer* optimizer_;
  };


  /** somehow this has to be in the header file, otherwise it doesn't link :( */
  template <class T>
  void DesignSpace::ExtractResults(shared_ptr<BaseResult> base_result)
  {
    // our results are up to now scalar!
    Result<T>& result = dynamic_cast<Result<T> &>(*base_result);

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

    if(ri->definedOn == ResultInfo::NODE)
      FillNodeResults(result, def);
    else
      FillElementResults(result, def);
  }

  template <class T>
  void DesignSpace::FillNodeResults(Result<T>& result, ResultDescription& descr)
  {
    Vector<T>& actSol = result.GetVector();
    actSol.Resize(result.GetEntityList()->GetSize());
    EntityIterator it = result.GetEntityList()->GetIterator();

    for(it.Begin(); !it.IsEnd(); it++ )
    {
      unsigned int node = it.GetNode();

      actSol[it.GetPos()] = GetNodalValue(node, descr.value);
    }
  }

  template <class T>
  void DesignSpace::FillElementResults(Result<T>& result, ResultDescription& descr)
  {
    Vector<T>& result_data = result.GetVector();

    // this is our entitiy result, a scalar or a vector of dim 2/3
    unsigned int dofs = result.GetResultInfo()->dofNames.GetSize();
    assert(dofs >= 1 && dofs <= 3);
    StdVector<double> result_value(dofs);

    // set the result as we need it
    result_data.Resize(elements_ * dofs);

    // search where in data we are
    int base = FindDesign(descr.design);

    // loop over elements from result. We have to do it this way as the the connection
    // of design element and result element is the element(->elemeNum) but we cannot
    // search in the result for an element.
    EntityIterator it = result.GetEntityList()->GetIterator();
    for ( it.Begin(); !it.IsEnd(); it++ )
    {
      // for elements not in the design region we set to one
      for(unsigned int i = 0; i < dofs; i++)
        result_value[i] = 1.0;

      if(FindRegion(it.GetElem()->regionId) >= 0)
      {
        // note that the index is from the first design set!
        unsigned int base_index = Find(it.GetElem()->elemNum);
        // base=0 is first!
        unsigned int data_index = (base * elements_) + base_index;
        DesignElement& de = data[data_index];
        de.GetValue(descr, result_value, dofs);

#ifdef CHECK_INDEX
        if(de.elem->elemNum != it.GetElem()->elemNum)
          EXCEPTION("mixed up indices:" << de.elem->elemNum << "!=" << it.GetElem()->elemNum
              << " base_index=" << base_index << " data_index=" << data_index << " it.Pos()=" << it.GetPos());
#endif
      }
      for(unsigned int i = 0; i < dofs; i++)
        result_data[it.GetPos() * dofs + i] = result_value[i];
    }
  }



} // end of namespace

#endif /*DESIGN_SPACE_HH_*/

