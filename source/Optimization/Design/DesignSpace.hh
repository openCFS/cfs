#ifndef DESIGN_SPACE_HH_
#define DESIGN_SPACE_HH_

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Optimization/TransferFunction.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/ErsatzMaterial.hh"
#include "Optimization/Design/DesignMaterial.hh"
#include "Forms/baseForm.hh"
#include "General/Enum.hh"

// we need it for the template implementation
#include "Domain/resultInfo.hh"
#include "Utils/result.hh"

namespace CoupledField
{
  template <class TYPE> class StdVector;
  class SinglePDE;
  class Elem;
  class BaseResult;
  class ResultInfo;
  class BaseOptimizer;

  /** This is the container of DesingElements which also holds the transferFunctions.
   * It can be initialized by Optimization of can contain the ersatz material stuff. */
  class DesignSpace
  {
    public:
     /** Constructor for SIMP type Optimization - there we lay on a region which contains also n# elements
      * @param result the result description list  */
     DesignSpace(StdVector<RegionIdType>& regionIds, ParamNodeList& design, ParamNodeList& transfer, ParamNodeList& result,
         ErsatzMaterial::Method method = ErsatzMaterial::NO_METHOD);

     virtual ~DesignSpace();
    
     /** creates the corresponding DesignSpace object depending on the method */
     static DesignSpace* CreateInstance(StdVector<RegionIdType> regionIds, ParamNodeList& design, ParamNodeList& transfer, ParamNodeList& result,
              ErsatzMaterial::Method method = ErsatzMaterial::NO_METHOD);

     /** PostInit as usual when not all can be stuffed into the constructor
      * @param objectives the number of objectives
      * @param constraints the number of constraints to initialize constraintGradients */
     void PostInit(int objectives, int constraints);

     /** Consist all regions of the design of a regular grid?.
      * In the derived design space we assume a non-regular grid for SHAPE_OPT and SHAPE_PARAM_MAT */
     virtual bool IsRegular() const
     {
       return all_regions_regular_;
     }

     /** Set the DesignMaterial this is only used in parametric material optimization and therefore not in constructor
      * @param dm ParamNode in XML */
     void SetDesignMaterial(PtrParamNode dm);

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
     bool HasErsatzMaterialTensor() const
     {
       return designMaterial != NULL;
     }
     
     /** Returns true if optimization does provide a mass, currently density is not handled by this */
     bool HasErsatzMaterialMass(){
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
     
     /** Calculates the corresponding Mass for the given element, this is usually tensor trace
      * @param elem Element
      * @param direction if !=NO_DERIVATIVE calculate the derivative instead of value
      */
     double GetErsatzMaterialMass(const Elem* elem, DesignElement::Type direction);

     /** This gets back a uniquely defined transfer function.
      * @param throw_exception if false NULL is returned when nothing is found! */
     TransferFunction* GetTransferFunction(DesignElement::Type design, Optimization::Application application, bool throw_exception = true);

     /** Try to determine the transfer function from the design element uniquely */
     TransferFunction* GetTransferFunction(DesignElement* de);

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

     /** Similar but more general as WriteDesignToExtern().
      * @param out if it has a window writes to the window of the vector! */
     void WriteGradientToExtern(StdVector<double>& out, DesignElement::ValueSpecifier vs,
                                DesignElement::Access access, Condition* g = NULL, bool scaling = true) const
     {
       if(g == NULL || g->HasDenseJacobian())
         WriteDenseGradientToExtern(out, vs, access, g, scaling); // is virtual!
       else
         WriteSparseGradientToExtern(out, vs, access, g, scaling);
     }

     /** provide the upper and lower bounds on the design variables to the optimizer */
     virtual void WriteBoundsToExtern(double* x_l, double* x_u) const;
     
     /** Sets the value of the described design element to 0
      * @param vs what values to set. Not all make sense -> exception
      * @param design with design elements to set, DEFAULT applies for all design types */
     virtual void Reset(DesignElement::ValueSpecifier vs, DesignElement::Type design = DesignElement::DEFAULT);
     
     /** This disables the transfer functions -> sets them to NO_TYPE. This is used
      * in SIMP to calculate the original stiffness matrices.
      * The setting from the XML file is stored -> to be undone with EnableTranferFunctions() */
     void DisableTransferFunctions();

     /** Enables the transfer functions -> sets again to the xml settings after
      * temporarily disabled via DisableTransferFunctions() */
     void EnableTransferFunctions();

     /** throws also in release mode an exception if there is more than one design */
     void AssertOneDesignOnly();

     /** Service method to find our index in the design vector
      * @return -1 if not throw_exception and not found */
     int FindDesign(DesignElement::Type dt, bool throw_exception = true);

     /** Service method to find a specific design element by element number and design type */
     DesignElement* Find(unsigned int elemNum, DesignElement::Type dt, bool throw_exception = true);

     /** Searches for the element idx.
      * @param elem checks region of elem or region of vol elemens if elem is a SurfElem
      * @param throw_region_exception if no region matches returns -1 or throws exception
      * @return either the element index for GetErsatzMaterialFactor() or -1 if no region matches
      * @exception throw_region_exception suppresses only exceptions on non-matching regions! */
     int Find(const Elem* elem, bool throw_region_exception);

     /** finds the index of the design element in design.data for the element.
      * Is very fast O(1) */
     int Find(unsigned int elemNum, bool throw_exception = true);


     /** Finds the index of this design element within the data set.
      * Is rather fast O(1) for few design types.
      * @return -1 if nothing found and not throw_exception */
     int Find(const DesignElement* de, bool throw_exception = true);

     /** When we have more design types this is a divisor of data.GetSieze() */
     unsigned int GetNumberOfElements() { return elements; }

     /** The number of optimization variables over all regions. Counts every time.
      * Takes constant regions into account - hence can be smaller than the number of elements even for
      * multiple regions. */
     virtual unsigned int GetNumberOfVariables() const;
     
     /** the tensor exists only if specified as SIMP-option */
     void SetBiMatTensor(const Matrix<double>& t) { bimattensor_ = t; }
     const Matrix<double>& GetBiMatTensor() const { return bimattensor_; }
     
     /** This is our real design data, a set of DesignElements.
      * Size is design.GetSize() * elements */
     StdVector<DesignElement> data;

     /** Our transfer functions */
     StdVector<TransferFunction> transfer;

     /** Here we store the design we have. Check with Find() for the element */
     StdVector<DesignElement::Type> design;

     /** Reference to DesignMaterial */
     DesignMaterial* designMaterial;

     /** Here we store result descriptions as defined in DesignElement.hh */
     StdVector<ResultDescription> resultDescriptions;

     /** Might be nonsense if our constructor is no simp or ersatz material one! */
     int GetRegionId()
     {
       assert(regions.GetSize() == 1);
       return regions[0].regionId;
     }

     /** We can define results in the optimization part: <pre>
      * <result id="optResult_2" design="density" access="plain" value="costGradient"  />
      * <result id="optResult_3" design="density" access="plain" value="objective" /></pre>
      * @return 0, 1, 2 is the index (optResult_x+1) for DesignElement::specialResult_[]. -1 if no
      * result is specified in XML */
     int GetSpecialResultIndex(DesignElement::Type design, DesignElement::ValueSpecifier value,
                                  DesignElement::Detail detail = DesignElement::NONE,
                                  DesignElement::Access access = DesignElement::PLAIN);

     /** Dumps the design space */
     std::string ToString();

     /** Writes summary information about design variables and transfer functions into the node */
     void ToInfo(PtrParamNode in);
     
     struct DesignRegion{
       RegionIdType regionId;
       unsigned int base;
       unsigned int elements;
       bool constant;
     };
     
     StdVector<DesignRegion> regions;

     /** save parameters for scaling the design to [0..1] in the optimizer: 
      * our design = scaling * optimizer_design + translation 
      * given for every design, for every region */
     StdVector<StdVector<double> > scale_design;
     StdVector<StdVector<double> > translate_design;

     /** for SIMP type constructor we have a number of elements,
      * data size = num of design * num region elements */
     unsigned int elements;

    protected:


     /** handles design and region reordering */
     virtual void WriteDenseGradientToExtern(StdVector<double>& out, DesignElement::ValueSpecifier vs,
                                DesignElement::Access access, Condition* g = NULL, bool scaling = true) const;

     /** This number identifies the design space. It is always incremented if ReadDesignFromExtern() reads
      * a different design */
     int design_id;

     /** Sets all Material Parameters in designMaterial for given element
      * @param elem the element to be considered
      */
     bool CollectMaterialParametersForElement(const Elem* elem);

   private:
     /** Extracts a nodal value */
     double GetNodalValue(unsigned int nodeNumber, DesignElement::ValueSpecifier vs);

     template <class T>
     void ExtractResults(shared_ptr<BaseResult> result);

     /** Helper method for ExtraxtResults() */
     template <class T>
     void FillElementResults(Result<T>& result, ResultDescription& descr);

     template <class T>
     void FillNodeResults(Result<T>& result, ResultDescription& descr);

     /** Sets a ResultInfo object by a result type optResult_1/2/3 where the detailed
      * description is in the optimization/SIMP element in the xml file
      * @param solutionType switch to typedef when environment.hh is gone
      * @return null if this solution type was not given in xml or a filled version.
      *         The deletion of the object is the duty of the caller (or shared pointer) */
     ResultInfo* GetResultInfo(ResultDescription& rd);


     /** as regionIds_ does not exist as StdVector anymore, a simple replacement for regionIds_.Find() */
     int FindRegion(RegionIdType regionId);


     /** can handle the sparse slope constraint but no reordering as the dense version */
     void WriteSparseGradientToExtern(StdVector<double>& out, DesignElement::ValueSpecifier vs,
                                DesignElement::Access access, Condition* g = NULL, bool scaling = true) const;
     
     /** second material tensor for bimaterial optimization, i. e. we do not use material+void
      *  but material+material2 
      *  tensor is set by SIMP class */
     Matrix<double> bimattensor_;
     
     /** We afford a large element number to design index mapping.
      * sorted by the elemNum the design index stored.
      * For multiple designs only for the first (index < elements) is stored.
      * -1 for element numbers with no associated design */
     StdVector<int> elemToDesign;

     /** This transforms FormName (Integrator) to Application -> so the enum is actually Application 
      * but here int because of circular includes :( */
     Enum<int> applicationForm;

     /** We have to know the level set method to map to nodal values */
     BaseOptimizer* optimizer_;

     /** are all regions regular.
      * Note, that in the derived design space a irregular grid is assumed! */
     bool all_regions_regular_;
  };




} // end of namespace

#endif /*DESIGN_SPACE_HH_*/

