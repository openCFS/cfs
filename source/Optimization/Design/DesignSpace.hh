#ifndef DESIGN_SPACE_HH_
#define DESIGN_SPACE_HH_

#include <stddef.h>
#include <complex>
#include <string>
#include <utility>

#include "DataInOut/ParamHandling/ParamNode.hh"
// we need it for the template implementation
#include "Forms/baseForm.hh"
#include "General/Enum.hh"
#include "General/environment.hh"
#include "Optimization/Condition.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/DesignMaterial.hh"
#include "Optimization/ErsatzMaterial.hh"
#include "Optimization/Optimization.hh"
#include "Utils/StdVector.hh"

namespace CoupledField {
class TransferFunction;
template <class TYPE> class Matrix;
template <class TYPE> class Result;
}  // namespace CoupledField

namespace CoupledField
{
  class BaseMaterial;
  class BaseOptimizer;
  class BaseResult;
  class SinglePDE;
  struct Elem;
  struct ResultInfo;

  /** This is the container of DesingElements which also holds the transferFunctions.
   * It can be initialized by Optimization of can contain the ersatz material stuff. */
  class DesignSpace
  {
    public:
     /** Constructor for SIMP type Optimization - there we lay on a region which contains also n# elements
      * @param pn we search for design, transferFunction, result and pamping  */
     DesignSpace(StdVector<RegionIdType>& regions, PtrParamNode pn, ErsatzMaterial::Method method = ErsatzMaterial::NO_METHOD);

     virtual ~DesignSpace();
    
     /** Creates the corresponding DesignSpace object depending on the method */
     static DesignSpace* CreateInstance(StdVector<RegionIdType> regions, PtrParamNode pn, ErsatzMaterial::Method method = ErsatzMaterial::NO_METHOD);

     /** Create a clone for projection method
      * @return you have to take care of the pointer! */
     DesignSpace* Clone();

     /** PostInit as usual when not all can be stuffed into the constructor
      * @param objectives the number of objectives
      * @param constraints the number of constraints to initialize constraintGradients */
     virtual void PostInit(int objectives, int constraints);

     /** Consist all regions of the design of a regular grid?.
      * In the derived design space we assume a non-regular grid for SHAPE_OPT and SHAPE_PARAM_MAT */
     virtual bool IsRegular() const
     {
       return all_regions_regular_;
     }

     /** Set the DesignMaterial this is only used in parametric material optimization and therefore not in constructor
      * @param dm ParamNode in XML */
     void SetDesignMaterial(PtrParamNode dm, OptimizationMaterial::System material);

     /** Set the optimizer, required for level set give the level set values as nodal values.
      * Otherwise not required to be called */
     void SetOptimizer(BaseOptimizer* bo) { optimizer_ = bo; }

     /** This gives the ersatz material factor for an element.
      *  This fulfills the trick, that there might be more transfer function for
      *  a single element -> as in the coupling term of Piezo SIMP.
      * @param design_index use Find(Elem*, bool) to find your index -> is complicated, check it!
      * @param applic finds the real transfer function, see  GetErsatzMaterialFactor(unsigned int, const BaseForm*)
      * @return a good factor or an exception is thrown */
     double GetErsatzMaterialFactor(unsigned int design_index, Optimization::Application applic, bool forBimaterial = false);

     /** assigns the pamping matrix: pamping_ * rho * (1-rho) * M_0. (Sigmund; Morphology; 2007)
      * The mesh is assumed irregular as we have not the ErsatzMaterial::OptimizatioMaterial.
      * This method is only to be used via domain. ErsatzMaterial has its own implementation in AddMassToStiffness() */
     bool GetErsatzMaterialPamping(const Elem* elem, Matrix<double>& elemMat);

     /** Convenience version
      * @see  GetErsatzMaterialFactor(unsigned int, Optimization::Application) */
     double GetErsatzMaterialFactor(unsigned int design_index, const BaseForm* form, bool forBimaterial = false)
     {
       return GetErsatzMaterialFactor(design_index, (Optimization::Application) applicationForm.Parse(form->GetName()), forBimaterial);
     }

     /** do we have a slack variable and such an AuxDesign? */
     virtual bool HasSlackVariable() const { return false; }

     /** returns the slack variable if present or throws an exception */
     virtual double GetSlackVariable() const { assert(false); return -1; }

     /** Returns true if optimization does provide a complete tensor, not just a density */
     bool HasErsatzMaterialTensor() const { return designMaterial != NULL; }
     
     /** Returns true if optimization does provide a mass, currently density is not handled by this */
     bool HasErsatzMaterialMass() const { return designMaterial != NULL; }
     
     /** Returns true if optimization also provides damping parameters for Rayleigh-Damping (alpha, beta) */
     bool HasErsatzMaterialDamping() {
       return(designMaterial != NULL && designMaterial->DampingIsDesign());
     }

     bool HasPiezoCouplingTensor() const { return designMaterial != NULL; }

     bool HasDielecTensor() const { return designMaterial != NULL; }

     /** gives either elasticity tensor, dielec tensor or piezo coupling tensor
      * @param type TENSOR_TRACE, ELAST_ALL, DIELEC_TRACE, DIELEC_ALL, PIEZO_ALL. Allways the complete tensor!
      * @see GetErsatzMaterialTensor() */
     bool GetTensor(Matrix<double>& t, DesignElement::Type type, SubTensorType subTensor, const Elem* elem, DesignElement::Type direction, DesignMaterial::Notation notation = DesignMaterial::VOIGT);

     /** Calculates the corresponding ErsatzMaterialTensor for the given element
      * @param t holds the resulting MaterialTensor
      * @param subTensor classifies the kind of Tensor needed
      * @param elem Element
      * @param design_index index of designElement
      * @param direction if !=DEFAULT calculate derivative of Tensor instead of Tensor 
      * @returns whether the given element is subject to optimization and the tensor therefore could be retrieved */
     bool GetErsatzMaterialTensor(Matrix<double>& t, SubTensorType subTensor, const Elem* elem, DesignElement::Type direction, DesignMaterial::Notation notation = DesignMaterial::VOIGT);
     
     bool GetDielecTensor(Matrix<double>& t, const Elem* elem, DesignElement::Type direction);

     bool GetPiezoCouplingTensor(Matrix<double>& t, const Elem* elem, DesignElement::Type direction);

     /** Calculates the corresponding Mass for the given element, this is usually tensor trace
      * @param elem Element
      * @param direction if !=NO_DERIVATIVE calculate the derivative instead of value
      */
     double GetErsatzMaterialMass(const Elem* elem, DesignElement::Type direction);
     
     /** Get the ErsatzMaterialDampingParameters
      * @param alpha Damping Parameter alpha
      * @param beta Damping Parameter beta
      * @param elem the Element for which the parameters should be returned
      * @param direction if given return derivative in that direction
      * @return whether DampingParameters are optimized at all  */
     bool GetErsatzMaterialDamping(double& alpha, double& beta, const Elem* elem,
         DesignElement::Type direction = DesignElement::NO_DERIVATIVE);
     
     /** Get the correct Damping parameter, alpha for Mass, beta for Stiffness */
     bool GetErsatzMaterialDampingParameterForIntegrator(const Elem* elem, BaseForm* integrator, double& param);

     /** This gets back a uniquely defined transfer function.
      * @param throw_exception if false NULL is returned when nothing is found!
      * @param use_single when there is only one transfer function, use this one and ignore design and application */
     TransferFunction* GetTransferFunction(DesignElement::Type design, Optimization::Application application, bool throw_exception = true, bool use_single = false);

     /** Try to determine the transfer function from the design element uniquely */
     TransferFunction* GetTransferFunction(const DesignElement* de);

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
     
     /** Compare the design with the present. Does not change anything!
      * @return true if the designs are equal and ReadDesignFromExtern() would give the old design id */
     virtual bool CompareDesign(const double* space_in);
           
     /** gives the initial guess (for the design space) 
      * @param space_out to this array of GetDesignSpaceSize() the initial guess is written to.
      * @param scaling false to return the unscaled design variables (for logging), 
      * true to return the variables as scaled for the optimizer 
      * @return the internal design_id as calculated by ReadDesignFromExtern()
      * @see SetDesignSpace() */
     virtual int WriteDesignToExtern(double* space_out, bool scaling = true) const;
     int WriteDesignToExtern(StdVector<double>& space_out, bool scaling = true) const;
     int WriteDesignToExtern(Vector<double>& space_out, bool scaling = true) const;

     /** Similar but more general as WriteDesignToExtern().
      * @param out if it has a window writes to the window of the vector! */
     virtual void WriteGradientToExtern(StdVector<double>& out, DesignElement::ValueSpecifier vs,
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
     int FindDesign(DesignElement::Type dt, bool throw_exception = true) const;

     /** gives a design element by idx. Handles als AuxDesign */
     virtual BaseDesignElement* GetDesignElement(unsigned int idx);

     /** Service method to find a specific design element by element number and design type */
     DesignElement* Find(unsigned int elemNum, DesignElement::Type dt, bool throw_exception = true, bool include_pseudo_designs = false);

     /** Searches for the element idx.
      * @param elem checks region of elem or region of vol elemens if elem is a SurfElem
      * @param throw_region_exception if no region matches returns -1 or throws exception
      * @return either the element index for GetErsatzMaterialFactor() or -1 if no region matches
      * @exception throw_region_exception suppresses only exceptions on non-matching regions! */
     int Find(const Elem* elem, bool throw_region_exception);

     /** finds the index of the design element in design.data for the element.
      * Is very fast O(1) */
     int Find(unsigned int elemNum, bool throw_exception = true, bool include_pseudo_designs = false);

     /** When we have more design types this is a divisor of data.GetSieze() */
     unsigned int GetNumberOfElements() { return elements; }

     /** The number of optimization variables over all regions. Counts every time.
      * Takes constant regions into account - hence can be smaller than the number of elements even for
      * multiple regions. */
     virtual unsigned int GetNumberOfVariables() const;

     /** this is the number of ersatz material variables. Only below GetNumberOfVariables()
      * if AuxDesign or ShapeDesign is used*/
     unsigned int GetNumberOfErsatzMaterialVariables() const { return DesignSpace::GetNumberOfVariables(); }
     
     /** Find the element with the largest Filter neighborhood, if no filter is used or if the
      * value is not unique (what should be the case) any suitable is returned.
      * We do not cache the result, and search all, so use with care. */
     DesignElement* FindElementWithLargesFilter();

     /** Get Pamping value (e.g. Sigmund; Morpology; 2007)
      * Extend to regions if necessary!
      * @return 0 if not set. */
     double GetPampingValue() const { return pamping_; }

     /** Do we do non_design_vicinity for larger filter/slopes */
     bool DoNonDesignVicinity() const { return non_design_vicinity_; }

     /** This is our real design data, a set of DesignElements.
      * Size is design.GetSize() * elements
      * @see pseudoDesigns_
      * @see totalElements_*/
     StdVector<DesignElement> data;

     /** This is the total set of design variables, including aux designs.
      * @see the difference to totalElements_ */
     StdVector<BaseDesignElement*> full_data;

     /** Our transfer functions */
     StdVector<TransferFunction> transfer;

     /** Here we store the design we have. Check with Find() for the element */
     StdVector<DesignElement::Type> design;

     /** Reference to DesignMaterial */
     DesignMaterial* designMaterial;

     /** Here we store result descriptions as defined in DesignElement.hh */
     StdVector<ResultDescription> resultDescriptions;

     /** Might be nonsense if our constructor is no simp or ersatz material one! */
     int GetRegionId() const
     {
       // in the case of homTracking with shape optimization we do not have a 
       // design region, but GetRegionId is still called
       // so we cannot assert regions.GetSize() == 1 as before
       if(regions.GetSize() == 0) return -1;
       return regions[0][0].regionId;
     }

     /** returns the current design id */
     int GetCurrentDesignId() const { return design_id; }

     /** We can define results in the optimization part: <pre>
      * <result id="optResult_2" design="density" access="plain" value="costGradient"  />
      * <result id="optResult_3" design="density" access="plain" value="objective" /></pre>
      * @return 0, 1, 2 is the index (optResult_x+1) for DesignElement::specialResult_[]. -1 if no
      * result is specified in XML */
     int GetSpecialResultIndex(DesignElement::Type design, DesignElement::ValueSpecifier value,
                                  DesignElement::Detail detail = DesignElement::NONE,
                                  DesignElement::Access access = DesignElement::PLAIN,
                                  const std::string& excitation = "");

     /** Dumps the design space */
     std::string ToString();

     /** Writes summary information about design variables and transfer functions into the node */
     virtual void ToInfo(PtrParamNode in);
     
     typedef enum { VARIABLE, CONSTANT_PER_REGION, CONSTANT_ON_ALL_REGIONS, FIXED } DesignConstant;
     
     static Enum<DesignConstant> designConstant;
     
     /** This holds information about a region, valid for one design.
      * save parameters for scaling the design to [0..1] in the optimizer: 
      * our design = scaling * optimizer_design + translation 
      * given for every design, for every region */
     class DesignRegion
     {
     public:
       /** Default constructor as C++ has no defaults :( */
       DesignRegion();

       DesignElement::Type design;
       RegionIdType regionId;
       unsigned int base;
       unsigned int elements;
       DesignConstant constant;
       double scale_design;
       double translate_design;

       void SetBiMaterial(const std::string& material) { bimaterial_ = material; }

       bool HasBiMaterial() const;
       /** the material is PDE dependent therefore we create and cache it on the fly. This makes it
        * easy to be also simple for load ersatz material */
       BaseMaterial* GetBiMaterial(const MaterialClass mc);

       std::string ToString() const;

       void ToInfo(PtrParamNode node) const;
     private:
       std::string bimaterial_;
       StdVector<std::pair<BaseMaterial*, MaterialClass> > materials_;
     };
     
     /** trivial find */
     DesignRegion* GetRegion(RegionIdType id, bool throw_exception = true);

     /** Convenience function */
     BaseMaterial* GetBiMaterial(RegionIdType reg, Optimization::Application app, bool throw_exception = true);

     /** This now is a vector of design and region regions[design][region] */
     StdVector<StdVector<DesignRegion> > regions;

     /** it is convenient to have such a vector for some functions. Taken from regions! */
     StdVector<RegionIdType>& GetRegionIds() { return regionIds_; }

     /** Check if a region is within the (pseudo) design domain
      * @param include_pseudo also consider pseudoDesigns_. This might result in always true ?! */
     bool Contains(const RegionIdType reg, bool include_pseduo = false) const;

     /** Add a pseudo design region. Either for off-design optimization (stress) for off-design vicinity elements.
      * If the region/ type exists nothing is added.
      * @param write_out export the designs elements. either already existing or newly added
      * @return if the region was added or existed already */
     bool RegisterPseudoDesignRegion(RegionIdType region, DesignElement::Type dt, StdVector<DesignElement*>* write_out = NULL);

     /** ErsatzMaterial::Solution::Read() needs the pseudo design elements. */
     StdVector<StdVector<DesignElement> >& GetPseudoDesignRegions() { return pseudoDesigns_; }

     /** Gives all used elements. This are DesignElement/FEM elements. For all designs.
      * The read design and pseudo design elements if there are some */
     StdVector<DesignElement*>& GetTotalElements() { return totalElements_; }

     /** Calculates the sum of the registered elements.
      * Necessary for the DesignElement constructor of pseudo elements (to be registered later) such that the virtual element
      * index can be set for extended pde vector element solution storage. */
     unsigned int CalcPseudoDesignElements() const;

     /** for SIMP type constructor we have a number of elements,
      * data size = num of design * num region elements */
     unsigned int elements;

    protected:


     /** handles design and region reordering */
     void WriteDenseGradientToExtern(StdVector<double>& out, DesignElement::ValueSpecifier vs,
                                DesignElement::Access access, Condition* g = NULL, bool scaling = true) const;

     /** can handle the sparse slope constraint but no reordering as the dense version */
     void WriteSparseGradientToExtern(StdVector<double>& out, DesignElement::ValueSpecifier vs,
                                DesignElement::Access access, Condition* g = NULL, bool scaling = true) const;


     /** This number identifies the design space. It is always incremented if ReadDesignFromExtern() reads
      * a different design */
     int design_id;

     /** Sets all Material Parameters in designMaterial for given element
      * @param elem the element to be considered
      */
     bool CollectMaterialParametersForElement(const Elem* elem);
     
   private:

     /** Helper for the constructor.
      * @param tf for tanh and heaviside and in the physcical case!! scaling and offset in tf is set!!
      * @return 'lower' or what is defined by the physical lower (e.g. by physical_lower or adapt_lower)  */
     double DetermineLowerBound(PtrParamNode pn, TransferFunction* tf);

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

     
     /** We afford a large element number to design index mapping.
      * sorted by the elemNum the design index is stored.
      * For real designs, only for the first design design (index < elements) the index is stored. @see Find()!!
      * For pseudo designs the all designs for the same element have the same index within their sub pseudoDesign_ vector
      * The boolean is design type (true is standard design, false is pseudo design.
      * -1 for element numbers with no associated design */
     StdVector<std::pair<int, bool> > elemToDesign;

     /** This transforms FormName (Integrator) to Application -> so the enum is actually Application 
      * but here int because of circular includes :( */
     Enum<int> applicationForm;

     /** We have to know the level set method to map to nodal values */
     BaseOptimizer* optimizer_;

     /** are all regions regular.
      * Note, that in the derived design space a irregular grid is assumed! */
     bool all_regions_regular_;

     /** just a cache from regions */
     StdVector<RegionIdType> regionIds_;

     /** This are all non-design elements. Either for non-design optimization (stress constraints) or
      * as vicinity design elements for non_design_vicinity="true" in <ersatzMaterial>. They are identified
      * by the region-id and the design-type of the first element of each vector.
      * They are stored here, such that we can report their special results via FillElementResults().
      * The array is grouped by regions and design type! but completely unsorted. One has to check the regionId and design of an
      * arbitrary element.
      * @see RegisterPseudoDesignRegion() */
     StdVector<StdVector<DesignElement> > pseudoDesigns_;

     /** This is a sequential list of all elements. design and pseudo-design. Multiple designs, more elements :)
      * It does not include aux design as in full_data!
      * Updated eventually by RegisterPseudoDesignRegion(). */
     StdVector<DesignElement*> totalElements_;


     /** Do we want to include the non-design vicinity? For Heaviside filters and slopes */
     bool non_design_vicinity_;

     /** the pamping parameter. Extend to region on request :) */
     double pamping_;

     /** Here we save the constructing param nodes to allow to create a clone for the projection method */
     PtrParamNode pn_;

  };




} // end of namespace

#endif /*DESIGN_SPACE_HH_*/

