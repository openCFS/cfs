#ifndef DESIGN_SPACE_HH_
#define DESIGN_SPACE_HH_

#include <stddef.h>
#include <complex>
#include <string>
#include <utility>

#include "DataInOut/ParamHandling/ParamNode.hh"
// we need it for the template implementation
#include "General/Enum.hh"
#include "General/Environment.hh"
#include "Forms/BiLinForms/BiLinearForm.hh"
#include "Optimization/Condition.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/DesignMaterial.hh"
#include "Optimization/ErsatzMaterial.hh"
#include "Optimization/Optimization.hh"
#include "Optimization/Transform.hh"
#include "Utils/StdVector.hh"
#include "MatVec/CRS_Matrix.hh"

namespace CoupledField
{
  template <class TYPE> class Matrix;
  template <class TYPE> class Result;
  struct Elem;
  struct ResultInfo;
  class CoefFunctionOpt;
  struct LocPointMapped;
  class BaseMaterial;
  class BaseOptimizer;
  class BaseResult;
  class SinglePDE;
  struct MultiMaterial;
  class Context;
  class LocalElementCache;
  struct DensityFilterMat;


  struct DensityFilterMat
    {
      Vector<double> inv_weighted_sum;
      Vector<double> filtered_vec;
      CRS_Matrix<double> filter_mat;
      DesignElement::Type designType;
      void AssembleFilterMatrix(StdVector<DesignElement>&data, int sum_neighbour, int filter_idx, unsigned int start = 0, unsigned int end = 0);
      void CacheDensityFilteredValue(const Vector<double>& design_vec);
      void ExportDensityFilterMatrix(std::string filename = "filter_matrix.mtx");

    };
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

     /** PostInit as usual when not all can be stuffed into the constructor
      * @param objectives the number of objectives
      * @param constraints the number of constraints to initialize constraintGradients */
     virtual void PostInit(int objectives, int constraints);

     /** Consist all regions of the design of a regular grid?.
      * In the derived design space we assume a non-regular grid for SHAPE_OPT and SHAPE_PARAM_MAT
      * the L-mesh of the stress constraint benchmark is meshed by gid with different positions of
      * element nodes, such that one cannot use the same element matrix, even if the grid is regular
      * therefore the attribute designSpace/enforce_unstructured
      * Regular means: All elements have same size but not necessarily that the domain is square. */
     bool IsRegular() const { return is_regular_; }

     /**
      * Is the design the whole mesh consisting of regular elements and completely filled by elements?
      * Is not true for a sparse mesh? */
     bool IsCubic() const { return is_cubic_; }

     /** Set the DesignMaterial this is only used in parametric material optimization and therefore not in constructor
      * @param dm ParamNode in XML
      * @param material optimization material object */
     void SetDesignMaterial(PtrParamNode dm, OptimizationMaterial::System material);

     /** returns the type of the DesignMaterial of the Ersatzmaterial **/
     DesignMaterial::Type GetDesignMaterialType()
     {
       assert(Optimization::context->dm != NULL);
       return Optimization::context->dm->GetType();
     }

     /** @see ErsaszMaterial::IsParamMat() for a weaker test */
     ErsatzMaterial::Method GetMethod() const { return this->method_; };

     /** Do we do multiscale FEM, where we model not the tensor as in FEM but the local element matrix */
     bool DoMSFEM() const { return Optimization::context->dm != NULL && Optimization::context->dm->GetType() == DesignMaterial::MSFEM_C1; }

     /** Set the optimizer, required for level set give the level set values as nodal values.
      * Otherwise not required to be called */
     void SetOptimizer(BaseOptimizer* bo) { optimizer_ = bo; }

     /** Check if a region is subject to optimization (in principle, might be more complicated for piezo, ... */
     int FindRegion(RegionIdType regionId) const;

     /** Is based on LocalElementCache and applies optimization if indicated by the coef function of the form.
      * Does only work for SIMP type calls
      * @return true if we could set retMat, even if no optimization (application of design) was applied */
     template <class T>
     bool ApplyPhysicalDesignElementMatrix(BiLinearForm* form, Matrix<T>& retMat, const Elem* elem);

     /** Handles ParamMat, including MS-FEM and fallback for SIMP for ApplyPhysicalDesignElementMatrix() if disabled
      * @return true if design and retMat is set */
     template <class T>
     bool ApplyPhysicalDesign(const CoefFunctionOpt* coef, Matrix<T>& retMat, const LocPointMapped* lpm);

     /** Performs the optimization for the scalar case. This
      * @return true if design and retScal is set */
     template <class T>
     bool ApplyPhysicalDesign(const CoefFunctionOpt* coef, T& retScal, const LocPointMapped* lpm);

     /** Performs the optimization for the scalar case. This
      * @return true if design and retScal is set */
     template <class T>
     bool ApplyPhysicalDesign(const CoefFunctionOpt* coef, Vector<T>& retVec, const LocPointMapped* lpm);

     /** Checks if tensor is positive definite. */
     template <class T>
     bool TestTensorPosDef(Matrix<T>& retMat, const LocPointMapped* lpm, DesignElement::Type direction);

     /** This gives the ersatz material factor for an element.
      *  This fulfills the trick, that there might be more transfer function for
      *  a single element -> as in the coupling term of Piezo SIMP.
      * @param design_index use Find(Elem*, bool) to find your index -> is complicated, check it!
      * @param applic finds the real transfer function, see  GetErsatzMaterialFactor(unsigned int, const BaseForm*)
      * @return a good factor or an exception is thrown */
     double GetErsatzMaterialFactor(unsigned int design_index, App::Type applic, bool forBimaterial = false);

     /** @return 1.0 if there could be no proper transfer function found and save_transfer_function is true. Otherwise exception */
     double GetErsatzMaterialFactor(DesignElement* de, App::Type applic, bool forBimaterial = false, bool save_transfer_function = false);

     /** assigns the pamping matrix: pamping_ * rho * (1-rho) * M_0. (Sigmund; Morphology; 2007)
      * The mesh is assumed irregular as we have not the ErsatzMaterial::OptimizatioMaterial.
      * This method is only to be used via domain. ErsatzMaterial has its own implementation in AddMassToStiffness() */
     bool GetErsatzMaterialPamping(const Elem* elem, Matrix<double>& elemMat);

     /** do we have a slack variable and such an AuxDesign? */
     virtual bool HasSlackVariable() const { return false; }
     virtual bool HasAlphaVariable() const { return false; }

     /** returns the slack variable if present or throws an exception */
     virtual double GetSlackVariable() const { assert(false); return -1; }
     virtual double GetAlphaVariable() const { assert(false); return -1; }
/*

     /** Returns true if optimization also provides damping parameters for Rayleigh-Damping (alpha, beta) */
     //bool HasErsatzMaterialDamping() { return(designMaterial != NULL && designMaterial->DampingIsDesign()); }


     /** Calculates the corresponding ErsatzElementMatrix for the given element
      * @param t holds the resulting Element Matrix
      * @param elem Elem pointer
      * @param direction if !=DEFAULT calculate derivative of Element matrix instead of element matrix
      * @returns whether the given element is subject to optimization and the element matrix therefore could be retrieved */
     bool GetErsatzElementMatrix(Matrix<double>& t, const Elem* elem, DesignElement::Type direction);

     /** Get the ErsatzMaterialDampingParameters
      * @param alpha Damping Parameter alpha
      * @param beta Damping Parameter beta
      * @param elem the Element for which the parameters should be returned
      * @param direction if given return derivative in that direction
      * @return whether DampingParameters are optimized at all  */
     //bool GetErsatzMaterialDamping(double& alpha, double& beta, const Elem* elem, DesignElement::Type direction = DesignElement::NO_DERIVATIVE);

     /** Get the correct Damping parameter, alpha for Mass, beta for Stiffness */
     //bool GetErsatzMaterialDampingParameterForIntegrator(const Elem* elem, /* FIXME BaseForm* integrator, */double& param);

     bool HasMultiMaterial() const { return !multimaterial.IsEmpty();}

     StdVector<MultiMaterial>& GetMultiMaterials() { return multimaterial; }

     /** Gives an assembled multimaterial tensor.
      * @param elem for the element number
      * @param tf transfer function, if not given searches by itself
      * @param mc to find the right one in the piezo case
      * @param derivative if given it contains the the index of the propert material */
     bool GetMultiMaterialTensor(Matrix<double>& t, const Elem* elem, TransferFunction* tf = NULL, SubTensorType subTensor = NO_TENSOR, MaterialClass mc = MECHANIC, const DesignElement* derivative = NULL);

     /** This gets back a uniquely defined transfer function.
      * @param throw_exception if false NULL is returned when nothing is found!
      * @param use_single when there is only one transfer function, use this one and ignore design and application */
     TransferFunction* GetTransferFunction(DesignElement::Type design, App::Type application, bool throw_exception = true, bool use_single = false);

     /** Try to determine the transfer function from the design element uniquely */
     TransferFunction* GetTransferFunction(const DesignElement* de);

     /** Apply the transformations if they shall be. The transformation is identified by the excitation of the context if not explicitly given.
      * @param fallback to be returned if transformation does not apply. E.g. again the de parameter
      * @param trans optionally give the transform, such it does not come from context
      * @return null if it did not apply or transformation was out of space (e.g. when rotating) */
      DesignElement* ApplyTransformations(const DesignElement* de, DesignElement* fallback = NULL, Transform* trans = NULL) const;

     /**<p>check the optResult_1/2/3/... from the optimization/simp/result elements against
      * element results in the pde and conditionally add it as store results to the pde.</p>
      * @param pde where to checke for store results
      * @param warn shall an warning be printed if the result is not referred in the pde */
     void AppendOptimizationResults(SinglePDE* pde, bool warn);

     /** <p>Copies the relevant data from the design element to the result such that it can be written
      * to the output (e.g. gid). This is called from the CalcResults() methods of the relevant pdes
      * when the result is defined in the xml file.</p>
      * Note, that the result might come from the 'predefined' MECH_PSEUDO_DENSITY or ELEC_PSEUDO_POLARIZATION
      * or from the OPT_RESULT_1/2/3 where the actual result is defined in a xml element under
      * optimization/SIMP</p> */
     template <class T>
     void ExtractResults(shared_ptr<BaseResult> result);

     /** <p>Take the design space from an external optimizer and apply it on the problem.
      * Do not solve the problem, but call this method before CalcCostFunction() :)</p>
      * <p>The background is, that external optimizers might have a own design space array
      * (only doubles) where we have a StdVector of the complex DesignElement.</p>
      * @param space_in the design space (in variable). Size is GetDesignSpaceSize()
      * @return the design_id which is the old one if space_in did not change the design. */
     virtual int ReadDesignFromExtern(const double* space_in, bool setAndWriteCurrent = true);
     int ReadDesignFromExtern(const StdVector<double>& space, bool setAndWriteCurrent = true);
     int ReadDesignFromExtern(const Vector<double>& space, bool setAndWriteCurrent = true);

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

     /** Similar to virtual WriteDesignToExtern: but we want to have data for a specific design type
      * @param space_out to this array of GetDesignSpaceSize() the initial guess is written to.
      * @param scaling false to return the unscaled design variables (for logging),
      * true to return the variables as scaled for the optimizer
      * @param type: allow to pick particular design type if we don't want to output all available designs
      * @return the internal design_id as calculated by ReadDesignFromExtern()
      * @see SetDesignSpace() */
     int WriteDesignToExtern(double* space_out, DesignElement::Type type, bool scaling = true) const;

     int WriteDesignToExtern(StdVector<double>& space_out, bool scaling = true, DesignElement::Type type = DesignElement::ALL_DESIGNS) const;
     int WriteDesignToExtern(Vector<double>& space_out, bool scaling = true, DesignElement::Type type = DesignElement::ALL_DESIGNS) const;

     /** Similar but more general as WriteDesignToExtern().
      * @param out if it has a window writes to the window of the vector! */
     virtual void WriteGradientToExtern(StdVector<double>& out, DesignElement::ValueSpecifier vs, DesignElement::Access access, Function* f, bool scaling = true)
     {
       if(f == NULL || f->HasDenseJacobian())
         WriteDenseGradientToExtern(out, vs, access, f, scaling); // is virtual!
       else
         WriteSparseGradientToExtern(out, vs, access, f, scaling);
     }

     /** provide the upper and lower bounds on the design variables to the optimizer */
     virtual void WriteBoundsToExtern(StdVector<double>& x_l, StdVector<double>& x_u) const;

     /** allow specification of design type */
     virtual void WriteBoundsToExtern(StdVector<double>& x_l, StdVector<double>& x_u, DesignElement::Type type) const;

     virtual void WriteBoundsToExtern(double* x_l, double* x_u) const;

     /** Sets the value of the described design element to 0
      * @param vs what values to set. Not all make sense -> exception
      * @param design with design elements to set, DEFAULT applies for all design types */
     virtual void Reset(DesignElement::ValueSpecifier vs, DesignElement::Type design = DesignElement::DEFAULT);

     /** creates a gnuplot file for the current iteration. To be triggered by gradplot for FeaturedDesign stuff and only implemented shapeMap and spaghetti.
      * Called for every iteration
      * @see gradplot_ */
     virtual void WriteGradientFile() {} ;

     /** This disables the transfer functions -> sets them to NO_TYPE. This is used
      * in SIMP to calculate the original stiffness matrices.
      * The setting from the XML file is stored -> to be undone with EnableTranferFunctions() */
     void DisableTransferFunctions();

     /** Enables the transfer functions -> sets again to the xml settings after
      * temporarily disabled via DisableTransferFunctions() */
     void EnableTransferFunctions();

     /** throws also in release mode an exception if there is more than one design */
     void AssertOneDesignOnly();

     /** Service method to find our index in the design vector.
      * @return -1 if not throw_exception and not found
      * @see double context for ShapeMapDesign::FindDesign()! */
     virtual int FindDesign(DesignElement::Type dt, bool throw_exception = true) const;

     /** Service wrapper for FindDesign(), doesn't throw exception if dt cannot be find  */
     bool HasDesign(DesignElement::Type dt) const { return FindDesign(dt, false) != -1; }

     /** gives a design element by idx. Handles als AuxDesign */
     virtual BaseDesignElement* GetDesignElement(unsigned int idx);

     /** Service method to find a specific design element by element number and design type
      * @param mm_index multmaterial index. -1 for none*/
     DesignElement* Find(unsigned int elemNum, DesignElement::Type dt, bool throw_exception = true, bool include_pseudo_designs = false, int mm_index = -1);

     /** finds the index of the design element in design.data for the element.
      * Is very fast O(1) */
     int Find(unsigned int elemNum, bool throw_exception = true, bool include_pseudo_designs = false)
     {
       // LOG_DBG3(designSpace) << "Find e=" << elemNum << " ipd=" << include_pseudo_designs << " idx=" << elemToDesign[elemNum].first << " sec=" << elemToDesign[elemNum].second;
       int idx = elemToDesign[elemNum].first;
       // reset pseudo designs when we don't look for them explicitly
       if(idx != -1 && !include_pseudo_designs && elemToDesign[elemNum].second == false)
         idx = -1;
       if(idx == -1 && throw_exception)
         EXCEPTION("could not find element " << elemNum << " in our (pseudo) design space");
       return idx;
     }

     /** Searches for the element idx.
      * @param elem checks region of elem or region of vol elemens if elem is a SurfElem
      * @param throw_region_exception if no region matches returns -1 or throws exception
      * @return either the element index for GetErsatzMaterialFactor() or -1 if no region matches
      * @exception throw_region_exception suppresses only exceptions on non-matching regions! */
     int Find(const Elem* elem, bool throw_exception)
     {
       // no extensions for pseudo designs implemented, yet!
       if(FindRegion(elem->regionId) >= 0)
         return Find(elem->elemNum, throw_exception);
       // we might have surface element and it is pointing to a design element
       const SurfElem* se = dynamic_cast<const SurfElem*>(elem);
       // no chance, we are wrong
       if(se == NULL) {
         if(!throw_exception) return -1;
         EXCEPTION("element " << elem->ToString() << " not in design regions" );
       }
       else
         for(unsigned int i = 0; i < se->ptVolElems.size(); i++)
           if(se->ptVolElems[i] != NULL && FindRegion(se->ptVolElems[i]->regionId) >= 0)
             return Find(se->ptVolElems[i]->elemNum);

       if(!throw_exception)
         return -1;
       EXCEPTION("element " << elem->ToString() << " has no volume element in design region");
     }

     /** helper class to for form tpye depending on material. Exception if nothing known.  */
     static const string ToForm(MaterialClass mc, MaterialType mt);

     /** When we have more design types this is a divisor of data.GetSize() */
     unsigned int GetNumberOfElements() { return elements; }

     /** number of designs*/
     unsigned int GetNumberOfDesigns() const {return design.GetSize(); }

     /** The number of optimization variables over all regions. Counts every time.
      * Takes constant regions into account - hence can be smaller than the number of elements even for
      * multiple regions. */
     virtual unsigned int GetNumberOfVariables() const;

     /** this is the number of ersatz material variables. Only below GetNumberOfVariables()
      * if AuxDesign or ShapeDesign is used*/
     unsigned int GetNumberOfErsatzMaterialVariables() const { return DesignSpace::GetNumberOfVariables(); }

     /** this is the number of Aux/Shape variables */
     virtual int GetNumberOfAuxParameters() const { return 0; }

     /** this is the number of shape mapping param variables. > 0 means we do shape mapping */
     virtual int GetNumberOfFeatureMappingVariables() const { return 0; }

     /** Get Pamping value (e.g. Sigmund; Morpology; 2007)
      * Extend to regions if necessary!
      * @return 0 if not set. */
     double GetPampingValue() const { return pamping_; }

     /** Do we do non_design_vicinity for larger filter/slopes */
     bool DoNonDesignVicinity() const { return non_design_vicinity_; }

     PtrParamNode GetInfo() { return info_; }

     /** allows to set a line in the header of the density.xml via DesnsityFile::Create(). Done by SpaghettiDensity */
     virtual void AddToDensityHeader(PtrParamNode pn) {};

     /** SpaghettiDesign has a own python module (script) which can be selected for python functions as script=design */
     virtual pyObject* GetPythonModule() { return NULL; }

     /** implement python function get_opt_filter_values()
      * give python all GlobalFilter properties of filter as string/string dict
      * @param args which filter - default is 0 for the first one. Check 'total_filters' in result for range info */
     pyObject* PythonGetFilterProperties(pyObject* args) const;

     /** implement python function set_opt_filter_values()
      * From a dict, accept "beta" and "eta" and optionally "non_lin_scale" and "non_lin_offset".
      * If not both non_lin_ is given, it is computed automatically
      * @param args tuple with filter idx and string/string dict */
     void PythonSetFilterProperties(pyObject* args);

     /** the global LocalElementCache instance, not necessary enabled. Only a pointer for include reasons */
     LocalElementCache* elementCache = NULL;

     /** This is our real design data, a set of DesignElements.
      * Size is design.GetSize() * elements
      * @see pseudoDesigns_
      * @see totalElements_*/
     StdVector<DesignElement> data;

     /** This is the total set of design variables, including aux designs. Only used for own optimizers like FeasPP
      * @see the difference to totalElements_ */
     StdVector<BaseDesignElement*> full_data;

     /** Our transfer functions */
     StdVector<TransferFunction> transfer;

     /** Our tranformations */
     StdVector<Transform> transform;

     /** Here we store the designs we have. Check with Find() for the element.
      * Note, that multimaterial_density is not unique here! */
     StdVector<DesignID> design;

     /** Reference to DesignMaterial, allow for multiple design materials e.g. multisequence optimization */
     StdVector<DesignMaterial*> designMaterials;

     /** Here we store result descriptions as defined in DesignElement.hh */
     StdVector<ResultDescription> resultDescriptions;

     /** Here we store the GlobalFilter part of the element Filter data.
      * By region/design/excitation */
     StdVector<GlobalFilter> filter;

     /** What is the filter type we have? We do not support mixing density and sensitivity filter.
      * Note that a filter can have too small radius which means no effect.
      * @Return Filter::NO_FILTERING, Filter::DENSITY, Filter::SENSITIVITY */
     Filter::Type GetFilterType() const { return filter_type_; }

     void SetFilterType(Filter::Type ft) { filter_type_ = ft; }

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

     /** Currently only used by SGP::OuterToDesign() */
     void IncrementDesignId() { design_id++; }

     /** We can define results in the optimization part: <pre>
      * <result id="optResult_2" design="density" access="plain" value="costGradient"  />
      * <result id="optResult_3" design="density" access="plain" value="objective" /></pre>
      * @return 0, 1, 2 is the index (optResult_x+1) for DesignElement::specialResult_[]. -1 if no
      * result is specified in XML */
     int GetSpecialResultIndex(DesignElement::Type design, DesignElement::ValueSpecifier value,
                                  DesignElement::Detail detail = DesignElement::NONE,
                                  DesignElement::Access access = DesignElement::PLAIN,
                                  const std::string& excitation = "");

     /** results of special type generic where the specification generic is created dynamically and written somewhere to .info.xml.
      * E.g. for Python spaghetti
      * @param value shall be GENERIC_ELEM */
     int GetSpecialResultIndex(DesignElement::ValueSpecifier value, const std::string& generic);

     /** collect all generic results. GENERIC_ELEM and possibly late the nodes also */
     StdVector<const ResultDescription*> GetGenericResults() const;

     /** Dumps the design space
      * @param level 0 many data, 1 only design values */
     std::string ToString(int level = 0);

     /** Writes summary information about design variables and transfer functions into the node
      * @param em might be NULL if called from read ersatz material */
     virtual void ToInfo(ErsatzMaterial* em);

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

       DesignElement::Type design = DesignElement::NO_TYPE;
       RegionIdType regionId = -1;
       unsigned int base = 0;
       unsigned int elements = 0;
       DesignConstant constant = VARIABLE;
       double scale_design = 0.0;
       double translate_design = 0.0;

       /** points to the multimaterial array or is NULL */
       MultiMaterial* multimaterial = NULL;

       /** rho^3 * E_0 + (1-rho)^3 * E_bimat. Not concurrently with HasGroundMaterial()  */
       bool HasBiMaterial() const { return has_bimat; }

       /** E_ground + rho^3 * E_0. Not concurrently with HasBiMaterial() */
       bool HasGroundMaterial() const { return has_grndmat; }

       /** convenience function */
       bool HasScndMaterial() const { return has_bimat || has_grndmat; }

       /** the material is PDE dependent therefore we create and cache it on the fly. This makes it
        * easy to be also simple for load ersatz material. Shall be thread-save
        * @pde for some material we need the SubTensorType. If not given, the pde is obtained slightly more slowely*/
       PtrCoefFct GetScndMaterial(MaterialClass mc, MaterialType mt, SinglePDE* pde = nullptr);

       /** returns existing data, As the other GetBiMaterial() is called, data may be created on the fly. */
       PtrCoefFct GetSncdMaterial(const string& integrator);

       PtrCoefFct GetScndMaterial(BiLinearForm* form) { return GetSncdMaterial(form->GetName()); }

       std::string ToString() const;

       void ToInfo(PtrParamNode node) const;

       void SetBiMaterial(const std::string& material);
       void SetGroundMaterial(const std::string& material);

       /** Here we cache the lower end material class. Complicated because of piezo and stiffness, density */
       std::map<MaterialClass, std::map<MaterialType, PtrCoefFct> > scnd_materials;

       /** the label for the info.xml. bimaterial or ground matrial. Cannot be concurrently */
       std::string scnd_material;
     private:
       /* bitmaterial is rho^3 * E_0 + (1-rho)^3 * E_bimat. Not concurrentl with ground material  */
       bool has_bimat = false;

       /* ground material is simplified bimaterial E_ground + rho^3 * E_0 where the gradient does not change */
       bool has_grndmat = false;
     };

     /** Get DesignRegion.  */
     DesignRegion* GetRegion(RegionIdType id, DesignElement::Type dt = DesignElement::NO_TYPE, int multimaterial_index = -1, bool throw_exception = true);

     DesignRegion* GetRegion(RegionIdType id, MultiMaterial* mm, bool throw_exception = true);

     DesignRegion* GetRegion(RegionIdType id, MaterialClass mc, MaterialType mt, bool throw_exception = true);

     /** flattens the content of regions which means ordered by designs all regions.
      * Always creates the data, so don't call it too often */
     StdVector<DesignRegion*> GetRegions(DesignElement::Type dt = DesignElement::ALL_DESIGNS);

     /** This is a vector of design and nested regions: regions[design][region].
      Design is here the unique design. */
     StdVector<StdVector<DesignRegion> > regions;

     /** return a debug string with regions information */
     std::string DumpRegions() const;

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

     /** Solution::Read() needs the pseudo design elements. */
     StdVector<StdVector<DesignElement> >& GetPseudoDesignRegions() { return pseudoDesigns_; }

     /** Gives all used elements. This are DesignElement/FEM elements. For all designs.
      * The read design and pseudo design elements if there are some */
     StdVector<DesignElement*>& GetTotalElements() { return totalElements_; }

     /** Calculates the sum of the registered elements.
      * Necessary for the DesignElement constructor of pseudo elements (to be registered later) such that the virtual element
      * index can be set for extended pde vector element solution storage. */
     unsigned int CalcPseudoDesignElements() const;

     /** Helper function for state tracking (nodeId ist 1-based). Function: 4 * s * (1 - s), where s is the interpolated density at node with nodeId, elemId is necessary for derivative */
     double EvalInterfaceFunction(int nodeId, bool derivative = false);

     /** Helper function for state tracking; elemId is necessary for derivative */
     double CalcAverageDensityAtNode(int nodeId, bool derivative = false);

     /** for heat optimization with interface function: this = load * temperate; load is normed to 1 */
     double CalcTemperatureAtInterface(int nodeId);

     /** for SIMP type constructor we have a number of elements,
      * data size = num of design * num region elements */
     unsigned int elements;


     // A vector that holds the filter weights Matrix and the filtered Vec for all the filters
     StdVector<DensityFilterMat> density_filter;

     // If set filtering is done by matrix Vector operation by assembling a filter mat.
     // this internal bool is set to true in default mode if the number of different design is one
     bool is_matrix_filt;

     // If set filtering matrix is written to MatrixMarket format file
     // this internal bool is set to false in default mode
     bool write_matrix_filt;

  protected:


     /** handles design and region reordering */
     void WriteDenseGradientToExtern(StdVector<double>& out, DesignElement::ValueSpecifier vs,
                                DesignElement::Access access, Function* f, bool scaling = true) const;

     /** can handle the sparse slope constraint but no reordering as the dense version */
     void WriteSparseGradientToExtern(StdVector<double>& out, DesignElement::ValueSpecifier vs,
                                DesignElement::Access access, Function* f, bool scaling = true) const;

     /** Initialized the LocalElementCache. Called by PostInit() */
     void SetupLocalElementCache();

     /** very few special results don't write there content to DesignElement at once, e.g. generic results for spaghetti.
      * Then there is the option to implement this function and provide the data.
      * Called right at the beginning of ExtractResults() */
     virtual void PrepareSpecialResults() {}

     /** This number identifies the design space. It is always incremented if ReadDesignFromExtern() reads
      * a different design */
     int design_id;

     /** This is the design space info node */
     PtrParamNode info_;

     /** the timer for design setup for Constructor and PostInit(). */
     shared_ptr<Timer> setup_timer_;

     /** file for WriteGradientFile() */
     std::ofstream gradplot_;

  private:

     /** Helper for the constructor.
      * @param tf for tanh and heaviside and in the physical case! scaling and offset in tf is set!
      * @return 'lower'/'upper' or what is defined by the physical lower/physical upper  */
     double DetermineBound(PtrParamNode pn, TransferFunction* tf, const string& bound);

     /** Extracts a nodal value */
     double GetNodalValue(unsigned int nodeNumber, DesignElement::ValueSpecifier vs);

     /** Helper method for ExtraxtResults() */
     template <class T>
     void FillElementResults(Result<T>& result, ResultDescription& descr);

     template <class T>
     void FillNodeResults(Result<T>& result, ResultDescription& descr);

     /** Creates a ResultInfo object by a result type optResult_* where the detailed
      * description is in the optimization/SIMP element in the xml file.
      * stuff like definedOn, label, units, ...
      * @param solutionType switch to typedef when Environment.hh is gone
      * @return null if this solution type was not given in xml or a filled version.
      *         The deletion of the object is the duty of the caller (or shared pointer) */
     shared_ptr<ResultInfo> GenerateResultInfo(ResultDescription& rd);

     /** Setup the multimaterial vector and do sanity checks */
     void SetupMultiMaterial(ParamNodeList design_list);

     /** Here we store the multimaterial data */
     StdVector<MultiMaterial> multimaterial;

     /** We afford a large element number to design index mapping.
      * sorted by the elemNum the design index is stored.
      * For real designs, only for the first design design (index < elements) the index is stored. @see Find()!!
      * For pseudo designs the all designs for the same element have the same index within their sub pseudoDesign_ vector
      * The boolean is design type (true is standard design, false is pseudo design.
      * -1 for element numbers with no associated design */
     StdVector<std::pair<int, bool> > elemToDesign;

     /** This transforms FormName (Integrator) to App::Type -> so the enum is actually App::Type
      * but here int because of circular includes :( */
     Enum<int> applicationForm;

     /** We have to know the level set method to map to nodal values */
     BaseOptimizer* optimizer_;

     /** are all regions regular and also not designSpace/enforce_unstructured is set.
      * Note, that in the derived design space a irregular grid is assumed! */
     bool is_regular_;

     /** generally we want local element caching but for debug purpose we might want to switch it off */
     bool local_element_caching_;

     /** if the full regular space is filled with design */
     bool is_cubic_;

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

     ErsatzMaterial::Method method_;

     /** shortcut for filter type. Determined by DesignStructure::GetCommonFilterType() */
     Filter::Type filter_type_ = Filter::NO_FILTERING;
  };


  struct MultiMaterial
  {
    /** the material is PDE dependent therefore we create and cache it on the fly. This makes it
     * easy to be also simple for load ersatz material */
    BaseMaterial* GetMultiMaterial(const MaterialClass mc);

    /** for all material classes */
    void ToInfo(PtrParamNode in);

    /** material name, to be allow creation on the fly. */
    std::string name;

    /** redundant multimaterial index, starting from 0 */
    int index;

    /** this is a list of materials. One for elasticity and three for piezoelectricity */
    StdVector<std::pair<BaseMaterial*, MaterialClass> > material;
  };


} // end of namespace

#endif /*DESIGN_SPACE_HH_*/
