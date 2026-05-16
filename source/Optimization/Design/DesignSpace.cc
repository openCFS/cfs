#include "def_use_embedded_python.hh"

#include <assert.h>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include "Utils/mathParser/mathParser.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ParamHandling/MaterialHandler.hh"
#include "Domain/CoefFunction/CoefFunctionCompound.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "Domain/CoefFunction/CoefFunctionOpt.hh"
#include "Domain/CoefFunction/CoefFunctionConst.hh"
#include "Domain/CoefFunction/CoefXpr.hh"
#include "Domain/Domain.hh"
#include "Domain/Mesh/Grid.hh"
#include "Domain/Results/ResultInfo.hh"
#include "Domain/ElemMapping/Elem.hh"
#include "Domain/ElemMapping/EntityLists.hh"
#include "Domain/ElemMapping/ElemShapeMap.hh"
#include "Domain/ElemMapping/SurfElem.hh"
#include "Forms/BiLinForms/BiLinearForm.hh"
#include "Forms/LinForms/LinearForm.hh"
#include "Forms/BiLinForms/BDBInt.hh"
#include "Driver/Assemble.hh"
#include "Driver/BaseDriver.hh"
#include "General/Enum.hh"
#include "General/Exception.hh"
#include "MatVec/Matrix.hh"
#include "Optimization/Condition.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/LocalElementCache.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Design/ShapeDesign.hh"
#include "Optimization/Design/ShapeMapDesign.hh"
#include "Optimization/Design/FeatureMappingDesign.hh"
#ifdef USE_EMBEDDED_PYTHON // currently only the Python version
  #include "Optimization/Design/SpaghettiDesign.hh"
#endif
#include "Optimization/Design/SplineBoxDesign.hh"
#include "Optimization/Design/DensityFile.hh"
#include "Optimization/Design/MaterialTensor.hh"
#include "Optimization/Function.hh"
#include "Optimization/LevelSet.hh"
#include "Optimization/OptimizationMaterial.hh"
#include "Optimization/Optimizer/BaseOptimizer.hh"
#include "Optimization/Optimizer/ShapeOptimizer.hh"
#include "Optimization/TransferFunction.hh"
#include "Optimization/ErsatzMaterial.hh"
#include "Optimization/MagSIMP.hh"
#include "Optimization/Excitation.hh"
#include "Optimization/Context.hh"
#include "PDE/SinglePDE.hh"
#include "Utils/StdVector.hh"
#include "boost/lexical_cast.hpp"
#include <iomanip>

using namespace CoupledField;
using DE = DesignElement; // shortcut
namespace CoupledField {
template <class TYPE> class Vector;
}

using std::complex;
using std::string;
using boost::lexical_cast;

DEFINE_LOG(designSpace, "designSpace")
DEFINE_LOG(ersatz, "ersatzMaterialFactor")

DesignSpace::DesignSpace(StdVector<RegionIdType>& reg_data, PtrParamNode pn, ErsatzMaterial::Method method)
{
  LOG_DBG(designSpace) << "DesignSpace for regions=" << reg_data;

  method_ = method;
  pn_ = pn;
  info_ = domain->GetInfoRoot()->Get("optimization")->Get(ParamNode::HEADER)->Get("designSpace");

  setup_timer_ = info_->Get("setup_design/timer")->AsTimer();
  setup_timer_->Start();

  // check not only for pn->Has("designSpace") to hande the load_ersatzmatrial case
  non_design_vicinity_ = pn->Has("designSpace/non_design_vicinity") ? pn->Get("designSpace/non_design_vicinity")->As<bool>() : false;
  is_regular_ = (pn->Has("designSpace/enforce_unstructured") && pn->Get("designSpace/enforce_unstructured")->As<bool>()) ? false : domain->GetGrid()->IsRegionRegular(reg_data);
  local_element_caching_ = pn->Has("designSpace/local_element_cache") ? pn->Get("designSpace/local_element_cache")->As<bool>() : true;

  // make sure we have a context, even when we have no optimization
  if(!Optimization::manager.IsInitialized())
    Optimization::manager.Init();

  // for convenience
  assert(Optimization::manager.context.GetSize() > 0);
  // better not call SpaghettiDesign::IsMaterial() from this constructor
  if(ErsatzMaterial::IsParamMat(method))
    designMaterials.Resize(Optimization::manager.context.GetSize(), nullptr); // //initialize with null to avoid undefined behavior when destructor tries to delete uninitialized data

  regionIds_ = reg_data;
  design_id = 0;
  optimizer_ = NULL;
  applicationForm.SetName("DesignSpace::ApplicationForm");
  applicationForm.Add(App::ELEC, "linGradBDBInt");
  applicationForm.Add(App::MECH, "LinElastInt");
  // We follow for the stress, strain calculation the transfer functions of mech
  applicationForm.Add(App::MECH, "MechStressStrain", false);
  applicationForm.Add(App::MECH, "PiezoStressStrain", false);
  applicationForm.Add(App::BUCKLING, "PreStressInt", false);
  applicationForm.Add(App::HEAT, "HeatConductivity", false);
  applicationForm.Add(App::HEAT, "ConvectiveStiffInt", false);
  applicationForm.Add(App::MAG, "CurlCurlIntegrator", false);
  applicationForm.Add(App::MAG, "CurlCurlIntegrator-NL", false);
  applicationForm.Add(App::MAG, "CoilIntegrator", false);
  applicationForm.Add(App::PIEZO_COUPLING, "linPiezoCoupling");
  applicationForm.Add(App::CHARGE_DENSITY, "LinNeumannInt");
  applicationForm.Add(App::PRESSURE, "PressureLinForm");
  applicationForm.Add(App::MASS, "MassInt");
  // acoustic and heat
  applicationForm.Add(App::ACOUSTIC, "LaplaceIntegrator", false);
  applicationForm.Add(App::ACOUSTIC, "MassIntegrator", false);
  // read the elements
  elements = domain->GetGrid()->GetNumElems(reg_data);

  pamping_ = pn->Has("pamping") ? pn->Get("pamping/value")->As<double>() : 0.0;

  // store the CFS element (number) to design element mapping.
  // Used by Find() and the filter and vicinity neighbors
  elemToDesign.Resize(domain->GetGrid()->GetNumElems() + 1, std::make_pair(-1, true)); // 1 based.
  // if only shape opt is done, ignore the rest here, the DesignSpace is "empty"
  if(method == ErsatzMaterial::SHAPE_OPT){
    return;
  }

  // setup designs
  ParamNodeList pn_design = pn->GetList("design");

  // preprocess multimaterial - does not know regions yet
  SetupMultiMaterial(pn_design);

  int mm_count = 0;

  // number of different designs, where multimaterial design is special
  for(unsigned int d = 0; d < pn_design.GetSize(); d++)
  {
    DesignElement::Type dt = DesignElement::type.Parse(pn_design[d]->Get("name")->As<string>());
    // actually we want the designs unique, but we make an exception for MULTIMATERIAL_DENSITY
    // don't add slack, it is not really a design but triggers AuxDesign
    if(dt == DesignElement::SLACK || dt == DesignElement::ALPHA)
      continue;

    if(dt == DesignElement::MULTIMATERIAL)
    {
      DesignID id(dt, &multimaterial[mm_count]);
      design.Push_back(id);
      mm_count++;
    }
    else if(FindDesign(dt, false) < 0)
    {

      double rb = pn_design[d]->Has("relative_bound") ? pn_design[d]->Get("relative_bound")->As<double>(): -1.0;
      bool   eb = pn_design[d]->Has("enforce_bounds") ? pn_design[d]->Get("enforce_bounds")->As<bool>(): false;
      design.Push_back(DesignID(dt, NULL, rb, eb));
    }
    // tolerate non unique designs - e.g. for different regions
  }

  is_cubic_ = is_regular_ && design.GetSize() == Product(domain->GetGrid()->CalcRegulardGridDiscretization());

  // now read the transfer functions
  ParamNodeList trans_in = pn->GetList("transferFunction");

  if(!ErsatzMaterial::IsParamMat(method))
  {
    if(trans_in.GetSize() == 0)
      throw Exception("no transferFunctions given");
    transfer.Reserve(trans_in.GetSize());
    if(trans_in.GetSize() < design.GetSize() && ! HasMultiMaterial())
      throw Exception("less transferFunctions than design variable types is infeasible");
  }
  else
  {
    transfer.Reserve(trans_in.GetSize() + 1); // We reserve space for all given TransferFunctions plus the fallback IDENTITY transfer function
    transfer.Push_back(TransferFunction()); // add fallback IDENTITY transfer function at transfer[0] for parameters with no given TransferFunction
  }

  for(PtrParamNode in : trans_in)
    transfer.Push_back(TransferFunction(in, design.GetSize() == 1 ? design[0].design : DesignElement::NO_TYPE));
    // check for mass if we have harmonic and density in PostInit() before the pde's are not ready
    // we cannot call TransferFunction::RegisterTune() as we have no Optimization finished yet

  // read the optional transformations
  if(pn->Has("transform"))
  {
    ParamNodeList tr_in = pn->Get("transform")->GetChildren();
    for(unsigned int i = 0; i < tr_in.GetSize(); i++) {
      transform.Push_back(Transform(tr_in[i], this));
      transform.Last().index = i;
      if(lexical_cast<string>(i) != transform.Last().excitation_str)
        EXCEPTION("The " << (i+1) << ".transformation has excitation '" << transform.Last().excitation_str << "' but should have '" << i << "'");
    }
  }

  if(elements == 0 || design.IsEmpty())
  { // this may happen in shape optimization
    if(method !=ErsatzMaterial::SHAPE_GRAD)
    {
      if(elements == 0) throw Exception("empty regions");
      if(design.IsEmpty()) throw Exception("no designs types given.");
    }
  }
  else // 'standard' SIMP case
  {
    // seed initialization for random design
    srand(time(NULL));

    // set our own structure with is element times design parameters
    data.Reserve(elements * design.GetSize());
    totalElements_.Reserve(elements  * design.GetSize()); // the quick access copy which also combines pseudo design elements

    const unsigned int nr = reg_data.GetSize();
    // check whether all regions have all design variables and all design variables are given in all regions
    StdVector<bool> region_design;
    region_design.Resize(design.GetSize() * nr, true);
    StdVector<string> regNames;
    domain->GetGrid()->GetRegionNames(regNames);
    // only temporary
    StdVector<Elem*> elems;
    // set size of regions
    regions.Resize(design.GetSize());
    // iterate over all designs, we have to always have the same order of designs not the order from the xml file
    for(unsigned int dti = 0; dti < design.GetSize(); dti++)
    {
      assert(regions[dti].IsEmpty());
      regions[dti].Reserve(reg_data.GetSize()); // we use push back as the order of r is unknown!!
      DesignElement::Type dt = design[dti].design;
      // dt is unique or one of many multimaterials

      // find the proper design in pn_design
      int mm_count = -1; // multimaterial only
      for(unsigned int d = 0; d < pn_design.GetSize(); d++)
      {
        PtrParamNode curr_design_pn = pn_design[d];
        if(DesignElement::type.Parse(curr_design_pn->Get("name")->As<string>()) != dt)
          continue;
        if(dt == DesignElement::MULTIMATERIAL)
        {
          mm_count++; // now zero for the first
          if(design[dti].multimaterial->index != mm_count)
            continue;
        }
        string design_reg = curr_design_pn->Get("region")->As<string>();
        string design_bim = curr_design_pn->Has("bimaterial") ? curr_design_pn->Get("bimaterial")->As<string>() : "";
        string design_grm = curr_design_pn->Has("groundmaterial") ? curr_design_pn->Get("groundmaterial")->As<string>() : "";

        if(design_bim != "" && design_grm != "")
          throw Exception("'bimaterial' and 'groundmaterial' not concurrently possible");

        for(unsigned int r = 0; r < nr; r++)
        {
          domain->GetGrid()->GetElems(elems, reg_data[r]);
          unsigned int n = elems.GetSize();
          string reg = regNames[reg_data[r]];

          bool design_all = design_reg == "all";
          if(design_all || design_reg == reg)
          {
            if(!region_design[r*design.GetSize() + dti])
              throw Exception("Design/Region combination given twice!");
            region_design[r*design.GetSize() + dti] = false;

            // this is now done for all designs per region, this is called for every design and every region here
            DesignRegion* prev = regions[dti].IsEmpty() ? NULL : &(regions[dti].Last());
            DesignRegion tmp;
            regions[dti].Push_back(tmp);
            DesignRegion& dr   = regions[dti].Last();
            dr.design = dt;
            dr.multimaterial = dt == DesignElement::MULTIMATERIAL ? &(multimaterial[mm_count]) : NULL;
            dr.regionId = reg_data[r];
            dr.base = prev == NULL ? dti * elements : prev->base + prev->elements;
            LOG_DBG2(designSpace) << "dti=" << dti << " d=" << d << " r=" << r << " el=" << elements << " size=" << regions[dti].GetSize()
                                  << " old_base=" << (prev != NULL ? (int) prev->base : -1)
                                  << " old_el=" << (prev != NULL ? (int) prev->elements : -1);
            dr.elements = n;
            dr.constant = VARIABLE;
            if(curr_design_pn->Has("constant") && curr_design_pn->Get("constant")->As<bool>())
              dr.constant = design_all ? CONSTANT_ON_ALL_REGIONS : CONSTANT_PER_REGION; // we have a constant design-value on that region
            if(curr_design_pn->Has("fixed") && curr_design_pn->Get("fixed")->As<bool>()) // for load ersatz material only we have no default value from the schema!
              dr.constant = FIXED; // fixed overwrites all other settings

            dr.scale_design = 1.0;
            dr.translate_design = 0.0;
            if(design_bim != "")
              dr.SetBiMaterial(design_bim);
            if(design_grm != "")
              dr.SetGroundMaterial(design_grm);

            // for tanh and heaviside scaling and offset is set in the physical case
            TransferFunction* tf = GetTransferFunction(dt, App::MECH, false); // assume mech - otherwise we normally don't penalize - change if you need it
            double lower = DetermineBound(curr_design_pn, tf, "lower");
            double upper = DetermineBound(curr_design_pn, tf, "upper");

            // this is strange stuff from Bastian, nobody else is using
            if(curr_design_pn->Has("scale") && curr_design_pn->Get("scale")->As<bool>()){
              dr.scale_design = (upper - lower);
              dr.translate_design = lower;
            }
            bool random = curr_design_pn->Get("initial")->As<string>() == "random";
            double initial = -1;

            MathParser* mp = domain->GetMathParser();
            unsigned int mHandle = 4711;
            string expr = curr_design_pn->Get("initial")->As<string>();
            bool initDependsOnSpace = CoefFunction::ExprDependsOnSpace(mp,expr);
            if (initDependsOnSpace) {
              mHandle = mp->GetNewHandle(true);
              mp->SetExpr(mHandle,expr);
              domain->GetGrid()->SetElementBarycenters(reg_data[r], true);
            }
            else
              initial = random ? -1.0 : curr_design_pn->Get("initial")->As<double>();

            if(!random && (initial < lower || initial > upper)) {
              info_->Get(ParamNode::HEADER)->SetWarning("Initial value for design " + DesignElement::type.ToString(dt) + " not within bounds");
              if (initDependsOnSpace)
                info_->Get(ParamNode::HEADER)->SetWarning("Set initial value for design " + DesignElement::type.ToString(dt) + " to next valid value");
            }

            for(unsigned int e = 0; e < n; e++)
            {
              DesignElement de(dt, lower, upper, elems[e], data.GetSize(), dr.multimaterial);

              if (initDependsOnSpace)
              {
                mp->SetCoordinates(mHandle, *(domain->GetCoordSystem()), de.elem->extended->barycenter.GetCoordVector());
                initial = mp->Eval(mHandle);
                if (initial < lower || initial > upper)
                  initial = std::max(std::min(initial, upper), lower);
              }

              de.SetDesign(random ? (((double) rand()/(double) RAND_MAX) * (upper - lower) + lower) : initial);

              data.Push_back(de);
              totalElements_.Push_back(&data.Last());
              // append rucksack :)
              if(method == ErsatzMaterial::SIMP_METHOD
                  || method == ErsatzMaterial::PARAM_MAT
                  || method == ErsatzMaterial::SHAPE_MAP
                  || method == ErsatzMaterial::SPAGHETTI
                  || method == ErsatzMaterial::SPAGHETTI_PARAM_MAT
                  || method == ErsatzMaterial::FEATURE_MAPPING
                  || method == ErsatzMaterial::FEATURE_MAPPING_PARAM_MAT
                  || method == ErsatzMaterial::SPLINE_BOX)
              {
                DesignElement* ptr = &(data.Last());
                assert(ptr->simp == NULL);
                ptr->simp = new SIMPElement(ptr);
              }
              // store the element mapping only for the first design
              int idx = (int) data.GetSize() - 1;
              if(idx < (int) elements)
              {
                elemToDesign[de.elem->elemNum].first = idx;
                elemToDesign[de.elem->elemNum].second = true; // real designs here!
              }
            }
            if(mHandle != 4711)
              mp->ReleaseHandle(mHandle);
          }
        }
      }
    }
    for(unsigned int rd=0; rd < design.GetSize() * nr; rd++){
      if(region_design[rd]){
        throw Exception("not all designs given for all regions");
      }
    }
    LOG_DBG(designSpace) << "data size: " << elements << "*" << pn_design.GetSize() << "=" << data.GetSize();
  } // no design elements given

  // determine regions domain bounds
  domainBounds = domain->GetGrid()->CalcRegionsBoundingBox(regionIds_);
  LOG_DBG(designSpace) << "DS: domainBounds=" << domainBounds.ToString();

  // copy to be extended by aux design
  full_data.Resize(data.GetSize());
  for(unsigned int i = 0; i < data.GetSize(); i++)
    full_data[i] = dynamic_cast<BaseDesignElement*>(&(data[i]));

  // set the result descriptions which identify the solution types
  ParamNodeList result = pn->GetList("result");
  resultDescriptions.Reserve(result.GetSize());
  for(unsigned int i = 0; i < result.GetSize(); i++)
    resultDescriptions.Push_back(ResultDescription(result[i]));

  // reserve for the worst case. non_design_vicinity and off-design optimization
  pseudoDesigns_.Reserve(domain->GetGrid()->GetNumRegions() * design.GetSize());
  if(non_design_vicinity_)
  {
    for(unsigned int d = 0; d < design.GetSize(); d++)
      for(unsigned int r = 0; r < domain->GetGrid()->regionData.GetSize(); r++)
        if(!Contains(domain->GetGrid()->regionData[r].id))
          RegisterPseudoDesignRegion(domain->GetGrid()->regionData[r].id, design[d].design);
  }

   bool is_mat_possible = false;

   // If any of the region is constant case we get the design vector of only one per constant region and lot of functions
   // work on that basis. So the normal assembling of filter matrix or the other functions
   // which give out design vector needs to be modified for making mat vec based filtering work for constant region
   StdVector<DesignRegion>& cur_des = regions.Last();
   // Assume if one design element has the attribute constant_on_all_region , all other design element should have this.
   // matrix filter is required for SGP solver, but the member optimizer_ of Optimization is not set yet -> thus, take it from xml file
   if ((cur_des.Last().constant == VARIABLE && design.GetSize() == 1) || domain->GetParamRoot()->Get("optimization/optimizer/type")->As<string>() == "sgp")
     is_mat_possible = true;
   // won't work with non-indentiy mat_filter transfer function -> is also checked in DesginStructure::SetFilter()
   if(pn->Has("filters/filter/density/type") && pn->Get("filters/filter/density/type")->As<string>() == "material")
     is_mat_possible = false;
   // Check if matrix filtering is enabled by the user, there is no default value in the schema file
   // we cannot use the matrix yet for multiple design types, yet this is easy to extend!
   is_matrix_filt = is_mat_possible;
   if(pn->Has("filters/use_mat_filt"))
   {
     is_matrix_filt = pn->Get("filters/use_mat_filt")->As<bool>();
     if(is_matrix_filt && !is_mat_possible)
       EXCEPTION("use_mat_filter is implemented only for non constant region and single design type");
   }

   // Check if write filter matrix is enabled by the user
   write_matrix_filt = false;
   if(pn->Has("filters/write_mat_filt"))
     write_matrix_filt = pn->Get("filters/write_mat_filt")->As<bool>();

   setup_timer_->Stop();
}

DesignSpace::~DesignSpace(){
  for (UInt i = 0; i < designMaterials.GetSize();i++) {
    delete designMaterials[i];
    designMaterials[i] = NULL;
  }
  delete elementCache;
  elementCache = NULL;
}

double DesignSpace::DetermineBound(PtrParamNode pn, TransferFunction* tf, const string& bound)
{
  if(pn->Get("fixed")->As<bool>())
    return pn->Get("initial")->As<double>();

  bool pl = pn->Has("physical_" + bound);

  // find the proper lower/upper value. We have a design lower/upper and a physical lower/upper. The user can decide what to set.
  if(pn->Has(bound) && pl)
    throw Exception("In 'design' the attributes '" + bound + "' and 'physical_" + bound + "' must not be given concurrently.");
  if(!pl)
  {
    if(!pn->Has(bound))
      throw Exception("In 'design' give '" + bound + "' or 'physical_" + bound + "'.");
    return pn->Get(bound)->As<double>();
  }

  // we have to find the lower/upper bound by the transfer function.
  assert(tf != NULL);

  double physical = pl ? pn->Get("physical_" + bound)->As<double>() : pn->Get(bound)->As<double>();

  switch(tf->GetType())
  {
  case TransferFunction::SIMP_TYPE:
    return std::pow(physical, 1.0/tf->GetParam());

  case TransferFunction::RAMP:
    return (physical + tf->GetParam() * physical ) / (1 + tf->GetParam() * physical);

  case TransferFunction::HASHIN_SHTRIKMAN:
    return 3.0 / (1/physical + 5);

  case TransferFunction::NO_TYPE:
  case TransferFunction::FULL:
  case TransferFunction::FIXED:
    throw Exception("Invalid transfer function type for physical " + bound + " bound in design.");

  case TransferFunction::IDENTITY:
    return physical;

  case TransferFunction::HEAVISIDE:
  case TransferFunction::TANH:
  case TransferFunction::EXPRESSION:
  {
    // disable scaling
    tf->SetScaling(1.0);
    tf->SetOffset(0.0);
    // we perform the scaling by a complex way: sc * func + of with sc=scaling and of=offset
    // original func(1) = uv, func(lower) = lv
    // sc * lv + of = physical
    // sc * uv + of = 1
    // lower/upper is set to physical!
    // example. For beta=5 and eta=0.5 tanh(>= 0) >= 0.0066. A negative design is not feasible!
    // therefore we have to set scaling and lower=physical
    double lv = tf->Transform(physical);
    double uv = tf->Transform(1.0); // we hope this is always true!!!
    tf->SetScaling((1.0 - (physical - lv/uv)/(1-lv/uv))/uv);
    tf->SetOffset((physical - lv/uv) / (1 - lv/uv));
    return physical;
  }
  }
  assert(false);
  return -1;
}


DesignSpace* DesignSpace::CreateInstance(StdVector<RegionIdType> reg_data, PtrParamNode pn, ErsatzMaterial::Method method)
{
  switch(method)
  {
  case ErsatzMaterial::SHAPE_OPT:
  case ErsatzMaterial::SHAPE_PARAM_MAT:
    return new ShapeDesign(reg_data, pn, method);
  case ErsatzMaterial::SHAPE_MAP:
    return new ShapeMapDesign(reg_data, pn, method);
  case ErsatzMaterial::FEATURE_MAPPING:
  case ErsatzMaterial::FEATURE_MAPPING_PARAM_MAT:
    return new FeatureMappingDesign(reg_data, pn, method); // set anisotropic
  case ErsatzMaterial::SPAGHETTI:
  case ErsatzMaterial::SPAGHETTI_PARAM_MAT:
#ifdef USE_EMBEDDED_PYTHON
    return new SpaghettiDesign(reg_data,pn,method);
#else
    EXCEPTION("currently spaghetti optimization requires embedded Python")
#endif
  case ErsatzMaterial::SPLINE_BOX:
    return new SplineBoxDesign(reg_data, pn, method);
  default:
    if(pn->HasByVal("design", "name", "slack") ||  pn->HasByVal("design", "name", "alpha"))
      return new AuxDesign(reg_data, pn, method); // slack variable and eventually also alpha
    else
      return new DesignSpace(reg_data, pn, method);
  }
}

void DesignSpace::PostInit(int objectives, int constraints)
{
  setup_timer_->Start();

  for(unsigned int i = 0; i < transfer.GetSize(); i++)
  {
    TransferFunction* tf = &(transfer[i]);
    tf->RegisterTune(domain->GetOptimization());
  }

  // probably not the smartest way to omit a transferFunction check - Fabian
  if(ErsatzMaterial::IsParamMat(method_))
  {
    assert(Optimization::manager.IsInitialized());
    for(unsigned int i = 0; i < Optimization::manager.context.GetSize(); i++)
    {
      Context& ctxt = Optimization::manager.context[i];
      // FIXME there might be an multi sequence issue
      if(ctxt.IsComplex() && FindDesign(DesignElement::DENSITY, false) >= 0)
      {
        TransferFunction* tf = GetTransferFunction(DesignElement::DENSITY, App::MASS, false); // silent
        if(tf == NULL && ctxt.pde->GetName() != "electrostatic")
          info_->Get(ParamNode::HEADER)->Get("transferFunctions")->SetWarning("no transfer function 'mass' given for harmonic model");
      }
    }
  }

  // Log to info if we use matrix filtering
  info_->Get("filters")->Get("use_mat_filt")->SetValue(is_matrix_filt);

  LOG_DBG(designSpace) << "# objectives = " << objectives << ", # constraints = " << constraints;
  DesignElement::SetDesignSpace(this);
  for(unsigned int i = 0, n = totalElements_.GetSize(); i < n; i++)
    totalElements_[i]->PostInit(objectives, constraints);


  // we don't do caching for pure load ersatz material
  if(local_element_caching_ && domain->GetOptimization() != NULL)
  {
    // set the whole LocalElementCache to do SIMP element based and cache also ParamMat gradients
    elementCache = new LocalElementCache(this); // not yet activated

    // this is a virtual Function.
    SetupLocalElementCache();
  }

  setup_timer_->Stop();
}

void DesignSpace::SetupLocalElementCache()
{
  // DesignRegion::scnd_materials_ is not filled yet, this is intended to be done on the fly

  // iterate over all context. LocalElementCache sees only the current context
  assert(Optimization::manager.IsInitialized());
  assert(Optimization::context->context_idx == 0);
  for(unsigned int i = 0; i < Optimization::manager.context.GetSize(); i++)
  {
    Optimization::manager.SwitchContext(i);
    Context* ctxt = Optimization::context;   // load elements
    assert((ctxt->DoBloch() && ctxt->num_bloch_wave_vectors > 0) || (!ctxt->DoBloch() && ctxt->num_bloch_wave_vectors == 0));
    assert(!ctxt->DoBloch() || (ctxt->GetEigenFrequencyDriver()->GetCurrentWaveVectorIndex() == 0)); // to switch back

    // loop over the bloch wave numbers of we have them otherwise at least once w/o bloch but standard
    for(unsigned int w = 0; w < std::max((unsigned int) 1, ctxt->num_bloch_wave_vectors); w++)
    {
      if(ctxt->DoBloch())
        ctxt->GetEigenFrequencyDriver()->SetCurrentWaveVector(w);

      // we can cache material derivatives only for FMO
      // we don't need this stuff for  SGP
      if(Optimization::context->dm != NULL && domain->GetOptimization()->GetOptimizerType() != Optimization::SGP_SOLVER)
        elementCache->InitMechMatDeriv(regionIds_); // no init org! and no piezo stuff
      else
      {
        // The SIMP case with bimat and, when implemented, multimaterial
        elementCache->InitOrg();

        // regions is a vector of design with vectors of regions.
        // Example: first density with elasticity-tensor for BDBInt and then mass for MassInt
        for(unsigned int d = 0; d < regions.GetSize(); d++) // design
          for(unsigned int r = 0; r < regions[0].GetSize(); r++) // region
            if(regions[d][r].HasScndMaterial())
              elementCache->InitShadow(&regions[d][r]);

        // TODO add Multimaterial
      } // end of non-designMaterial case
    } // end of wave vector loop
    // reset wave vector to the first one
    if(ctxt->DoBloch())
      ctxt->GetEigenFrequencyDriver()->SetCurrentWaveVector(0);
  } // end of context loop

  Optimization::manager.SwitchContext(0); // go back to first which is what we expect.
}



bool DesignSpace::Contains(const RegionIdType reg, bool include_pseduo) const
{
  for(unsigned int i = 0, n = regions[0].GetSize(); i < n; i++)
    if(regions[0][i].regionId == reg)
      return true;
  if(!include_pseduo)
    return false;
  for(unsigned int i = 0, n = pseudoDesigns_.GetSize(); i < n; i++)
    if(pseudoDesigns_[i][0].elem->regionId == reg)
      return true;
  return false;
}
bool DesignSpace::RegisterPseudoDesignRegion(RegionIdType region, DesignElement::Type dt, StdVector<DesignElement*>* write_out)
{
  bool added = false;
  // check if the stuff already exists
  StdVector<DesignElement>* ptr = NULL;
  for(unsigned int i = 0; i < pseudoDesigns_.GetSize(); i++)
    if(pseudoDesigns_[i][0].elem->regionId == region && pseudoDesigns_[i][0].GetType() == dt)
      ptr = &(pseudoDesigns_[i]);
  // add if the stuff does not exist
  if(ptr == NULL)
  {
    domain->GetGrid()->SetElementBarycenters(region, true);
    StdVector<Elem*> elems;
    domain->GetGrid()->GetElems(elems, region);
    StdVector<DesignElement> tmp;
    // the Push_back must not resize!!!!
    assert(pseudoDesigns_.GetCapacity() >= pseudoDesigns_.GetSize() + 1);
    pseudoDesigns_.Push_back(tmp);
    StdVector<DesignElement>& vec = pseudoDesigns_.Last();
    // construct pseudo design elements
    vec.Reserve(elems.GetSize());
    // see GetPseudoElementIndex()
    int pei_base = GetNumberOfElements() + CalcPseudoDesignElements();
    // also store in totalElements_
    totalElements_.Reserve(totalElements_.GetSize() + elems.GetSize());
    for(unsigned int i = 0; i < elems.GetSize(); i++)
    {
      DesignElement del(elems[i], dt, i, pei_base + i); // the index i is very critical, what happens here??
      vec.Push_back(del); // copy constructor
      DesignElement* de = &(vec.Last());
      totalElements_.Push_back(de);
      de->SetDesign(1.0); // fixed, should work, if not extend
      assert(de->simp == NULL);
      de->simp = new SIMPElement(de);
      // store properly in elemToDesign
      assert(!(elemToDesign[de->elem->elemNum].first != -1 && elemToDesign[de->elem->elemNum].second == false)); // don't overwrite real design
      assert(!(elemToDesign[de->elem->elemNum].first != -1 && elemToDesign[de->elem->elemNum].first != (int) i));      // within pseudo designs als indices for the designs are the same
      elemToDesign[de->elem->elemNum].first = i;
      elemToDesign[de->elem->elemNum].second = false; // we have a pseudo design!
    }
    ptr = &(pseudoDesigns_.Last());
    assert(ptr->GetSize() == vec.GetSize());
    added = true;
  }
  // export new or old stuff?
  if(write_out != NULL)
  {
    write_out->Resize(ptr->GetSize());
    for(unsigned int i = 0, n = ptr->GetSize(); i < n; i++)
      (*write_out)[i] = &((*ptr)[i]);
  }
  return added;
}

unsigned int DesignSpace::CalcPseudoDesignElements() const
{
  unsigned int sum = 0;
  for(unsigned int i = 0; i < pseudoDesigns_.GetSize(); i++)
    sum += pseudoDesigns_[i].GetSize();
  return sum;
}

void DesignSpace::SetDesignMaterial(PtrParamNode dm, OptimizationMaterial::System material)
{
  designMaterials[Optimization::context->context_idx] = new DesignMaterial(dm, material, design, this);
  Optimization::context->dm = designMaterials[Optimization::context->context_idx];
}

void DesignSpace::SetDesignMaterial(DesignMaterial::Type type)
{
  designMaterials[Optimization::context->context_idx] = new DesignMaterial(type, this);
  Optimization::context->dm = designMaterials[Optimization::context->context_idx];
}

void DesignSpace::AppendOptimizationResults(SinglePDE* pde, bool warn)
{
  // set the result descriptions which identify the solution types
  for(unsigned int i = 0; i < resultDescriptions.GetSize(); i++)
  {
    // this stuff is from the optimization results in the xml file
    ResultDescription& rd = resultDescriptions[i];
    // generate ResultInfo objects with the names, ... generated from the description
    shared_ptr<ResultInfo> opt_res = GenerateResultInfo(rd);
    // this also adds the result as available result
    pde->DefineFieldResult(shared_ptr<FeFunction<double> >(new FeFunction<double>(NULL)), opt_res);
    // this compares the result with storeResults in the pde and activates it.
    bool added = pde->CheckStoreResult(opt_res);
    if(warn && !added) {
      assert(opt_res->definedOn == ResultInfo::NODE || opt_res->definedOn == ResultInfo::ELEMENT);
      std::stringstream ss;
      ss << "'" << SolutionTypeEnum.ToString(rd.solutionType) 
          << "' defined as 'result' in optimization, but not referenced as " 
          << (opt_res->definedOn == ResultInfo::NODE ? "node" : "element")
          << "Result";
      info_->SetWarning(ss.str());
    }
  }
}

double DesignSpace::EvalInterfaceFunction(int nodeId, bool derivative)
{
  double dens = CalcAverageDensityAtNode(nodeId,false);
  // with shape mapping density might be slightly larger one for tanh_sum or much larger with sum
  assert(dens < 1.01);
  dens = std::min(1.0, dens); // not very smooth but otherwise we open hell :(

  if (derivative)
    return 4.0 * CalcAverageDensityAtNode(nodeId,true) * (1.0 - 2.0 * dens);
  else
    return 4.0 *  dens * (1.0 - dens);
}

double DesignSpace::CalcAverageDensityAtNode(int nodeId, bool derivative)
{
  StdVector<const Elem*> elems;
  domain->GetGrid()->GetElemsNextToNode(elems, nodeId);
  double tmp = 0;
  int found = 0;
  double lower = 0.0;
  //FIXME Assume design elements are all of the same type and application is HEAT
  TransferFunction* tf = GetTransferFunction(data[0].GetType(), App::HEAT);
  lower = tf->Transform(data[0].GetLowerBound());
  double den = 1.0 / (1.0 - lower);

  for (unsigned int index = 0; index < elems.GetSize(); index++)
  {
    // s_i = 1/N_i \sum_{e \in N_i} [(rho_e - rho_min) / (1 - rho_min)]
    int design_index = Find(elems[index],false);
    if(design_index >= 0)
    {
      DesignElement& de = data[design_index];

      tmp += (de.GetPhysicalDesign(Optimization::context) - lower) * den;
      found++;

      LOG_DBG3(designSpace) << "EIF el="  << elems[index]->elemNum << " f=" << (de.GetPhysicalDesign(Optimization::context) - lower) * den;
    }
  }

  if(found == 0)
    EXCEPTION("CADAN: Node has no neighbor elements!!")

  if (derivative) {
    return 1.0 / (double) found * den;
  }
  else
    return tmp / (double) found;
}

double DesignSpace::GetNodalValue(unsigned int nodeNumber, DesignElement::ValueSpecifier vs, ResultDescription& descr)
{
  ShapeOptimizer* shopt = dynamic_cast<ShapeOptimizer*>(optimizer_);
  // Commented out for state tracking values at nodes
  // FIXME maybe throw an Exception? This should not be called without a levelset
  if(shopt != nullptr) 
  {
    if(shopt->ptrLS_ == nullptr)
      return 0.0;
    assert(shopt->ptrLS_->GetNodePointer(nodeNumber) != nullptr);
  }

  switch(vs)
  {
  case DesignElement::LEVEL_SET_VALUE:
    return shopt->ptrLS_->GetNodePointer(nodeNumber)->value;
  case DesignElement::LEVEL_SET_STATE:
    return shopt->ptrLS_->GetNodePointer(nodeNumber)->state;
  case DesignElement::SHAPEGRAD_NODE_VALUE:
    return shopt->ptrLS_->GetNodePointer(nodeNumber)->shapegrad;
  case DesignElement::LEVEL_SET_GRAD_XP:
    return shopt->ptrLS_->GetGradientAtNode(nodeNumber, 0);
  case DesignElement::LEVEL_SET_GRAD_XN:
    return shopt->ptrLS_->GetGradientAtNode(nodeNumber, 1);
  case DesignElement::LEVEL_SET_GRAD_YP:
    return shopt->ptrLS_->GetGradientAtNode(nodeNumber, 2);
  case DesignElement::LEVEL_SET_GRAD_YN:
    return shopt->ptrLS_->GetGradientAtNode(nodeNumber, 3);
  case DesignElement::LEVEL_SET_GRAD_ZP:
    return shopt->ptrLS_->GetGradientAtNode(nodeNumber, 4);
  case DesignElement::LEVEL_SET_GRAD_ZN:
    return shopt->ptrLS_->GetGradientAtNode(nodeNumber, 5);
  case DesignElement::HEAT_NODAL_TRACK_VAL:
    return dynamic_cast<ErsatzMaterial*>(domain->GetOptimization())->CalcStateTrackingAtNode(nodeNumber);
  case DesignElement::TEMP_AT_INTERFACE:
    return dynamic_cast<ErsatzMaterial*>(domain->GetOptimization())->CalcTempAtInterface(nodeNumber);
  case DesignElement::FEATURE_DISTANCE:
  case DesignElement::FEATURE_PROJECTED:
  {
    FeatureMappingDesign* fmd = dynamic_cast<FeatureMappingDesign*>(this);
    if(fmd == nullptr)
      throw Exception("no feature mapping design but feature mapping result requested: " + DesignElement::valueSpecifier.ToString(vs));
    return fmd->GetNodalSpecialResult(nodeNumber, descr);
  }  

  default:
    EXCEPTION("case not implemented")
  }
}

shared_ptr<ResultInfo> DesignSpace::GenerateResultInfo(ResultDescription& rd)
{
  // <result id="optResult_1" design="density" access="plain" value="costGradient"/>
  shared_ptr<ResultInfo> ri(new ResultInfo);

  // I hate it!!! :(
  ri->resultType = (SolutionType) rd.solutionType;
  // no space and brackets to have no problems with info.xml and no problems with the paraview calculator
  if(rd.value == DE::GENERIC_ELEM)
    ri->resultName = rd.generic;
  else
  {
    string access = "_" + DE::access.ToString(rd.access);
    if(rd.design == DE::FEATURE || rd.design == DE::ALL_FEATURES || rd.value == DE::FEATURE_GRAD || rd.detail == DE::GRAD_DISTANCE)
    {
      assert(rd.access == DE::PLAIN);
      access = ""; // skip irrelevant data _plain
    }
    string feature_id = "";
    if(rd.design == DE::FEATURE || FeatureVariable::IsFeatureVariable(rd.design))
      feature_id = "_" + rd.generic; // in the DE::FEATURE case it is even feature_2
    ri->resultName = DE::valueSpecifier.ToString(rd.value) + "_"
                   + (rd.detail != DE::NONE ? (DE::detail.ToString(rd.detail) + "_") : "")
                   + DE::type.ToString(rd.design) 
                   + feature_id
                   + access
                   + (rd.excitation >= 0 ? ("_ex_" + std::to_string(rd.excitation)) : "");
  }
  ri->unit = "";
  ri->entryType = ResultInfo::SCALAR;
  ri->dofNames = "";
  ri->fromOptimization = true;

  // in most cases we are on elements,
  switch(rd.value)
  {
  case DE::LEVEL_SET_VALUE:
  case DE::LEVEL_SET_STATE:
  case DE::SHAPEGRAD_NODE_VALUE:
  case DE::LEVEL_SET_GRAD_XP:
  case DE::LEVEL_SET_GRAD_XN:
  case DE::LEVEL_SET_GRAD_YP:
  case DE::LEVEL_SET_GRAD_YN:
  case DE::LEVEL_SET_GRAD_ZP:
  case DE::LEVEL_SET_GRAD_ZN:
  case DE::HEAT_NODAL_TRACK_VAL:
  case DE::TEMP_AT_INTERFACE:
  case DE::FEATURE_DISTANCE:
  case DE::FEATURE_PROJECTED:
    ri->definedOn = ResultInfo::NODE;
    break;
  default:
    ri->definedOn = ResultInfo::ELEMENT;
  }
  // FIXME (really needed?) ri->fctType = shared_ptr<ConstFct>(new ConstFct() );
  // let the caller or a shared pointer delete it finally
  return ri;
}

int DesignSpace::GetSpecialResultIndex(DesignElement::Type design, DesignElement::ValueSpecifier value,
                                       DesignElement::Detail detail, DesignElement::Access access, const string& excitation)
{
  assert(design != DesignElement::NO_TYPE); // this cannot be set in xml. DEFAULT can also only be set by omitting the attribute

  for(unsigned int i = 0; i < resultDescriptions.GetSize(); i++)
  {
    const ResultDescription& rd = resultDescriptions[i];
    // if either rd.desgin from xml or the given design is DEFAULT, we do NOT compare both
    if(rd.design != DesignElement::DEFAULT && design != DesignElement::DEFAULT && rd.design != design)
      continue;

    if(rd.value != value || rd.detail != detail || rd.access != access)
      continue;

    if(rd.excitation >= 0 && lexical_cast<string>(rd.excitation) != excitation)
      continue;

    if(rd.solutionType < OPT_RESULT_1 || rd.solutionType >= OPT_RESULT_BOUND)
      throw Exception("invalid solution type");

    return rd.solutionType - OPT_RESULT_1;
  }
  return -1; // the specified triple was not specified such in xml
}

int DesignSpace::GetSpecialResultIndex(DesignElement::ValueSpecifier value, const string& generic)
{
  assert(value == DesignElement::GENERIC_ELEM);

  for(const ResultDescription& rd : resultDescriptions)
    if(rd.value == value && rd.generic == generic)
      return rd.solutionType - OPT_RESULT_1;

  return -1;
}

StdVector<const ResultDescription*> DesignSpace::GetGenericResults() const
{
  StdVector<const ResultDescription*> res;

  for(const ResultDescription& rd : resultDescriptions)
    if(rd.value == DesignElement::GENERIC_ELEM)
      res.Push_back(&rd);

  return res;
}


void DesignSpace::AssertOneDesignOnly()
{
  if(design.GetSize() != 1)
    throw Exception("A feature relies on a single design only!");
}

int DesignSpace::FindDesign(DesignElement::Type dt, bool throw_exception) const
{
  // do a fallback for NO_TYPE and DEFAULT
  if(design.GetSize() == 1 && (dt == DesignElement::NO_TYPE || dt == DesignElement::DEFAULT || dt == DesignElement::ALL_DESIGNS))
    return 0;
  // this is not a real type of design, but volume constraint can operate on it, if optimization returns a complete tensor
  if(dt == DesignElement::MECH_TRACE && Optimization::context->dm != NULL)
    return 0;
  // feature mapping designs do not apply - return 0 for special results
  if(dt == DE::FEATURE || FeatureVariable::IsFeatureVariable(dt))
    return 0;

  // search where in data we are
  int base = -1;
  for(unsigned int i = 0; i < design.GetSize(); i++)
    if(design[i].design == dt)
      base = i;
  if(base == -1 && throw_exception)
    EXCEPTION("Design " << DesignElement::type.ToString(dt) << " not within " << design.GetSize() << " actual designs.");
  return base;
}

template <class T>
bool DesignSpace::ApplyPhysicalDesignElementMatrix(BiLinearForm* form, Matrix<T>& retMat, const Elem* elem)
{
  // We can do this only for SIMP type optimization.
  // ParamMat modifies the tensor which needs to be applied to the BDB form and MS-FEM has also nothing to do
  // with local element caching. For these call ApplyPhysicalDesign()
  if(Optimization::context->dm != NULL)
    return false;

  if(elementCache == NULL)
    return false;

  // for buckling we must not use local element caching as the local element matrices depend on the current stresses
  if(Optimization::context->DoBuckling())
    return false;

  // load the element matrix to apply optimization to it. If true, retMat is set with org material local element matrix
  if(!elementCache->CachedOrgElement<T>(retMat, form, elem))
    return false;

  // now we have the element matrix set. It is not necessarily a design element but then it helps to
  // speed up assembly over the iterations
  assert(retMat.GetNumRows() >= 3);

  // we cannot check for the region here, if form is a linear form (e.g.
  // pressure) but the design variable comes from elements one dimension higher
  int idx = Find(elem->elemNum, false); // This is very fast, just a lookup in an array
  if(idx == -1)
    return true; // we have the material but cannot proceed, hence retMat contains already the material for simulation, set in CachedOrgElement()

  // just validate that we have indeed the optimization case
  BaseBDBInt* bdb = dynamic_cast<BaseBDBInt*>(form);
  assert(bdb != NULL);
  assert(bdb->GetCoef());
  shared_ptr<CoefFunctionOpt> coef = dynamic_pointer_cast<CoefFunctionOpt>(bdb->GetCoef());
  assert(coef);
  assert(coef->GetState() == CoefFunctionOpt::OPT);

  App::Type app = (App::Type) applicationForm.Parse(form->GetName());
  double factor = GetErsatzMaterialFactor(idx, app, false); // this is not the bimat case

  if(app == App::MAG)
  {
    // in the mag case we cannot simply multiply the K_0 with f(rho).
    // K_0 contains v_0*v_r as material, we need v_0*(1+v_r*f(rho)-f(rho)).
    // This is done by multiplying K_0 in the non-grad case by v_0*(1+v_r*f(rho)-f(rho))/(v_0*v_r)
    // see SIMP::SetElementK() for the derivative (only done there)

    // in the Optimization case nu_r is cached in the linear case. Here we do the general nonlinear case
    MagSIMP* ms = domain->GetOptimization() != NULL ? dynamic_cast<MagSIMP*>(domain->GetOptimization()) : NULL;
    double nu_r = ms != NULL ? ms->GetRelactivity(elem, domain->GetDim()) : MagSIMP::ExtractRelactivity(coef->orgMat.get());
    double nu_0 = MagSIMP::nu_0;
    assert(nu_r > 0 && nu_0 > 0);

    double f_rho = factor; // penalized rho
    factor = nu_0*(1+nu_r*f_rho-f_rho)/(nu_0*nu_r);
    LOG_DBG2(designSpace) << "APDEM el="  << elem->elemNum << " mag reg=" << elem->regionId << " nu_r=" << nu_r << " nu_0=" << nu_0 << " f_rho=" << f_rho << " -> " << factor;
  }

  retMat *= factor;

  // check bimat : TODO handle multimaterial
  DesignRegion* dr = GetRegion(elem->regionId);
  double bimat_factor = -1.0;
  if(dr->HasBiMaterial())
  {
    const Matrix<T>& other = elementCache->CachedShadowElement<T>(form->GetName(), elem, dr->GetScndMaterial(form));
    bimat_factor = GetErsatzMaterialFactor(idx, app, true); // this is the bimat case
    retMat.Add(bimat_factor, other); // rho^p * E_l + (1-rho)^p * E_u
  }
  else if(dr->HasGroundMaterial())
  {
    const Matrix<T>& other = elementCache->CachedShadowElement<T>(form->GetName(), elem, dr->GetScndMaterial(form));
    retMat.Add(1.0, other); // E_g + rho^p * E_0
  }

  LOG_DBG2(designSpace) << "APDEM el="  << elem->elemNum << " mt=" << MaterialTypeEnum.ToString(coef->GetMaterialType()) << " f=" << factor << " bf=" << bimat_factor;
  LOG_DBG3(designSpace) << "APDEM el="  << elem->elemNum << " -> " << retMat.ToString();
  return true;
}


/** Performs the optimization.
 * @return true if design and coefMat is set */
template <class T>
bool DesignSpace::ApplyPhysicalDesign(const CoefFunctionOpt* coef, Matrix<T>& retMat, const LocPointMapped* lpm)
{
  // we cannot check for the region here, if form is a linear form (e.g.
  // pressure) but the design variable comes from elements one dimension higher
  // For the magnetic or bimat case we might have more design elements for a FE-Elem, hence idx is arbitrary!
  int idx = Find(lpm->ptEl->elemNum, false); // This is very fast, just a lookup in an array
  if(idx == -1)
    return false;

  CoefFunctionOpt::State old_state_inner = CoefFunctionOpt::OPT;

  //LOG_DBG3(designSpace) << "APD: form name=" << coef->GetForm()->GetName() << ", dir=" << DesignElement::type.ToString(coef->GetMaterialDerivative());
  // check if we shall perform param-mat -> construct the tensor by ourselves instead of multiplying it with the mat tensor
  if(Optimization::context->dm != NULL) // easy to extend to piezo and other stuff!
  {
    if(DoMSFEM())
    {
      assert(IsRegular());
      return Optimization::context->dm->GetErsatzElementMatrixMSFEM(dynamic_cast<Matrix<double>& >(retMat),lpm->ptEl,coef->GetMaterialDerivative());
    }

    // if we do prestressing (e.g. in buckling analysis), the tensor is never constructed but instead generated from stresses
    if(coef->GetForm()->GetName() != "PreStressInt")
    {
      MaterialTensor<T> retTensor(VOIGT, &retMat, false);
      // this here is the standard case for design material!!!
      bool ret = Optimization::context->dm->GetMechTensor(retTensor, coef->subTensor, lpm->ptEl, coef->GetMaterialDerivative(), false, coef);
      // anisotropic spaghetti have bias (region material), the anisotropic tensor is given as param mat fixed values
      // we add the bias here as we have no coef in DesignMaterial::GetMechTensor()
      if(ret && Optimization::context->dm->HasBias() && coef->GetMaterialDerivative() == DesignElement::NO_DERIVATIVE)
      {
        const CoefFunctionConst<T>* org = dynamic_cast<const CoefFunctionConst<T>*>(coef->orgMat.get());
        assert(org != nullptr);
        retMat.Add(1.0, org->GetTensor());
      }
      // for anisotropic feature mapping the org material is the main material, but we may add ground material
      DesignRegion* dr = GetRegion(lpm->ptEl->regionId);
      assert(!dr->HasBiMaterial()); // not implemented
      if(ret && dr->HasGroundMaterial() && coef->GetMaterialDerivative() == DesignElement::NO_DERIVATIVE)
      {
       retMat.Add(1.0, dr->GetCachedTensor<T>(MECHANIC,coef->GetMaterialType()));
       LOG_DBG3(designSpace) << "APD(M) el="  << lpm->ptEl->elemNum << " GM=" << dr->GetCachedTensor<T>(MECHANIC,coef->GetMaterialType()).ToString();
      }

      LOG_DBG2(designSpace) << "APD: e=" << lpm->ptEl->elemNum << " dm=1 " << " d=" << DesignElement::type.ToString(coef->GetMaterialDerivative())
                            <<  " r=" << ret << " rm=" << retMat.ToString();

      return ret;
      //LOG_DBG3(designSpace) << "APD: form name=" << coef->GetForm()->GetName() << ", dir=" << DesignElement::type.ToString(coef->GetMaterialDerivative());
    }

    if(coef->GetMaterialDerivative() != DesignElement::NO_DERIVATIVE)
    {
      assert(coef->GetForm()->GetName() == "PreStressInt");
      // We have "PreStressInt" and param-mat -> we have to set the direction for the stress computation
      // stress = D' * B * u, i.e. stress is a nested coeffunction
      // In OptimizationMaterial::ComputeElementMatrix the outer coeffunction is set to material derivative
      // and we do it here for the inner one (resulting in D')
      // Only then, `coef->orgMat->GetTensor(retMat, *lpm)` will use the correct tensor D' for retMat = stress
      CoefFunctionCompound<Double>* stressTens = dynamic_cast<CoefFunctionCompound<Double>*>(coef->orgMat.get());
      assert(stressTens != NULL);

      std::map<string, PtrCoefFct>& map = stressTens->GetCoefFcts();
      assert(map.find("a") != map.end());

      CoefFunctionFlux<Complex>* stressVec = dynamic_cast<CoefFunctionFlux<Complex>*>(map["a"].get());
      assert(stressVec != NULL);

      std::map<RegionIdType, BaseBDBInt* > forms = stressVec->GetForms();
      BaseBDBInt* bdb = forms[lpm->ptEl->regionId];
      assert(bdb != NULL);
      assert(bdb->GetCoef());

      CoefFunctionOpt* innerCoef = dynamic_cast<CoefFunctionOpt*>(bdb->GetCoef().get());
      old_state_inner = innerCoef->GetState();
      innerCoef->SetToMaterialDerivative(coef->GetMaterialDerivative());
    }
  } // end of design material

  // this is legacy stuff, most times ApplyPhysicalDesignElementMatrix() shall be used
  if (coef->GetForm()->GetName() == "PreStressInt")
    assert(retMat.GetNumCols() <= (domain->GetDim() == 2 ? 4 : 9));
  else
    assert(retMat.GetNumCols() <= (domain->GetDim() == 2 ? 3 : 6));
  assert(coef->GetForm() != NULL); // needs to be set manually via CoefFunctionOpt::SetForm()
  double factor = -4711; // set below
  double bimat_factor = -1.0;

  App::Type app = (App::Type) applicationForm.Parse(coef->GetForm()->GetName());

  Optimization* opt = domain->GetOptimization();
  if(app == App::BUCKLING && opt && opt->me && opt->me->DoHomogenization())
    dynamic_cast<ErsatzMaterial*>(opt)->CalcStressesForBucklingHomogenization(dynamic_cast<Matrix<Double>&>(retMat), lpm);
  else
    // we store the original material tensor in retMat
    coef->orgMat->GetTensor(retMat, *lpm);

  if(Optimization::context->dm != NULL && coef->GetMaterialDerivative() != DesignElement::NO_DERIVATIVE)
  {
    assert(coef->GetForm()->GetName() == "PreStressInt");
    // Reset the inner CoefFunction
    CoefFunctionCompound<Double>* stressTens = dynamic_cast<CoefFunctionCompound<Double>*>(coef->orgMat.get());
    std::map<string, PtrCoefFct>& map = stressTens->GetCoefFcts();
    CoefFunctionFlux<Complex>* stressVec = dynamic_cast<CoefFunctionFlux<Complex>*>(map["a"].get());
    std::map<RegionIdType, BaseBDBInt* > forms = stressVec->GetForms();
    BaseBDBInt* bdb = forms[lpm->ptEl->regionId];
    CoefFunctionOpt* innerCoef = dynamic_cast<CoefFunctionOpt*>(bdb->GetCoef().get());
    if(old_state_inner == CoefFunctionOpt::ORG)
      innerCoef->SetToOrgMaterial();
    else
      innerCoef->SetToOptimization();
  } // end of design material derivative

  if(app == App::MAG)
  {
    // retMat = nu_0 * nu_r
    // we assume the org mat to be a dim x dim diagonal matrix
    assert(retMat.GetNumRows() == domain->GetDim());
    assert(retMat.GetNumCols() == retMat.GetNumRows());
    assert(retMat[0][1] == 0.0); // shall be a diagonal matrix
    assert(retMat.IsSymmetric());

    const double nu_0 = 1/(4*M_PI*1e-7);

    // be sure not use RHS_DENSITY
    DesignElement* de = Find(lpm->ptEl->elemNum, DesignElement::DENSITY, true);
    factor = GetErsatzMaterialFactor(de, app, false);

    for(unsigned int i = 0; i < retMat.GetNumRows(); i++)
      retMat[i][i] = (retMat[i][i] * factor) + (1-factor) * nu_0;
  }
  else if(app == App::BUCKLING)
  {
    // we already applied the ErsatzMaterialFactor in the calculation of stresses
    // @see CoefFunctionFlux::GetVector
    factor = 1.0;

    ParamNodeList sequenceSteps = domain->GetParamRoot()->GetList("sequenceStep");
    double stressFilter = -1;
    for(unsigned int i=0; i<sequenceSteps.GetSize(); ++i)
      if(sequenceSteps[i]->Get("analysis")->Has("buckling"))
        stressFilter = sequenceSteps[i]->Get("analysis")->Get("buckling")->Get("stressFilter")->As<double>();

    if(stressFilter > 0.0)
    {
      // stress filtering
      // for each element, where density < stressFilter, we set the stress to 10^-10
      DesignElement* de = Find(lpm->ptEl->elemNum, DesignElement::DENSITY, false);
      if(!de)
        de = Find(lpm->ptEl->elemNum, DesignElement::STIFF1, false);
      if(de)
      {
        double density = de->GetDesign(BaseDesignElement::PLAIN);
        if(density < stressFilter)
          retMat.InitValue(pow(10,-10));
      }
    }
  }
  else // else this the plain common stuff
  {
    factor = GetErsatzMaterialFactor(idx, app, false); // this is not the bimat case
    LOG_DBG3(designSpace) << "APD(M) el="  << lpm->ptEl->elemNum << " org=" << retMat.ToString() << " factor=" << factor;
	  retMat *= factor; // true for mech and almost all other stuff
	  assert(!retMat.ContainsNaN());
  }
  assert(factor != -4711);

  DesignRegion* dr = GetRegion(lpm->ptEl->regionId);
  if(dr->HasBiMaterial() && app != App::BUCKLING)
  {
    const Matrix<T>& tmp = dr->GetCachedTensor<T>(app == App::HEAT ? THERMIC : MECHANIC, coef->GetMaterialType());
    bimat_factor = GetErsatzMaterialFactor(idx, app, true); // this is the bimat case
    LOG_DBG3(designSpace) << "APD(M) el="  << lpm->ptEl->elemNum << " before=" << retMat.ToString() << " bimat " << tmp.ToString() << " bf=" << bimat_factor;
    retMat.Add(bimat_factor,tmp); // rho^p * E_o + (1-rho)^p * E_b
  }
  else if(dr->HasGroundMaterial())
  {
    const Matrix<T>& tmp = dr->GetCachedTensor<T>(app == App::HEAT ? THERMIC : MECHANIC, coef->GetMaterialType());
    retMat.Add(1.0, tmp);
    LOG_DBG3(designSpace) << "APD(M) el="  << lpm->ptEl->elemNum << " GM tmp=" << tmp.ToString();
  }

  LOG_DBG2(designSpace) << "APD(M) el="  << lpm->ptEl->elemNum << " mt=" << MaterialTypeEnum.ToString(coef->GetMaterialType()) << " f=" << factor << " bf=" << bimat_factor;
  LOG_DBG3(designSpace) << "APD(M) el="  << lpm->ptEl->elemNum << " -> " << retMat.ToString();
  return true;
}


/** Performs the optimization for scalar material as for the mass
 * @return true if it is a design  variable and retScal is set */
template <class T>
bool DesignSpace::ApplyPhysicalDesign(const CoefFunctionOpt* coef, T& retScal, const LocPointMapped* lpm)
{
  // we cannot check for the region here, if form is a linear form (e.g.
  // pressure) but the design variable comes from elements one dimension higher
  int idx = Find(lpm->ptEl->elemNum, false);
  if(idx == -1)
    return false;

  // check for param mat -> e.g. scalar mass
  if(Optimization::context->dm != NULL)
  {
    assert(!DoMSFEM());
    retScal = Optimization::context->dm->GetMechMass(lpm->ptEl, coef->GetMaterialDerivative());
    LOG_DBG3(designSpace) << "APD(s) el="  << lpm->ptEl->elemNum << " d=" << DesignElement::type.ToString(coef->GetMaterialDerivative()) << " -> " << retScal;
    return true; // note that we have no plausibility check in GetMechMass()
  }

  double bimat_factor = -1.0;

  assert(coef->GetForm() != nullptr);
  App::Type app = (App::Type) applicationForm.Parse(coef->GetForm()->GetName());

  // factor is the pseudo density case, in case it has the penalty parameter applied
  double factor = -4711;

  // retScal becomes the original value
  coef->orgMat->GetScalar(retScal, *lpm);
  LOG_DBG3(designSpace) << "APD(s) el="  << lpm->ptEl->elemNum << " org=" << retScal << " app=" << applicationForm.ToString(app);

  // we need it in MAG, if we load density file, then we have a scalar value
  if(app == App::MAG)
  {
    // retScal = nu_0 * nu_r
    // be sure not use RHS_DENSITY
    DesignElement* de = Find(lpm->ptEl->elemNum, DesignElement::DENSITY, true);
    factor = GetErsatzMaterialFactor(de, app, false);
    const double nu_0 = 1/(4*M_PI*1e-7);
    retScal = (retScal * factor) + (1-factor) * nu_0;
  }
  else
  {
    factor = GetErsatzMaterialFactor(idx, app, false); // Not the bimat case
    retScal *= factor;
  }
  assert(factor != -4711);

  DesignRegion* dr = GetRegion(lpm->ptEl->regionId);
  assert(!dr->HasGroundMaterial()); // implement!
  if(dr->HasBiMaterial())
  {
    T bimat;
    if(app == App::HEAT)
    {
      assert(coef->GetMaterialType() == NO_MATERIAL); // DENSITY * HEAT_CAPACITY
      T dens = dr->GetCachedScalar<T>(THERMIC, DENSITY);
      bimat  = dr->GetCachedScalar<T>(THERMIC, HEAT_CAPACITY);
      LOG_DBG3(designSpace) << "APD(s) el="  << lpm->ptEl->elemNum << " bimat heat dens=" << dens << " cond=" << bimat;
      bimat *= dens;
    }
    else if(app == App::ACOUSTIC)
    {
      assert(coef->GetMaterialType() == DENSITY || coef->GetMaterialType() == ACOU_BULK_MODULUS);
      // we either scale 1/rho (for the stiffness matrix) or 1/K (for the mass matrix)
      dr->GetScndMaterial(ACOUSTIC, coef->GetMaterialType(), coef->GetPDE())->GetScalar(bimat, *lpm);
      // TODO: check why this fails! std::bad_variant_access: bad_variant_access
      //bimat  = dr->GetCachedScalar<T>(ACOUSTIC, coef->GetMaterialType());
    }
    else
    {
      // this is the default (MECHANIC) case
      bimat = dr->GetCachedScalar<T>(MECHANIC, DENSITY);
      // dr->GetScndMaterial(MECHANIC, DENSITY, coef->GetPDE())->GetScalar(bimat, *lpm);
    }

    bimat_factor = GetErsatzMaterialFactor(idx, app, true); // this is the bimat case
    retScal += bimat_factor * bimat; // rho^p * E_l + (1-rho)^p * E_u
  }
  LOG_DBG3(designSpace) << "APD(s) el="  << lpm->ptEl->elemNum << " mt=" << MaterialTypeEnum.ToString(coef->GetMaterialType()) << " f=" << factor << " bf=" << bimat_factor << " -> " << retScal;
  return true;
}

/** Performs the optimization for scalar material as for the mass
 * @return true if it is a design  variable and retVec is set */
template <class T>
bool DesignSpace::ApplyPhysicalDesign(const CoefFunctionOpt* coef, Vector<T>& retVec, const LocPointMapped* lpm)
{
  assert(Optimization::context->pde != NULL);

  LOG_DBG2(designSpace) << "APD (vec) el="  << lpm->ptEl->elemNum << " mt=" << MaterialTypeEnum.ToString(coef->GetMaterialType());

  if (Optimization::context->pde->GetParamNode()->Has("bcsAndLoads/designDependentHeatSource"))
  {

    assert(Optimization::context->pde->GetParamNode()->Has("bcsAndLoads/designDependentHeatSource"));
    //StdVector<Elem*> elems = domain->GetGrid()->GetElemsByNode(lpm->lp.number);

    coef->orgMat->GetVector(retVec, *lpm);

    retVec[0] *=  EvalInterfaceFunction(lpm->lp.number) / (double) data.GetSize();
    // TODO: what is if we are not in the design???
  }
  else
  {
    App::Type app = (App::Type) applicationForm.Parse(coef->GetFormL()->GetName());

    if(app == App::MAG)
    {
      // including coil optimization
      coef->orgMat->GetVector(retVec, *lpm);
      assert(retVec.GetSize() != 0);
      DesignElement* de = Find(lpm->ptEl->elemNum, DesignElement::RHS_DENSITY, false);

      LOG_DBG3(designSpace) << "APD(V): elem=" << lpm->ptEl->elemNum << " de=" << (de != NULL ? de->ToString() : "NULL") << " org retVec= " << retVec.ToString();
      if(de == NULL)
        return false;
      double factor = GetErsatzMaterialFactor(de, app, false); // Not the bimat case
      retVec *= factor;
      LOG_DBG3(designSpace) << "APD(V): factor=" << factor << " -> retVec= " << retVec.ToString();
    }
    else
      assert(false);
  }

  return true;
}

template <class T>
bool DesignSpace::TestTensorPosDef(Matrix<T>& retMat, const LocPointMapped* lpm, DesignElement::Type direction) {
  Vector<Double> lp_w;
  assert(retMat.GetNumRows() == retMat.GetNumCols());
  lp_w.Resize(retMat.GetNumRows());

  //lp_w.Init();
  retMat.eigenvaluesWithLapack(lp_w);
  if (direction == BaseDesignElement::NO_DERIVATIVE) {
    for (unsigned int i = 0; i < lp_w.GetSize();i++) {
      if (lp_w[i] < EPS) {
        throw Exception("The material tensor of element '" + lexical_cast<string>(lpm->ptEl->elemNum) + "' is not positive definite! '" + "'The tensor is given by '" + retMat.ToString());
        return false;
      }
    }
  }
  return true;
}

double DesignSpace::GetErsatzMaterialFactor(unsigned int design_index, App::Type applic, bool forBimaterial)
{
  // now do the trick, that the piezo coupling factor might be a product of the
  // density transfer function and the polarization transfer function
  assert(applic != App::MAG); // we don't want to multiply DENS with RHS_DENSITY

  double result = 1.0;
  // go over all design elements we have (one for design only, with polarization
  // it is two
  for(unsigned int index = design_index; index < data.GetSize(); index += elements)
  {
    // note that this loop with loop normally once or twice (piezo)
    DesignElement* de = &data[index];
    result *= GetErsatzMaterialFactor(de, applic, forBimaterial, true); // 1 if tf not found
  }
  return result;
}

double DesignSpace::GetErsatzMaterialFactor(DesignElement* de, App::Type applic, bool forBimaterial, bool save_transfer_function)
{
  // The design of the current element
  DesignElement::Type dt = de->GetType();
  // find the transfer function for our form and application.
  // There is not necessary a transfer function -> e.g. polarization
  // is for the piezo only defined on the coupling
  TransferFunction* tf = GetTransferFunction(dt, applic, false);
  if(tf == NULL)
  {
    if(save_transfer_function)
      return 1.0;
    else
      EXCEPTION("no transfer function found for dt=" << dt << " and applic=" << applic);
  }

  LOG_DBG3(designSpace) << "GEMF: dt=" << DesignElement::type.ToString(dt) << " app=" << Optimization::application.ToString(applic) << " tf found=" << (tf != NULL);
  // when we have a transformation we want the physical value for the source design
  DesignElement* trans = ApplyTransformations(de);
  DesignElement* use = trans != NULL ? trans : de;

  double transformed = tf->Transform(use, DesignElement::SMART, forBimaterial); // handles design filtering
  LOG_DBG3(designSpace) << "GEMF: ErsatzMaterial for " << de->elem->elemNum
      << " trans to " << DesignElement::ToString(trans,true)
      << "/" << Optimization::application.ToString(applic) << " for "
      << DesignElement::type.ToString(dt) << ": "
      << TransferFunction::type.ToString(tf->GetType()) << "("
      << use->GetDesign(DesignElement::PLAIN) << ") = " << transformed
      << " ex=" << (domain->GetOptimization() != NULL ? Optimization::context->GetExcitation()->index : -1)
      << " -> " << transformed;
  return transformed;
}


bool DesignSpace::GetErsatzMaterialPamping(const Elem* elem, Matrix<double>& elemMat)
{
  // see also implementation ErsatzMaterial::AddMassToStiffness() for match!!!
  static MechMat mm = MechMat(this); // Assumes irregular mesh :(
  // pamping at all -> see Sigmund; Morphology; 2007
  assert(GetPampingValue() >= 0);
  // have design?
  DesignElement* de = Find(elem->elemNum, DesignElement::DENSITY, false);
  if(de == NULL)
    return false;
  // we use the physical design variable to match better
  TransferFunction* tf = GetTransferFunction(de->GetType(), App::MASS);
  double tv = tf->Transform(de, DesignElement::SMART); // be consistent with ErsatzMaterial::AddMassToStiffness()
  // now the original mass matrix
  const Matrix<double>& mass = dynamic_cast<const Matrix<double>&>(mm.Mass(de->elem)); // FIXME might be complex!
  LOG_DBG3(designSpace) << "GEMP e=" << elem->elemNum << " mass=" << mass.ToString();
  elemMat.Resize(mass.GetNumRows(), mass.GetNumCols());
  elemMat.Assign(mass, tv * (1.0-tv) * GetPampingValue());
  LOG_DBG3(designSpace) << "GEMP e=" << elem->elemNum << "rv=" << tv << " p=" << GetPampingValue() << " -> " << (tv * (1.0-tv) * GetPampingValue());
  LOG_DBG3(designSpace) << "GEMP e=" << elem->elemNum << " ->" << elemMat.ToString();
  return true;
}

/*
bool DesignSpace::GetErsatzMaterialDamping(double& alpha, double& beta, const Elem* elem, DesignElement::Type direction){
  if(CollectMaterialParametersForElement(elem)){
    return(designMaterial->GetMaterialDamping(alpha, beta, direction));
  }
  return(false);
}

*/
bool DesignSpace::GetErsatzElementMatrix(Matrix<double>& t, const Elem* elem, DesignElement::Type direction){
  if(Optimization::context->dm != NULL){
    Optimization::context->dm->GetErsatzElementMatrixMSFEM(t, elem, direction);
    return(true);
  }
  return(false);
}

//bool DesignSpace::GetErsatzMaterialDampingParameterForIntegrator(const Elem* elem, /* FIXME BaseForm* form*, */ double& param)
//{
//  assert(false);
//  return false;

  /* FIXME
  if(CollectMaterialParametersForElement(elem)){
    double dummy = 0.0;
    if(form->GetName() == "MassInt") return(designMaterial->GetMaterialDamping(param, dummy));
    if(form->GetName() == "LinElastInt") return(designMaterial->GetMaterialDamping(dummy, param));
  }
  return(false);
  */
//}

bool DesignSpace::GetMultiMaterialTensor(Matrix<double>& t, const Elem* elem, TransferFunction* tf, SubTensorType stt, MaterialClass mc, const DesignElement* derivative)
{
  if(multimaterial.IsEmpty())
    return false;

  t.Init(); // even if we don't know the size, otherwise we sum up

  if(tf == NULL && derivative != NULL)
    tf = GetTransferFunction(derivative);
  if(tf == NULL)
    tf = GetTransferFunction(DesignElement::MULTIMATERIAL, App::MECH);

  if(derivative != NULL)
  {
    assert(derivative->multimaterial != NULL);
    assert(elem->elemNum == derivative->elem->elemNum);

    BaseMaterial* bm = derivative->multimaterial->GetMultiMaterial(mc);
    bm->GetTensor(t, BaseMaterial::ConvertMaterialClass(mc), Global::REAL, stt); // up to now no complex stuff

    t *= tf->Derivative(derivative, DesignElement::SMART);
  }
  else
  {
    int idx = Find(elem, false);

    if(idx < 0)
      return false;

    static Matrix<double> tmp;

    for(unsigned int d = 0; d < design.GetSize(); d++)
    {
      DesignElement& de = data[elements * d + idx];
      if(de.GetType() == DesignElement::MULTIMATERIAL)
      {
        assert(de.multimaterial != NULL);
        BaseMaterial* bm = de.multimaterial->GetMultiMaterial(mc);
        bm->GetTensor(tmp, BaseMaterial::ConvertMaterialClass(mc), Global::REAL, stt); // up to now no complex stuff

        // initial tensor initialization
        if(t.GetNumRows() == 0)
        {
          t.Resize(tmp.GetNumRows(), tmp.GetNumCols());
          t.Init();
        }

        t.Add(tf->Transform(&de, DesignElement::SMART), tmp);

        LOG_DBG3(designSpace) << "GMMT e=" << elem->elemNum << " des=" << d << " pl=" << de.GetDesign(DesignElement::PLAIN) << " sm=" << de.GetDesign(DesignElement::PLAIN) << " tf=" << tf->Transform(&de, DesignElement::SMART) << " tmp=" << tmp.ToString();
      }
    }
  }

  LOG_DBG3(designSpace) << "GMMT e=" << elem->elemNum << " d=" << (derivative == NULL ? -1 : derivative->multimaterial->index) << " -> " << t.ToString();

  return true;
}


TransferFunction* DesignSpace::GetTransferFunction(const DesignElement* de)
{
  TransferFunction* res = NULL;
  for(unsigned int i = 0; i < transfer.GetSize(); i++)
  {
    if(transfer[i].GetDesign() == de->GetType() || transfer[i].GetDesign() == DesignElement::DEFAULT)
    {
      if(res != NULL)
        EXCEPTION("Cannot determine unique transfer function for " << de->ToString());
      res = &transfer[i];
    }
  }
  if(res == NULL)
    EXCEPTION("None of the " << transfer.GetSize() << " transfer functions matches " << de->ToString());
  return res;
}

TransferFunction* DesignSpace::GetTransferFunction(DesignElement::Type design, App::Type application, bool throw_exception, bool use_single)
{
  //if(HasNonDensityDesignMaterial())
  //  return &transfer[0]; // this will always point to an identity transfer function, so CalcU1KU2 in ErsatzMaterial will simply work for parametric material optimization

  if(use_single && transfer.GetSize() == 1)
    return &transfer[0];

  for(unsigned int i = 0; i < transfer.GetSize(); i++)
  {
    TransferFunction* tf = &transfer[i];
    // it would be better to call CalcTrivialVolume() not unconditionally with app = App::MECH (fails for App::LBM)
    if(transfer.GetSize() > 1 && tf->GetApplication() != application) // be only that sensitive when we have more than one transfer function
      continue;
    if(tf->GetDesign() == design)
      return tf;
    if(this->design.GetSize() == 1 && design == DesignElement::DEFAULT)
      return tf;
  }
  if(Optimization::context->dm != NULL)
    return &transfer[0]; // this will always point to an identity transfer function, so CalcU1KU2 in ErsatzMaterial will simply work for parametric material optimization with no / not all TransferFunctions given

  if(throw_exception)
   throw Exception("the desired transfer function for design '" + DesignElement::type.ToString(design)
                    + "' in application '" + Optimization::application.ToString(application) + "' is not contained");
  return NULL;
}

DesignElement* DesignSpace::ApplyTransformations(const DesignElement* de, DesignElement* fallback, Transform* trans) const
{
  DesignElement* found = NULL;
  Context* ctxt = Optimization::context;

  if(transform.IsEmpty() && trans == NULL)
    return fallback;

  if(trans != NULL)
  {
    found = trans->FindSource(de);
  }
  else
  {
    Excitation* excite = ctxt->GetExcitation();
    assert(!(excite->transform == NULL && transform.GetSize() > 1));

    if(excite->transform != NULL)
    {
      found = excite->transform->FindSource(de);
      LOG_DBG2(designSpace) << "AT: de=" << de->ToString() << " ce=" << excite->label << " a=" << excite->transform->ToString() << " -> " << DesignElement::ToString(found);
    }
    else
    {
      assert(transform.GetSize() == 1);
      found = transform[0].FindSource(de);
    }
  }

  return found == NULL ? fallback : found;
}


int DesignSpace::ReadDesignFromExtern(const double* space, bool setAndWriteCurrent)
{
  bool new_design = false;
  const unsigned int nd = design.GetSize();
  unsigned int s = 0;
  for(unsigned int des = 0; des < nd; des++)
  {
    StdVector<DesignRegion>& cur_des = regions[des];
    const unsigned int nr = regions[des].GetSize();
    for(unsigned int r = 0; r < nr; r++)
    {
      DesignRegion& cur_reg = cur_des[r];
      const double scaling = cur_reg.scale_design;
      const double translation = cur_reg.translate_design;
      const unsigned int u = cur_reg.base + cur_reg.elements;
      if(cur_reg.constant == VARIABLE)
      {
        for(unsigned int d = cur_reg.base; d < u; d++)
        {
          const double v = space[s] * scaling + translation;
          if(!new_design && data[d].GetDesign(DesignElement::PLAIN) != v)
            new_design = true;

          LOG_DBG3(designSpace) << "ReadDesignFromExtern: setting design: data[" << d << "] = " << v << " = in[" << s << "]";
          data[d].SetDesign(v);
          s++; // advance in every step
        } // for d
      }
      else if(cur_reg.constant == CONSTANT_PER_REGION || cur_reg.constant == CONSTANT_ON_ALL_REGIONS)
      { // in FIXED case, nothing is done
        const double v = space[s] * scaling + translation;
        for(unsigned int d = cur_reg.base; d < u; d++)
        {
          if(!new_design && data[d].GetDesign(DesignElement::PLAIN) != v)
            new_design = true;

          LOG_DBG3(designSpace) << "ReadDesignFromExtern: setting design (constant region): data[" << d << "] = " << v << " = in[" << s << "]";
          data[d].SetDesign(v);
        } // for d
        if(cur_reg.constant == CONSTANT_PER_REGION || (cur_reg.constant == CONSTANT_ON_ALL_REGIONS && (r == nr-1) ) )
          s++; // only advance after having set all element of this region (or even all regions) to the corresponding value
      } // if/else constant
    } // for r
  } // for des
  assert(s == DesignSpace::GetNumberOfVariables());

  if(new_design)
    design_id++;

  // filter the current design
  // this has to be done before we write to the density file
  if (pn_->Has("filters") && is_matrix_filt) {
    Vector<double> des_vec;
    des_vec.Replace(const_cast<double*>(space), DesignSpace::GetNumberOfVariables());
    for(DensityFilterMat& filt : density_filter)
      filt.CacheDensityFilteredValue(des_vec);
  }

  // for cases where computation of an iteration fails (e.g. Bloch) we have the design which causes the error.
  // the design will again be written in Optimization::CommitIteration()
  if(new_design && domain->GetOptimization())
  {
    DensityFile* df = dynamic_cast<ErsatzMaterial*>(domain->GetOptimization())->GetDensityFile();
    if(df && setAndWriteCurrent)
      df->SetAndWriteCurrent(domain->GetOptimization()->GetCurrentIteration());
  }
  // allow for mutiple filters
  Vector<double> des_vec;
  des_vec.Replace(const_cast<double*>(space), DesignSpace::GetNumberOfVariables());
  if (density_filter.GetSize() > 0 && is_matrix_filt){
     for(unsigned int i = 0; i < density_filter.GetSize(); i++)
       density_filter[i].CacheDensityFilteredValue(des_vec);
  }
  return design_id;
}

int DesignSpace::ReadDesignFromExtern(const StdVector<double>& space, bool setAndWriteCurrent)
{
  return ReadDesignFromExtern(space.GetPointer(), setAndWriteCurrent);
}

int DesignSpace::ReadDesignFromExtern(const Vector<double>& space, bool setAndWriteCurrent)
{
  return ReadDesignFromExtern(space.GetPointer(), setAndWriteCurrent);
}


bool DesignSpace::CompareDesign(const double* space)
{
  const unsigned int nd = design.GetSize();
  unsigned int s = 0;
  for(unsigned int des = 0; des < nd; des++){
    StdVector<DesignRegion>& cur_des = regions[des];
    const unsigned int nr = regions[des].GetSize();
    for(unsigned int r = 0; r < nr; r++){
      DesignRegion& cur_reg = cur_des[r];
      const double scaling = cur_reg.scale_design;
      const double translation = cur_reg.translate_design;
      const unsigned int u = cur_reg.base + cur_reg.elements;
      if(cur_reg.constant == VARIABLE) {
        for(unsigned int d = cur_reg.base; d < u; d++){
          const double v = space[s] * scaling + translation;
          if(data[d].GetDesign(DesignElement::PLAIN) != v) {
            return(false);
          }
          LOG_DBG3(designSpace) << "ReadDesignFromExtern: setting design: data[" << d << "] = " << v << " = in[" << s << "]";
          data[d].SetDesign(v);
          s++; // advance in every step
        } // for d
      }else if(cur_reg.constant == CONSTANT_PER_REGION || cur_reg.constant == CONSTANT_ON_ALL_REGIONS){ // in FIXED case, nothing is done
        const double v = space[s] * scaling + translation;
        for(unsigned int d = cur_reg.base; d < u; d++){
          if(data[d].GetDesign(DesignElement::PLAIN) != v) {
            return(false);
          }
          LOG_DBG3(designSpace) << "ReadDesignFromExtern: setting design (constant region): data[" << d << "] = " << v << " = in[" << s << "]";
          data[d].SetDesign(v);
        } // for d
        if(cur_reg.constant == CONSTANT_PER_REGION || (cur_reg.constant == CONSTANT_ON_ALL_REGIONS && (r == nr-1) ) ){
          s++; // only advance after having set all element of this region (or even all regions) to the corresponding value
        }
      } // if/else constant
    } // for r
  } // for des
  assert(s == DesignSpace::GetNumberOfVariables());
  return(true);
}

int DesignSpace::WriteDesignToExtern(double* space, bool scaling) const
{
  const unsigned int nd = design.GetSize();
  unsigned int d = 0;
  for(unsigned int des = 0; des < nd; des++){
    const StdVector<DesignRegion>& cur_des = regions[des];
    const unsigned int nr = cur_des.GetSize();
    for(unsigned int r = 0; r < nr; r++){
      const DesignRegion& cur_reg = cur_des[r];
      LOG_DBG2(designSpace) << "WDTE: dr=" << cur_reg.ToString();
      const double rscaling = scaling ? 1.0 /cur_reg.scale_design : 1.0;
      const double translation = scaling ? cur_reg.translate_design : 0.0;
      if(cur_reg.constant == VARIABLE){
        const unsigned int u = cur_reg.base + cur_reg.elements;
        for(unsigned int s = cur_reg.base; s < u; s++){
          LOG_DBG3(designSpace) << "WriteDesignToExtern: non-constant region " << r << ": out[" << d << "] = design[" << s << "]=" << data[s].GetDesign(DesignElement::PLAIN);
          space[d++] = (data[s].GetDesign(DesignElement::PLAIN) - translation) * rscaling;
        }
      }else if(cur_reg.constant == CONSTANT_PER_REGION || cur_reg.constant == CONSTANT_ON_ALL_REGIONS) { // in FIXED case nothing is done
        LOG_DBG3(designSpace) << "WriteDesignToExtern: constant " << (cur_reg.constant == CONSTANT_PER_REGION ? "region" : "design") << r << ": out[" << d << "] = design[" << cur_reg.base << "]=" << data[cur_reg.base].GetDesign(DesignElement::PLAIN);
        space[d++] = (data[cur_reg.base].GetDesign(DesignElement::PLAIN) - translation) * rscaling;
        if(cur_reg.constant == CONSTANT_ON_ALL_REGIONS) break; // the other regions are ignored
      } // if/else constant
    } // for r
  } // for des
  assert(d == DesignSpace::GetNumberOfVariables());

  return design_id;
}

int DesignSpace::WriteDesignToExtern(double* space, DesignElement::Type type, bool scaling) const
{
  unsigned int d = 0;
  int did = FindDesign(type, true);
  const StdVector<DesignRegion>& cur_des = regions[did];
  const unsigned int nr = cur_des.GetSize();
  for(unsigned int r = 0; r < nr; r++){
    const DesignRegion& cur_reg = cur_des[r];
    LOG_DBG2(designSpace) << "WDTE: dr=" << cur_reg.ToString();
    const double rscaling = scaling ? 1.0 /cur_reg.scale_design : 1.0;
    const double translation = scaling ? cur_reg.translate_design : 0.0;
    if(cur_reg.constant == VARIABLE){
      const unsigned int u = cur_reg.base + cur_reg.elements;
      for(unsigned int s = cur_reg.base; s < u; s++){
        LOG_DBG3(designSpace) << "WriteDesignToExtern: non-constant region " << r << ": out[" << d << "] = design[" << s << "]=" << data[s].GetDesign(DesignElement::PLAIN);

        space[d++] = (data[s].GetDesign(DesignElement::PLAIN) - translation) * rscaling;
      }
    }else if(cur_reg.constant == CONSTANT_PER_REGION || cur_reg.constant == CONSTANT_ON_ALL_REGIONS) { // in FIXED case nothing is done
      LOG_DBG3(designSpace) << "WriteDesignToExtern: constant " << (cur_reg.constant == CONSTANT_PER_REGION ? "region" : "design") << r << ": out[" << d << "] = design[" << cur_reg.base << "]=" << data[cur_reg.base].GetDesign(DesignElement::PLAIN);
      space[d++] = (data[cur_reg.base].GetDesign(DesignElement::PLAIN) - translation) * rscaling;
      if(cur_reg.constant == CONSTANT_ON_ALL_REGIONS) break; // the other regions are ignored
    } // if/else constant
  } // for r

  return design_id;
}

int DesignSpace::WriteDesignToExtern(StdVector<double>& space_out, bool scaling, DesignElement::Type type) const
{
  assert(GetNumberOfVariables() > 0);
  space_out.Reserve(GetNumberOfVariables());
  if (type == DesignElement::ALL_DESIGNS)
    return WriteDesignToExtern(space_out.GetPointer(), scaling);
  else
    return WriteDesignToExtern(space_out.GetPointer(), type, scaling);
}

int DesignSpace::WriteDesignToExtern(Vector<double>& space_out, bool scaling, DesignElement::Type type) const
{
  space_out.Resize(GetNumberOfVariables());
  if (type == DesignElement::ALL_DESIGNS)
    return WriteDesignToExtern(space_out.GetPointer(), scaling);
  else
    return WriteDesignToExtern(space_out.GetPointer(), type, scaling);
}

void DesignSpace::WriteBoundsToExtern(StdVector<double>& x_l, StdVector<double>& x_u) const
{
  x_l.Resize(GetNumberOfVariables());
  x_u.Resize(GetNumberOfVariables());
  // TODO remove the ugly pointer variant!
  WriteBoundsToExtern(x_l.GetPointer(), x_u.GetPointer());
}

void DesignSpace::WriteBoundsToExtern(StdVector<double>& x_l, StdVector<double>& x_u, DesignElement::Type type) const {
  // make sure push_back does not hurt
  assert(x_l.empty() && x_u.empty());
  int did = FindDesign(type, true);
  const StdVector<DesignRegion>& cur_des = regions[did];
  const unsigned int nr = cur_des.GetSize();

  for(unsigned int r = 0; r < nr; r++){
    const DesignRegion& cur_reg = cur_des[r];
    const double rscaling = 1.0 / cur_reg.scale_design;
    const double translation = cur_reg.translate_design;
    if(cur_reg.constant == VARIABLE){
      const unsigned int u = cur_reg.base + cur_reg.elements;
      for(unsigned int s = cur_reg.base; s < u; s++){
        x_l.push_back((data[s].GetLowerBound() - translation) * rscaling);
        x_u.push_back((data[s].GetUpperBound() - translation) * rscaling);
      }
    }else if(cur_reg.constant == CONSTANT_PER_REGION || cur_reg.constant == CONSTANT_ON_ALL_REGIONS) { // in FIXED case nothing is done
      x_l.push_back((data[cur_reg.base].GetLowerBound() - translation) * rscaling);
      x_u.push_back((data[cur_reg.base].GetUpperBound() - translation) * rscaling);
      if(cur_reg.constant == CONSTANT_ON_ALL_REGIONS) break;
    }
  }
}

void DesignSpace::WriteBoundsToExtern(double* x_l, double* x_u) const {
  const unsigned int nd = design.GetSize();
  const unsigned int nr = regions[0].GetSize();
  unsigned int d = 0;
  for(unsigned int des = 0; des < nd; des++){
    const StdVector<DesignRegion>& cur_des = regions[des];
    for(unsigned int r = 0; r < nr; r++){
      const DesignRegion& cur_reg = cur_des[r];
      const double rscaling = 1.0 / cur_reg.scale_design;
      const double translation = cur_reg.translate_design;
      if(cur_reg.constant == VARIABLE){
        const unsigned int u = cur_reg.base + cur_reg.elements;
        for(unsigned int s = cur_reg.base; s < u; s++){
          x_l[d] = (data[s].GetLowerBound() - translation) * rscaling;
          x_u[d++] = (data[s].GetUpperBound() - translation) * rscaling;
        }
      }else if(cur_reg.constant == CONSTANT_PER_REGION || cur_reg.constant == CONSTANT_ON_ALL_REGIONS) { // in FIXED case nothing is done
        x_l[d] = (data[cur_reg.base].GetLowerBound() - translation) * rscaling;
        x_u[d++] = (data[cur_reg.base].GetUpperBound() - translation) * rscaling;
        if(cur_reg.constant == CONSTANT_ON_ALL_REGIONS) break;
      }
    }
  }
  assert(d == DesignSpace::GetNumberOfVariables());
}
void DesignSpace::WriteSparseGradientToExtern(StdVector<double>& out, DesignElement::ValueSpecifier vs, DesignElement::Access access, Function* f, bool use_scaling) const
{
  // Bastian did some complicated reordering stuff. For the only case of sparse Jacobians (slope constraints)
  // we'll have only the simple standard situation .. if this changes you have at least a test case :) Fabian
  // This should work now as long there is only one region. Jannis
  // had to weaken this condition for DESIGN_TRACKING in debug mode
  assert((regions[0].GetSize() == 1) || (f->GetType() != Function::DESIGN_TRACKING));
  assert(f != NULL); // only constraints can have sparse Jacobians
  assert(!f->IsObjective()); // only constraints can have sparse Jacobians
  
  unsigned int data_size = DesignSpace::GetNumberOfVariables(); // do not take aux variables

  StdVector<unsigned int>& sparsity = f->GetSparsityPattern();

  assert(out.window.GetSize() == sparsity.GetSize());
  unsigned int base = out.window.GetStart();
  for(unsigned int i = 0; i < sparsity.GetSize(); i++)
  {
    unsigned int s = sparsity[i];
    if(s <= data_size){ // else we have parts of the sparsity pattern on the aux design
      assert(out.InWindow(base + i));
      double scaling = use_scaling ? regions[FindDesign(data[s].GetType())][0].scale_design : 1.0;
      out[base + i] = data[sparsity[i]].GetValue(vs, access, f) * scaling;
    }
  }
}
void DesignSpace::WriteDenseGradientToExtern(StdVector<double>& out, DesignElement::ValueSpecifier vs, DesignElement::Access access, Function* f, bool use_scaling) const
{
  // this does now do reordering as gradients are reordered in the optimizer
  // must be set in the constructor! might be trivial volume fraction or from file!!
  assert(f != NULL);
  assert(!(vs == DesignElement::COST_GRADIENT && !f->IsObjective()));
  unsigned int n0 = out.window.GetStart(); // to grow up to the total number of design variables
  unsigned int n = n0;
  const unsigned int nd = design.GetSize();

  f->GetExcitation()->Apply(false); // this takes the proper gradient for robustness and transformation

  for(unsigned int des = 0; des < nd; des++)
  {
    const StdVector<DesignRegion>& cur_des = regions[des];
    const unsigned int nr = regions[0].GetSize();
    for(unsigned int r = 0; r < nr; r++)
    {
      const DesignRegion& cur_reg = cur_des[r];
      const double scaling = cur_reg.scale_design;
      const unsigned int u = cur_reg.base + cur_reg.elements;
      if(cur_reg.constant == VARIABLE)
      {
        for(unsigned int s = cur_reg.base; s < u; s++)
        {
          LOG_DBG3(designSpace) << "DS:WDGtE: non-constant region r=" << r << " rid=" << cur_reg.regionId << " out[" << n << "] = design[" << s << "]=" << data[s].GetValue(vs, access, f);
          assert(out.InWindow(n));
          out[n++] = data[s].GetValue(vs, access, f) * scaling;
        }
      }
      else if(cur_reg.constant == CONSTANT_PER_REGION || cur_reg.constant == CONSTANT_ON_ALL_REGIONS) // in FIXED case nothing is done
      {
        if(cur_reg.constant == CONSTANT_PER_REGION || r == 0){
          out[n] = 0;
        }
        for(unsigned int s = cur_reg.base; s < u; s++)
        {
          assert(out.InWindow(n));
          out[n] += data[s].GetValue(vs, access, f) * scaling;
          LOG_DBG3(designSpace) << "WriteDenseGradientToExtern: constant region " << r << ": out[" << n << "] += design[" << s << "]=" << data[s].GetValue(vs, access, f);
        }
        LOG_DBG3(designSpace) << "WriteDenseGradientToExtern: constant region " << r << ": sum = " << out[n] / scaling;
        if(cur_reg.constant == CONSTANT_PER_REGION || (cur_reg.constant == CONSTANT_ON_ALL_REGIONS && (r == nr-1) ) )
        {
          n++;
        }
      }
    }
  }
  assert(n - n0 == DesignSpace::GetNumberOfVariables());
}
void DesignSpace::Reset(DesignElement::ValueSpecifier vs, DesignElement::Type design)
{
  unsigned int start = (design == DesignElement::DEFAULT || design == DesignElement::MECH_TRACE) ? 0 : FindDesign(design) * elements;
  unsigned int end   = (design == DesignElement::DEFAULT || design == DesignElement::MECH_TRACE) ? data.GetSize() : start + elements;
  LOG_DBG3(designSpace) << "Reset: vs=" << DesignElement::valueSpecifier.ToString(vs) << " design="
                        << DesignElement::type.ToString(design) << " from " << start << " to " << end;

  // speed up by repeating loops
  switch(vs)
  {
  case DesignElement::DESIGN:
    for(unsigned int i = start; i < end; i++)
      data[i].SetDesign(0.0);
    break;
  case DesignElement::CONSTRAINT_GRADIENT:
  case DesignElement::COST_GRADIENT:
  case DesignElement::FUNCTION_GRADIENT:
    for(unsigned int i = start; i < end; i++)
      data[i].Reset(vs);
    break;
  default:
    if(end-start > 0)
      throw Exception("value specifier not handled");
  }
}

void DesignSpace::ResetGradient(Function* f)
{
  for(DesignElement& de : data)
    de.Reset(DesignElement::FUNCTION_GRADIENT, f);
}

inline
BaseDesignElement* DesignSpace::GetDesignElement(unsigned int idx)
{
  assert(idx < data.GetSize());
  return dynamic_cast<BaseDesignElement*>(&data[idx]);
}

DesignElement* DesignSpace::Find(unsigned int elemNum, DesignElement::Type dt, bool throw_exception, bool include_pseudo_designs, int mm_index)
{
  int idx = Find(elemNum, throw_exception, include_pseudo_designs);
  if(idx == -1)
    return NULL; // an exception was already thrown if desired
  // check for real design or pseudo design
  if(elemToDesign[elemNum].second == true)
  {
    for(unsigned int d = 0, nd = design.GetSize(); d < nd; d++)
    {
      DesignElement& de = data[elements * d + idx];
      if(de.GetType() == dt)
      {
        assert(de.elem->elemNum == elemNum);
        assert(mm_index == -1 || (mm_index >= 0 && de.multimaterial != NULL));

        if(mm_index < 0 || de.multimaterial == NULL || mm_index == de.multimaterial->index)
          return &de;
      }
    }
  }
  else
  {
    for(unsigned int i = 0, n = pseudoDesigns_.GetSize(); i < n; i++)
    {
      // LOG_DBG3(designSpace) << "Find e=" << elemNum << " pd! idx=" << idx << " dt=" << dt << " i=" << i << " test_e=" << pseudoDesigns_[i][idx].elem->elemNum << " test_d=" << pseudoDesigns_[i][idx].GetType();
      // there might be different regions within different pseudo designs!
      if(idx < (int) pseudoDesigns_[i].GetSize() && pseudoDesigns_[i][idx].elem->elemNum == elemNum && pseudoDesigns_[i][idx].GetType() == dt)
        return &(pseudoDesigns_[i][idx]);
    }
  }
  if(throw_exception)
    throw Exception("design type not in design or pseudo design region problem");
  return NULL;
}

void DesignSpace::ToInfo(ErsatzMaterial* em)
{
  // em == null for the case we create this design space only for loading ersatz material within a simulation
  PtrParamNode in = em != NULL ? info_ : domain->GetInfoRoot()->Get("loadErsatzMaterial");

  PtrParamNode tf = in->Get("transferFunctions");
  for(unsigned int i = 0; i < transfer.GetSize(); i++)
    transfer[i].ToInfo(tf->Get("transferFunction", ParamNode::APPEND));

  if(!transform.IsEmpty())
  {
    PtrParamNode t = in->Get("transform");
    for(unsigned int i = 0; i < transform.GetSize(); i++)
      transform[i].ToInfo(t);
  }

  PtrParamNode dv = in->Get("designVariables");
  dv->Get("opt_vars")->SetValue(GetNumberOfVariables());
  // TODO @Bastian - add you shape stuff if you like
  dv->Get("field_elems")->SetValue(data.GetSize()); // density, angle, ... -> no variables for featured designs
  dv->Get("field_fe_elemes")->SetValue(GetNumberOfElements()); // fe-elems the field elems are assigned to
  for(unsigned int i = 0; i < design.GetSize(); i++)
  {
    DesignElement& de = data[i * elements];
    // FIXME an arbitrary transfer function is nonsense!
    if(design[this->FindDesign(de.GetType())].relative_bound > 0.) {
      dv->Get("design", ParamNode::APPEND)->Get("relative_bound")->SetValue(design[this->FindDesign(de.GetType())].relative_bound);
    }
    de.ToInfo(dv->Get("design", ParamNode::APPEND), GetTransferFunction(de.GetType(), App::MECH, false), em); // silent!
  }

  in->Get("pamping")->SetValue(pamping_);
  in->Get("regular")->SetValue(IsRegular());

  if(elementCache != NULL)
    elementCache->ToInfo(in->Get("localElementCache"));
  else
    in->Get("localElementCache")->SetValue("disabled");

  if(regions.GetSize() > 0)
  {
    PtrParamNode rs = in->Get("regions");
    for(unsigned int i = 0; i < regions[0].GetSize(); i++)
      regions[0][i].ToInfo(rs->Get("region", ParamNode::APPEND));
  }

  for(const DesignMaterial* dm : designMaterials)
    dm->ToInfo(in->Get("designMaterial", ParamNode::APPEND));
}

string DesignSpace::ToString(int level)
{
  std::stringstream ss;

  assert(level == 0 || level == 1);
  if(level == 0)
  {
    ss << "design[";
    for(unsigned int i = 0; i < data.GetSize(); i++)
    {
      DesignElement* de = &data[i];
      ss << i << ":elem=" << de->elem->elemNum;
      if(de->vicinity != NULL) ss << " " << de->vicinity->ToString();
      ss << " ";
    }
    ss << "]";
  }
  if(level == 1)
  {
    for(unsigned int i = 0; i < data.GetSize(); i++)
      ss << data[i].GetPlainDesignValue() << ", ";
  }
  return ss.str();
}
void DesignSpace::DisableTransferFunctions()
{
  for(unsigned int i = 0; i < transfer.GetSize(); i++) transfer[i].Enable(false);
}
 /** Enables the transfer functions -> sets again to the xml settings after
  * temporarily disabled via DisableTransferFunctions() */
void DesignSpace::EnableTransferFunctions()
{
  for(unsigned int i = 0; i < transfer.GetSize(); i++) transfer[i].Enable(true);
}

unsigned int DesignSpace::GetNumberOfVariables() const
{
  const unsigned int nd = design.GetSize();
  const unsigned int nr = regions.GetSize() > 0 ? regions[0].GetSize() : 0; // e.g. for pure shape optimization
  unsigned int n = 0;
  for(unsigned int des = 0; des < nd; des++){
    for(unsigned int r = 0; r < nr; r++){
      const DesignRegion& cur_reg = regions[des][r];
      switch(cur_reg.constant){
      case VARIABLE:
        n += cur_reg.elements;
        break;
      case CONSTANT_PER_REGION:
        n++;
        break;
      case CONSTANT_ON_ALL_REGIONS:
        if(r == 0) n++;
        break;
      case FIXED:
        break;
      }
    }
  }
  return n;
}

string DesignSpace::DumpRegions() const
{
  std::stringstream ss;
  for(unsigned int d = 0; d < regions.GetSize(); d++)
    for(unsigned int r = 0; r < regions[d].GetSize(); r++)
      ss << "des=" << d << " reg=" << r << ": " << regions[d][r].ToString() << " ";
  return ss.str();
}

int DesignSpace::FindRegion(RegionIdType regionId) const
{
  if(regions.GetSize() > 0)
  {
    const StdVector<DesignRegion>& regs = regions[0];
    const unsigned int nr = regs.GetSize();
    for(unsigned int r = 0; r < nr; r++)
    {
      if(regs[r].regionId == regionId)
        return r;
    }
  }
  return -1;
}

template <class T>
void DesignSpace::ExtractResults(shared_ptr<BaseResult> base_result)
{
  // in the load ersatz material case the context is not set
  assert(Optimization::context);
  Optimization::context->Update();

  // only for rare generic cases
  PrepareSpecialResults();

  // our results are up to now scalar!
  Result<T>& result = dynamic_cast<Result<T> &>(*base_result);
  // the description of the result
  shared_ptr<ResultInfo> ri = result.GetResultInfo();
  // Work with a result description. This is either a result description from the
  // xml file or when using the "predefined" *_PSEUDO_* we set it here.
  ResultDescription def;
  // set the defaults to be maybe replaced by a resultDescription
  def.solutionType = ri->resultType;
  // this is clearly nonsense if the result/solution type is OPT_RESULT_*
  switch(ri->resultType)
  {
  case PSEUDO_DENSITY:
  case MECH_PSEUDO_DENSITY:
  case PHYSICAL_PSEUDO_DENSITY:
  case ELEC_PHYSICAL_PSEUDO_DENSITY:
    def.design = DesignElement::DENSITY;
    break;
  case ELEC_PSEUDO_POLARIZATION:
    def.design = DesignElement::POLARIZATION;
    break;
  case RHS_PSEUDO_DENSITY:
  case PHYSICAL_RHS_PSEUDO_DENSITY:
    def.design = DesignElement::RHS_DENSITY;
    break;
  default:
    // to be overwritten by the ResultDescription
    def.design = DesignElement::DENSITY;
    break;
  }
  // somehow critical! but only for density filtering, if at all.
  def.access = DesignElement::IsPhysical(ri->resultType) ? DesignElement::SMART : DesignElement::PLAIN;
  def.value  = DesignElement::DESIGN;
  // ignore defaults if there is a result description for the OPT_RESULT_* case
  for(unsigned int i = 0; i < resultDescriptions.GetSize(); i++)
    if(resultDescriptions[i].solutionType == ri->resultType)
      def = resultDescriptions[i];

  LOG_DBG(designSpace) << "ER: def=" << def.ToString();

  // this enables excitation specific physical designs (robust, transformation)
  if(def.excitation >= 0 && domain->GetOptimization() != NULL)
  {
    StdVector<Excitation>& mex = domain->GetOptimization()->GetMultipleExcitation()->excitations;
    if(def.excitation >= (int) mex.GetSize())
      EXCEPTION("'result' has too large 'excitation' index " << def.excitation << " for only " << mex.GetSize() << " excitations");

    mex[def.excitation].Apply(false);
    LOG_DBG(designSpace) << "ER: apply excitation " << mex[def.excitation].GetFullLabel();
  }

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
    actSol[it.GetPos()] = GetNodalValue(node, descr.value, descr);
  }
}

template <class T>
void DesignSpace::FillElementResults(Result<T>& result, ResultDescription& descr)
{
  Vector<T>& result_data = result.GetVector();
  // this is our entity result, a scalar or a vector of dim 2/3
  unsigned int dofs = result.GetResultInfo()->dofNames.GetSize();
  assert(dofs >= 1 && dofs <= 3);
  StdVector<double> result_value(dofs);
  // loop over elements from result. We have to do it this way as the the connection
  // of design element and result element is the element(->elemNum) but we cannot
  // search in the result for an element.
  EntityIterator it = result.GetEntityList()->GetIterator();
  // set the result as we need it
  result_data.Resize(result.GetEntityList()->GetSize() * dofs);
  SolutionType st = result.GetResultInfo()->resultType;

  // the value when we are not in a design domain
  // the default value is 0.0 but 1 for densities
  double none = 0.0;
  switch(st)
  {
  case MECH_PSEUDO_DENSITY:
  case PSEUDO_DENSITY:
  case PHYSICAL_PSEUDO_DENSITY:
  case ELEC_PSEUDO_POLARIZATION:
  case ELEC_PHYSICAL_PSEUDO_DENSITY:
  case RHS_PSEUDO_DENSITY:
  case PHYSICAL_RHS_PSEUDO_DENSITY:
    none = 1.0;
    break;
  default:
    break;
  }
  // search where in data we are. -1 when not found
  int base = (st == MECH_ELEM_VOL || st == MECH_ELEM_POROSITY) ? 0 : FindDesign(descr.design, false);
  Excitation* ex = domain->GetOptimization() != NULL ? Optimization::context->GetExcitation() : NULL;

  for (it.Begin(); !it.IsEnd(); it++)
  {
    // for elements not in the design region we set to the default value
    for(unsigned int i = 0; i < dofs; i++)
      result_value[i] = none;
    if(base >= 0 && FindRegion(it.GetElem()->regionId) >= 0)
    {
      // note that the index is from the first design set!
      unsigned int base_index = Find(it.GetElem()->elemNum, true, true); // exception and pseudo designs (?)
      // base=0 is first!
      unsigned int data_index = (base * elements) + base_index;
      DesignElement* org = &data[data_index];

      // we need to transform manually only for smart design with excitation given. The physicalPseudoDensity has it by itself
      if(descr.solutionType >= OPT_RESULT_1 && descr.solutionType < OPT_RESULT_BOUND && descr.access == DesignElement::SMART && descr.excitation >= 0 && ex != NULL && ex->transform != NULL)
      {
        DesignElement* trans = ApplyTransformations(org, org, NULL);
        trans->GetValue(descr, result_value, dofs);
      } else if (st == MECH_ELEM_VOL) {
        result_value = org->CalcVolume();
      } else if (st == MECH_ELEM_POROSITY) {
        result_value = org->GetElemPorosity();
      }
      else
        org->GetValue(descr, result_value, dofs);

      #ifdef CHECK_INDEX
        if(org->elem->elemNum != it.GetElem()->elemNum)
          EXCEPTION("mixed up indices:" << org->elem->elemNum << "!=" << it.GetElem()->elemNum
              << " base_index=" << base_index << " data_index=" << data_index << " it.Pos()=" << it.GetPos());
      #endif
    }
    else
    {
      // there might be the case that the function is not within the design space, e.g. stress
      // is only implemented for OPT_RESULT_x
      int ori = DesignElement::GetOptResultIndex(descr.solutionType);
      for(unsigned int f = 0; ori != -1 && f < pseudoDesigns_.GetSize(); f++)
      {
        StdVector<DesignElement>& data = pseudoDesigns_[f];
        // has the pseudo design the right special result?
        if(it.GetElem()->regionId == data[0].elem->regionId)
        {
          // search it slowly and add it up -> it will be a very special case anyway
          for(unsigned int e = 0; e < data.GetSize(); e++)
            if(data[e].elem == it.GetElem())
            {
              // make sure the result description is unique and we don't overwrite
              assert(result_value[0] == 0.0);
              data[e].GetValue(descr, result_value, dofs);
            }
        }
      }
    }
    for(unsigned int i = 0; i < dofs; i++)
      result_data[it.GetPos() * dofs + i] = result_value[i];
  }
}

DesignSpace::DesignRegion* DesignSpace::GetRegion(RegionIdType id, MultiMaterial* mm, bool throw_exception)
{
  assert(mm != NULL);
  return GetRegion(id, DesignElement::MULTIMATERIAL, mm->index, throw_exception);
}


DesignSpace::DesignRegion* DesignSpace::GetRegion(RegionIdType id, DesignElement::Type dt, int multimaterial_index, bool throw_exception)
{
  assert(!((dt == DesignElement::MULTIMATERIAL) && multimaterial_index < 0));

  for(unsigned int d = 0, dn = regions.GetSize(); d < dn; d++)
  {
    StdVector<DesignRegion>& regs = regions[d];
    if((dt != DesignElement::NO_TYPE && dt != DesignElement::ALL_DESIGNS) && regs[0].design != dt)
      continue;
    if((dt == DesignElement::MULTIMATERIAL) && regs[0].design == dt && regs[0].multimaterial->index != multimaterial_index)
      continue;
    for(unsigned r = 0, rn = regs.GetSize(); r < rn; r++)
    {
      if(regs[r].regionId == id)
        return &regs[r];
    }
  }

  if(throw_exception)
    EXCEPTION("cannot find design region");
  return NULL;
}


DesignSpace::DesignRegion* DesignSpace::GetRegion(RegionIdType id, MaterialClass mc, MaterialType mt, bool throw_exception)
{
  for(unsigned int d = 0, dn = regions.GetSize(); d < dn; d++)
    for(unsigned r = 0, rn = regions[d].GetSize(); r < rn; r++)
    {
      DesignRegion& dr = regions[d][r];
      if(dr.regionId == id)
      {
        LOG_DBG2(designSpace) << "GR: mc=" << mc << " mt=" << mt << " dr_sm=" << dr.scnd_materials.count(mc) <<  
        " dr_sm_mc=" <<(dr.scnd_materials.count(mc) > 0 ? dr.scnd_materials[mc].count(mt) : -1);

        for(auto mymc: dr.scnd_materials)
          for(auto mymt: mymc.second)
            {LOG_DBG2(designSpace) << "GR: mymc=" << mymc.first << " mymt=" << mymt.first << " mymt.name=" << MaterialTypeEnum.ToString(mymt.first);}

        if(dr.scnd_materials.count(mc) > 0 && dr.scnd_materials[mc].count(mt) > 0)
          return &dr;
      }
    }

  if(throw_exception)
    throw Exception("Cannot find design region " + std::to_string(id) + " for MaterialType::" + MaterialTypeEnum.ToString(mt));
  return NULL;
}

StdVector<DesignSpace::DesignRegion*> DesignSpace::GetRegions(DesignElement::Type dt)
{
  StdVector<DesignRegion*> res;
  res.Reserve(regions.GetSize() * regions[0].GetSize());
  for(StdVector<DesignRegion>& des : regions)
    if(dt == DesignElement::ALL_DESIGNS || des[0].design == dt)
      for(DesignRegion& reg : des)
        res.Push_back(&reg);

  return res;
}


DesignSpace::DesignRegion::DesignRegion()
{
  regionId = -1;
  multimaterial = NULL;
}

string DesignSpace::DesignRegion::ToString() const
{
  std::stringstream ss;
  ss << " d=" << DesignElement::type.ToString(design) << " reg=" << regionId << " base=" << base;
  ss << " elem=" << elements << " sd=" << scale_design << " td=" << translate_design << " scnd_material=" << scnd_material;
  //ss << " dc=" << DesignSpace::designConstant.ToString(constant);
  return ss.str();
}

void DesignSpace::DesignRegion::SetBiMaterial(const std::string& material)
{
  assert(!has_grndmat);
  scnd_material = material;
  has_bimat = true;
}

void DesignSpace::DesignRegion::SetGroundMaterial(const std::string& material)
{
  assert(!has_bimat);
  scnd_material = material;
  has_grndmat = true;
}

const DesignSpace::DesignRegion::CoefValue& DesignSpace::DesignRegion::GetCachedValue(MaterialClass mc, MaterialType mt) const
{
  assert(scnd_material_cached.size() > 0);
  assert(scnd_material_cached.find(mc) != scnd_material_cached.end());
  const auto& mtv = scnd_material_cached.at(mc);
  
  assert(mtv.size() > 0);
  assert(mtv.find(mt) != mtv.end());
  return mtv.at(mt);
}

void DesignSpace::SetupMultiMaterial(ParamNodeList design_list)
{
  assert(multimaterial.IsEmpty());

  for(unsigned int d = 0; d < design_list.GetSize(); d++)
  {
    PtrParamNode pn = design_list[d];
    DesignElement::Type dt = DesignElement::type.Parse(pn->Get("name")->As<string>());
    if(dt == DesignElement::MULTIMATERIAL)
    {
      if(!pn->Has("material"))
        throw Exception("mutlimaterial designs require the 'material' attribute");
      string material = pn->Get("material")->As<string>();

      // check for material
      for(unsigned int m = 0; m < multimaterial.GetSize(); m++)
      {
        if(multimaterial[m].name == material)
          throw Exception("multimaterial design " + material + " not unique");
      }

      if(pn->Get("region")->As<string>() != "all")
        throw Exception("multimaterial not yet compatible with multiregion");

      multimaterial.Push_back(MultiMaterial());
      MultiMaterial& mm = multimaterial.Last();
      mm.name = material;
      mm.index = multimaterial.GetSize() - 1;
      // the real material id set on the fly by GetMultiMaterial()
    }
    else if(pn->Has("material"))
      throw Exception("the 'design' attribute 'material' is only for multimaterial designs");
  }
}

PtrCoefFct DesignSpace::DesignRegion::GetScndMaterial(MaterialClass mc, MaterialType mt, SinglePDE* pde)
{
  assert(scnd_material != ""); // check with HasBiMaterial()!

  if(scnd_materials.count(mc) == 0 || scnd_materials[mc].count(mt) == 0)
  {
    // This part is only called once, because then we have saved it in the map
    #pragma omp critical (DR_GMB)
    {
      // apparently first run
      MaterialHandler* matLoader = domain->GetMaterialHandler();
      BaseMaterial* mat = matLoader->LoadMaterial(scnd_material, mc);
      MathParser* mp = domain->GetMathParser();

      switch(mt)
      {
      case DENSITY:
        if (mc != ACOUSTIC){
          scnd_materials[mc][mt] = mat->GetScalCoefFnc(DENSITY,Global::REAL);
        } else {
          // because in acoustics the density factor is 1/density
          scnd_materials[mc][mt] = CoefFunction::Generate(mp, Global::COMPLEX, CoefXprBinOp(mp, "1.0", mat->GetScalCoefFnc(DENSITY,Global::COMPLEX), CoefXpr::OP_DIV));
        }
        break;

      case ACOU_BULK_MODULUS:
      {
        // because in acoustics the bulk modulus factor is 1/bulk
        scnd_materials[mc][mt] = CoefFunction::Generate(mp, Global::COMPLEX, CoefXprBinOp(mp, "1.0", mat->GetScalCoefFnc(ACOU_BULK_MODULUS,Global::COMPLEX), CoefXpr::OP_DIV));
        break;
      }

      case MECH_STIFFNESS_TENSOR:
      {
        SinglePDE* sp = pde ? pde : domain->GetSinglePDE("mechanic");
        if(Optimization::context->DoBloch() || sp->HasComplexMatData(regionId))
          scnd_materials[mc][mt] = mat->GetTensorCoefFnc(MECH_STIFFNESS_TENSOR, sp->GetSubTensorType(), Global::COMPLEX);
        else
          scnd_materials[mc][mt] = mat->GetTensorCoefFnc(MECH_STIFFNESS_TENSOR, sp->GetSubTensorType(), Global::REAL);
        break;
      }

      case HEAT_CONDUCTIVITY_TENSOR:
      {
        SubTensorType stt = pde ? pde->GetSubTensorType() : domain->GetSinglePDE("heatConduction")->GetSubTensorType();
        scnd_materials[mc][mt] = mat->GetTensorCoefFnc(HEAT_CONDUCTIVITY_TENSOR, stt, Global::REAL);
        break;
      }

      case HEAT_CAPACITY:
        scnd_materials[mc][mt] = mat->GetScalCoefFnc(HEAT_CAPACITY,Global::REAL);
        break;

      default:
        assert(false);
      }
    } // omp critical
  }
  LOG_DBG3(designSpace) << "DS:DR:GBM: mc=" << mc << " mt=" << MaterialTypeEnum.ToString(mt) << " -> " << scnd_materials[mc][mt]->ToString();
  return scnd_materials[mc][mt];
}

PtrCoefFct DesignSpace::DesignRegion::GetSncdMaterial(const string& integrator)
{
  if(integrator == "LinElastInt")
    return GetScndMaterial(MECHANIC, MECH_STIFFNESS_TENSOR);
  if(integrator == "MassInt")
    return GetScndMaterial(MECHANIC, DENSITY);
  if(integrator == "HeatConductivity")
    return GetScndMaterial(THERMIC, HEAT_CONDUCTIVITY_TENSOR);
  if(integrator == "LaplaceIntegrator")
    return GetScndMaterial(ACOUSTIC, DENSITY);
  if(integrator == "MassIntegrator")
    return GetScndMaterial(ACOUSTIC, ACOU_BULK_MODULUS);
  assert(false);
  return PtrCoefFct();
}


void DesignSpace::DesignRegion::ToInfo(PtrParamNode node) const
{
  node->Get("name")->SetValue(domain->GetGrid()->GetRegion().ToString(regionId));
  node->Get("elements")->SetValue(elements);
  assert(!(HasBiMaterial() && HasGroundMaterial()));
  node->Get("bimaterial")->SetValue(HasBiMaterial() ? scnd_material : "-");
  node->Get("groundmaterial")->SetValue(HasGroundMaterial() ? scnd_material : "-");
}

const string DesignSpace::ToForm(MaterialClass mc, MaterialType mt)
{
  switch(mc)
  {
  case MECHANIC:
    switch(mt)
    {
    case MECH_STIFFNESS_TENSOR:
      return "LinElastInt";
    case DENSITY:
      return "MassInt";
    default:
      assert(false);
    }
    break;

  case THERMIC:
    switch(mt)
    {
    case HEAT_CONDUCTIVITY_TENSOR:
      return "HeatConductivity";
    default:
      assert(false);
    }
    break;
  default:
    assert(false);
  }
  return ""; // shall not happen
}

void MultiMaterial::ToInfo(PtrParamNode in)
{
  Matrix<double> E;
  BaseMaterial* bm = NULL;

  SubTensorType stt = Optimization::context->pde->GetSubTensorType();

  assert(!Optimization::context->DoMultiSequence());

  switch(Optimization::context->mat->GetSystem())
  {
  case OptimizationMaterial::PIEZOCOUPLING:
    bm = GetMultiMaterial(ELECTROSTATIC);
    bm->GetTensor(E, BaseMaterial::ConvertMaterialClass(ELECTROSTATIC), Global::REAL, stt);
    in->Get("electrostatic")->SetValue(E);

    bm = GetMultiMaterial(PIEZO);
    bm->GetTensor(E, BaseMaterial::ConvertMaterialClass(PIEZO), Global::REAL, stt);
    in->Get("piezo")->SetValue(E);

    // no break by intention
  case OptimizationMaterial::MECH:
    bm = GetMultiMaterial(MECHANIC);
    bm->GetTensor(E, BaseMaterial::ConvertMaterialClass(MECHANIC), Global::REAL, stt);
    in->Get("mechanic")->SetValue(E);
    break;

  default:
    assert(false);
    break;
  }
}

BaseMaterial* MultiMaterial::GetMultiMaterial(const MaterialClass mc)
{
  for(unsigned int m = 0; m < material.GetSize(); m++)
    if(material[m].second == mc)
      return material[m].first;

  // apparently first run
  MaterialHandler* matLoader = domain->GetMaterialHandler();
  BaseMaterial* mat = matLoader->LoadMaterial(name, mc);
  material.Push_back(std::make_pair(mat, mc));
  return mat;
}



void DensityFilterMat::AssembleFilterMatrix(StdVector<DesignElement>&data, int sum_neighbours,int filter_idx, unsigned int start, unsigned int end){

  // We just get all the design elements and for each filter create a sparse matrix
  // For the sparse matrix we require row_index(element number) , column index(neighbour idx), and weights array
  // Implementing this above in the neigbhor search will require use of critical sections. So lets just stick to looping over all elements and extracting

    int num_elem = data.GetSize();
    if (end -start > 0) {
      num_elem = end - start;
    }
    int nnz = (sum_neighbours+num_elem);
    this->filter_mat.SetSize(num_elem,num_elem,nnz);

    UInt *colPointer=this->filter_mat.GetColPointer();
    UInt *rowPointer=this->filter_mat.GetRowPointer();
    double *dataPtr=this->filter_mat.GetDataPointer();

    this->filtered_vec.Resize(num_elem);
    this->inv_weighted_sum.Resize(num_elem);

    int lastIndex=0;
    rowPointer[0]=lastIndex;
    if (!(end -start > 0)) {
      start = 0;
      end = data.GetSize();
    }
    for (UInt i=start;i < end; i++){

      auto neighbours = data[i].simp->filter[filter_idx].neighborhood;
      // Set this weight sum so that we don't recalculate it
      data[i].simp->filter[filter_idx].weight_sum = (data[i].simp->filter[filter_idx].CalcWeightSum(true));
      this->inv_weighted_sum[i-start] = (1/ data[i].simp->filter[filter_idx].weight_sum);
      colPointer[lastIndex]=i-start;
      dataPtr[lastIndex]= data[i].simp->filter[filter_idx].weight  * this->inv_weighted_sum[i-start];
      for (UInt j=0;j<neighbours.GetSize();j++){
        colPointer[lastIndex+j+1]=neighbours[j].neighbour->GetIndex()-start;
        dataPtr[lastIndex+j+1]=(neighbours[j].weight)* this->inv_weighted_sum[i-start];
      }
      lastIndex +=(neighbours.GetSize()+1); // Since Neighbours doesn't include the own element
      rowPointer[i+1-start]=lastIndex;
    }
}



void DensityFilterMat::CacheDensityFilteredValue(const Vector<double>& design_vec)
{
  this->filter_mat.Mult(design_vec, this->filtered_vec);
}

void DensityFilterMat::ExportDensityFilterMatrix(string filename){
  this->filter_mat.ExportMatrixMarket(filename,"filter_matrix");
}


// explicit template instantiation
template void DesignSpace::ExtractResults<double>(shared_ptr<BaseResult> base_result);
template void DesignSpace::ExtractResults<complex<double> >(shared_ptr<BaseResult> base_result);
template void DesignSpace::FillNodeResults<double>(Result<double>& result, ResultDescription& descr);
template void DesignSpace::FillNodeResults<complex<double> >(Result<complex<double> >& result, ResultDescription& descr);
template void DesignSpace::FillElementResults<double>(Result<double>& result, ResultDescription& descr);
template void DesignSpace::FillElementResults<complex<double> >(Result<complex<double> >& result, ResultDescription& descr);
template bool DesignSpace::ApplyPhysicalDesign<double>(const CoefFunctionOpt* coef, Matrix<double>& retMat, const LocPointMapped* lpm);
template bool DesignSpace::ApplyPhysicalDesign<double>(const CoefFunctionOpt* coef, Vector<double>& retVEc, const LocPointMapped* lpm);
template bool DesignSpace::ApplyPhysicalDesign<complex<double> >(const CoefFunctionOpt* coef, Matrix<complex<double> >& retMat, const LocPointMapped* lpm);
template bool DesignSpace::ApplyPhysicalDesign<double>(const CoefFunctionOpt* coef, double& retScal, const LocPointMapped* lpm);
template bool DesignSpace::ApplyPhysicalDesign<complex<double> >(const CoefFunctionOpt* coef, complex<double>& retScal, const LocPointMapped* lpm);
template bool DesignSpace::ApplyPhysicalDesignElementMatrix<double>(BiLinearForm* form, Matrix<double>& retMat, const Elem* elem);
template bool DesignSpace::ApplyPhysicalDesignElementMatrix<complex<double> >(BiLinearForm* form, Matrix<complex<double> >& retMat, const Elem* elem);
template bool DesignSpace::TestTensorPosDef<double>(Matrix<double>& retMat, const LocPointMapped* lpm, DesignElement::Type direction);
template bool DesignSpace::TestTensorPosDef<complex<double> >(Matrix<complex<double> >& retMat, const LocPointMapped* lpm, DesignElement::Type direction);
