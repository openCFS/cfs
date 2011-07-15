#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/ShapeDesign.hh"
#include "Optimization/TransferFunction.hh"
#include "Optimization/Condition.hh"
#include "Optimization/Optimizer/BaseOptimizer.hh"
#include "Optimization/Optimizer/ShapeOptimizer.hh"
#include "Optimization/LevelSet.hh"
#include "Optimization/OptimizationMaterial.hh"
#include "General/exception.hh"
#include "General/Enum.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"
#include "Domain/resultInfo.hh"
#include "DataInOut/MaterialHandler.hh"
#include "Utils/result.hh"
#include "Utils/StdVector.hh"
#include "PDE/SinglePDE.hh"
#include "DataInOut/Logging/cfslog.hh"

#include <boost/lexical_cast.hpp>

using namespace CoupledField;

using std::complex;

// declare class specific logging stream
DECLARE_LOG(designSpace)
DEFINE_LOG(designSpace, "designSpace")

// declare class specific logging stream
DECLARE_LOG(ersatz)
DEFINE_LOG(ersatz, "ersatzMaterialFactor")



DesignSpace::DesignSpace(StdVector<RegionIdType>& reg_data, PtrParamNode pn, ErsatzMaterial::Method method)
{
  LOG_TRACE(designSpace) << "DesignSpace for regions=" << reg_data;
  all_regions_regular_ = domain->GetGrid()->IsRegionRegular(reg_data);
  
  // for convenience
  regionIds_ = reg_data;

  design_id = 0;
  optimizer_ = NULL;

  designMaterial = NULL;

  applicationForm.SetName("DesignSpace::ApplicationForm");
  applicationForm.Add(Optimization::ELEC, "linGradBDBInt");
  applicationForm.Add(Optimization::MECH, "linElastInt");
  // We follow for the stress, strain calculation the transfer functions of mech
  applicationForm.Add(Optimization::MECH, "MechStressStrain", false);
  applicationForm.Add(Optimization::PIEZO_COUPLING, "linPiezoCoupling");
  applicationForm.Add(Optimization::CHARGE_DENSITY, "LinNeumannInt");
  applicationForm.Add(Optimization::PRESSURE, "PressureLinForm");
  applicationForm.Add(Optimization::MASS, "MassInt");
  // acoustic and heat
  applicationForm.Add(Optimization::LAPLACE, "LaplaceInt");

  // read the elements
  elements = domain->GetGrid()->GetNumElems(reg_data);
  
  pamping_ = pn->Has("pamping") ? pn->Get("pamping/value")->As<double>() : 0.0;

  ParamNodeList trans_in = pn->GetList("transferFunction");

  // should we adapt the lower bound of the design variable to the penalty exponent of the transfer function? 
  bool adapt_lower(false);

  // check for non-design-vicinity
  non_design_vicinity_ = pn->Has("designSpace") ? pn->Get("designSpace/non_design_vicinity")->As<bool>() : false;

  // store the CFS element (number) to design element mapping.
  // Used by Find() and the filter and vicinity neighbors
  elemToDesign.Resize(domain->GetGrid()->GetNumElems() + 1, std::make_pair(-1, true)); // 1 based.

  // if only shape opt is done, ignore the rest here, the DesignSpace is "empty"
  if(method == ErsatzMaterial::SHAPE_OPT){
    return;
  }
  
  // number of different designs
  ParamNodeList pn_design = pn->GetList("design");
  unsigned int nd = 0;
  for(unsigned int d = 0; d < pn_design.GetSize(); d++){
    DesignElement::Type dt = DesignElement::type.Parse(pn_design[d]->Get("name")->As<std::string>());
    if(design.Find(dt) < 0){
      design.Push_back(dt);
      nd++;
    }
  }
  
  if(elements == 0 || nd == 0)
  { // this may happen in shape optimization
    if(method !=ErsatzMaterial::SHAPE_GRAD)
    {
      if(elements == 0) throw Exception("empty regions");
      if(nd == 0) throw Exception("no designs types given.");
    }
  }
  else // 'standard' SIMP case
  {
    // set our own structure with is element times design parameters
    data.Reserve(elements * nd);
    totalElements_.Reserve(elements  * nd); // the quick access copy which also combines pseudo design elements
    
    const unsigned int nr = reg_data.GetSize();

    // check whether all regions have all design variables and all design variables are given in all regions
    StdVector<bool> region_design;
    region_design.Resize(nd * nr, true);

    StdVector<std::string> regNames;
    domain->GetGrid()->GetRegionNames(regNames);

    // only temporary
    StdVector<Elem*> elems;

    // set size of regions
    regions.Resize(nd);

    // iterate over all designs, we have to always have the same order of designs not the order from the xml file
    for(unsigned int dti = 0; dti < nd; dti++)
    {
      regions[dti].Resize(reg_data.GetSize());
      
      DesignElement::Type dt = design[dti];
      
      for(unsigned int d = 0; d < pn_design.GetSize(); d++)
      {
        PtrParamNode curr_design_pn = pn_design[d];
        if(DesignElement::type.Parse(curr_design_pn->Get("name")->As<std::string>()) != dt) continue;

        std::string design_reg = curr_design_pn->Get("region")->As<std::string>();
        std::string design_bim = curr_design_pn->Has("bimaterial") ? curr_design_pn->Get("bimaterial")->As<std::string>() : "";
        // there is no default for adapt_lower in the load ersatz material case
        adapt_lower = false;
        if(curr_design_pn->Has("adapt_lower"))
          adapt_lower = curr_design_pn->Get("adapt_lower")->As<bool>();
        
        for(unsigned int r = 0; r < nr; r++)
        {
          domain->GetGrid()->GetElems(elems, reg_data[r]);
          unsigned int n = elems.GetSize();
          std::string reg = regNames[reg_data[r]];
          
          bool design_all = design_reg == "all";
          if(design_all || design_reg == reg){
            if(!region_design[r*nd + dti]){
              throw Exception("Design/Region combination given twice!");
            }
            region_design[r*nd + dti] = false;
            
            // this is now done for all designs per region, this is called for every design and every region here
            regions[dti][r].design = dt;
            regions[dti][r].regionId = reg_data[r];
            regions[dti][r].base = r == 0 ? dti * elements : regions[dti][r-1].base + regions[dti][r-1].elements;
            regions[dti][r].elements = n;

            regions[dti][r].constant = VARIABLE;
            if(curr_design_pn->Has("constant") && curr_design_pn->Get("constant")->As<bool>()){
              // we have a constant densign-value on that region
              regions[dti][r].constant = design_all ? CONSTANT_ON_ALL_REGIONS : CONSTANT_PER_REGION;
            }
            if(curr_design_pn->Get("fixed")->As<bool>()){ // fixed overwrites all other settings
              regions[dti][r].constant = FIXED;
            }
            
            regions[dti][r].scale_design = 1.0;
            regions[dti][r].translate_design = 0.0;
            if(design_bim != "")
              regions[dti][r].SetBiMaterial(design_bim);

            if(curr_design_pn->Has("scale") && curr_design_pn->Get("scale")->As<bool>()){
              double upper = curr_design_pn->Get("upper")->As<double>();
              double lower = curr_design_pn->Get("lower")->As<double>();
              regions[dti][r].scale_design = (upper - lower);
              regions[dti][r].translate_design = lower;
            }

            double initial = curr_design_pn->Get("initial")->As<double>();
            LOG_DBG2(designSpace) << "add design " << dt << ":" << DesignElement::type.ToString(dt)  
                                  << " initial=" << boost::lexical_cast<std::string>(initial);
                        


            for(unsigned int e = 0; e < n; e++)
            {
              DesignElement de(curr_design_pn, elems[e], data.GetSize());
              // simple extension: if the passed value is smaller
              // than 0.0 (which makes no sense as the element contribution
              // is always positive), we put a random number -> random initial design
              de.SetDesign(initial < 0.0 ? ((rand()%100 + 10)/110.0) : initial);
              
              data.Push_back(de);
              totalElements_.Push_back(&data.Last());

              // append rucksack :)
              if(method == ErsatzMaterial::SIMP_METHOD)
              {
                DesignElement* ptr = &(data.Last());
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
          }
        }
      }
    }

    for(unsigned int rd=0; rd < nd * nr; rd++){
      if(region_design[rd]){
        throw Exception("not all designs given for all regions");
      }
    }

    LOG_DBG(designSpace) << "data size: " << elements << "*" << pn_design.GetSize() << "=" << data.GetSize();


    // now read the transfer functions
    if(method != ErsatzMaterial::PARAM_MAT && method != ErsatzMaterial::SHAPE_PARAM_MAT)
    {
      if(trans_in.GetSize() == 0)
        throw Exception("no transferFunctions given");

      transfer.Reserve(trans_in.GetSize());
      if(trans_in.GetSize() < nd)
        throw Exception("less transferFunctions than design variable types is infeasible");

      for(unsigned int i = 0; i < trans_in.GetSize(); i++)
        transfer.Push_back(TransferFunction(trans_in[i], design.GetSize() == 1 ? design[0] : DesignElement::NO_TYPE));
    }

  } // no design elements given

  // set the result descriptions which identify the solution types
  ParamNodeList result = pn->GetList("result");
  resultDescriptions.Reserve(result.GetSize());
  for(unsigned int i = 0; i < result.GetSize(); i++)
    resultDescriptions.Push_back(ResultDescription(result[i]));
  
  if(adapt_lower)
  {
    // get the param from the transfer function
    const double par(1.0/GetTransferFunction(DesignElement::DENSITY, Optimization::MECH)->GetParam());
    
    // adapt lower bound of density
    for(unsigned int e = 0, n = data.GetSize(); e < n; ++e)
    {
      DesignElement &curr_de = data[e];
      curr_de.SetLowerBound(std::pow(curr_de.GetLowerBound(), par));
    }
  }

  // reserve for the worst case. non_design_vicinity and off-design optimization
  pseudoDesigns_.Reserve(domain->GetGrid()->GetNumRegions() * design.GetSize());

  if(non_design_vicinity_)
  {
    for(unsigned int d = 0; d < design.GetSize(); d++)
      for(unsigned int r = 0; r < domain->GetGrid()->regionData.GetSize(); r++)
        if(!Contains(domain->GetGrid()->regionData[r].id))
          RegisterPseudoDesignRegion(domain->GetGrid()->regionData[r].id, design[d]);
  }
}

DesignSpace::~DesignSpace(){
  if(designMaterial != NULL){
    delete designMaterial;
    designMaterial = NULL;
  }
}

DesignSpace* DesignSpace::CreateInstance(StdVector<RegionIdType> reg_data, PtrParamNode pn, ErsatzMaterial::Method method){
  switch(method){
  case ErsatzMaterial::SHAPE_OPT:
  case ErsatzMaterial::SHAPE_PARAM_MAT:
    return new ShapeDesign(reg_data, pn, method);
  default:
    return new DesignSpace(reg_data, pn, method);
  }
}

void DesignSpace::PostInit(int objectives, int constraints)
{
  LOG_DBG(designSpace) << "# objectives = " << objectives << ", # constraints = " << constraints;
  DesignElement::SetDesignSpace(this);

  for(unsigned int i = 0, n = totalElements_.GetSize(); i < n; i++)
    totalElements_[i]->PostInit(objectives, constraints);
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
    StdVector<Elem*> elems;
    domain->GetGrid()->GetElems(elems, region);

    StdVector<DesignElement> tmp;
    // the Push_back must not resize!!!!
    assert(pseudoDesigns_.Capacity() >= pseudoDesigns_.GetSize() + 1);
    pseudoDesigns_.Push_back(tmp);

    StdVector<DesignElement>& vec = pseudoDesigns_.Last();

    // construct pesudo design elements
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
      elemToDesign[de->elem->elemNum].second = false; // we have a pesudo design!
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

void DesignSpace::SetDesignMaterial(PtrParamNode dm){
  if(transfer.GetSize() > 0)
    throw Exception("designmaterial can not be given when using transferFunctions");
  transfer.Push_back(TransferFunction()); // create an identity transfer function
  DesignMaterial::SetEnums();
  designMaterial = new DesignMaterial(dm, design);  
}

void DesignSpace::AppendOptimizationResults(SinglePDE* pde)
{
  // set the result descriptions which identify the solution types
  for(unsigned int i = 0; i < resultDescriptions.GetSize(); i++)
  {
	  ResultDescription& rd = resultDescriptions[i];
    ResultInfo* ri = GetResultInfo(rd);
    shared_ptr<ResultInfo> opt_res(ri);
    bool added = pde->CheckStoreResult(opt_res);
    if(!added)
    {
      std::ostringstream os;
      os << "'optimization defines' '" << SolutionTypeEnum.ToString(rd.solutionType) << "' but no PDE references it in it's 'storeResults'";

      WARN(os.str().c_str());
    }
  }
}

double DesignSpace::GetNodalValue(unsigned int nodeNumber, DesignElement::ValueSpecifier vs)
{
  ShapeOptimizer* shopt = dynamic_cast<ShapeOptimizer*>(optimizer_);
  if(shopt == NULL) EXCEPTION("No level set optimizer activated");
  // FIXME maybe throw an Exception? This should not be called without a levelset
  if(shopt->ptrLS_ == NULL) return 0.0;
  assert(shopt->ptrLS_->GetNodePointer(nodeNumber) != NULL);
  
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
  default:
    EXCEPTION("case not implemented")
  }
}

ResultInfo* DesignSpace::GetResultInfo(ResultDescription& rd)
{
  // <result id="optResult_1" design="density" access="plain" value="costGradient"/>

  ResultInfo* ri = new ResultInfo();
  // I hate it!!! :(
  ri->resultType = (SolutionType) rd.solutionType;

  // no space and brackets to have no problems with info.xml and no problems with the paraview calculator
  ri->resultName = DesignElement::valueSpecifier.ToString(rd.value) + "_"
                   + (rd.detail != DesignElement::NONE ? (DesignElement::detail.ToString(rd.detail) + "_") : "")
                   + DesignElement::type.ToString(rd.design) + "_"
                   + DesignElement::access.ToString(rd.access)
                   + (rd.excitation != "" ? ("_ex_" + rd.excitation) : "");

  ri->unit = "";

  ri->entryType = ResultInfo::SCALAR;
  ri->dofNames = "";

  // in most cases we are on elements,
  switch(rd.value)
  {
  case DesignElement::LEVEL_SET_VALUE:
  case DesignElement::LEVEL_SET_STATE:
  case DesignElement::SHAPEGRAD_NODE_VALUE:
  case DesignElement::LEVEL_SET_GRAD_XP:
  case DesignElement::LEVEL_SET_GRAD_XN:
  case DesignElement::LEVEL_SET_GRAD_YP:
  case DesignElement::LEVEL_SET_GRAD_YN:
  case DesignElement::LEVEL_SET_GRAD_ZP:
  case DesignElement::LEVEL_SET_GRAD_ZN:
    ri->definedOn = ResultInfo::NODE;
    break;
  default:
    ri->definedOn = ResultInfo::ELEMENT;
  }
  ri->fctType = shared_ptr<ConstFct>(new ConstFct() );

  // let the caller or a shared pointer delete it finally
  return ri;
}

int DesignSpace::GetSpecialResultIndex(DesignElement::Type design, DesignElement::ValueSpecifier value,
                                       DesignElement::Detail detail, DesignElement::Access access, const std::string& excitation)
{
  for(unsigned int i = 0; i < resultDescriptions.GetSize(); i++)
  {
    const ResultDescription& rd = resultDescriptions[i];
    // two step check
    if(rd.design != design || rd.value != value || rd.detail != detail || rd.access != access) continue;
    // second check
    if(rd.excitation != "" && rd.excitation != excitation) continue;

    // we are right.
    switch(rd.solutionType)
    {
      case OPT_RESULT_1: return 0;
      case OPT_RESULT_2: return 1;
      case OPT_RESULT_3: return 2;
      case OPT_RESULT_4: return 3;
      case OPT_RESULT_5: return 4;
      case OPT_RESULT_6: return 5;
      case OPT_RESULT_7: return 6;
      case OPT_RESULT_8: return 7;
      case OPT_RESULT_9: return 8;
      default: throw Exception("invalid solution type");
    }
  }

  return -1; // the specified triple was not specified such in xml
}

void DesignSpace::AssertOneDesignOnly()
{
  if(design.GetSize() != 1)
    throw Exception("A feature relies on a single design only!");
}


int DesignSpace::FindDesign(DesignElement::Type dt, bool throw_exception)
{
  // do a fallback for NO_TYPE and DEFAULT
  if(design.GetSize() == 1 && (dt == DesignElement::NO_TYPE || dt == DesignElement::DEFAULT))
    return 0;

  if(dt == DesignElement::TENSOR_TRACE && HasErsatzMaterialTensor()) // this is not a real type of design, but volume constraint can operate on it, if optimization returns a complete tensor
    return 0;

  // search where in data we are
  int base = -1;
  for(unsigned int i = 0; i < design.GetSize(); i++)
    if(design[i] == dt) base = i;

  if(base == -1 && throw_exception)
    EXCEPTION("Design " << DesignElement::type.ToString(dt) << " not within " << design.GetSize() << " actual designs.");

  return base;
}


double DesignSpace::GetErsatzMaterialFactor(unsigned int design_index, Optimization::Application applic)
{
  // now do the trick, that the piezo coupling factor might be a product of the
  // density transfer function and the polarization transfer function
  double result = 1.0;

  // go over all design elements we have (one for design only, with polarization
  // it is two
  for(unsigned int index = design_index; index < data.GetSize(); index += elements)
  {
    // note that this loop with loop normally once or twice (piezo)
    DesignElement* de = &data[index];
    // The design of the current element
    DesignElement::Type dt = de->GetType();

    // find the transfer function for our form and application.
    // There is not necessary a transfer function -> e.g. polarization
    // is for the piezo only defined on the coupling
    TransferFunction* tf = GetTransferFunction(dt, applic, false);
    LOG_DBG3(designSpace) << "GEMF: dt=" << DesignElement::type.ToString(dt) << " app=" << Optimization::application.ToString(applic) << " tf found=" << (tf != NULL);
    // multiply our transfer function
    if(tf != NULL) {
      double transformed = tf->Transform(de, DesignElement::SMART); // handles design filtering
      LOG_DBG3(designSpace) << "GEMF: ErsatzMaterial for " << de->elem->elemNum << "/"
                       << Optimization::application.ToString(applic) << " for "
                       << DesignElement::type.ToString(dt) << ": "
                       << TransferFunction::type.ToString(tf->GetType()) << "("
                       << de->GetDesign(DesignElement::PLAIN) << ") = " << transformed
                       << " -> * " << result << " = " << (result * transformed);
      result *= transformed;

    }
  }

  return result;
}

bool DesignSpace::GetErsatzMaterialPamping(const Elem* elem, Matrix<double>& elemMat)
{
  // see also implementation SIMP::AddMassToStiffness() for match!!!
  static MechMat mm = MechMat(this); // Assumes irregular mesh :(

  // pamping at all -> see Sigmund; Morphology; 2007
  assert(GetPampingValue() >= 0);
  // have design?
  DesignElement* de = Find(elem->elemNum, DesignElement::DENSITY, false);
  if(de == NULL)
    return false;

  // we use the physical design variable to match better
  TransferFunction* tf = GetTransferFunction(de->GetType(), Optimization::MASS);
  double tv = tf->Transform(de, DesignElement::SMART); // be consistent with SIMP::AddMassToStiffness()
  // now the original mass matrix
  const Matrix<double>& mass = mm.Mass(de->elem);
  LOG_DBG3(designSpace) << "GEMP e=" << elem->elemNum << " mass=" << mass.ToString();
  elemMat.Resize(mass.GetNumRows(), mass.GetNumCols());
  elemMat.Assign(mass, tv * (1.0-tv) * GetPampingValue());

  LOG_DBG3(designSpace) << "GEMP e=" << elem->elemNum << "rv=" << tv << " p=" << GetPampingValue() << " -> " << (tv * (1.0-tv) * GetPampingValue());
  LOG_DBG3(designSpace) << "GEMP e=" << elem->elemNum << " ->" << elemMat.ToString();
  return true;
}


bool DesignSpace::CollectMaterialParametersForElement(const Elem* elem){
  int base = Find(elem, false);
  if(base < 0){
    return(false);
  }
  
  for(unsigned int index = base; index < data.GetSize(); index += elements){
    DesignElement* de = &data[index];
    designMaterial->SetParameter(de->GetType(), de->GetDesign(DesignElement::PLAIN));
  }
  return(true);
}

bool DesignSpace::GetErsatzMaterialTensor(Matrix<double>& t, SubTensorType subTensor, const Elem* elem, DesignElement::Type direction){
  // collect all parameters
  if(CollectMaterialParametersForElement(elem)){
    designMaterial->GetMaterialTensor(t, subTensor, direction);
    return(true);
  }
  return(false);
}

double DesignSpace::GetErsatzMaterialMass(const Elem* elem, DesignElement::Type direction){
  // collect all parameters
  if(CollectMaterialParametersForElement(elem)){
    return(designMaterial->GetMaterialMass(direction));
  }
  return(1.0);
}

bool DesignSpace::GetErsatzMaterialDamping(double& alpha, double& beta, const Elem* elem, DesignElement::Type direction){
  if(CollectMaterialParametersForElement(elem)){
    return(designMaterial->GetMaterialDamping(alpha, beta, direction));
  }
  return(false);
}

bool DesignSpace::GetErsatzMaterialDampingParameterForIntegrator(const Elem* elem, BaseForm* form, double& param){
  if(CollectMaterialParametersForElement(elem)){
    double dummy = 0.0;
    if(form->GetName() == "MassInt") return(designMaterial->GetMaterialDamping(param, dummy));
    if(form->GetName() == "linElastInt") return(designMaterial->GetMaterialDamping(dummy, param));
  }
  return(false);
}

TransferFunction* DesignSpace::GetTransferFunction(DesignElement* de)
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


TransferFunction* DesignSpace::GetTransferFunction(DesignElement::Type design, Optimization::Application application, bool throw_exception)
{
  if(HasErsatzMaterialTensor()){
    return &transfer[0]; // this will always point to an identity transfer function, so CalcU1KU2 in ErsatzMaterial will simply work for parametric material optimization
  }

  for(unsigned int i = 0; i < transfer.GetSize(); i++)
  {
    TransferFunction* tf = &transfer[i];
    if(tf->GetApplication() != application)
      continue;
    if(tf->GetDesign() == design)
      return tf;
    if(this->design.GetSize() == 1 && design == DesignElement::DEFAULT)
      return tf;
  }

  if(!throw_exception) return NULL;
  else throw Exception("the desired transfer function for design '" + DesignElement::type.ToString(design)
                        + "' in application '" + Optimization::application.ToString(application)
                        + "' is not contained");
}


int DesignSpace::ReadDesignFromExtern(const double* space)
{
  bool new_design = false;
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
          if(!new_design && data[d].GetDesign(DesignElement::PLAIN) != v) {
            new_design = true;
          }
          LOG_DBG3(designSpace) << "ReadDesignFromExtern: setting design: data[" << d << "] = " << v << " = in[" << s << "]";
          data[d].SetDesign(v);
          s++; // advance in every step
        } // for d
      }else if(cur_reg.constant == CONSTANT_PER_REGION || cur_reg.constant == CONSTANT_ON_ALL_REGIONS){ // in FIXED case, nothing is done
        const double v = space[s] * scaling + translation;
        for(unsigned int d = cur_reg.base; d < u; d++){
          if(!new_design && data[d].GetDesign(DesignElement::PLAIN) != v) {
            new_design = true;
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
  if(new_design) design_id++;
  return design_id;
}

int DesignSpace::ReadDesignFromExtern(const StdVector<double>& space)
{
  return ReadDesignFromExtern(space.GetPointer());
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

int DesignSpace::WriteDesignToExtern(StdVector<double>& space_out, bool scaling) const
{
  space_out.Reserve(GetNumberOfVariables());
  return WriteDesignToExtern(space_out.GetPointer(), scaling);
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

void DesignSpace::WriteSparseGradientToExtern(StdVector<double>& out, DesignElement::ValueSpecifier vs, DesignElement::Access access, Condition* g, bool use_scaling) const
{
  // Bastian did some complicated reordering stuff. For the only case of sparse Jacobians (slope constraints)
  // we'll have only the simple standard situation .. if this changes you have at least a test case :) Fabian
  
  assert(regions.GetSize() == 1 && regions[0].GetSize() == 1); // only one region with one design
  assert(g != NULL); // only constraints can have sparse Jacobians 

  const double scaling = use_scaling ? regions[0][0].scale_design : 1.0;
  
  StdVector<unsigned int>& sparsity = g->GetSparsityPattern();
  
  assert(out.window.GetSize() == sparsity.GetSize());
  unsigned int base = out.window.GetStart();
  for(unsigned int i = 0; i < sparsity.GetSize(); i++)
  {
    assert(out.InWindow(base + i));
    out[base + i] = data[sparsity[i]].GetValue(vs, access, g) * scaling;
  }
}

void DesignSpace::WriteDenseGradientToExtern(StdVector<double>& out, DesignElement::ValueSpecifier vs, DesignElement::Access access, Condition* g, bool use_scaling) const
{
  // this does now do reordering as gradients are reordered in the optimizer
  // must be set in the constructor! might be trivial volume fraction or from file!!
  assert(!(vs == DesignElement::COST_GRADIENT && g != NULL));

  unsigned int n0 = out.window.GetStart(); // to grow up to the total number of design variables
  unsigned int n = n0;

  const unsigned int nd = design.GetSize();

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
          assert(out.InWindow(n));
          LOG_DBG3(designSpace) << "WriteDenseGradientToExtern: non-constant region " << r << ": out[" << n << "] = design[" << s << "]=" << data[s].GetValue(vs, access, g);
          out[n++] = data[s].GetValue(vs, access, g) * scaling;
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
          out[n] += data[s].GetValue(vs, access, g) * scaling;
          LOG_DBG3(designSpace) << "WriteDenseGradientToExtern: constant region " << r << ": out[" << n << "] += design[" << s << "]=" << data[s].GetValue(vs, access, g);
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
  unsigned int start = design == DesignElement::DEFAULT || DesignElement::TENSOR_TRACE ? 0 : FindDesign(design) * elements;
  unsigned int end   = design == DesignElement::DEFAULT || DesignElement::TENSOR_TRACE ? data.GetSize() : start + elements;

  LOG_DBG3(designSpace) << "Reset: vs=" << DesignElement::valueSpecifier.ToString(vs) << " design="
                        << DesignElement::type.ToString(design) << " from " << start << " to " << end;

  for(unsigned int i = start; i < end; i++)
  {
    DesignElement& de = data[i];

    switch(vs)
    {
      case DesignElement::DESIGN:
           de.SetDesign(0.0);
           break;

      case DesignElement::CONSTRAINT_GRADIENT:
      case DesignElement::COST_GRADIENT:
           de.Reset(vs);
           break;

      default: throw Exception("value specifier not handled");
    }
  }
}
DesignElement* DesignSpace::Find(unsigned int elemNum, DesignElement::Type dt, bool throw_exception, bool include_pseudo_designs)
{
  int idx = Find(elemNum, throw_exception, include_pseudo_designs);

  if(idx == -1) return NULL; // an exception was already thrown if desired

  // check for real design or pseudo design
  if(elemToDesign[elemNum].second == true)
  {
    for(unsigned int d = 0, nd = design.GetSize(); d < nd; d++)
    {
      unsigned int pos = elements * d + idx;
      if(data[pos].GetType() == dt)
      {
        if(data[pos].elem->elemNum != elemNum) throw Exception("index mixed up");
        return &data[pos];
      }
    }
  }
  else
  {
    for(unsigned int i = 0, n = pseudoDesigns_.GetSize(); i < n; i++)
    {
       assert(pseudoDesigns_[i][idx].elem->elemNum == elemNum);
       if(pseudoDesigns_[i][idx].GetType() == dt)
         return &(pseudoDesigns_[i][idx]);
    }
  }

  if(throw_exception) throw Exception("design type not in design");
  return NULL;
}


inline
int DesignSpace::Find(unsigned int elemNum, bool throw_exception, bool include_pseudo_designs)
{
  int idx = elemToDesign[elemNum].first;

  // reset pseudo designs when we don't look for them explicitly
  if(idx != -1 && !include_pseudo_designs && elemToDesign[elemNum].second == false)
    idx = -1;

  if(idx == -1 && throw_exception)
    EXCEPTION("could not find element " << elemNum << " in our (pseudo) design space");

  return idx;
}

int DesignSpace::Find(const Elem* elem, bool throw_exception)
{
  // no extensions for pseudo designs implemented, yet!
  if(FindRegion(elem->regionId) >= 0) return Find(elem->elemNum, throw_exception);

  // we might have surface element and it is pointing to a design element
  const SurfElem* se = dynamic_cast<const SurfElem*>(elem);

  // no chance, we are wrong
  if(se == NULL)
  {
    if(!throw_exception) return -1;
    EXCEPTION("element " << elem->ToString() << " not in design regions" );
  }

  if(se->ptVolElem1 != NULL && FindRegion(se->ptVolElem1->regionId) >= 0)
    return Find(se->ptVolElem1->elemNum);

  if(se->ptVolElem2 != NULL && FindRegion(se->ptVolElem2->regionId) >= 0)
    return Find(se->ptVolElem2->elemNum);

  if(!throw_exception) return -1;

  EXCEPTION("element " << elem->ToString() << " has no volume element in design region");
}

DesignElement* DesignSpace::FindElementWithLargesFilter()
{
  int max = 0;
  DesignElement* res = NULL;

  for(unsigned int i = 0, n = data.GetSize(); i < n; i++)
  {
    DesignElement* de = &data[i];
    if(de->simp != NULL && (int) de->simp->neighborhood.GetSize() > max)
      res = de;
  }

  return res;
}

void DesignSpace::ToInfo(PtrParamNode in)
{
  PtrParamNode tf = in->Get("transferFunctions");
  for(unsigned int i = 0; i < transfer.GetSize(); i++)
    transfer[i].ToInfo(tf->Get("transferFunction", ParamNode::APPEND));

  PtrParamNode dv = in->Get("designVariables");
  dv->Get("optimizationVariables")->SetValue(GetNumberOfVariables());
  // TODO @Bastian - add you shape stuff if you like
  dv->Get("modelVariables")->SetValue(data.GetSize());
  dv->Get("modelElements")->SetValue(GetNumberOfElements());
  for(unsigned int i = 0; i < design.GetSize(); i++)
    data[i * elements].ToInfo(dv->Get("design", ParamNode::APPEND));

  in->Get("pamping")->SetValue(pamping_);
  if(regions.GetSize() > 0)
  {
    PtrParamNode rs = in->Get("regions");
    for(unsigned int i = 0; i < regions[0].GetSize(); i++)
      regions[0][i].ToInfo(rs->Get("region", ParamNode::APPEND));
  }

}


std::string DesignSpace::ToString()
{
  std::stringstream ss;
  ss << "design[";
  const unsigned int data_size(data.GetSize());
  for(unsigned int i = 0; i < data_size; i++)
  {
    DesignElement* de = &data[i];
    ss << i << ":elem=" << de->elem->elemNum;
    if(de->vicinity != NULL) ss << " " << de->vicinity->ToString();
    ss << " ";
  }
  ss << "]";
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
  const unsigned int nr = regions[0].GetSize();
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

int DesignSpace::FindRegion(RegionIdType regionId){
  if(regions.GetSize() > 0){
    const StdVector<DesignRegion>& regs = regions[0];
    const unsigned int nr = regs.GetSize();  
    for(unsigned int r = 0; r < nr; r++){
      if(regs[r].regionId == regionId){
        return r;
      }
    }
  }
  return -1;
}


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
  // this is clearly nonsense if the result/solution type is OPT_RESULT_*
  switch(ri->resultType)
  {
  case MECH_PSEUDO_DENSITY:
  case PHYSICAL_PSEUDO_DENSITY:
    def.design = DesignElement::DENSITY;
    break;
  case ELEC_PSEUDO_POLARIZATION:
    def.design = DesignElement::POLARIZATION;
    break;
  case ACOU_PSEUDO_DENSITY:
    def.design = DesignElement::ACOU_DENSITY;
    break;
  default:
    // to be overwritten by the ResultDescription
    def.design = DesignElement::DENSITY;
  }
  // somehow critical! but only for density filtering, if at all.
  def.access = ri->resultType == PHYSICAL_PSEUDO_DENSITY ? DesignElement::SMART : DesignElement::PLAIN;
  def.value  = DesignElement::DESIGN;

  ResultDescription& descr = def;
  // ignore defaults if there is a result description for the OPT_RESULT_* case
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

  // this is our entity result, a scalar or a vector of dim 2/3
  unsigned int dofs = result.GetResultInfo()->dofNames.GetSize();
  assert(dofs >= 1 && dofs <= 3);
  StdVector<double> result_value(dofs);


  // search where in data we are
  int base = FindDesign(descr.design);

  // loop over elements from result. We have to do it this way as the the connection
  // of design element and result element is the element(->elemeNum) but we cannot
  // search in the result for an element.
  EntityIterator it = result.GetEntityList()->GetIterator();

  // set the result as we need it
  result_data.Resize(result.GetEntityList()->GetSize() * dofs);

  // the default value is 0.0 but 1 for densities
  SolutionType st = result.GetResultInfo()->resultType;
  double none = st == MECH_PSEUDO_DENSITY || st == PHYSICAL_PSEUDO_DENSITY || st == ELEC_PSEUDO_POLARIZATION ? 1.0 : 0.0; 
  
  for ( it.Begin(); !it.IsEnd(); it++ )
  {
    // for elements not in the design region we set to to the default value
    for(unsigned int i = 0; i < dofs; i++)
      result_value[i] = none;

    if(FindRegion(it.GetElem()->regionId) >= 0)
    {
      // note that the index is from the first design set!
      unsigned int base_index = Find(it.GetElem()->elemNum, true, true); // exception and pseudo designs (?)
      // base=0 is first!
      unsigned int data_index = (base * elements) + base_index;
      DesignElement& de = data[data_index];
      de.GetValue(descr, result_value, dofs);

      #ifdef CHECK_INDEX
        if(de.elem->elemNum != it.GetElem()->elemNum)
          EXCEPTION("mixed up indices:" << de.elem->elemNum << "!=" << it.GetElem()->elemNum
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

BaseMaterial* DesignSpace::GetBiMaterial(RegionIdType reg, Optimization::Application app, bool throw_exception)
{
  DesignRegion* dr = GetRegion(reg, false); // tolerant for off-design stress constraints
  if(dr == NULL || !dr->HasBiMaterial())
  {
    if(!throw_exception)
      return NULL;
    else
      EXCEPTION("cannot find bimaterial for region" << reg << " and application " << app);
  }

  switch(app)
  {
  case Optimization::MECH:
    return dr->GetBiMaterial(MECHANIC);
  case Optimization::ELEC:
    return dr->GetBiMaterial(ELECTROSTATIC);
  case Optimization::PIEZO_COUPLING:
    return dr->GetBiMaterial(PIEZO);
  default:
    assert(false); // implement what you need!
    return NULL;
  }
}

DesignSpace::DesignRegion* DesignSpace::GetRegion(RegionIdType id, bool throw_exception)
{
  for(unsigned int i = 0, n = regions[0].GetSize(); i < n; i++)
    if(regions[0][i].regionId == id)
      return &regions[0][i];
  if(!throw_exception)
    return NULL;
  else
    EXCEPTION("invalid region id");
}

DesignSpace::DesignRegion::DesignRegion()
{
  regionId = -1;
}

bool DesignSpace::DesignRegion::HasBiMaterial() const
{
  return bimaterial_ != "";
}


BaseMaterial* DesignSpace::DesignRegion::GetBiMaterial(const MaterialClass mc)
{
  assert(bimaterial_ != ""); // check with HasBiMaterial()!

  for(unsigned int i = 0; i < materials_.GetSize(); i++)
    if(materials_[i].second == mc)
      return materials_[i].first;

  // apparently first run
  MaterialHandler* matLoader = domain->GetMaterialHandler();
  BaseMaterial* mat = matLoader->LoadMaterial(bimaterial_, mc);
  materials_.Push_back(std::make_pair(mat, mc));

  return mat;
}

void DesignSpace::DesignRegion::ToInfo(PtrParamNode node) const
{
  node->Get("name")->SetValue(domain->GetGrid()->GetRegion().ToString(regionId));
  node->Get("elements")->SetValue(elements);
  node->Get("bimaterial")->SetValue(HasBiMaterial() ? bimaterial_ : "-");
}

// explicit template instantiation for GCC compiler
#ifdef __GNUC__


template
void DesignSpace::ExtractResults<double>(shared_ptr<BaseResult> base_result);

template
void DesignSpace::ExtractResults<complex<double> >(shared_ptr<BaseResult> base_result);

template
void DesignSpace::FillNodeResults<double>(Result<double>& result, ResultDescription& descr);

template
void DesignSpace::FillNodeResults<complex<double> >(Result<complex<double> >& result, ResultDescription& descr);

template
void DesignSpace::FillElementResults<double>(Result<double>& result, ResultDescription& descr);

template
void DesignSpace::FillElementResults<complex<double> >(Result<complex<double> >& result, ResultDescription& descr);

#endif
