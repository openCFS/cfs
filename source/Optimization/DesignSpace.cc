#include "Optimization/DesignSpace.hh"
#include "Optimization/DesignElement.hh"
#include "Optimization/TransferFunction.hh"
#include "Optimization/Condition.hh"
#include "Optimization/BaseOptimizer.hh"
#include "Optimization/LevelSet.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/InfoNode.hh"
#include "General/exception.hh"
#include "General/Enum.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"
#include "Domain/resultInfo.hh"
#include "Utils/result.hh"
#include "Utils/StdVector.hh"
#include "PDE/SinglePDE.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "Optimization/ShapeDesign.hh"

using namespace CoupledField;

// declare class specific logging stream
DECLARE_LOG(designSpace)
DEFINE_LOG(designSpace, "designSpace")

// declare class specific logging stream
DECLARE_LOG(ersatz)
DEFINE_LOG(ersatz, "ersatzMaterialFactor")

DesignSpace::DesignSpace(StdVector<RegionIdType>& regionIds, StdVector<ParamNode*>& pn_design, StdVector<ParamNode*>& trans_in, StdVector<ParamNode*>& result, ErsatzMaterial::Method method)
{
  LOG_TRACE(designSpace) << "DesignSpace for regions=" << regionIds << " #designs=" << pn_design.GetSize()
                         << " #transferFunctions=" << trans_in.GetSize() << " #results=" << result.GetSize();
  last_find_index_ = 0;

  regions_.Resize(regionIds.GetSize());

  design_id = 0;
  optimizer_ = NULL;

  designMaterial = NULL;

  // only temporary
  StdVector<Elem*> elems;

  applicationForm.SetName("DesignSpace::ApplicationForm");
  applicationForm.Add(Optimization::ELEC, "linElecInt");
  applicationForm.Add(Optimization::MECH, "linElastInt");
  // We follow for the stress, strain calculation the transfer functions of mech
  applicationForm.Add(Optimization::MECH, "MechStressStrain", false);
  applicationForm.Add(Optimization::PIEZO_COUPLING, "linPiezoCoupling");
  applicationForm.Add(Optimization::CHARGE_DENSITY, "LinNeumannInt");
  applicationForm.Add(Optimization::PRESSURE, "PressureLinForm");
  applicationForm.Add(Optimization::MASS, "MassInt");

  // read the elements
  elements_ = domain->GetGrid()->GetNumElems(regionIds);

  // number of different designs
  unsigned int nd = 0;
  for(unsigned int d = 0; d < pn_design.GetSize(); d++){
    DesignElement::Type dt = DesignElement::type.Parse(pn_design[d]->Get("name"));
    if(design.Find(dt) < 0){
      design.Push_back(dt);
      nd++;
    }
  }
  
  if(elements_ == 0 || nd == 0){ // this may happen in shape optimization

    if(method != ErsatzMaterial::SHAPE_OPT && method != ErsatzMaterial::SHAPE_PARAM_MAT){
      if(elements_ == 0) throw Exception("empty regions");
      if(nd == 0) throw Exception("no designs given");
    }

  }else{

  // set our own structure with is element times design parameters
  data.Reserve(elements_ * nd);

  // check whether all regions have all design variables and all design variables are given in all regions
  StdVector<bool> region_design;
  region_design.Resize(nd * regionIds.GetSize(), true);
  
  // scaling
  scale_design.Resize(nd);
  translate_design.Resize(nd);
  for(unsigned int d = 0; d < nd; d++){
    scale_design[d].Resize(regionIds.GetSize(), 1.0);
    translate_design[d].Resize(regionIds.GetSize(), 0.0);
  }

  StdVector<std::string> regNames;
  domain->GetGrid()->GetRegionNames(regNames);
  
  // iterate over all designs, we have to always have the same order of designs not the order from the xml file
  for(unsigned int dti = 0; dti < design.GetSize(); dti++){
    DesignElement::Type dt = design[dti];
    for(unsigned int d = 0; d < pn_design.GetSize(); d++) {
      if(DesignElement::type.Parse(pn_design[d]->Get("name")) != dt){
        continue;
      }
      std::string design_reg = pn_design[d]->Get("region")->AsString();

      for(unsigned int r = 0; r < regionIds.GetSize(); r++){
        domain->GetGrid()->GetElems(elems, regionIds[r]);
        unsigned int n = elems.GetSize();
        std::string reg = regNames[regionIds[r]];

        if(d == 0){
          regions_[r].regionId = regionIds[r];
          regions_[r].base = r == 0 ? 0 : regions_[r-1].base + regions_[r-1].elements;
          regions_[r].elements = n;
          regions_[r].constant = false;
        }

        if(design_reg == "all" || design_reg == reg){
          int di = design.Find(dt);
          if(!region_design[r*nd + di]){
            throw Exception("Design/Region combination given twice!");
          }
          region_design[r*nd + di] = false;
          if(pn_design[d]->Get("constant")->AsBool()){ // we have a constant densign-value on that region
            regions_[r].constant = true;
          }
          if(pn_design[d]->Get("scale")->AsBool()){
            double upper = pn_design[d]->Get("upper")->AsDouble();
            double lower = pn_design[d]->Get("lower")->AsDouble();
            scale_design[di][r] = (upper - lower);
            translate_design[di][r] = lower;
          }

          double initial = pn_design[d]->Get("initial")->AsDouble();
          LOG_DBG2(designSpace) << "add design " << dt << ":" << DesignElement::type.ToString(dt)
          << " initial=" << initial;

          for(unsigned int e = 0; e < n; e++)
          {
            DesignElement de(pn_design[d], elems[e]);
            de.SetDesign(initial);
            data.Push_back(de);
            switch(method)
            {
            case ErsatzMaterial::SIMP_METHOD:
              data.Last().simp = new SIMPElement(&data.Last());
              break;
            case ErsatzMaterial::PARAM_MAT:
            case ErsatzMaterial::NO_METHOD:
            case ErsatzMaterial::SHAPE_PARAM_MAT:
              break;
            default: assert(false);
            }
          }
        }
      }
    }
  }

  for(unsigned int rd=0; rd < nd * regionIds.GetSize(); rd++){
    if(region_design[rd]){
      throw Exception("not all designs given for all regions");
    }
  }

  LOG_DBG(designSpace) << "data size: " << elements_ << "*" << pn_design.GetSize() << "=" << data.GetSize();


  // now read the transfer functions
  if(method != ErsatzMaterial::PARAM_MAT && method != ErsatzMaterial::SHAPE_PARAM_MAT)
  {
    if(trans_in.GetSize() == 0)
      throw Exception("no transferFunctions given");

    transfer.Reserve(trans_in.GetSize());
    if(trans_in.GetSize() < nd)
      throw Exception("less transferFunctions than design variable types is inveasible");

    for(unsigned int i = 0; i < trans_in.GetSize(); i++)
      transfer.Push_back(TransferFunction(trans_in[i], design.GetSize() == 1 ? design[0] : DesignElement::NO_TYPE));
  }
  
  } // no designelements given

  // set the result descriptions which identify the solution types
  resultDescriptions.Reserve(result.GetSize());
  for(unsigned int i = 0; i < result.GetSize(); i++)
    resultDescriptions.Push_back(ResultDescription(result[i]));
}

DesignSpace::~DesignSpace(){
  if(designMaterial != NULL){
    delete designMaterial;
    designMaterial = NULL;
  }
}

DesignSpace* DesignSpace::CreateInstance(StdVector<RegionIdType> regionIds, StdVector<ParamNode*>& design, StdVector<ParamNode*>& transfer, StdVector<ParamNode*>& result, ErsatzMaterial::Method method){
  switch(method){
  case ErsatzMaterial::SHAPE_OPT:
  case ErsatzMaterial::SHAPE_PARAM_MAT:
    return new ShapeDesign(regionIds, design, transfer, result, method);
  default:
    return new DesignSpace(regionIds, design, transfer, result, method);
  }
}

void DesignSpace::PostInit(int constraints)
{
  for(unsigned int i = 0; i < data.GetSize(); i++)
    data[i].PostInit(constraints);
}

void DesignSpace::SetDesignMaterial(ParamNode* dm){
  if(transfer.GetSize() > 0)
    throw Exception("designmaterial can not be given when using transferFunctions");
  transfer.Push_back(TransferFunction()); // create an identity transfer function
  DesignMaterial::SetEnums();
  designMaterial = new DesignMaterial(dm);
  if(designMaterial->RequiredDesigns() != design.GetSize())
    throw Exception("number of required designs not correct");
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
      os << "<optimization> defines '" << SolutionTypeEnum.ToString(rd.solutionType) << "' but no PDE references it in it's <storeResults>";
      Warning(os.str().c_str());
    }
  }
}


double DesignSpace::GetNodalValue(unsigned int nodeNumber, DesignElement::ValueSpecifier vs)
{
  switch(vs)
  {
  case DesignElement::LEVEL_SET_VALUE:
  {
    LevelSet* ls = dynamic_cast<LevelSet*>(optimizer_);
    if(ls == NULL) throw Exception("No level set optimizer activated");
    return ls->Find(nodeNumber)->GetCenterValue();
  }

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

  ri->resultName = DesignElement::valueSpecifier.ToString(rd.value) + "_"
                   + (rd.detail != DesignElement::NONE ? (DesignElement::detail.ToString(rd.detail) + "_") : "")
                   + DesignElement::type.ToString(rd.design) + " ("
                   + DesignElement::access.ToString(rd.access) + ")";


  ri->unit = "";
  // in most cases we have a scalar result
  switch(rd.value)
  {
  case DesignElement::LEVEL_SET_NORMAL:
    ri->entryType = ResultInfo::VECTOR;
    ri->dofNames.Resize(domain->GetGrid()->GetDim());
    ri->dofNames[0] = "x";
    ri->dofNames[1] = "y";
    if(ri->dofNames.GetSize() == 3) ri->dofNames[2] = "z";
    break;
  default:
    ri->entryType = ResultInfo::SCALAR;
    ri->dofNames = "";
  }

  // in most cases we are on elements,
  switch(rd.value)
  {
  case DesignElement::LEVEL_SET_VALUE:
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
                                       DesignElement::Detail detail, DesignElement::Access access)
{
  for(unsigned int i = 0; i < resultDescriptions.GetSize(); i++)
  {
    const ResultDescription& rd = resultDescriptions[i];
    if(rd.design != design || rd.value != value || rd.detail != detail || rd.access != access) continue;

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

  return -1; // the specified tripel was not specified such in xml
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
  for(unsigned int index = design_index; index < data.GetSize(); index += elements_)
  {
    // note that this loop with loop normaly once or twice (piezo)
    DesignElement* de = &data[index];
    // The design of the current element
    DesignElement::Type dt = de->GetType();

    // find the transfer function for our form and application.
    // There is not necessary a transfer function -> e.g. polarization
    // is for the piezo only defined on the coupling
    TransferFunction* tf = GetTransferFunction(dt, applic, false);
    // multiply our transfer function
    if(tf != NULL) {
      double transformed = tf->Transform(de);
      LOG_DBG3(ersatz) << "ErsatzMaterial for " << de->elem->elemNum << "/"
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

bool DesignSpace::GetErsatzMaterialTensor(Matrix<double>& t, SubTensorType subTensor, const Elem* elem, DesignElement::Type direction){
  int base = Find(elem, false);
  if(base < 0){
    return(false);
  }
  // collect all parameters
  for(unsigned int index = base; index < data.GetSize(); index += elements_){
    DesignElement* de = &data[index];
    designMaterial->SetParameter(de->GetType(), de->GetDesign(DesignElement::PLAIN));
  }
  // get the Material-Tensor
  designMaterial->GetMaterialTensor(t, subTensor, direction);
  return(true);
}

TransferFunction* DesignSpace::GetTransferFunction(DesignElement::Type design, Optimization::Application application, bool throw_exception)
{
  if(HasErsatzMaterialTensor()){
    return &transfer[0]; // this will always point to an identity transfer function, so CalcU1KU2 in ErsatzMaterial will simply work for parametric material optimization
  }

  for(unsigned int i = 0; i < transfer.GetSize(); i++)
  {
    TransferFunction* tf = &transfer[i];
    if(tf->GetDesign() == design && tf->GetApplication() == application)
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
  unsigned int s = 0;
  for(unsigned int des = 0; des < design.GetSize(); des++){
    const unsigned int base = des * elements_;
    for(unsigned int r = 0; r < regions_.GetSize(); r++){
      const double scaling = scale_design[des][r];
      const double translation = translate_design[des][r];
      const unsigned int u = base + regions_[r].base + regions_[r].elements;
      if(regions_[r].constant){
        const double v = space[s] * scaling + translation;
        for(unsigned int d = base + regions_[r].base; d < u; d++){
          if(!new_design && data[d].GetDesign(DesignElement::PLAIN) != v)
            new_design = true;
          data[d].SetDesign(v);
        } // for d
        s++; // only advance after having set all element of this region to the corresponding value
      }else{
        for(unsigned int d = base + regions_[r].base; d < u; d++){
          double v = space[s] * scaling + translation;
          if(!new_design && data[d].GetDesign(DesignElement::PLAIN) != v)
            new_design = true;
          data[d].SetDesign(v);
          s++; // advance in every step
        } // for d
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

int DesignSpace::WriteDesignToExtern(double* space, bool scaling) const
{
  unsigned int d = 0;
  for(unsigned int des = 0; des < design.GetSize(); des++){
    const unsigned int base = des * elements_;
    for(unsigned int r = 0; r < regions_.GetSize(); r++){
      const double rscaling = scaling ? 1.0 / scale_design[des][r] : 1.0;
      const double translation = scaling ? translate_design[des][r] : 0.0;
      if(regions_[r].constant){
        space[d++] = (data[base + regions_[r].base].GetDesign(DesignElement::PLAIN) - translation) * rscaling;
      }else{
        const unsigned int u = base + regions_[r].base + regions_[r].elements;
        for(unsigned int s = base + regions_[r].base; s < u; s++){
          space[d++] = (data[s].GetDesign(DesignElement::PLAIN) - translation) * rscaling;
        }
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
  unsigned int d = 0;
  for(unsigned int des = 0; des < design.GetSize(); des++){
    const unsigned int base = des * elements_;
    for(unsigned int r = 0; r < regions_.GetSize(); r++){
      const double rscaling = 1.0 / scale_design[des][r];
      const double translation = translate_design[des][r];
      if(regions_[r].constant){
        x_l[d] = (data[base + regions_[r].base].GetLowerBound() - translation) * rscaling;
        x_u[d++] = (data[base + regions_[r].base].GetUpperBound() - translation) * rscaling;
      }else{
        const unsigned int u = base + regions_[r].base + regions_[r].elements;
        for(unsigned int s = base + regions_[r].base; s < u; s++){
          x_l[d] = (data[s].GetLowerBound() - translation) * rscaling;
          x_u[d++] = (data[s].GetUpperBound() - translation) * rscaling;
        }
      }
    }
  }
  assert(d == DesignSpace::GetNumberOfVariables());
}

void DesignSpace::WriteGradientToExtern(double* out, DesignElement::ValueSpecifier vs, DesignElement::Access access, Condition* constraint, bool use_scaling) const
{
  // this does now do reordering as gradients are reordered in the optimizer
  // must be set in the constructor! might be trivial volume fraction or from file!!
  unsigned int d = 0;
  for(unsigned int des = 0; des < design.GetSize(); des++){
    const unsigned int base = des * elements_;
    for(unsigned int r = 0; r < regions_.GetSize(); r++){
      const double scaling = use_scaling ? scale_design[des][r] : 1.0;
      const unsigned int u = base + regions_[r].base + regions_[r].elements;
      if(regions_[r].constant){
        out[d] = 0;
        for(unsigned int s = base + regions_[r].base; s < u; s++){
          out[d] += (vs == DesignElement::COST_GRADIENT ? data[s].GetValue(vs, access) : data[s].GetConstraintGradient(constraint)) * scaling;
        }
        d++;
      }else{
        for(unsigned int s = base + regions_[r].base; s < u; s++){
          out[d++] = (vs == DesignElement::COST_GRADIENT ? data[s].GetValue(vs, access) : data[s].GetConstraintGradient(constraint)) * scaling;
        }
      }
    }
  }
  assert(d == DesignSpace::GetNumberOfVariables()); // this does only work on the non-shape part of the gradient, the rest is done in ShapeDesign
}

void DesignSpace::Reset(DesignElement::ValueSpecifier vs, DesignElement::Type design)
{
  unsigned int start = design == DesignElement::DEFAULT ? 0 : FindDesign(design) * elements_;
  unsigned int end   = design == DesignElement::DEFAULT ? data.GetSize() : start + elements_;

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

      case DesignElement::COST_GRADIENT:
           de.SetObjectiveGradient(0.0);
           break;

      default: throw Exception("value specifier not handled");
    }
  }
}
DesignElement* DesignSpace::Find(unsigned int elemNum, DesignElement::Type dt, bool throw_exception)
{
  int idx = Find(elemNum, throw_exception);
  if(idx == -1) return NULL;
  for(unsigned int d = 0; d < design.GetSize(); d++)
  {
    unsigned int pos = elements_ * d + idx;
    if(data[pos].GetType() == dt)
    {
      if(data[pos].elem->elemNum != elemNum) throw Exception("index mixed up");
      return &data[pos];
    }
  }

  if(throw_exception) throw Exception("design type not in design");
  return NULL;
}


int DesignSpace::Find(unsigned int elemNum, bool throw_exception)
{
  // to easy find optimizing. Check if we are at the same element as
  // last time
  if(data[last_find_index_].elem->elemNum == elemNum)
    return last_find_index_;

  // check if we are at the next element than last time.
  // if we are at the boundary the loop does the wrap around
  if(++last_find_index_ < elements_ && data[last_find_index_].elem->elemNum == elemNum)
    return last_find_index_;

  // search all
  for(unsigned int i = 0; i < elements_; i++)
  {
    if(data[i].elem->elemNum == elemNum)
    {
      last_find_index_ = i;
      return i;
    }
  }

  if(throw_exception) EXCEPTION("could not find element " << elemNum << " in our design space");
  return -1;
}

int DesignSpace::Find(const Elem* elem, bool throw_exception)
{
  if(FindRegion(elem->regionId) >= 0) return Find(elem->elemNum);

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

void DesignSpace::ToInfo(InfoNode* in)
{
  InfoNode* tf = in->Get("transferFunctions");
  for(unsigned int i = 0; i < transfer.GetSize(); i++)
    transfer[i].ToInfo(tf->Get("transferFunction", InfoNode::APPEND));

  InfoNode* dv = in->Get("designVariables");
  dv->Get("totalDesignVariables")->SetValue(data.GetSize());
  for(unsigned int i = 0; i < design.GetSize(); i++)
    data[i * elements_].ToInfo(dv->Get("design", InfoNode::APPEND));

  InfoNode* rs = in->Get("regions");
  for(unsigned int i = 0; i < regions_.GetSize(); i++)
  {
    InfoNode* r = rs->Get("region", InfoNode::APPEND);
    r->Get("name")->SetValue(domain->GetGrid()->RegionIdToName(regions_[i].regionId));
    r->Get("elements")->SetValue(regions_[i].elements);
  }
}


std::string DesignSpace::ToString()
{
  std::stringstream ss;
  ss << "design[";
  for(unsigned int i = 0; i < data.GetSize(); i++)
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

unsigned int DesignSpace::GetNumberOfVariables() const {
  unsigned int n = 0;
  for(unsigned int r = 0; r < regions_.GetSize(); r++){
    n += regions_[r].constant ? 1 : regions_[r].elements;
  }
  return design.GetSize() * n;
}

int DesignSpace::FindRegion(RegionIdType regionId){
  for(unsigned int r = 0; r < regions_.GetSize(); r++){
    if(regions_[r].regionId == regionId){
      return r;
    }
  }
  return -1;
}
