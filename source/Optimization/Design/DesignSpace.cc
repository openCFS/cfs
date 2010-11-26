#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/ShapeDesign.hh"
#include "Optimization/TransferFunction.hh"
#include "Optimization/Condition.hh"
#include "Optimization/Optimizer/BaseOptimizer.hh"
#include "Optimization/Optimizer/ShapeOptimizer.hh"
#include "Optimization/LevelSet.hh"
#include "General/exception.hh"
#include "General/Enum.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"
#include "Domain/resultInfo.hh"
#include "Utils/result.hh"
#include "Utils/StdVector.hh"
#include "PDE/SinglePDE.hh"
#include "DataInOut/Logging/cfslog.hh"

using namespace CoupledField;

using std::complex;

// declare class specific logging stream
DECLARE_LOG(designSpace)
DEFINE_LOG(designSpace, "designSpace")

// declare class specific logging stream
DECLARE_LOG(ersatz)
DEFINE_LOG(ersatz, "ersatzMaterialFactor")



DesignSpace::DesignSpace(StdVector<RegionIdType>& regionIds, ParamNodeList &pn_design, ParamNodeList &trans_in, ParamNodeList &result, ErsatzMaterial::Method method)
{
  LOG_TRACE(designSpace) << "DesignSpace for regions=" << regionIds << " #designs=" << pn_design.GetSize()
                         << " #transferFunctions=" << trans_in.GetSize() << " #results=" << result.GetSize();
  regions.Resize(regionIds.GetSize());
  all_regions_regular_ = domain->GetGrid()->IsRegionRegular(regionIds);

  design_id = 0;
  optimizer_ = NULL;

  designMaterial = NULL;

  applicationForm.SetName("DesignSpace::ApplicationForm");
  applicationForm.Add(Optimization::ELEC, "linElecInt");
  applicationForm.Add(Optimization::MECH, "linElastInt");
  // We follow for the stress, strain calculation the transfer functions of mech
  applicationForm.Add(Optimization::MECH, "MechStressStrain", false);
  applicationForm.Add(Optimization::PIEZO_COUPLING, "linPiezoCoupling");
  applicationForm.Add(Optimization::CHARGE_DENSITY, "LinNeumannInt");
  applicationForm.Add(Optimization::PRESSURE, "PressureLinForm");
  applicationForm.Add(Optimization::MASS, "MassInt");
  applicationForm.Add(Optimization::HEAT, "LaplaceInt");

  // read the elements
  elements = domain->GetGrid()->GetNumElems(regionIds);
  
  // should we adapt the lower bound of the design variable to the penalty exponent of the transfer function? 
  bool adapt_lower(false);

  // number of different designs
  unsigned int nd = 0;
  for(unsigned int d = 0; d < pn_design.GetSize(); d++){
    DesignElement::Type dt = DesignElement::type.Parse(pn_design[d]->Get("name")->As<std::string>());
    if(design.Find(dt) < 0){
      design.Push_back(dt);
      nd++;
    }
  }
  
  // store the CFS element (number) to design element mapping.
  // Used by Find() and the filter and vicinity neighbors
  elemToDesign.Resize(domain->GetGrid()->GetNumElems() + 1, -1); // 1 based.

  if(elements == 0 || nd == 0)
  { // this may happen in shape optimization
    if(method != ErsatzMaterial::SHAPE_OPT && method != ErsatzMaterial::SHAPE_PARAM_MAT
        && method !=ErsatzMaterial::SHAPE_GRAD)
    {
      if(elements == 0) throw Exception("empty regions");
      if(nd == 0) throw Exception("no designs types given.");
    }
  }
  else // 'standard' SIMP case
  {
    // set our own structure with is element times design parameters
    data.Reserve(elements * nd);

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

    // only temporary
    StdVector<Elem*> elems;

    // iterate over all designs, we have to always have the same order of designs not the order from the xml file
    for(unsigned int dti = 0; dti < design.GetSize(); dti++)
    {
      DesignElement::Type dt = design[dti];
      for(unsigned int d = 0; d < pn_design.GetSize(); d++)
      {
        PtrParamNode curr_design_pn = pn_design[d];
        if(DesignElement::type.Parse(curr_design_pn->Get("name")->As<std::string>()) != dt) continue;

        std::string design_reg = curr_design_pn->Get("region")->As<std::string>();

        // there is no default for adapt_lower in the load ersatz material case
        adapt_lower = false;
        if(curr_design_pn->Has("adapt_lower"))
          adapt_lower = curr_design_pn->Get("adapt_lower")->As<bool>();
        
        for(unsigned int r = 0; r < regionIds.GetSize(); r++)
        {
          domain->GetGrid()->GetElems(elems, regionIds[r]);
          unsigned int n = elems.GetSize();
          std::string reg = regNames[regionIds[r]];

          // do this only for the first design per region (?)
          if(d == 0){
            regions[r].regionId = regionIds[r];
            regions[r].base = r == 0 ? 0 : regions[r-1].base + regions[r-1].elements;
            regions[r].elements = n;
            regions[r].constant = false;
          }

          if(design_reg == "all" || design_reg == reg)
          {
            int di = design.Find(dt);
            if(!region_design[r*nd + di]){
              throw Exception("Design/Region combination given twice!");
            }
            region_design[r*nd + di] = false;
            if(curr_design_pn->Has("constant") && curr_design_pn->Get("constant")->As<bool>()){
              // we have a constant densign-value on that region
              regions[r].constant = true;
            }
            
            if(curr_design_pn->Has("scale") && curr_design_pn->Get("scale")->As<bool>()){
              double upper = curr_design_pn->Get("upper")->As<double>();
              double lower = curr_design_pn->Get("lower")->As<double>();
              scale_design[di][r] = (upper - lower);
              translate_design[di][r] = lower;
            }

            double initial = curr_design_pn->Get("initial")->As<double>();
            LOG_DBG2(designSpace) << "add design " << dt << ":" << DesignElement::type.ToString(dt)
              << " initial=" << initial;
                        
            for(unsigned int e = 0; e < n; e++)
            {
              DesignElement de(curr_design_pn, elems[e]);
              de.SetDesign(initial);
              
              data.Push_back(de);

              // append rucksack :)
              if(method == ErsatzMaterial::SIMP_METHOD)
              {
                DesignElement* ptr = &(data.Last());
                ptr->simp = new SIMPElement(ptr);
              }

              // store the element mapping only for the first design
              int idx = (int) data.GetSize() - 1;
              if(idx < (int) elements)
                elemToDesign[de.elem->elemNum] = idx;
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
  resultDescriptions.Reserve(result.GetSize());
  for(unsigned int i = 0; i < result.GetSize(); i++)
    resultDescriptions.Push_back(ResultDescription(result[i]));
  
  bimattensor_ = Matrix<double>();
  
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
}

DesignSpace::~DesignSpace(){
  if(designMaterial != NULL){
    delete designMaterial;
    designMaterial = NULL;
  }
}

DesignSpace* DesignSpace::CreateInstance(StdVector<RegionIdType> regionIds, ParamNodeList &design, ParamNodeList& transfer, ParamNodeList& result, ErsatzMaterial::Method method){
  switch(method){
  case ErsatzMaterial::SHAPE_OPT:
  case ErsatzMaterial::SHAPE_PARAM_MAT:
    return new ShapeDesign(regionIds, design, transfer, result, method);
  default:
    return new DesignSpace(regionIds, design, transfer, result, method);
  }
}

void DesignSpace::PostInit(int objectives, int constraints)
{
  LOG_DBG(designSpace) << "# objectives = " << objectives << ", # constraints = " << constraints;
  DesignElement::SetDesignSpace(this);

  for(unsigned int i = 0; i < data.GetSize(); i++)
    data[i].PostInit(objectives, constraints);
}

bool DesignSpace::Contains(const RegionIdType reg) const
{
  for(unsigned int i = 0, n = regions.GetSize(); i < n; i++)
    if(regions[i].regionId == reg)
      return true;

  return false;
}

void DesignSpace::SetDesignMaterial(PtrParamNode dm){
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

  // in info.xml result output the " (" is replaced by "_("!
  ri->resultName = DesignElement::valueSpecifier.ToString(rd.value) + "_"
                   + (rd.detail != DesignElement::NONE ? (DesignElement::detail.ToString(rd.detail) + "_") : "")
                   + DesignElement::type.ToString(rd.design) + " ("
                   + DesignElement::access.ToString(rd.access) + ")"
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
      double transformed = tf->Transform(de, DesignElement::SMART); // handles design filtering
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
  unsigned int s = 0;
  for(unsigned int des = 0; des < design.GetSize(); des++){
    const unsigned int base = des * elements;
    for(unsigned int r = 0; r < regions.GetSize(); r++){
      const double scaling = scale_design[des][r];
      const double translation = translate_design[des][r];
      const unsigned int u = base + regions[r].base + regions[r].elements;
      if(regions[r].constant){
        const double v = space[s] * scaling + translation;
        for(unsigned int d = base + regions[r].base; d < u; d++){
          if(!new_design && data[d].GetDesign(DesignElement::PLAIN) != v) {
            new_design = true;
          }
          LOG_DBG3(designSpace) << "ReadDesignFromExtern: setting design (constant region): data[" << d << "] = " << v;
          data[d].SetDesign(v);
        } // for d
        s++; // only advance after having set all element of this region to the corresponding value
      }else{
        for(unsigned int d = base + regions[r].base; d < u; d++){
          double v = space[s] * scaling + translation;
          if(!new_design && data[d].GetDesign(DesignElement::PLAIN) != v) {
            new_design = true;
          }
          LOG_DBG3(designSpace) << "ReadDesignFromExtern: setting design: data[" << d << "] = " << v;
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


bool DesignSpace::CompareDesign(const double* space)
{
  unsigned int s = 0;
  for(unsigned int des = 0; des < design.GetSize(); des++){
    const unsigned int base = des * elements;
    for(unsigned int r = 0; r < regions.GetSize(); r++){
      const double scaling = scale_design[des][r];
      const double translation = translate_design[des][r];
      const unsigned int u = base + regions[r].base + regions[r].elements;
      if(regions[r].constant){
        const double v = space[s] * scaling + translation;
        for(unsigned int d = base + regions[r].base; d < u; d++){
          if(data[d].GetDesign(DesignElement::PLAIN) != v)
            return false;
        } // for d
        s++; // only advance after having set all element of this region to the corresponding value
      }else{
        for(unsigned int d = base + regions[r].base; d < u; d++){
          double v = space[s] * scaling + translation;
          if(data[d].GetDesign(DesignElement::PLAIN) != v)
            return false;
          s++; // advance in every step
        } // for d
      } // if/else constant
    } // for r
  } // for des
  assert(s == DesignSpace::GetNumberOfVariables());
  return true;
}

int DesignSpace::WriteDesignToExtern(double* space, bool scaling) const
{
  unsigned int d = 0;
  for(unsigned int des = 0; des < design.GetSize(); des++){
    const unsigned int base = des * elements;
    for(unsigned int r = 0; r < regions.GetSize(); r++){
      const double rscaling = scaling ? 1.0 / scale_design[des][r] : 1.0;
      const double translation = scaling ? translate_design[des][r] : 0.0;
      if(regions[r].constant){
        LOG_DBG3(designSpace) << "WriteDesignToExtern: constant region " << r << " design[" << base + regions[r].base << "]=" << data[base + regions[r].base].GetDesign(DesignElement::PLAIN);
        space[d++] = (data[base + regions[r].base].GetDesign(DesignElement::PLAIN) - translation) * rscaling;
      }else{
        const unsigned int u = base + regions[r].base + regions[r].elements;
        for(unsigned int s = base + regions[r].base; s < u; s++){
          LOG_DBG3(designSpace) << "WriteDesignToExtern: non-constant region " << r << " design[" << s << "]=" << data[s].GetDesign(DesignElement::PLAIN);
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
    const unsigned int base = des * elements;
    for(unsigned int r = 0; r < regions.GetSize(); r++){
      const double rscaling = 1.0 / scale_design[des][r];
      const double translation = translate_design[des][r];
      if(regions[r].constant){
        x_l[d] = (data[base + regions[r].base].GetLowerBound() - translation) * rscaling;
        x_u[d++] = (data[base + regions[r].base].GetUpperBound() - translation) * rscaling;
      }else{
        const unsigned int u = base + regions[r].base + regions[r].elements;
        for(unsigned int s = base + regions[r].base; s < u; s++){
          x_l[d] = (data[s].GetLowerBound() - translation) * rscaling;
          x_u[d++] = (data[s].GetUpperBound() - translation) * rscaling;
        }
      }
    }
  }
  assert(d == DesignSpace::GetNumberOfVariables());
}

void DesignSpace::WriteSparseGradientToExtern(StdVector<double>& out, DesignElement::ValueSpecifier vs, DesignElement::Access access, Condition* g, bool use_scaling) const
{
  // Bastian did some complicated reordering stuff. For the only case of sparse Jacobians (slope constraints)
  // we'll have only the simple standard situtation .. if this changes you have at least a test case :) Fabian
  
  assert(regions.GetSize() == 1);
  assert(design.GetSize() == 1);
  assert(g != NULL); // only constraints can have sparse Jacobians 

  const double scaling = use_scaling ? scale_design[0][0] : 1.0;
  
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

  unsigned int n = out.window.GetStart(); // to grow up to the total number of design variables

  for(unsigned int des = 0; des < design.GetSize(); des++)
  {
    const unsigned int base = des * elements;

    for(unsigned int r = 0; r < regions.GetSize(); r++)
    {
      const double scaling = use_scaling ? scale_design[des][r] : 1.0;
      const unsigned int region_elements = base + regions[r].base + regions[r].elements;

      if(regions[r].constant)
      {
        out[n] = 0;
        for(unsigned int s = base + regions[r].base; s < region_elements; s++)
        {
          assert(out.InWindow(n));
          out[n] += data[s].GetValue(vs, access, g) * scaling;
          LOG_DBG3(designSpace) << "WriteDenseGradientToExtern: constant region " << r << " design[" << s << "]=" << data[s].GetValue(vs, access, g);
        }
        LOG_DBG3(designSpace) << "WriteDenseGradientToExtern: constant region " << r << " sum = " << out[n] / scaling;
        n++;
      }
      else
      {
        for(unsigned int s = base + regions[r].base; s < region_elements; s++)
        {
          assert(out.InWindow(n));
          LOG_DBG3(designSpace) << "WriteDenseGradientToExtern: non-constant region " << r << " design[" << s << "]=" << data[s].GetValue(vs, access, g);
          out[n++] = data[s].GetValue(vs, access, g) * scaling;
        }
      }
    }
  }
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
DesignElement* DesignSpace::Find(unsigned int elemNum, DesignElement::Type dt, bool throw_exception)
{
  int idx = Find(elemNum, throw_exception);
  if(idx == -1) return NULL; // an exception was already thrown if desired
  for(unsigned int d = 0; d < design.GetSize(); d++)
  {
    unsigned int pos = elements * d + idx;
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
  int idx = elemToDesign[elemNum];

  if(idx == -1 && throw_exception)
    EXCEPTION("could not find element " << elemNum << " in our design space");

  return idx;
}

int DesignSpace::Find(const DesignElement* de, bool throw_exception)
{
  int idx = Find(de->elem->elemNum, throw_exception);

  for(unsigned int d = 0; idx > -1 && d < design.GetSize(); d++)
    if(data[elements * d + idx].GetType() == de->GetType())
      return idx;

  // had wrong design type
  if(throw_exception)
    EXCEPTION("requested design element num has index " << idx << " but invalid design type " << de->GetType());

  return -1;
}

int DesignSpace::Find(const Elem* elem, bool throw_exception)
{
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

  PtrParamNode rs = in->Get("regions");
  for(unsigned int i = 0; i < regions.GetSize(); i++)
  {
    PtrParamNode r = rs->Get("region", ParamNode::APPEND);
    r->Get("name")->SetValue(domain->GetGrid()->GetRegion().ToString(regions[i].regionId));
    r->Get("elements")->SetValue(regions[i].elements);
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
  unsigned int n = 0;
  for(unsigned int r = 0; r < regions.GetSize(); r++){
    n += regions[r].constant ? 1 : regions[r].elements;
  }
  return design.GetSize() * n;
}

int DesignSpace::FindRegion(RegionIdType regionId){
  for(unsigned int r = 0; r < regions.GetSize(); r++){
    if(regions[r].regionId == regionId){
      return r;
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
  def.design = ri->resultType == ELEC_PSEUDO_POLARIZATION ? DesignElement::POLARIZATION : DesignElement::DENSITY;
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
      unsigned int data_index = (base * elements) + base_index;
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
