#include "Optimization/DesignSpace.hh"
#include "Optimization/DesignElement.hh"
#include "Optimization/TransferFunction.hh"
#include "Optimization/Condition.hh"
#include "Optimization/BaseOptimizer.hh"
#include "Optimization/ShapeOptimizer.hh"
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
#include "Optimization/ShapeDesign.hh"

using namespace CoupledField;

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
  regions_.Resize(regionIds.GetSize());
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
  elements_ = domain->GetGrid()->GetNumElems(regionIds);

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

  if(elements_ == 0 || nd == 0)
  { // this may happen in shape optimization
    if(method != ErsatzMaterial::SHAPE_OPT && method != ErsatzMaterial::SHAPE_PARAM_MAT
        && method !=ErsatzMaterial::SHAPE_GRAD)
    {
      if(elements_ == 0) throw Exception("empty regions");
      if(nd == 0) throw Exception("no designs given");
    }
  }
  else // 'standard' SIMP case
  {
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

        for(unsigned int r = 0; r < regionIds.GetSize(); r++)
        {
          domain->GetGrid()->GetElems(elems, regionIds[r]);
          unsigned int n = elems.GetSize();
          std::string reg = regNames[regionIds[r]];

          // do this only for the first design per region (?)
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
            if(curr_design_pn->Get("constant")->As<bool>()){ // we have a constant densign-value on that region
              regions_[r].constant = true;
            }
            if(curr_design_pn->Get("scale")->As<bool>()){
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
              if(idx < (int) elements_)
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

  } // no design elements given

  // set the result descriptions which identify the solution types
  resultDescriptions.Reserve(result.GetSize());
  for(unsigned int i = 0; i < result.GetSize(); i++)
    resultDescriptions.Push_back(ResultDescription(result[i]));
  
  bimattensor_ = Matrix<double>();
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
  DesignElement::SetDesignSpace(this);

  for(unsigned int i = 0; i < data.GetSize(); i++)
    data[i].PostInit(objectives, constraints);
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

  ri->resultName = DesignElement::valueSpecifier.ToString(rd.value) + "_"
                   + (rd.detail != DesignElement::NONE ? (DesignElement::detail.ToString(rd.detail) + "_") : "")
                   + DesignElement::type.ToString(rd.design) + " ("
                   + DesignElement::access.ToString(rd.access) + ")";


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

  // nothing found, another chance
  if(transfer.GetSize() == 1 && design == DesignElement::DEFAULT) // kind of dangerous without application check!!
    return &transfer[0];

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

void DesignSpace::WriteSparseGradientToExtern(StdVector<double>& out, DesignElement::ValueSpecifier vs, DesignElement::Access access, Condition* g, bool use_scaling) const
{
  // Bastian did some complicated reordering stuff. For the only case of sparse Jacobians (slope constraints)
  // we'll have only the simple standard situtation .. if this changes you have at least a test case :) Fabian
  
  assert(regions_.GetSize() == 1);
  assert(design.GetSize() == 1);
  assert(g != NULL); // only constraints can have sparse Jacobians 

  const double scaling = use_scaling ? scale_design[0][0] : 1.0;
  
  StdVector<unsigned int>& sparsity = g->GetSparsityPattern();
  
  assert(out.window.GetSize() == sparsity.GetSize());
  unsigned int base = out.window.GetStart();
  for(unsigned int i = 0; i < sparsity.GetSize(); i++)
  {
    assert(out.InWindow(base + i));
    out[base + i] = data[sparsity[i]].GetGradient(NULL, g) * scaling;
  }
}

void DesignSpace::WriteDenseGradientToExtern(StdVector<double>& out, DesignElement::ValueSpecifier vs, DesignElement::Access access, Condition* g, bool use_scaling) const
{
  // this does now do reordering as gradients are reordered in the optimizer
  // must be set in the constructor! might be trivial volume fraction or from file!!
  bool cost_grad = vs == DesignElement::COST_GRADIENT;
  assert(!(cost_grad && g != NULL));

  unsigned int n = out.window.GetStart(); // to grow up to the total number of design variables

  for(unsigned int des = 0; des < design.GetSize(); des++)
  {
    const unsigned int base = des * elements_;

    for(unsigned int r = 0; r < regions_.GetSize(); r++)
    {
      const double scaling = use_scaling ? scale_design[des][r] : 1.0;
      const unsigned int region_elements = base + regions_[r].base + regions_[r].elements;

      if(regions_[r].constant)
      {
        out[n] = 0;
        for(unsigned int s = base + regions_[r].base; s < region_elements; s++)
        {
          assert(out.InWindow(n));
          out[n] += (cost_grad ? data[s].GetValue(vs, access) : data[s].GetGradient(NULL, g)) * scaling;
        }
        n++;
      }
      else
      {
        for(unsigned int s = base + regions_[r].base; s < region_elements; s++)
        {
          assert(out.InWindow(n));
          out[n++] = (cost_grad ? data[s].GetValue(vs, access) : data[s].GetGradient(NULL, g)) * scaling;
        }
      }
    }
  }
}

void DesignSpace::Reset(DesignElement::ValueSpecifier vs, DesignElement::Type design)
{
  unsigned int start = design == DesignElement::DEFAULT || DesignElement::TENSOR_TRACE ? 0 : FindDesign(design) * elements_;
  unsigned int end   = design == DesignElement::DEFAULT || DesignElement::TENSOR_TRACE ? data.GetSize() : start + elements_;

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
  int idx = elemToDesign[elemNum];

  if(idx == -1 && throw_exception)
    EXCEPTION("could not find element " << elemNum << " in our design space");

  return idx;
}

int DesignSpace::Find(const DesignElement* de, bool throw_exception)
{
  int idx = Find(de->elem->elemNum, throw_exception);

  for(unsigned int d = 0; idx > -1 && d < design.GetSize(); d++)
    if(data[elements_ * d + idx].GetType() == de->GetType())
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
    data[i * elements_].ToInfo(dv->Get("design", ParamNode::APPEND));

  PtrParamNode rs = in->Get("regions");
  for(unsigned int i = 0; i < regions_.GetSize(); i++)
  {
    PtrParamNode r = rs->Get("region", ParamNode::APPEND);
    r->Get("name")->SetValue(domain->GetGrid()->GetRegion().ToString(regions_[i].regionId));
    r->Get("elements")->SetValue(regions_[i].elements);
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
