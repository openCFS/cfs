#include "Optimization/DesignSpace.hh"
#include "Optimization/DesignElement.hh"
#include "Optimization/TransferFunction.hh"
#include "Optimization/Condition.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
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

// declare class specific logging stream
DECLARE_LOG(designSpace)
DEFINE_LOG(designSpace, "designSpace")

// declare class specific logging stream
DECLARE_LOG(ersatz)
DEFINE_LOG(ersatz, "ersatzMaterial")

DesignSpace::DesignSpace(int regionId, StdVector<ParamNode*>& pn_design, StdVector<ParamNode*>& trans_in, StdVector<ParamNode*>& result)
{
  LOG_TRACE(designSpace) << "DesignSpace for region=" << regionId << " #designs=" << pn_design.GetSize() 
                         << " #transferFunctions=" << trans_in.GetSize() << " #results=" << result.GetSize();
  last_find_index_ = 0;
  
  regionId_ = regionId;

  // only temporary 
  StdVector<Elem*> elems;

  applicationForm.SetName("DesignSpace::ApplicationForm");
  applicationForm.Add(Optimization::ELEC, "linElecInt");
  applicationForm.Add(Optimization::MECH, "linElastInt");
  applicationForm.Add(Optimization::PIEZO_COUPLING, "linPiezoCoupling");
  applicationForm.Add(Optimization::CHARGE_DENSITY, "LinNeumannInt");
  applicationForm.Add(Optimization::PRESSURE, "PressureLinForm");  

  // read the elements
  domain->GetGrid()->GetElems(elems, regionId);
  elements_ = elems.GetSize();
  if(elements_ == 0) throw Exception("empty region");

  // set our own structure with is element times design parameters
  data.Reserve(elements_ * pn_design.GetSize());
  LOG_DBG(designSpace) << "data size: " << elements_ << "*" << pn_design.GetSize() << "=" << elements_ * pn_design.GetSize(); 
   
  // initialize for all designs for all elements with type, inital and element
  design.Resize(pn_design.GetSize());
  for(unsigned int d = 0; d < pn_design.GetSize(); d++)
  { 
    DesignElement::Type dt = (DesignElement::Type) DesignElement::type.Parse(pn_design[d]->Get("name")->AsString());
    double initial = pn_design[d]->Get("initial")->AsDouble();
    design[d]=dt;
    LOG_DBG2(designSpace) << "add design " << dt << ":" << DesignElement::type.ToString(dt);
    
    for(unsigned int e = 0; e < elements_; e++) 
    {
      DesignElement de;
      data.Push_back(de);
      data.Last().SetType(dt);
      data.Last().elem = elems[e];
      data.Last().SetDesign(initial);
    }
  } 
  
  // now read the transfer functions
  transfer.Reserve(trans_in.GetSize());
  if(trans_in.GetSize() < pn_design.GetSize()) 
    throw Exception("less transferFunctions than design variable types is inveasible");
  
  for(unsigned int i = 0; i < trans_in.GetSize(); i++)
  {
    transfer.Push_back(TransferFunction(trans_in[i]));
  }
  
  // set the result descriptions which identify the solution types
  resultDescriptions_.Reserve(result.GetSize());
  for(unsigned int i = 0; i < result.GetSize(); i++)
    resultDescriptions_.Push_back(ResultDescription(result[i]));
}

void DesignSpace::PostInit(int constraints)
{
  for(unsigned int i = 0; i < data.GetSize(); i++)
    data[i].constraintGradient.Resize(constraints);
}

void DesignSpace::AppendOptimizationResults(SinglePDE* pde)
{
  // set the result descriptions which identify the solution types
  for(unsigned int i = 0; i < resultDescriptions_.GetSize(); i++)
  {
    ResultInfo* ri = GetResultInfo(resultDescriptions_[i]);
    shared_ptr<ResultInfo> opt_res(ri);    
    bool added = pde->CheckStoreResult(opt_res);
    if(!added) Warning("Warning: optimization result not as pde result declared"); 
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

  ri->dofNames = "";
  ri->unit = "";
  ri->entryType = ResultInfo::SCALAR;
  ri->definedOn = ResultInfo::ELEMENT;
  ri->fctType = shared_ptr<ConstFct>(new ConstFct() );
  
  // let the caller or a shared pointer delete it finally
  return ri;
}

int DesignSpace::GetSpecialResultIndex(DesignElement::Type design, DesignElement::ValueSpecifier value, DesignElement::Detail detail)
{
  for(unsigned int i = 0; i < resultDescriptions_.GetSize(); i++)
  {
    const ResultDescription& rd = resultDescriptions_[i];
    if(rd.design != design || rd.value != value || rd.detail != detail) continue;
    
    // we are right. 
    switch(rd.solutionType)
    {
      case OPT_RESULT_1: return 0;
      case OPT_RESULT_2: return 1;
      case OPT_RESULT_3: return 2;
      default: throw Exception("invalid solution type");
    }
  }

  return -1; // the specified tripel was not specified such in xml
}

void DesignSpace::ExtractResults(shared_ptr<BaseResult> base_result)
{
  // our results are up to now scalar!
  Result<double>& result = dynamic_cast<Result<double>&>(*base_result);
  Vector<double>& result_data = result.GetVector();
  
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
  for(unsigned int i = 0; i < resultDescriptions_.GetSize(); i++)
    if(resultDescriptions_[i].solutionType == ri->resultType) 
      descr = resultDescriptions_[i];
  
  // search where in data we are
  int base = FindDesign(descr.design);

  // loop over elements from result. We have to do it this way as the the connection
  // of design element and result element is the element(->elemeNum) but we cannot
  // search in the result for an element.
  EntityIterator it = result.GetEntityList()->GetIterator();
  for ( it.Begin(); !it.IsEnd(); it++ ) 
  {
    // note that the index is from the first design set!
    unsigned int base_index = Find(it.GetElem()->elemNum);
    // base=0 is first!
    unsigned int data_index = (base * elements_) + base_index; 
    DesignElement& de = data[data_index];
  
    #ifdef CHECK_INDEX
    if(de.elem->elemNum != it.GetElem()->elemNum) {
      EXCEPTION("mixed up indices:" << de.elem->elemNum << "!=" << it.GetElem()->elemNum
                << " base_index=" << base_index << " data_index=" << data_index << " it.Pos()=" << it.GetPos());
    }
    #endif    
    
    result_data[it.GetPos()] = de.GetValue(&descr);
  }
}

int DesignSpace::FindDesign(DesignElement::Type dt, bool throw_exception)
{
  // do a fallback for NO_TYPE and DEFAULT
  if(design.GetSize() == 1 && (dt == DesignElement::NO_TYPE || dt == DesignElement::DEFAULT))
    return 0; 
  
  // search where in data we are
  int base = -1;
  for(unsigned int i = 0; i < design.GetSize(); i++)
    if(design[i] == dt) base = i;
  
  if(base == -1 && throw_exception) 
    EXCEPTION("Design " << DesignElement::type.ToString(dt) << " not within " << design.GetSize() << " actual designs.");
                   
  return base;                    
} 


double DesignSpace::GetErsatzMaterialFactor(const Elem* elem, const BaseForm* form)
{
  // are we in the relevant region at all?
  if(elem->regionId != regionId_) 
    throw Exception("the element has the wrong region id");
  

  // the application we have (MECH, ELEC, PIEZO_COUPLING)
  Optimization::Application applic = (Optimization::Application) applicationForm.Parse(form->GetName());
  // Start with finding the base DE
  unsigned int design_index = Find(elem->elemNum);
  
  return GetErsatzMaterialFactor(design_index, applic);
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

TransferFunction* DesignSpace::GetTransferFunction(DesignElement::Type design, Optimization::Application application, bool throw_exception)
{
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

void DesignSpace::ReadDesignFromExtern(const double* space)
{
   for(unsigned int i = 0; i < data.GetSize(); i++) 
       data[i].SetDesign(space[i]);
}  

void DesignSpace::WriteDesignToExtern(double* space)
{
    // must be set in the constructor! might be trivial volume fraction or from file!!
    for(unsigned int i = 0; i < data.GetSize(); i++)
       space[i] = data[i].GetDesign(DesignElement::PLAIN);
}  


DesignElement* DesignSpace::Find(unsigned int elemNum, DesignElement::Type dt)
{
  unsigned int idx = Find(elemNum);
  for(unsigned int d = 0; d < design.GetSize(); d++)
  {
    unsigned int pos = elements_ * d + idx; 
    if(data[pos].GetType() == dt)
    {
      if(data[pos].elem->elemNum != elemNum) throw Exception("index mixed up");
      return &data[pos]; 
    }
  }
  
  throw Exception("design type not in design");
}

unsigned int DesignSpace::Find(unsigned int elemNum)
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
    
  throw Exception("could not find element in our design space");  
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

  
       
