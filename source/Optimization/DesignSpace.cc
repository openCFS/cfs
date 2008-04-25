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
DEFINE_LOG(ersatz, "ersatzMaterialFactor")

DesignSpace::DesignSpace(int regionId, StdVector<ParamNode*>& pn_design, StdVector<ParamNode*>& trans_in, StdVector<ParamNode*>& result, ErsatzMaterial::Method method)
{
  LOG_TRACE(designSpace) << "DesignSpace for region=" << regionId << " #designs=" << pn_design.GetSize() 
                         << " #transferFunctions=" << trans_in.GetSize() << " #results=" << result.GetSize();
  last_find_index_ = 0;
  
  regionId_ = regionId;
  design_id = 0;

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
  applicationForm.Add(Optimization::SURFACE_NORMAL, "SurfaceNormalInt");

  // read the elements
  domain->GetGrid()->GetElems(elems, regionId);
  elements_ = elems.GetSize();
  if(elements_ == 0) throw Exception("empty region");

  // set our own structure with is element times design parameters
  data.Reserve(elements_ * pn_design.GetSize());
  
  // initialize for all designs for all elements with type, inital and element
  design.Resize(pn_design.GetSize());
  for(unsigned int d = 0; d < pn_design.GetSize(); d++)
  { 
    DesignElement::Type dt = DesignElement::type.Parse(pn_design[d]->Get("name"));
    design[d]=dt;
    double initial = pn_design[d]->Get("initial")->AsDouble();
    LOG_DBG2(designSpace) << "add design " << dt << ":" << DesignElement::type.ToString(dt)
                          << " initial=" << initial;
    
    for(unsigned int e = 0; e < elements_; e++) 
    {
      DesignElement de(pn_design[d]);
      de.elem = elems[e];
      de.SetDesign(initial);
      data.Push_back(de);
      switch(method)
      {
      case ErsatzMaterial::SIMP_METHOD:
        data.Last().simp = new SIMPElement(&data.Last());
        break;
      case ErsatzMaterial::NO_METHOD:
        break;
      default: assert(false);
      }
    }
  } 

  LOG_DBG(designSpace) << "data size: " << elements_ << "*" << pn_design.GetSize() << "=" << data.GetSize(); 
  
  
  // now read the transfer functions
  transfer.Reserve(trans_in.GetSize());
  if(trans_in.GetSize() < pn_design.GetSize()) 
    throw Exception("less transferFunctions than design variable types is inveasible");
  
  for(unsigned int i = 0; i < trans_in.GetSize(); i++)
  {
    transfer.Push_back(TransferFunction(trans_in[i]));
  }
  
  // set the result descriptions which identify the solution types
  resultDescriptions.Reserve(result.GetSize());
  for(unsigned int i = 0; i < result.GetSize(); i++)
    resultDescriptions.Push_back(ResultDescription(result[i]));
}

void DesignSpace::PostInit(int constraints)
{
  for(unsigned int i = 0; i < data.GetSize(); i++)
    data[i].constraintGradient.Resize(constraints);
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
      os << "Warning: optimization result '" << rd.ToString() << "' not as pde result declared";
      Warning(os.str().c_str());
    
    }		     
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

int DesignSpace::ReadDesignFromExtern(const double* space)
{
  bool new_design = false;
  for(unsigned int i = 0; i < data.GetSize(); i++)
  {
    // set new_design -> save some time if already new
    if(!new_design && data[i].GetDesign(DesignElement::PLAIN) != space[i]) 
      new_design = true;
    data[i].SetDesign(space[i]);
  }
  if(new_design) design_id++;
  return design_id;
}  

int DesignSpace::WriteDesignToExtern(double* space) const
{
    // must be set in the constructor! might be trivial volume fraction or from file!!
    for(unsigned int i = 0; i < data.GetSize(); i++)
       space[i] = data[i].GetDesign(DesignElement::PLAIN);
    
    return design_id;
}  

void DesignSpace::WriteGradientToExtern(double* out, DesignElement::Type det,
                 DesignElement::ValueSpecifier vs, DesignElement::Access access) const
{
  // must be set in the constructor! might be trivial volume fraction or from file!!
  for(unsigned int i = 0; i < data.GetSize(); i++)
  {
    const DesignElement& de = data[i];
    
    bool relevant = det == DesignElement::DEFAULT 
                 || det == DesignElement::NO_TYPE
                 || det == de.GetType();

    out[i] = relevant ? data[i].GetValue(vs, access) : 0.0;
  }
}
  
void DesignSpace::Reset(DesignElement::Type design, DesignElement::ValueSpecifier vs)
{
  unsigned int base = FindDesign(design) * elements_;

  for(unsigned int i = 0; i < elements_; i++)
  {
    DesignElement& de = data[base + i];
    if(de.GetType() != design) continue;

    switch(vs)
    {
      case DesignElement::DESIGN: 
           de.SetDesign(0.0);
           break;
      
      case DesignElement::COST_GRADIENT: 
           de.SetObjectiveGradient(0.0);
           break;
           
      default: throw Exception("design not handled");     
    }
  }
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
 
  EXCEPTION("could not find element " << elemNum << " in our design space");
}

int DesignSpace::Find(const Elem* elem, bool throw_exception)
{
  if(elem->regionId == regionId_) return Find(elem->elemNum);
  
  // we might have surface element and it is pointing to a design element
  const SurfElem* se = dynamic_cast<const SurfElem*>(elem);

  // no chance, we are wrong
  if(se == NULL) 
  {
    if(!throw_exception) return -1;
    EXCEPTION("element " << elem->ToString() << " not in design region " << regionId_);
  }

  if(se->ptVolElem1 != NULL && se->ptVolElem1->regionId == regionId_)
    return Find(se->ptVolElem1->elemNum); 

  if(se->ptVolElem2 != NULL && se->ptVolElem2->regionId == regionId_)
    return Find(se->ptVolElem2->elemNum);     

  if(!throw_exception) return -1;
  
  EXCEPTION("element " << elem->ToString() << " has no volume element in design region");
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

  
       
