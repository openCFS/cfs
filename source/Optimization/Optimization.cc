#include <def_use_ipopt.hh>
#include <def_use_scpip.hh>

#include "Optimization/Optimization.hh"
#include "Optimization/SIMP.hh"
#include "Optimization/PiezoSIMP.hh"
#include "Optimization/DesignElement.hh"
#include "Optimization/DesignSpace.hh"
#include "Optimization/OptimalityCondition.hh"
#include "Optimization/EvaluateOnly.hh"
#include "Driver/basedriver.hh"
#include "Driver/singleDriver.hh"
#include "PDE/StdPDE.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/resultHandler.hh"
#include "General/exception.hh"

// IPOPT and SCPIP are not necessarily linked
#ifdef USE_IPOPT
  #include "Optimization/IPOPTHolder.hh"
#endif
#ifdef USE_SCPIP
  #include "Optimization/SCPIP.hh"
#endif

using namespace CoupledField;

DECLARE_LOG(opt)
DEFINE_LOG(opt, "opt")

// instantiation of the static elements
Enum<Optimization::ObjectiveType>    Optimization::objectiveType;
Enum<Optimization::Optimizer>        Optimization::optimizer;
Enum<Optimization::Application>      Optimization::application;


Optimization::Objective::Objective(ParamNode* pn)
{
  // the current value -> check <Get/Set>Value() when altering the presets!
  this->value_ = -1.0;
  
  type = objectiveType.Parse(pn->Get("type"));

  if(pn->Get("stoppingRule")->AsString() != "relative")
      throw Exception("stopping rule not implemented yet");

  task = pn->Get("task")->AsString() == "minimize" ? MINIMIZE : MAXIMIZE; 

  stop.value = pn->Get("value")->AsDouble();
    
  stop.queue = pn->Get("queue")->AsUInt();  
  if(stop.queue < 1) throw Exception("minimal queue value for stopping rule is 1");
}
  
double Optimization::Objective::GetValue()
{
  return value_;
}
 
void Optimization::Objective::SetValue(double val)
{
  value_ = val;   
}  

  

Optimization::Optimization()
{
  this->logFile_ = NULL;
  this->design = NULL;
  this->baseOptimizer_ = NULL;
  this->harmonic = false;
  this->currentIteration = 1; // a 1 or 0 can make a lot of difference! 0 is initial design!
  this->problemSolvedCounter = 0;
  this->problemWithinIteration = 0;
  this->logFileHeader = "iteration\tcost\tchange\tproblems"; // constraints to be added later

  // inject the driver and tell him that we do optimization
  BaseDriver* driver = domain->GetDriver();
  if(driver->GetDriverClass() != BaseDriver::SINGLE_DRIVER) 
    throw Exception("optimization not implemented for driver " + driver->GetDriverClass());

  ParamNode* pn = param->Get("optimization");       

  // the tool to solve the optimization problem 
  optimizer_ = optimizer.Parse(pn->Get("optimizer")->Get("type"));
  maxIterations = pn->Get("optimizer")->Get("maxIterations")->AsInt();

  // the cost function is mandatory     
  cost = new Objective(pn->Get("costFunction"));
  
  // the constraints are optional and might not be real constraints!
  StdVector<ParamNode*> list = pn->GetList("constraint");
  for(unsigned int i = 0; i < list.GetSize(); i++)
  {
     // the index is the current constraints size, it such works also if there are outputs
     Condition g(list[i], constraints.GetSize());
     if(g.active) constraints.Push_back(g);
             else outputs.Push_back(g);
  }           
  // first the constraints then the outputs!
  for(unsigned int i = 0; i < constraints.GetSize(); i++)
    this->logFileHeader += "\t" + constraints[i].ToString();
  for(unsigned int i = 0; i < outputs.GetSize(); i++)
    this->logFileHeader += "\t" + outputs[i].ToString();


  // make a log file for gnuplot?
  if(pn->Has("log")) 
  {
      logFile_ = new std::ofstream(pn->Get("log")->AsString().c_str());
      if(logFile_ == NULL) 
        throw Exception("cannot open log file " + pn->Get("log")->AsString() + " for writing");
  }
}

Optimization::~Optimization()
{
    // if write to file close it
    if(logFile_ != NULL)
    {
        logFile_->close();
        delete logFile_;
        logFile_ = NULL;
    }
  
    if(cost != NULL) { delete cost; cost = NULL; }
    if(design != NULL) { delete design; design = NULL; }
    if(baseOptimizer_ != NULL) { delete baseOptimizer_; baseOptimizer_ = NULL; }
    
}  

void Optimization::PostInit()
{
  ParamNode* opt = param->Get("optimization")->Get("optimizer");
  
  switch(optimizer_)
  {
    case IPOPT_SOLVER:
         #ifdef USE_IPOPT  
           baseOptimizer_ = new IPOPTHolder(this, opt);
         #else
           throw Exception("CFS++ was compiled w/o IPOPT!");
         #endif
         break;
         
    case SCPIP_SOLVER:
         #ifdef USE_SCPIP
           baseOptimizer_ = new SCPIP(this, opt);
         #else
           throw Exception("CFS++ was compiled w/o SCPIP");
         #endif
         break;  
    
    case OPTIMALITY_CONDITION:
         baseOptimizer_ = new OptimalityCondition(this, opt);
         break; 
         
    case EVALUATE_INITIAL_DESIGN:
         baseOptimizer_ = new EvaluateOnly(this, opt);
         break;
         
    default: throw Exception("optimizer not implemented");     
  }
  // add plot logging of the optimizer
  this->logFileHeader += baseOptimizer_->LogFileHeader();
}


void Optimization::SetEnums()
{
  objectiveType.SetName("Optimization::ObjectiveType");
  objectiveType.Add(COMPLIANCE, "compliance");
  objectiveType.Add(OUTPUT, "output");
  objectiveType.Add(CONJUGATE_OUTPUT, "conjugateOutput");
  objectiveType.Add(GLOBAL_DYNAMIC_COMPLIANCE, "globalDynamicCompliance");
  objectiveType.Add(RADIATION, "radiation");
  
  Condition::name.SetName("Contraint::Name");
  Condition::name.Add(Condition::VOLUME, "volume");
  Condition::name.Add(Condition::GREYNESS, "greyness");
  Condition::name.Add(Condition::GAUSS_GREYNESS, "gaussGreyness");  

  Condition::type.SetName("Contraint::Type");
  Condition::type.Add(Condition::EQUAL, "equal");
  Condition::type.Add(Condition::LOWER_BOUND, "lowerBound");
  Condition::type.Add(Condition::UPPER_BOUND, "upperBound");
  
  optimizer.SetName("Optimization::Optimizer");
  optimizer.Add(OPTIMALITY_CONDITION, "optimalityCondition");
  optimizer.Add(IPOPT_SOLVER, "ipopt");
  optimizer.Add(SCPIP_SOLVER, "scpip");  
  optimizer.Add(EVALUATE_INITIAL_DESIGN, "evaluateInitialDesign");  

  ErsatzMaterial::method.SetName("ErsatzMaterial::Method");
  ErsatzMaterial::method.Add(ErsatzMaterial::SIMP_METHOD, "simp");
  ErsatzMaterial::method.Add(ErsatzMaterial::FREE_MAT, "freeMat");
  ErsatzMaterial::method.Add(ErsatzMaterial::SHAPE_GRAD, "shapeGrad");
  
  SIMP::system.SetName("SIMP::System");
  SIMP::system.Add(SIMP::PIEZO, "piezo");
  SIMP::system.Add(SIMP::MECHANIC, "mechanic");
  
  application.SetName("Optimization::Application");
  application.Add(NO_APP, "no_app");
  application.Add(MECH, "mech");
  application.Add(MASS, "mass");
  application.Add(ELEC, "elec");
  application.Add(PIEZO_COUPLING, "piezoCoupling");
  application.Add(PRESSURE, "pressure");  
  application.Add(CHARGE_DENSITY, "chargeDensity");
  application.Add(SURFACE_NORMAL, "surfaceNormal");
}    



bool Optimization::IsMinimumReached()
{
    // this currently only implements relative stopping rule

    // we need a minimum number of itations to be sure we are in a minimum
    if(cost->history.GetSize() <= cost->stop.queue) return false;
    
    for(unsigned int i = cost->history.GetSize()-1; i >= (cost->history.GetSize() - cost->stop.queue); i--) 
    {
        double delta = std::abs(cost->history[i] - cost->history[i-1]);
        double rel = std::abs(delta / cost->history[i]);
        if(rel > cost->stop.value) return false; 
    }  
    
    // the relative values for the whole queue are smaller than the requirement -> we are done! :)
    return true;
}
 

/** read only the very basic stuff */
Optimization* Optimization::CreateInstance()
{
  // set the enums we need
  Optimization::SetEnums(); // sets also ErsatzMaterial::Method
  DesignElement::SetEnums();

  if(!param->Has("optimization")) return NULL;
  
  // we assume ersatz material, currently there is nothing else.
  // note, we read method again in the ersat material constructor.
  ParamNode* em = param->Get("optimization")->Get("ersatzMaterial"); 
  ErsatzMaterial::Method method = ErsatzMaterial::method.Parse(em->Get("method"));
  switch(method)
  {
  case ErsatzMaterial::SIMP_METHOD:
  {
    // which simp type?
    SIMP::System system = SIMP::system.Parse(em->Get("SIMP")->Get("system"));
    SIMP* simp = (system == SIMP::MECHANIC) ? new SIMP(): new PiezoSIMP();
    simp->PostInit();
    return simp;
  }
  // FreeMat, ShapeGrad, ...
  default: throw Exception("Optimization not implemented");
  }
} 

void Optimization::SolveProblem()
{
  // one driver is one multisequence step. We do this stuff here
  // and call the driver->StoreResults() multiple times f
  
  ResultHandler* rh = domain->GetResultHandler();
  unsigned int mss = domain->GetDriver()->GetActSequenceStep();
  rh->BeginMultiSequenceStep(mss, BasePDE::TRANSIENT, 9999); // max steps is high
  baseOptimizer_->SolveProblem();
  rh->FinishMultiSequenceStep();
  rh->Finalize();
}

std::string Optimization::GetSolveComment()
{
  std::ostringstream os;
  os << "iter"<< currentIteration;
  if(problemWithinIteration > 0) os << "_cnt" << problemWithinIteration;
  return os.str();
}

void Optimization::SolveStateProblem()
{
    BaseDriver* driver = domain->GetDriver();

    // we solve, give a part of the filename in case we use export linear system
    // but do not store the results. This is to be done in CommitIteration
    driver->SolveProblem(false, GetSolveComment());
      
    problemSolvedCounter++;  
    problemWithinIteration++;
} 


double Optimization::CalcSymmetry(DesignElement::Type de, DesignElement::ValueSpecifier vs, DesignElement::Access access)
{
  // the symmetry works only for squared models with a horizontal symmetry axis
  
  // our special result index
  int res_idx = design->GetSpecialResultIndex(de, vs, DesignElement::SYMMETRY, access);
  
  // plausibility check for squared
  int edge = (int) std::sqrt(design->data.GetSize());
  if(edge * edge != (int) design->data.GetSize()) 
    throw Exception("Symmetry not possible as models seems to be not squared");
  
  int max = (int) design->data.GetSize();
  double sum = 0;
  // we assume the first element (1) be in the lower left corner and then a lexicaly
  // ordering to right and then the rows up.
  for(int i = 0; i < edge/2; i++)
  {
    for(int j = 0; j < edge; j++)
    {
      // our data index
      int idx = i*edge + j;
      double idx_val = design->data[idx].GetValue(vs, access);
      
      // our counterpart
      int cntr = max - (i*edge)-(edge-j);
      double cntr_val = design->data[cntr].GetValue(vs, access);
      
      double err = (idx_val-cntr_val) / idx_val;
      LOG_DBG3(opt) << "Symmetry " << DesignElement::valueSpecifier.ToString(vs) 
                    << "(" << DesignElement::access.ToString(access) << "): "
                    << idx << " (" << idx_val << ") : " << cntr << " (" << cntr_val << ") -> " << err;

      if(res_idx >= 0) design->data[idx].specialResult[res_idx] = err;
      sum += err;
    }
  }
  
  return sum / (double) max;
}

void Optimization::EvaluateSpecialResults()
{
  for(unsigned int i = 0; i < design->resultDescriptions.GetSize(); i++)
  {
    const ResultDescription& rd = design->resultDescriptions[i];

    if(rd.detail == DesignElement::SYMMETRY)
      CalcSymmetry(rd.design, rd.value, rd.access);
  }
}

void Optimization::StoreResults(double step_val)
{
  // For PiezoSIMP we can do storing there and this method is overwritten
  // and might do nothing
  
  // this will write the CFS result and history file
  domain->GetDriver()->StoreResults(step_val == -1 ? currentIteration : step_val);
}

void Optimization::CommitIteration()
{
  LOG_TRACE2(opt) << "CommitIteration " << currentIteration << " ojective=" << cost->GetValue();

  // store the real cost -> not a scaled one
  cost->history.Push_back(cost->GetValue());

  // eventually set special result
  EvaluateSpecialResults();

  // this writes the most current solved forward problem
  // via the driver to gid or whatever  
  StoreResults();

  // save this iteration
  design->WriteDesignToExtern(last_iteration.GetPointer());

  currentIteration++;
  problemWithinIteration = 0;

  if(logFile_)
  {
    // write the header only once. Note that we start with one!
    if(currentIteration == 2) *logFile_ << logFileHeader << std::endl; 
    LogFileLine(logFile_);
    baseOptimizer_->LogFileLine(logFile_);
    *logFile_ << std::endl;
  } 

  // IPOPT does own logging -> otherwise show the user we are alive    
  if(optimizer_ != IPOPT_SOLVER)
    std::cout << "iteration " << (currentIteration-1) << " -> cost = " 
    << cost->history.Last() << std::endl;
}

void Optimization::LogFileLine(std::ofstream* out)
{
  // calculate the relative cost change
  double change = cost->history.GetSize() < 2 ? 0.0 : (cost->history.Last() 
         - cost->history[cost->history.GetSize() - 2]) / cost->history.Last();
  
  *out << (currentIteration-1) << "\t" << cost->history.Last() << "\t" << change
       << "\t" << problemSolvedCounter;
       
  for(unsigned int i = 0; i < constraints.GetSize(); i++)    
    *out << "\t" << CalcConstraint(&constraints[i]);

  for(unsigned int i = 0; i < outputs.GetSize(); i++)    
    *out << "\t" << CalcConstraint(&outputs[i]);
  
  out->flush();
}


Condition& Optimization::GetConstraint(Condition::Name name, DesignElement::Type design)
{
  // be save and check for uniqueness! 
  int count = 0;
  
  for(unsigned int i = 0; i < constraints.GetSize(); i++)
    if(constraints[i].GetName() == name &&
       (design != DesignElement::NO_TYPE ? constraints[i].design == design : true)) 
          count++;

  for(unsigned int i = 0; i < outputs.GetSize(); i++)
    if(outputs[i].GetName() == name &&
       (design != DesignElement::NO_TYPE ? outputs[i].design == design : true)) 
          count++;

  if(count > 1) 
    throw Exception("constraint " + Condition::name.ToString(name) + "is not unique"); 
  
  for(unsigned int i = 0; i < constraints.GetSize(); i++)
    if(constraints[i].GetName() == name) return constraints[i];

  for(unsigned int i = 0; i < outputs.GetSize(); i++)
    if(outputs[i].GetName() == name) return outputs[i];
    
  throw Exception("no constraint " + Condition::name.ToString(name) + " found"); 
}


