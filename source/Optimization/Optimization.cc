#include "Optimization/Optimization.hh"
#include "Optimization/SIMP.hh"
#include "Optimization/PiezoSIMP.hh"
#include "Optimization/DesignElement.hh"
#include "Optimization/DesignSpace.hh"
#include "Driver/basedriver.hh"
#include "Driver/singleDriver.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/exception.hh"


using namespace CoupledField;


// instantiation of the static elements
Enum Optimization::optimizationType;
Enum Optimization::objectiveType;
Enum Optimization::optimizerType;
Enum Optimization::application;
Enum PiezoSIMP::storage;


Optimization::Objective::Objective(ParamNode* pn)
{
  // the current value -> check <Get/Set>Value() when altering the presets!
  this->value_ = -1.0;
  
  type = static_cast<ObjectiveType>(objectiveType.Parse(pn->Get("type")->AsString()));

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
  this->scpip_ = NULL;
  this->logFile_ = NULL;
  this->design = NULL;
  this->currentIteration = 1;
  this->problemSolvedCounter = 0;
  this->problemWithinIteration = 0;
  this->logFileHeader = "iteration\tcost\tchange\tproblems"; // constraints to be added later

  // inject the driver and tell him that we do optimization
  BaseDriver* driver = domain->GetDriver();
  if(driver->GetDriverClass() != BaseDriver::SINGLE_DRIVER) 
    throw Exception("optimization not implemented for driver " + driver->GetDriverClass());
  driver->SetOptimization(true);

  ParamNode* pn = param->Get("optimization");       

  // the optimization problem
  optimization = (OptimizationType) optimizationType.Parse(pn->Get("type")->AsString());
 
  // the tool to solve the optimization problem 
  optimizer = (OptimizerType) optimizerType.Parse(pn->Get("optimizer")->Get("type")->AsString());
  maxIterations = pn->Get("optimizer")->Get("maxIterations")->AsInt();

  // if optimalityCondition is in XML we read and don't care if this is 
  // our optimizer
  if(pn->Has("optimalityCondition")) {
    this->move_limit = pn->Get("optimalityCondition")->Get("move_limit")->AsDouble();
    this->oc_damping = pn->Get("optimalityCondition")->Get("damping")->AsDouble();
  } else {
    // the following values are standard in mech SIMP -> see e.g. the 99 lines paper
    this->move_limit = 0.2;
    this->oc_damping = 0.5;
  }
  
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
    
    if(scpip_ != NULL) { delete scpip_; scpip_ = NULL; }
}  

void Optimization::PostInit()
{
   // construct IPOPT
   if(optimizer == IPOPT_SOLVER) 
   {
     #ifdef USE_IPOPT  
       ipopt_ = new IPOPT(this, param->Get("optimization")->Get("optimizer")->Get("ipopt", false));
     #else
       throw Exception("CFS++ was compiled w/o IPOPT!");
     #endif      
   }
   if(optimizer == SCPIP_SOLVER)
   {
     #ifdef USE_SCPIP
       scpip_ = new SCPIP(this, param->Get("optimization")->Get("optimizer")->Get("scpip", false));
     #else
       throw Exception("CFS++ was compiled w/o SCPIP");
     #endif  
   }
}

void Optimization::SetEnums()
{
  optimizationType.SetName("Optimization::OptimizationType");
  optimizationType.Add(SIMP_TYPE, "SIMP");
  
  objectiveType.SetName("Optimization::ObjectiveType");
  objectiveType.Add(COMPLIANCE, "compliance");
  objectiveType.Add(TRANSDUCTION, "transduction");
  
  Condition::name.SetName("Contraint::Name");
  Condition::name.Add(Condition::VOLUME, "volume");
  Condition::name.Add(Condition::GREYNESS, "greyness");
  Condition::name.Add(Condition::GAUSS_GREYNESS, "gaussGreyness");  

  Condition::type.SetName("Contraint::Type");
  Condition::type.Add(Condition::EQUAL, "equal");
  Condition::type.Add(Condition::LOWER_BOUND, "lowerBound");
  Condition::type.Add(Condition::UPPER_BOUND, "upperBound");
  
  optimizerType.SetName("Optimization::OptimizerType");
  optimizerType.Add(OPTIMALITY_CONDITION, "optimalityCondition");
  optimizerType.Add(IPOPT_SOLVER, "ipopt");
  optimizerType.Add(SCPIP_SOLVER, "scpip");  
  optimizerType.Add(EVALUATE_INITIAL_DESIGN, "evaluateInitialDesign");  
  
  SIMP::system.SetName("SIMP::System");
  SIMP::system.Add(SIMP::PIEZO, "piezo");
  SIMP::system.Add(SIMP::MECHANIC, "mechanic");
  
  application.SetName("Optimization::Application");
  application.Add(MECH, "mech");
  application.Add(ELEC, "elec");
  application.Add(PIEZO_COUPLING, "piezoCoupling");
  application.Add(PRESSURE, "pressure");  
  application.Add(CHARGE_DENSITY, "chargeDensity");  
  
  
  PiezoSIMP::storage.SetName("PiezoSIMP::Storage");
  PiezoSIMP::storage.Add(PiezoSIMP::COMMIT, "commit");
  PiezoSIMP::storage.Add(PiezoSIMP::ELEC_CASE, "elec_case");
  PiezoSIMP::storage.Add(PiezoSIMP::MECH_CASE, "mech_case");
  PiezoSIMP::storage.Add(PiezoSIMP::BOTH_CASES, "both_cases");
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
  Optimization::SetEnums();
  DesignElement::SetEnums();

  if(!param->Has("optimization")) return NULL;

  std::string string = param->Get("optimization")->Get("type")->AsString();

  // the actual parameters are read in the constructors
  switch(static_cast<OptimizationType>(optimizationType.Parse(string)))
  {
    case SIMP_TYPE: 
    {
      string = param->Get("optimization")->Get("SIMP")->Get("system")->AsString();
      // determine subtype - this is again done in the constructor
      SIMP::System system = (SIMP::System) SIMP::system.Parse(string);
      SIMP* simp = (system == SIMP::MECHANIC) ? new SIMP(): new PiezoSIMP();
      simp->PostInit();
      return simp;
    }   
    
    default: throw Exception("Optimization " + string + " not implemented");
  }
} 

void Optimization::SolveProblem()
{
   switch(optimizer)
   {
      case OPTIMALITY_CONDITION: 
           SolveProblemManually();
           break;
      
      case IPOPT_SOLVER:
           SolveStateProblem(); // ipopt starts with the gradient ...
           #ifdef USE_IPOPT
             ipopt_->SolveProblem();
           #else
             throw Exception("CFS++ was compiled w/o IPOPT");
           #endif
           break;

      case SCPIP_SOLVER:
           // scpip starts with the function evaluation
           #ifdef USE_SCPIP 
             scpip_->SolveProblem();
           #else
             throw Exception("CFS++ was compiled w/o SCPIP");
           #endif  
           break;
           
      case EVALUATE_INITIAL_DESIGN:
           // no iterations but only the evaluation of the initial guesses.
           EvaluateInitialDesign();
           break;     

           
      default: throw Exception("optimizer not implemented");     
   }
}


void Optimization::EvaluateInitialDesign()
{
  // solve the state problem with the initial guess.
  std::cout << "Evaluate state problem for initial guess ..." << std::endl;
  // note that when PiezoSIMP and storage is not PiezoSIMP::COMMIT we might have wrong 
  // special results before CalcObjective()
  SolveStateProblem();
  std::cout << "objective: " << CalcObjective() << std::endl;
  // see note above!
  if(optimization == SIMP_TYPE 
     && (dynamic_cast<SIMP*>(this))->GetSystem() == SIMP::PIEZO
     && (dynamic_cast<PiezoSIMP*>(this))->GetStorage() != PiezoSIMP::COMMIT)
    SolveStateProblem();
  
  
  // calc gradients, they might be sored in store results!
  CalcObjectiveGradient(NULL);

  for(unsigned int i = 0; i < constraints.GetSize(); i++)
  { 
     std::cout << "constraint " << constraints[i].ToString() << ": " 
               << CalcConstraint(&constraints[i]) << std::endl;
     CalcConstraintGradient(&constraints[i], NULL);          
  }

  for(unsigned int i = 0; i < outputs.GetSize(); i++)
  { 
     std::cout << "output_only " << outputs[i].ToString() << ": "
               << CalcConstraint(&outputs[i]) << std::endl;
  }
  CommitIteration();
}


void Optimization::SolveStateProblem()
{
    BaseDriver* driver = domain->GetDriver();

    // we alter the driver in a way that SolveProblem() can behave unique in the first
    // call - if we are a single driver!
    if(problemSolvedCounter > 0) {
       SingleDriver* sd = dynamic_cast<SingleDriver*>(driver);
       sd->SetConsecutiveRun(true);
    }

    driver->SolveProblem(currentIteration);
      
    problemSolvedCounter++;  
    problemWithinIteration++;
} 

bool Optimization::NeedObjectiveEval(const double* space_in)
{
   if(design->data.GetSize() != last_evaluation.GetSize())
     throw Exception("index mixed up");
     
   for(unsigned int i = 0; i < last_evaluation.GetSize(); i++)
     if(space_in[i] != last_evaluation[i]) return true;

   return false; // same same   
} 

void Optimization::StoreResults(double step_val)
{
  // For PiezoSIMP we can do storing there and this method is overwritten
  // and might do nothing
  
  // this will write the CFS result and history file
  domain->GetDriver()->StoreResults(step_val);
}

void Optimization::CommitIteration()
{
    // store the real cost -> not a scaled one
    cost->history.Push_back(cost->GetValue());

    // this writes the most corrent solved forward problem
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
        *logFile_ << std::endl;
      } 
 
      // IPOPT does own logging -> otherwise show the user we are alive    
    if(optimizer != IPOPT_SOLVER)
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


void Optimization::SolveProblemManually()
{
  bool first = true; 
  
  while(!IsMinimumReached() && currentIteration <= maxIterations)
  {
      // adjust the design parameters but not if first iteration 
      if(!first) ExecuteOptimizationStep();

      // solve the state problem
      SolveStateProblem();
      
      // every state problem is an iteration
      if(!first) CommitIteration();
      first = false; 
  }
  
  if(currentIteration >= maxIterations-1) {
     std::cout << " max iterations reached" << std::endl;
  }
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


