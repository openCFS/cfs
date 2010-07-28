#include <def_use_ipopt.hh>
#include <def_use_scpip.hh>
#include <def_use_snopt.hh>
#include <map>

#include "Optimization/Optimization.hh"
#include "Optimization/SIMP.hh"
#include "Optimization/PiezoSIMP.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Design/DesignStructure.hh"
#include "Optimization/Optimizer/OptimalityCondition.hh"
#include "Optimization/LevelSet.hh"
#include "Optimization/Optimizer/EvaluateOnly.hh"
#include "Optimization/Optimizer/ShapeOptimizer.hh"
#include "Optimization/ShapeGrad.hh"
#include "Optimization/Optimizer/GradientCheck.hh"
#include "Optimization/ParamMat.hh"
#include "Optimization/OptimizationMaterial.hh"
#include "Driver/assemble.hh"
#include "Driver/basedriver.hh"
#include "Driver/harmonicDriver.hh"
#include "Driver/singleDriver.hh"
#include "Driver/stdSolveStep.hh"
#include "PDE/StdPDE.hh"
#include "PDE/SinglePDE.hh"
#include "PDE/mechPDE.hh"
#include "PDE/elecPDE.hh" // for polarization matrix, see class TopGrad
#include "Domain/domain.hh"
#include "Domain/grid.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/resultHandler.hh"
#include "DataInOut/programOptions.hh"
#include "General/exception.hh"
#include <boost/filesystem.hpp>
#include "Optimization/ShapeOpt.hh"

// IPOPT, SCPIP and SnOpt are not necessarily linked
#ifdef USE_IPOPT
  #include "Optimization/Optimizer/IPOPTHolder.hh"
#endif
#ifdef USE_SCPIP
  #include "Optimization/Optimizer/SCPIP.hh"
#endif
#ifdef USE_SNOPT
  #include "Optimization/Optimizer/SnOpt.hh"
#endif

using namespace CoupledField;
using namespace std;
namespace fs = boost::filesystem;



DECLARE_LOG(opt)
DEFINE_LOG(opt, "opt")

// instantiation of the static elements
Enum<Optimization::Optimizer>        Optimization::optimizer;
Enum<Optimization::Application>      Optimization::application;
Enum<Optimization::CommitMode>       Optimization::commitMode;

Enum<Optimization::MultipleExcitation::Type>  Optimization::MultipleExcitation::type;

Optimization::Optimization()
{
  this->lastStoredResult_ = -1;
  this->design = NULL;
  this->baseOptimizer_ = NULL;
  this->harmonic = BasePDE::IsComplex(domain->GetDriver()->GetAnalysisType());
  this->currentIteration = 0; // a 1 or 0 can make a lot of difference! 0 is initial design!
  this->problemSolvedCounter = 0;
  this->problemWithinIteration = 0;
  this->grid = domain->GetGrid();

  // inject the driver and tell him that we do optimization
  BaseDriver* driver = domain->GetDriver();
  if(driver->GetDriverClass() != BaseDriver::SINGLE_DRIVER)
    throw Exception("optimization not implemented for driver " + driver->GetDriverClass());

  optInfoNode = info->Get("optimization");   // store our info results here
  PtrParamNode pn = param->Get("optimization"); // read our parameters from the xml file
  
  // in transient optimization one can specify the initial value as a solution to a static problem and a weight for it (just in tracking)
  firstStepStatic = pn->Has("firstStepStatic");
  if(firstStepStatic){
    otherStepWeight = 1.0 - pn->Get("firstStepStatic")->Get("weight")->As<Double>();
  }else{
    otherStepWeight = 1.0;
  }

  // the tool to solve the optimization problem
  optimizer_ = optimizer.Parse(pn->Get("optimizer")->Get("type")->As<std::string>());
  maxIterations = pn->Get("optimizer")->Get("maxIterations")->As<Integer>();

  // might read a multiObjective problem
  objectives.Read(pn->Get("costFunction"));
  objectives.ToInfo(optInfoNode->Get(ParamNode::HEADER)->Get("objective"));

  // constraints to be added later -- it is so much easier with the ParamNodes
  this->log.fileHeader = harmonic ? "#iter\tfreq" : "#iter";
  for(unsigned int i = 0; i < objectives.data.GetSize(); i++)
    this->log.fileHeader += "\t" + objectives.data[i]->GetName();
  this->log.fileHeader += "\tchange\tproblems";

  // multiple excitations are are toggled via attribute. Only if enabled we read the optional element
  // actually part of costFunction - but we store in Optimization itself!
  bool me = pn->Get("costFunction")->Get("multiple_excitation")->As<bool>();
  this->multiple_excitation = new MultipleExcitation(me, me ? pn->Get("costFunction")->Get("multipleExcitation", ParamNode::PASS) : PtrParamNode());
  if(me) this->multiple_excitation->ToInfo(optInfoNode->Get(ParamNode::HEADER)->Get("multipleExcitations"));

  // slope constraints to be processed in SIMP -> Constraints::PostProc
  ParamNodeList list = pn->GetList("constraint");
  constraints.Read(list);
  PtrParamNode in = optInfoNode->Get(ParamNode::HEADER)->Get("constraints");

  // call this again after PostProc()
  constraints.ToInfo(optInfoNode->Get(ParamNode::HEADER)->Get("constraints"));

  for(unsigned int i = 0; i < constraints.all.GetSize(); i++)
  {
    Condition* g = constraints.all[i];
    if(!g->IsLocalCondition())
      this->log.fileHeader += "\t" + g->ToString();
    else
    {
      this->log.fileHeader += "\tmax_" + g->ToString();
      this->log.fileHeader += "\tmean_" + g->ToString();
    }
  }

  log.Init(pn->Get("log")->As<std::string>(), pn->Get("logging", ParamNode::PASS)); // is fail save

  // the commit stuff
  string cm = pn->Has("commit") ? pn->Get("commit")->Get("mode")->As<std::string>() : "forward";
  this->commitMode_ = commitMode.Parse(cm);
  this->commitStride = pn->Has("commit") ? pn->Get("commit")->Get("stride")->As<Integer>() : 1;
  optInfoNode->Get("commit")->Get("mode")->SetValue(cm);
  optInfoNode->Get("commit")->Get("stride")->SetValue(commitStride);
  
  // write out the directory where the HALTOPT file will be searched for
  optInfoNode->Get("haltopt_directory")->SetValue(fs::current_path().directory_string());

  // remove a stop file, if found
  if(fs::exists("HALTOPT"))
  {
    bool good = fs::remove("HALTOPT");
    if(!good) throw new Exception("Could not remove file 'HALTOPT' after detection");
  }
}

Optimization::~Optimization()
{
  delete design; design = NULL;
  delete baseOptimizer_; baseOptimizer_ = NULL;
  delete multiple_excitation; multiple_excitation = NULL;
}


void Optimization::PostInitSecond()
{
  PtrParamNode opt = param->Get("optimization")->Get("optimizer");

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
         
    case SNOPT_SOLVER:
         #ifdef USE_SNOPT
           baseOptimizer_ = new SnOpt(this, opt);
         #else
           throw Exception("CFS++ was compiled w/o SnOpt");
         #endif
      break;

    case OPTIMALITY_CONDITION:
         baseOptimizer_ = new OptimalityCondition(this, opt);
         break;

    case SHAPE_SOLVER:
         baseOptimizer_ = new ShapeOptimizer(this, opt);
         break;

    case EVALUATE_INITIAL_DESIGN:
         baseOptimizer_ = new EvaluateOnly(this, opt);
         break;
         
    case GRADIENT_CHECK:
         baseOptimizer_ = new GradientCheck(this, opt);
         break;

    default: throw Exception("optimizer not implemented");
  }
  if (this->log.design) {
    for (unsigned int i = 0; i < design->GetNumberOfVariables(); i++) {
      this->log.fileHeader += "\t";
      this->log.fileHeader += "design";
    }
  }
  if (this->log.designGradient) {
    for (unsigned int i = 0; i < design->GetNumberOfVariables(); i++) {
      this->log.fileHeader += "\t";
      this->log.fileHeader += "designGradient";
    }
  }
  design->SetOptimizer(baseOptimizer_);
  // add plot logging of the optimizer
  this->log.fileHeader += baseOptimizer_->LogFileHeader();
}


void Optimization::SetEnums()
{
  Function::type.SetName("Function::Type");
  Function::type.Add(Function::MULTI_OBJECTIVE, "multiObjective");
  Function::type.Add(Function::COMPLIANCE, "compliance");
  Function::type.Add(Function::OUTPUT, "output");
  Function::type.Add(Function::DYNAMIC_OUTPUT, "dynamicOutput");
  Function::type.Add(Function::ABS_DYN_OUTPUT_SQUARED, "absDynamicOutputSquared");
  Function::type.Add(Function::GLOBAL_DYNAMIC_COMPLIANCE, "globalDynamicCompliance");
  Function::type.Add(Function::CONJUGATE_COMPLIANCE, "conjugateCompliance");
  Function::type.Add(Function::VOLUME, "volume");
  Function::type.Add(Function::PENALIZED_VOLUME, "penalizedVolume");
  Function::type.Add(Function::GAP, "gap");
  Function::type.Add(Function::REALVOLUME, "realvolume");
  Function::type.Add(Function::TRACKING, "tracking");
  Function::type.Add(Function::ELEC_ENERGY, "elecEnergy");
  Function::type.Add(Function::ENERGY_FLUX, "energyFlux");
  Function::type.Add(Function::HOMOGENIZATION_TENSOR, "homTensor");
  Function::type.Add(Function::HOMOGENIZATION_TRACKING, "homTracking");
  Function::type.Add(Function::POISSONS_RATIO, "poissonsRatio");
  Function::type.Add(Function::YOUNGS_MODULUS, "youngsModulus");
  Function::type.Add(Function::TYCHONOFF, "tychonoff");
  Function::type.Add(Function::TEMPERATURE, "temperature");
  Function::type.Add(Function::GREYNESS, "greyness");
  Function::type.Add(Function::ISOTROPY, "isotropy");
  Function::type.Add(Function::SLOPE, "slope");
  Function::type.Add(Function::GLOBAL_SLOPE, "globalSlope");
  Function::type.Add(Function::CHECKERBOARD, "checkerboard");
  Function::type.Add(Function::GLOBAL_CHECKERBOARD, "globalCheckerboard");
  Function::type.Add(Function::MOLE, "mole");
  Function::type.Add(Function::GLOBAL_MOLE, "globalMole");
  Function::type.Add(Function::OSCILLATION, "oscillation");
  Function::type.Add(Function::GLOBAL_OSCILLATION, "globalOscillation");

  Function::Local::locality.SetName("Function::Local::Locality");
  Function::Local::locality.Add(Function::Local::DEFAULT, "default");
  Function::Local::locality.Add(Function::Local::NEXT, "next");
  Function::Local::locality.Add(Function::Local::NEXT_AND_REVERSE, "next_and_reverse");
  Function::Local::locality.Add(Function::Local::PREV_NEXT_AND_REVERSE, "prev_next_and_reverse");
  Function::Local::locality.Add(Function::Local::DEG_45_STAR, "45_deg_star");
  Function::Local::locality.Add(Function::Local::DEG_45_STAR_AND_REVERSE, "45_deg_star_and_reverse");

  Condition::bound.SetName("Condition::Bound");
  Condition::bound.Add(Condition::EQUAL, "equal");
  Condition::bound.Add(Condition::LOWER_BOUND, "lowerBound");
  Condition::bound.Add(Condition::UPPER_BOUND, "upperBound");

  DesignStructure::filterSpace.SetName("DesignStructure::FilterSpace");
  DesignStructure::filterSpace.Add(DesignStructure::RADIUS, "radius");
  DesignStructure::filterSpace.Add(DesignStructure::VOLUME_RADIUS, "volumeRadius");
  DesignStructure::filterSpace.Add(DesignStructure::MAX_EDGE, "maxEdge");

  optimizer.SetName("Optimization::Optimizer");
  optimizer.Add(OPTIMALITY_CONDITION, "optimalityCondition");
  optimizer.Add(IPOPT_SOLVER, "ipopt");
  optimizer.Add(SCPIP_SOLVER, "scpip");
  optimizer.Add(SNOPT_SOLVER, "snopt");
  optimizer.Add(SHAPE_SOLVER, "shapeOpt");
  optimizer.Add(EVALUATE_INITIAL_DESIGN, "evaluateInitialDesign");
  optimizer.Add(GRADIENT_CHECK, "gradientCheck");

  ErsatzMaterial::method.SetName("ErsatzMaterial::Method");
  ErsatzMaterial::method.Add(ErsatzMaterial::SIMP_METHOD, "simp");
  ErsatzMaterial::method.Add(ErsatzMaterial::PARAM_MAT, "paramMat");
  ErsatzMaterial::method.Add(ErsatzMaterial::SHAPE_GRAD, "shapeGrad");
  ErsatzMaterial::method.Add(ErsatzMaterial::SHAPE_OPT, "shapeOpt");
  ErsatzMaterial::method.Add(ErsatzMaterial::SHAPE_PARAM_MAT, "shapeParamMat");
  
  ErsatzMaterial::commitMode.SetName("ErsatzMaterial::CommitMode");
  ErsatzMaterial::commitMode.Add(ErsatzMaterial::FORWARD, "forward");
  ErsatzMaterial::commitMode.Add(ErsatzMaterial::ADJOINT, "adjoint");
  ErsatzMaterial::commitMode.Add(ErsatzMaterial::BOTH, "both_cases");

  OptimizationMaterial::system.SetName("OptimizationMaterial::System");
  OptimizationMaterial::system.Add(OptimizationMaterial::PIEZO, "piezo");
  OptimizationMaterial::system.Add(OptimizationMaterial::MECHANIC, "mechanic");
  OptimizationMaterial::system.Add(OptimizationMaterial::HEAT, "heat");

  application.SetName("Optimization::Application");
  application.Add(NO_APP, "no_app");
  application.Add(MECH, "mech");
  application.Add(MASS, "mass");
  application.Add(ELEC, "elec");
  application.Add(PIEZO_COUPLING, "piezoCoupling");
  application.Add(PRESSURE, "pressure");
  application.Add(CHARGE_DENSITY, "chargeDensity");

  LevelSet::Action::type.SetName("LevelSet::Action::Type");
  LevelSet::Action::type.Add(LevelSet::Action::SIGNED_DISTANCE_FIELD, "signedDistanceField");
  LevelSet::Action::type.Add(LevelSet::Action::TRIVIAL_HOLE, "trivialHole");
  LevelSet::Action::type.Add(LevelSet::Action::DO_SHAPE_STEP, "shapeStep");

  MultipleExcitation::type.SetName("Optimization::MultipleExcitation::Type");
  MultipleExcitation::type.Add(MultipleExcitation::NO_TYPE, "no_type");
  MultipleExcitation::type.Add(MultipleExcitation::FIXED_WEIGHT, "fixed_weights");
  MultipleExcitation::type.Add(MultipleExcitation::META_OBJECTIVE, "meta_objective");
  MultipleExcitation::type.Add(MultipleExcitation::HOMOGENIZATION_TEST_STRAINS, "homogenizationTestStrains");
  MultipleExcitation::type.Add(MultipleExcitation::POLARIZATION_MATRIX, "polarizationMatrix");
}

bool Optimization::IsTransient() const{
  return(domain->GetDriver()->GetAnalysisType() == BasePDE::TRANSIENT);
}

double Optimization::GetStepWeight(unsigned int ts) const{
  unsigned int nts = domain->GetDriver()->GetNumSteps();
  if(IsFirstTransientStepStatic()){
    if(ts == 0){
      return((1.0 - otherStepWeight));
    }else{
      return(otherStepWeight / (nts-1));
    }
  }else{
    return(1.0 / nts);
  }
}

bool Optimization::DoStopOptimization()
{
  PtrParamNode in = optInfoNode->Get(ParamNode::SUMMARY)->Get("break");
  // check if the HALTOPT file exists
  if(fs::exists("HALTOPT"))
  {
    bool good = fs::remove("HALTOPT");
    if(!good) throw new Exception("Could not remove file 'HALTOPT' after detection");
    in->Get("converged")->SetValue("no");
    in->Get("reason/msg")->SetValue("Detected file 'HALTOPT'");
    return true;
  }
  
  // this currently only implements relative stopping rule
  Objective* cost = objectives.data[0];

  // we need a minimum number of iterations to be sure we are in a minimum
  unsigned int hs = objectives.GetHistorySize();
  if(hs <= cost->stop.queue) return false;

  for(unsigned int i = hs-1; i >= (hs - cost->stop.queue); i--)
  {
    double delta = objectives.GetHistoryValue(true, i) - objectives.GetHistoryValue(true, i-1);
    double rel = abs(delta / objectives.GetHistoryValue(true, i));
    if(rel > cost->stop.value) return false;
  }

  // the relative values for the whole queue are smaller than the requirement -> we are done! :)
  in->Get("converged")->SetValue("practically");
  in->Get("reason/msg")->SetValue("Too small change in objective function");
  in->Get("reason/queue")->SetValue(cost->stop.queue);
  in->Get("reason/relative")->SetValue(cost->stop.value);
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
  // note, we read method again in the ersatz material constructor.
  PtrParamNode em = param->Get("optimization")->Get("ersatzMaterial");
  ErsatzMaterial::Method method = 
      ErsatzMaterial::method.Parse(em->Get("method")->As<std::string>());
  
  Optimization* opt = NULL;
  
  switch(method)
  {
  case ErsatzMaterial::SIMP_METHOD:
    switch(OptimizationMaterial::system.Parse(em->Get("material")->As<std::string>()))
    {
    case OptimizationMaterial::MECHANIC:
      opt = new SIMP();
      break;
      
    case OptimizationMaterial::PIEZO:
      opt = new PiezoSIMP();
      break;
      
    default: assert(false);
    }
    break;
    
  // FreeMat, ShapeGrad, ...
  case ErsatzMaterial::PARAM_MAT:
    opt = new ParamMat();
    break;
  case ErsatzMaterial::SHAPE_OPT:
  case ErsatzMaterial::SHAPE_PARAM_MAT:
    opt = new ShapeOpt();
    break;
  case ErsatzMaterial::SHAPE_GRAD:
  {
    opt = new ShapeGrad();
    break;
  }
  default: throw Exception("Optimization not implemented");
  }
  
  // we have to do this, as PostInitSecond does already run CalcObjective/Gradient
  domain->SetOptimization(opt);
  
  // second initialization phase, constructs material
  opt->PostInit();
  // third initialization phase, constructs optimizer
  opt->PostInitSecond();

  return opt;
}

void Optimization::SolveProblem()
{
  // one driver is one multisequence step. We do this stuff here
  // and call the driver->StoreResults() multiple times f
  
  ResultHandler* rh = NULL;

  if(!IsTransient()){ // transient optimization saves results in a different way
    rh = domain->GetResultHandler();
    unsigned int mss = domain->GetDriver()->GetActSequenceStep();
    rh->BeginMultiSequenceStep(mss, BasePDE::TRANSIENT, 999); // max steps is high
  }

  Exception* e = NULL;
  try
  {
    baseOptimizer_->SolveOptimizationProblem();
  }
  catch(Exception& ex)
  {
    e = new Exception(ex);
  }

  if(!IsTransient()){ // transient optimization saves results in a different way
    // do the finally - try to write results even if the optimizer broke down
    FinalizeStoreResults(); // when we have strides the results are written
    rh->FinishMultiSequenceStep();
    rh->Finalize();
  }
  if(e != NULL) throw *e;
  delete e;
}
PtrParamNode Optimization::CreateAdjointAnalysisIdNode(std::string child_name)
{
  BaseDriver* driver = domain->GetDriver();
  PtrParamNode base = driver->GetAnalysisId();
  PtrParamNode in = driver->CreateAnalysisIdChild(base, child_name);
  
  return in;
}


void Optimization::SolveStateProblem(Excitation* excite)
{
  BaseDriver* driver = domain->GetDriver();

  // this is the forward problem. Store the analysis_id in the driver such that
  // we can aquire it for the adjoint problem
  PtrParamNode analysis_id = excite == NULL ? driver->CreateAnalysisId("iter", currentIteration)
                                         : driver->CreateAnalysisId("iter", currentIteration, "excite", excite->index);
  
  if(IsTransient() && problemSolvedCounter > 0){ // transient optimization always has a mech pde
    SinglePDE* mech = domain->GetSinglePDE("mechanic");
    mech->ReReadResults();
    design->AppendOptimizationResults(mech);
    mech->GetSolveStep()->ReInit();
  }
                                         
  // Do not store the results. This is to be done in CommitIteration
  if(!harmonic || excite == NULL) 
    driver->SolveProblem(IsTransient(), analysis_id, false); // static and transient optimization
  else
    dynamic_cast<HarmonicDriver*>(driver)->ComputeFrequencyStep(excite->f_link->step, analysis_id);

  problemSolvedCounter++;
  problemWithinIteration++;
}

void Optimization::SolveAdjointProblem(Excitation* excite, Objective* cost){
  // does almost the same as SolveStateProblem now, but passing, that we want the adjoint to be solved
  BaseDriver* driver = domain->GetDriver();
  
  AdjointParameters adjointParams(cost, excite);
  
  if(IsTransient()){
    SinglePDE* mech = domain->GetSinglePDE("mechanic");
    mech->GetSolveStep()->ReInit();
  }

  // Do not store the results. This is adjoint.
  if(!harmonic) 
    driver->SolveProblem(false, CreateAdjointAnalysisIdNode(), &adjointParams); // static and transient optimization
  else
      EXCEPTION("Harmonic adjoint not implemented!");
}

void Optimization::SolveAdjointProblems(Excitation* excite){
  for(unsigned int o = 0; o < objectives.data.GetSize(); ++o){
    SolveAdjointProblem(excite, objectives.data[o]);
  }
}

double Optimization::CalcSymmetry(DesignElement::Type de, DesignElement::ValueSpecifier vs, DesignElement::Access access)
{
  // the symmetry works only for squared models with a horizontal symmetry axis

  // our special result index
  int res_idx = design->GetSpecialResultIndex(de, vs, DesignElement::SYMMETRY, access);

  // plausibility check for squared
  int edge = (int) sqrt(design->data.GetSize());
  if(edge * edge != (int) design->data.GetSize())
    throw Exception("Symmetry not possible as models seems to be not squared");

  int max = (int) design->data.GetSize();
  double sum = 0;
  // we assume the first element (1) be in the lower left corner and then a lexical
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
  if(!IsTransient()){ // transient optimization saves results in a different way
    domain->GetDriver()->StoreResults(step_val == -1 ? currentIteration : step_val);
  }
}

void Optimization::FinalizeStoreResults()
{
  // after the last CommitIteration the iteration counter was incremented
  bool store = currentIteration-1 != lastStoredResult_ && currentIteration > 1;
  LOG_DBG(opt) << "CheckFinalStoreResults: currentIteration=" << currentIteration << " lastStoredResult="
               << lastStoredResult_ << " store=" << store;
  if(store)
    StoreResults(currentIteration-1);
}


PtrParamNode Optimization::CommitIteration(bool keep_iteration_number)
{
  // store the real cost -> not a scaled one
  objectives.PushBackHistory();

  // eventually set special result
  EvaluateSpecialResults();

  // this writes the most current solved forward problem via the driver to gid or whatever
  bool store = currentIteration == 0 || commitStride == 1 || (commitStride > 0 && currentIteration % commitStride == 0);
  LOG_TRACE2(opt) << "CommitIteration " << currentIteration << " ojective=" << objectives.GetHistoryValue() << " store=" << store;
  if(store)
  {
    StoreResults();
    lastStoredResult_ = currentIteration;
    // see FinalizeStoreResults() !
  }

  // save this iteration
  design->WriteDesignToExtern(last_iteration.GetPointer());

  // also log to info node, append the iteration
  PtrParamNode iteration = optInfoNode->Get(ParamNode::PROCESS)->Get("iteration", ParamNode::APPEND);

  // write the header only once - we might keep the iteration number
  if(log.file) if(objectives.GetHistorySize() == 1) *log.file << log.fileHeader << endl;
  LogFileLine(log.file, iteration); // also ParamNode is to be written
  baseOptimizer_->LogFileLine(log.file, iteration);
  if(log.file) *log.file << endl;

  // IPOPT does own logging -> otherwise show the user we are alive
  std::string f = GetIterationFrequency();
  if(optimizer_ != IPOPT_SOLVER && optimizer_ != SNOPT_SOLVER)
  {
    cout << "iteration " << (currentIteration);
    if(f != "") cout << " f = " << f << "Hz";
    cout << " -> cost = " << objectives.GetHistoryValue() << endl;
  }

  if(!keep_iteration_number)
  {
    currentIteration++;
    problemWithinIteration = 0;
  }

  return iteration;
}

void Optimization::LogFileLine(ofstream* out, PtrParamNode iteration)
{
  // calculate the relative cost change
  int hs = objectives.GetHistorySize();
  double change = 0.0;
  if(hs >= 2)
  {
    double last = objectives.GetHistoryValue(true, hs-1);
    change = (last - objectives.GetHistoryValue(true, hs-2)) / last;
  }

  if(out)
  {
    *out << currentIteration;
    if(harmonic) *out << "\t" << GetIterationFrequency();

    for(unsigned int i = 0; i < objectives.data.GetSize(); i++)
      *out << "\t" << objectives.data[i]->GetValue();

    *out << "\t" << change << "\t" << problemSolvedCounter;
  }

  iteration->Get("number")->SetValue(currentIteration);
  if(harmonic) iteration->Get("frequency")->SetValue(GetIterationFrequency());

  for(unsigned int i = 0; i < objectives.data.GetSize(); i++)
    iteration->Get(Objective::type.ToString(objectives.data[i]->GetType()))->SetValue(objectives.data[i]->GetValue());

  iteration->Get("change")->SetValue(change);
  iteration->Get("problemsSolved")->SetValue(problemSolvedCounter);

  // For iteration 0 we want also the constraint values but they were not evaluated.
  // For any iteration we need to evaluate the observe constraints
  // A problem are the slope constraints, they need to be evaluated in local mode
  // and Done() forms the global result
  for(unsigned int i = 0, m = constraints.view->GetNumberOfTotalConstraints(); i < m; i++)
  {
    Condition* g = constraints.view->Get(i); // traverse in local mode
    if(g->GetValue() == -1.0 || g->IsObservation())
      CalcConstraint(g);
  }
  constraints.view->Done();

  for(unsigned int i = 0; i < constraints.all.GetSize(); i++)
  {
    Condition* g = constraints.all[i]; // Now traverse in global mode

    if(g->IsLocalCondition())
    {
      LocalCondition* local = dynamic_cast<LocalCondition*>(g);
      double max  = local->CalcMaxValue();
      double mean = local->CalcMeanValue();
      if(out)
        *out << "\tmax_" << local->ToString() << max << "\tmean_" << local->ToString() << mean;
      iteration->Get("max_" + g->ToString())->SetValue(max);
      iteration->Get("mean_" + g->ToString())->SetValue(mean);
    }
    else
    {
      double value = g->GetValue();
      if(g->delta_logging) value = value - g->GetBoundValue();
      if(out) *out << "\t" << value;
      iteration->Get(g->ToString())->SetValue(value);
    }
  }

  if(out && log.design){
    StdVector<double> d;
    d.Resize(design->GetNumberOfVariables());
    design->WriteDesignToExtern(d.GetPointer(), false);
    for(unsigned int i = 0; i < design->GetNumberOfVariables(); i++){
      *out << "\t" << d[i];
    }
  }
  
  if(out && log.designGradient){
    StdVector<double> d;
    d.Resize(design->GetNumberOfVariables());
    d.window.Set(d);
    design->WriteGradientToExtern(d, DesignElement::COST_GRADIENT, DesignElement::PLAIN, NULL, false);
    for(unsigned int i = 0; i < design->GetNumberOfVariables(); i++){
      *out << "\t" << d[i];
    }
  }

  if(out) out->flush();
}


Optimization::MultipleExcitation::MultipleExcitation(bool multiple, PtrParamNode pn)
{
  // set defaults
  this->stride = 1;
  this->damping = 1.0;
  this->max_gain = 1e4;
  this->multiple_excitation_ = multiple;
  this->type_ = NO_TYPE; // to be eventually overwritten soon

  // if disabled, we don't read anything
  if(!multiple || pn == NULL) return;

  this->type_ = type.Parse(pn->Get("type")->As<std::string>());

  // adjust defaults
  if(pn->Has("adjustWeights")) {
    this->stride   = pn->Get("adjustWeights")->Get("stride")->As<Integer>();
    this->max_gain = pn->Get("adjustWeights")->Get("max_gain")->As<Double>();
    this->damping  = pn->Get("adjustWeights")->Get("damping")->As<Double>();
  }
}

void Optimization::MultipleExcitation::ToInfo(PtrParamNode in) const
{
  if(!multiple_excitation_) return;

  if(type_ == META_OBJECTIVE)  {
    PtrParamNode in_ = in->Get("metaObjective");
    in_->Get("type")->SetValue("best_value");
    in_->Get("damping")->SetValue(damping);
    in_->Get("stride")->SetValue(stride);
    in_->Get("max_gain")->SetValue(max_gain);
  } else
    in->Get("type")->SetValue(type.ToString(type_));
}


Excitation::Excitation()
{
  this->index = -1; // must be updated
  this->frequency = -1.0;
  this->f_link = NULL;
  this->weight = 1.0;
  this->normalized_weight = 1.0;
  linForms = new StdVector<LinearFormContext*>();  
  AddLinFormsFromAssemble();
}

void Excitation::AddLinFormsFromAssemble(){
  // all already set linear Forms
  StdVector<LinearFormContext*>& assLinForms = domain->GetBasePDE()->getPDE_assemble()->GetLinForms();
  for(unsigned int i = 0; i < assLinForms.GetSize(); ++i){
    linForms->Push_back(assLinForms[i]);
  }
}

void Excitation::Apply()
{
  domain->GetOptimization()->applied_excitation = this;
  Assemble* a = domain->GetBasePDE()->getPDE_assemble();
  a->SetLinForms(linForms); // this works always, if not used just a copy of assembles linearForms
  if(! loads.IsEmpty()){
    a->SetLoads(loads);
  }
  // a frequency cannot really be applied but has to be used as parameter
  // in the driver call
}

double Excitation::GetOmega()
{
  if(frequency < 0.0)
    EXCEPTION("No frequency given");
  return 2 * M_PI * frequency;
}

double Excitation::GetFactor(Function* cost)
{
  double factor = 1.0; // default
  
  if(cost->FactorOmegaOmega())
  {
    if(frequency < 0.0) EXCEPTION("No frequency given");
    factor *= 4 * M_PI * M_PI * frequency * frequency;
  }
  
  return factor;
}

double Excitation::GetWeightedFactor(Function* cost)
{
  return normalized_weight * GetFactor(cost);
}

void Excitation::ReadTrackings(PtrParamNode ts){
  StdPDE* mech = domain->GetStdPDE("mechanic");
  trackings.Clear();
  ParamNodeList tracking_list = ts->GetChildren();
  for(unsigned int i = 0; i < tracking_list.GetSize(); i++){
    std::string name, dof, value;
    dof = "";
    tracking_list[i]->GetValue("name", name);
    tracking_list[i]->GetValue("dof", dof, ParamNode::PASS);
    tracking_list[i]->GetValue("value", value);
    shared_ptr<LoadBc> actLoad( new LoadBc );
    shared_ptr<EntityList> actList = domain->GetGrid()->GetEntityList( EntityList::NODE_LIST, name, EntityList::NAMED_NODES );
    actLoad->entities = actList;
    actLoad->eqnMap = mech->GetEqnMap();
    if ( dof == "" ) {
      actLoad->dof = 1;
    } else {
      actLoad->dof = mech->GetResultInfo(MECH_DISPLACEMENT)->GetDofIndex(dof);
    }
    actLoad->value = value;
    trackings.Push_back(actLoad);
  }
}

void Excitation::ReadLoads(PtrParamNode ls)
{
  // loads
  loads.Clear();
  MechPDE* mech = (MechPDE*)domain->GetSinglePDE("mechanic");
  mech->ReadLoads(ls->GetList("load"), loads);

  // pressures
  StdVector<shared_ptr<EntityList> > pressSurf;
  StdVector<std::string> pressVals;
  StdVector<std::string> pressPhase;
  mech->ReadPressureLoadsFromXML(ls, pressSurf, pressVals, pressPhase);
  mech->DefinePressureIntegrators(pressSurf, pressVals, pressPhase, linForms);
  
  // regionLoads
  std::map<RegionIdType, SinglePDE::RegionLoad> regionLoads;
  mech->ReadRegionLoadsFromXML(ls, regionLoads);
  mech->DefineRegionLoadIntegrators(regionLoads, linForms);
}

void Excitation::ReadTestStrain(const Vector<double>& vec)
{
  this->test_strain = vec;

  loads.Clear();
  MechPDE* mech = dynamic_cast<MechPDE*>(domain->GetSinglePDE("mechanic"));
  mech->DefineTestStrainIntegrators(vec, linForms);
}

void Excitation::SetPolarizationMatrixRHS(const Vector<double>& mechp,
    const Vector<double>& elecp, const int num)
{
  loads.Clear();
  MechPDE* mech = dynamic_cast<MechPDE*>(domain->GetSinglePDE("mechanic"));
  mech->DefinePolarizationMatrixIntegrators(mechp, linForms, num);
  
  ElecPDE* elec = dynamic_cast<ElecPDE*>(domain->GetSinglePDE("electrostatic"));
  elec->DefinePolarizationMatrixIntegrators(elecp, linForms, num);
}

Optimization::Log::Log()
{
  this->design = false;
  this->designGradient = false;
  this->file = NULL;
  this->fileHeader = "";
}

void Optimization::Log::Init(const string& log_name, PtrParamNode pn_log)
{
  if(log_name != "false")
  {
    string name = log_name == "[problem]" || log_name == "true" ? (progOpts->GetSimName() + ".plot.dat") : log_name;
    file = new ofstream(name.c_str());
    if(file == NULL)
      throw Exception("cannot open log file " + name + " for writing");

    if(pn_log != NULL)
    {
      design = pn_log->Get("design")->As<bool>();
      designGradient = pn_log->Get("designGradient")->As<bool>();
    }
  }
}

 Optimization::Log::~Log()
 {
   // if write to file close it
   if(file != NULL)
   {
     file->close();
     delete file;
     file = NULL;
   }
 }
