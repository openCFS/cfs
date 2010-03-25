#include <string>
#include <cmath>
#include <limits>
#include <iomanip>

#include "Optimization/SIMP.hh"
#include "Optimization/DesignSpace.hh"
#include "Optimization/DesignElement.hh"
#include "Optimization/DesignStructure.hh"
#include "Optimization/OptimizationMaterial.hh"
#include "Domain/domain.hh"
#include "Domain/surfElem.hh"
#include "General/exception.hh"
#include "Utils/basenodestoresol.hh"
#include "PDE/SinglePDE.hh"
#include "PDE/mechPDE.hh"
#include "PDE/heatCondPDE.hh"
#include "PDE/acousticPDE.hh"
#include "Forms/baseForm.hh"
#include "Forms/linSurfForm.hh"
#include "Elements/basefe.hh"
#include "CoupledPDE/DirectCoupledPDE.hh"
#include "Utils/StdVector.hh"
#include "MatVec/vector.hh"
#include "Utils/result.hh"
#include "Utils/mathParser/mathParser.hh"
#include "Driver/assemble.hh"
#include "Driver/baseSolveStep.hh"
#include "Driver/harmonicDriver.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/programOptions.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "OLAS/algsys/basesystem.hh"
#include "OLAS/algsys/baseentrymanipulator.hh"
#include "MatVec/basematrix.hh"
#include "MatVec/stdmatrix.hh"
#include "Materials/mechanicMaterial.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"

using namespace CoupledField;
using namespace std;

DECLARE_LOG(conditions)
DEFINE_LOG(conditions, "conditions")

DECLARE_LOG(em)
DEFINE_LOG(em, "ersatzMaterial")

Enum<ErsatzMaterial::Method> ErsatzMaterial::method;

ErsatzMaterial::ErsatzMaterial() :
  Optimization(),
  material(NULL), // to be set in PostInit()
  dim(domain->GetGrid()->GetDim())
{
  /** We store here the solution */
  tracking_ = NULL;
  volume_fraction_ = 0.0;
  structure_ = NULL;
  pde = NULL;
  acou_near_field_elems_ = shared_ptr<BaseResult>(new Result<complex<double> >());

  pn = param->Get("optimization")->Get("ersatzMaterial");

  method_ = method.Parse(pn->Get("method")->As<std::string>());

  homogenization_ = false;
  for(unsigned int i = 0; i < objectives.data.GetSize(); i++)
    if(objectives.data[i]->IsHomogenization()) homogenization_ = true;

  for(unsigned int i = 0; i < constraints.all.GetSize(); i++)
    if(constraints.all[i]->IsHomogenization()) homogenization_ = true;

  optInfoNode->Get(ParamNode::HEADER)->Get("homogenization")->SetValue(homogenization_);

  homogenizedTensor.Resize(dim == 2 ? 3 : 6, dim == 2 ? 3 : 6);
  homogenizedTensor.Init();

  // region stuff
  ParamNodeList region_list = pn->GetList("region");
  if(region_list.IsEmpty() && (method_ != SHAPE_OPT && method_ != SHAPE_PARAM_MAT))
    EXCEPTION("no region given!");

  for(unsigned int i=0; i < region_list.GetSize(); i++){
    std::string reg = region_list[i]->Has("name") ? region_list[i]->Get("name")->As<std::string>() : region_list[i]->As<std::string>();
    if(!domain->GetGrid()->GetRegion().IsValid(reg))
      throw Exception("region given in ersatzMaterial is invalid");
    regionIds.Push_back(domain->GetGrid()->GetRegion().Parse(reg));
  }


  // set up the design space elements, note PiezoSIMP have only POLARIZATION
  // this includes the transfer functions!
  ParamNodeList design_list = pn->GetList("design");
  ParamNodeList transfer_list = pn->GetList("transferFunction");
  ParamNodeList result = pn->GetList("result");
  design = DesignSpace::CreateInstance(regionIds, design_list, transfer_list, result, method_);
  // make basic loggings
  design->ToInfo(optInfoNode->Get(ParamNode::HEADER)->Get("designSpace"));

  assume_constant_element_matrices_ = design->IsRegular() && method_ != ErsatzMaterial::PARAM_MAT;

  // optionally write the densities to an xml file
  if(pn->Has("export"))
  {
    std::string name = pn->Get("export/file")->As<std::string>();
    if(name == "[problem]") name = progOpts->GetSimName() + ".density.xml";
    exportDesign = CreateExportDesign(name, design_list, transfer_list);
    exportDesignAllIterations = pn->Get("export/save")->As<std::string>() == "all";
    exportDesignFinallyOnly   = pn->Get("export/write")->As<std::string>() == "finally";
  }

  // we assume always to have either a mechanic or heatcond pde.
  SetPDEs();

  harmonic = BasePDE::IsComplex(pde->GetAnalysisType());

  if(homogenization_ && !pde->HasPeriodicBC())
    EXCEPTION("homogenization requires periodic boundary conditions");

  // Get the assemble class
  assemble_ = pde->getPDE_assemble();

  // manipulate the PDE (the integrator of its bilinear form) such that
  // it is set up with the new density values in the optimization loops
  for(unsigned int r=0; r < regionIds.GetSize(); r++){
    if(pde->GetName() == "mechanic")
    {
      GetForm(regionIds[r], pde, pde, "linElastInt")->SetSolDependent(true);
      if(harmonic)
        GetForm(regionIds[r], pde, pde, "MassInt")->SetSolDependent(true);
    }
    if(pde->GetName() == "heatConduction")
      GetForm(regionIds[r], pde, pde, "LaplaceInt")->SetSolDependent(true);
  }

  // check our constraints, the shall have only valid designs
  for(unsigned int i = 0; i < constraints.all.GetSize();  i++)
  {
    Condition* g = constraints.all[i];
    if(g->design == DesignElement::UNITY && (method_ == SHAPE_OPT || method_ == SHAPE_PARAM_MAT))
      continue;

    if(g->design == DesignElement::TENSOR_TRACE && (method_ == PARAM_MAT || method_ == SHAPE_PARAM_MAT))
      continue;

    if(design->FindDesign(g->design, false) == -1)
      throw Exception("constraint " + g->ToString() + " operates on invalid design variable");

  }


  // give the domain this data, s.th. the ersatz material approach is applied
  domain->SetErsatzMaterial(design);

  // add optimization results to the pde
  design->AppendOptimizationResults(pde);

  // forward and adjoint are initialized in PostInit()

  // it's cheap as well -> no multiple load cases yet!
  tracking_ = new Solution(this);

  // check for multiple loadcases (might be frequencies)
  PrepareMultipleExcitations();
}

ErsatzMaterial::~ErsatzMaterial()
{
  delete tracking_;

  // if write to file close the xml envelope and the file
  if(exportDesign != NULL)
  {
    if(exportDesignFinallyOnly) { // in CommitIteration or here?
      std::string name = pn->Get("export/file")->As<std::string>();
      if(name == "[problem]") name = progOpts->GetSimName() + ".density.xml";
      exportDesign->ToFile( name); 
    }
    exportDesign.reset();
  }

  delete material;

  delete structure_;

  // "remove" the ersatzMaterial (=data) from the domain
  domain->SetErsatzMaterial(NULL);
}

void ErsatzMaterial::PostInit()
{
  // might be constructed in SIMP::PostInit()
  if(structure_ == NULL)
    structure_ = new DesignStructure(this);

  // the constraints size is known and the shapeDesign constructor is finished -> PostInit design
  design->PostInit(objectives.data.GetSize(), constraints.all.GetSize());
  // post init slope constraints when the design is there
  constraints.PostProc(design, structure_);

  // number of design variables in shape optimization is only known here
  // make a copy of the old iteration to calculate the move
  last_iteration.Resize(design->GetNumberOfVariables());

  // note the difference between function evaluations (line search) and iterations!
  last_evaluation.Resize(design->GetNumberOfVariables());

  // for simplicity we always store forward and adjoint and have a multiple object
  forward.Init(this);
  adjoint.Init(this); // it's cheap - even if not used


  for(unsigned int i = 0; i < objectives.data.GetSize(); i++)
  {
    Objective* cost = objectives.data[i];

    std::string objective = "'" + cost->type.ToString(cost->GetType()) + "'";

    // checks
    switch(cost->GetType())
    {
      case Objective::COMPLIANCE:
      case Objective::OUTPUT: // it would work but is saver not to allow
      case Objective::TRACKING:
      case Objective::HOMOGENIZATION_TENSOR:
      case Objective::HOMOGENIZATION_TRACKING:
      case Objective::POISSONS_RATIO:
      case Objective::YOUNGS_MODULUS:
      case Objective::ELEC_ENERGY: // it simply does not work yet in the harmonics
        if(harmonic) throw Exception(objective + " is only for static state problems");
        break;

      case Objective::DYNAMIC_OUTPUT:
      case Objective::GLOBAL_DYNAMIC_COMPLIANCE:
      case Objective::ABS_DYN_OUTPUT_SQUARED:
        if(!harmonic) throw Exception(objective + " is only for harmonic state problems");
        break;

      default:
        break;
    }

    // actions
    switch(cost->GetType())
    {
      case Objective::OUTPUT:
      case Objective::DYNAMIC_OUTPUT:
      case Objective::ABS_DYN_OUTPUT_SQUARED:
      {
        PtrParamNode output = cost->pn->Get("output");

        if(output->Has("displacement"))
          pde->ReadLoads(output->GetList("displacement"), output_nodes_);
        if(output->Has("elecPotential"))
          domain->GetSinglePDE("electrostatic")->ReadLoads(output->GetList("elecPotential"), output_nodes_);
        if(output->Has("acoustic"))
          domain->GetSinglePDE("acoustic")->ReadLoads(output->GetList("acoustic"), output_nodes_);

        if(output_nodes_.GetSize() == 0)
          throw Exception("no output optimization targets given");
        break;
      }
      case Objective::ACOU_NEAR_FIELD:
      {
        Grid* grid = domain->GetGrid();
        acou_near_field_elems_->SetResultInfo(ToPDE(ACOUSTIC)->GetResultInfo(ACOU_POTENTIAL));

        // id could be also named elements but we use a surface element region (ONE surfeRegion by schema)
        string region = cost->pn->Get("output/acousticNearField/surfRegion/name")->As<std::string>();
        shared_ptr<EntityList> list = grid->GetEntityList(EntityList::SURF_ELEM_LIST, region, EntityList::REGION);
        acou_near_field_elems_->SetEntityList(list);
      }
      default:
        break;
    }
  }

  // create Material class
  if(pde != NULL)
  {
    OptimizationMaterial::System system = 
        OptimizationMaterial::system.Parse(pn->Get("material")->As<std::string>());
    switch(system)
    {
    case OptimizationMaterial::PIEZO:
      material = new OptPiezoMat(this);
      break;

    case OptimizationMaterial::MECHANIC:
      material = new OptMechMat(this);
      break;

    case OptimizationMaterial::HEAT:
      material = new HeatMat(this);
      break;

    default: assert(false);
    }
    optInfoNode->Get(ParamNode::HEADER)->Get("material")->SetValue(OptimizationMaterial::system.ToString(material->GetSystem()));
  }

  Optimization::PostInit();
}


void ErsatzMaterial::PrepareMultipleExcitations()
{
  // by definition, we always have at least one excitation.
  // This does not necessarily need to ba a load (but might be voltage, pressure, ...)
  excitations.Resize(1);
  excitations[0].index = 0;

  PtrParamNode pncf = param->Get("optimization")->Get("costFunction");

  // the actual multipleExcitation description is read in Optimization as part of
  // objective function block

  int num_freq  = !harmonic ? 0 : dynamic_cast<HarmonicDriver*>(domain->GetDriver())->freqs.GetSize();
  int num_loads = assemble_->GetLoads().GetSize(); // to be faked later for homogenization test strains

  ParamNodeList pnexcitations;

  double weight_sum = 0.0;

  // plausibility check for homogenization
  if(homogenization_ && (!multiple_excitation->IsEnabled() || !multiple_excitation->DoHomogenization()))
    throw Exception("A homogenization objective/constraint is set but no homogenization test strain excitation");
  if(multiple_excitation->IsEnabled() && multiple_excitation->DoHomogenization() && !homogenization_)
    throw Exception("No homogenization objective/constraint for homogenization test strain excitation");

  // initialize data and do simple plausibilty check. Note that also 1 is "multiple"
  if(multiple_excitation->IsEnabled())
  {
    // either every single load is an excitation (Fabian's way) (done before)
    // or allow combinations of loads, pressures, regionLoads and trackings in one excitation (Bastian) (only non-harmonic) (this is done here)
    if(!harmonic && pncf->Has("multipleExcitation/excitations"))
    {
      pnexcitations = pncf->Get("multipleExcitation/excitations")->GetChildren();
      num_loads = pnexcitations.GetSize();
    }

    // decide if we have multiple loads or multiple frequencies
    // we cannot do both!
    if(num_freq > 1 && num_loads > 1)
      throw Exception("Cannot to concurrent multiple excitations for multiple loads and multiple frequencies");

    // sets and resizes excitations with strain loads
    if(multiple_excitation->DoHomogenization())
    {
      num_loads = SetHomogenizationTestStrains();
      weight_sum = 1; // all 0 but the last 1
    }

    // sets and resizes excitations with polarization matrix loads
    if(multiple_excitation->DoPolarizationMatrix())
    {
      // FIXME: maybe more sanity checks needed
      assert(!multiple_excitation->DoHomogenization()); // cannot do both
      
      num_loads = SetPolarizationMatrixExcitations();
      weight_sum = 1; // all 0 but the last 1
    }
    
    // the following is validated by above and 1 frequency and 0 loads is harmless
    if(num_freq > 1)  excitations.Resize(num_freq);
    if(num_loads > 1) excitations.Resize(num_loads); // redundant for homogenization

    // we average the solutions(s) only for output.
    // In the calculations we average the individual calculations
  }


  // read in tracking parameters from XML, for the first and only load
  if(pncf->Has("trackings")){
    excitations[0].ReadTrackings(pncf->Get("trackings"));
  }

  // this sets the first and only excitation even when we have multiple harmonic forward case
  // but not multiple excitations. Then only the first frquency is called.
  // Fills the excitations list with the data provided in the xml file as problem description
  if(harmonic)
  {
    HarmonicDriver* hd = dynamic_cast<HarmonicDriver*>(domain->GetDriver());

    for(unsigned int i = 0; i < excitations.GetSize(); i++)
    {
      excitations[i].index = i;
      excitations[i].frequency = hd->freqs[i].freq;
      excitations[i].weight    = hd->freqs[i].weight;
      excitations[i].f_link    = &hd->freqs[i];

      weight_sum += excitations[i].weight;
    }
  }
  if(!harmonic && multiple_excitation->IsEnabled()) // multiple loads case
  {
    MathParser * parser = domain->GetMathParser();
    unsigned int handle = parser->GetNewHandle();

    if(pnexcitations.GetSize() > 0)
    { //
      for(unsigned int i = 0; i < excitations.GetSize(); i++)
      {
        parser->SetExpr(handle, pnexcitations[i]->Get("weight")->As<std::string>());
        const double weight = parser->Eval(handle);
        excitations[i].weight = weight;
        weight_sum += weight;

        if(pnexcitations[i]->Has("trackings")){
          excitations[i].ReadTrackings(pnexcitations[i]->Get("trackings"));
        }
        excitations[i].ReadLoads(pnexcitations[i]->Get("loads"));
      }

    }
    if(pnexcitations.GetSize() == 0 
        && !multiple_excitation->DoHomogenization()
        && !multiple_excitation->DoPolarizationMatrix())
    {
      LoadList loads = assemble_->GetLoads();

      for(unsigned int i = 0; i < excitations.GetSize(); i++)
      {
        if(num_loads > (int) i)
        {
          excitations[i].loads.Push_back(loads[i]);

          parser->SetExpr(handle, loads[i]->weight);
          const double weight = parser->Eval(handle);
          excitations[i].weight = weight;

          weight_sum += weight;
        }
      }
    }

    parser->ReleaseHandle(handle);
  }

  // output summary
  // -------------------

  // calc the inital normalized_weight and print info.

  // for many concurrent mechanical loads do no print information
  if(multiple_excitation->IsEnabled() || harmonic)
  {
    PtrParamNode in = optInfoNode->Get(ParamNode::HEADER)->Get("excitations");

    if(!multiple_excitation->IsEnabled() && num_freq > 1 && optimizer_ != EVALUATE_INITIAL_DESIGN)
    {
      stringstream ss;
      ss << "Solve only for 1. frequency (" << excitations[0].frequency << "Hz) of "
         << num_freq << " as multiple excitations are disabled";
      in->Get(ParamNode::WARNING)->SetValue(ss.str());
    }

    // communicate what we have but also normalize the weights!
    for(unsigned int i = 0; i < excitations.GetSize(); i++)
    {
      Excitation& ex = excitations[i];
      ex.index = i;

      ex.normalized_weight = ex.weight / weight_sum;

      PtrParamNode exin = in->Get("excitation", ParamNode::APPEND);
      exin->Get("index")->SetValue(ex.index);
      if(! ex.loads.IsEmpty())
        ex.loads[0]->ToInfo(exin->Get("load"));
      if(ex.frequency >= 0.0)
        exin->Get("frequency")->SetValue(ex.frequency);
      if(ex.test_strain.GetSize() > 0)
        exin->Get("testStrain")->SetValue(ex.test_strain.ToString());
    }
  }
}


int ErsatzMaterial::SetHomogenizationTestStrains()
{
  assert(homogenization_);

  double ts[6][6] =  { {1.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                        {0.0, 1.0, 0.0, 0.0, 0.0, 0.0 },
                        {0.0, 0.0, 1.0, 0.0, 0.0, 0.0 },
                        {0.0, 0.0, 0.0, 1.0, 0.0, 0.0 },
                        {0.0, 0.0, 0.0, 0.0, 1.0, 0.0 },
                        {0.0, 0.0, 0.0, 0.0, 0.0, 1.0 } };

  Vector<double> vec;

  int cases = dim == 2 ? 3 : 6;
  excitations.Resize(cases);

  for(int i = 0, cnt = 0; i < 6; ++i)
  {
    // in 2D only 0, 1 and 5
    if(dim == 2 && (i == 2 || i == 3 || i == 4)) continue;

    vec.Fill(ts[i], 6);

    excitations[cnt].ReadTestStrain(vec);
    // The homogenized tensor can only be evaluated for the last excitation!
    excitations[cnt].weight = i < cases-1 ? 0.0 : 1.0;
    ++cnt;
  }

  return excitations.GetSize();
}

int ErsatzMaterial::SetPolarizationMatrixExcitations()
{
  // collect all necessary data
  int dim = domain->GetGrid()->GetDim();
  int dim1 = dim == 3 ? 6 : 3;
  // int dim2 = dim == 3 ? 3 : 2;
  
  // read material information from xml-file
  // check for tensor element
  PtrParamNode pol = param->Get("polarizationMatrix");
  Matrix<double> mechmatrix;
  Matrix<double> elecmatrix;
  Matrix<double> couplmatrix;
  
  mechmatrix.Resize(dim1, dim1);
  elecmatrix.Resize(dim, dim);
  couplmatrix.Resize(dim1, dim);
  
  ParamTools::AsTensor<double>(pol->Get("mechTensor/real"), dim1, dim1, mechmatrix);
  ParamTools::AsTensor<double>(pol->Get("elecTensor/real"), dim, dim, elecmatrix);
  ParamTools::AsTensor<double>(pol->Get("couplingTensor/real"), dim1, dim, couplmatrix);

  const int cases(dim == 2 ? 5 : 9);
  excitations.Resize(cases);

  // decompose the vector with the material parameters into the
  // part for mech and the one for elec
  Vector<double> mechpart(dim1, 0.0);
  Vector<double> elecpart(dim, 0.0);
  
  for(int i = 0; i < cases; ++i)
  {
    // this must be the entry of A^0
    if(i < cases - dim)
    {
      for(int n = 0; n < dim1; ++n)
        mechpart[n] = mechmatrix[n][i];
      for(int n = 0; n < dim; ++n)
        elecpart[n] = couplmatrix[i][n];
    }
    else
    {
      for(int n = 0; n < dim1; ++n)
        mechpart[n] = couplmatrix[n][i];
      for(int n = 0; n < dim; ++n)
        elecpart[n] = elecmatrix[i][n];
    }

    excitations[i].SetPolarizationMatrixRHS(mechpart, elecpart, i);
    // The homogenized tensor can only be evaluated for the last excitation!
    excitations[i].weight = i < cases-1 ? 0.0 : 1.0;
  }

  return cases;
}


void ErsatzMaterial::NormalizeMultipleExcitations()
{
  if(!multiple_excitation->IsEnabled()) return;

  const unsigned int exsi(excitations.GetSize());
  if(multiple_excitation->DoAdjustWeights())
  {
    // find best cost
    double best;
    if(objectives.DoMaximize())
    {
      best = numeric_limits<double>::min();
      for(unsigned int i = 0; i < exsi; i++) best = max(best, excitations[i].cost);
    }
    else
    {
      best = numeric_limits<double>::max();
      for(unsigned int i = 0; i < exsi; i++) best = min(best, excitations[i].cost);
    }

    for(unsigned int i = 0; i < exsi; i++)
    {
      // this is for "best_value": We push (weight the gradient) the worst result most to equalize/over-/undercompensate
      // the algorithm is the foolowing: w_k^p J_k = const -> we define w_k' for J_k = J_best to be one.
      // By this we can compute all w_k'. The we sum up all w_k' to sum_w_k' and get w_k = w_k'/sum_w_k'
      double t = objectives.DoMaximize() ? best / excitations[i].cost : excitations[i].cost / best;
      t = std::pow(t, 1.0/multiple_excitation->damping); // fails for negative quotient with damping < 0
      excitations[i].weight = min(t, multiple_excitation->max_gain);;
      LOG_DBG(em) << "NormalizeMultipleExcitations: excitation[" << i << "] best: "
      << best << " gain: " << t << " weight: " << excitations[i].weight;
    }
  }

  // normalize
  double weight_sum = 0.0;

  for(unsigned int i = 0; i < exsi; i++)
    weight_sum += excitations[i].weight;

  for(unsigned int i = 0; i < exsi; i++)
  {
    excitations[i].normalized_weight = excitations[i].weight / weight_sum;
    LOG_DBG(em) << "NormalizeMultipleExcitations: excitation[" << i << "].normalized_weight <- " << excitations[i].normalized_weight;
  }
}

PtrParamNode ErsatzMaterial::CreateExportDesign(const std::string& filename, ParamNodeList& des, ParamNodeList& tfs)
{
   PtrParamNode in = PtrParamNode(new ParamNode(ParamNode::INSERT));
   in->SetName("cfsErsatzMaterial");

   // write header
   PtrParamNode in_ = in->Get("header");
   for(unsigned int i = 0; i < des.GetSize(); i++)
     in_->Get("dummy", ParamNode::APPEND)->SetValue(des[i]); // name is overwritten

   for(unsigned int i = 0; i < tfs.GetSize(); i++)
     in_->Get("dummy", ParamNode::APPEND)->SetValue(tfs[i]);

   return in;

}

void ErsatzMaterial::StoreResults(double step_val)
{
  // check if we have to play with the results, forward is default mode!
  if(commitMode_ == FORWARD || commitMode_ == BOTH)
  {
    // in case of FORWARD this is redundant but it hanldes multiple exitations ans is cheap

    // check for multiple to find out what is actually to be restores
    for(map<Application, SinglePDE*>::iterator it = pdes.begin(); it != pdes.end(); ++it)
    {
      if(multiple_excitation->IsEnabled())
        Solution::Write(it->second, forward.GetMultiple()); // todo: later a pde dependcy is required
      else
        forward.Get(0)->Write(it->second); // todo: currently always the same
    }
    // call real implementation in Optimization
    Optimization::StoreResults(step_val);
  }
  if(commitMode_ == ADJOINT || commitMode_ == BOTH)
  {
    // "forward" is not possible -> hence it is both_cases or adjoint
    if(adjoint.Get(0)->GetVector(Solution::RAW_VECTOR) == NULL) EXCEPTION("Adjoint solution cannot be written on 'commit' because there is none");

    for(map<Application, SinglePDE*>::iterator it = pdes.begin(); it != pdes.end(); ++it)
    {
      if(multiple_excitation->IsEnabled())
        Solution::Write(it->second, adjoint.GetMultiple());
      else
        adjoint.Get(0)->Write(it->second);
    }

    step_val = step_val == -1 ? currentIteration + 0.5 : step_val + 0.5;
    Optimization::StoreResults(step_val);
  }
}


PtrParamNode ErsatzMaterial::CommitIteration(bool keep_iteration_number)
{
  // will write the cfs results and the log file
  PtrParamNode iter = Optimization::CommitIteration(keep_iteration_number);
  // add our multiple excitation stuff here (only in info.xml, this would be to complex for dat
  if(multiple_excitation->IsEnabled())
  {
    for(unsigned int i = 0; i < excitations.GetSize(); i++)
    {
      PtrParamNode e = iter->Get("excitation", ParamNode::APPEND);
      e->Get("index")->SetValue(excitations[i].index);
      e->Get("objective")->SetValue(excitations[i].cost);
      e->Get("normalized_weight")->SetValue(excitations[i].normalized_weight);
    }
  }

  if(homogenization_)
  {
    PtrParamNode in = iter->Get("homogenizedTensor");
    in->SetValue(new Matrix<double>(homogenizedTensor));
    in->Get("norm_L2")->SetValue(homogenizedTensor.NormL2());
    in = in->Get("isotropy");
    in->Get("err")->SetValue(MechanicMaterial::CalcIsotropyError(homogenizedTensor, false));
    in->Get("err_normed")->SetValue(MechanicMaterial::CalcIsotropyError(homogenizedTensor, true));
    in->Get("poissons_ratio")->SetValue(MechanicMaterial::CalcIsotropicPoissonsRatio(homogenizedTensor));
    in->Get("E")->SetValue(MechanicMaterial::CalcIsotropicYoungsModulus(homogenizedTensor));
  }

  if(exportDesign != NULL)
  {
    SetCurrentExportDesign(); // appends or replaces
    if(!exportDesignFinallyOnly) { // here or on destructor
      std::string name = pn->Get("export/file")->As<std::string>();
      if(name == "[problem]") name = progOpts->GetSimName() + ".density.xml";
      exportDesign->ToFile(name);
    }
  }

  return iter;
}

void ErsatzMaterial::SetCurrentExportDesign()
{
  PtrParamNode in = exportDesign->Get("set", exportDesignAllIterations ? ParamNode::APPEND : ParamNode::INSERT);
  // add the entry, note that the iteration counter was incremented in base implementation
  in->Get("id")->SetValue(currentIteration-1);

  // the loop can be quite long. When it is fille already, ParamNode::Get() is cached
  // but when it is filled it would search. Do it faster!

  // even for !exportDesignAllIterations in is empty in the first commit

  ParamNodeList& list = in->GetChildren();

  const unsigned int bias = 1; // there one "id" element as first element

  assert(list.GetSize() == bias || list.GetSize() == design->data.GetSize() + bias);
  assert(list[0]->GetName() == "id");

  bool append = list.GetSize() == bias;
  assert(!(exportDesignAllIterations && ! append));

  if(append)
    list.Resize(design->data.GetSize() + 1);

  for(unsigned int i = 0, s = design->data.GetSize(); i < s; ++i)
  {
    DesignElement* de = &design->data[i];
    PtrParamNode el = append ? in->SetNewChild("element", i + bias) : list[i + bias];
    el->Get("nr")->SetValue(de->elem->elemNum);
    el->Get("type")->SetValue(DesignElement::type.ToString(de->GetType()));
    el->Get("design")->SetValue(de->GetDesign(DesignElement::PLAIN));
    el->Get("gradient")->SetValue(de->GetObjectiveGradient(DesignElement::PLAIN));
  }
}

string ErsatzMaterial::GetIterationFrequency()
{
  if(!harmonic) return "";

  // make clear, when doing *real* multiple excitations, that this is no single
  // frequency result
  if(multiple_excitation->IsEnabled() && excitations.GetSize() > 1) return "(mult)";

  double frequency = dynamic_cast<HarmonicDriver*>(domain->GetDriver())->GetActFreq();
  // as we control the fractional digits, we do not use lexical_cast<string>
  stringstream ss;
  ss << fixed << std::setprecision(1) << frequency;
  return ss.str();
}


BaseForm* ErsatzMaterial::GetForm(const RegionIdType regionId, StdPDE* pde1, StdPDE* pde2, const std::string& integrator)
{
  return assemble_->GetBiLinForm(regionId, pde1, pde2, integrator)->GetIntegrator();
}

double ErsatzMaterial::CalcObjective()
{
  design->WriteDesignToExtern(last_evaluation.GetPointer());

  // in objective.value_ we store the sum over all excitations w/o penalty but with normalization
  // in excitation.cost we store the sum over all objectives with penalty but w/o normalization

  // reset the objective values such that we can sum up normalized but unpenalized values
  for(unsigned int o = 0; o < objectives.data.GetSize(); o++)
    objectives.data[o]->ResetValue();

  double result = 0.0;

  // the multiple excitation case is a special case - for all other cases this is executed once
  for(unsigned int e = 0; e < excitations.GetSize(); e++)
  {
    Excitation& excite = excitations[e];
    excite.cost = 0.0;

    for(unsigned int o = 0; o < objectives.data.GetSize(); o++)
    {
      Objective* f = objectives.data[o];

      // some objectives are only to be evaluated for the last excitation
      bool last = excite.index == (int) excitations.GetSize() - 1;
      if(f->DoEvaluateOnce() && !last) continue;

      double ov = CalcObjective(excite, f); // this is virtual!
      excite.cost += ov * f->GetPenalty();

      // we ignore the weight if the evaluation happens only once!
      double weight = f->DoEvaluateOnce() ? 1.0 : excite.normalized_weight;

      f->AddValue(ov * weight);

      result += ov * f->GetPenalty() * weight;
      LOG_DBG(em) << "CalcObjective: ex=" << e << " obj=" << f->type.ToString(f->GetType()) << " ov=" << ov
                  << " penalty" << f->GetPenalty() << " ex.cost=" << excite.cost << " nw=" << excite.normalized_weight
                  << " wei=" << weight << " f->val=" << f->GetValue() << " result=" << result;
    }
  }

  return result;
}

double ErsatzMaterial::CalcObjective(Excitation& excite, Objective* cost)
{
  if(harmonic) return CalcObjective<std::complex<double> >(excite, cost);
          else return CalcObjective<double>(excite, cost);
}


template <class T>
double ErsatzMaterial::CalcObjective(Excitation& excite, Objective* cost)
{
  assert(!(excite.index < (int) excitations.GetSize() - 1 && cost->DoEvaluateOnce()));

  switch(cost->GetType())
  {
  case Objective::COMPLIANCE:
    return CalcCompliance(excite, cost, NULL, false);

  case Objective::TRACKING:
    return CalcTracking(excite, cost, NULL, false);

  case Objective::TYCHONOFF:
    return IntegrateDesignVariable(cost, NULL, false, DesignElement::NO_TYPE, true, true, 2.0);

  case Objective::VOLUME:
  case Objective::PENALIZED_VOLUME: // the exponent is as parameter in const
  case Objective::GAP: // volume - penalized volume
    return CalcVolume(cost, NULL, false, true);

  case Objective::GLOBAL_DYNAMIC_COMPLIANCE:
    return CalcGlobalDynamicCompliance(excite, cost);

  case Objective::OUTPUT:
  case Objective::DYNAMIC_OUTPUT:
  case Objective::CONJUGATE_COMPLIANCE:
  case Objective::ABS_DYN_OUTPUT_SQUARED:
    return CalcOutputObjective<T>(excite, cost);

  case Objective::HOMOGENIZATION_TENSOR:
  {
    Matrix<double> hom_tensor = CalcHomogenizedTensor();
    if(cost->HasHomogenizationEntry())
    {
      return hom_tensor[cost->coord.first -1][cost->coord.second - 1];
    }
    else
    {
      std::cout << "Homogenized Tensor: " << std::endl << hom_tensor.ToString(0, true);
      return hom_tensor.NormL2();
    }
  }

  case Objective::HOMOGENIZATION_TRACKING:
  {
    Matrix<double> hom_tensor = CalcHomogenizedTensor();
    return 0.5 * cost->GetTensor().DiffNormL2(hom_tensor) * cost->GetTensor().DiffNormL2(hom_tensor);
  }

  case Objective::POISSONS_RATIO:
  case Objective::YOUNGS_MODULUS:
    return CalcPoissonsRatioAndYoungsModulus(cost, NULL, false);

  case Objective::ACOU_NEAR_FIELD:
    return CalcAcouNearField(excite, cost);

  case Objective::TEMPERATURE:
    return 1.0; // FIXMEHEAT

  default: throw Exception("objective no handled");
  }
}

void ErsatzMaterial::CalcObjectiveGradient(StdVector<double>* grad_out)
{
  // reset the cost gradients in the design elements and sum them up in a weighted way
  // to perform multiple loads
  design->Reset(DesignElement::COST_GRADIENT);

  for(unsigned int obj = 0; obj < objectives.data.GetSize(); obj++)
  {
    Objective* cost = objectives.data[obj];
    // the multiple excitation case is a special case - for all other cases this is executed once
    for(unsigned int idx = 0; idx < excitations.GetSize(); idx++)
    {
      Excitation& excite = excitations[idx];
      // some objectives are only to be evaluated for the last excitation
      bool last = excite.index == (int) excitations.GetSize() - 1;
      if(!last && cost->DoEvaluateOnce()) continue;

      switch(cost->GetType())
      {
        // Note, that in SIMP case we handle compliance already  in SIMO
        case Objective::COMPLIANCE:
          // the multiple load cases implementation for the gradient is in SIMP
          CalcCompliance(excite, cost, NULL, true);
          break;

        case Objective::TRACKING:
          CalcTracking(excite, cost, NULL, true);
          break;

        case Objective::VOLUME:
        case Objective::PENALIZED_VOLUME: // penalization is in the parameter
        case Objective::GAP: // volume - penalized volume
          CalcVolume(cost, NULL, true, true);
          break;

        case Objective::TYCHONOFF:
          IntegrateDesignVariable(cost, NULL, true, DesignElement::NO_TYPE, true, true, 2.0);
          break;

        case Objective::HOMOGENIZATION_TRACKING:
          CalcHomogenizedTrackingGradient(cost->GetTensor(), CalcHomogenizedTensor(), cost, NULL);
          break;

        case Objective::HOMOGENIZATION_TENSOR:
          // if there s no "coord" set it is only meant for evaluateInitialDesign for forward homogenization
          if(cost->HasHomogenizationEntry())
          {

            StdVector<double> tmp;
            CalcHomogenizedTensorEntry(cost->coord, true, tmp);
            for(unsigned int e = 0, ne = design->GetNumberOfElements(); e < ne; e++)
              design->data[e].AddGradient(cost, NULL, tmp[e]);
          }
          break;
        
        case Objective::POISSONS_RATIO:
        case Objective::YOUNGS_MODULUS:
          CalcPoissonsRatioAndYoungsModulus(cost, NULL, true);
          break;

        default:
          CalcObjectiveGradient(excite, cost);
      }
    }
  }

  if(grad_out != NULL)
    design->WriteGradientToExtern(*grad_out, DesignElement::COST_GRADIENT, DesignElement::SMART);

}

int ErsatzMaterial::GetSpecialResultIndex(Application app1, Application app2, CalcMode calcMode, Condition* constraint)
{
  stringstream label;
  label << application.ToString(app1) << "_" << application.ToString(app2);
  if(calcMode == CONJ_QUAD) label << "_quad";

  DesignElement& de = design->data[0];

  int index = -1;
  DesignElement::Detail detail = de.NONE;
  if(de.detail.IsValid(label.str()))
  {
    detail = de.detail.Parse(label.str());
    DesignElement::ValueSpecifier vs = constraint != NULL ? de.CONSTRAINT_GRADIENT : de.COST_GRADIENT;
    index = design->GetSpecialResultIndex(de.DEFAULT, vs, detail, de.PLAIN);
  }
  LOG_DBG2(em) << "GetSpecialResultIndex: label='" << label.str() << "' detail=" << de.detail.ToString(detail) << " index=" << index;
  return index;
}

template <class T>
double ErsatzMaterial::CalcU1KU2(TransferFunction* tf, StdVector<SingleVector*>& u1,
                       Application app, StdVector<SingleVector*>& u2,
                       SurfaceRef* rhs, double factor, CalcMode calcMode, Objective* f, Condition* g, int res_idx)
{
  LOG_DBG2(em) << "CalcU1KU2(): tf=" << (tf ? tf->ToString() : "NULL") << " #u1=" << u1.GetSize()
                 << " app=" << application.ToString(app) << " #u2=" << u2.GetSize() << " calcMode="
                 << calcMode << " factor=" << factor << " rhs=" << (rhs == NULL ? "NULL" : rhs->ToString(1));

  // This solves <l,K'*u-f'> or <u1, K' * u2 - f'> for all elements and adds it up to the element gradients

  assert(u1.GetSize() != 0);
  assert(u1.GetSize() == u2.GetSize());

  double sum = 0.0;
  // the dimensions of our matrix is determined by u1_vec and u2_vec.
  // mat will be filled by SetElementK where also the derivative form most cases is built in
  Matrix<T> mat(u1[0]->GetSize(), u2[0]->GetSize());
  Vector<T> mat_vec(u1[0]->GetSize());

  TransferFunction* rtf = rhs != NULL && rhs->valid ? design->GetTransferFunction(tf->GetDesign(), rhs->app) : NULL;

  // traverse over our elements
  // in ErsatzMaterialTensor case we loop over all elements, else only over the elements belonging to this design
  int elements = design->GetNumberOfElements();
  int base_lower = 0;
  int base_upper = 0;
  if(domain->HasErsatzMaterialTensor()){
    base_upper = design->data.GetSize();
  }else{
    base_lower = design->FindDesign(tf->GetDesign()) * elements;
    base_upper = base_lower + elements;
  }
  LOG_DBG2(em) << "elements=" << elements << " base=" << base_lower << " base_upper=" << base_upper;

  // create an element list to gain the iterator in the loop
  ElemList elemList(domain->GetGrid());

  // for ParamMat we need the derivative w.r.t. every designvariable, else the base loop is only run once
  for(int base = base_lower; base < base_upper; base += elements)
  {
    for(int e = 0 ; e < elements; e++)
    {
      Vector<T>& u1_vec = dynamic_cast<Vector<T>& >(*u1[e]);
      Vector<T>& u2_vec = dynamic_cast<Vector<T>& >(*u2[e]);

      DesignElement* de = &design->data[e + base];

      LOG_DBG3(em) << "u1:" << e << ": " << u1_vec.ToString();
      LOG_DBG3(em) << "u2:" << e << ": " << u2_vec.ToString();

      // u1^T K' u2 - f' -> find "K'"
      SetElementK(de, app, dynamic_cast<DenseMatrix*>(&mat), calcMode);
      LOG_DBG3(em) << "mat: " << mat.ToString();

      // We generelly solve u1^T K' u2 as (K' u1)^T u2 
      // u1^T K' u2 - f' -> calc K' u2"
      mat_vec = mat * u2_vec;
      LOG_DBG3(em) << "mat * u2: " << mat_vec.ToString();

      // u1^T K' u2 - f' -> calc "- f'"
      assert(!(calcMode == CONJ_QUAD && rtf != NULL)); // no sensitive rhs here!
      if(rtf != NULL) SubstractGradSurfaceRHS(de, rtf, rhs, mat_vec);
      LOG_DBG3(em) << "-f': " << mat_vec.ToString();

      // u1^T(K' u2 - f') -> calc "u1^T *" or <u1, *> 
      // the difference is the conjugate complex in the harmonic inner product case!
      T sp;
      if(calcMode == CONJ_QUAD) mat_vec.Inner(u1_vec, sp); // u1 = u2 = u! 
                           else sp = mat_vec * u1_vec ;

      // when doing complex Jensen 22.07.07 shows that we always have 2 * Re(lamda * grad S * u)
      // the factor gives the negative sign
      // in real case it is simple value = factor * sp.
      // factor shall be +/- 1!
      double this_value = factor;
      if(calcMode == STANDARD && harmonic) this_value *= 2 * ((complex<double>) sp).real(); // 2 * Re{...}
                                      else this_value *= ((complex<double>) sp).real(); // CONJ_QUAD or real STANDARD

      de->AddGradient(f, g, this_value);

      LOG_DBG3(em) << "<l,K'*u-f'> = " << sp << " -> " << this_value << " sum = " << de->GetGradient(f,g);

      sum += this_value;

      if(res_idx != -1) de->specialResult[res_idx] = this_value;
    }
  }

  return sum;
}

// template instantiation stuff
template double ErsatzMaterial::CalcU1KU2<double>(TransferFunction* tf, StdVector<SingleVector*>& u1,
                       Application app, StdVector<SingleVector*>& u2,
                       SurfaceRef* rhs, double factor, CalcMode calcMode, Objective* f, Condition* g, int res_idx);

template double ErsatzMaterial::CalcU1KU2<std::complex<double> >(TransferFunction* tf, StdVector<SingleVector*>& u1,
                       Application app, StdVector<SingleVector*>& u2,
                       SurfaceRef* rhs, double factor, CalcMode calcMode, Objective* f, Condition* g,  int res_idx);

template <class T>
void ErsatzMaterial::SubstractGradSurfaceRHS(DesignElement* de, TransferFunction* tf, SurfaceRef* ref, Vector<T>& in_out)
{
  // we have to find the nodes which are common between de->elem
  // and the surface element which is one dimension smaller

  // not all elements do necessary lay on a surface and then not all nodes
  assert(ref != NULL && ref->valid);


  // nodes (numbers) of our design element
  StdVector<unsigned int>& de_nodes = de->elem->connect;
  Vector<T>& rhs = dynamic_cast<Vector<T>& >(*(ref->vec));

  // in_out is scalar (potential) or vectorial (x,y,(z))
  assert(in_out.GetSize() >= de_nodes.GetSize());

  int dof = in_out.GetSize() / de_nodes.GetSize();
  assert(dof == 1 || dof == 2 || dof == 3);
  assert(dof == (int) rhs.GetSize());

  // all node numbers of the surface are in a set
  std::set<unsigned int>::iterator it;
  // compare with the node numbers of our design element
  for(unsigned int n = 0; n < de_nodes.GetSize(); n++)
  {
    it = ref->nodes.find(de_nodes[n]);
    if(it != ref->nodes.end())
    {
      LOG_DBG3(em) << "SubstractGradSurfaceRHS : node " << n << " is common with elem "
                     << *it << " in surface: K'u = " << in_out.ToString();

      // find the the sensitivity of the rhs w.r.t the design volume element!
      double factor = tf->Derivative(de);

      // we do not really construct a rhs vector (with some/many zeros) but substract
      // for all design nodes common with the surface directly the entries
      for(int d = 0; d < dof; d++)
        in_out[n*dof + d] -= factor * rhs[d];
      LOG_DBG3(em) << "... -" << factor << "*" << rhs.ToString()
                     << " -> " << in_out.ToString();
    }
  }
}

void ErsatzMaterial::CalcConstraintGradient(Condition* g, StdVector<double>* grad_out)
{
  design->Reset(DesignElement::CONSTRAINT_GRADIENT, g != NULL ? g->design : DesignElement::DEFAULT);
  CalcConstraint(g, true, grad_out);
}

double ErsatzMaterial::CalcConstraint(Condition* constraint)
{
  return CalcConstraint(constraint, false); // no gradient
}

double ErsatzMaterial::CalcConstraint(Condition* g, bool derivative, StdVector<double>* grad_out)
{
  // Just for access of enums
  DesignElement de;

  // handle default parameter
  if(g == NULL)
  {
    if(constraints.active.GetSize() != 1)
      throw Exception("default constraint only valid with exactly one constraint");
    g = constraints.active[0];
  }
  double result = 0.0;

  switch(g->GetType())
  {
    case Condition::VOLUME:
    case Condition::PENALIZED_VOLUME:
    case Condition::GAP:
         result = CalcVolume(NULL, g, derivative, true);
         break;

    case Condition::REALVOLUME:
         result = CalcVolume(NULL, g, derivative, false);
         break;
         
    case Condition::COMPLIANCE:
         assert(excitations.GetSize() == 1);
         result = CalcCompliance(excitations[0], NULL, g, derivative);
         break;

    case Condition::GREYNESS:
         result = CalcGreyness(g, derivative);
         break;

    case Condition::SLOPE:
         result = CalcSlopeConstraint(g, derivative);
         break;

    case Condition::CHECKERBOARD:
         result = CalcCheckerboard(g, derivative);
         break;

    case Condition::HOMOGENIZATION_TENSOR:
         result = CalcHomogenizedTensorConstraint(g, derivative);
         break;

    case Condition::HOMOGENIZATION_TRACKING:
    {
      double diff = g->GetTensor().DiffNormL2(CalcHomogenizedTensor());
      result = 0.5 * diff * diff;
      break;
    }

    case Objective::POISSONS_RATIO:
    case Objective::YOUNGS_MODULUS:
          result = CalcPoissonsRatioAndYoungsModulus(NULL, g, derivative);
          break;

    default:
      throw Exception("Constraint not implemented", __FILE__, __LINE__);
  }

  // there is no single scalar value for the gradient
  if(!derivative)
    g->SetValue(result);

  // copies from the design element gradient data to a memory array for external optimizers
  if(grad_out != NULL)
    design->WriteGradientToExtern(*grad_out, de.CONSTRAINT_GRADIENT, de.SMART, g);

  // if there is a <result ... value="constraintGradient" detail="penalizedVolume/*"
  if(derivative && g->special_result_idx != -1)
  {
    int base = design->FindDesign(g->design);
    int n    = design->GetNumberOfElements();
    for(int i = n * base; i < n * (base + 1); i++)
      design->data[i].specialResult[g->special_result_idx] = design->data[i].GetGradient(NULL, g);
  }

  LOG_DBG2(em) << "CalcConstraint " << g->ToString() << " -> " << (derivative ? "derivative" : lexical_cast<std::string>(result));
  return result;
}

double ErsatzMaterial::IntegrateDesignVariable(Objective* f, Condition* g, bool derivative,
    DesignElement::Type dtype, bool normalized, bool scale, double exponent)
{
  // this replaces and enhances calculation of volume, it is used by regularization
  // when not assuming a regular grid, computation of Volume is not as simple
  // we also consider working only on a given region, when used as constrain
  // use dtype == NO_TYPE to iterate over all designs, but do not calculate tensor trace even if available

  // do we want the physical value? Don't make GTF() fault tolerant as we assume the physcial value!
  TransferFunction* tf = Function::GetFunction(f, g)->IsPhysical() ? design->GetTransferFunction(dtype, MECH) : NULL;

  Matrix<Double> cornerCoords;
  Grid* grid = domain->GetGrid();
  double fraction = f != NULL ? volume_fraction_ : g->volume_fraction;
  bool allDesignsRelevant = dtype == DesignElement::TENSOR_TRACE || dtype == DesignElement::DEFAULT || dtype == DesignElement::NO_TYPE;
  bool calculateTensorTrace = domain->HasErsatzMaterialTensor() && (dtype == DesignElement::TENSOR_TRACE || dtype == DesignElement::DEFAULT); // tensor trace is calculated if dtype == DEFAULT or TENSOR_TRACE and a tensor available
  if(calculateTensorTrace && scale){
    throw Exception("Cannot calculate Tensor Trace Volume on scaled design variables!");
  }
  // check whether we have to calculate the full volume
  if(normalized){
    if(g != NULL && g->design == DesignElement::UNITY){ // this will always return 1, does only make sense for unnormalized (real) volume
      // the gradient is not set, it is really 0 on every element but should not be used
      //TODO: Do We need a warning here? 
      return(derivative ? 0.0 : 1.0);
    }
  }else{
    if(design->IsRegular()){ // we use 1/(the volume of the first element) as fraction
      grid->GetElemNodesCoord(cornerCoords, design->data[0].elem->connect, true);
      fraction = 1.0 / design->data[0].elem->ptElem->CalcVolume(cornerCoords, false);
    }else{
      fraction = 1.0; // if no normalization needed, fraction is simply 1.0
    }
  }
  if( fraction == 0.0 ){
    for(unsigned int d = 0; d < design->design.GetSize(); d++){
      if(allDesignsRelevant || design->design[d] == dtype){
        const unsigned int base = d * design->elements_;
        for(unsigned int r = 0; r < design->regions_.GetSize(); r++){
          if(f != NULL || g->IsForRegion(design->regions_[r].regionId)){
            if(design->IsRegular()){
              fraction += design->regions_[r].elements;
            }else{
              const unsigned int u = base + design->regions_[r].base + design->regions_[r].elements;
              for(unsigned int i = base + design->regions_[r].base; i < u; i++){
                DesignElement* de = &design->data[i];
                grid->GetElemNodesCoord(cornerCoords, de->elem->connect, true);
                fraction += de->elem->ptElem->CalcVolume(cornerCoords, false);
              }
            }
          }
        }
      }
    }
    fraction = 1.0 / fraction;
    if(g == NULL){
      volume_fraction_ = fraction;
    }else{
      g->volume_fraction = fraction;
    }
  }

  double sum = 0.0;

  for(unsigned int d = 0; d < design->design.GetSize(); d++){
    if(allDesignsRelevant || design->design[d] == dtype){
      const unsigned int base = d * design->elements_;
      for(unsigned int r = 0; r < design->regions_.GetSize(); r++){
        if(f != NULL || g->IsForRegion(design->regions_[r].regionId)){
          const double scaling = scale ? design->scale_design[d][r] : 1.0;
          const double rscaling = 1.0 / scaling;
          const double translation = scale ? design->translate_design[d][r] : 0.0;
          const unsigned int u = base + design->regions_[r].base + design->regions_[r].elements;
          for(unsigned int i = base + design->regions_[r].base; i < u; i++)
          {
            DesignElement* de = &design->data[i];
            // standard or derivative case?
            if(derivative)
            {
              double val = 0.0;
              if(calculateTensorTrace){
                Matrix<double> material;
                GetErsatzMaterialTensor(material, de->elem, de->GetType());
                material.Trace(val);
                if(exponent != 1.0){
                  double des;
                  GetErsatzMaterialTensor(material, de->elem);
                  material.Trace(des);
                  val *= exponent * std::pow(des, exponent - 1.0);
                }
              }else{
                if(exponent != 1.0){
                  // by tf we mark if we want the physical stuff.
                  double des_val = tf != NULL ? tf->Derivative(de) : de->GetDesign(DesignElement::PLAIN);
                  val = exponent * std::pow(((des_val - translation) * rscaling), exponent - 1.0);
                }
                else
                  val = 1.0;
              }
              val *= rscaling;  // the gradient will be multiplied by scaling later again, but if it is calculated on the scaled designs here (as in regularization) this should not happen
              val *= fraction;
              if(!design->IsRegular()){
                grid->GetElemNodesCoord(cornerCoords, de->elem->connect, true);
                const double vol = de->elem->ptElem->CalcVolume(cornerCoords, false);
                val *= vol;
              }
              de->AddGradient(f, g, val);
            }
            else
            {
              // no derivative
              double des;
              if(calculateTensorTrace){ // use the trace of the stiffness Tensor as "volume"
                Matrix<double> material;
                GetErsatzMaterialTensor(material, de->elem);
                material.Trace(des);
              }
              else
              {
                // tf marks the physical function
                double des_val = tf != NULL ? tf->Transform(de) : de->GetDesign(DesignElement::PLAIN);
                des = (des_val - translation) * rscaling;
              }
              des = std::pow(des, exponent);
              if(design->IsRegular()){
                sum += des;
              }else{
                grid->GetElemNodesCoord(cornerCoords, de->elem->connect, true);
                const double vol = de->elem->ptElem->CalcVolume(cornerCoords, false);
                sum += des * vol;
              }
            } // if derivative
          } // for element (i)
        } // right region
      } // for region
      if(calculateTensorTrace && !derivative){  // the derivative has to be calculated for all designs but the tensor itself only once
        break;
      }
    } // if relevant
  } // for design
  return derivative ? -1.0 : sum * fraction;
}

double ErsatzMaterial::CalcVolume(Objective* f, Condition* g, bool derivative, bool normalized)
{
  // obly a constraint has a design type set (if it has)
  DesignElement::Type des  = g == NULL ? DesignElement::DEFAULT : g->design;

  // parameter is form the base function
  Function::Type func = f != NULL ? f->GetType() : g->GetType();
  double        exp  = f != NULL ? f->GetParameter() : g->GetParameter();

  switch(func)
  {
  case Function::VOLUME:
    return IntegrateDesignVariable(f, g, derivative, des, normalized, false, 1.0); // no scaling, exponent=1

  case Function::PENALIZED_VOLUME:
  case Function::REALVOLUME:
    return IntegrateDesignVariable(f, g, derivative, des, normalized, false, exp);

  case Function::GAP:
  {
    if(!derivative)
    {
      double vol = IntegrateDesignVariable(f, g, false, des, normalized, false, 1.0);
      double pen = IntegrateDesignVariable(f, g, false, des, normalized, false, exp);
      return vol - pen;
    }
    else
    {
      if(!design->IsRegular())
        throw Exception(Function::type.ToString(func) + " only implemented for regular grids");

      CalcRegularGapConstraint(f, g, des);
      return -1.0;
    }
  }

  default: assert(false);
  }
  return -1.0; // cannot happen due to assert
}

void ErsatzMaterial::CalcRegularGapConstraint(Objective* f, Condition* g, DesignElement::Type dt)
{
  assert(design->IsRegular());

  unsigned int des = design->FindDesign(dt);
  unsigned int ele = design->GetNumberOfElements();

  // exponent for penalized volume
  const double exp = f != NULL ? g->GetParameter() : g->GetParameter();

  for(unsigned int i = des * ele; i < (des + 1) * ele; i++)
  {
    DesignElement& de = design->data[i];

    // derivative of vol is 1
    // derivative of pen (x^p) is p*x^(p-1)
    double pen_grad = exp * std::pow(de.GetDesign(DesignElement::SMART), exp - 1.0);
    // normalize
    double grad = (1.0 - pen_grad) / (double) ele;

    de.AddGradient(f, g, grad);
  }
}


double ErsatzMaterial::CalcGlobalDynamicCompliance(Excitation& excite, Objective* f)
{
  //GLOBAL_DYNAMIC_COMPLIANCE
  // c = u^T conj(u) -> "A note on sensitivity analysis of linear dynamic systems with
  //                          harmonic excitation"; Jakob S. Jensen; June 22, 2007
  Vector<complex<double> >& u = forward.Get(excite)->GetComplexVector(Solution::RAW_VECTOR);
  assert(u.GetSize() != 0);

  complex<double> csp;
  u.Inner(u, csp);  // the inner product is sum over u_i * conj(u_i);
  assert(csp.imag() == 0);
  return csp.real() * excite.GetFactor(f);
}

double ErsatzMaterial::CalcCompliance(Excitation& excite, Objective* cost, Condition* g, bool derivative)
{
  // note for the derivative case gradients are summed up (with weights)
  double result = 0.0;
  if(derivative)
  {
    // calculate the compliance which is according to
    // "A 99 line topology optimization code written in Matlab"; O.Sigmund, 2001
    // -> dc/dx_e = -p * x_e ^(p-1) u_e^T k_0 u_e
    TransferFunction* tf = design->GetTransferFunction(DesignElement::DENSITY, MECH, true);
    double factor = -1.0 * excite.GetWeightedFactor(cost);

    CalcU1KU2(tf, forward.Get(excite)->elem[MECH], MECH, forward.Get(excite)->elem[MECH], NULL, factor, STANDARD, cost, g);
  }
  else
  {
    // compliance is easier computed using f^T u on nodes with force
    // to avoid any work for assembling force again, we simply calculate solution times rhs from the system
    Vector<double>& u = forward.Get(excite)->GetRealVector(Solution::RAW_VECTOR);
    Vector<double>& f = forward.Get(excite)->GetRealVector(Solution::RHS_VECTOR);
    double sp = 0.0;
    u.Inner(f, sp);
    result = sp * excite.GetFactor(cost);
    LOG_DBG(em) << "CalcCompliance(): result=" << result << " sp=" << sp << " u=" << u.ToString() << " f=" << f.ToString(1);
  }
  return result;
}

template <class T>
double ErsatzMaterial::CalcOutputObjective(Excitation& excite, Objective* cost)
{
  // The output is <l,u> where l is -1* rhs of the adjoint pde
  // Here the rhs of the adoint pde is not -1 -> hence we do -1*<l,u>
  //
  // We always take l from rhs and don't store it explicitly.
  // Note that one has to use the algsys RHS! The PDE RHS is still from the
  // forward simulation!
  Vector<T>& u = dynamic_cast<Vector<T> & >(*(forward.Get(excite)->GetVector(Solution::RAW_VECTOR)));
  Vector<T>& l = dynamic_cast<Vector<T> & >(*(adjoint.Get(excite)->GetVector(Solution::SEL_VECTOR)));
  assert(u.GetSize() == l.GetSize());

  LOG_DBG2(em) << "OUPTUT: adjoint sel (l): " << l.ToString(1);
  LOG_DBG2(em) << "OUPTUT: forward sol (u): " << u.ToString(0);

  double result = 0.0;

  switch(cost->GetType())
  {
    case Objective::OUTPUT:
    {
      // this is <l, u> which is for complex not really defined as it might be non-real
      T inner = u.Inner(l);
      result = ((complex<double>) inner).real();
      result *= excite.GetFactor(cost);
      LOG_DBG2(em) << "output <l,u>: " << inner << " * " << excite.GetFactor(cost) << " -> " << result;
      break;
    }

    case Objective::ABS_DYN_OUTPUT_SQUARED:
    {
      // |<u,l>|^2 = Re{<u,l>}^2 + Im{<u,l>}^2 
      T ul = u.Inner(l);
      double re_ul = ((complex<double>) ul).real();
      double im_ul = ((complex<double>) ul).imag();
      result = re_ul * re_ul + im_ul * im_ul;
      result *= excite.GetFactor(cost);
      LOG_DBG2(em) << "output |<u,l>|^2 = Re{<u,l>}^2 + Im{<u,l>}^2: " << re_ul << "^2  + " << im_ul << "^2 -> " << result;
      break;
    }


    case Objective::DYNAMIC_OUTPUT:
    case Objective::CONJUGATE_COMPLIANCE:
    {
      // this is <u,L conj(u)> and only defined for the harmonic case!
      if(!harmonic) throw Exception("'" + cost->type.ToString(cost->GetType()) + "' is only defined for harmonic!");

      // we loop over the vectors and do the scalar product by hand as we have
      // no diagonal matrix version of l
      for(unsigned int i = 0; i < u.GetSize(); i++ )
      {
        if(l[i] == 0.0) continue; // we skip this so we can make output for the real cases

        complex<double> u_val = (complex<double>) u[i];
        // make sure we have no penalization stuff!
        assert(std::abs(u_val) < 1e15);

        double sp = std::real(std::conj(u_val) * l[i] * u_val);
        LOG_DBG2(em) << "CalcObjective: " << std::conj(u_val) << " * " << l[i] << " * " << u_val << " -> " << sp;
        result += sp;
      }
      LOG_DBG2(em) << "output <u,L u*>: " << result << " * " << excite.GetFactor(cost) << " -> " << result * excite.GetFactor(cost);
      result *= excite.GetFactor(cost);
      break;
    }

    default: EXCEPTION("Not handled");
  }

  return result;
}

double ErsatzMaterial::CalcAcouNearField(Excitation& excite, Objective* f)
{
  // see CalcOutputObjective() for more details
  // calc -1*Re{j*w*psi^T L grad_n psi^*}
  //    =    Im{  w*psi^T L grad_n psi^*}

  if(!harmonic) throw Exception("'acousticNearField' is only defined for harmonic!");

  Vector<complex<double> >&      psi = forward.Get(excite, Solutions::STANDARD)->GetComplexVector(Solution::RAW_VECTOR);
  Vector<complex<double> >& grad_psi = forward.Get(excite, Solutions::GRAD)->GetComplexVector(Solution::RAW_VECTOR);
  Vector<complex<double> >&      sel = adjoint.Get(excite, Solutions::STANDARD)->GetComplexVector(Solution::SEL_VECTOR);

  assert(psi.GetSize() == grad_psi.GetSize());
  assert(psi.GetSize() == sel.GetSize());

  LOG_DBG2(em) << "CANF: forward sol (psi): " << psi.ToString(1);
  LOG_DBG2(em) << "CANF: forward sol (grad_psi): " << grad_psi.ToString(1);
  LOG_DBG2(em) << "CANF: adjoint rhs (l): " << sel.ToString(1);

  double result = 0.0;

  // we loop over the vectors and do the scalar product by hand
  for(unsigned int i = 0; i < psi.GetSize(); i++ )
  {
    if(sel[i] == 0.0) continue; // is solution part of gamma_opt?

    double sp = std::imag(excite.GetOmega() * psi[i] * sel[i] * std::conj(grad_psi[i]));

    LOG_DBG2(em) << "CANF: " << i << ": Im{" << excite.GetOmega() << " * " << psi[i] << " * "
                 << sel[i] << " * " << std::conj(grad_psi[i]) << "} -> " << sp;
    result += sp;
  }

  LOG_DBG2(em) << "CANF: " << result << " * " << excite.GetFactor(f) << " -> " << result * excite.GetFactor(f);

  result *= excite.GetFactor(f);

  return result;
}

void ErsatzMaterial::SolveTrackingProblem(Excitation& excite, bool designelem, bool gridelem)
{
  Vector<double>& u = forward.Get(excite)->GetRealVector(Solution::RAW_VECTOR); // Tracking is only implemented for non-harmonic
  LOG_DBG3(em) << "SolveTrackingProblem: displacement vector: (" << u.GetSize() << ") " << u.ToString();

  shared_ptr<EqnMap> eqnMap = pde->GetEqnMap();
  Grid * grd = domain->GetGrid();
  MathParser * parser = domain->GetMathParser();
  unsigned int mathParserHandle = parser->GetNewHandle();
  CoordSystem * coosy = domain->GetCoordSystem();

  Vector<Double> rhs;
  assemble_->GetAlgSys()->GetRHSVal(rhs);

  // set rhs to 0
  rhs.Init();

  // set rhs for tracking nodes
  Vector<double> GlobCoord;
  for(unsigned int l=0; l < excite.trackings.GetSize(); l++){
    TrackingBc const & actTrack = *(excite.trackings[l]);
    EntityIterator it = actTrack.entities->GetIterator();
    const int dof = actTrack.dof;
    for(it.Begin(); !it.IsEnd(); it++){
      grd->GetNodeCoordinate(GlobCoord, it.GetNode());
      parser->SetCoordinates(mathParserHandle, *coosy, GlobCoord);
      parser->SetExpr(mathParserHandle, actTrack.value);
      const double uref = parser->Eval(mathParserHandle);
      const int eqnr = eqnMap->GetNodeEqn(it.GetNode(), dof) - 1; // equation numbers are 1 based but vector u is 0 based
      if(eqnr>=0){
        rhs[eqnr] = u[eqnr] - uref;
        LOG_DBG2(em) << "SolveTrackingProblem: tracking setting RHS equation " << eqnr << " (Node: " << it.GetNode() << ", dof: " << dof << ") to " << rhs[eqnr];
      }
    }
  }

  // solve system
  assemble_->GetAlgSys()->InitRHS(rhs);
  // skipped "_tracking"
  assemble_->GetAlgSys()->Solve(CreateAdjointAnalysisIdNode());


  // save solution to tracking_ and restore original solution
  Vector<Double> tmpSol;
  assemble_->GetAlgSys()->GetSolutionVal(tmpSol);
  assert(tmpSol.GetSize() > 0);
  ((StdPDE*)pde)->SaveSolution(tmpSol.GetPointer(), tmpSol.GetSize());
  if(designelem)
    tracking_->Read(Solution::ELEMENT_VECTORS, pde, MECH);
  if(gridelem)
    tracking_->Read(Solution::GRIDELEM_VECTORS, pde, MECH);
  assert(tmpSol.GetSize() == forward.Get(excite)->GetRealVector(Solution::RAW_VECTOR).GetSize());
  forward.Get(excite)->Write(pde);
}

double ErsatzMaterial::CalcTracking(Excitation& excite, Objective* f, Condition* g, bool derivative, bool solveproblem)
{
  if(derivative){
    // calculate the tracking functional gradient, which is z^T k_i u,
    // where Kz = ut
    // where ut = M^T M (u-u0) = I_\Gamma (u-u0)
    if(solveproblem){
      SolveTrackingProblem(excite);
    }

    // calculate gradient z^T k_i u
    // note that in multiple excitations case (this is always now), we do sum this up
    TransferFunction* tf = design->GetTransferFunction(DesignElement::DENSITY, MECH, false);
    double val = CalcU1KU2(tf, tracking_->elem[MECH], MECH, forward.Get(excite)->elem[MECH], NULL, -1.0*excite.normalized_weight, STANDARD, f, g);
    return val;
  }else{
    // prepare for calculation
    Vector<double>& u = forward.Get(excite)->GetRealVector(Solution::RAW_VECTOR); // Tracking is only implemented for non-harmonic
    LOG_DBG3(em) << "CalcTracking: displacement vector: (" << u.GetSize() << ") " << u.ToString();

    shared_ptr<EqnMap> eqnMap = pde->GetEqnMap();
    Grid * grd = domain->GetGrid();
    MathParser * parser = domain->GetMathParser();
    unsigned int mathParserHandle = parser->GetNewHandle();
    CoordSystem * coosy = domain->GetCoordSystem();

    // the tracking functional is ut^T ut (ut as above)
    double val = 0;
    Vector<double> GlobCoord;
    for(unsigned int l=0; l < excite.trackings.GetSize(); l++){
      TrackingBc const & actTrack = *(excite.trackings[l]);
      EntityIterator it = actTrack.entities->GetIterator();
      const int dof = actTrack.dof;
      for(it.Begin(); !it.IsEnd(); it++){
        grd->GetNodeCoordinate(GlobCoord, it.GetNode(), true); // get Updated Coordinate
        parser->SetCoordinates(mathParserHandle, *coosy, GlobCoord);
        parser->SetExpr(mathParserHandle, actTrack.value);
        const double uref = parser->Eval(mathParserHandle);
        const int eqnr = eqnMap->GetNodeEqn(it.GetNode(), dof) - 1; // equation numbers are 1 based, the vector u is 0 based
        const double uact = eqnr>=0 ? u[eqnr] : 0.0;
        val += (uact - uref) * (uact - uref);
      }
    }
    return val/2;
  }
}

double ErsatzMaterial::CalcPoissonsRatioAndYoungsModulus(Objective* cost, Condition* g, bool derivative)
{
  assert(cost != NULL || g != NULL); // at least one of them must be given!
  
  bool poisson(false);
  if(cost != NULL)
  {
    assert(cost->GetType() == Objective::POISSONS_RATIO || cost->GetType() == Objective::YOUNGS_MODULUS);
    poisson = cost->GetType() == Objective::POISSONS_RATIO;
  }
  else
  {
    assert(g->GetType() == Condition::POISSONS_RATIO || g->GetType() == Condition::YOUNGS_MODULUS);
    poisson = g->GetType() == Condition::POISSONS_RATIO;
  }

  Matrix<double> hom_tensor = CalcHomogenizedTensor();

  double result = 0.0;

  if(derivative)
  {
    StdVector<double> dE11;
    CalcHomogenizedTensorEntry(make_pair(1,1), true, dE11);
    StdVector<double> dE12;
    CalcHomogenizedTensorEntry(make_pair(1,2), true, dE12);

    const double E11 = hom_tensor[0][0];
    const double E12 = hom_tensor[0][1];
    double grad(0.0);
    const double denom = (E11 + E12) * (E11 + E12);
    
    for(unsigned int e = 0, ne = design->GetNumberOfElements(); e < ne; e++)
    {      
      if(poisson) 
        grad = (dE12[e] * E11 - E12 * dE11[e]);
      else
      {
        grad  = (E11 * E11 + 2.0 * E11 * E12 + 3.0 * E12 * E12) * dE11[e];
        grad -= (4.0 * E11 * E12 + 2.0 * E12 * E12) * dE12[e];
      }
      
      grad /= denom;

      design->data[e].AddGradient(cost, g, grad);
    }
  }
  else
  {
    if(poisson)
      result = MechanicMaterial::CalcIsotropicPoissonsRatio(hom_tensor);
    else
      result = MechanicMaterial::CalcIsotropicYoungsModulus(hom_tensor);
  }

  return result;
}


void ErsatzMaterial::SetTestStrainMatrix(Matrix<double> &matrix, const Vector<double> &vec)
{
  assert(matrix.GetNumCols() == dim);
  assert(matrix.GetNumRows() == dim);

/*
 * suggested by Fabian S.
 *   matrix[0][0] = vec[0];
  matrix[1][1] = vec[1];
  matrix[0][1] = 0.5 * vec[5]; // Voigt notation!
  matrix[1][0] = 0.5 * vec[5]; // this is 0 because of the symmetry
  if(dim == 3)
  {
    matrix[2][2] = vec[2];
    matrix[1][2] = 0.5 * vec[3];
    matrix[2][1] = 0.5 * vec[3]; // symmetry again

    matrix[0][2] = 0.5 * vec[4];
    matrix[2][0] = 0.5 * vec[4]; // symmetry again
  }
*/

  matrix[0][0] = vec[0];
  matrix[1][1] = vec[1];
  matrix[0][1] = vec[5]; // voigt notation!
  matrix[1][0] = 0.0; // because of symmetry we need factor 0.5

  if(dim == 3)
  {
    matrix[2][2] = vec[2];
    matrix[1][2] = vec[3];
    matrix[2][1] = 0.0; // symmetry again

    matrix[0][2] = vec[4];
    matrix[2][0] = 0.0; // symmetry again
  }
}

Matrix<double> ErsatzMaterial::CalcHomogenizedTensor()
{
  const double cube_vol(domain->GetGrid()->CalcVolumeSpannedByNamedNodes());

  unsigned int ex_size = excitations.GetSize();
  assert((dim == 2 && ex_size == 3) || (dim == 3 && ex_size == 6));

  Matrix<double> test_strain_matrix_ij(dim, dim);
  Matrix<double> test_strain_matrix_kl(dim, dim);

  Matrix<double> result(ex_size, ex_size);
  result.Init();

  for(unsigned int ij = 0; ij < ex_size; ++ij)
  {
    SetTestStrainMatrix(test_strain_matrix_ij, excitations[ij].test_strain);
    StdVector<SingleVector*> &u1 = forward.Get(ij)->elem[MECH]; // equal to \chi^{ij}

    for(unsigned int kl = 0; kl < ex_size; ++kl)
    {
      if(ij > kl) // already computed this entry!
      {
        // by construction, the resulting matrix must be symmetric
        // so we can skip the calculation and do a simple assignment instead
        result[ij][kl] = result[kl][ij];
        continue;
      }

      SetTestStrainMatrix(test_strain_matrix_kl, excitations[kl].test_strain);
      StdVector<SingleVector*> &u2 = forward.Get(kl)->elem[MECH]; // equal to \chi^{kl}

      // loop over elements. In the gradient case not summed up
      for(int e = 0, ne = design->GetNumberOfElements(); e < ne; ++e)
      {
        DesignElement* de = &design->data[e];
        Vector<double>& u1_vec = dynamic_cast<Vector<double>& >(*u1[e]);
        Vector<double>& u2_vec = dynamic_cast<Vector<double>& >(*u2[e]);

        // prepare for calculation
        double p = CalcHomogenizedElementProduct(this, de, false, u1_vec, u2_vec, test_strain_matrix_ij, test_strain_matrix_kl);
        result[ij][kl] += p/cube_vol; // normalize for volume
      }
    } // end of kl loop
  } // end of ij loop

  // save e.g. for CommitIteration()
  homogenizedTensor.Assign(result, 1.0);

  return result;
}

void ErsatzMaterial::CalcHomogenizedTrackingGradient(const Matrix<double>& target, const Matrix<double>& hom, Objective* f, Condition* g)
{
  const double cube_vol(domain->GetGrid()->CalcVolumeSpannedByNamedNodes());

  // TODO Would be E^* - E^H if expression templates would work
  Matrix<double> diff_tensor;
  diff_tensor = target - hom;

  const unsigned int ex_size(excitations.GetSize());
  assert((dim == 2 && ex_size == 3) || (dim == 3 && ex_size == 6));

  Matrix<double> test_strain_matrix_ij(dim, dim);
  Matrix<double> test_strain_matrix_kl(dim, dim);

  // our derivative tensor
  Matrix<double> hom_tensor_deriv(ex_size, ex_size);
  hom_tensor_deriv.Init(); // we set and do not add - hence one init is enough

  // loop over elements.
  for(int e = 0, ne = design->GetNumberOfElements(); e < ne; ++e)
  {
    DesignElement* de = &design->data[e];

    for(unsigned int ij = 0; ij < ex_size; ++ij)
    {
      SetTestStrainMatrix(test_strain_matrix_ij, excitations[ij].test_strain);
      StdVector<SingleVector*> &u1 = forward.Get(ij)->elem[MECH]; // equal to \chi^{ij}
      Vector<double>& u1_vec = dynamic_cast<Vector<double>& >(*u1[e]);

      for(unsigned int kl = 0; kl < ex_size; ++kl)
      {
        if(ij > kl) // already computed this entry!
        {
          // by construction, the resulting matrix must be symmetric
          // so we can skip the calculation and do a simple assignment instead
          hom_tensor_deriv[ij][kl] = hom_tensor_deriv[kl][ij];
          continue;
        }

        SetTestStrainMatrix(test_strain_matrix_kl, excitations[kl].test_strain);
        StdVector<SingleVector*> &u2 = forward.Get(kl)->elem[MECH]; // equal to \chi^{kl}
        Vector<double>& u2_vec = dynamic_cast<Vector<double>& >(*u2[e]);

        // prepare for calculation
        double p = CalcHomogenizedElementProduct(this, de, true, u1_vec, u2_vec, test_strain_matrix_ij, test_strain_matrix_kl);
        hom_tensor_deriv[ij][kl] = p/cube_vol; // normalize for volume
      } // end of kl loop
    } // end of ij loop

    // hom_tensor_deriv is completely set.
    // (E^* - E^H) * - d(E^H)/d(rho_e) -> therefore the minus !
    double grad = -1.0 * diff_tensor.ScalarProduct(hom_tensor_deriv);

    de->AddGradient(f, g, grad);
  } // element loop
}


double ErsatzMaterial::CalcHomogenizedTensorConstraint(Condition* g, bool derivative)
{
  // me make use of the multi-purpose CalcHomogenizedTensorEntry()
  StdVector<double> grad;

  double result = 0.0;

  // we have the form Exx or Exx-Eyy or Exx-Eyy-Ezz-Ezz for isotropy constraints

  for(unsigned int i = 0; i < g->coords.GetSize(); i++)
  {
    double t = CalcHomogenizedTensorEntry(g->coords[i], derivative, grad);

    if(derivative)
    {
      for(int e = 0, ne = design->GetNumberOfElements(); e < ne; ++e)
        design->data[e].AddGradient(NULL, g, (i == 0 ? 1.0 : -1.0) * grad[e]);
    }
    else
    {
      result = (i == 0 ? t : result-t);

      homogenizedTensor[g->coords[i].first - 1][g->coords[i].second - 1] = t;
      // all tensors are symmetric. Makes reading easier!
      homogenizedTensor[g->coords[i].second - 1][g->coords[i].first - 1] = t;

      //LOG_DBG(em) << "CHTC: g=" << g->ToString() << " coord=" << i << " ["
      //    << g->coords[i].first << "-1][" << g->coords[i].second << "-1] = " << t;
    }
  }

  return result;
}

double ErsatzMaterial::CalcHomogenizedTensorEntry(const std::pair<int, int> entry, bool derivative, StdVector<double>& grad_out)
{
  const double cube_vol(domain->GetGrid()->CalcVolumeSpannedByNamedNodes());

  assert((dim == 2 && excitations.GetSize() == 3) || (dim == 3 && excitations.GetSize() == 6));

  Matrix<double> test_strain_matrix_ij(dim, dim);
  Matrix<double> test_strain_matrix_kl(dim, dim);

  const unsigned int ij = entry.first - 1;
  const unsigned int kl = entry.second - 1;

  SetTestStrainMatrix(test_strain_matrix_ij, excitations[ij].test_strain);
  StdVector<SingleVector*> &u1 = forward.Get(ij)->elem[MECH]; // equal to \chi^{ij}

  SetTestStrainMatrix(test_strain_matrix_kl, excitations[kl].test_strain);
  StdVector<SingleVector*> &u2 = forward.Get(kl)->elem[MECH]; // equal to \chi^{kl}

  double result = 0.0;

  if(derivative) grad_out.Resize(design->GetNumberOfElements());

  // loop over elements. In the gradient case not summed up
  for(int e = 0, ne = design->GetNumberOfElements(); e < ne; ++e)
  {
    DesignElement* de = &design->data[e];
    Vector<double>& u1_vec = dynamic_cast<Vector<double>& >(*u1[e]);
    Vector<double>& u2_vec = dynamic_cast<Vector<double>& >(*u2[e]);

    // prepare for calculation
    double p = CalcHomogenizedElementProduct(this, de, derivative, u1_vec, u2_vec, test_strain_matrix_ij, test_strain_matrix_kl);
    result += p/cube_vol; // normalize for volume

    if(derivative)
    {
      grad_out[e] = result;
      result = 0.0; // reset such that we do not sum up for the next case!
    }
  }

  return result; // in the non-derivative case this is the sum.
}


double ErsatzMaterial::CalcHomogenizedElementProduct(ErsatzMaterial* obj, DesignElement* de, bool derivative,
    Vector<double>& u1_vec, Vector<double>& u2_vec,
    Matrix<double>& test_strain_matrix_ij, Matrix<double>& test_strain_matrix_kl)
{
  assert(u1_vec.GetSize() == u2_vec.GetSize());

  const unsigned int dim_ = obj->dim;

  LOG_DBG3(em) << "CHEP: de=" << de->ToString() << " u1_vec=" << u1_vec.ToString() << " u2_vec=" << u2_vec.ToString();

  // TODO too much temporary matrices and vectors!

  // from the coordinates of this element we build a "test displacement" vector
  // u1(,2)_0. it contains linear strains which are assembled in the following lines
  // these strains are not unique! an arbitrary constant can be added without change

  // coordinates of "this" element
  Matrix<double> tmp_mat;
  // coordinates of current element, not updated lagrangian
  domain->GetGrid()->GetElemNodesCoord(tmp_mat, de->elem->connect, false);

  Matrix<double> u1_tmp;
  u1_tmp = test_strain_matrix_ij * tmp_mat;
  Matrix<double> u2_tmp;
  u2_tmp = test_strain_matrix_kl * tmp_mat;

  assert(u1_tmp.GetNumCols() == u2_tmp.GetNumCols());
  assert(u1_tmp.GetNumRows() == u2_tmp.GetNumRows());
  assert(u1_tmp.GetNumRows() == dim_);
  assert(u1_tmp.GetNumRows() * u1_tmp.GetNumCols() == u1_vec.GetSize());

  // u1_tmp, u2_tmp have to be transformed into vectors
  // 2D: from 2x4 to 8x1 on quad elems, 2x3 to 6x1 on triangles
  Vector<double> u1_0(u1_vec.GetSize());
  Vector<double> u2_0(u2_vec.GetSize()); // u1 and u2 have the same size
  for(unsigned int out = 0, cols = u1_tmp.GetNumCols(); out < cols; ++out)
  {
    for(unsigned int in = 0; in < dim_; ++in)
    {
      u1_0[out*dim_ + in] = u1_tmp[in][out]; // equal to \chi^{0(ij)}
      u2_0[out*dim_ + in] = u2_tmp[in][out]; // equal to \chi^{0(kl)}
    }
  }

  u1_0 -= u1_vec;
  u2_0 -= u2_vec;

  // reuse tmp_mat as elementK-Matrix
  // Matrix<double> k_mat;
  obj->SetElementK(de, MECH, &tmp_mat, STANDARD, derivative);

  Vector<double> mat_vec;
  mat_vec = tmp_mat * u1_0;
  double result = mat_vec * u2_0;

  LOG_DBG3(em) << "CHEP de=" << de->ToString() << " result=" << result;

  return result;
}

double ErsatzMaterial::CalcGreyness(Condition* g, bool derivative)
{
  double greyness = 0.0; // element greyness
  int counter = 0; // to make it sure for different design variables!

  double lb, ub, span, org_value, value, grad, eval;
  lb = ub = value = grad = eval = span = 0.0;

  // we have to divide the gradients by their relative volume = fraction
  double fraction = g->design == DesignElement::DEFAULT ?
                    design->data.GetSize() : design->GetNumberOfElements();

  // do we want the physical value?
  TransferFunction* tf = g->IsPhysical() ? design->GetTransferFunction(g->design, MECH) : NULL;

  // go over the complete design space to set gradients of other types to 0
  for(unsigned int i = 0; i < design->data.GetSize(); i++)
  {
    DesignElement* de = &design->data[i];
    bool relevant = g->design == DesignElement::DEFAULT || g->design == de->GetType();
    if(relevant)
    {
      lb = de->GetLowerBound();
      ub = de->GetUpperBound();
      span = ub-lb;
      if(g->IsPhysical())
        org_value = derivative ? tf->Derivative(de) : tf->Transform(de);
      else
        org_value = de->GetDesign(DesignElement::PLAIN);

      // We normalize for a design variable from [0;1]
      // this has minor effect on density [0.001;1] but is important
      // for polarization[-1;1]
      value = (org_value - lb) / span;
    }

    if(derivative)
    {
      if(relevant)
      {
        // standard greyness without parameters! times 4 so we have 1 for maximum greyness
        // Note that we transformed to [0,1]
        // f(x)=4*(1-x)x * 4 = 4*(-x^2+x+1)
        // f'(x)= 4*(-2x+1)
        grad = (-8.0 * value + 4.0)/span;
      }
      // divide by fraction
      // the gradient of non-relevant is 0
      grad = relevant ? grad / fraction : 0.0;

      // set 0.0 if not relevant
      design->data[i].AddGradient(NULL, g, grad);
    }
    else // not derivative but function evaluation
    {
       // we normalized the greyness to value from [0;1]
       if(relevant)
       {
         eval = 4.0* (1-value) * (value);
         greyness += eval;
         counter++;
       }
    }
    LOG_DBG3(conditions) << Condition::type.ToString(g->GetType())
                         << " derive=" << derivative << " relevant=" << relevant
                         << " elem " << de->elem->elemNum << " des_value: " << de->GetDesign(DesignElement::PLAIN)
                         << " value = " << org_value
                         << " -> " << value << " grad=" << grad << " eval=" << eval
                         << " fraction=" << fraction << " counter=" << counter;

  }

  return greyness / (double) counter;
}

double ErsatzMaterial::CalcSlopeConstraint(Condition* g, bool derivative)
{
  SlopeCondition* slope = dynamic_cast<SlopeCondition*>(g);

  // take care, similar logic in SlopeCondition::GetSparsityPattern() !
  
  assert(slope->IsLocal());

  // as in Petersson and Sigmund, 1998 we have a constraint only for full X_P, Y_P (, Z_P)
  // neighborhood.

  unsigned int  de_idx     = slope->GetCurrentVirtualElement();
  DesignElement* de         = &(design->data[de_idx]);
  int            neigh_num  = slope->GetCurrentVirtualNeighbor(); // 0, 2 (,4)
  DesignElement* neigh      = de->vicinity->GetNeighbour((VicinityElement::Neighbour) neigh_num);

  // the abs(slope) is done by two inequality constraints. Therefore the 2 * dim
  // a) x_i+1 - x_i <= c*h
  // b) -x_i+1 + x_i <= c*h was x_i+1 - x_i >= -c*h
  bool case_a = slope->GetCurrentVirtualSign() == 1;

  if(derivative)
  {
    // the slope constraint is linear so the gradient is constant
    // furthermore it is dense.
    // We set it to the position where it belongs to

    // there are two entries per virtual constraint

    // the index of the neighbor element
    unsigned int neigh_idx = design->Find(neigh);

    design->data[neigh_idx].AddGradient(NULL, g, case_a ? 1.0 : -1.0);
    design->data[de_idx].AddGradient(NULL, g, case_a ? -1.0 : 1.0);

    LOG_DBG2(em) << "CSC: grad ce=" << slope->GetCurrentVirtualElement()
                 << " cn=" << slope->GetCurrentVirtualNeighbor()
                 << " cs=" << slope->GetCurrentVirtualSign()
                 << " ov=" << (case_a ? -1.0 : 1.0) << " nv=" << (case_a ? 1.0 : -1.0);
  }
  else
  {
    double own_val   = de->GetValue(DesignElement::DESIGN, DesignElement::PLAIN);
    double other_val = neigh->GetValue(DesignElement::DESIGN, DesignElement::PLAIN);

    double res = case_a ? other_val - own_val : -1.0 * other_val + own_val;

    LOG_DBG2(em) << "CSC: ce=" << slope->GetCurrentVirtualElement()
                 << " cn=" << slope->GetCurrentVirtualNeighbor()
                 << " cs=" << slope->GetCurrentVirtualSign()
                 << " ov=" << own_val << " nv=" << other_val << " -> " << res;
    return res;
  }
  return -1.0; // derivative
}


double ErsatzMaterial::CalcCheckerboard(Condition* g, bool derivative)
{
  assert(!derivative);

  // global result
  double sum = 0.0;
  int    div = 0; // divider to normalize sum

  // store as special result?
  int idx = design->GetSpecialResultIndex(DesignElement::DEFAULT, DesignElement::CHECKERBOARD);

  // loop over all elements
  unsigned int base  = design->FindDesign(g->design);
  unsigned int elems = design->GetNumberOfElements();
  for(unsigned int e = base * elems; e < (base + 1) * elems; e++)
  {
    DesignElement&   de  = design->data[e];
    VicinityElement* vic = de.vicinity;

    // this value
    double own = de.GetValue(DesignElement::DESIGN);
    double elem_max = 0.0;

    for(unsigned int d = 0; d < dim; d++)
    {
      VicinityElement::Neighbour left_idx  = (VicinityElement::Neighbour) (d*2 + 1); // VicinityElement::X_N, Y_N(, Z_N)
      VicinityElement::Neighbour right_idx = (VicinityElement::Neighbour) (d*2);    // VicinityElement::X_P, Y_P(, Z_P)

      if(vic->HasNeighbor(left_idx) && vic->HasNeighbor(right_idx))
      {
        double left  = vic->GetNeighbour(left_idx)->GetValue(DesignElement::DESIGN);
        double right = vic->GetNeighbour(right_idx)->GetValue(DesignElement::DESIGN);

        // "Heaviside"(rho_i - max( rho_i-1, rho_i+1)
        double smaller = std::max(0.0, own - std::max(left, right));

        // "Heaviside"(min( rho_i-1, rho_i+1) - rho_i)
        double larger = std::max(0.0, std::min(left, right) - own);

        elem_max = std::max(elem_max, larger + smaller);

        LOG_DBG2(em) << "CalcCheckerboard: e=" << de.elem << " dim=" << d << " own=" << own
                     << " left=" << left << " right=" << right << " smaller=" << smaller
                     << " larger=" << larger << " elem_max=" << elem_max
                     << " sum=" << sum;
      }

      if(idx > 0)
        design->data[e].specialResult[idx] = elem_max;

      sum += elem_max;
      div++;
    }
  }

  return sum / div;
}


void ErsatzMaterial::SolveStateProblem(Excitation* ev_only_exite)
{
  forward.GetMultiple()->Init();
  adjoint.GetMultiple()->Init();

  // if ev_only_exite is set we use the given excitation
  // -> it shall not coincide
  assert(!(ev_only_exite != NULL && multiple_excitation->IsEnabled()));

  // shall we normalize afterwards?
  bool normalize = false;

  for(unsigned int e = 0; e < excitations.GetSize(); e++)
  {
    Excitation* excite = ev_only_exite != NULL ? ev_only_exite : &excitations[e];
    excite->Apply();

    // this is true for all problem types
    Optimization::SolveStateProblem(excite);
    StorePDESolution(*excite, forward, true, false, "forward");

    // check our objectives if there is an adoint problem to solve
    for(unsigned int o = 0; o < objectives.data.GetSize(); o++)
    {
      switch(objectives.data[o]->GetType())
      {
        case Objective::COMPLIANCE:
        case Objective::TRACKING:
        case Objective::HOMOGENIZATION_TENSOR:
        case Objective::HOMOGENIZATION_TRACKING:
        case Objective::POISSONS_RATIO:
        case Objective::YOUNGS_MODULUS:
        case Objective::TYCHONOFF:
        case Objective::VOLUME:
        case Objective::PENALIZED_VOLUME:
        case Objective::GAP:
        case Objective::TEMPERATURE:
          break; // already done before switch

        case Objective::OUTPUT:
        case Objective::CONJUGATE_COMPLIANCE:
        case Objective::ABS_DYN_OUTPUT_SQUARED:
        case Objective::GLOBAL_DYNAMIC_COMPLIANCE:
        case Objective::DYNAMIC_OUTPUT:
        case Objective::ELEC_ENERGY:
          SolveAdjointProblem(*excite, objectives.data[o]);
          break;

        case Objective::ACOU_NEAR_FIELD:
          // sets forward(GRAD) RAW vector
          // and  adjoint(STANDARD) SEL vector
          PostProcAcouNearFieldSolution(*excite,
                                        forward.Get(*excite, Solutions::GRAD)->GetComplexVector(Solution::RAW_VECTOR),
                                        adjoint.Get(*excite)->GetComplexVector(Solution::SEL_VECTOR));
          SolveAcouNearFieldAdjointProblems(*excite);
          break;

        default:
          assert(false);
      }

      // when we do multiple excitations with adjusted weights we calculate the objective here
      // to find the best weights. CalcObjective is so cheap, it is done later again by
      // the optimizer but then the weights are set
      if(multiple_excitation->DoAdjustWeights())
      {
        // in the first iteration we adjust the weights
        // stride = 0 is only the first time
        if((multiple_excitation->stride < 1 && GetCurrentIteration() == 0) ||
            ((GetCurrentIteration() % multiple_excitation->stride) == 0))
        {
          excite->cost = CalcObjective(*excite, objectives.data[o]); // to be eventually normalized
          normalize = true;
        }
      }
    }
  }

  // we need only to normalize when the properties have changed. This also reflects the stride
  if(normalize)
    NormalizeMultipleExcitations();
}

void ErsatzMaterial::StorePDESolution(Excitation &excite, Solutions& solutions, bool read_rhs, bool save_sol, const std::string& comment, Solutions::Key key)
{
  Solution& sol = *(solutions.Get(excite, key));

  SingleVector* raw = NULL; // last result = every result for multiple solutions

  // store solution element wise for gradient and raw vector for objective.
  // This is redundant as currently the solution is the global one!
  for(map<Application, SinglePDE*>::iterator it = pdes.begin(); it != pdes.end(); ++it)
  {
    stringstream ss;
    ss << "StorePDESolution: prob=" << comment << " excite=" << excite.index << " pde: " << it->first;

    sol.Read(Solution::ELEMENT_VECTORS, it->second, it->first, save_sol);
    raw = sol.Read(Solution::RAW_VECTOR, it->second, it->first, save_sol);

    LOG_DBG2(em) << ss.str() << " sol: " << raw->ToString();

    if(read_rhs)
    {
      sol.Read(Solution::RHS_VECTOR, it->second, it->first, save_sol);
      LOG_DBG2(em) << ss.str() << " rhs: " << sol.GetVector(Solution::RHS_VECTOR)->ToString(1);
    }
  }

  // currently all pdes have the same raw vector, therefore we use the last one here.
  if(multiple_excitation->IsEnabled())
  {
    // check for very first call, then multiple has no size yet
    if(solutions.GetMultiple()->GetSize() == 0)
    {
      solutions.GetMultiple()->Resize(raw->GetSize());
      solutions.GetMultiple()->Init();
    }

    // add the current solution in a weighted sense to the multiple solution vector
    if(harmonic)
    {
      complex<double> w = complex<double>(excite.normalized_weight);
      dynamic_cast<Vector<Complex>*>(solutions.GetMultiple())->Add(w, dynamic_cast<Vector<Complex>&>(*raw));
    }
    else
    {
      double w = excite.normalized_weight;
      dynamic_cast<Vector<double>*>(solutions.GetMultiple())->Add(w, dynamic_cast<Vector<double>&>(*raw));
    }

    // TODO: This does not work any more since SingleVector::Add(Double Complex a, const SingleVector &vec) is no child of Vector<T>::Add(T a, const SingleVector &vec)
    //    solutions.multiple->Add(excite.normalized_weight, *raw);
    LOG_DBG2(em) << "StorePDESolution: multiple= " << solutions.GetMultiple()->ToString();
  }
}

void ErsatzMaterial::PostProcAcouNearFieldSolution(Excitation &excite, Vector<complex<double> >& grad_n, Vector<complex<double> >& sel)
{
  AcousticPDE* acou_pde = dynamic_cast<AcousticPDE*>(ToPDE(ACOUSTIC));

  // for all surface regions we get grad psi for the center points of the surface elements
  // the following has to be done:
  // -> calculate grad_n = grad * normal
  // -> find all assigned nodes
  // -> store these nodes for selection
  // -> average the grad_n values at these nodes from the neighbor elements

  StdVector<Vector<complex<double> > > grad;
  // apply surface normal to grad if vol2 is used.
  acou_pde->CalcGradSurfaceElement<complex<double> >(acou_near_field_elems_, ACOU_POTENTIAL, true, grad);

  // grad_n  will take the nodal grad_n psi as postproc forward result
  // temporary the se-node values are summed up and divided by the number of sums
  //Vector<complex<double> >& grad_n = forward.Get(excite, Solutions::GRAD)->GetComplexVector(Solution::RAW_VECTOR);
  grad_n.Resize(forward.Get(excite)->GetVector(Solution::RAW_VECTOR)->GetSize());
  grad_n.Init(0.0);

  // The selection vector for the objective. Not from output nodes but via surface regions
  // temporary we count the nodes from the surface elements to normalize grad_n. Set to ones in the end
  // Vector<complex<double> >& sel = adjoint.Get(excite)->GetComplexVector(Solution::SEL_VECTOR);
  // sel might be just created, then the size needs to be set
  sel.Resize(grad_n.GetSize());
  sel.Init(0.0);

  Vector<double> normal(dim);

  Result<complex<double> >& res = dynamic_cast<Result<complex<double> >& >(*acou_near_field_elems_);
  EntityIterator it = res.GetEntityList()->GetIterator();

  // loop over all surface elements
  for(it.Begin(); !it.IsEnd(); it++)
  {
    // determine normal
    const SurfElem* se = it.GetSurfElem();
    domain->GetGrid()->CalcSurfNormal(normal, *se);

    complex<double> grad_n_val = normal * grad[it.GetPos()];

    // also have the normal sign correction from CalcGradSurfaceElement()
    grad_n_val *= se->normalSign;

    const StdVector<unsigned int>& nodes = se->connect;

    for(unsigned int n = 0; n < nodes.GetSize(); n++)
    {
      const unsigned int node = nodes[n];
      grad_n[node] += grad_n_val;
      sel[node] += 1.0;
      LOG_DBG2(em) << "PPANFS: se=" << se->elemNum << " node=" << node << " sum=" << grad_n[node] << " count=" << sel[node] << " avg=" << (grad_n[node] / sel[node]);
    }
    LOG_DBG2(em) << "PPANFS: se=" << se->elemNum << " val=" << grad_n_val << " nodes=" << nodes.ToString();
  }

  // normalize and reset selection
  for(unsigned int i = 0; i < sel.GetSize(); i++)
  {
    assert((sel[i] == 0.0 && grad_n[i] == 0.0) || (sel[i] != 0.0 && grad_n[i] != 0.0));
    if(sel[i] != 0.0)
    {
      assert(sel[i].imag() == 0.0);
      grad_n[i] /= sel[i];
      sel[i] = 1.0; // normalize back
    }
  }

  LOG_DBG(em) << "PPANFS: sel=" << sel.ToString(1);
  LOG_DBG(em) << "PPANFS: grad=" << grad_n.ToString(1);
}



StdVector<pair<SinglePDE*, IdBcList> > ErsatzMaterial::SetHDBC()
{
  // IDBC values in the forward problem are homogeneous in the adjoint PDE. Consider elec and mech
  // store org value as idbc_list per pde
  StdVector<pair<SinglePDE*, IdBcList> > org_idbc;
  org_idbc.Reserve(pdes.size());

  for(map<Application, SinglePDE*>::iterator it = pdes.begin(); it != pdes.end(); ++it)
  {
    // we need a deep copy. Create new list for the PDE
    org_idbc.Push_back(make_pair(it->second, IdBcList()));

    // get the idbc list of the current pde
    IdBcList  idbc_list = it->second->GetIDBCList();
    for(unsigned int bc = 0; bc < idbc_list.GetSize(); bc++)
    {
      // store original value -> we have to do a deep copy!
      org_idbc.Last().second.Push_back(shared_ptr<InhomDirichletBc>(new InhomDirichletBc(*(idbc_list[bc]))));
      // make homogeneous values! -> this is only stored in the PDE, we have to call SetBCs() later!
      (*idbc_list[bc]).value = "0.0";
      (*idbc_list[bc]).phase = "0.0";
      LOG_DBG(em) << "Set IDBC to HDBC (" << it->second->GetName() << ") -> value "
                    << (*(org_idbc.Last().second[bc])).value << " -> " << (*(it->second->GetIDBCList()[bc])).value;
    }
  }

  return org_idbc;
}

void ErsatzMaterial::ResetHDBC(StdVector<pair<SinglePDE*, IdBcList> >& org_idbc)
{
  // reset the original IDBC which were HDBC for the adjoint PDE
  for(unsigned int p = 0; p < org_idbc.GetSize(); p++)
  {
    // get the idbc list of the current pde
    IdBcList idbc_list = org_idbc[p].first->GetIDBCList();
    IdBcList org_list  = org_idbc[p].second;
    assert(idbc_list.GetSize() == org_list.GetSize());
    for(unsigned int bc = 0; bc < idbc_list.GetSize(); bc++) {
      (*idbc_list[bc]).value = (*org_list[bc]).value;
      (*idbc_list[bc]).phase = (*org_list[bc]).phase;
      LOG_DBG(em) << "Reset HDBC (" << org_idbc[p].first->GetName() << ") -> " << (*(org_idbc[p].first->GetIDBCList()[bc])).value;
    }
    org_idbc[p].first->SetBCs();
  }
}

template <class T>
void ErsatzMaterial::SolveAdjointProblem(Excitation& excite, Objective* cost)
{
  // the forward problem was already solved and stored !!

  // Set the rhs and solve for it
  SetAndSolveAdjointRHS<T>(excite, cost);

  // store the stuff -> no rhs but special handling of element results
  StorePDESolution(excite, adjoint, false, true, "adjoint");

  // write back the solution s.th. CommitIteraion() makes StoreResults() properly.
  for(map<Application, SinglePDE*>::iterator it = pdes.begin(); it != pdes.end(); ++it)
    forward.Get(excite)->Write(it->second);
}

void ErsatzMaterial::SolveAcouNearFieldAdjointProblems(Excitation& excite)
{
  Vector<complex<double> >& psi      = forward.Get(excite)->GetComplexVector(Solution::RAW_VECTOR);
  Vector<complex<double> >& grad_psi = forward.Get(excite, Solutions::GRAD)->GetComplexVector(Solution::RAW_VECTOR);
  Vector<complex<double> >& sel      = adjoint.Get(excite)->GetComplexVector(Solution::SEL_VECTOR);
  assert(sel.NormMax() != 0.0); // must be set

  // first solve S lmbd_grad_n = - omega/(2*j) L psi
  // create first olas rhs vector
  Vector<complex<double> >& rhs_lmbd_grad = adjoint.Get(excite, Solutions::GRAD)->GetComplexVector(Solution::RHS_VECTOR);
  rhs_lmbd_grad.Resize(psi.GetSize()); // not set in the first call

  // any adjoint PDE has HDBC instead of IDBC -> will be undone later ResetHDBC()
  StdVector<pair<SinglePDE*, IdBcList> > org_idbc = SetHDBC();

  complex<double> factor = -1.0 * excite.GetOmega() / complex<double>(0.0, 2.0); // - omega/(2*j)
  for(unsigned int i = 0; i < rhs_lmbd_grad.GetSize(); i++)
  {
    rhs_lmbd_grad[i] = factor * sel[i] * psi[i];
  }

  LOG_DBG2(em) << "SACNFAP: rhs_lmbd_grad=" << rhs_lmbd_grad.ToString(1);

  // solve the problem. We get lmbd_grad_n
  assemble_->GetAlgSys()->InitRHS(rhs_lmbd_grad);
  assemble_->GetAlgSys()->Solve(CreateAdjointAnalysisIdNode()); // todo

  // store the stuff -> no rhs but special handling of element results
  StorePDESolution(excite, adjoint, false, true, "grad_lmbd_adjoint", Solutions::GRAD);

  // now calculate the second adjoint problem.
  // the second adjoint problem solves for standard lambda but requires the solution of the first
  // adjoint problem
  // S lmbd = omega/j L grad_n(psi) + S grad_n(lmbd_grad_n)

  // lmbd_grad_n is the solution of the first adjoint pde, grad_n(lmbd_grad_n) is found like grad_n(psi)
  Vector<complex<double> >& grad_lmbd_grad = adjoint.Get(excite, Solutions::GRAD)->GetComplexVector(Solution::RAW_VECTOR);
  Vector<complex<double> >  dummy_sel; // should become sel!
  PostProcAcouNearFieldSolution(excite, grad_lmbd_grad, dummy_sel);

  Vector<complex<double> >& rhs_lmbd = adjoint.Get(excite, Solutions::STANDARD)->GetComplexVector(Solution::RHS_VECTOR);

  factor = excite.GetOmega() / complex<double>(0.0, 1.0); // omega/j
  for(unsigned int i = 0; i < rhs_lmbd_grad.GetSize(); i++)
  {
    // TODO!!!! + S grad_n(lmbd_grad_n) is missing!
    rhs_lmbd[i] = factor * sel[i] * grad_psi[i];
  }

  // solve the problem. We get lmbd
  assemble_->GetAlgSys()->InitRHS(rhs_lmbd);
  assemble_->GetAlgSys()->Solve(CreateAdjointAnalysisIdNode()); // todo

  // store the stuff -> no rhs but special handling of element results
  StorePDESolution(excite, adjoint, false, true, "lmbd_adjoint", Solutions::STANDARD);

  ResetHDBC(org_idbc);
}

template <class T>
void ErsatzMaterial::SetAndSolveAdjointRHS(Excitation& excite, Objective* cost)
{
  // the adjoint RHS might be an output stuff, then the loads are changed.
  // save and restore them in any case.
  LoadList org_loads = assemble_->GetLoads();

  // set pseudo loads (if there are output nodes)
  ConstructSelection(excite);

  // any adjoint PDE has HDBC instead of IDBC -> will be undone later ResetHDBC()
  StdVector<pair<SinglePDE*, IdBcList> > org_idbc = SetHDBC();

  // set the adjoint rhs
  ConstructAdjointRHS(excite, cost);

  // make IDBC to HDBC -> moves the stuff to OLAS
  for(map<Application, SinglePDE*>::iterator it = pdes.begin(); it != pdes.end(); ++it)
    it->second->SetBCs();

  // calculate adjoint problem
  assemble_->GetAlgSys()->Solve(CreateAdjointAnalysisIdNode());

  ResetHDBC(org_idbc);

  // reset the original loads, the have been changed in the output case
  assemble_->SetLoads(org_loads);
}

void ErsatzMaterial::ConstructSelection(Excitation& excite)
{
  // overwrite the assemble loads with "pesudo loads"s loads
  assemble_->SetLoads(output_nodes_);

  // set our own RHS but delete first as Assemble adds
  assemble_->GetAlgSys()->InitRHS();

  // assemble the output nodes
  assemble_->AssembleRHSLoads();

  // save the "pseudo loading" which is the selection as the rhs for the adjoint
  // This is exactly what has been constructed. Not that for an adjoint RHS it needs
  // post processing.
  adjoint.Get(excite)->Read(Solution::SEL_VECTOR, pde);

  LOG_DBG2(em) << "ConstructSelection: excite=" << excite.index << " sel=" << adjoint.Get(excite)->GetVector(Solution::SEL_VECTOR)->ToString(1);
}


void ErsatzMaterial::ConstructRealAdjointRHS(Excitation& excite, Objective* cost)
{
  // this can only handle static output!
  assert(cost->GetType() == Objective::OUTPUT);

  Vector<double>& l = adjoint.Get(excite)->GetRealVector(Solution::SEL_VECTOR);
  // create a own OLAS vector

  Vector<double> rhs(l.GetSize());

  rhs = l * -1.0;

  assemble_->GetAlgSys()->InitRHS(rhs);
  assert(rhs.NormMax() != 0.0);
  LOG_DBG2(em) << "CARHS<double>: rhs before solving: " << rhs.ToString(1);
}

void ErsatzMaterial::ConstructComplexAdjointRHS(Excitation& excite, Objective* cost)
{
  // we handle only complex cases
  Vector<complex<double> >& u = forward.Get(excite)->GetComplexVector(Solution::RAW_VECTOR);
  Vector<complex<double> >& l = adjoint.Get(excite)->GetComplexVector(Solution::SEL_VECTOR);

  // create a OLAS vector
  Vector<complex<double> > rhs(u.GetSize());

  switch(cost->GetType())
  {
  case Objective::DYNAMIC_OUTPUT:       // rhs is from "output loads" and set in adjoint...rhs
    // the correct conjugate_output case is -L * u*, always complex!
    for(unsigned int i = 0; i < rhs.GetSize(); i++)
      rhs[i] = -1.0 * l[i] * std::conj(u[i]);
    break;

  case Objective::ABS_DYN_OUTPUT_SQUARED: // J = Re{<u,l>}^2 + Im{<u,l>}^2 .
  {
    // RHS = -0.5 * (2*<u_R,l> l - j 2*<u_I,l> l)
    //     = - <u_R,l> l + j <u_I,l> l
    //     = - Re{<u,l>} l + j Im{<u,l>} l
    // find the complex scalar product ul <u,l>
    complex<double> ul = u.Inner(l);
    LOG_DBG2(em) << "AdjustComplexAdjointRHS: <u,l> = " << ul;

    for(unsigned int i = 0; i < rhs.GetSize(); i++)
      rhs[i] = complex<double>(-1.0 * ul.real() * l[i].real(), ul.imag() * l[i].real());

  }
  case Objective::CONJUGATE_COMPLIANCE: // rhs is from original excitation, we stored it in forward...rhs
  {
    forward.Get(excite)->Read(Solution::RHS_VECTOR, pde); // set
    Vector<complex<double> >& org_rhs = forward.Get(excite)->GetComplexVector(Solution::RHS_VECTOR); // read
    assert(org_rhs.GetSize() == u.GetSize());

    // the actual rhs for the adjoint pde is org_rhs * conj(u) -> this is only stored in OLAS!!!
    for(unsigned int i = 0, n = rhs.GetSize(); i < n; i++)
      rhs[i] = -1.0 * org_rhs[i] * std::conj(u[i]);
    break;
  }


  case Objective::GLOBAL_DYNAMIC_COMPLIANCE:
    // S lambda = -conj(u);
    for(int i = 0, n = u.GetSize(); i < n; ++i)
      rhs[i] = -1.0 * std::conj(u[i]);
    break;

  default:
    assert(true); // e.g. for ELEC_ENERGY the rhs is set in PiezoSIMP::ConstructAdjointRHS()
  }

  assemble_->GetAlgSys()->InitRHS(rhs);
  assert(rhs.NormMax() != 0.0);
  LOG_DBG2(em) << "CARHS<complex>: rhs before solving: " << rhs.ToString(1);
}


void ErsatzMaterial::ConstructAdjointRHS(Excitation& excite, Objective* cost)
{
  // cannot be inlined due to linker problems
  if(harmonic) ConstructComplexAdjointRHS(excite, cost);
          else ConstructRealAdjointRHS(excite, cost);
}

void ErsatzMaterial::SetPDEs()
{
  // we either set MECH (almost always) or HEAT -> set the common pde attribute!
  pde = domain->GetSinglePDE("mechanic", false);
  if(pde != NULL)
    pdes[MECH] = pde;
  else
  {
    pde = domain->GetSinglePDE("heatConduction", true);
    pdes[HEAT] = pde;
  }

  // optionally also ACOUSTIC
  if(objectives.Has(Function::ACOU_NEAR_FIELD))
    pdes[ACOUSTIC] = domain->GetSinglePDE("acoustic", true);

  // ELEC is set in PiezoSIMP()
}



SinglePDE* ErsatzMaterial::ToPDE(Application app, bool throw_exception) const
{
  map<Application, SinglePDE*>::const_iterator it = pdes.find(app);
  if(it != pdes.end())
    return it->second;

  // nothing found
  if(throw_exception)
      EXCEPTION("No PDE '" << app << "' stored");

  return NULL;
}


Optimization::Application ErsatzMaterial::ToApp(SinglePDE* pde) const
{
  if(pde->GetName() == "electrostatic") return ELEC;
  if(pde->GetName() == "mechanic") return MECH;
  if(pde->GetName() == "heatConduction") return HEAT;
  if(pde->GetName() == "acoustic") return ACOUSTIC;

  throw Exception("invalid");
}


void ErsatzMaterial::GetErsatzMaterialTensor(Matrix<double>& mat, Elem* elem, DesignElement::Type direction){
  linElastInt* form = (linElastInt*)GetForm(elem->regionId, pde, pde, "linElastInt");
  form->calcDMat(mat, elem, direction);
}

ErsatzMaterial::Solutions::Solutions()
{
  this->em_ = NULL;
}

ErsatzMaterial::Solutions::~Solutions()
{
  for(map<Key, Unit*>::iterator it = data_.begin(); it != data_.end(); ++it)
    delete it->second;
}


void ErsatzMaterial::Solutions::Init(ErsatzMaterial* em)
{
  this->em_ = em;
}

ErsatzMaterial::Solutions::Unit::Unit()
{
  multiple = NULL;
}

ErsatzMaterial::Solutions::Unit::Unit(ErsatzMaterial* em)
{
  data.Resize(em->excitations.GetSize());
  for(unsigned int i=0; i < data.GetSize(); i++)
    data[i] = new Solution(em);

  if(em->IsHarmonic())
    multiple = new Vector<complex<double> >();
  else
    multiple = new Vector<double>();
}

ErsatzMaterial::Solutions::Unit::~Unit()
{
  for(unsigned int i = 0; i < data.GetSize(); i++)
  {
    delete data[i];
    data[i] = NULL;
  }

  delete multiple;
  multiple = NULL;
}

ErsatzMaterial::Solution* ErsatzMaterial::Solutions::Get(Excitation& excitation, Key key)
{
  return Get(excitation.index, key);
}

ErsatzMaterial::Solution* ErsatzMaterial::Solutions::Get(int excitation_index, Key key)
{
  map<Key, Unit*>::iterator it = data_.find(key);
  if(it != data_.end())
    return it->second->data[excitation_index];

  // we need to create the data for this key first
  data_[key] = new Unit(em_);
  return data_[key]->data[excitation_index];
}

SingleVector* ErsatzMaterial::Solutions::GetMultiple(Key key)
{
  map<Key, Unit*>::iterator it = data_.find(key);
  if(it != data_.end())
    return it->second->multiple;

  data_[key] = new Unit(em_);
  return data_[key]->multiple;
}


ErsatzMaterial::Solution::Solution(ErsatzMaterial* em)
{
  this->em_ = em;
  this->raw = NULL;
  this->rhs = NULL;
  this->select = NULL;
}

ErsatzMaterial::Solution::~Solution()
{
  delete raw;
  delete rhs;
  delete select;

  std::map<Application, StdVector<SingleVector* > >::iterator elem_iter;
  for(elem_iter = elem.begin(); elem_iter != elem.end(); elem_iter++)
  {
    StdVector<SingleVector* >& data = elem_iter->second;
    for(unsigned int i = 0; i < data.GetSize(); i++)
    {
      delete data[i];
    }
  }
}

template <class T>
void ErsatzMaterial::Solution::Write(StdPDE* pde)
{
  T* ptr = NULL;
  Vector<T>* v = static_cast<Vector<T>*>(raw);
  ptr = v->GetPointer();
  assert(ptr != NULL);
  assert(raw->GetSize() != 0);
  pde->SaveSolution(ptr, raw->GetSize());
}

void ErsatzMaterial::Solution::Write(StdPDE* pde, SingleVector* vec)
{
  // we are static
  bool harmonic = BasePDE::IsComplex(domain->GetDriver()->GetAnalysisType());

  if(harmonic)
  {
    Vector<complex<double> >& sol = dynamic_cast<Vector<complex<double> >& >(*vec);
    pde->SaveSolution(sol.GetPointer(), sol.GetSize());
  }
  else
  {
    Vector<double>& sol = dynamic_cast<Vector<double>& >(*vec);
    pde->SaveSolution(sol.GetPointer(), sol.GetSize());
  }
}

template <class T>
SingleVector* ErsatzMaterial::Solution::GetVector(StorageType st, bool create)
{
  switch(st)
  {
  case RAW_VECTOR:
    assert(raw != NULL || create);
    if(raw == NULL && create) // this will crash for !create in release
      raw = new Vector<T>();
    return raw;

  case RHS_VECTOR:
    assert(rhs != NULL || create);
    if(rhs == NULL && create)
      rhs = new Vector<T>();
    return rhs;

  case SEL_VECTOR:
    assert(select != NULL || create);
    if(select == NULL && create)
      select = new Vector<T>();
    return select;

  default:
    assert(false);
  }

  EXCEPTION("false");
}

SingleVector* ErsatzMaterial::Solution::GetVector(StorageType st)
{
  // as create is false, it makes no difference of we use the double or complex variant
  return GetVector<double>(st, false);
}

Vector<double>& ErsatzMaterial::Solution::GetRealVector(StorageType st)
{
  // we know what we want - hence we can create the result on the fly
  return dynamic_cast<Vector<double>& >(*GetVector<double>(st, true));
}

Vector<complex<double> >& ErsatzMaterial::Solution::GetComplexVector(StorageType st)
{
  return dynamic_cast<Vector<complex<double> >& >(*GetVector<complex<double> >(st, true));
}

template <class T>
SingleVector* ErsatzMaterial::Solution::Read(StorageType st, StdPDE* pde, Application app, bool save_sol)
{
  SolutionType solt;
  switch(app)
  {
  case NO_APP: // up to now
    solt = MECH_DISPLACEMENT;
    break;
  case MECH:
    solt = MECH_DISPLACEMENT;
    break;
  case ELEC:
    solt = ELEC_POTENTIAL;
    break;
  case HEAT:
    solt = HEAT_TEMPERATURE;
    break;
  case ACOUSTIC:
    solt = ACOU_POTENTIAL;
    break;
  default:
    EXCEPTION("Solution type not implemented");
  }

  shared_ptr<ResultInfo> resinfo = pde->GetResultInfo(solt);
  switch(st)
  {
    case GRIDELEM_VECTORS:
    {
      StdVector<SingleVector*>& elem_vec = gridelem[app];
      Grid* grd = domain->GetGrid();
      int n = grd->GetNumElems();
      if(elem_vec.GetSize() == 0){
        elem_vec.Resize(n);
        for(int ve = 0; ve < n; ve++){
          elem_vec[ve] = new Vector<T>;
        }
      }
      ElemList elemList(grd);
      for(int e = 0; e < n; e++){
        elemList.SetElement(grd->GetElem(e+1)); // GetElem is 1-based
        const EntityIterator& it = elemList.GetIterator();
        pde->GetSolVecOfElement((Vector<T>&) *elem_vec[e], it, resinfo);
      }
      return NULL;
    }
    case ELEMENT_VECTORS:
    {
      if(save_sol)
      {
        // store the solution in the PDE! This is necessary to read the element vector
        // in the adjoint case!
        Vector<T> tmpSol;
        em_->assemble_->GetAlgSys()->GetSolutionVal(tmpSol);
        pde->SaveSolution(tmpSol.GetPointer(), tmpSol.GetSize() );
      }

      // we save the element vectors in elem_vec. Might be empty the first call
      StdVector<SingleVector*>& elem_vec = elem[app];
      int n = em_->design->GetNumberOfElements();

      // check for first call
      if(elem_vec.GetSize() == 0)
      {
        elem_vec.Resize(n);
        for(int ve = 0; ve < n; ve++)
          elem_vec[ve] = new Vector<T>;
      }

      // create an element list to gain the iterator in the loop
      ElemList elemList(domain->GetGrid());

      // store the results in our own structure
      for(int e = 0; e < n; e++)
      {
        DesignElement* de = &em_->design->data[e];

        elemList.SetElement(de->elem);
        const EntityIterator& it = elemList.GetIterator();

        pde->GetSolVecOfElement((Vector<T>&) *elem_vec[e], it, resinfo);
      }
      return NULL;
    }
    case SEL_VECTOR: // interpreted as RHS_VECTOR
    case RAW_VECTOR:
    case RHS_VECTOR:
    {
      // It is best to get the RHS from the algebraic system (OLAS), the PDE
      // might not check about a further adjont calculation

      // get access to solution
      BaseSystem* bs = em_->assemble_->GetAlgSys();
      Vector<T> tmp;

      if(st == RAW_VECTOR)
        bs->GetSolutionVal(tmp);
      else
        bs->GetRHSVal(tmp); // RHS_VECTOR and SEL_VECTOR

      SingleVector* ptr = NULL;
      switch(st)
      {
      case RAW_VECTOR:
        if(raw == NULL) raw = new Vector<T>(tmp.GetSize());
        ptr = raw;
        break;
      case RHS_VECTOR:
        if(rhs == NULL) rhs = new Vector<T>(tmp.GetSize());
        ptr = rhs;
        break;
      case SEL_VECTOR:
        if(select == NULL) select = new Vector<T>(tmp.GetSize());
        ptr = select;
        break;
      default:
        assert(false);
      }

      *ptr = tmp; // copy constructor

      return ptr;
    }
  }

  throw Exception("false");
}

