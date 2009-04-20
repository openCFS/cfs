#include <string>
#include <cmath>
#include <limits>
#include <iomanip>

#include "Optimization/SIMP.hh"
#include "Optimization/DesignSpace.hh"
#include "Optimization/DesignElement.hh"
#include "Optimization/OptimizationMaterial.hh"
#include "Domain/domain.hh"
#include "Domain/surfElem.hh"
#include "General/exception.hh"
#include "Utils/basenodestoresol.hh"
#include "PDE/SinglePDE.hh"
#include "PDE/mechPDE.hh"
#include "Forms/baseForm.hh"
#include "Forms/linSurfForm.hh"
#include "CoupledPDE/DirectCoupledPDE.hh"
#include "Utils/StdVector.hh"
#include "MatVec/vector.hh"
#include "Utils/result.hh"
#include "Utils/mathParser/mathParser.hh"
#include "Driver/assemble.hh"
#include "Driver/baseSolveStep.hh"
#include "Driver/harmonicDriver.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/InfoNode.hh"
#include "DataInOut/programOptions.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "OLAS/algsys/basesystem.hh"
#include "OLAS/algsys/baseentrymanipulator.hh"
#include "MatVec/basematrix.hh"
#include "MatVec/stdmatrix.hh"


using namespace CoupledField;
using namespace std;

DECLARE_LOG(conditions)
DEFINE_LOG(conditions, "conditions")

DECLARE_LOG(em)
DEFINE_LOG(em, "ersatzMaterial")

Enum<ErsatzMaterial::Method> ErsatzMaterial::method;

ErsatzMaterial::ErsatzMaterial() : Optimization()
{
  /** We store here the solution */
  tracking_ = NULL;
  volume_fraction_ = 0.0;

  material = NULL; // to be set in PostInit()

  pn = param->Get("optimization")->Get("ersatzMaterial");

  method_ = method.Parse(pn->Get("method"));

  // region stuff
  StdVector<ParamNode*> region_list = pn->GetList("region");
  for(unsigned int i=0; i < region_list.GetSize(); i++){
    std::string reg = region_list[i]->Has("name") ? region_list[i]->Get("name")->AsString() : region_list[i]->AsString();
    if(!domain->GetGrid()->HasRegion(reg))
      throw Exception("region given in ersatzMaterial is invalid");
    regionIds.Push_back(domain->GetGrid()->RegionNameToId(reg));
  }

  // grid regular?
  assume_regular_grid_ = pn->Get("grid_regular")->AsBool();
  assume_constant_element_matrices_ = assume_regular_grid_ && method_ != ErsatzMaterial::PARAM_MAT;

  // we assume always to have a mechanic pde. in pde[] we might have more pdes.
  // historically "mech" is a pointer to pde[0]
  // when there is no more mech, then "mech" might be replaced by pde[0]
  mech = dynamic_cast<MechPDE*>(domain->GetSinglePDE("mechanic"));
  assert(pde.GetSize() == 0);
  pde.Push_back(mech);

  harmonic = BasePDE::IsComplex(mech->GetAnalysisType());

  // Get the assemble class
  assemble_ = mech->getPDE_assemble();

  // manipulate the PDE (the Integrator of its bilinear form) such that
  // it is set up with the new densitiy values in the optimization loops
  for(unsigned int r=0; r < regionIds.GetSize(); r++){
    GetForm(regionIds[r], mech, mech, "linElastInt")->SetSolDependent(true);
    if(harmonic)
      GetForm(regionIds[r], mech, mech, "MassInt")->SetSolDependent(true);
  }

  // set up the design space elements, note PiezoSIMP have only POLARIZATION
  // this includes the transfer functions!
  StdVector<ParamNode*> design_list = pn->GetList("design");
  StdVector<ParamNode*> transfer_list = pn->GetList("transferFunction");
  StdVector<ParamNode*> result = pn->GetList("result");
  design = new DesignSpace(regionIds, design_list, transfer_list, result, method_);
  // make basic loggings
  design->ToInfo(optInfoNode->Get(InfoNode::HEADER)->Get("designSpace"));

  // check our constraints, the shall have only valid designs
  for(unsigned int i = 0; i < constraints.GetSize();  i++)
    if(design->FindDesign(constraints[i].design, false) == -1)
      throw Exception("constraint " + constraints[i].ToString() + " operates on invalid design variable");
  for(unsigned int i = 0; i < outputs.GetSize();  i++)
    if(design->FindDesign(outputs[i].design, false) == -1)
      throw Exception("output constraint " + outputs[i].ToString() + " operates on invalid design variable");

  // now we know the constraints size -> PostInit design
  design->PostInit(constraints.GetSize() + outputs.GetSize());

  // give the domain this data, s.th. the ersatz material approch is appield
  domain->SetErsatzMaterial(design);

  // add optimization results to the pde
  design->AppendOptimizationResults(mech);

  // optionally write the densities to an xml file
  if(pn->Has("export"))
  {
    std::string name = pn->Get("export")->Get("file")->AsString();
    if(name == "[problem]") name = progOpts->GetSimName() + ".density.xml";
    exportDesign = CreateExportDesign(name, design_list, transfer_list);
    exportDesignAllIterations = pn->Get("export")->Get("save")->AsString() == "all";
  }
  else exportDesign = NULL;


  // read in tracking parameters from XML
  StdVector<ParamNode*> tracking_list = pn->GetList("tracking");
  for(unsigned int i = 0; i < tracking_list.GetSize(); i++){
    std::string name, dof, value;
    dof = "";
    tracking_list[i]->Get("name", name);
    tracking_list[i]->Get("dof", dof, false);
    tracking_list[i]->Get("value", value);
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
    tracks_.Push_back(actLoad);
  }

  // check for multiple loadcases (might be frequencies)
  PrepareMultipleExcitations();

  // forward and adjoint are initialized in PostInit()

  // it's cheap as well -> no multiple load cases yet!
  tracking_ = new Solution(this);

  // make a copy of the old iteration to calculate the move
  last_iteration.Resize(design->data.GetSize());

  // note the difference between function evaluations (line search) and iterations!
  last_evaluation.Resize(design->data.GetSize());
}

ErsatzMaterial::~ErsatzMaterial()
{
  delete tracking_;

  // if write to file close the xml envelope and the file
  if(exportDesign != NULL)
  {
    // a ToFile() shall not be necessary
    delete exportDesign; 
    exportDesign = NULL;
  }

  // "remove" the ersatzMaterial (=data) from the domain
  domain->SetErsatzMaterial(NULL);
}

void ErsatzMaterial::PostInit()
{
  // for simplicity we always store forward and adjoint and have a multiple object
  forward.Init(excitations.GetSize(), this);
  adjoint.Init(excitations.GetSize(), this); // it's cheap - even if not used

  std::string objective = "'" + objectiveType.ToString(cost->type) + "'";

  // checks
  switch(cost->type)
  {
  case OUTPUT:
  case COMPLIANCE:
  case TRACKING:
  case VOLUME:
    if(harmonic) throw Exception(objective + " is only for static state problems");
    break;

  case DYNAMIC_OUTPUT:
  case GLOBAL_DYNAMIC_COMPLIANCE:
    if(!harmonic) throw Exception(objective + " is only for harmonic state problems");
    break;

  default:
    break;
  }

  // actions
  // note, that Radiation has its bilinear form already set in MechPDE::DefineIntegrators()
  switch(cost->type)
  {
  case OUTPUT:
  case DYNAMIC_OUTPUT:
  {
    // optimizing for nodes is mechanic, optimizing for loads is electrostatic, even if we SIMP mechanic!
    ParamNode* output = cost->pn->Get("output");

    if(output->Has("nodes") && output->Has("loads"))
      throw Exception("current implementation cannot optimize nodes and loads concurrently");

    if(output->Has("displacement"))
      mech->ReadLoads(output->GetList("displacement"), output_nodes_);
    if(output->Has("elecPotential"))
      domain->GetSinglePDE("electrostatic")->ReadLoads(output->GetList("elecPotential"), output_nodes_);
    if(output->Has("acoustic"))
      domain->GetSinglePDE("acoustic")->ReadLoads(output->GetList("acoustic"), output_nodes_);

    if(output_nodes_.GetSize() == 0)
      throw Exception("no output optimization nodes/loads given");
    break;
  }
  // radiation is done in SIMP!
  default:
    break;
  }

  // create Material class
  OptimizationMaterial::System system = OptimizationMaterial::system.Parse(pn->Get("material"));
  switch(system)
  {
  case OptimizationMaterial::PIEZO:
    material = new OptPiezoMat(this);
    break;
    
  case OptimizationMaterial::MECHANIC:
    material = new OptMechMat(this);
    break;
  
  default: assert(false);
  }
  optInfoNode->Get(InfoNode::HEADER)->Get("material")->SetValue(OptimizationMaterial::system.ToString(material->GetSystem()));
  
  Optimization::PostInit();
}


void ErsatzMaterial::PrepareMultipleExcitations()
{
  // by definition, we always have at least one excitation.
  // This does not necessarily need to ba a load (but might be voltage, pressure, ...)
  excitations.Resize(1);
  excitations[0].index = 0;

  // the actual multipleExcitation description is read in Optimization as part of
  // objective function block

  int num_freq  = !harmonic ? 0 : dynamic_cast<HarmonicDriver*>(domain->GetDriver())->freqs.GetSize();
  int num_loads = assemble_->GetLoads().GetSize();

  // initialize data and do simple plausibilty check. Note that also 1 is "multiple"
  if(multiple_excitation->IsEnabled())
  {
    // we average the solutions(s) only for output.
    // In the calculations we average the individual calculations

    // decide if we have multiple loads or multiple frequencies
    // we connot do both!
    if(num_freq > 1 && num_loads > 1)
      throw Exception("Cannot to concurrent multiple excitations for multiple loads and multiple frequencies");

    // the following is validated by above and 1 frequency and 0 loads is harmless
    if(num_freq > 1)
      excitations.Resize(num_freq);
    if(num_loads > 1)
      excitations.Resize(num_loads);
  }

  double weight_sum = 0.0;


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
    LoadList loads = assemble_->GetLoads();

    MathParser * parser = domain->GetMathParser();
    unsigned int handle = parser->GetNewHandle();

    for(unsigned int i = 0; i < excitations.GetSize(); i++)
    {
      excitations[i].index = i;

      if(num_loads > (int) i)
      {
        excitations[i].load = loads[i];

        parser->SetExpr(handle, excitations[i].load->weight == "" ? "1.0" : excitations[i].load->weight);
        const double weight = parser->Eval(handle);
        excitations[i].weight = weight;

        weight_sum += weight;
      }
    }

    parser->ReleaseHandle(handle);
  }

  // calc the inital normalized_weight and print info. T

  // for many concurrent mechanical loads do no print information
  if(multiple_excitation->IsEnabled() || harmonic)
  {
    InfoNode* in = optInfoNode->Get(InfoNode::HEADER)->Get("exitations");

    if(!multiple_excitation->IsEnabled() && num_freq > 1 && optimizer_ != EVALUATE_INITIAL_DESIGN)
    {
      stringstream ss;
      ss << "Solve only for 1. frequency (" << excitations[0].frequency << "Hz) of "
         << num_freq << " as multiple excitations are disabled";
      in->Get(InfoNode::WARNING)->SetValue(ss.str());
    }

    // communicate what we have but also normalize the weights!
    for(unsigned int i = 0; i < excitations.GetSize(); i++)
    {
      Excitation& ex = excitations[i];

      ex.normalized_weight = ex.weight / weight_sum;

      InfoNode* exin = in->Get("excitation", InfoNode::APPEND);
      exin->Get("index")->SetValue(ex.index);
      if(ex.load)
        ex.load->ToInfo(exin->Get("load"));
      if(ex.frequency >= 0.0)
        exin->Get("frequency")->SetValue(ex.frequency);
    }
  }
}


void ErsatzMaterial::NormalizeMultipleExcitations()
{
  if(multiple_excitation->IsEnabled())
  {
    if(multiple_excitation->DoAdjustWeights())
    {
      // find best cost
      double best;
      if(cost->task == MAXIMIZE)
      {
        best = numeric_limits<double>::min();
        for(unsigned int i = 0; i < excitations.GetSize(); i++) best = max(best, excitations[i].cost);
      }
      else
      {
        best = numeric_limits<double>::max();
        for(unsigned int i = 0; i < excitations.GetSize(); i++) best = min(best, excitations[i].cost);
      }

      for(unsigned int i = 0; i < excitations.GetSize(); i++)
      {
        // this is for "best_value": We push (weight the gradient) the worst result most to equalize/over-/undercompensate
        // the algorithm is the foolowing: w_k^p J_k = const -> we define w_k' for J_k = J_best to be one.
        // By this we can compute all w_k'. The we sum up all w_k' to sum_w_k' and get w_k = w_k'/sum_w_k'
        double t = cost->task == MAXIMIZE ? best / excitations[i].cost : excitations[i].cost / best;
        t = std::pow(t, 1.0/multiple_excitation->damping); // fails for negative quotient with damping < 0
        excitations[i].weight = min(t, multiple_excitation->max_gain);;
        LOG_DBG(em) << "NormalizeMultipleExcitations: excitation[" << i << "] best: "
                    << best << " gain: " << t << " weight: " << excitations[i].weight;
      }
    }

    // normalize
    double weight_sum = 0.0;

    for(unsigned int i = 0; i < excitations.GetSize(); i++)
      weight_sum += excitations[i].weight;

    for(unsigned int i = 0; i < excitations.GetSize(); i++)
    {
      excitations[i].normalized_weight = excitations[i].weight / weight_sum;
      LOG_DBG(em) << "NormalizeMultipleExcitations: excitation[" << i << "].normalized_weight <- " << excitations[i].normalized_weight;
    }
  }

}


InfoNode* ErsatzMaterial::CreateExportDesign(const std::string& filename, StdVector<ParamNode*>& des, StdVector<ParamNode*>& tfs)
{
   InfoNode* in = new InfoNode(filename, "<?xml version=\"1.0\"?>");
   in->SetName("cfsErsatzMaterial");
   
   // write header
   InfoNode* in_ = in->Get("header");
   for(unsigned int i = 0; i < des.GetSize(); i++)
     in_->Get("dummy", InfoNode::APPEND)->SetValue(des[i]); // name is overwritten

   for(unsigned int i = 0; i < tfs.GetSize(); i++)
     in_->Get("dummy", InfoNode::APPEND)->SetValue(tfs[i]);

   return in;
   
}

void ErsatzMaterial::StoreResults(double step_val)
{
  // check if we have to play with the results, forward is default mode!
  if(commitMode_ == FORWARD || commitMode_ == BOTH)
  {
    // in case of FORWARD this is redundant but it hanldes multiple exitations ans is cheap

    // check for multiple to find out what is actually to be restores
    for(unsigned int i = 0; i < pde.GetSize(); i++)
    {
      if(multiple_excitation->IsEnabled())
        Solution::Write(pde[i], forward.multiple); // todo: later a pde dependcy is required
      else
        forward.data[0]->Write(pde[i], ToApp(pde[i])); // todo: currently always the same
    }
    // call real implementation in Optimization
    Optimization::StoreResults(step_val);
  }
  if(commitMode_ == ADJOINT || commitMode_ == BOTH)
  {
    // "forward" is not possible -> hence it is both_cases or adjoint
    if(adjoint.data[0]->raw.empty()) EXCEPTION("Ajdoint solution cannot be written on 'commit' because there is none");

    for(unsigned int i = 0; i < pde.GetSize(); i++)
    {
      if(multiple_excitation->IsEnabled())
        Solution::Write(pde[i], adjoint.multiple);
      else
        adjoint.data[0]->Write(pde[i], ToApp(pde[i]));
    }

    step_val = step_val == -1 ? currentIteration + 0.5 : step_val + 0.5;
    Optimization::StoreResults(step_val);
  }
}


InfoNode* ErsatzMaterial::CommitIteration(bool keep_iteration_number)
{
  // will write the cfs results and the log file
  InfoNode* iter = Optimization::CommitIteration(keep_iteration_number);
  // add our multiple excitation stuff here (only in info.xml, this would be to complex for dat
  if(multiple_excitation->IsEnabled())
  {
    for(unsigned int i = 0; i < excitations.GetSize(); i++)
    {
      InfoNode* e = iter->Get("excitation", InfoNode::APPEND);
      e->Get("index")->SetValue(excitations[i].index);
      e->Get("objective")->SetValue(excitations[i].cost);
      e->Get("normalized_weight")->SetValue(excitations[i].normalized_weight);
    }
  }

  if(exportDesign != NULL)
  {
    InfoNode* in = exportDesign->Get("set", exportDesignAllIterations ? InfoNode::APPEND : InfoNode::USE_EXISTING);
    // add the entry, note that the iteratation counter was incremented in base implementation
    in->Get("id")->SetValue(currentIteration-1);

    for(unsigned int i = 0; i < design->data.GetSize(); i++)
    {
      DesignElement* de = &design->data[i];
      InfoNode* el = in->Get("element", "nr", lexical_cast<string>(de->elem->elemNum));
      el->Get("type")->SetValue(DesignElement::type.ToString(de->GetType()));
      el->Get("design")->SetValue(de->GetDesign(DesignElement::PLAIN));
      el->Get("gradient")->SetValue(de->GetObjectiveGradient(DesignElement::PLAIN));
      el->Get("filt_grad")->SetValue(de->GetObjectiveGradient(DesignElement::SMART));
    }


    exportDesign->ToFile();
  }

  return iter;
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


BaseForm* ErsatzMaterial::GetForm(RegionIdType regionId, StdPDE* pde1, StdPDE* pde2, const std::string& integrator)
{
  return assemble_->GetBiLinForm(regionId, pde1, pde2, integrator)->GetIntegrator();
}

double ErsatzMaterial::CalcObjective()
{
  design->WriteDesignToExtern(last_evaluation.GetPointer()); // todo: what for? - Fabian

  double result = 0.0;

  // the multiple excitation case is a special case - for all other cases this is executed once
  for(unsigned int i = 0; i < excitations.GetSize(); i++)
  {
    Excitation& excite = excitations[i];

    excite.cost = CalcObjective(excite); // this is virtual!

    assert(excite.normalized_weight > 0.0);
    result += excite.cost * excite.normalized_weight;
  }

  cost->SetValue(result);
  return result;
}

double ErsatzMaterial::CalcObjective(Excitation& excite)
{
  if(harmonic) return CalcObjective<std::complex<double> >(excite);
          else return CalcObjective<double>(excite);
}


template <class T>
double ErsatzMaterial::CalcObjective(Excitation& excite)
{
  switch(cost->type)
  {
  case COMPLIANCE:
    return CalcCompliance(excite, false, NULL);

  case TRACKING:
    return CalcTracking(false, NULL);
    break;

  case VOLUME:
    return CalcVolume(false, NULL);

  case GLOBAL_DYNAMIC_COMPLIANCE:
    return CalcGlobalDynamicCompliance(excite);

  case OUTPUT:
  case DYNAMIC_OUTPUT:
  case CONJUGATE_COMPLIANCE:
    return CalcOutputObjective<T>(excite);

  default: throw Exception("objective no handled");
  }
}

void ErsatzMaterial::CalcObjectiveGradient(double* grad_out)
{
  // reset the cost gradients in the design elements and sum them up in a weighted way
  // to perform multiple loads
  design->Reset(DesignElement::COST_GRADIENT);

  // the multiple excitation case is a special case - for all other cases this is executed once
  for(unsigned int idx = 0; idx < excitations.GetSize(); idx++)
  {
    switch(cost->type)
    {
    // Note, that in SIMP case we handle compliance already  in SIMO
    case COMPLIANCE:
      // the multiple load cases implementation for the gradient is in SIMP
      CalcCompliance(excitations[idx], true, NULL);
      break;

    case TRACKING:
      // no multiple load cases for tracking implemented
      assert(excitations.GetSize() == 1);
      CalcTracking(true, NULL);
      break;

    case VOLUME:
      // does not depend on multipl load cases
      CalcVolume(true, NULL);
      break;

    default:
      CalcObjectiveGradient(excitations[idx]);
    }
  }

  if(grad_out != NULL)
    design->WriteGradientToExtern(grad_out, DesignElement::COST_GRADIENT, DesignElement::SMART);
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
                       SurfaceRef* rhs, double factor, CalcMode calcMode, Condition* constraint, int res_idx)
{
  LOG_DBG2(em) << "CalcU1KU2(): tf=" << tf->ToString() << " #u1=" << u1.GetSize()
                 << " app=" << application.ToString(app) << " #u2=" << u2.GetSize() << " calcMode="
                 << calcMode << " factor=" << factor << " rhs=" << (rhs == NULL ? "NULL" : rhs->ToString(1));

  // This solves <l,K'*u-f'> or <u1, K' * u2 - f'> for all elements and adds it up to the element gradients

  assert(u1.GetSize() != 0);
  assert(u1.GetSize() == u2.GetSize());

  double sum = 0.0;
  // the dimenions of our matrix is determinded by u1_vec and u2_vec.
  // mat will be filled by SetElementK where also the derivative form most cases is built in
  Matrix<T> mat(u1[0]->GetSize(), u2[0]->GetSize());
  Vector<T> mat_vec(u1[0]->GetSize());

  TransferFunction* rtf = rhs != NULL && rhs->valid ? design->GetTransferFunction(tf->GetDesign(), rhs->app) : NULL;

  // traverse over our lements
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

      double added_value = 0.0;
      if(constraint)
      { // do we calculate this for a constraint or objective derivative
        added_value = de->GetConstraintGradient(constraint) + this_value;
        de->SetConstraintGradient(constraint, added_value);
      }
      else
      {
        added_value = de->GetObjectiveGradient(DesignElement::PLAIN) + this_value;
        de->SetObjectiveGradient(added_value);
      }

      LOG_DBG3(em) << "<l,K'*u-f'> = " << sp << " -> " << this_value << " + " << (added_value - this_value) << " = " << added_value;

      sum += this_value;
      
      if(res_idx != -1) de->specialResult[res_idx] = this_value; 
    }
  }

  return sum;
}

template double ErsatzMaterial::CalcU1KU2<double>(TransferFunction* tf, StdVector<SingleVector*>& u1,
                       Application app, StdVector<SingleVector*>& u2,
                       SurfaceRef* rhs, double factor, CalcMode calcMode, Condition* constraint, int res_idx);

template double ErsatzMaterial::CalcU1KU2<std::complex<double> >(TransferFunction* tf, StdVector<SingleVector*>& u1,
                       Application app, StdVector<SingleVector*>& u2,
                       SurfaceRef* rhs, double factor, CalcMode calcMode, Condition* constraint,  int res_idx);

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

void ErsatzMaterial::CalcConstraintGradient(Condition* constraint, double* grad_out)
{
  CalcConstraint(constraint, true, grad_out);
}

double ErsatzMaterial::CalcConstraint(Condition* constraint)  
{
  return CalcConstraint(constraint, false); // no gradient
}

double ErsatzMaterial::CalcConstraint(Condition* constraint, bool derivative, double* grad_out)
{
  // Just for access of enums
  DesignElement de;

  // handle default parameter
  if(constraint == NULL)
  {
    if(constraints.GetSize() != 1)
      throw Exception("default constraint only valid with exactly one constraint");
    constraint = &constraints[0];
  }
  double result = 0.0;

  switch(constraint->GetName())
  {
    case Condition::VOLUME:
         result = CalcVolume(derivative, constraint);
         break;

    case Condition::COMPLIANCE:
         assert(excitations.GetSize() == 1);
         result = CalcCompliance(excitations[0], derivative, constraint);
         break;

    case Condition::GREYNESS:
         result = CalcGreyness(derivative, constraint);
         break;

    case Condition::GAUSS_GREYNESS:
         result = CalcGreyness(derivative, constraint);
         break;

    default: throw Exception("Constraint notimplemented", __FILE__, __LINE__);
  }

  // TODO: Hi Bastian, do you really need this ?
  // method_ == PARAM_MAT ? de.DEFAULT : DesignElement::DENSITY
  if(grad_out != NULL)
    design->WriteGradientToExtern(grad_out, de.CONSTRAINT_GRADIENT, de.SMART, constraint);


  LOG_TRACE2(em) << "CalcConstraint " << constraint->ToString()
                  << " derivative=" << derivative << " -> " << result;
  return result;
}

double ErsatzMaterial::CalcVolume(bool derivative, Condition* constraint)
{
  // when not assuming a regular grid, computation of Volume is not as simple
  // we also consider working only on a given region, when used as constraint
  Matrix<Double> cornerCoords;
  Grid* grid = domain->GetGrid();
  bool isObjective = constraint == NULL;
  double fraction = isObjective?volume_fraction_:constraint->volume_fraction_;
  // if we use optimization with DesignMaterial then
  // if no design is explicitly specified, the constraint works with the trace of the tensor
  // else it works with exactly the specified variable integrating it
  bool allDesignsRelevant = constraint == NULL || constraint->design == DesignElement::TENSOR_TRACE || constraint->design == DesignElement::DEFAULT;
  bool ersatzMaterialTensor = domain->HasErsatzMaterialTensor() && allDesignsRelevant;
  unsigned int upper = ersatzMaterialTensor ? design->GetNumberOfElements() : design->data.GetSize();
  // check whether we have to calculate the full volume
  if( fraction == 0.0 ){
    if(assume_regular_grid_){
     for(unsigned int i=0; i < upper; i++){
        DesignElement* de = &design->data[i];
        // if relevant
        if( (allDesignsRelevant || constraint->design == de->GetType()) && (isObjective || constraint->IsForRegion(de->elem->regionId)) ){
          fraction += 1.0;
        }
      }
      fraction = 1.0 / fraction;
    }else{
      for(unsigned int i=0; i < upper; i++){
        DesignElement* de = &design->data[i];
        // if relevant
        if( (allDesignsRelevant || constraint->design == de->GetType()) && (isObjective || constraint->IsForRegion(de->elem->regionId)) ){
          grid->GetElemNodesCoord(cornerCoords, de->elem->connect, true );
          fraction += de->elem->ptElem->CalcVolume(cornerCoords, false);
        }
      }
      fraction = 1.0 / fraction;
    }
    if(constraint == NULL){
      volume_fraction_ = fraction;
    }else{
      constraint->volume_fraction_ = fraction;
    }
  }

  double sum = 0.0;

  // loop over all elements to set 0.0 for not relevant gradients
  for(unsigned int i = 0; i < (derivative ? design->data.GetSize() : upper); i++)
  {
    DesignElement* de = &design->data[i];
    bool relevant = (allDesignsRelevant || constraint->design == de->GetType()) && (isObjective || constraint->IsForRegion(de->elem->regionId));

    if(derivative)
    {
      double val = 0.0;
      if(relevant){
        if(ersatzMaterialTensor){
          Matrix<double> material;
          GetErsatzMaterialTensor(material, de->elem, de->GetType());
          assert(material.GetNumRows() == material.GetNumCols());
          for(unsigned int i=0; i < material.GetNumRows(); i++){
            val += material(i,i);
          }
        }else{
          val = 1.0;
        }
        val *= fraction;
        if(!assume_regular_grid_){
          grid->GetElemNodesCoord(cornerCoords, de->elem->connect, true);
          const double vol = de->elem->ptElem->CalcVolume(cornerCoords, false);
          val *= vol;
        }
      }
      if(constraint == NULL){
        de->SetObjectiveGradient(val);
      }else{
        de->SetConstraintGradient(constraint, val);
      }
    }
    else
    {
      if(relevant)
      {
        double des = 0;
        if(ersatzMaterialTensor){ // use the trace of the stiffness Tensor as "volume"
          Matrix<double> material;
          GetErsatzMaterialTensor(material, de->elem);
          assert(material.GetNumRows() == material.GetNumCols());
          for(unsigned int i=0; i < material.GetNumRows(); i++){
            des += material(i,i);
          }
        }else{
          des = de->GetDesign(DesignElement::PLAIN);
        }
        if(assume_regular_grid_){
          sum += des;
        }else{
          grid->GetElemNodesCoord(cornerCoords, de->elem->connect, true);
          const double vol = de->elem->ptElem->CalcVolume(cornerCoords, false);
          sum += des * vol;
        }
      }
    }
  }

  return !derivative ? sum * fraction : -1.0;
}

double ErsatzMaterial::CalcGlobalDynamicCompliance(Excitation& excite)
{
  //GLOBAL_DYNAMIC_COMPLIANCE
  // c = u^T conj(u) -> "A note on sensitivity analysis of linear dynamic systems with
  //                          harmonic excitation"; Jakob S. Jensen; June 22, 2007
  SingleVector* cfsvec = forward.data[excite.index]->raw[MECH];
  assert(cfsvec != NULL);
  assert(cfsvec->GetSize() != 0);
  Vector<complex<double> >& u = dynamic_cast<Vector<complex<double> >& >(*cfsvec);
  assert(u.GetSize() != 0);

  complex<double> csp;
  u.Inner(u, csp);  // the inner product is sum over u_i * conj(u_i);
  assert(csp.imag() == 0);
  return csp.real() * excite.GetFactor(cost);
}

double ErsatzMaterial::CalcCompliance(Excitation& excite, bool derivative, Condition* constraint)
{
  // for the derivate case there is a another implementation in SIMP which handles multiple load cases.
  int idx = excite.index;
  double result = 0.0;
  if(derivative)
  {
    // calculate the compliance which is according to
    // "A 99 line topology optimization code written in Matlab"; O.Sigmund, 2001
    // -> dc/dx_e = -p * x_e ^(p-1) u_e^T k_0 u_e
    TransferFunction* tf = design->GetTransferFunction(DesignElement::DENSITY, MECH, true);
    double factor = -1.0 * excite.GetWeightedFactor(cost);
    
    CalcU1KU2(tf, forward.data[idx]->elem[MECH], MECH, forward.data[idx]->elem[MECH], NULL, factor, STANDARD, constraint);
  }
  else
  {
    // compliance is easier computed using f^T u on nodes with force
    // to avoid any work for assembling force again, we simply calculate solution times rhs from the system
    Vector<double>& u = dynamic_cast< Vector<double>& > (*forward.data[idx]->raw[MECH]);
    Vector<double>& f = dynamic_cast< Vector<double>& > (*forward.data[idx]->rhs[MECH]);
    double sp = 0.0;
    u.Inner(f, sp);
    result = sp * excite.GetFactor(cost);
    LOG_DBG(em) << "CalcCompliance(): result=" << result << " sp=" << sp << " u=" << u.ToString() << " f=" << f.ToString();
  }
  return result;
}

template <class T>
double ErsatzMaterial::CalcOutputObjective(Excitation& excite)
{
  int idx = excite.index;

  // The output is <l,u> where l is -1* rhs of the adjoint pde
  // Here the rhs of the adoint pde is not -1 -> hence we do -1*<l,u>
  //
  // We always take l from rhs and don't store it explicitly.
  // Note that one has to use the algsys RHS! The PDE RHS is still from the
  // forward simulation!
  Vector<T>& u = dynamic_cast<Vector<T> & >(*(forward.data[idx]->raw[MECH]));
  Vector<T>& l = dynamic_cast<Vector<T> & >(*(adjoint.data[idx]->rhs[MECH]));
  assert(u.GetSize() == l.GetSize());

  LOG_DBG2(em) << "OUPTUT: adjoint rhs (l): " << l.ToString();
  LOG_DBG2(em) << "OUPTUT: forward sol (u): " << u.ToString();

  double result = 0.0;

  if(cost->type == OUTPUT)
  {
    // this is <l, u> which is for complex not really defined as it might be non-real
    T inner;
    u.Inner(l, inner);
    result = ((complex<double>) inner).real();
    result *= excite.GetFactor(cost);
    LOG_DBG2(em) << "output <l,u>: " << inner << " * " << excite.GetFactor(cost) << " -> " << result;
  }
  else
  {
    // this is <u,L conj(u)> and only defined for the harmonic case!
    if(!harmonic) throw Exception("'conjugateOutput' is only defined for harmonic!");

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
  }

  return result;
}

double ErsatzMaterial::CalcTracking(bool derivative, Condition* constraint)
{
  // prepare for calculation
  Vector<double>& u = dynamic_cast< Vector<double>& > (*forward.data[0]->raw[MECH]); // Tracking is only implemented for non-harmonic
  LOG_DBG3(em) << "CalcTracking: displacement vector: (" << u.GetSize() << ") " << u.ToString();

  shared_ptr<EqnMap> eqnMap = mech->GetEqnMap();
  Grid * grd = domain->GetGrid();
  MathParser * parser = domain->GetMathParser();
  unsigned int mathParserHandle = parser->GetNewHandle();
  CoordSystem * coosy = domain->GetCoordSystem();

  if(derivative){
    // calculate the tracking functional gradient, which is z^T k_i u,
    // where Kz = ut
    // where ut = M^T M (u-u0) = I_\Gamma (u-u0)
    Vector<Double> rhs;
    assemble_->GetAlgSys()->GetRHSVal(rhs);

    // set rhs to 0
//    for(int i=0; i < size; i++){
//      rhs[i] = 0;
//    }
    
    rhs.Init();

    // set rhs for tracking nodes
    Vector<double> GlobCoord;
    for(unsigned int l=0; l < tracks_.GetSize(); l++){
      TrackingBc const & actTrack = *tracks_[l];
      EntityIterator it = actTrack.entities->GetIterator();
      const int dof = actTrack.dof;
      for(it.Begin(); !it.IsEnd(); it++){
        grd->GetNodeCoordinate(GlobCoord, it.GetNode());
        parser->SetCoordinates(mathParserHandle, *coosy, GlobCoord);
        parser->SetExpr(mathParserHandle, actTrack.value);
        const double uref = parser->Eval(mathParserHandle);
        const int eqnr = eqnMap->GetNodeEqn(it.GetNode(), dof) - 1; // equation numbers are 1 based but vector u is 0 based
        rhs[eqnr] = u[eqnr] - uref;
        LOG_DBG2(em) << "CalcTracking: tracking setting RHS equation " << eqnr << " to " << rhs[eqnr];
      }
    }

    // solve system
    assemble_->GetAlgSys()->InitRHS(rhs);
    assemble_->GetAlgSys()->Solve(GetSolveComment() + "_tracking");

    // save solution to tracking_ and restore original solution
    Vector<Double> tmpSol;
    assemble_->GetAlgSys()->GetSolutionVal(tmpSol);
    assert(tmpSol.GetSize() > 0);
    ((StdPDE*)mech)->SaveSolution(tmpSol.GetPointer(), tmpSol.GetSize());
    tracking_->Read(Solution::ELEMENT_VECTORS, mech, MECH);
    assert(tmpSol.GetSize()== forward.data[0]->raw[MECH]->GetSize());
    forward.data[0]->Write(mech, MECH);

    // Hi Bastian: U read ELEMENT_VECTORS but want to print RAW - Fabian
    // LOG_DBG3(em) << "CalcTracking: solution of subproblem " << tracking_->raw[MECH]->ToString();

    // calculate gradient z^T k_i u
    TransferFunction* tf = design->GetTransferFunction(DesignElement::DENSITY, MECH, false);
    //double val = CalcU1KU2(tf, tracking_->elem[MECH], MECH, forward.data[0]->elem[MECH], NULL, false, -1.0, constraint);
    double val = CalcU1KU2(tf, tracking_->elem[MECH], MECH, forward.data[0]->elem[MECH], NULL, -1.0, STANDARD, constraint);
    return val;
  }else{
    // the tracking functional is ut^T ut (ut as above)
    double val = 0;
    Vector<double> GlobCoord;
    for(unsigned int l=0; l < tracks_.GetSize(); l++){
      TrackingBc const & actTrack = *tracks_[l];
      EntityIterator it = actTrack.entities->GetIterator();
      const int dof = actTrack.dof;
      for(it.Begin(); !it.IsEnd(); it++){
        grd->GetNodeCoordinate(GlobCoord, it.GetNode());
        parser->SetCoordinates(mathParserHandle, *coosy, GlobCoord);
        parser->SetExpr(mathParserHandle, actTrack.value);
        const double uref = parser->Eval(mathParserHandle);
        const int eqnr = eqnMap->GetNodeEqn(it.GetNode(), dof) - 1; // equation numbers are 1 based, the vector u is 0 based
        val += (u[eqnr] - uref) * (u[eqnr] - uref);
      }
    }
    return val/2;
  }
}

double ErsatzMaterial::CalcGreyness(bool derivative, Condition* constraint)
{
  double greyness = 0.0; // element greyness
  int counter = 0; // to make it sure for different design variables!

  double lb, ub, span, org_value, value, grad, eval;
  lb = ub = value = grad = eval = span = 0.0;
  bool gauss = constraint->GetName() == Condition::GAUSS_GREYNESS;
  // only used in gauss case
  double h = constraint->parameter;

  // we have to divide the gradients by their relalive volume = fraction
  double fraction = constraint->design == DesignElement::DEFAULT ?
                    design->data.GetSize() : design->GetNumberOfElements();

  // go over the complemte design space to set gradients of other types to 0
  for(unsigned int i = 0; i < design->data.GetSize(); i++)
  {
    DesignElement* de = &design->data[i];
    bool relevant = constraint->design == DesignElement::DEFAULT || constraint->design == de->GetType();
    if(relevant)
    {
      lb = de->GetLowerBound();
      ub = de->GetUpperBound();
      span = ub-lb;
      org_value = design->data[i].GetDesign(DesignElement::PLAIN);

      if(gauss)
      {
        // We normalize for a design variable to [-1;1] !
        // nothing to do for polarization [-1:1] but an issue for density [0.001:1]
        value = (org_value - lb) * (2.0/span) - 1.0;
      }
      else
      {
        // We normalize for a design variable from [0;1]
        // this has minor effect on density [0.001;1] but is important
        // for polarization[-1;1]
        value = (org_value - lb) / span;
      }
    }

    if(derivative)
    {
      if(relevant && gauss)
      {
        // a gauss type greyness (own idea, but most probably not unique :))
        // y = is a normalized version of x for the range [-1:1] (see above!)
        // f(y)=(e^(-1 * (y^2)/h) - e^(-1*1/h))/(1-e^(-1*1/h))
        // f'(y)=-2*y/h * e^(-1 * (y^2)/h)/(1-e^(-1*1/h))
        grad = (2.0/span) * -2.0 * (value/h) * exp(-1.0 * value * value / h) / (1-exp(-1.0/h));
      }
      if(relevant && !gauss)
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
      design->data[i].SetConstraintGradient(constraint, grad);
    }
    else // not derivative but function evaluation
    {
       // we normalized the greyness to value from [0;1]
       if(relevant)
       {
         if(gauss) eval = (exp(-1.0 * value * value / h) - exp(-1.0/h))
                         / (1-exp(-1.0/h));
              else eval = 4.0* (1-value) * (value);
         greyness += eval;
         counter++;
       }
    }
    LOG_DBG3(conditions) << Condition::name.ToString(constraint->GetName())
                         << " derive=" << derivative << " relevant=" << relevant
                         << " elem " << de->elem->elemNum << " value: " << org_value
                         << " -> " << value << " grad=" << grad << " eval=" << eval
                         << " fraction=" << fraction << " counter=" << counter
                         << " h=" << h;

  }

  return greyness / (double) counter;
}


void ErsatzMaterial::SolveStateProblem(Excitation* ev_only_exite)
{
  forward.multiple->Init();
  adjoint.multiple->Init();

  // if ev_only_exite is set we use the given excitation
  // -> it shall not coincide
  assert(!(ev_only_exite != NULL && multiple_excitation->IsEnabled()));

  // shall we normalize afterwards?
  bool normalize = false;

  for(unsigned int i = 0; i < excitations.GetSize(); i++)
  {
    Excitation* excite = ev_only_exite != NULL ? ev_only_exite : &excitations[i];
    excite->Apply();

    // this is true for all problem types
    Optimization::SolveStateProblem(excite);
    StorePDESolution(*excite, forward, true, false, "forward");

    switch(cost->type)
    {
    case COMPLIANCE:
    case TRACKING:
      break; // already done before switch

    case OUTPUT:
    case CONJUGATE_COMPLIANCE:
    case GLOBAL_DYNAMIC_COMPLIANCE:
    case DYNAMIC_OUTPUT:
    case ELEC_ENERGY:
      SolveAdjointProblem(*excite);
      break;

    default:
      assert(false);
    }

    // when we do multiple excitations with adjustes weights we calculate the objective here
    // to find the best weights. CalcObjective is so cheap, it is done later again by
    // the optimizer but then the weights are set
    if(multiple_excitation->DoAdjustWeights())
    {
      // in the first iteration we adjust the weights
      // stride = 0 is only the first time
      if((multiple_excitation->stride < 1 && GetCurrentIteration() == 0) ||
          ((GetCurrentIteration() % multiple_excitation->stride) == 0))
      {
        excite->cost = CalcObjective(*excite); // to be eventually normalized 
        normalize = true;
      }
    }
  }

  // we need only to normalize when the properties have changed. This also reflects the stride
  if(normalize)
    NormalizeMultipleExcitations();
}

void ErsatzMaterial::StorePDESolution(Excitation &excite,  Solutions& solutions, bool read_rhs, bool save_sol, const std::string& comment)
{
  Solution& sol = *(solutions.data[excite.index]);

  SingleVector* raw = NULL; // last result = every result for multiple solutions

  // store solution element wise for gradient and raw vector for objective.
  // This is redundant as currently the solution is the global one!
  for(unsigned int i=0; i < pde.GetSize(); i++)
  {
    Application app = ToApp(pde[i]);
    stringstream ss;
    ss << "StorePDESolution: prob=" << comment << " excite=" << excite.index << " pde: " << i;

    sol.Read(Solution::ELEMENT_VECTORS, pde[i], app, save_sol);
    raw = sol.Read(Solution::RAW_VECTOR, pde[i], app, save_sol);

    LOG_DBG2(em) << ss.str() << " sol: " << raw->ToString();

    if(read_rhs)
    {
      sol.Read(Solution::RHS_VECTOR, pde[i], app, save_sol);
      LOG_DBG2(em) << ss.str() << " rhs: " << sol.GetVector(Solution::RHS_VECTOR, app)->ToString();
    }
  }

  // currently all pdes have the same raw vector, therefore we use the last one here.
  if(multiple_excitation->IsEnabled())
  {
    // check for very first call, then multiple has no size yet
    if(solutions.multiple->GetSize() == 0)
    {
      solutions.multiple->Resize(raw->GetSize());
      solutions.multiple->Init();
    }

    // add the current solution in a weighted sense to the multiple solution vector
    if(harmonic)
    {
      complex<double> w = complex<double>(excite.normalized_weight);
      dynamic_cast<Vector<Complex>*>(solutions.multiple)->Add(w, dynamic_cast<Vector<Complex>&>(*raw));
    }
    else
    {
      double w = excite.normalized_weight;
      dynamic_cast<Vector<double>*>(solutions.multiple)->Add(w, dynamic_cast<Vector<double>&>(*raw));
    }
    
    // TODO: This does not work any more since SingleVector::Add(Double Complex a, const SingleVector &vec) is no child of Vector<T>::Add(T a, const SingleVector &vec)
    //    solutions.multiple->Add(excite.normalized_weight, *raw);
    LOG_DBG2(em) << "StorePDESolution: multiple= " << solutions.multiple->ToString();
  }
}


StdVector<IdBcList> ErsatzMaterial::SetHDBC()
{
  // IDBC values in the forward problem are homogeneous in the adjoint PDE. Consider elec and mech
  // store org value as idbc_list per pde
  StdVector<IdBcList> org_idbc(pde.GetSize());
  for(unsigned int p = 0; p < pde.GetSize(); p++) {
    // get the idbc list of the current pde
    IdBcList idbc_list = pde[p]->GetIDBCList();
    for(unsigned int bc = 0; bc < idbc_list.GetSize(); bc++) {
      // store original value -> we have to do a deep copy!
      org_idbc[p].Push_back(shared_ptr<InhomDirichletBc>(new InhomDirichletBc(*(idbc_list[bc]))));
      // make homogeneous values! -> this is only stored in the PDE, we have to call SetBCs() later!
      (*idbc_list[bc]).value = "0.0";
      (*idbc_list[bc]).phase = "0.0";
      LOG_DBG(em) << "Set IDBC to HDBC (" << pde[p]->GetName() << ") -> value "
                    << (*(org_idbc[p][bc])).value << " -> " << (*(pde[p]->GetIDBCList()[bc])).value;
    }
  }

  return org_idbc;
}

void ErsatzMaterial::ResetHDBC(StdVector<IdBcList> org_idbc)
{
  // reset the original IDBC which were HDBC for the adjoint PDE
  for(unsigned int p = 0; p < pde.GetSize(); p++) {
    // get the idbc list of the current pde
    IdBcList idbc_list = pde[p]->GetIDBCList();
    IdBcList org_list  = org_idbc[p];
    assert(idbc_list.GetSize() == org_list.GetSize());
    for(unsigned int bc = 0; bc < idbc_list.GetSize(); bc++) {
      (*idbc_list[bc]).value = (*org_list[bc]).value;
      (*idbc_list[bc]).phase = (*org_list[bc]).phase;
      LOG_DBG(em) << "Reset HDBC (" << pde[p]->GetName() << ") -> " << (*(pde[p]->GetIDBCList()[bc])).value;
    }
    pde[p]->SetBCs();
  }
}

template <class T>
void ErsatzMaterial::SolveAdjointProblem(Excitation& excite)
{
  int idx = excite.index; // 0 is the value if not doing multiple excitations

  // the forward problem was already solved and stored !!

  // Set the rhs and solve for it
  SetAndSolveAdjointRHS<T>(excite);

  // store the stuff -> no rhs but special handling of element results
  StorePDESolution(excite, adjoint, false, true, "adjoint");

  // write back the solution s.th. CommitIteraion() makes StoreResults() properly.
  for(unsigned int i=0; i < pde.GetSize(); i++)
    forward.data[idx]->Write( pde[i], ToApp(pde[i]));
}



template <class T>
void ErsatzMaterial::SetAndSolveAdjointRHS(Excitation& excite)
{
  // the adjoint RHS might be an output stuff, then the loads are changed.
  // save and restore them in any case.
  LoadList org_loads = assemble_->GetLoads();
  ConstructAdjointRHS(excite);
  
  // any adjoint PDE has HDBC instead of IDBC -> will be undone later ResetHDBC()
  StdVector<IdBcList> org_idbc = SetHDBC();

  // post process the RHS in the algebraic system
  if(harmonic)
    AdjustComplexAdjointRHS(excite);

  // make IDBC to HDBC -> moves the stuff to OLAS
  for(unsigned int p = 0; p < pde.GetSize(); p++)
    pde[p]->SetBCs();

  // calculate adjoint problem
  assemble_->GetAlgSys()->Solve(GetSolveComment(&excite) + "_adjoint");

  ResetHDBC(org_idbc);

  // reset the original loads, the have been changed in the output case
  assemble_->SetLoads(org_loads);
}


template <class T>
void ErsatzMaterial::ConstructAdjointRHS(Excitation& excite)
{
  // overwrite the assemble loads with our adjoint rhs loads
  // Note, that for the harmonic case, this is not yet the rhs for the adjoint PDE as
  // the multiplication with conj(u) is missing!

  assemble_->SetLoads(output_nodes_);

  // set our own RHS but delete first as Assemble adds
  assemble_->GetAlgSys()->InitRHS();

  // assemble the output nodes
  assemble_->AssembleRHSLoads();

  // multiply with -1 as the correct RHS is -dJ/du with J objective and u solution
  // TODO kick this fucking OLAS nonsense

  // first read the RHS with missing -1. it has no idbc and penalization set yet
  // we need this w/o -1 for the objective calculation.
  adjoint.data[excite.index]->Read(Solution::RHS_VECTOR, mech, MECH);  // up to now the vector is for all PDEs
  Vector<T>& t = dynamic_cast<Vector<T>& >(*(adjoint.data[excite.index]->rhs[MECH]));
  // create a own OLAS vector
  Vector<T> rhs(t.GetSize());
  // set the OLAS vector with -1 loads

  // hande the excite factor
  LOG_DBG2(em) << "ConstructOutputRHS: excite=" << excite.index;

  for(unsigned int i = 0; i < t.GetSize(); i++)
    rhs[i] = -1.0 * t[i];

  // write the new rhs to the algsys. Note, that adjoint[excite.index]->rhs[MECH]) misses the -1
  assemble_->GetAlgSys()->InitRHS(rhs);

  LOG_DBG2(em) << "ConstructOutputRHS: pure output vector w/o idbc: " << adjoint.data[excite.index]->rhs[MECH]->ToString();
  
  return;
}

void ErsatzMaterial::ConstructAdjointRHS(Excitation& excite)
{
  // here in ErsatzMaterial this is outout stuff
  if(harmonic) ConstructAdjointRHS<std::complex<double> >(excite);
          else ConstructAdjointRHS<double>(excite);
}

void ErsatzMaterial::AdjustComplexAdjointRHS(Excitation& excite)
{
  // we handle only complex cases
  Vector<complex<double> >& u = dynamic_cast<Vector<complex<double> >& >(*forward.data[excite.index]->raw[MECH]);

  // create a OLAS vector - note, that this is still fucking 1-based!! :)
  // todo: work directly on OLAS memory once the f*@!! OLAS buffer is kicked out of windows
  Vector<complex<double> > rhs(u.GetSize());

  int idx = excite.index;

  switch(cost->type)
  {
  case DYNAMIC_OUTPUT:       // rhs is from "output loads" and set in adjoint...rhs
  {
    // the correct conjugate_output case is -L * u* -- always complex!
    // we already have -1 in the
    Vector<complex<double> >& out = *(adjoint.data[idx]->GetComplexPointer(Solution::RHS_VECTOR, MECH));

    for(unsigned int i = 0; i < rhs.GetSize(); i++)
      rhs[i] = -1.0 * out[i] * std::conj(u[i]);
    break;
  }

  case CONJUGATE_COMPLIANCE: // rhs is from original excitation, we store it in forward...rhs
  {
    forward.data[idx]->Read(Solution::RHS_VECTOR, mech, MECH); // set
    Vector<complex<double> >& org_rhs = *(forward.data[idx]->GetComplexPointer(Solution::RHS_VECTOR, MECH)); // read
    assert(org_rhs.GetSize() == u.GetSize());

    // the actual rhs for the adjoint pde is org_rhs * conj(u) -> this is only stored in OLAS!!!
    for(unsigned int i = 0, n = rhs.GetSize(); i < n; i++)
      rhs[i] = -1.0 * org_rhs[i] * std::conj(u[i]);
    break;
  }


  case GLOBAL_DYNAMIC_COMPLIANCE:
    // S lambda = -conj(u);
    for(int i = 0, n = u.GetSize(); i < n; ++i)
      rhs[i] = -1.0 * std::conj(u[i]);
    break;

  default:
    assert(true); // e.g. for ELEC_ENERGY the rhs is set in PiezoSIMP::ConstructAdjointRHS()
  }

  assemble_->GetAlgSys()->InitRHS(rhs);
  LOG_DBG2(em) << "AdjustAdjointRHS: rhs before solving: " << rhs.ToString();
}



SinglePDE* ErsatzMaterial::ToPDE(Optimization::Application app, bool throw_exception) const
{
  switch(app)
  {
  case MECH: return pde[0];
  case ELEC: return pde[1];
  default: if(throw_exception) throw Exception("cannot map to pde");
                          else return NULL;
  }
}


Optimization::Application ErsatzMaterial::ToApp(SinglePDE* pde) const
{
  if(pde->GetName() == "electrostatic") return ELEC;
  if(pde->GetName() == "mechanic") return MECH;
  EXCEPTION("pdetype unknown;")
}



void ErsatzMaterial::GetErsatzMaterialTensor(Matrix<double>& mat, Elem* elem, DesignElement::Type direction){
  linElastInt* form = (linElastInt*)GetForm(elem->regionId, mech, mech, "linElastInt");
  form->calcDMat(mat, elem, direction);
}

ErsatzMaterial::Solutions::Solutions()
{
  this->multiple = NULL;
}

ErsatzMaterial::Solutions::~Solutions()
{
  for(unsigned int i = 0; i < data.GetSize(); i++)
  {
    delete data[i];
    data[i] = NULL;
  }

  delete multiple;
  multiple = NULL;
}

void ErsatzMaterial::Solutions::Init(int size, ErsatzMaterial* em)
{
  data.Resize(size);
  for(int i=0; i < size; i++)
    data[i] = new Solution(em);

  if(em->IsHarmonic())
    multiple = new Vector<complex<double> >();
  else
    multiple = new Vector<double>();
}

ErsatzMaterial::Solution::Solution(ErsatzMaterial* em)
{
  this->em_ = em;
}

ErsatzMaterial::Solution::~Solution()
{
  std::map<Application, SingleVector* >::iterator vec_iter;
  for(vec_iter = raw.begin(); vec_iter != raw.end(); ++vec_iter)
    delete vec_iter->second;

  for(vec_iter = rhs.begin(); vec_iter != rhs.end(); ++vec_iter)
   delete vec_iter->second;

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
void ErsatzMaterial::Solution::Write(StdPDE* pde, Application app)
{
  T* ptr = NULL;
  Vector<T>* v = static_cast<Vector<T>*>(raw[app]);  
  ptr = v->GetPointer();
  assert(ptr != NULL);
  assert(raw[app]->GetSize() != 0);
  pde->SaveSolution(ptr, raw[app]->GetSize());


  //  void SinglePDE::SaveSolution( const Complex * ptSol, UInt size ) {
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


SingleVector* ErsatzMaterial::Solution::GetVector(StorageType st, Application app)
{
  SingleVector* ptr = st == RAW_VECTOR ? raw[app] : rhs[app];
  assert(ptr != NULL);

  return ptr;
}


Vector<double>* ErsatzMaterial::Solution::GetRealPointer(StorageType st, Application app)
{
  SingleVector* ptr = st == RAW_VECTOR ? raw[app] : rhs[app];
  assert(ptr != NULL);

  Vector<double>* result = dynamic_cast<Vector<double>*>(ptr);
  assert(result != NULL);

  return result;
}

Vector<complex<double> >* ErsatzMaterial::Solution::GetComplexPointer(StorageType st, Application app)
{
  SingleVector* ptr = st == RAW_VECTOR ? raw[app] : rhs[app];
  assert(ptr != NULL);

  Vector<complex<double> >* result = dynamic_cast<Vector<complex<double> >*>(ptr);
  assert(result != NULL);

  return result;
}

template <class T>
SingleVector* ErsatzMaterial::Solution::Read(StorageType st, StdPDE* pde, Application app, bool save_sol)
{
  switch(st)
  {
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

        SolutionType solt = app == MECH ? MECH_DISPLACEMENT : ELEC_POTENTIAL;
        pde->GetSolVecOfElement((Vector<T>&) *elem_vec[e], it, pde->GetResultInfo(solt));
      }

      return NULL;
    }
    case RAW_VECTOR:
    case RHS_VECTOR:
    {
      // It is best to get the RHS from the algebraic system (OLAS), the PDE
      // might not check about a further adjont calculation

      // get access to solution
      BaseSystem* bs = em_->assemble_->GetAlgSys();
      Vector<T> tmpSol;

      if( st == RAW_VECTOR ) 
      {
        bs->GetSolutionVal(tmpSol);
      } else {
        bs->GetRHSVal(tmpSol);
      }
      UInt size = tmpSol.GetSize();
      std::map<Application, SingleVector* >& map_ = st == RAW_VECTOR ? raw : rhs;
      // check for first call
      if(map_.find(app) == map_.end())
        map_[app] = new Vector<T>(size);

      Vector<T>* raw_vec = dynamic_cast<Vector<T>* >(map_[app]);

      // todo: do better when we interchange nicely with OLAS and not any further by pointers!
      for(UInt i = 0; i < size; i++)
        (*raw_vec)[i] = tmpSol[i];

      return raw_vec;
    }
    default:
      assert(false);
      return NULL;
  }
}

