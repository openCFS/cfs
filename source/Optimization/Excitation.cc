#include <assert.h>
#include <stddef.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <ostream>

#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/Logging/log.hpp"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"
#include "Domain/Domain.hh"
#include "Domain/ElemMapping/EntityLists.hh"
#include "Domain/Mesh/Grid.hh"
#include "Domain/Results/ResultInfo.hh"
#include "Driver/Assemble.hh"
#include "Driver/BaseDriver.hh"
#include "General/defs.hh"
#include "General/Environment.hh"
#include "General/Exception.hh"
#include "MatVec/Matrix.hh"
#include "Optimization/Excitation.hh"
#include "Optimization/Function.hh"
#include "Optimization/Objective.hh"
#include "Optimization/Optimization.hh"
#include "PDE/SinglePDE.hh"
#include "PDE/StdPDE.hh"
#include "PDE/BasePDE.hh"
#include "PDE/ElecPDE.hh"
#include "Utils/mathParser/mathParser.hh"

using namespace CoupledField;
using namespace std;

DECLARE_LOG(exlog)
DEFINE_LOG(exlog, "excite")

Enum<MultipleExcitation::Type>  MultipleExcitation::type;

MultipleExcitation::MultipleExcitation(bool multiple, PtrParamNode pn)
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

void MultipleExcitation::ToInfo(PtrParamNode in) const
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

Excitation* MultipleExcitation::GetExcitation(const std::string& label, bool quiet)
{
  for(unsigned int i = 0; i < excitations.GetSize(); i++)
    if(excitations[i].label == label)
      return &(excitations[i]);

  if(quiet)
    return NULL;
  else
    throw Exception("None of the " + lexical_cast<string>(excitations.GetSize()) + " excitations has a label '" + label + "'");
}


void MultipleExcitation::PrepareMultipleExcitations(SinglePDE* pde, PtrParamNode optInfoNode, bool harmonic, bool eval_inital_design)
{
  // by definition, we always have at least one excitation.
  // This does not necessarily need to be a load (but might be voltage, pressure, ...)
  excitations.Resize(1);
  excitations[0].index = 0;

  PtrParamNode pncf = domain->GetParamRoot()->Get("optimization/costFunction");

  // the actual multipleExcitation description is read in Optimization as part of
  // objective function block

  int num_freq  = !harmonic ? 0 : dynamic_cast<HarmonicDriver*>(domain->GetDriver())->freqs.GetSize();

  pde->GetAssemble()->GetLinForms().GetSize();

  int num_loads = pde->GetAssemble()->GetLinForms().GetSize();  // to be faked later for homogenization test strains

  ParamNodeList pnexcitations;

  double weight_sum = 0.0;

  // initialize data and do simple plausibility check. Note that also 1 is "multiple"
  if(IsEnabled())
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
    if(DoHomogenization())
    {
      num_loads = SetHomogenizationTestStrains();
      weight_sum = 1; // all 0 but the first 1
    }

    // sets and resizes excitations with test charges
    if(DoMaxwellHomogenization())
    {
      num_loads = SetMaxwellHomogenizationTestCharges();
      weight_sum = 1; // all 0 but the first 1
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
      Excitation& ex = excitations[i];
      ex.index = i;
      if (DoMaxwellHomogenization()) {
        ex.frequency = hd->freqs[0].freq;
        ex.label = lexical_cast<string>(ex.frequency);
        //ex.weight    = hd->freqs[0].weight;
        ex.f_link    = &hd->freqs[0];
      }
      else {
        ex.frequency = hd->freqs[i].freq;
        ex.label = lexical_cast<string>(ex.frequency);
        ex.weight    = hd->freqs[i].weight;
        ex.f_link    = &hd->freqs[i];

        weight_sum += ex.weight;
      }

    }
  }
  if(!harmonic && IsEnabled()) // multiple loads case
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
        && !DoHomogenization()
        && !DoMaxwellHomogenization())
    {
      assert(false);
      /* FIXME
      LoadList loads = pde->GetAssemble()->GetLoads();

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
      */
    }

    parser->ReleaseHandle(handle);
  }

  // output summary
  // -------------------

  // finally verify that the labels are set
  for(unsigned int i = 0; i < excitations.GetSize(); i++)
    if(excitations[i].label == "")
      excitations[i].label = lexical_cast<string>(i);

  // calc the inital normalized_weight and print info.

  // for many concurrent mechanical loads do no print information
  if(IsEnabled() || harmonic)
  {
    PtrParamNode in = optInfoNode->Get(ParamNode::HEADER)->Get("excitations");

    if(!IsEnabled() && num_freq > 1 && !eval_inital_design)
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
      exin->Get("label")->SetValue(ex.label);
      assert(false);
/* FIXME
      if(! ex.loads.IsEmpty())
        ex.loads[0]->ToInfo(exin->Get("load"));
        */
      if(ex.frequency >= 0.0)
        exin->Get("frequency")->SetValue(ex.frequency);
      if(ex.test_strain.GetSize() > 0)
        exin->Get("testStrain")->SetValue(ex.test_strain.ToString());
      if(ex.test_charge.GetSize() > 0)
          exin->Get("testCharge")->SetValue(ex.test_charge.ToString());
    }
  }

}


int MultipleExcitation::SetHomogenizationTestStrains()
{
  // excitations.Resize(1);
  // Excitation& ex = excitations[0];
  // vec.Fill(ts[0], 6);
  // ex.ReadTestStrain(vec);
  // ex.weight = 1.0;
  // if(true) return 1;

  unsigned int dim = domain->GetGrid()->GetDim();

  int cases = dim == 2 ? 3 : 6;
  assert(excitations.GetSize() == 1);
  excitations.Resize(cases);

  assert((int) MechPDE::X == 0);
  assert((int) MechPDE::XY == 5);

  for(int i = 0, cnt = 0; i < 6; i++)
  {
    MechPDE::TestStrain ts = (MechPDE::TestStrain) i;
    Excitation& ex = excitations[cnt];

    // in 2D only 0, 1 and 5
    if(dim == 2 && (i == 2 || i == 3 || i == 4)) continue;

    ex.label = MechPDE::testStrain.ToString(ts);

    ex.ReadTestStrain(ts);
    // The homogenized tensor can only be evaluated for the last excitation!
    ex.weight = i < cases-1 ? 0.0 : 1.0;
    ++cnt;
  }

  return excitations.GetSize();
}

int MultipleExcitation::SetMaxwellHomogenizationTestCharges()
{
  double ts[3][3] =  { {1.0, 0.0, 0.0 },
                        {0.0, 1.0, 0.0 },
                        {0.0, 0.0, 1.0 } };

  Vector<double> vec;

  unsigned int dim = domain->GetGrid()->GetDim();

  int cases = dim == 2 ? 2 : 3;
  excitations.Resize(cases);

  for(int i = 0, cnt = 0; i < cases; ++i)
  {
    Excitation& ex = excitations[cnt];

    // in 2D only 0, and 1
    if(dim == 2 && (i == 2 )) continue;

    if(i == 0) ex.label = "x";
    if(i == 1) ex.label = "y";
    if(i == 2) ex.label = "z";

    vec.Fill(ts[i], dim);
 //   std::cout << vec.ToString() << std::endl;

    ex.ReadTestCharges(vec);
    // The homogenized tensor can only be evaluated for the last excitation!
    ex.weight = i < cases-1 ? 0.0 : 1.0;
    ++cnt;
  }

  return excitations.GetSize();
}



void MultipleExcitation::NormalizeMultipleExcitations(ObjectiveContainer* objectives)
{
  if(!IsEnabled()) return;

  const unsigned int exsi(excitations.GetSize());
  if(DoAdjustWeights())
  {
    // find best cost
    double best;
    if(objectives->DoMaximize())
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
      double t = objectives->DoMaximize() ? best / excitations[i].cost : excitations[i].cost / best;
      t = std::pow(t, 1.0/damping); // fails for negative quotient with damping < 0
      excitations[i].weight = min(t, max_gain);;
      LOG_DBG(exlog) << "NormalizeMultipleExcitations: excitation[" << i << "] best: "
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
    LOG_DBG(exlog) << "NormalizeMultipleExcitations: excitation[" << i << "].normalized_weight <- " << excitations[i].normalized_weight;
  }
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
  StdVector<LinearFormContext*>& assLinForms = domain->GetBasePDE()->GetAssemble()->GetLinForms();
  for(unsigned int i = 0; i < assLinForms.GetSize(); ++i){
    linForms->Push_back(assLinForms[i]);
  }
}

void Excitation::Apply()
{
  domain->GetOptimization()->applied_excitation = this;
  Assemble* a = domain->GetBasePDE()->GetAssemble();
  a->SetLinForms(linForms); // this works always, if not used just a copy of assembles linearForms
  if(! loads.IsEmpty()){
    assert(false);
    // FIXME
    // a->SetLoads(loads);
  }
  // a frequency cannot really be applied but has to be used as parameter
  // in the driver call
}

double Excitation::GetOmega() const
{
  if(frequency < 0.0)
    EXCEPTION("No frequency given");
  return 2 * M_PI * frequency;
}

double Excitation::GetFactor(Function* cost) const
{
  double factor = 1.0; // default

  if(cost->FactorOmegaOmega())
  {
    if(frequency < 0.0) EXCEPTION("No frequency given");
    factor *= 4 * M_PI * M_PI * frequency * frequency;
  }

  return factor;
}

double Excitation::GetWeightedFactor(Function* f) const
{
  if(f->DoEvaluateAlways())
   return normalized_weight * GetFactor(f);
  else
    return GetFactor(f);
}

void Excitation::ReadTrackings(PtrParamNode ts){
  assert(false);;
  /* FIXME
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
  */
}

void Excitation::ReadLoads(PtrParamNode ls)
{
  loads.Clear();

  SinglePDE* mech = domain->GetSinglePDE("mechanic");
  StdVector<std::string> comp;
  comp.Push_back("x");
  comp.Push_back("y");
  comp.Push_back("z");
  StdVector<PtrCoefFct > coef;
  bool update_geo;

  mech->ReadRhsExcitation("force", comp, ResultInfo::VECTOR, mech->IsComplex(), loads, coef, update_geo, ls);

  /* FIXME
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
  */
}

void Excitation::ReadTestStrain(MechPDE::TestStrain ts)
{
  assert(false);
  /* FIXME
  MechPDE* mech = dynamic_cast<MechPDE*>(domain->GetSinglePDE("mechanic"));

  loads.Clear();

  this->test_strain = mech->CalcTestStrainVector(ts);
  mech->DefineTestStrainIntegrator(ts, linForms);
  */
}

void Excitation::ReadTestCharges(const Vector<double>& vec)
{
  assert(false);
  /** FIXME
  this->test_charge = vec;

  loads.Clear();
  ElecPDE* elec = dynamic_cast<ElecPDE*>(domain->GetSinglePDE("electrostatic"));
  elec->SetRegionCharges(vec);
  elec->DefineMaxwellHomIntegrators(linForms);
  */
}


