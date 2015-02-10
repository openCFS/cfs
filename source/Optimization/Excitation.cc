#include <assert.h>
#include <stddef.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <ostream>
#include <boost/lexical_cast.hpp>

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
#include "Forms/LinForms/LinearForm.hh"
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


void MultipleExcitation::PrepareMultipleExcitations(Optimization* opt, bool eval_inital_design)
{
  Assemble* ass = domain->GetBasePDE()->GetAssemble();

  // by definition, we always have at least one excitation.
  // This does not necessarily need to be a load (but might be voltage, pressure, ...)
  excitations.Resize(1);
  excitations[0].index = 0;

  PtrParamNode pn = domain->GetParamRoot()->Get("optimization/costFunction");

  // the actual multipleExcitation description is read in Optimization as part of
  // objective function block

  int num_freq  = opt->IsHarmonic() ? dynamic_cast<HarmonicDriver*>(domain->GetDriver())->freqs.GetSize() : 0;

  int num_loads = ass->GetLinForms().GetSize();  // to be faked later for homogenization test strains
  LOG_DBG(exlog) << "PME: linForms from assemble: " << num_loads;

  // this might become explicitly given excite in the optimization xml description
  ParamNodeList pn_ex;

  double weight_sum = 0.0;

  // initialize data and do simple plausibility check. Note that also 1 is "multiple"
  if(IsEnabled())
  {
    // either every single load from bcsAndLoads is an excitation or allow combinations of loads, pressures, regionLoads
    // and trackings in one excitation when specified in multipleExcitation (only non-harmonic) (this is done here)
    if(!opt->IsHarmonic() && pn->Has("multipleExcitation/excitations"))
    {
      pn_ex = pn->Get("multipleExcitation/excitations")->GetChildren();
      num_loads = pn_ex.GetSize();
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

    // the following is validated by above and 1 frequency and 0 loads is harmless
    if(num_freq > 1)  excitations.Resize(num_freq);
    if(num_loads > 1) excitations.Resize(num_loads); // redundant for homogenization

    // we average the solutions(s) only for output.
    // In the calculations we average the individual calculations
  }


  // read in tracking parameters from XML, for the first and only load
  if(pn->Has("trackings")){
    excitations[0].ReadTrackings(pn->Get("trackings"));
  }

  // this sets the first and only excitation even when we have multiple harmonic forward case
  // but not multiple excitations. Then only the first frquency is called.
  // Fills the excitations list with the data provided in the xml file as problem description
  if(opt->IsHarmonic())
  {
    HarmonicDriver* hd = dynamic_cast<HarmonicDriver*>(domain->GetDriver());


    for(unsigned int i = 0; i < excitations.GetSize(); i++)
    {
      Excitation& ex = excitations[i];
      ex.index = i;
      ex.frequency = hd->freqs[i].freq;
      ex.label = lexical_cast<string>(ex.frequency);
      ex.weight    = hd->freqs[i].weight;
      ex.f_link    = &hd->freqs[i];

      weight_sum += ex.weight;
    }
  }
  if(!opt->IsHarmonic() && IsEnabled()) // multiple loads case
  {
    // when the loads are given in the optimization section of the xml file
    if(pn_ex.GetSize() > 0)
    {
      // this allows flexible weights, currently only when the loads are given in the optimization part and not for the simulation
      MathParser* parser = domain->GetMathParser();
      unsigned int handle = parser->GetNewHandle();

      // don't mix optimization excitations with bcsAndLoads stuff
      if(ass->GetLinForms().GetSize() > 0)
        opt->optInfoNode->Get(ParamNode::HEADER)->Get(ParamNode::WARNING)->
            SetValue("Excitations are given in optimization mulitload. "
            + boost::lexical_cast<std::string>(ass->GetLinForms().GetSize())
            + " loads from the bcsAndLoads section will be deleted.");

      ass->GetLinForms().Clear(); // the forms are gone!

      assert(excitations.GetSize() == pn_ex.GetSize());
      for(unsigned int i = 0; i < excitations.GetSize(); i++)
      {
        parser->SetExpr(handle, pn_ex[i]->Get("weight")->As<std::string>());
        const double weight = parser->Eval(handle);
        excitations[i].weight = weight;
        weight_sum += weight;

        if(pn_ex[i]->Has("trackings")){
          excitations[i].ReadTrackings(pn_ex[i]->Get("trackings"));
        }

        excitations[i].ReadLoads(pn_ex[i]); // possibly multiple forms for one excitation has it's own Resize(0)
      }

      parser->ReleaseHandle(handle);
    }
    // take the loads from the bcsAndLoads section. Store them here and reduce them to one.
    // Apply() will then change the entities for this loads
    if(pn_ex.GetSize() == 0 && !DoHomogenization())
    {
      // a force is a linarForm with integrator NodalForceInt
      StdVector<LinearFormContext*>& forms = ass->GetLinForms(true); // take the memory ownership

      for(unsigned int i = 0; i < excitations.GetSize(); i++) // via num_loads or num_freq
      {
        // don't do anything for multiple excitation enabled with only one load
        if(num_loads > (int) i)
        {
          // static load case
          excitations[i].forms.Push_back(forms[i]); // one form for one excitation
          LOG_DBG(exlog) << "PME: form " << i << " = " << forms[i]->GetIntegrator()->GetName() << " entities=" << forms[i]->GetEntities()->GetName();

          // we do not support to give a weight in the force section of the simulation
          excitations[i].weight = 1.0;
          weight_sum += excitations[i].weight;
        }
      }

      // "remove" the loads from the simulation. From now on we Apply() ist
      forms.Resize(0); // won't delete content but set the internal size_ counter
    }
  }

  // output summary
  // -------------------

  // finally verify that the labels are set
  for(unsigned int i = 0; i < excitations.GetSize(); i++)
    if(excitations[i].label == "")
      excitations[i].label = lexical_cast<string>(i);

  // calc the inital normalized_weight and print info.

  if(IsEnabled() || opt->IsHarmonic())
  {
    PtrParamNode in = opt->optInfoNode->Get(ParamNode::HEADER)->Get("excitations");

    if(!IsEnabled() && num_freq > 1 && !eval_inital_design)
    {
      stringstream ss;
      ss << "Solve only for 1. frequency (" << excitations[0].frequency << "Hz) of " << num_freq << " as multiple excitations are disabled";
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
      exin->Get("weight")->SetValue(ex.weight);
      exin->Get("normalized_weight")->SetValue(ex.normalized_weight);
      if(ex.forms.GetSize() == 1)
        exin->Get("load")->SetValue(ex.forms[0]->GetEntities()->GetName());
      if(ex.forms.GetSize() > 1)
      {
        if(ex.forms.GetSize() > 5)
          exin->Get("loads")->SetValue(ex.forms.GetSize());
        else
        {
          for(unsigned l = 0; l < ex.forms.GetSize(); l++)
            exin->Get(ex.forms[l]->GetIntegrator()->GetName(), ParamNode::APPEND)->Get("entities")->SetValue(ex.forms[l]->GetEntities()->GetName());
        }
      }
      if(ex.frequency >= 0.0)
        exin->Get("frequency")->SetValue(ex.frequency);
      if(ex.test_strain.GetSize() > 0)
        exin->Get("testStrain")->SetValue(ex.test_strain.ToString());
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

      LOG_DBG(exlog) << "NormalizeMultipleExcitations: excitation[" << i << "] best: " << best << " gain: " << t << " weight: " << excitations[i].weight;
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
  this->assemble = domain->GetBasePDE()->GetAssemble();
}

Excitation::~Excitation()
{
  // we have the ownerhip of the loads (form) with multiple excitation. We took if from Assemble::linForms_
  for(unsigned int i = 0; i < forms.GetSize(); i++)
    delete forms[i];
}

void Excitation::Apply()
{
  domain->GetOptimization()->applied_excitation = this;

  if(forms.GetSize() > 0)
  {
    assemble->GetLinForms() = forms; // let the copy constructor doe the stuff

    assert(assemble->GetLinForms().GetSize() == forms.GetSize());
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
  // the original loads need to be cleared before
  assert(assemble->GetLinForms().GetSize() == 0);

  // reads all loads and adds them to Assemble::linForms_
  dynamic_cast<SinglePDE*>(domain->GetBasePDE())->DefineRhsLoadIntegrators(ls);

  // own vector with the pointers in assemble and we have to delete the content
  forms = assemble->GetLinForms(true); // take ownership. copy constructor!
  assert(forms.GetSize() > 0);

  // delete the stuff from assemble to prepare for the next ReadLoads()
  assemble->GetLinForms().Resize(0); // the elements must not be destructed
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



