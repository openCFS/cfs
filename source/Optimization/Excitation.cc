/////#include <assert.h>
#include <stddef.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <ostream>
#include <boost/lexical_cast.hpp>

#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"
#include "Domain/CoefFunction/CoefFunctionCompound.hh"
#include "Domain/Domain.hh"
#include "Domain/ElemMapping/EntityLists.hh"
#include "Domain/Mesh/Grid.hh"
#include "Domain/Results/ResultInfo.hh"
#include "Driver/Assemble.hh"
#include "Driver/BaseDriver.hh"
#include "Driver/EigenFrequencyDriver.hh"
#include "Forms/BiLinForms/BDBInt.hh"
#include "General/defs.hh"
#include "General/Environment.hh"
#include "General/Exception.hh"
#include "MatVec/Matrix.hh"
#include "Optimization/Excitation.hh"
#include "Optimization/Function.hh"
#include "Optimization/Objective.hh"
#include "Optimization/Optimization.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "PDE/SinglePDE.hh"
#include "PDE/StdPDE.hh"
#include "PDE/BasePDE.hh"
#include "PDE/ElecPDE.hh"
#include "PDE/HeatPDE.hh"
#include "Forms/LinForms/LinearForm.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Utils/mathParser/mathParser.hh"

using namespace CoupledField;
using namespace std;

DEFINE_LOG(exlog, "excite")

Enum<MultipleExcitation::Type>  MultipleExcitation::type;

MultipleExcitation::MultipleExcitation(bool multiple, PtrParamNode pn)
{
  // set defaults
  this->stride = 1;
  this->damping = 1.0;
  this->max_gain = 1e4;
  this->multiple_excitation_ = multiple;
  this->type_ = NO_TYPE; // to be possibly overwritten soon
  this->num_trans_  = -1; // to be set later
  this->sequence_ = 1; // default
  this->principle_ = 0;

  // if disabled, we don't read anything
  if(multiple && pn != NULL)
  {
    this->type_ = type.Parse(pn->Get("type")->As<std::string>());
    this->sequence_ = pn->Get("sequence")->As<int>();
    // validate sequence
    assert(Optimization::manager.IsInitialized());
    if(sequence_ > (int) Optimization::manager.context.GetSize()) // 1-based!
      throw Exception("In 'multipleExcitation' the 1-based 'sequence' attribute has a too high value.");

    // adjust defaults
    if(pn->Has("adjustWeights")) {
      this->stride   = pn->Get("adjustWeights")->Get("stride")->As<Integer>();
      this->max_gain = pn->Get("adjustWeights")->Get("max_gain")->As<Double>();
      this->damping  = pn->Get("adjustWeights")->Get("damping")->As<Double>();
    }
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
  LOG_DBG(exlog) << "ME:GE: label: " << label;
  for(unsigned int i = 0; i < excitations.GetSize(); i++) {
    LOG_DBG(exlog) << "ME:GE: candidate: " << excitations[i].label;
    if(excitations[i].label == label)
      return &(excitations[i]);
  }
  if(quiet)
    return NULL;
  else
    throw Exception("None of the " + lexical_cast<string>(excitations.GetSize()) + " excitations has a label '" + label + "'");
}

unsigned int MultipleExcitation::GetNumberHomogenization(App::Type app) const
{
  // when we have only one sequence total_base_ must be the number of strains
  if(!DoHomogenization())
    return 0;
  else{
    UInt dim = domain->GetGrid()->GetDim();
    switch(app)
    {
      case App::MECH:
      case App::BUCKLING:
        return dim == 2 ? 3 : 6;
      case App::HEAT:
        return dim;
      default:
        throw Exception("Homogenization only implemented for heat and mech PDE! Requested: " + Optimization::context->ToPDE(app)->GetName());
    }
  }
}

bool MultipleExcitation::IsEnabled(int sequence) const
{
  assert(sequence == -1 || (sequence >= 1 && sequence <= (int) Optimization::manager.context.GetSize()));

  if(sequence == -1)
    return multiple_excitation_;
  // for bloch we do not set multiple excitation, it is implicit and we might use ME for homogenization
  if(Optimization::manager.context[sequence-1].DoBloch())
    return true;

  if(Optimization::manager.context[sequence-1].DoBuckling())
    return true;

  return multiple_excitation_ && sequence == sequence_;
}

void MultipleExcitation::WriteInInfo(PtrParamNode in)
{
  // communicate what we have but also normalize the weights!
  for(unsigned int i = 0; i < excitations.GetSize(); i++)
  {
    Excitation& ex = excitations[i];

    PtrParamNode exin = in->Get("excitation", ParamNode::APPEND);
    exin->Get("index")->SetValue(ex.index);
    exin->Get("meta_idx")->SetValue(ex.meta_index);
    if(Optimization::context->DoMultiSequence())
      exin->Get("sequence")->SetValue(ex.sequence);
    if(ex.transform != NULL)
      exin->Get("transform")->SetValue(ex.transform->ToString(1));
    if(ex.robust)
      exin->Get("robust_idx")->SetValue(ex.robust_filter_idx);

    exin->Get("label")->SetValue(ex.label);
    if(ex.transform || ex.robust)
      exin->Get("full_label")->SetValue(ex.GetFullLabel());

    exin->Get("weight")->SetValue(ex.weight);
    exin->Get("normalized_weight")->SetValue(ex.normalized_weight);
    if(ex.forms.GetSize() == 1)
      exin->Get("load")->SetValue(ex.forms[0]->GetEntities()->GetName());
    if(DoMetaExcitation(ex.sequence))
      exin->Get("reassemble")->SetValue(ex.reassemble);

    if(ex.forms.GetSize() > 1)
    {
      if (ex.forms.GetSize() > 5)
        exin->Get("loads")->SetValue(ex.forms.GetSize());
      else
      {
        for (unsigned l = 0; l < ex.forms.GetSize(); l++)
          exin->Get(ex.forms[l]->GetIntegrator()->GetName(), ParamNode::APPEND)->Get("entities")->SetValue(ex.forms[l]->GetEntities()->GetName());
      }
    }
    if(ex.frequency >= 0.0)
      exin->Get("frequency")->SetValue(ex.frequency);
  }

  domain->GetInfoRoot()->ToFile();
}

void MultipleExcitation::SetHarmonic(Context* ctxt, unsigned int base, int num_freq)
{
  // this sets the first and only excitation even when we have multiple harmonic forward case
  // but not multiple excitations. Then only the first frequency is called.
  // Fills the excitations list with the data provided in the xml file as problem description

  assert(excitations.GetCapacity() >= base + num_freq);
  excitations.Resize(base + num_freq);

  for (int i = 0; i < num_freq; i++)
  {
    Excitation& ex = excitations[base + i];
    SetHarmonicExcitation(ctxt, ex, i);
    ex.reassemble = true;
  }
}

void MultipleExcitation::SetHarmonicExcitation(Context* ctxt, Excitation& ex, int freq_idx)
{
  HarmonicDriver* hd = ctxt->GetHarmonicDriver();
  assert(freq_idx < (int) hd->freqs.GetSize());

  // note that we might have the case, that the excitation comes from a harmonic single freq multiple load case

  ex.frequency = hd->freqs[freq_idx].freq;
  assert(!(ex.label != "" && freq_idx > 0));
  if(ex.label == "") // don't overwrite the single frequenc multiple load case
    ex.label = lexical_cast<string>(freq_idx);
  ex.weight = hd->freqs[freq_idx].weight;
  ex.f_link = &hd->freqs[freq_idx];
}


void MultipleExcitation::SetBlochWaves(Context* ctxt, unsigned int base, int num_wave)
{
  const StdVector<Vector<double> >& wv =  ctxt->GetEigenFrequencyDriver()->wave_vectors;
  assert(excitations.GetCapacity() >= base + num_wave);
  excitations.Resize(base + num_wave);

  for (int i = 0; i < num_wave; i++)
  {
    Excitation& ex = excitations[base + i];
    ex.weight = i < num_wave - 1 ? 0 : 1; // we assume the slack variable as objective!
    ex.wave_vector = wv[i];
    ex.label = "(" + ex.wave_vector.ToString(TS_PLAIN,",") + ")";
  }
}


void MultipleExcitation::InitializeMultipleExcitations(Optimization* opt, ContextManager* manager)
{
  // this is to be called once prior to any PrepareMultipleExcitations() call

  // we have no "do_transform" either there is transformation or not
  num_trans_  = opt->GetDesign()->transform.GetSize();
  assert(!opt->GetDesign()->data.IsEmpty() && opt->GetDesign()->data[0].simp != NULL);
  int num_robust =  opt->GetDesign()->data[0].simp->filter.GetSize();

  robust_.Resize(manager->context.GetSize());

  // when robust and multi sequence we need to define what we want in the xml file
  if(num_robust > 1 && manager->context.GetSize() > 1)
  {
    if(!opt->optParamNode->Has("costFunction/multipleExcitation/robust"))
      throw Exception("when robust filters and multi sequence 'multipleExcitation/robust' elements are mandatory.");
    ParamNodeList list = opt->optParamNode->Get("costFunction/multipleExcitation")->GetList("robust");
    if(list.GetSize() != manager->context.GetSize())
      throw Exception("when robust filters and multi sequence you need a 'robust' in 'multipleExcitation' for every sequence.");
    int enable_count = 0;
    for(unsigned int i = 0; i < list.GetSize(); i++)
    {
      if(list[i]->Get("sequence")->As<unsigned int>() != i+1)
        throw Exception("the 'robust' elements in 'multipleExcitation' need ascending 'sequence'");
      if(list[i]->Get("enable")->As<bool>()) {
        robust_[i].num_robust = num_robust;
        enable_count++;
        if(list[i]->Has("filter"))
          throw Exception("in 'multipleExcitation/robust' use 'filter' only when not enabled");
      }
      else {
        if(!(list[i]->Has("filter")))
          throw Exception("in 'multipleExcitation/robust' when not enabled 'filter' is mandatory.");
        robust_[i].alt_filter = list[i]->Get("filter")->As<int>();
        if(robust_[i].alt_filter >= num_robust)
          throw Exception("'filter' value in 'multipleExcitation/robust' too large.");
      }
    }
    if(enable_count == 0)
      throw Exception("enable at least one sequence in 'multipleExcitation/robust'.");
  }
  // if robust and not multi sequence it is easy
  if(num_robust > 1 && manager->context.GetSize() == 1)
  {
    robust_[0].num_robust = num_robust;
    robust_[0].alt_filter = -1; // does not apply
    if(opt->optParamNode->Has("costFunction/multipleExcitation/robust"))
      opt->optInfoNode->Get(ParamNode::HEADER)->SetWarning("ignore 'multipleExcitation/robust' as we have only one sequence");
  }
  // if not robust or not multi sequence
  if(num_robust <= 1 && opt->optParamNode->Has("costFunction/multipleExcitation/robust"))
    opt->optInfoNode->Get(ParamNode::HEADER)->SetWarning("ignore 'multipleExcitation/robust' as we have no robust filters");

  // We do not want to have a copy of Excitations on excitation.Resize() therefore find an upper bound
  unsigned int upper_bound = 1;

  for(unsigned int i = 0; i < manager->context.GetSize(); i++)
  {
    Context& ctxt = manager->context[i];

    // no C++11 yet :(
    unsigned int tmp = std::max(ctxt.num_harm_freq, ctxt.num_bloch_wave_vectors);
    // we don't know the real number of loads, hence make an estimate. It's just to reserve!
    if(ctxt.analysis == BasePDE::STATIC)
      tmp = std::max(tmp, (unsigned int) 10); // FIXME we don't have Context::num_static_loads yet
    upper_bound += tmp;
  }

  upper_bound *= std::max(1,num_trans_);
  upper_bound *= GetNumberRobust(NULL, true); // any context, >= 1
  excitations.Reserve(upper_bound);
}

void MultipleExcitation::PrepareMultipleExcitations(Optimization* opt, Context* ctxt)
{
  // we assume InitializeMultipleExcitations() to be called
  assert(num_trans_ >= 0 && robust_[ctxt->context_idx].num_robust >= 0);

  Assemble* ass = ctxt->pde->GetAssemble();

  // the already existing excitations from a prior context
  unsigned int context_base = excitations.GetSize();
  assert((ctxt->context_idx == 0 &&  context_base == 0) || (ctxt->context_idx > 0 && context_base > 0)); // at least one per context

  // by definition, we always have at least one excitation.
  // This does not necessarily need to be a load (but might be voltage, pressure, ...)
  assert(excitations.GetCapacity() >= context_base + 1);
  excitations.Resize(context_base + 1);
  excitations[context_base].index = context_base;

  PtrParamNode pn = opt->optParamNode->Get("costFunction");

  // the actual multipleExcitation description is read in optimization as part
  // of objective function block
  int num_freq = ctxt->num_harm_freq; // might be 1 for multiple complex loads, eg. magnetic coupling

  // bloch mode analysis wave vectors
  int num_wave = ctxt->num_bloch_wave_vectors;

  // to be overwritten when we have multiple excitations in optimization or to be faked later for homogenization test strains
  int num_loads = ass->GetLinForms().GetSize();

  LOG_DBG(exlog) << "PME: linForms from assemble: " << num_loads << " trans: " << num_trans_;

  // this might become explicitly given excite in the optimization xml description
  ParamNodeList pn_ex;

  LOG_DBG(exlog) << "PME: cs=" << ctxt->sequence << " en=" << IsEnabled() << " seq=" << sequence_ << " cb=" << context_base << " cap=" << excitations.GetCapacity();

  // initialize data and do simple plausibility check. Note that also 1 is "multiple"
  // only for our sequence. We assume only one multipleExcitations element in xml even for multiple sequence optimization
  if(IsEnabled(ctxt->sequence))
  {
    if(sequence_ == ctxt->sequence)
    {
      // either every single load from bcsAndLoads is an excitation or allow combinations of loads, pressures, regionLoads
      // and trackings in one excitation when specified in multipleExcitation (only non-harmonic) (this is done here)
      if(!ctxt->IsHarmonic() && pn->Has("multipleExcitation/excitations"))
      {
        pn_ex = pn->Get("multipleExcitation/excitations")->GetChildren();
        num_loads = pn_ex.GetSize(); // in the magnetic coil case, we overwrite the num_loads from bcs
      }

      // decide if we have multiple loads or multiple frequencies
      // we cannot do both!

      if(num_freq > 1 && num_loads > 1)
        throw Exception("Cannot do concurrent multiple excitations for multiple loads and multiple frequencies");

      assert(!(num_wave > 0 && (num_freq > 0 || num_loads > 0  || num_trans_ > 0)));

      // sets and resizes excitations with strain loads
      if(DoHomogenization())
        num_loads = SetHomogenizationTestStrains(context_base, ctxt);
    }

    // usually we provide a macrostress in the xml, so we have only one excitation
    // one could change that to more excitations and use unit vectors as test stresses
    if(DoHomogenization() && ctxt->DoBuckling())
      //num_loads += SetHomogenizationTestStrains(context_base, ctxt);
      num_loads += 1;

    // we average the solutions(s) only for output.
    // In the calculations we average the individual calculations
  } // end IsEnabled()
  else
  {
    if(DoRobust(ctxt) && Optimization::manager.context.GetSize() == 1)
      throw Exception("robust filters are defined but 'multiple_excitation' is not enabled in 'costFunction");
    if(num_trans_ > 1)
      throw Exception("transformations are defined but 'multiple_excitation' is not enabled in 'costFunction");
    if(pn->Has("multipleExcitation/excitations"))
      opt->optInfoNode->Get(ParamNode::HEADER)->SetWarning("'multiple_excitations' set to false but 'multipleExcitation/excitations' given");
  }


  // read in tracking parameters from XML, for the first and only load
  if(pn->Has("trackings"))
    excitations[context_base].ReadTrackings(pn->Get("trackings"));

  if(num_wave > 0)
    SetBlochWaves(ctxt, context_base, num_wave);

  // when we do magnetic coupling with 1 frequency but two loads
  if(ctxt->IsHarmonic() && (num_freq > 1 || num_loads <= 1))
    SetHarmonic(ctxt, context_base, num_freq);

  // SetLoadCases() recognizes a harmonic load
  if( IsEnabled(ctxt->sequence) &&
      ( ( (!ctxt->IsComplex() || num_freq == 1) && !ctxt->DoBloch() ) // multiple loads case
      || ( DoHomogenization() && ctxt->DoBuckling() ) ) ) // homogenization with buckling -> SetLoadCases only resizes 'exitations'
    SetLoadCases(ctxt, context_base, pn_ex, num_loads, opt); // when the loads are given in the optimization section of the xml file

  assert(ctxt->sequence >= 1);
  assert((int) excitations.GetSize() >= ctxt->sequence); // at least one excitation per sequence
  for(unsigned int i = context_base; i < excitations.GetSize(); i++)
    excitations[i].sequence = ctxt->sequence;
}

void MultipleExcitation::FinalizeMultipleExcitations(Optimization* opt, ContextManager* manager, bool eval_inital_design)
{
  // for bloch stiffness we might have robust only for stiffness
  for(unsigned int c = 0; c < manager->context.GetSize(); c++)
    if(DoRobust(&(manager->context[c])))
      ApplyRobust(&(manager->context[c])); // multiply by the robust filters

  if(DoTransform())
    for(unsigned int c = 0; c < manager->context.GetSize(); c++)
      ApplyTransformations(&(manager->context[c]), opt->GetDesign()); // multiply the existing transformation

  PtrParamNode in = opt->optInfoNode->Get(ParamNode::HEADER)->Get("excitations");

  for(unsigned int c = 0; c < manager->context.GetSize(); c++)
  {
    Context& ctxt = manager->context[c];
    unsigned int basic_excitations = 0;

    // assign context specific excitations
    for(unsigned int e= 0; e < excitations.GetSize(); e++)
    {
      Excitation& ex = excitations[e];
      if(ex.sequence == ctxt.sequence)
      {
        ctxt.excitations.Push_back(&ex); // assume sufficient small number such that we don't have to deal with reserve
        if(ex.meta_index == 0)
          basic_excitations++; // filter out robust and transformation
      }
    }
    ctxt.SetMultiExcitations(this, basic_excitations);

    // finally verify that the labels are set
    double weight_sum = 0.0;
    for(unsigned int i = 0; i < ctxt.excitations.GetSize(); i++)
    {
      Excitation* ex = ctxt.excitations[i];
      weight_sum += ex->weight;
      if(ex->label == "")
        ex->label = lexical_cast<string>(i);
      if(ctxt.DoMultiSequence()) // relabel to capture sequence
        ex->label = "s_" + lexical_cast<string>(ctxt.sequence) + "-" + ex->label;
    }
    // set excitation index, calculate the initial normalized_weight and print info.
    if(IsEnabled(ctxt.sequence) || ctxt.IsHarmonic())
    {
      // we need to set the global excitation index!
      unsigned int base = 0;
      for(int t = 0; t < ctxt.context_idx; t++) // smaller by intention to not count ourselves
        base += manager->context[t].excitations.GetSize(); // already set up to this point
      for(unsigned int i = 0; i < ctxt.excitations.GetSize(); i++)
      {
        Excitation* ex = ctxt.excitations[i];
        ex->index = base + i;
        // fixes bug for pressure loads
        if (std::abs(weight_sum) < std::numeric_limits<float>::epsilon()) {
          // multiple pressure loads not implemented
          assert(ctxt.excitations.GetSize() == 1);
          ex->normalized_weight = 1.;
        } else {
          ex->normalized_weight = ex->weight / weight_sum;
        }

        LOG_DBG2(exlog) << "FME: i=" << i << " l=" << ex->label << " w=" << ex->weight << " nw=" << ex->normalized_weight << " ws=" << weight_sum << " idx=" << ex->index;
      }
    }

    if(!IsEnabled() && ctxt.num_harm_freq > 1 && !eval_inital_design)
    {
      stringstream ss;
      ss << "Solve only for 1. frequency (" << excitations[0].frequency
          << "Hz) of " << ctxt.num_harm_freq << " as multiple excitations are disabled";
      in->SetWarning(ss.str());
    }
  } // end ctxt loop

  WriteInInfo(in);
}


int MultipleExcitation::SetHomogenizationTestStrains(unsigned int base, Context* ctxt)
{
  unsigned int dim = domain->GetGrid()->GetDim();
  int cases = GetNumberHomogenization(ctxt->ToApp());

  assert(excitations.GetCapacity() >= base + cases * GetNumberMeta(ctxt, true)); // so we need no copy constructor in ApplyTransformation()
  excitations.Resize(base + cases);
  assert((int) MechPDE::X == 0);
  assert((int) MechPDE::Y == 1);
  assert((int) MechPDE::Z == 2);
  assert((int) MechPDE::XY == 5);
  assert((int) MechPDE::X == (int) HeatPDE::X);
  assert((int) MechPDE::Y == (int) HeatPDE::Y);
  assert((int) MechPDE::Z == (int) HeatPDE::Z);

  StdVector<int> testStrains(cases);
  // assumes order of TestStrain Enum to be X=0, Y=1, Z=2, YZ=3, XZ=4, XY=5
  for (int i = 0; i < cases; i++)
    testStrains[i] = i;

  App::Type app = ctxt->ToApp();
  assert(app == App::MECH || app == App::HEAT || app == App::BUCKLING);

  if ((app == App::MECH || app == App::BUCKLING) && dim == 2){
    testStrains.Clear(cases);
    testStrains.Push_back(MechPDE::X);
    testStrains.Push_back(MechPDE::Y);
    testStrains.Push_back(MechPDE::XY);
  }

  for(int i = 0; i < cases; i++)
  {
    MechPDE::TestStrain ts = (MechPDE::TestStrain) testStrains[i];
    Excitation& ex = excitations[base + i];

    //MechPDE has more test strains than HeatPDE, naming is the same
    ex.label = MechPDE::testStrain.ToString(ts);

    if(app == App::MECH || app == App::HEAT)
    {
      ex.ReadTestStrain(ctxt, ts);
      // The homogenized tensor can only be evaluated for the last excitation!
      ex.weight = i < cases-1 ? 0.0 : 1.0;
      LOG_DBG3(exlog) << "SHTS: i=" << i << " f=" << ex.forms.GetSize() << " i=" << (ex.forms.First()->GetIntegrator() == NULL ? "NULL" : ex.forms.First()->GetIntegrator()->GetName());
    }
    LOG_DBG3(exlog) << "SHTS: i=" << i << " base=" << base << " label:" << ex.label << " testStrain: " << ex.test_strain.ToString();
  }

  return cases;
}

unsigned int MultipleExcitation::CountExcitations(const Context* ctxt) const
{
  unsigned int c = 0;
  for(unsigned int i = 0; i < excitations.GetSize(); i++)
    if(excitations[i].sequence == ctxt->sequence)
      c++;
  return c;
}

StdVector<unsigned int> MultipleExcitation::GetExcitations(const Context* ctxt) const
{
  StdVector<unsigned int> res;
  for(unsigned int i = 0; i < excitations.GetSize(); i++)
    if(excitations[i].sequence == ctxt->sequence)
      res.Push_back(i);
  return res;
}


void MultipleExcitation::ApplyRobust(const Context* ctxt)
{
  // we have the case for non multi sequence, then we simply multiply the current excitations
  // for multi sequence we might have bloch where we don't want to do robust and then robust homogenization
  // for the later case we need alt_filter for bloch

  // even continue if we have non robust as we might have to set alt_filter (see description above)
  Robust& robust = robust_[ctxt->context_idx];

  StdVector<unsigned int> org_ex = GetExcitations(ctxt); // the current excitations for this context

  unsigned int principle = excitations.GetSize(); // actually Context::basic_excitations_
  unsigned int n_ex_ctxt = org_ex.GetSize();

  assert(robust.num_robust >= 1 && excitations.GetSize() >= 1); // num_robust_ might be 0 or 1 if we don't do robust

  assert(excitations.GetCapacity() >= (principle - n_ex_ctxt) + n_ex_ctxt * robust.num_robust); // the *additional* excitations only!
  excitations.Resize((principle - n_ex_ctxt) + n_ex_ctxt * robust.num_robust);

  LOG_DBG2(exlog) << "AR: c=" << ctxt->context_idx << " principle=" << principle << " n_ex_ctxt=" << n_ex_ctxt
      << " num_robust=" << robust.num_robust << " alt_filter=" << robust.alt_filter;

  // we now loop over excitations where the later part might is uninitialized in the robust case after resize
  // we start from 0 to get also the existing ones
  for(int r = 0; r < robust.num_robust; r++)
  {
    // traverse over the original excitations
    for(unsigned int b = 0; b < org_ex.GetSize(); b++)
    {
      unsigned org_idx = org_ex[b];
      Excitation& base = excitations[org_ex[b]]; // independent from r the original stuff
      assert(base.sequence == ctxt->sequence);

      // the new index is org_idx if r==0 or it is in the additional resized excitations space >= principle
      unsigned int new_idx = r == 0 ? org_idx : principle + (r-1) * n_ex_ctxt + b;
      Excitation& ex   = excitations[new_idx]; // ex == base for the first transform

      if(r > 0) // from the Resize() above only the first block is set. Copy the test strains, loads or frequencies
        ex = base; // default copy constructors are cool

      if(b == 0) // the original test strains and loads did not need to reassemble. We need it for a new type of transform only
        ex.reassemble = true;

      ex.robust = true;
      ex.meta_index = r; // this might be changed when we do transformation
      // for the multi sequence case when not all is robust we have to set alt_filter!!
      if(robust.num_robust > 1)
        ex.robust_filter_idx = r;
      else
        ex.robust_filter_idx = robust.alt_filter;

      // only set a label if it was not set already
      assert(!(ex.label == "" && principle > 1));
      if(ex.label == "" || principle == 1)
        ex.label = boost::lexical_cast<std::string>(r);

      // TODO up to now we do not handle weights for the transformations!
      // If we have total_base_ > 1 we assume the weights were properly set and are copied
      if(principle == 1)
        ex.weight = 1.0;

      LOG_DBG2(exlog) << "AR: r=" << r << " b=" << b << " org_idx=" << org_idx << " new_idx=" << new_idx << " rfi=" << ex.robust_filter_idx
          << " f=" << ex.forms.GetSize() << " i=" << (ex.forms.First()->GetIntegrator() == NULL ? "NULL" : ex.forms.First()->GetIntegrator()->GetName())
          << " m=" << ex.meta_index << " ra=" << ex.reassemble << " w=" << ex.weight;
    }
  }
}


void MultipleExcitation::ApplyTransformations(const Context* ctxt, DesignSpace* space)
{
  assert(!ctxt->DoMultiSequence()); // just not implemented. robustness needs to apply within sequence?
  //assert(excitations.GetSize() == principle_); // first robust, then transformation
  unsigned int principle = DoRobust(ctxt) ? excitations.GetSize() / GetNumberRobust(ctxt) : excitations.GetSize(); // actually Context::basic_excitations_

  StdVector<Transform>& trans = space->transform;

  assert(num_trans_ == (int) trans.GetSize());

  // robust comes before transformation!
  assert((!DoRobust(ctxt) && excitations.GetSize() == principle) || (DoRobust(ctxt) && excitations.GetSize() == principle * GetNumberRobust(ctxt)));
  assert(num_trans_ >= 1 && excitations.GetSize() >= 1);
  unsigned int old_base = excitations.GetSize();

  // multiply excitations. Robust comes first, the transformation
  assert(excitations.GetCapacity() >= principle * GetNumberMeta(ctxt, true));
  excitations.Resize(principle * GetNumberMeta(ctxt, true));

  for(unsigned int t = 0; t < trans.GetSize(); t++)
  {
    Transform* tr = &trans[t];

    for(unsigned int b = 0; b < old_base; b++)
    {
      Excitation& base = excitations[b];
      Excitation& ex   = excitations[old_base * t + b]; // ex == base for the first transform

      if(t > 0) // from the Resize() above only the first block is set. Copy the test strains, loads or frequencies
        ex = base; // default copy constructors are cool

      if(b == 0) // the original test strains and loads did not need to reassemble. We need it for a new type of transform only
        ex.reassemble = true;

      ex.transform = tr; // the transformations are block wise 0,0,0,0,1,1,1,1,2,2,2,2,....
      // the meta_index handles robust and transformation
      if(!DoRobust(ctxt))
        ex.meta_index = t;
      else
        ex.meta_index = t * GetNumberRobust(ctxt) + ex.robust_filter_idx;

      // only set a label if it was not set already
      assert(!(ex.label == "" && principle > 1));
      if(ex.label == "" || principle == 1)
        ex.label = boost::lexical_cast<std::string>(ex.meta_index);

      // TODO up to now we do not handle weights for the transformations!
      // If we total_base_ > 1 we assume the weights were properly set and are copied
      if(principle == 1)
        ex.weight = 1.0;

      LOG_DBG3(exlog) << "AT: t=" << t << " b=" << b << " f=" << ex.forms.GetSize() << " i=" << (ex.forms.IsEmpty() || ex.forms.First()->GetIntegrator() == NULL ? "NULL" : ex.forms.First()->GetIntegrator()->GetName())
                          << " m=" << ex.meta_index << " ra=" << ex.reassemble << " w=" << ex.weight;
    }
  }
}

LinearFormContext* MultipleExcitation::SearchFormByCoilId(StdVector<LinearFormContext*>& forms, const string& query)
{
  for(unsigned int i = 0; i < forms.GetSize(); i++)
  {
    LinearFormContext* fc = forms[i];
    LinearForm* lf = fc->GetIntegrator();
    LOG_DBG3(exlog) << "ME::SFBCI: Name of form number " << i << " is " << lf->GetName();
    if(lf->GetName() != "CoilIntegrator"){
      LOG_DBG3(exlog) << "ME::SFBCI: Name of form number " << i << " is " << lf->GetName();
      continue;
    }
    // is there a difference 2D/3D??
    string id = lf->IsComplex() ? dynamic_cast<BUIntegrator<Complex>*>(lf)->GetId() : dynamic_cast<BUIntegrator<double>*>(lf)->GetId();
    LOG_DBG3(exlog) << "ME::SFBCI: query: " << query << " id: " << id;
    if(query == id){
      LOG_DBG2(exlog) << "ME::SFBCI: Found form number " << i << " corresponding to coil " << query;
      return fc;
    }
  }
  // nothing found
  throw Exception("could not find CurlIntegration with coil id " + query + " in " + to_string(forms.GetSize()) + " forms");
}

void MultipleExcitation::SetCoils(unsigned int base, Assemble* ass, const ParamNodeList& pn_ex, int num_loads, MathParser* parser, unsigned int handle) {
  // take the loads from the bcsAndLoads section. Store them here and reduce them to one.
  // Apply() will then change the entities for this loads
  // Therefore we do not call ReadLoads() which would read the loads now. This does not work as nicely for Magnetic, therefore we use the approach described above
  StdVector<LinearFormContext*>& forms = ass->GetLinForms(true); // take the memory ownership

  ParamNodeList exc = ParamNode::GetListByChild(pn_ex, "excitation");
  for(int l = 0; l < num_loads; l++)
  {
    // get weight
    PtrParamNode pnex = pn_ex[l];
    parser->SetExpr(handle, pnex->Get("weight")->As<std::string>());
    const double weight = parser->Eval(handle);
    Excitation& ex = excitations[base + l];
    // get coils for this excitation
    // exc is an excitation which is @weight and coils; we want to extract the coils
    ParamNodeList children = exc[l]->GetChildren();
    ParamNodeList coils = ParamNode::GetListByChild(children, "coil");
    LOG_DBG(exlog) << "ME:SLC: Number of coils in excitation number "<< l << " : " << coils.GetSize();
    // push all coils
    for(unsigned int c = 0; c < coils.GetSize(); c++) {
      string query = coils[c]->GetChild()->ToString(2);
      LOG_DBG(exlog) << "ME:SLC: adding form of coil with id " << query;
      LinearFormContext* excForm = SearchFormByCoilId(forms, query); // throws when not found
      ex.forms.Push_back(excForm); // add form to excitation
      LOG_DBG(exlog)<< "ME:SLC: form " << excForm->GetIntegrator()->GetName() << " entities=" << excForm->GetEntities()->GetName();
    }

    // here we support to give a weight in the optimization section
    ex.weight = weight;
    LOG_DBG(exlog)<< "ME:SLC: Setting weight " << weight << " for excitation number " << l;
  }

  // "remove" the loads from the simulation. From now on we Apply() it
  forms.Resize(0); // won't delete content but set the internal size_ counter
}

void MultipleExcitation::SetLoadCases(Context* ctxt, unsigned int base, const ParamNodeList& pn_ex, int num_loads, Optimization* opt)
{
  Assemble* ass = ctxt->pde->GetAssemble();

  assert((ctxt->sequence == 1 && base == 0) || ((int) base >= ctxt->sequence-1));
  assert(excitations.GetCapacity() >= base + num_loads * GetNumberMeta(ctxt, true));
  excitations.Resize(base + num_loads);

  LOG_DBG(exlog)<< "ME:SLC base=" << base << " pn=" << pn_ex.GetSize() << " nl=" << num_loads << " a-lf size=" << ass->GetLinForms().GetSize();

  // when the loads are given in the optimization section of the xml file
  if(pn_ex.GetSize() > 0)
  {
    // this allows flexible weights, currently only when the loads are given in the optimization part and not for the simulation
    MathParser* parser = domain->GetMathParser();
    unsigned int handle = parser->GetNewHandle();

    // the coils are a special case, as we need here the definition in pde/bcsAndLoads AND also the references in optimization/multipleExcitation
    if(ParamNode::GetListByGrandChild(pn_ex , "coil").GetSize() > 0) {
      SetCoils(base, ass, pn_ex, num_loads, parser, handle);
    } else {
      // don't mix optimization excitations with bcsAndLoads stuff
      if(ass->GetLinForms().GetSize() > 0) {
        std::string msg = "Excitations are given in optimization mulitload. " + to_string(ass->GetLinForms().GetSize()) + " loads from the bcsAndLoads section will be deleted.";
        opt->optInfoNode->Get(ParamNode::HEADER)->SetWarning(msg);
      }

      ass->GetLinForms().Clear(); // the forms eventually read by cfs are ignored, we use our own from optimization

      assert(excitations.GetSize() - base == pn_ex.GetSize());

      for(int i = 0; i < num_loads; i++)
      {
        PtrParamNode pnex = pn_ex[i];
        parser->SetExpr(handle, pnex->Get("weight")->As<std::string>());
        // we should not land here in coil case
        assert(!pnex->Has("coil"));
        const double weight = parser->Eval(handle);
        Excitation& ex = excitations[base + i];
        ex.weight = weight;

        if (pnex->Has("trackings"))
          ex.ReadTrackings(pnex->Get("trackings"));

        ex.ReadLoads(ctxt, pnex); // possibly multiple forms for one excitation has its own Resize(0)
      }
    }

    parser->ReleaseHandle(handle);
  }
  // take the loads from the bcsAndLoads section. Store them here and reduce them to one.
  // Apply() will then change the entities for this loads
  if(pn_ex.GetSize() == 0 && !DoHomogenization() && !ctxt->DoBloch())
  {
    // a force is a linarForm with integrator NodalForceInt
    StdVector<LinearFormContext*>& forms = ass->GetLinForms(true); // take the memory ownership
    assert((int) forms.GetSize() == num_loads);

    for(int i = 0; i < num_loads; i++) // via num_loads or num_freq
    {
      // static load case
      Excitation& ex = excitations[base + i];
      ex.forms.Push_back(forms[i]); // one form for one excitation
      LOG_DBG(exlog)<< "PME: form " << i << " = " << forms[i]->GetIntegrator()->GetName() << " entities=" << forms[i]->GetEntities()->GetName();

      // we do not support to give a weight in the force section of the simulation
      ex.weight = 1.0;
    }

    // "remove" the loads from the simulation. From now on we Apply() ist
    forms.Resize(0); // won't delete content but set the internal size_ counter
  }

  // check for the harmonic single frequency case, then we also store this information in the excitation
  assert(!(ctxt->IsHarmonic() && ctxt->num_harm_freq != 1)); // might also be extended for multiple frequencies
  if(ctxt->IsHarmonic())
  {
    for(int i = 0; i < num_loads; i++) // via num_loads or num_freq
    {
      Excitation& ex = excitations[base + i];
      SetHarmonicExcitation(ctxt, ex, 0); // we have exactly one frequency!
      ex.label = ""; // we don't want the frequency as label
    }
  }
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

bool MultipleExcitation::DoRobust(const Context* ctxt) const
{
  assert(ctxt != NULL);
  return robust_[ctxt->context_idx].num_robust > 1;
}

unsigned int MultipleExcitation::GetNumberRobust(const Context* ctxt, bool mininum_one) const
{
  if(ctxt != NULL)
    return mininum_one ? std::max(robust_[ctxt->context_idx].num_robust, 1) : robust_[ctxt->context_idx].num_robust;

  for(unsigned int i = 0; i < robust_.GetSize(); i++)
    if(robust_[i].num_robust > 1)
      return robust_[i].num_robust;

  // nothing found
  return mininum_one ? 1 : 0;
}

unsigned int MultipleExcitation::GetNumberMeta(const Context* ctxt, bool minimum_one) const
{
  return GetNumberTransform(minimum_one) * GetNumberRobust(ctxt, minimum_one); // GNR accepts NULL for ctxt
}

unsigned int MultipleExcitation::GetUniqueFrequencies() const
{
  StdVector<double> freq;
  for(unsigned int i = 0; i < excitations.GetSize(); i++)
  {
    const Excitation& ex = excitations[i];
    if(ex.frequency >= 0 && !freq.Contains(ex.frequency))
      freq.Push_back(ex.frequency);
  }
  return freq.GetSize();
}


bool MultipleExcitation::DoMetaExcitation(const Context* ctxt) const
{
  if(DoTransform())
    return true;

  if(ctxt != NULL)
    return DoRobust(ctxt);

  for(unsigned int i = 0; i < robust_.GetSize(); i++)
    if(robust_[i].num_robust > 1)
      return true;

  return false;
}

bool MultipleExcitation::DoMetaExcitation(int sequence) const
{
  assert(sequence >= 1);
  const Context& ctxt = Optimization::manager.context[sequence-1];
  return DoMetaExcitation(&ctxt);
}

Excitation::Excitation()
{
  this->index = -1; // must be updated
  this->sequence = -1; // set correctly in prepare
  this->meta_index = 0;
  this->frequency = -1.0;
  this->f_link = NULL;
  this->cost = -1.0;
  this->transform = NULL;
  this->robust = false;
  this->robust_filter_idx = 0; // a good default value
  this->weight = 1.0;
  this->normalized_weight = 1.0;
  this->reassemble = false; // normally
}

Excitation::~Excitation()
{
  // we have the ownerhip of the loads (form) with multiple excitation. We took if from Assemble::linForms_
  // Coils are a special case, they are owned by the PDE
  for(unsigned int i = 0; i < forms.GetSize(); i++)
    if(meta_index == 0)
      if(forms[i]->GetIntegrator()->GetName() != "CoilIntegrator")
         delete forms[i];
}

bool Excitation::Apply(bool switch_context)
{
  LOG_DBG(exlog) << "A: sc=" << switch_context << " curr_ex=" << Optimization::context->GetExcitation()->index << " new_ex=" << index;
  bool switched = false;
  if(switch_context && Optimization::context->sequence != this->sequence)
  {
    Optimization::manager.SwitchContext(this); // also sets the excitation
    switched = true;
    assert(Optimization::context->sequence == this->sequence);
    LOG_DBG(exlog) << "A: switched context to sequence " << sequence;
  }
  Optimization::context->SetExcitation(this);

  assert(!switch_context || Optimization::context->sequence == sequence);

  if(forms.GetSize() > 0)
  {
    Context& ctxt = Optimization::manager.GetContext(this);
    Assemble* ass = ctxt.pde->GetAssemble();
    ass->GetLinForms() = forms; // let the copy constructor do the stuff

    assert(ass->GetLinForms().GetSize() == forms.GetSize());
    LOG_DBG(exlog) << "A: set form " << forms[0]->ToString();
    assert(ctxt.pde->GetAnalysisType() == forms[0]->GetPde()->GetAnalysisType());
  }

  // a frequency cannot really be applied but has to be used as parameter
  // in the driver call
  // the same holds for the bloch mode analysis
  return switched;
}

std::map<RegionIdType, PtrCoefFct> Excitation::GetStressCoefFctFromExcitation(UInt index)
{
  std::map<RegionIdType, PtrCoefFct> stress_map;

  Context& ctxt = Optimization::manager.GetContext(this);
  assert(ctxt.GetBucklingDriver());
  assert(domain->GetSingleDriver()->GetAnalysisType() == BasePDE::BUCKLING);

  RegionIdType actRegion;
  SinglePDE* pde = ctxt.pde; // buckling PDE (complex!)

  //  Loop over all regions
  std::map<RegionIdType, BaseMaterial*>& materials = pde->GetMaterialData();
  for(auto mat : materials)
  {
    // Set current region
    actRegion = mat.first;

    // get excitation at index
    Optimization* opt = domain->GetOptimization();
    assert(opt->GetMultipleExcitation()->excitations.GetSize() > 1);
    Excitation& excite = opt->GetMultipleExcitation()->excitations[index];

    // get the solution (displacement) of this excitation
    ErsatzMaterial* em = dynamic_cast<ErsatzMaterial*>(opt);
    StateSolution* ss = em->forward.Get(excite, NULL, -1);
    Vector<Double> mech_displ;
    mech_displ = ss->GetRealVector(StateSolution::RAW_VECTOR);

    // convert to complex, as buckling PDE is complex
    Vector<Complex> complex_mech_displ;
    complex_mech_displ.Resize(mech_displ.GetSize());
    complex_mech_displ.Init();
    complex_mech_displ.SetPart(Global::REAL, mech_displ);

    // inject solution into buckling pde
    ss->Write(pde, &complex_mech_displ);

    // get coef function for stresses
    PtrCoefFct stress = pde->GetCoefFct(MECH_STRESS);
    assert(stress && stress.get());

    stress_map[actRegion] = stress;
  }
  return stress_map;
}

// in buckling replace the old stresses, i.e. previous iteration
// (or dummy for first iteration), by the new ones
void Excitation::SetStressCoefFct(std::map<RegionIdType, PtrCoefFct> stress_map)
{
  Context& ctxt = Optimization::manager.GetContext(this);
  assert(ctxt.GetBucklingDriver());
  assert(domain->GetSingleDriver()->GetAnalysisType() == BasePDE::BUCKLING);

  RegionIdType actRegion;
  SinglePDE* pde = ctxt.pde; // buckling PDE (complex!)

  //  Loop over all regions
  std::map<RegionIdType, BaseMaterial*>& materials = pde->GetMaterialData();
  assert(materials.size() == stress_map.size());
  for(auto mat : materials)
  {
    // Set current region
    actRegion = mat.first;
    // we want to change the stress integrator of the geometric stiffness matrix
    BiLinFormContext* blfc = pde->GetAssemble()->GetBiLinForm("PreStressInt", actRegion, pde, pde, false);

    // the coef function of the integrator is a compound, i.e. a bundle of multiple coeffunctions
    BDBInt<Double,Double>* bdb = dynamic_cast<BDBInt<Double,Double>*>(blfc->GetIntegrator());
    CoefFunctionCompound<Double>* stressTens = dynamic_cast<CoefFunctionCompound<Double>*>(ctxt.mat->GetOrgMatCoef(bdb));
    assert(stressTens != NULL);

    // get map of coeffunctions
    std::map<std::string, PtrCoefFct>& map = stressTens->GetCoefFcts();

    // set the coeffunction for stresses. key "a" was chosen in MechPDE::CreatePreStressFct
    assert(map.find("a") != map.end());
    map["a"] = stress_map[actRegion];
  }

  this->reassemble = true;
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
  if(f->DoEvaluateAlways(sequence)) // sequence matches always
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

void Excitation::ReadLoads(Context* ctxt, PtrParamNode ls)
{
  // the original loads need to be cleared before
  Assemble* ass = ctxt->pde->GetAssemble();
  assert(ass->GetLinForms().GetSize() == 0);
  LOG_DBG(exlog) << "RL: paramnode ls=" << ls->GetChildren().GetSize() << " ls=" << ls->GetName();
  // reads all loads and adds them to Assemble::linForms_
  ctxt->pde->DefineRhsLoadIntegrators(ls);

  LOG_DBG(exlog) << "RL: paramnode ls: " << ls << " ctxt exc size: " << ctxt->excitations.GetSize();
  // own vector with the pointers in assemble and we have to delete the content
  forms = ass->GetLinForms(true); // take ownership. copy constructor!
  assert(forms.GetSize() > 0);

  // delete the stuff from assemble to prepare for the next (()
  ass->GetLinForms().Resize(0); // the elements must not be destructed
}

void Excitation::ReadTestStrain(Context* ctxt, MechPDE::TestStrain ts)
{
  // the original loads need to be cleared before
  Assemble* ass = ctxt->pde->GetAssemble();
  forms = ass->GetLinForms(true); // take ownership for our destructor which is common for ReadLoads(). copy constructor!
  assert(ass->GetLinForms().GetSize() == 0);
  assert(ctxt->ToApp() == App::MECH || ctxt->ToApp() == App::HEAT);
  switch(ctxt->ToApp())
  {
    case App::MECH:
    {
      MechPDE* mech = dynamic_cast<MechPDE*>(ctxt->pde);
      assert(mech->GetAnalysisType() == BasePDE::STATIC);
      mech->DefineTestStrainIntegrator(ts, &forms);
      break;
    }
    case App::HEAT:
    {
      HeatPDE* heat = dynamic_cast<HeatPDE*>(ctxt->pde);
      assert(heat->GetAnalysisType() == BasePDE::STATIC);
      heat->DefineTestStrainIntegrator((HeatPDE::TestStrain) ts, &forms);
      break;
    }
    default:
      throw Exception("ReadTestStrain() only implemented for MECH and HEAT!");
      break;
  }

  test_strain.Resize(6, 0.0); // always full dimensions
  test_strain[ts] = 1.0;

  LOG_DBG3(exlog) << test_strain.ToString();
}

std::string Excitation::GetMetaLabel() const
{
  if(transform == NULL && !robust)
    return "";

  std::stringstream ss;
  if(robust)
    ss << "f_" << robust_filter_idx;
  if(transform != NULL)
    ss << (robust ? "_" : "") << "t_" << transform->ToString(0);

  return ss.str();
}

std::string Excitation::GetFullLabel() const
{
  if(transform == NULL && !robust)
    return label;
  else
    return label + "_" + GetMetaLabel();
}

bool Excitation::DoBloch() const
{
  return Optimization::manager.GetContext(this).DoBloch();
}


int Excitation::GetWaveNumber() const
{
  assert(Optimization::manager.GetContext(this).DoBloch());
  // the wave number is for the first sequence the excitation index, for higher sequence we need to make it relative
  return index - Optimization::manager.GetContext(this).excitations[0]->index;
}


