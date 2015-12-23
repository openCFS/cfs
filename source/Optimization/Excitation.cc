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
#include "Driver/EigenFrequencyDriver.hh"
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
  this->type_ = NO_TYPE; // to be possibly overwritten soon
  this->num_trans_  = -1; // to be set later
  this->num_robust_ = -1; // to be set later
  this->sequence_ = 1; // default

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
  for(unsigned int i = 0; i < excitations.GetSize(); i++)
    if(excitations[i].label == label)
      return &(excitations[i]);

  if(quiet)
    return NULL;
  else
    throw Exception("None of the " + lexical_cast<string>(excitations.GetSize()) + " excitations has a label '" + label + "'");
}



unsigned int MultipleExcitation::GetNumberHomogenization() const
{
  // when we have only one sequence total_base_ must be the number of strains
  if(!DoHomogenization())
    return 0;
  else
    return domain->GetGrid()->GetDim() == 2 ? 3 : 6;
}

bool MultipleExcitation::IsEnabled(int sequence) const
{
  assert(sequence == -1 || (sequence >= 1 && sequence <= (int) Optimization::manager.context.GetSize()));

  if(sequence == -1)
    return multiple_excitation_;
  else
    return multiple_excitation_ && sequence == sequence_;
}

void MultipleExcitation::WriteInInfo(int num_freq, bool eval_inital_design,  double weight_sum, Optimization* opt)
{
  PtrParamNode in = opt->optInfoNode->Get(ParamNode::HEADER)->Get("excitations");

  if(!IsEnabled() && num_freq > 1 && !eval_inital_design)
  {
    stringstream ss;
    ss << "Solve only for 1. frequency (" << excitations[0].frequency
        << "Hz) of " << num_freq << " as multiple excitations are disabled";
    in->Get(ParamNode::WARNING)->SetValue(ss.str());
  }

  // communicate what we have but also normalize the weights!
  for(unsigned int i = 0; i < excitations.GetSize(); i++)
  {
    Excitation& ex = excitations[i];

    PtrParamNode exin = in->Get("excitation", ParamNode::APPEND);
    exin->Get("index")->SetValue(ex.index);
    exin->Get("label")->SetValue(ex.label);
    exin->Get("meta_idx")->SetValue(ex.meta_index);
    if(Optimization::context->DoMultiSequence())
      exin->Get("sequence")->SetValue(ex.sequence);
    if(ex.transform != NULL)
      exin->Get("transform")->SetValue(ex.transform->ToString(1));
    if(ex.robust)
      exin->Get("robust_filter_idx")->SetValue(ex.robust_filter_idx);
    if(ex.transform || ex.robust)
      exin->Get("full_label")->SetValue(ex.GetFullLabel());
    exin->Get("weight")->SetValue(ex.weight);
    exin->Get("normalized_weight")->SetValue(ex.normalized_weight);
    if(ex.forms.GetSize() == 1)
      exin->Get("load")->SetValue(ex.forms[0]->GetEntities()->GetName());
    if(DoMetaExcitation())
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

  HarmonicDriver* hd = ctxt->GetHarmonicDriver();

  assert(excitations.Capacity() >= base + num_freq);
  excitations.Resize(base + num_freq);

  for (unsigned int i = 0; i < excitations.GetSize(); i++)
  {
    Excitation& ex = excitations[base + i];
    ex.frequency = hd->freqs[i].freq;
    ex.label = lexical_cast<string>(ex.frequency);
    ex.weight = hd->freqs[i].weight;
    ex.f_link = &hd->freqs[i];
    ex.reassemble = true;
  }
}

void MultipleExcitation::SetBlochWaves(Context* ctxt, unsigned int base, int num_wave)
{
  const StdVector<Vector<double> >& wv =  ctxt->GetEigenFrequencyDriver()->wave_vectors;

  assert(excitations.Capacity() >= base + num_wave);
  excitations.Resize(base + num_wave);

  for (unsigned int i = 0; i < excitations.GetSize(); i++)
  {
    Excitation& ex = excitations[base + i];
    ex.weight = i < excitations.GetSize() - 1 ? 0 : 1; // we assume the slack variable as objective!
    ex.wave_vector = wv[i];
    ex.label = "(" + ex.wave_vector.ToString(0, ',') + ")";
  }
}


void MultipleExcitation::InitializeMultipleExcitations(Optimization* opt, ContextManager* manager)
{
  // this is to be called once prior to any PrepareMultipleExcitations() call

  // we have no "do_transform" either there is transformation or not
  num_trans_  = opt->GetDesign()->transform.GetSize();
  num_robust_ = opt->GetDesign()->data[0].simp->filter.GetSize();

  // We do not want to have a copy of Excitations on excitation.Resize() therefore find an upper bound
  unsigned int upper_bound = 1;

  for(unsigned int i = 0; i < manager->context.GetSize(); i++)
  {
    Context& ctxt = manager->context[i];

    // no C++11 yet :(
    unsigned int tmp = std::max(ctxt.num_harm_freq, ctxt.num_bloch_wave_vectors);
    // we don't know determine the real number of loads, hence make an estimate. It's just to reserve!
    if(ctxt.analysis == BasePDE::STATIC)
      tmp = std::max(tmp, (unsigned int) 10); // FIXME we don't have Context::num_static_loads yet
    upper_bound += tmp;
  }

  upper_bound *= std::max(1,num_trans_) * std::max(1,num_robust_);
  excitations.Reserve(upper_bound);
}

void MultipleExcitation::PrepareMultipleExcitations(Optimization* opt, Context* ctxt)
{
  // we assume InitializeMultipleExcitations() to be called
  assert(num_trans_ >= 0 && num_robust_ >= 0);

  Assemble* ass = ctxt->pde->GetAssemble();

  // the already existing excitations from prior a context
  unsigned int context_base = excitations.GetSize();
  assert((ctxt->context_idx == 0 &&  context_base == 0) || (ctxt->context_idx > 0 && context_base > 0)); // at least one per context

  // by definition, we always have at least one excitation.
  // This does not necessarily need to be a load (but might be voltage, pressure, ...)
  assert(excitations.Capacity() >= context_base + 1); // avoid moving by resize like hell
  excitations.Resize(context_base + 1);
  excitations[context_base].index = context_base;

  PtrParamNode pn = opt->optParamNode->Get("costFunction");

  // the actual multipleExcitation description is read in Optimization as part of
  // objective function block
  int num_freq = ctxt->num_harm_freq;

  // bloch mode analysis wave vectors
  int num_wave = ctxt->num_bloch_wave_vectors;

  int num_loads = ass->GetLinForms().GetSize();  // to be faked later for homogenization test strains

  LOG_DBG(exlog) << "PME: linForms from assemble: " << num_loads << " trans: " << num_trans_;

  // this might become explicitly given excite in the optimization xml description
  ParamNodeList pn_ex;

  LOG_DBG(exlog) << "PME: cs=" << ctxt->sequence << " en=" << IsEnabled() << " seq=" << sequence_ << " cb=" << context_base << " cap=" << excitations.Capacity();

  // initialize data and do simple plausibility check. Note that also 1 is "multiple"
  // only for our sequence. We assume only one multipleExcitations element in xml even for multiple sequence optimzation
  if(IsEnabled(ctxt->sequence) && sequence_ == ctxt->sequence)
  {
    // either every single load from bcsAndLoads is an excitation or allow combinations of loads, pressures, regionLoads
    // and trackings in one excitation when specified in multipleExcitation (only non-harmonic) (this is done here)
    if(!ctxt->IsHarmonic() && pn->Has("multipleExcitation/excitations"))
    {
      pn_ex = pn->Get("multipleExcitation/excitations")->GetChildren();
      num_loads = pn_ex.GetSize();
    }

    // decide if we have multiple loads or multiple frequencies
    // we cannot do both!

    if(num_freq > 1 && num_loads > 1)
      throw Exception("Cannot to concurrent multiple excitations for multiple loads and multiple frequencies");

    assert(!(num_wave > 0 && (num_freq > 0 || num_loads > 0  || num_trans_ > 0)));

    // sets and resizes excitations with strain loads
    if(DoHomogenization())
      num_loads = SetHomogenizationTestStrains(context_base, ctxt);

    // we average the solutions(s) only for output.
    // In the calculations we average the individual calculations
  } // end IsEnabled()
  else
  {
    if(DoRobust())
      throw Exception("robust filters are defined but 'multiple_excitation' is not enabled in 'costFunction");
    if(num_trans_ > 1)
      throw Exception("transformations are defined but 'multiple_excitation' is not enabled in 'costFunction");
    if(pn->Has("multipleExcitation/excitations"))
      opt->optInfoNode->Get(ParamNode::HEADER)->Get(ParamNode::WARNING)->SetValue("'multiple_excitations' set to false but 'multipleExcitation/excitations' given");
  }


  // read in tracking parameters from XML, for the first and only load
  if(pn->Has("trackings"))
    excitations[context_base].ReadTrackings(pn->Get("trackings"));

  if(num_wave > 0)
    SetBlochWaves(ctxt, context_base, num_wave);

  if(ctxt->IsHarmonic())
    SetHarmonic(ctxt, context_base, num_freq);

  if(!ctxt->IsComplex() && IsEnabled(ctxt->sequence) && !ctxt->DoBloch()) // multiple loads case
    SetLoadCases(ctxt, context_base, pn_ex, num_loads, opt); // when the loads are given in the optimization section of the xml file

  assert(ctxt->sequence >= 1);
  assert((int) excitations.GetSize() >= ctxt->sequence); // at least one excitation per sequence
  for(unsigned int i = context_base; i < excitations.GetSize(); i++)
    excitations[i].sequence = ctxt->sequence;
}

void MultipleExcitation::FinalizeMultipleExcitations(Optimization* opt, ContextManager* manager, bool eval_inital_design)
{
  if(DoRobust())
    ApplyRobust(opt->GetDesign()); // multiply by the robust filters

  if(DoTransform())
    ApplyTransformations(opt->GetDesign()); // multiply the existing transformation

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
      for(unsigned int t = 0; t < ctxt.context_idx; t++) // smaller by intention to not count ourselves
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

        LOG_DBG3(exlog) << "PME: i=" << i << " l=" << ex->label << " w=" << ex->weight << " nw=" << ex->normalized_weight << " ws=" << weight_sum;
      }
      WriteInInfo(ctxt.num_harm_freq, eval_inital_design, weight_sum, opt);
    }
  }
}



int MultipleExcitation::SetHomogenizationTestStrains(unsigned int base, Context* ctxt)
{
  unsigned int dim = domain->GetGrid()->GetDim();

  int cases = dim == 2 ? 3 : 6;
  assert(excitations.Capacity() >= base + cases * GetNumberMeta(true)); // so we need no copy constructor in ApplyTransformation()
  excitations.Resize(base + cases);

  assert((int) MechPDE::X == 0);
  assert((int) MechPDE::XY == 5);

  for(int i = 0, cnt = 0; i < 6; i++)
  {
    MechPDE::TestStrain ts = (MechPDE::TestStrain) i;
    Excitation& ex = excitations[base + cnt];

    // in 2D only 0, 1 and 5
    if(dim == 2 && (i == 2 || i == 3 || i == 4)) continue;

    ex.label = MechPDE::testStrain.ToString(ts);

    ex.ReadTestStrain(ctxt, ts);
    // The homogenized tensor can only be evaluated for the last excitation!
    ex.weight = i < cases-1 ? 0.0 : 1.0;
    ++cnt;

    LOG_DBG3(exlog) << "SHTS: i=" << i << " f=" << ex.forms.GetSize() << " i=" << (ex.forms.First()->GetIntegrator() == NULL ? "NULL" : ex.forms.First()->GetIntegrator()->GetName());
  }

  return cases;
}

void MultipleExcitation::ApplyRobust(DesignSpace* space)
{
  assert(!Optimization::context->DoMultiSequence()); // just not implemented. robustness needs to apply within sequence?
  //assert(excitations.GetSize() == principle_); // first robust, then transformation
  unsigned int principle = excitations.GetSize(); // actually Context::basic_excitations_
  assert(num_robust_ >= 1 && excitations.GetSize() >= 1); // num_robust_ might be 0 or 1 if we don't do robust

  // multiply excitations
  assert(excitations.Capacity() >= principle * num_robust_);
  excitations.Resize(principle * num_robust_);

  for(int r = 0; r < num_robust_; r++)
  {
    for(unsigned int b = 0; b < principle; b++)
    {
      Excitation& base = excitations[b];
      Excitation& ex   = excitations[principle * r + b]; // ex == base for the first transform

      if(r > 0) // from the Resize() above only the first block is set. Copy the test strains, loads or frequencies
        ex = base; // default copy constructors are cool

      if(b == 0) // the original test strains and loads did not need to reassemble. We need it for a new type of transform only
        ex.reassemble = true;

      ex.robust = true;
      ex.meta_index = r; // this might be changed when we do transformation
      ex.robust_filter_idx = r;

      // only set a label if it was not set already
      assert(!(ex.label == "" && principle > 1));
      if(ex.label == "" || principle == 1)
        ex.label = boost::lexical_cast<std::string>(r);

      // TODO up to now we do not handle weights for the transformations!
      // If we have total_base_ > 1 we assume the weights were properly set and are copied
      if(principle == 1)
        ex.weight = 1.0;

      LOG_DBG3(exlog) << "AR: r=" << r << " b=" << b << " f=" << ex.forms.GetSize() << " i=" << (ex.forms.First()->GetIntegrator() == NULL ? "NULL" : ex.forms.First()->GetIntegrator()->GetName())
                      << " m=" << ex.meta_index << " ra=" << ex.reassemble << " w=" << ex.weight;

    }
  }
}


void MultipleExcitation::ApplyTransformations(DesignSpace* space)
{
  assert(!Optimization::context->DoMultiSequence()); // just not implemented. robustness needs to apply within sequence?
  //assert(excitations.GetSize() == principle_); // first robust, then transformation
  unsigned int principle = DoRobust() ? excitations.GetSize() / num_robust_ : excitations.GetSize(); // actually Context::basic_excitations_

  StdVector<Transform>& trans = space->transform;

  assert(num_trans_ == (int) trans.GetSize());

  // robust comes before transformation!
  assert((!DoRobust() && excitations.GetSize() == principle) || (DoRobust() && excitations.GetSize() == principle * num_robust_));
  assert(num_trans_ >= 1 && excitations.GetSize() >= 1);
  unsigned int old_base = excitations.GetSize();

  // multiply excitations. Robust comes first, the transformation
  assert(excitations.Capacity() >= principle * GetNumberMeta(true));
  excitations.Resize(principle * GetNumberMeta(true));

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
      if(!DoRobust())
        ex.meta_index = t;
      else
        ex.meta_index = t * GetNumberRobust() + ex.robust_filter_idx;

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

void MultipleExcitation::SetLoadCases(Context* ctxt, unsigned int base, const ParamNodeList& pn_ex, int num_loads, Optimization* opt)
{
  Assemble* ass = ctxt->pde->GetAssemble();

  assert((ctxt->sequence == 1 && base == 0) || ((int) base >= ctxt->sequence-1));
  assert(excitations.Capacity() >= base + num_loads * GetNumberMeta(true));
  excitations.Resize(base + num_loads);

  // when the loads are given in the optimization section of the xml file
  if(pn_ex.GetSize() > 0)
  {
    // this allows flexible weights, currently only when the loads are given in the optimization part and not for the simulation
    MathParser* parser = domain->GetMathParser();
    unsigned int handle = parser->GetNewHandle();

    // don't mix optimization excitations with bcsAndLoads stuff
    if (ass->GetLinForms().GetSize() > 0) {
      std::string msg = "Excitations are given in optimization mulitload. "
          + boost::lexical_cast<std::string>(ass->GetLinForms().GetSize())
          + " loads from the bcsAndLoads section will be deleted.";
      opt->optInfoNode->Get(ParamNode::HEADER)->Get(ParamNode::WARNING)->SetValue(msg);
    }

    ass->GetLinForms().Clear(); // the forms are gone!

    assert(excitations.GetSize() - base == pn_ex.GetSize());

    for(int i = 0; i < num_loads; i++)
    {
      parser->SetExpr(handle, pn_ex[i]->Get("weight")->As<std::string>());
      const double weight = parser->Eval(handle);
      Excitation& ex = excitations[base + i];
      ex.weight = weight;

      if (pn_ex[i]->Has("trackings"))
        ex.ReadTrackings(pn_ex[i]->Get("trackings"));

      ex.ReadLoads(ctxt, pn_ex[i]); // possibly multiple forms for one excitation has it's own Resize(0)
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
      // don't do anything for multiple excitation enabled with only one load
      if (num_loads > (int) i)
      {
        // static load case
        Excitation& ex = excitations[base + i];
        ex.forms.Push_back(forms[i]); // one form for one excitation
        LOG_DBG(exlog)<< "PME: form " << i << " = " << forms[i]->GetIntegrator()->GetName() << " entities=" << forms[i]->GetEntities()->GetName();

        // we do not support to give a weight in the force section of the simulation
        ex.weight = 1.0;
      }
    }

    // "remove" the loads from the simulation. From now on we Apply() ist
    forms.Resize(0); // won't delete content but set the internal size_ counter
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
  for(unsigned int i = 0; i < forms.GetSize(); i++)
    if(meta_index == 0)
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

  // reads all loads and adds them to Assemble::linForms_
  ctxt->pde->DefineRhsLoadIntegrators(ls);


  // own vector with the pointers in assemble and we have to delete the content
  forms = ass->GetLinForms(true); // take ownership. copy constructor!
  assert(forms.GetSize() > 0);

  // delete the stuff from assemble to prepare for the next ReadLoads()
  ass->GetLinForms().Resize(0); // the elements must not be destructed
}

void Excitation::ReadTestStrain(Context* ctxt, MechPDE::TestStrain ts)
{
  // the original loads need to be cleared before
  Assemble* ass = ctxt->pde->GetAssemble();
  forms = ass->GetLinForms(true); // take ownership for our destructor which is common for ReadLoads(). copy constructor!
  assert(ass->GetLinForms().GetSize() == 0);

  MechPDE* mech = dynamic_cast<MechPDE*>(ctxt->ToPDE(App::MECH));
  assert(mech->GetAnalysisType() == BasePDE::STATIC);

  mech->DefineTestStrainIntegrator(ts, &forms);

  test_strain.Resize(6, 0.0); // always full dimensions
  test_strain[ts] = 1.0;
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



