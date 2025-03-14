#include "Optimization/Tune.hh"
#include "Optimization/Optimization.hh"
#include "Optimization/Condition.hh"
#include "Optimization/Design/Filter.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

namespace CoupledField {

DEFINE_LOG(tune, "tune")

Enum<Tune::Method> Tune::method;
Enum<Tune::Usage>  Tune::usage;

void Tune::Init(PtrParamNode pn, Usage use)
{
  if(pn == nullptr)
    throw Exception("A 'tune' element is expected but not given.");

  if(method.map.size() == 0)
  {
    method.SetName("Tune::Method");
    method.Add(OBJ, "obj");
    method.Add(MULT, "mult");
    method.Add(ADD, "add");

    usage.SetName("Tune::Usage");
    usage.Add(BETA, "beta");
    usage.Add(PENALTY, "penalty");
    usage.Add(FUNC_SCALE, "fscale");  
  }

  this->usage_ = use;

  // <tune param="beta" method="obj/mult/add" start="1" end="256" grow="1e-4" obj_max_grow="0.2" stride="1" stopping_greyness="true" />
  method_ = method.Parse(pn->Get("method")->As<string>());
  start = pn->Get("start")->As<double>();
  end  = pn->Has("end")  ? pn->Get("end")->As<double>()  : (use == BETA ? 256 : 6);
  minimal  = pn->Has("minimal")  ? pn->Get("minimal")->As<double>() : OFF;
  grow = pn->Has("grow") ? pn->Get("grow")->As<double>() : (method_ == OBJ ? 1e-4 : (method_ == MULT ? 2 : 0.1));
  max_grow_rate = pn->Get("obj_max_grow")->As<double>();
  stride =  pn->Get("stride")->As<unsigned int>();
  stopping_greyness_ = pn->Get("stopping_greyness")->As<bool>();

  // if we start from +/- 0.* we cannot grow (in magnitude) by multiplications
  one_scale = method_ != ADD && start > -1 && start < 1.0;

  LOG_DBG(tune) << "I " << usage.ToString(usage_) << " opt=" << opt << " method=" << method.ToString(method_) << " value=" << value << " start=" << start << " one_scale=" << one_scale;
}

void Tune::Register(double* value, Optimization* opt, GlobalFilter* f)
{
  assert(!(f != nullptr && usage_ != BETA));
  assert(value != nullptr && opt != nullptr);
  assert(!IsRegistered());
  this->value = value;
  this->opt = opt;
  if(f != nullptr)
    this->gf.Push_back(f);

  // tries to find a greyness stopping rule, if not found, a warning is issued
  if(stopping_greyness_)
    FindGraynessStoppingRule();

  // register ourself to get Update() called
  opt->tunes.Push_back(this);

  opt->RegisterAuxLogValue(usage.ToString(usage_), *value);

  // set initial value
  *value = start;

  LOG_DBG(tune) << "R: " << usage.ToString(usage_) << " greyness=" << stopping_greyness_ << " opt=" << opt << " value=" << value;
}

bool Tune::IsRegistered() const
{
  if(opt == nullptr || value == nullptr)
    return false;

  // check if really we are registered or an object we are copied from
  if(std::find(opt->tunes.begin(), opt->tunes.end(), this) != opt->tunes.end())
    return true;
  LOG_DBG(tune) << "IR tune not registered #tunes=" << opt->tunes.GetSize() << " opt " << opt << " value=" << value;
  return false;
}

void Tune::Append(double* value, GlobalFilter* f)
{
  LOG_DBG(tune) << "A " << value << " f=" << f << " e=" << external.ToString(TS_PYTHON) << " gf=" << gf.ToString(TS_PYTHON);

  assert(!external.Contains(value));
  external.Push_back(value);
  *value = *(this->value);

  assert(!gf.Contains(f));
  gf.Push_back(f);
}

void Tune::Remove(double* value, GlobalFilter* f)
{
  int vi = external.Find(value, true); // quiet - somehow debugging shows unknown beta adress in the robust case?!
  int fi = gf.Find(f, true);
  LOG_DBG(tune) << "R " << value << " f=" << f << " vi=" << vi << " fi=" << fi << " e=" << external.ToString(TS_PYTHON) << " gf=" << gf.ToString(TS_PYTHON);
  if(vi >= 0)
    external.Erase((size_t) vi);
  if(fi >= 0)
    gf.Erase((size_t) fi);
}

double Tune::CalcTransistionZone(double y) const
{
  assert(usage_ == BETA);
  double beta = *value;
  double eta = 0.5;

  // wolframalpha: solve ( y = (tanh(beta*eta)+tanh(beta*(x-eta)))/(tanh(beta*eta)+tanh(beta*(1-eta))) for x
  double tz= 1/beta * (std::atanh((y-1)*std::tanh(beta*eta) + y * std::tanh(beta-beta*eta)) + beta*eta);
  return tz;
}

void Tune::ToInfo(PtrParamNode info) const
{
  PtrParamNode in = info->Get("tune"); // extend when we actually have multiple tunes assigned (beta + eta, ...)
  in->Get("variable")->SetValue(usage.ToString(usage_));
  in->Get("method")->SetValue(method.ToString(method_));
  in->Get("start")->SetValue(start);
  in->Get("end")->SetValue(end);
  in->Get("minimal")->SetValue(minimal == OFF ? "-" : std::to_string(minimal));
  in->Get("grow")->SetValue(grow);
  if(method_ == OBJ)
    in->Get("max_grow_rate")->SetValue(max_grow_rate);
  in->Get("one_scale")->SetValue(one_scale);
  in->Get("stride")->SetValue(stride);
  in->Get("stopping")->SetValue(grayness ? grayness->function : "-");

  if(method_ == MULT && grow <= 1)
    in->SetWarning("when using tune method 'mult' the attribute 'grow' shall be greater 1");

}

void Tune::Update(unsigned int iter)
{
  assert(value != nullptr && opt != nullptr);

  // do nothing if we don't have the stride
  assert(stride >= 1);
  if(once_stopped_ || (iter % stride != 0))
    return;

  double cand = *value;
  if(one_scale)
    cand = (1-start) + ((end-(1-start))/end) * cand;

  switch(method_)
  {
  case NO_METHOD:
    assert(false);
    break;
  case ADD:
    cand += grow;
    break;
  case MULT:
    cand *= grow;
    break;
  case OBJ:
    if(opt->objectives.GetHistorySize() >= 2)
    {
      double f_p = opt->objectives.GetHistoryValue(true, -2);
      double f_k = opt->objectives.GetHistoryValue(true, -1);

      // Peter Dunning's formula. Don't divide by 0
      double test = f_p != f_k ? std::max((-grow / 2.0) * ((f_k + f_p) / (f_k - f_p)), 0.0) : 0.0;

      cand += std::min(test,max_grow_rate * cand);
      LOG_DBG(tune) << "U: iter=" << iter << " f_p=" << f_p << " f_k=" << f_k << " test=" << test;
    }
    break;
  }

  double rescaled = cand;
  if(one_scale)
    rescaled = (cand - (1-start))*end / (end-(1-start));


  LOG_DBG(tune) << "U: iter=" << iter << " m=" << method.ToString(method_) << " old=" << *value << " rescaled=" << rescaled << " end=" << end
                << " SG=" << SufficientlyGray() << " v=" << value << " cand=" << cand << " os=" << one_scale;

  if(rescaled > end)
  {
    once_stopped_ = true;

    LOG_DBG(tune) << "U: rescaled=" << rescaled << " end=" << end << " SG=" << SufficientlyGray() << " -> once_stopped_=true";
    return;
  }

  if(SufficientlyGray() && iter > 2) // small initial density can lead to too small grayness in the beginning
  {
    LOG_DBG(tune) << "U: rescaled=" << rescaled << " SG=" << SufficientlyGray() << " iter=" << iter << " -> once_stopped_=true";
    if(minimal == OFF || rescaled > minimal)
    {
      once_stopped_ = true;
      LOG_DBG(tune) << "U: SG=" << SufficientlyGray() << " minimal=" << minimal << " rescaled=" << rescaled << " -> once_stopped_=true";
      return;
    }
  }

  // do it: set the candidate to the value
  *value = rescaled;
  for(double* ptr : external)
    *ptr = rescaled;

  for(GlobalFilter* f : gf)
  {
    DesignSpace::DesignRegion* dr = opt->GetDesign()->GetRegion(f->region, f->design);
    f->SetNonLinCorrection(&(opt->GetDesign()->data[dr->base]));
    LOG_DBG(tune) << "U ref=" << dr->base << " b=" << f->beta << " e=" << f->eta << " SNLC -> scale=" << f->non_lin_scale << " offset=" << f->non_lin_offset << " f=" << f << " gf=" << gf.ToString(TS_PYTHON);
  }

  opt->SetAuxLogValue(usage.ToString(usage_), *value);
}

void Tune::FindGraynessStoppingRule()
{
  assert(opt != nullptr);

  // store the last found stopping rule containing greyness
  for(auto& rule : opt->objectives.stop)
    if(rule.function.find("greyness") != std::string::npos)
      grayness = &rule;

  LOG_DBG(tune) << "FGSR: grayness=" << (grayness ? grayness->function : "-");
}

bool Tune::SufficientlyGray()
{
  if(grayness == nullptr)
    return false;

  Condition* g = (Condition*) opt->constraints.Get(grayness->function, true); // throw exception
  double val = g->GetValue();
  double bound = grayness->value;
  assert(grayness->GetType() == grayness->BELOW_FUNCTION);

  LOG_DBG(tune) << "SG: g=" << g->ToString() << " v=" << val << " b=" << bound;
  return val <= bound;
}

} // end of namespace

