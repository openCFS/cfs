#include <assert.h>
#include <cmath>
#include <limits>
#include <ostream>
#include <sstream>

#include "DataInOut/Logging/LogConfigurator.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Design/Filter.hh"
#include "Optimization/Tune.hh"
#include "Utils/StdVector.hh"
#include "Utils/tools.hh"
#include "Utils/mathParser/mathParser.hh"

using namespace CoupledField;

DEFINE_LOG(flt, "filter")

void Filter::Dump() const
{
  double weight_sum = weight;
  double distance_avg = 0.0;
  for(unsigned int i = 0; i < neighborhood.GetSize(); i++)
  {
    weight_sum   += neighborhood[i].weight;
    distance_avg += neighborhood[i].distance;
  }
  distance_avg /= neighborhood.GetSize();

  std::cout << " weight sum " << weight_sum << " this weight " << weight <<" distance avg " << distance_avg << std::endl;
  for(unsigned int i = 0; i < neighborhood.GetSize(); i++)
    std::cout << "  n[" << i << "]: elem " << neighborhood[i].neighbour->elem->elemNum << " location "
              << neighborhood[i].neighbour->GetLocation()->ToString()
              << " dist=" << neighborhood[i].distance << " w=" << neighborhood[i].weight << std::endl;
}

GlobalFilter::GlobalFilter()
{
  LOG_DBG(flt) << "GF &beta=" << &beta << " tune=" << tune.IsSet() << " e_t=" << ext_tune;
}

GlobalFilter::~GlobalFilter()
{
  LOG_DBG(flt) << "~GF &beta=" << &beta << " tune=" << tune.IsSet() << " e_t=" << ext_tune;

  if(ext_tune != nullptr && ext_tune != (Tune*)(0) + UNSET_EXT_BETA_TUNE) // prevent conversion of different size
    ext_tune->Remove(&beta, this);
}

void GlobalFilter::PostCopy(bool register_tune)
{
  if(density == Filter::EXPRESSION)
    ReInitParser();

  if(register_tune && domain->GetOptimization() != nullptr) // handle pure load_ersatzmaterial
  {
    Tune* t = domain->GetOptimization()->SearchTune(Tune::BETA, true); // silent
    if(tune.IsSet())
    {
      assert(ext_tune == nullptr);
      if(t)
        throw Exception("beta tune already registered. In the robust case, the first filter shall have 'tune' beta, the other ones 'tuned'");
      tune.Register(&beta, domain->GetOptimization(), this);
    }
    else if(ext_tune == (Tune*)(0) + UNSET_EXT_BETA_TUNE) // prevent conversion of different size
    {
      if(t == nullptr)
        throw Exception("No beta tune registered. In the robust case, the first filter shall have 'tune' beta, the other ones 'tuned'");
      t->Append(&beta, this);
      ext_tune = t;
    }
  }
}

void GlobalFilter::SetNonLinCorrection(const DesignElement* ref)
{
  non_lin_scale = 1.0;
  non_lin_offset = 0;
  LOG_DBG(flt) << "GF:SNLC de=" << ref->ToString() << " this=" << ToString();

  if(type != Filter::DENSITY || density == Filter::STANDARD)
    return;

  assert(density != Filter::MATERIAL_PART);
  assert(density == Filter::MATERIAL || beta >= 0);
  assert(density != Filter::TANH || eta >= 0);

  // we ignore the bimat case - there we still have the tanh problem for small beta
  // we create a temporary global filter
  GlobalFilter tmp = *this;
  tmp.PostCopy(false); // we need no tune

  double lb = ref->GetLowerBound();
  double ub = ref->GetUpperBound();

  // in the material filter case we have three transfer functions involved and this is the filtered lower bound
  // in all other cases the filtered lower bound is simply the lower bound
  double flb = density == Filter::MATERIAL ? mat_filter.Transform(lb) : lb;

  // calc the pure values
  tmp.non_lin_scale = 1.0;
  tmp.non_lin_offset = 0;

  double org_u = SIMPElement::CalcNonLinFilter(ub, &tmp,ub);
  double org_l = SIMPElement::CalcNonLinFilter(flb, &tmp,lb);

  LOG_DBG(flt) << "GF:SNLC de=" << ref->ToString() << " ub=" << ub << " -> " << org_u << " lb=" << lb << " flb=" << flb << " -> " << org_l << " f=" << ToString();

  assert((org_u-org_l) >= 0.1);
  assert((org_u-org_l) <= 1.0);
  assert(org_u > org_l);

  // F is the filter. We scale with scale*F + offset
  // scale such that F(u)-F(l) == ub-lb
  non_lin_scale  = (ub-lb) / (org_u - org_l);
  non_lin_offset = lb - non_lin_scale * org_l;

  LOG_DBG(flt) << "GF:SNLC de=" << ref->ToString() << " f=" << Filter::density.ToString(density) << " d==" << (ub-lb) << " od=" << (org_u-org_l) << " -> s=" << non_lin_scale << " o=" << non_lin_offset;
  LOG_DBG(flt) << "GF:SNLC " << Filter::density.ToString(density) << "(" << ub << ")="<< SIMPElement::CalcNonLinFilter(ub, this, ub);
  LOG_DBG(flt) << "GF:SNLC " << Filter::density.ToString(density) << "(" << flb << ")="<< SIMPElement::CalcNonLinFilter(flb, this, lb);
  assert(close(SIMPElement::CalcNonLinFilter(ub, this, ub), ub));
  assert(close(SIMPElement::CalcNonLinFilter(flb, this, lb), lb));
}

std::string GlobalFilter::ToString() const
{
  std::stringstream ss;
  ss << "dt=" << DesignElement::type.ToString(design);
  ss << " reg=" << region << " rob=" << robust;
  ss << " t=" << Filter::type.ToString(type);
  if(type == Filter::DENSITY)
    ss << " df=" << Filter::density.ToString(density);
  else
    ss << " ss=" << Filter::sensitivity.ToString(sensitivity);
   ss << " v=" << value;
  return ss.str();
}

void GlobalFilter::ToInfo(PtrParamNode in)
{
  in->Get("type")->SetValue(Filter::filterSpace.ToString(filterspace));

  in->Get("value")->SetValue(value);
  in->Get("contribution")->SetValue(contribution == Filter::LINEAR ? "linear" : "constant");

  in->Get("design")->SetValue(DesignElement::type.ToString(design));

  if(type == Filter::SENSITIVITY)
    in->Get("sensitivity")->SetValue(Filter::sensitivity.ToString(sensitivity));

  if(type == Filter::DENSITY && density != Filter::STANDARD)
  {
    PtrParamNode in_ = in->Get(Filter::density.ToString(density));

    switch(density)
    {
    case Filter::EXPRESSION:
      in_->Get("func")->SetValue(parser_->GetExpr(function_handle_));
      in_->Get("grad")->SetValue(parser_->GetExpr(derivative_handle_));
      // intentionally no break
    case Filter::TANH:
      in_->Get("eta")->SetValue(eta);
      // intentionally no break
    case Filter::VOID_HEAVISIDE:
    case Filter::SOLID_HEAVISIDE:
      in_->Get("beta")->SetValue(tune.IsSet() ? "tune" : std::to_string(beta));
      if(tune.IsSet())
        tune.ToInfo(in_);
      break;
    case Filter::MATERIAL:
      mat_filter.ToInfo(in_->Get("filter"),true); // skip design and app
      mat_scale.ToInfo(in_->Get("scale"),true);
      mat_phase.ToInfo(in_->Get("phase"),true);
      break;
    case Filter::STANDARD:
    case Filter::MATERIAL_PART:
      break;
    }
    in_->Get("scaling")->SetValue(non_lin_scale);
    in_->Get("offset")->SetValue(non_lin_offset);
  }
  if(robust >= 0)
    in->Get("robust_excitation")->SetValue(robust);

  in->Get("elems")->SetValue(elements);
  in->Get("avg_radius")->SetValue(avg_radius);
  in->Get("avg_neighbors")->SetValue(avg_neigbor);
}

void GlobalFilter::ReInitParser()
{
  assert(parser_ != nullptr && function_handle_ > 0 && derivative_handle_ > 0);
  string func = parser_->GetExpr(function_handle_);
  string grad = parser_->GetExpr(derivative_handle_);
  assert(func.size() > 0 && grad.size() > 0);
  InitParser(func, grad); // intentionally do not delete old handles
}

void GlobalFilter::InitParser(const string& func_expr, const string& deriv_expr)
{
  // note that the handles can be alredy set when calling ReInitParser()
  parser_ = domain->GetMathParser();
  function_handle_  = parser_->GetNewHandle();
  derivative_handle_= parser_->GetNewHandle();
  parser_->SetValue(function_handle_, "r", value); // the filter value is the radius
  parser_->RegisterExternalVar(function_handle_, "b", &beta);
  parser_->RegisterExternalVar(function_handle_, "e", &eta); // the Euler number is _e
  parser_->RegisterExternalVar(function_handle_, "x", &expression_x_); // the density filtered input density
  parser_->SetValue(derivative_handle_, "r", value);
  parser_->RegisterExternalVar(derivative_handle_, "b", &beta);
  parser_->RegisterExternalVar(derivative_handle_, "e", &eta);
  parser_->RegisterExternalVar(derivative_handle_, "x", &expression_x_);
  parser_->SetExpr(function_handle_, func_expr);
  parser_->SetExpr(derivative_handle_, deriv_expr);
  LOG_DBG(flt) << "GF::IP fh=" << function_handle_ << " dh=" << derivative_handle_ << " vals: " << parser_->ToString(function_handle_) << " &ex=" << &expression_x_;
}

double GlobalFilter::EvalExpressionFunction(double x, bool scale) const
{
  assert(parser_ != nullptr);
  double res = 0;
  LOG_DBG3(flt) << "GF:EEF expression_x_=" << expression_x_ << " reg variables=" << parser_->GetRegisteredVariables(function_handle_);
  #pragma omp critical
  {
    expression_x_ = x;
    res = parser_->Eval(function_handle_);
  }
  if(scale)
    res = (non_lin_scale * res) + non_lin_offset;
  //result = function_parser_.Eval();
  LOG_DBG3(flt) << "GF:EEF x=" << x << " b=" << beta << " e=" << eta << " &ex=" << &expression_x_ << " fh=" << function_handle_ << " s=" << scale << " -> " << res;
  return res;
}

double GlobalFilter::EvalExpressionDerivative(double x, bool scale) const
{
  double res = 0;
  #pragma omp critical
  {
    expression_x_ = x;
    res = parser_->Eval(derivative_handle_);
  }
  if(scale)
    res *= non_lin_scale;
  //result = function_parser_.Eval();
  LOG_DBG3(flt) << "GF:EED x=" << x << " s=" << scale << " -> " << res;
  return res;
}

