#include <assert.h>
#include <cmath>
#include <limits>
#include <ostream>
#include <sstream>

#include "DataInOut/Logging/LogConfigurator.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Design/Filter.hh"
#include "Utils/StdVector.hh"
#include "Utils/tools.hh"

using namespace CoupledField;

EXTERN_LOG(ds)

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

void GlobalFilter::SetNonLinCorrection(const DesignElement* ref)
{
  non_lin_scale = 1.0;
  non_lin_offset = 0;
  LOG_DBG(ds) << "GF:SNLC de=" << ref->ToString() << " this=" << ToString();

  if(type != Filter::DENSITY || density == Filter::STANDARD)
    return;

  assert(density != Filter::MATERIAL_PART);
  assert(density == Filter::MATERIAL || beta >= 0);
  assert(density != Filter::TANH || eta >= 0);

  // we ignore the bimat case - there we still have the tanh problem for small beta

  // we create a temporary global filter
  GlobalFilter tmp = *this;

  double ub = ref->GetUpperBound();
  double lb = ref->GetLowerBound();

  // in the material filter case we have three transfer functions involved and this is the filtered lower bound
  // in all other cases the filtered lower bound is simply the lower bound
  double flb = density == Filter::MATERIAL ? mat_filter.Transform(lb) : lb;

  // calc the pure values
  tmp.non_lin_scale = 1.0;
  tmp.non_lin_offset = 0;

  double org_u = SIMPElement::CalcNonLinFilter(ub, &tmp,ub);
  double org_l = SIMPElement::CalcNonLinFilter(flb, &tmp,lb);

  LOG_DBG(ds) << "GF:SNLC de=" << ref->ToString() << " ub=" << ub << " -> " << org_u << " lb=" << lb << " flb=" << flb << " -> " << org_l << " f=" << ToString();

  assert((org_u-org_l) >= 0.1);
  assert((org_u-org_l) <= 1.0);
  assert(org_u > org_l);

  // F ist the filter. We scale with scale*F + offset
  // scale such that F(u)-F(l) == ub-lb
  non_lin_scale  = (ub-lb) / (org_u - org_l);
  non_lin_offset = lb - non_lin_scale * org_l;

  LOG_DBG(ds) << "GF:SNLC de=" << ref->ToString() << " f=" << Filter::density.ToString(density) << " d==" << (ub-lb) << " od=" << (org_u-org_l) << " -> s=" << non_lin_scale << " o=" << non_lin_offset;

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
    case Filter::TANH:
      in_->Get("eta")->SetValue(eta);
    case Filter::VOID_HEAVISIDE:
    case Filter::SOLID_HEAVISIDE:
      in_->Get("beta")->SetValue(beta);
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


