#include <assert.h>
#include <cmath>
#include <limits>
#include <ostream>

#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/Logging/log.hpp"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Design/Filter.hh"
#include "Utils/StdVector.hh"
#include "Utils/tools.hh"

DECLARE_LOG(ds)

using namespace CoupledField;

Filter::Filter()
{
  type_          = NO_FILTERING;
  sensitivity_   = SIGMUND;
  density_       = STANDARD;
  beta           = -1.0;
  eta            = -1.0;
  non_lin_scale  = -1.0;
  non_lin_offset = -1.0;
  explicit_lower_bound_ = std::numeric_limits<Double>::min();
  region         = NO_REGION_ID;
  weight         = 1.0;
  weight_sum     = -1.0;
}


Double Filter::CalcWeightSum(bool include_this) const
{
Double res = 0.0;

for(unsigned int i = 0, n = neighborhood.GetSize(); i < n; i++)
  res += neighborhood[i].weight;

if(include_this)
  res += this->weight;

return res;
}


Double Filter::GetLowerBound(const DesignElement* de) const
{
  if(explicit_lower_bound_ != std::numeric_limits<Double>::min())
    return explicit_lower_bound_;
  else
    return de->GetLowerBound();
}


void Filter::SetNonLinCorrection(const DesignElement* ref, unsigned int fix)
{
  if(type_ != DENSITY || !(density_ == SOLID_HEAVISIDE || density_ == VOID_HEAVISIDE || density_ == TANH))
    return;

  assert(beta >= 0);
  assert(density_ != TANH || eta >= 0);

  // we ignore the bimat case - the we still have the tahn problem for small beta

  DesignElement de = DesignElement();
  de.simp = new SIMPElement(&de);
  de.simp->filter = *this; // only copy data, changes don't distribute back
  de.elem = ref->elem; // otherwise the logings segfault

  Double ub = ref->GetUpperBound();
  Double lb = ref->GetLowerBound();

  // calc the pure values
  de.simp->filter[fix].non_lin_scale = 1.0;
  de.simp->filter[fix].non_lin_offset = 0;

  Double org_u = density_ == TANH ? de.simp->CalcTanh(ub, fix) : de.simp->CalcHeaviside(ub, fix);
  Double org_l = density_ == TANH ? de.simp->CalcTanh(lb, fix) : de.simp->CalcHeaviside(lb, fix);

  LOG_DBG(ds) << "SNLC de=" << ref->ToString() << " fix=" << fix << " f=" << density.ToString(density_) << " ub=" << ub << " -> " << org_u << " lb=" << lb << " -> " << org_l;

  assert((org_u-org_l) >= 0.1);
  assert((org_u-org_l) <= 1.0);
  assert(org_u > org_l);

  // F ist the filter. We scale with scale*F + offset
  // scale such that F(u)-F(l) == ub-lb
  this->non_lin_scale  = (ub-lb) / (org_u - org_l);
  this->non_lin_offset = lb - non_lin_scale * org_l;
  de.simp->filter = *this; // copy again for the asserts checking the results

  LOG_DBG(ds) << "SNLC de=" << ref->ToString() << " f=" << density.ToString(density_) << " d==" << (ub-lb) << " od=" << (org_u-org_l) << " -> s=" << non_lin_scale << " o=" << non_lin_offset;

  assert(close(density_ == TANH ? de.simp->CalcTanh(lb, fix) : de.simp->CalcHeaviside(lb, fix), lb));
  assert(close(density_ == TANH ? de.simp->CalcTanh(ub, fix) : de.simp->CalcHeaviside(ub, fix), ub));
}


void Filter::Dump() const
{
  Double weight_sum = weight;
  Double distance_avg = 0.0;
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
