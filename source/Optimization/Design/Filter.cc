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
  beta_          = -1.0;
  eta            = -1.0;
  non_lin_scale  = -1.0;
  non_lin_offset = -1.0;
  explicit_lower_bound_ = std::numeric_limits<double>::min();
}

double Filter::GetLowerBound(const DesignElement* de) const
{
  if(explicit_lower_bound_ != std::numeric_limits<double>::min())
    return explicit_lower_bound_;
  else
    return de->GetLowerBound();
}


void Filter::SetNonLinCorrection(const DesignElement* ref)
{
  if(type_ != DENSITY || !(density_ == SOLID_HEAVISIDE || density_ == VOID_HEAVISIDE || density_ == TANH))
    return;

  assert(beta_ >= 0);
  assert(density_ != TANH || eta >= 0);

  // we ignore the bimat case - the we still have the tahn problem for small beta

  DesignElement de = DesignElement();
  de.simp = new SIMPElement(&de);
  de.simp->filter = *this; // this is the important point!. Note that we copy the Filter, not the pointer. Any change on this will NOT be reflected in de.simp->filter!
  de.elem = ref->elem; // otherwise the logings segfault

  double ub = ref->GetUpperBound();
  double lb = ref->GetLowerBound();

  // calc the pure values
  de.simp->filter.non_lin_scale = 1.0;
  de.simp->filter.non_lin_offset = 0;

  double org_u = density_ == TANH ? de.simp->CalcTanh(ub) : de.simp->CalcHeaviside(ub);
  double org_l = density_ == TANH ? de.simp->CalcTanh(lb) : de.simp->CalcHeaviside(lb);

  LOG_DBG(ds) << "SNLC de=" << ref->ToString() << " f=" << density.ToString(density_) << " ub=" << ub << " -> " << org_u << " lb=" << lb << " -> " << org_l;

  assert((org_u-org_l) >= 0.1);
  assert((org_u-org_l) <= 1.0);
  assert(org_u > org_l);

  // F ist the filter. We scale with scale*F + offset
  // scale such that F(u)-F(l) == ub-lb
  non_lin_scale  = (ub-lb) / (org_u - org_l);
  non_lin_offset = lb - non_lin_scale * org_l;

  LOG_DBG(ds) << "SNLC de=" << ref->ToString() << " f=" << density.ToString(density_) << " d==" << (ub-lb) << " od=" << (org_u-org_l) << " -> s=" << non_lin_scale << " o=" << non_lin_offset;

  de.simp->filter = *this; // this is only for the asserts
  assert(close(density_ == TANH ? de.simp->CalcTanh(lb) : de.simp->CalcHeaviside(lb), lb));
  assert(close(density_ == TANH ? de.simp->CalcTanh(ub) : de.simp->CalcHeaviside(ub), ub));
}
