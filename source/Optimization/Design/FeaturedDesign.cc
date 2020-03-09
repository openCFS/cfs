#include "DataInOut/Logging/LogConfigurator.hh"
#include "Optimization/Design/FeaturedDesign.hh"

namespace CoupledField {

DEFINE_LOG(FD, "featuredDesign")

Enum<FeaturedDesign::IntStrategy> FeaturedDesign::intStrategy_;
unsigned int FeaturedDesign::dim_ = 99;

FeaturedDesign::FeaturedDesign(StdVector<RegionIdType>& regionIds, PtrParamNode pn, ErsatzMaterial::Method method)
: AuxDesign(regionIds, pn, method)
{
  intStrategy_.SetName("FeaturedDesign::IntStrategy");
  intStrategy_.Add(CONSTANT_FULL, "constant_full");
  intStrategy_.Add(FULL_OR_NOTHING, "full_or_nothing");

  this->dim_ = domain->GetGrid()->GetDim();
  this->export_fe_design_ = false; // we use the original design but don't communicate it via ReadDesignFromExtern(), ...
  this->tailing_aux_design_ = true; // we want our own design (e.g. SBD: param_, SMD: shape_param_ or better their opt_* versions) to take the role of DesignSpace::data

  this->mapping_timer_  = info_->Get("splineBox/mapping/timer")->AsTimer();
  this->mapping_timer_->SetLabel("splineBox_map");
  this->mapping_timer_->SetSub(); // already in eval_*
  this->gradient_timer_ = info_->Get("splineBox/gradient/timer")->AsTimer();
  this->gradient_timer_->SetLabel("splineBox_grad");
  this->gradient_timer_->SetSub(); // already in eval_*
}


double FeaturedDesign::Item::MaxDiffCornerValue() const
{
  assert(min_corner_value.GetSize() == max_corner_value.GetSize());
  assert(min_corner_value.GetSize() > 0);

  double diff = 0;

  for(unsigned int i = 0; i < min_corner_value.GetSize(); i++)
    diff = std::max(max_corner_value[i] - min_corner_value[i], diff);

  return diff;
}


int FeaturedDesign::Item::GetOrder(Vector<int>& order, const FeaturedDesign::NumInt& ni) const
{
  assert(order.GetSize() == min_corner_value.GetSize());
  assert(order.GetSize() == max_corner_value.GetSize());
  int max = 0;
  for(unsigned int s = 0; s < order.GetSize(); s++)
  {
    double min_val = min_corner_value[s];
    double max_val = max_corner_value[s];

    assert(min_val >= 0);
    assert(max_val >= min_val);
    assert(max_val <= 1);

    // prevent unused variable warning
    (void)(min_val);

    switch(ni.strategy_)
    {
    case FeaturedDesign::CONSTANT_FULL:
      order[s] = ni.max_order_;
      break;
    case FeaturedDesign::FULL_OR_NOTHING:
      order[s] = max_val < ni.sensitivity_ ? 0 : ni.max_order_;
      break;
    case FeaturedDesign::TAILORED:
      throw Exception("Not implemented");
      break;
    }

    max = std::max(max, order[s]);
  }
  return max;
}


void FeaturedDesign::NumInt::Init(FeaturedDesign* sbd, PtrParamNode pn, PtrParamNode info)
{
  this->sensitivity_ = pn->Get("sensitivity")->As<double>();
  assert(sensitivity_ > 0);
  this->max_order_ = pn->Get("integration_order")->As<int>();

  if(max_order_ < 1)
    info->SetWarning("minimal value for 'max_order' is '1'");

  this->strategy_ = intStrategy_.Parse(pn->Get("integration_strategy")->As<string>());

  cells_ = domain->GetGrid()->GetNumElems();
}


void FeaturedDesign::NumInt::ToInfo(PtrParamNode info) const
{
  info->Get("max_order")->SetValue(max_order_);
  info->Get("sensitivity")->SetValue(sensitivity_);
  info->Get("integration")->SetValue(FeaturedDesign::intStrategy_.ToString(strategy_));
  assert(cells_ > 0);
  PtrParamNode cells = info->Get("cells");
  cells->Get("integrate_fraction")->SetValue(int_cells_cnt_ / (double) cells_);
  cells->Get("avg_order")->SetValue(int_cells_order_sum_ / (double) int_cells_cnt_);
  cells->Get("total_int")->SetValue(std::pow(int_cells_order_sum_, domain->GetGrid()->GetDim()));
}


} // end of namespace


