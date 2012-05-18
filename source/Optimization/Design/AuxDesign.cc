#include <assert.h>
#include <ostream>
#include <string>

#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"
#include "General/exception.hh"
#include "Optimization/Design/AuxDesign.hh"

namespace CoupledField {

class Condition;
class Objective;

DECLARE_LOG(aux_des)
DEFINE_LOG(aux_des, "auxDesign")

AuxDesign::AuxDesign(StdVector<RegionIdType>& regions,  PtrParamNode pn, ErsatzMaterial::Method method, unsigned int naux)
  : DesignSpace(regions, pn, method)
{
  alsomatopt_ = true; // can be different in ShapeDesign
  scaling_ = 1.0;

  assert(naux == 0 || naux == 1);
  if(naux == 1)
  {
    slack_  = pn->GetByVal("design", "name", "slack");
    aux_design_.Reserve(1);
    BaseDesignElement de;
    de.SetLowerBound(slack_->Get("lower")->As<double>());
    de.SetUpperBound(slack_->Get("upper")->As<double>());
    de.SetDesign(slack_->Get("initial")->As<double>());
    // we need to PostInit!
    aux_design_.Push_back(de);
  }
}

void AuxDesign::PostInit(int objectives, int constraints)
{
  DesignSpace::PostInit(objectives, constraints);

  assert((slack_ != NULL && aux_design_.GetSize() == 1) || slack_ == NULL);

  if(slack_ != NULL)
  {
    LOG_DBG(aux_des) << "PI: #objectives = " << objectives << ", #constraints = " << constraints;
    DesignElement::SetDesignSpace(this);
    assert(aux_design_.GetSize() == 1);
    assert(objectives == 1);
    aux_design_[0].PostInit(objectives, constraints);
  }
}

void AuxDesign::ToInfo(PtrParamNode in)
{
  DesignSpace::ToInfo(in);

  if(slack_ != NULL)
  {
    PtrParamNode in_ = in->Get("designVariables")->Get("design", ParamNode::APPEND);
    in_->Get("type")->SetValue("slack");
    in_->Get("upperBound")->SetValue(aux_design_[0].GetUpperBound());
    in_->Get("lowerBound")->SetValue(aux_design_[0].GetLowerBound());
  }
}

int AuxDesign::ReadDesignFromExtern(const double* space_in)
{
  int old_design = design_id;

  if(alsomatopt_)
    DesignSpace::ReadDesignFromExtern(space_in);

  unsigned int offset = DesignSpace::GetNumberOfVariables(); // the size of the simp space - might be != data.GetSize()
  assert((alsomatopt_ && (offset <= elements * design.GetSize())) || (!alsomatopt_ && offset == 0));

  // might be set above
  bool new_design = old_design != design_id;

  for(unsigned int i=0; i < aux_design_.GetSize(); i++)
  {
    double v = space_in[offset + i] * scaling_;
    if(!new_design && v != aux_design_[i].GetDesign())
      new_design = true;

    aux_design_[i].SetDesign(v);
    LOG_DBG(aux_des) << "ReadDesignFromExtern: shapeparams_[i]=" << v;
  }
  if(new_design && design_id <= old_design) design_id++; // if new design and not already changed by DesignSpace
  return(design_id);
}


bool AuxDesign::CompareDesign(const double* space_in)
{
  if(!DesignSpace::CompareDesign(space_in))
    return false;

  unsigned int offset = DesignSpace::GetNumberOfVariables();

  for(unsigned int i=0; i < aux_design_.GetSize(); i++)
  {
    double v = space_in[offset + i] * scaling_;
    if(v != aux_design_[i].GetDesign())
      return false;
  }

  return true;
}

int AuxDesign::WriteDesignToExtern(double* space_out, bool scale) const
{
  if(alsomatopt_)
    DesignSpace::WriteDesignToExtern(space_out, scale);

  double rscaling = scale ? 1.0 / scaling_ : 1.0;
  unsigned int offset = DesignSpace::GetNumberOfVariables();

  for(unsigned int i=0; i < aux_design_.GetSize(); i++)
  {
    space_out[offset + i] = aux_design_[i].GetDesign() * rscaling;
    LOG_DBG(aux_des) << "WriteDesignToExtern: out[" << i << "]=" << space_out[offset +i];
  }
  return design_id;
}

void AuxDesign::WriteGradientToExtern(StdVector<double>& out, DesignElement::ValueSpecifier vs, DesignElement::Access access, Condition* g, bool scaling) const
{
  LOG_DBG(aux_des) << "WGTE: ad=" << aux_design_.GetSize() << " DS:GNOV=" << DesignSpace::GetNumberOfVariables() << " ows=" << out.window.GetSize();

  if(alsomatopt_)
  {
    // we call DesignSpace::WriteDenseGradientToExtern() for the ersatz material part.

    // in case the the gradient window covers DesignSpace design and aux design, we need to
    // modify the window before calling DesignSpace::Write*GradientToExtern().
    // this is e.g. not the case for sparse gradients not touching the aux design.
    StdVector<double>::Window org_window = out.window; // I like standard constructors :)

    if(out.window.GetStart() + out.window.GetSize() > data.GetSize())
      out.window.Set(out.window.GetStart(), out.window.GetSize() - aux_design_.GetSize());

    // The original window needs to be restored afterwards for assert() in BaseOptimization which
    // does not know about separation into shape and mat variables
    if(out.window.GetSize() != org_window.GetSize()) // restore original window
      out.window = org_window;
  }

  // makes use of the window within out even  if only a part of the window is used in the alsomatopt_ case
  WriteAuxGradientToExtern(out, g, scaling);
}


void AuxDesign::WriteAuxGradientToExtern(StdVector<double>& out, Condition* g, bool scale) const
{
  unsigned int base = out.window.GetStart()  + out.window.GetSize() - aux_design_.GetSize();

  LOG_DBG(aux_des) << "WAGTE: g=" << (g == NULL ? "null" : g->ToString()) << " ows=" << out.window.GetStart()
                    << " HS=" << HasSlackVariable() << " base=" << base;

  assert(aux_design_.GetSize() <= out.window.GetSize());
  double s = scale ? scaling_ : 1.0;
  for(unsigned int i=0; i < aux_design_.GetSize(); i++)
  {
    if(HasSlackVariable() && g != NULL && g->HasSlackBound())
      out[base + i] = -1.0;
    else
      out[base + i] = aux_design_[i].GetPlainGradient(NULL, g) * s;
    LOG_DBG3(aux_des) << "WAGTE: g=" << (g == NULL ? "null" : g->ToString()) << " out[" << base+i << "]=" << out[base+i]
                      << " slack case=" << (HasSlackVariable() && g != NULL && g->HasSlackBound());
  }
}

void AuxDesign::Reset(DesignElement::ValueSpecifier vs, DesignElement::Type design)
{
  for(unsigned int i=0; i < aux_design_.GetSize(); i++)
    aux_design_[i].Reset(vs);

  DesignSpace::Reset(vs, design);
}

void AuxDesign::WriteBoundsToExtern(double* x_l, double* x_u) const
{
  if(alsomatopt_)
    DesignSpace::WriteBoundsToExtern(x_l, x_u);

  unsigned int offset = DesignSpace::GetNumberOfVariables();
  assert(offset + aux_design_.GetSize() == GetNumberOfVariables());

  for(unsigned int i=0; i < aux_design_.GetSize(); i++)
  {
    x_l[offset + i] = aux_design_[i].GetLowerBound() / scaling_;
    x_u[offset + i] = aux_design_[i].GetUpperBound() / scaling_;
    LOG_DBG3(aux_des) << "WBTE: l[" << (offset + i) << "]=" << x_l[offset + i] << " u[" << (offset + i) << "]=" << x_u[offset + i];
  }
}

unsigned int AuxDesign::GetNumberOfVariables() const
{
  if(alsomatopt_){
    return(aux_design_.GetSize() + DesignSpace::GetNumberOfVariables());
  }else{
    return(aux_design_.GetSize());
  }
}

void AuxDesign::AddAuxDerivative(Function* f, unsigned int index, double value)
{
  aux_design_[index].AddGradient(f, value);
}

void AuxDesign::AddAuxDerivatives(Objective* f, Condition* g, StdVector<double>& d, double weight)
{
  assert(d.GetSize() == aux_design_.GetSize());
  for(unsigned int i = 0; i < aux_design_.GetSize(); i++){
    LOG_DBG3(aux_des) << "AAD[" << i << "]+=" << d[i] << "*" << weight << "=" << d[i]*weight;
    aux_design_[i].AddGradient(f, g, d[i]*weight);
  }
}


} // end of namespace
