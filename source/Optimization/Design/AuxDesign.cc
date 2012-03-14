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
DEFINE_LOG(aux_des, "aux_design")

AuxDesign::AuxDesign(StdVector<RegionIdType>& regions,  PtrParamNode pn, ErsatzMaterial::Method method)
  : DesignSpace(regions, pn, method)
{
  alsomatopt_ = true; // can be different in ShapeDesign
  scaling_ = 1.0;
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

void AuxDesign::WriteDenseGradientToExtern(StdVector<double>& out, DesignElement::ValueSpecifier vs, DesignElement::Access access, Condition* g, bool scaling) const
{
  assert(aux_design_.GetSize() + DesignSpace::GetNumberOfVariables() == out.window.GetSize());

  if(alsomatopt_)
  {
    // we call DesignSpace::WriteDenseGradientToExtern() for the ersatz material part.
    // therefore the window needs to be adjusted to a "subwindow" as WriteDenseGradientToExtern()
    // orientates itself on the given window.
    // The original window needs to be restored afterwards for assert() in BaseOptimization which
    // does not know about separation into shape and mat variables
    StdVector<double>::Window org_window = out.window; // I like standard constructors :)
    // shift the window to do the rest
    // replace the original window by a subwindow excluding the shape stuff
    out.window.Set(out.window.GetStart(), out.window.GetSize() - aux_design_.GetSize());
    DesignSpace::WriteDenseGradientToExtern(out, vs, access, g, scaling);

    // restore original window
    out.window = org_window;
  }

  // makes use of the window within out even  if only a part of the window is used in the alsomatopt_ case
  WriteAuxGradientToExtern(out, g, scaling);
}


void AuxDesign::WriteAuxGradientToExtern(StdVector<double>& out, Condition* g, bool scale) const
{
  unsigned int base = out.window.GetStart() + DesignSpace::GetNumberOfVariables();
  assert(aux_design_.GetSize() <= out.window.GetSize());
  double s = scale ? scaling_ : 1.0;
  for(unsigned int i=0; i < aux_design_.GetSize(); i++)
  {
    out[base + i] = aux_design_[i].GetPlainGradient(NULL, g) * s;
    LOG_DBG3(aux_des) << "WAGTE: out[" << base+i << "]=" << out[base+i];
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

void AuxDesign::AddAuxDerivatives(Objective* f, Condition* g, StdVector<double>& d, double weight)
{
  assert(d.GetSize() == aux_design_.GetSize());
  for(unsigned int i = 0; i < aux_design_.GetSize(); i++){
    LOG_DBG3(aux_des) << "AAD[" << i << "]+=" << d[i] << "*" << weight << "=" << d[i]*weight;
    aux_design_[i].AddGradient(f, g, d[i]*weight);
  }
}


} // end of namespace
