#include <assert.h>
#include <ostream>
#include <string>

#include "DataInOut/Logging/LogConfigurator.hh"
#include "General/Exception.hh"
#include "Optimization/Design/AuxDesign.hh"

namespace CoupledField {

class Condition;
class Objective;

DEFINE_LOG(aux_des, "auxDesign")

AuxDesign::AuxDesign(StdVector<RegionIdType>& regions,  PtrParamNode pn, ErsatzMaterial::Method method, int naux)
  : DesignSpace(regions, pn, method)
{
  setup_timer_->Start();

  export_fe_design_ = true; // can be different in ShapeDesign
  tailing_aux_design_ = true; // set false for shape optimization
  scaling_ = 1.0;

  if(naux == -1) { // check automatically
    naux = 0;
    if(pn->HasByVal("design", "name", "slack"))
      naux = 1;
    if(pn->HasByVal("design", "name", "alpha"))
      naux = 2;
  }
  assert(naux == 0 || naux == 1 || naux == 2);
  if(naux == 1 || naux == 2)
  {
    slack_  = pn->GetByVal("design", "name", "slack");
    aux_design_.Reserve(naux);
    BaseDesignElement de(BaseDesignElement::SLACK);
    de.SetLowerBound(slack_->Get("lower")->As<double>());
    de.SetUpperBound(slack_->Get("upper")->As<double>());
    de.SetDesign(slack_->Get("initial")->As<double>());
    de.SetIndex(0);  // We need to set the index to some value, because default is std::numeric_limits<unsigned int>::max()
    // we need to PostInit!
    aux_design_.Push_back(de);
  }

  if(naux == 2)
  {
    alpha_  = pn->GetByVal("design", "name", "alpha");
    BaseDesignElement de(BaseDesignElement::ALPHA);
    de.SetLowerBound(alpha_->Get("lower")->As<double>());
    de.SetUpperBound(alpha_->Get("upper")->As<double>());
    de.SetDesign(alpha_->Get("initial")->As<double>());
    aux_design_.Push_back(de);
  }

  setup_timer_->Stop();
}

void AuxDesign::PostInit(int objectives, int constraints)
{
  DesignSpace::PostInit(objectives, constraints);

  setup_timer_->Start();

  assert((slack_ != NULL && (aux_design_.GetSize() == 1 || aux_design_.GetSize() == 2)) || slack_ == NULL);
  assert((alpha_ != NULL && aux_design_.GetSize() == 2 && slack_ != NULL) || alpha_ == NULL);

  if(slack_ != NULL)
  {
    LOG_DBG(aux_des) << "PI: #objectives = " << objectives << ", #constraints = " << constraints;
    DesignElement::SetDesignSpace(this);
//    assert(objectives == 1);
    aux_design_[0].PostInit(objectives, constraints);
    if(alpha_ != NULL)
      aux_design_[1].PostInit(objectives, constraints);
  }

  // extend full_data
  unsigned int offset = AuxDesignOffset();
  full_data.Resize(offset + aux_design_.GetSize());
  for(unsigned int i = 0; i < aux_design_.GetSize(); i++)
    full_data[offset + i] = &(aux_design_[i]);

  setup_timer_->Stop();
}

void AuxDesign::ToInfo(ErsatzMaterial* em)
{
  DesignSpace::ToInfo(em);

  if(slack_ != NULL)
  {
    PtrParamNode in_ = info_->Get("designVariables")->Get("design", ParamNode::APPEND);
    in_->Get("type")->SetValue("slack");
    in_->Get("upperBound")->SetValue(aux_design_[0].GetUpperBound());
    in_->Get("lowerBound")->SetValue(aux_design_[0].GetLowerBound());
  }

  if(alpha_ != NULL)
  {
    PtrParamNode in_ = info_->Get("designVariables")->Get("design", ParamNode::APPEND);
    in_->Get("type")->SetValue("alpha");
    in_->Get("upperBound")->SetValue(aux_design_[1].GetUpperBound());
    in_->Get("lowerBound")->SetValue(aux_design_[1].GetLowerBound());
  }

}

int AuxDesign::ReadDesignFromExtern(const double* space_in, bool setAndWriteCurrent)
{
  int old_design = design_id;

  if(export_fe_design_)
    DesignSpace::ReadDesignFromExtern(space_in, setAndWriteCurrent);

  unsigned int offset = AuxDesignOffset(); // the size of the simp space - might be != data.GetSize()
  assert((export_fe_design_ && (offset <= elements * design.GetSize())) || !export_fe_design_);

  // design_id might be changed above in DesignSpace::ReadDesignFromExtern()
  bool new_design = old_design != design_id;

  for(unsigned int i=0; i < aux_design_.GetSize(); i++)
  {
    double v = space_in[offset + i] * scaling_;
    if(!new_design && v != aux_design_[i].GetDesign())
      new_design = true;

    aux_design_[i].SetDesign(v);
    LOG_DBG(aux_des) << "ReadDesignFromExtern: shapeparams_[i]=" << v;
  }
  if(new_design && design_id <= old_design) 
    design_id++; // if new design and not already changed by DesignSpace
  return design_id;
}

inline unsigned int AuxDesign::AuxDesignOffset() const
{
  if(export_fe_design_)
    return DesignSpace::GetNumberOfVariables();
  // no export_fe_design_ is shape optimization of shape mapping
  if(tailing_aux_design_)
    return GetNumberOfVariables() - aux_design_.GetSize(); // GetNumberOfVariables is virtual and covers in shape map design also aux_design_
  // we are first
  return 0;
}

bool AuxDesign::CompareDesign(const double* space_in)
{
  if(export_fe_design_ && !DesignSpace::CompareDesign(space_in))
    return false;

  unsigned int offset = AuxDesignOffset();

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
  if(export_fe_design_)
    DesignSpace::WriteDesignToExtern(space_out, scale);

  double rscaling = scale ? 1.0 / scaling_ : 1.0;
  unsigned int offset = AuxDesignOffset();

  for(unsigned int i=0; i < aux_design_.GetSize(); i++)
  {
    space_out[offset + i] = aux_design_[i].GetDesign() * rscaling;
    LOG_DBG(aux_des) << "WriteDesignToExtern: out[" << i << "]=" << space_out[offset +i];
  }
  return design_id;
}

void AuxDesign::WriteGradientToExtern(StdVector<double>& out, DesignElement::ValueSpecifier vs, DesignElement::Access access, Function* f, bool scaling)
{
  LOG_DBG(aux_des) << "WGTE: ad=" << aux_design_.GetSize() << " DS:GNOV=" << DesignSpace::GetNumberOfVariables() << " owst=" << out.window.GetStart() << " owsz=" << out.window.GetSize();

  bool write_aux = true;
  if(export_fe_design_ && ( f == NULL || f->GetType() != Function::SHAPE_INF) ) // SHAPE_INF, does have a sparse gradient, but no components of it are in the designspace, only in auxspace
  {
    // the number of DesignSpace variables is complicated because of constant region.
    unsigned int data_size = DesignSpace::GetNumberOfVariables(); 

    // we call DesignSpace::WriteDenseGradientToExtern() for the ersatz material part.

    // in case the the gradient window covers DesignSpace design and aux design, we need to
    // modify the window before calling DesignSpace::Write*GradientToExtern().
    // this is e.g. not the case for sparse gradients not touching the aux design.
    StdVector<double>::Window org_window = out.window; // I like standard constructors :)

    // FIXME! window != data!!!
    assert(out.window.GetSize() <= data_size + aux_design_.GetSize());
    if(out.window.GetSize() > data_size)
      out.window.Set(out.window.GetStart(), data_size);
    else
      write_aux = false;

    if(f == NULL || f->HasDenseJacobian())
      DesignSpace::WriteDenseGradientToExtern(out, vs, access, f, scaling);
    else
      DesignSpace::WriteSparseGradientToExtern(out, vs, access, f, scaling);

    // The original window needs to be restored afterwards for assert() in BaseOptimization which
    // does not know about separation into shape and mat variables
    if(out.window.GetSize() != org_window.GetSize()) // restore original window
      out.window = org_window;
  }

  // makes use of the window within out even  if only a part of the window is used in the export_fe_design_ case
  // check if there is something to write. E.g. for FeasPP out.size is the size sparsity size, don't overwrite
  // a single designBound value with 0 from aux_design_.
  if(write_aux) 
  {
    if(f == NULL || f->HasDenseJacobian()) 
      WriteAuxGradientToExtern(out, f, scaling);
    else
      WriteSparseAuxGradientToExtern(out, f, scaling);
  }
}


void AuxDesign::WriteAuxGradientToExtern(StdVector<double>& out, Function* f, bool scale) const
{
  Condition* g = dynamic_cast<Condition*>(f);

  assert(tailing_aux_design_); // otherwise the stuff below would fail ?!
  unsigned int base = out.window.GetStart()  + out.window.GetSize() - aux_design_.GetSize();

  LOG_DBG(aux_des) << "WAGTE: g=" << (f == NULL ? "null" : f->ToString()) << " ows=" << out.window.GetStart()
                    << " HS=" << HasSlackVariable() << " base=" << base;

  assert(aux_design_.GetSize() <= out.window.GetSize());
  double s = scale ? scaling_ : 1.0;
  for(unsigned int i=0; i < aux_design_.GetSize(); i++)
  {
    if(HasSlackVariable() && g != NULL && g->HasGeneralSlackBound())
    {
      assert(i == 0 || alpha_ != NULL);
      assert(aux_design_.GetSize() <= 2);
      // g = slack         -> g - slack = 0
      // g = alpha + slack -> g - alpha - slack = 0
      // g = alpha - slack -> g - alpha + slack = 0
      if(g->IsGeneralSlackBound(Condition::SLACK_VALUE))
      {
        if(i == 0)
          out[base + i] = -1.0;             // derivative w.r.t. slack
        if(i == 1)
          out[base + i] = 0.0;              // derivative w.r.t. alpha
      }
      else if(g->IsGeneralSlackBound(Condition::ALPHA_VALUE))
      {
        assert(HasSlackVariable());
        assert(aux_design_.GetSize() == 2);
        assert(alpha_ != NULL);

        if(i == 0)
          out[base + i] = 0.0;             // derivative w.r.t. slack
        if(i == 1)
          out[base + i] = -1.0;              // derivative w.r.t. alpha
      }
      else if(g->IsGeneralSlackBound(Condition::ALPHA_PLUS_SLACK_VALUE))
      {
        out[base + i] = -1.0;               // derivative w.r.t. slack and alpha is the same
      }
      else if(g->IsGeneralSlackBound(Condition::ALPHA_MINUS_SLACK_VALUE))
      {
        if(i == 0)
          out[base + i] = +1.0;             // derivative w.r.t. slack
        if(i == 1)
          out[base + i] = -1.0;             // derivative w.r.t. alpha
      }
    }
    else
    {
      if(g != NULL)
        out[base + i] = aux_design_[i].GetPlainGradient(g) * s;
      else
        out[base + i] = aux_design_[i].SumObjectiveGradient() * s;
    }

    LOG_DBG3(aux_des) << "WAGTE: g=" << (g == NULL ? "null" : g->ToString()) << " out[" << base+i << "]=" << out[base+i]
                      << " slack case=" << (HasSlackVariable() && g != NULL && g->HasGeneralSlackBound());
  }
}

void AuxDesign::WriteSparseAuxGradientToExtern(StdVector<double>& out, Function* f, bool scale) const {
  assert(f != NULL);
  assert(!f->IsObjective()); // only constraints can have sparse Jacobians

  Condition* g = dynamic_cast<Condition*>(f);
  assert(!HasSlackVariable() && !g->HasGeneralSlackBound());

  StdVector<unsigned int>& sparsity = g->GetSparsityPattern();

  unsigned int nonaux_size = DesignSpace::GetNumberOfVariables(); // do not take the virtual call

  assert(out.window.GetSize() == sparsity.GetSize());
  unsigned int base = out.window.GetStart();
  double s = scale ? scaling_ : 1.0;
  for(unsigned int i = 0; i < sparsity.GetSize(); i++)
  {
    assert(out.InWindow(base + i));
    out[base + i] = aux_design_[sparsity[i]-nonaux_size].GetPlainGradient(g) * s;
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
  if(export_fe_design_)
    DesignSpace::WriteBoundsToExtern(x_l, x_u);

  unsigned int offset = AuxDesignOffset();

  for(unsigned int i=0; i < aux_design_.GetSize(); i++)
  {
    x_l[offset + i] = aux_design_[i].GetLowerBound() / scaling_;
    x_u[offset + i] = aux_design_[i].GetUpperBound() / scaling_;
    LOG_DBG3(aux_des) << "WBTE: l[" << (offset + i) << "]=" << x_l[offset + i] << " u[" << (offset + i) << "]=" << x_u[offset + i];
  }
}

unsigned int AuxDesign::GetNumberOfVariables() const
{
  if(export_fe_design_){
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

BaseDesignElement* AuxDesign::GetAuxDesignElement(unsigned int idx){
 return(&aux_design_[idx]);
}

inline BaseDesignElement* AuxDesign::GetDesignElement(unsigned int idx)
{
  // FIXME: data.GetSize() != DesignSpace::GetNumberOfVariables()
  if(export_fe_design_ && idx < data.GetSize())
    return DesignSpace::GetDesignElement(idx);
  else
    return &aux_design_[idx - AuxDesignOffset()];
}


} // end of namespace
