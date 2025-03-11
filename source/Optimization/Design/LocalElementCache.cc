#include "Optimization/Design/LocalElementCache.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Domain/CoefFunction/CoefFunctionOpt.hh"
#include "Forms/BiLinForms/BDBInt.hh"
#include "Domain/Mesh/Grid.hh"
#include "Driver/Assemble.hh"
#include "Driver/FormsContexts.hh"
#include "Driver/EigenFrequencyDriver.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

DEFINE_LOG(lec, "localElementCache")

namespace CoupledField {

LocalElementCache::LocalElementCache(DesignSpace* space)
{
  /** reserve sufficient space to prevent expensive copying of full vectors on resize */
  this->space_ = space;
  this->regular_ = space->IsRegular();
  data_.Reserve(100); // we might have quite some wave vectors for bloch
}


void LocalElementCache::InitOrg()
{
  // if we dont't init, we switch back to org, otherwise to true
  bool next_state = active_;
  active_ = false;

  assert(Optimization::context->pdes.size() == 1);
  Assemble* assemble = Optimization::context->pde->GetAssemble();

  LOG_DBG(lec) << "IO: #forms=" << assemble->GetBiLinForms().size();

  // read the stuff from the current context
  for(std::set<BiLinFormContext*>::iterator it = assemble->GetBiLinForms().begin();
      it != assemble->GetBiLinForms().end(); it++ )
  {
    BiLinFormContext& context = **it;
    BiLinearForm*     form    = context.GetIntegrator();

    LOG_DBG(lec) << "IO: forms=" << form->GetName();
    LOG_DBG(lec) << "IO: pde1=" << context.GetFirstPde()->GetName();
    LOG_DBG(lec) << "IO: ent1=" << context.GetFirstEntities()->GetName();
    LOG_DBG(lec) << "IO: reg1=" << context.GetFirstEntities()->GetRegion();
    RegionIdType reg = context.GetFirstEntities()->GetRegion();

    // for buckling we must not use local element caching as the
    // local element matrices depend on the current stresses
    if (Optimization::context->DoBuckling() && form->GetName() == "PreStressInt")
      continue;

    // the data is not created when it exists for previous regions
    FormData& data = GetFormData(form, ORG, DesignElement::NO_DERIVATIVE, true);

    // sets the form temporarily to ORG
    next_state = FillFormData(data, form, reg) ? true : next_state; // if false we don't change the state
  }

  active_ = next_state;
}


void LocalElementCache::InitShadow(DesignSpace::DesignRegion* dr)
{
  assert(dr != NULL);
  assert(dr->HasScndMaterial());

  MaterialClass mc = NO_CLASS;
  MaterialType  mt = NO_MATERIAL;
  Context* ctxt = Optimization::context;

  switch(ctxt->mat->GetSystem())
  {
  case OptimizationMaterial::MECH:
  {
    mc = MECHANIC;
    mt = MECH_STIFFNESS_TENSOR;

    Init(space_->ToForm(mc, mt), dr->regionId, SHADOW, DesignElement::NO_DERIVATIVE, dr->GetScndMaterial(mc, mt));

    if(ctxt->pde->GetAssemble()->HasBiLinForm("MassInt", dr->regionId))
      Init(space_->ToForm(mc, DENSITY), dr->regionId, SHADOW, DesignElement::NO_DERIVATIVE, dr->GetScndMaterial(mc, DENSITY));
    break;
  }
  case OptimizationMaterial::HEAT:
  {
    mc = THERMIC;
    mt = HEAT_CONDUCTIVITY_TENSOR;
    Init(space_->ToForm(mc, mt), dr->regionId, SHADOW, DesignElement::NO_DERIVATIVE, dr->GetScndMaterial(mc, mt));
    break;
  }
  case OptimizationMaterial::MAG:
  {
    mc = ELECTROMAGNETIC;
    mt = MAG_RELUCTIVITY_TENSOR;
    Init(space_->ToForm(mc, mt), dr->regionId, SHADOW, DesignElement::NO_DERIVATIVE, dr->GetScndMaterial(mc, mt));
    break;
  }

  case OptimizationMaterial::PIEZOCOUPLING:
  case OptimizationMaterial::ACOUSTIC:
  case OptimizationMaterial::ELEC:
    assert(false);
    // implement. You might want to switch of element cache!
    break;
  case OptimizationMaterial::LBM:
  case OptimizationMaterial::NO_SYSTEM:
    assert(false);
    // there shall not be LBM bimaterial!
    break;
  } // end switch
}

void LocalElementCache::InitMechMatDeriv(StdVector<RegionIdType>& reg)
{
  // the design we test for. Note that FMO might be only for the diagonal entries, ...
  StdVector<DesignElement::Type> des;
  des.Push_back(DesignElement::MECH_11);
  des.Push_back(DesignElement::MECH_12);
  des.Push_back(DesignElement::MECH_13);
  des.Push_back(DesignElement::MECH_22);
  des.Push_back(DesignElement::MECH_23);
  des.Push_back(DesignElement::MECH_33);

  // when we have other param mat variables we cannot use LocalElementCache as these variables
  // are the only linear ones

  for(unsigned int d = 0; d < des.GetSize(); d++)
    for(unsigned int r = 0; r < reg.GetSize(); r++)
      if(space_->GetRegion(reg[r], des[d], -1, false) != NULL) // no exception to throw
        InitMatDeriv("LinElastInt", reg[r], des[d]);

  // TODO add Piezo ParamMat stuff with other integrators

}

bool LocalElementCache::InitMatDeriv(const string& integrator, RegionIdType reg, DesignElement::Type dir)
{
  switch(dir)
  {
  case DesignElement::MECH_11:
  case DesignElement::MECH_12:
  case DesignElement::MECH_13:
  case DesignElement::MECH_22:
  case DesignElement::MECH_23:
  case DesignElement::MECH_33:
    Init(integrator, reg, DIRECTION, dir, PtrCoefFct());
    return true;
  default:
    LOG_DBG(lec) << "IMD int=" << integrator << " reg=" << reg << " dir=" << DesignElement::type.ToString(dir) << " -> not to be cached";
    return false;
  }
}

void LocalElementCache::Init(const string& integrator, RegionIdType reg, Type type, DesignElement::Type dir, PtrCoefFct coef)
{
  active_ = false;

  assert(type == SHADOW || type == DIRECTION);

  SinglePDE* pde = Optimization::context->pde;
  assert(pde != NULL);

  BiLinFormContext* c = pde->GetAssemble()->GetBiLinForm(integrator, reg, pde, pde, false); // not silent
  BiLinearForm* form = c->GetIntegrator();

  FormData* data = GetFormData(integrator, type, dir, coef);
  if(data == NULL)
    data = AppendFormData(form, type, dir, coef);
  assert(data != NULL);

  FillFormData(*data, form, reg); // sets the form temporarily to SHADOW

  active_ = true;
}

void LocalElementCache::SetMatDeriv(const string& integrator, RegionIdType reg, DesignElement::Type dir, int design_id)
{
  assert(Optimization::context->DoBuckling());
  assert(space_->FindDesign(DesignElement::STIFF1, false) >= 0);

  if(space_->GetRegion(reg, DesignElement::STIFF1, -1, false) != NULL) // no exception to throw
  {
    active_ = false;

    SinglePDE* pde = Optimization::context->pde;
    assert(pde != NULL);

    BiLinFormContext* c = pde->GetAssemble()->GetBiLinForm(integrator, reg, pde, pde, false); // not silent
    BiLinearForm* form = c->GetIntegrator();
    assert(form->GetName() == "LinElastInt" || form->GetName() == "PreStressInt");

    LOG_DBG(lec) << "LEC:CC id=" << design_id << " region=" << reg << " int=" << form->GetName() << " dir=" << DesignElement::type.ToString(dir);

    PtrCoefFct coef = PtrCoefFct();

    FormData* data = GetFormData(integrator, DIRECTION, dir, coef, design_id);
    if(data == NULL)
      data = AppendFormData(form, DIRECTION, dir, coef, design_id);
    assert(data != NULL);

    BaseBDBInt* bdb = dynamic_cast<BaseBDBInt*>(form);
    shared_ptr<CoefFunctionOpt> opt = dynamic_pointer_cast<CoefFunctionOpt>(bdb->GetCoef());
    if(opt)
      opt->SetToOptimization();

    FillFormData(*data, form, reg); // sets the form temporarily to SHADOW

    active_ = true;
  }
}

bool LocalElementCache::FillFormData(FormData& data, BiLinearForm* form, RegionIdType reg)
{
  BaseBDBInt* bdb = dynamic_cast<BaseBDBInt*>(form);

  LOG_DBG(lec) << "FFD: form=" << form->GetName() << " bdb=" << (bdb != NULL) << " reg=" << reg;

  if(bdb == NULL || reg == NO_REGION_ID)
    return false;

  shared_ptr<CoefFunctionOpt> opt = dynamic_pointer_cast<CoefFunctionOpt>(bdb->GetCoef());
  if(opt != NULL)
  {
    assert(opt->GetState() == CoefFunctionOpt::OPT);
    switch(data.type)
    {
    case ORG:
      opt->SetToOrgMaterial();
      break;
    case SHADOW:
      assert(data.shadow);
      opt->SetToShadow(data.shadow);
      break;
    case DIRECTION:
      assert(data.dir != DesignElement::NO_DERIVATIVE);
      opt->SetToMaterialDerivative(data.dir);
      break;
    case NONE:
      assert(false);
      break; // please some compilers
    }
  }

  StdVector<Elem*> elems;
  ElemList elemList(domain->GetGrid());
  // region elements
  domain->GetGrid()->GetElems(elems, reg);
  EntityIterator ei = elemList.GetIterator();

  assert(CheckFormState(form, data.type));

  if(regular_ && Optimization::context->dm == NULL)
  {
    // the region case
    elemList.SetElement(elems[0]); // region's first element
    assert(ei.GetElem()->elemNum == elems[0]->elemNum);
    data.CalcElementMatrix(form, ei, true); // use the region part
  }
  else
  {
    // set every element of the region. Subject to parallelization as in Assemble :)
    for(unsigned int i = 0; i < elems.GetSize(); i++)
    {
      elemList.SetElement(elems[i]);
      assert(ei.GetElem()->elemNum == elems[i]->elemNum);
      data.CalcElementMatrix(form, ei, false); // no region!
    }
  }

  if(opt)
    opt->SetToOptimization();

  LOG_DBG(lec) << "FFD: opt=" << (opt != NULL) << " regular=" << regular_ << " elems=" << elems.GetSize();

  return true;
}


inline LocalElementCache::FormData& LocalElementCache::GetFormData(const BiLinearForm* form, Type type, DesignElement::Type dir, bool create)
{
  FormData* fd = GetFormData(form->GetName(), type, dir);

  if(fd == NULL)
  {
    if(!create) {
      EXCEPTION("cannot find cached element matrices for form " << form->GetName() << " of type " << type);
    } else
      fd = AppendFormData(form, type, dir);
  }
  assert(fd != NULL);
  assert(fd->Validate(form, type, dir));

  return *fd;
}

LocalElementCache::FormData* LocalElementCache::AppendFormData(const BiLinearForm* form, Type type, DesignElement::Type dir, PtrCoefFct shadow_coef, int design_id)
{
  switch(type)
  {
  case ORG:
    // make sure we fill first ORG, then SHADOW and DIRECTION at the end (usually then there is no SHADOW)
    assert(org_end_ == shadow_end_);
    assert(org_end_ == dir_end_);
    data_.Push_back(FormData());
    org_end_    = data_.GetSize();
    shadow_end_ = data_.GetSize();
    dir_end_    = data_.GetSize();
    break;
  case SHADOW:
    assert(shadow_end_ == dir_end_);
    assert(shadow_coef);
    data_.Push_back(FormData());
    data_.Last().shadow = shadow_coef;
    shadow_end_ = data_.GetSize();
    dir_end_    = data_.GetSize();
    break;
  case DIRECTION:
    assert(dir != DesignElement::NO_DERIVATIVE);
    data_.Push_back(FormData());
    data_.Last().dir = dir;
    data_.Last().designID = design_id;
    dir_end_    = data_.GetSize();
    break;
  default:
    assert(false);
  }

  FormData& fd = data_.Last();

  fd.Init(form, this->regular_);
  fd.idx = data_.GetSize()-1;
  fd.type = type;
  fd.contex = Optimization::context->context_idx;
  if(Optimization::context->DoBloch())
    fd.wave = Optimization::context->GetEigenFrequencyDriver()->GetCurrentWaveVector();

  assert(fd.integrator == form->GetName());
  assert(fd.isComplex == form->IsComplex());

  LOG_DBG(lec) << "ADF: " << fd.ToString();

  return &fd;
}

void LocalElementCache::ClearMatDeriv(const string& integrator, RegionIdType reg, DesignElement::Type dir, int design_id)
{
  assert(Optimization::context->DoBuckling());
  assert(space_->FindDesign(DesignElement::STIFF1, false) >= 0);

  if(space_->GetRegion(reg, DesignElement::STIFF1, -1, false) != NULL) // no exception to throw
  {
    PtrCoefFct coef = PtrCoefFct();

    FormData* data = GetFormData(integrator, DIRECTION, dir, coef, design_id);
    if(data != NULL)
    {
      data->region_cplx.Clear(false);
      data->region_real.Clear(false);
    }
  }
}

void LocalElementCache::FormData::Init(const BiLinearForm* form, bool structured)
{
  this->integrator = form->GetName();
  this->isComplex = form->IsComplex();

  Grid* grid = domain->GetGrid();
  if(structured && Optimization::context->dm == NULL)
  {
    int size = grid->regionData.GetSize();
    if(isComplex)
      region_cplx.Resize(size);
    else
      region_real.Resize(size);
  }
  else
  {
    int size = grid->GetNumElems(ALL_REGIONS) + 1; // 1-based
    if(isComplex)
      elem_cplx.Resize(size);
    else
      elem_real.Resize(size);
  }
}

inline void LocalElementCache::FormData::CalcElementMatrix(BiLinearForm* form, EntityIterator& entity, bool region)
{
  if(region)
  {
    // ensure that we don't double set values
    int reg = entity.GetElem()->regionId;
    assert(!(isComplex && region_cplx[reg].GetNumCols() != 0));
    assert(!(!isComplex && region_real[reg].GetNumCols() != 0));
    if(isComplex)
      form->CalcElementMatrix(region_cplx[reg], entity, entity);
    else
      form->CalcElementMatrix(region_real[reg], entity, entity);
  }
  else
  {
    unsigned int elem_num = entity.GetElem()->elemNum;
    assert(!(isComplex && elem_cplx[elem_num].GetNumCols() != 0));
    assert(!(!isComplex && elem_real[elem_num].GetNumCols() != 0));
    if(isComplex)
      form->CalcElementMatrix(elem_cplx[elem_num], entity, entity); // no 0-based conversion!
    else
      form->CalcElementMatrix(elem_real[elem_num], entity, entity);
  }
}

bool LocalElementCache::FormData::Validate(const BiLinearForm* form, Type type, DesignElement::Type dir)
{
  assert(this->integrator == form->GetName());
  assert(this->isComplex == form->IsComplex());
  assert(this->type == type);
  assert(this->dir == dir);
  return true;
}

double LocalElementCache::FormData::CalcMemory()
{
  int size = 0;

  for(unsigned int i = 0; i < region_real.GetSize(); i++)
    size += region_real[i].GetNumEntries() * sizeof(double);
  for(unsigned int i = 0; i < region_cplx.GetSize(); i++)
    size += region_cplx[i].GetNumEntries() * sizeof(Complex);
  for(unsigned int i = 0; i < elem_real.GetSize(); i++)
    size += elem_real[i].GetNumEntries() * sizeof(double);
  for(unsigned int i = 0; i < elem_cplx.GetSize(); i++)
    size += elem_cplx[i].GetNumEntries() * sizeof(Complex);

  return size / (1024*1024);
}

string LocalElementCache::FormData::ToString()
{
  std::stringstream ss;
  ss << "idx=" << idx << " form=" << integrator << " type=" << type
      << " cplx=" << this->isComplex << " ctxt=" << contex << " wv=" << wave.ToString();
  return ss.str();
}

template <>
const Matrix<double>& LocalElementCache::CachedElement<double>(const string& integrator, Type type, const Elem* elem, DesignElement::Type dir, PtrCoefFct shadow_coef, int design_id)
{
  FormData* data = GetFormData(integrator, type, dir, shadow_coef, design_id);
  assert(data != NULL);
  assert(data->type == type);
  assert(!data->isComplex);

  // TODO: bug if we cache not optimization data for simulation!!
  if(regular_ && Optimization::context->dm == NULL)
    return data->region_real[elem->regionId];
  else
    return data->elem_real[elem->elemNum];
}

template <>
const Matrix<complex<double> >& LocalElementCache::CachedElement<complex<double> >(const string& integrator, Type type, const Elem* elem, DesignElement::Type dir, PtrCoefFct shadow_coef, int design_id)
{
  FormData* data = GetFormData(integrator, type, dir, shadow_coef, design_id);
  assert(data != NULL);
  assert(data->type == type);
  assert(data->contex == (int) Optimization::context->context_idx);
  assert(data->isComplex);

  if(regular_ && Optimization::context->dm == NULL)
    return data->region_cplx[elem->regionId];
  else
   return data->elem_cplx[elem->elemNum];
}


template<class T>
bool LocalElementCache::CachedOrgElement(Matrix<T>& mat_out, BiLinearForm* form, const Elem* elem)
{
  if(!active_)
    return false;

  mat_out = CachedElement<T>(form->GetName(), ORG, elem); // copy constructor :)
  return true;
}


LocalElementCache::Type LocalElementCache::GetType(bool lower_bimat,  DesignElement::Type direction) const
{
  assert(!(lower_bimat && direction != DesignElement::NO_DERIVATIVE));

  if(lower_bimat)
    return SHADOW;
  if(direction != DesignElement::NO_DERIVATIVE)
    return DIRECTION;
  return ORG;
}

bool LocalElementCache::CheckFormState(BiLinearForm* form, Type type)
{
  BaseBDBInt* bdb = dynamic_cast<BaseBDBInt*>(form);
  if(bdb == NULL)
    return true;

  assert(bdb->GetCoef());

  shared_ptr<CoefFunctionOpt> opt = dynamic_pointer_cast<CoefFunctionOpt>(bdb->GetCoef());
  if(!opt)
    return true;

  switch(type)
  {
  case ORG:
    assert(opt->GetState() == CoefFunctionOpt::ORG);
    break;
  case SHADOW:
    assert(opt->GetState() == CoefFunctionOpt::SHADOW);
    break;
  case DIRECTION:
    assert(opt->GetState() == CoefFunctionOpt::DIRECTION);
    break;
  default:
    assert(false);
  }
  return true;
}

void LocalElementCache::ToInfo(PtrParamNode info)
{
  info->Get("org")->SetValue(org_end_);
  info->Get("shadow")->SetValue(shadow_end_ - org_end_);
  info->Get("mat_deriv")->SetValue(dir_end_ - shadow_end_);

  double mem = 0;
  for(unsigned int i = 0; i < data_.GetSize(); i++)
    mem += data_[i].CalcMemory();

  info->Get("megabytes")->SetValue(mem);

}

template bool LocalElementCache::CachedOrgElement<double>(Matrix<double>& out, BiLinearForm* form, const Elem* elem);
template bool LocalElementCache::CachedOrgElement<Complex>(Matrix<Complex>& out,  BiLinearForm* form, const Elem* elem);


} // end of namespace

