#include "Optimization/Design/LocalElementCache.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Domain/CoefFunction/CoefFunctionOpt.hh"
#include "Forms/BiLinForms/BDBInt.hh"
#include "Domain/Mesh/Grid.hh"
#include "Driver/Assemble.hh"
#include "Driver/FormsContexts.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/Logging/log.hpp"

DECLARE_LOG(lec)
DEFINE_LOG(lec, "localElementCache")

namespace CoupledField {

LocalElementCache::LocalElementCache(DesignSpace* space)
{
  /** reserve sufficient space to prevent expensive copying of full vectors on resize */
  this->space_ = space;
  this->regular_ = space->IsRegular(true); // also check "enforce_unstrucutred"
  data_.Reserve(100); // we might have quite some wave vektors for bloch
}


void LocalElementCache::InitOrg()
{
  // if we dont't init, we switch back to org, otherwise to true
  bool next_state = active_;
  active_ = false;

  assert(Optimization::manager.context.GetSize() == 1);
  assert(Optimization::manager.context[0].pdes.size() == 1);
  Assemble* assemble = Optimization::manager.context[0].pde->GetAssemble();

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

    // the data is not created when it exists for previous regions
    FormData& data = GetFormData(form, ORG, DesignElement::NO_DERIVATIVE, true);

    // sets the form temporarily to ORG
    next_state = FillFormData(data, form, reg) ? true : next_state; // if false we don't change the state
  }

  active_ = next_state;
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

  assert(Optimization::manager.context.GetSize() == 1);
  assert(Optimization::manager.context[0].pdes.size() == 1); // implement!
  SinglePDE* pde = Optimization::manager.context[0].pde;
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

  if(regular_)
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

LocalElementCache::FormData* LocalElementCache::AppendFormData(const BiLinearForm* form, Type type, DesignElement::Type dir, PtrCoefFct shadow_coef)
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
    dir_end_    = data_.GetSize();
    break;
  default:
    assert(false);
  }

  data_.Last().Init(form, this->regular_);
  data_.Last().type = type;

  assert(data_.Last().integrator == form->GetName());
  assert(data_.Last().isComplex == form->IsComplex());

  return &data_.Last();
}

void LocalElementCache::FormData::Init(const BiLinearForm* form, bool structured)
{
  this->integrator = form->GetName();
  this->isComplex = form->IsComplex();

  Grid* grid = domain->GetGrid();
  if(structured)
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
    assert(!(isComplex && region_cplx[elem_num].GetNumCols() != 0));
    assert(!(!isComplex && region_real[elem_num].GetNumCols() != 0));
    if(isComplex)
      form->CalcElementMatrix(elem_cplx[elem_num], entity, entity); // no 0-based conversion!
    else
      form->CalcElementMatrix(region_real[elem_num], entity, entity);
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

template <>
const Matrix<double>& LocalElementCache::CachedElement<double>(const string& integrator, Type type, const Elem* elem, DesignElement::Type dir, PtrCoefFct shadow_coef)
{
  FormData* data = GetFormData(integrator, type, dir, shadow_coef);
  assert(data != NULL);
  assert(data->type == type);
  assert(!data->isComplex);

  if(regular_)
    return data->region_real[elem->regionId];
  else
   return data->elem_real[elem->elemNum];
}

template <>
const Matrix<complex<double> >& LocalElementCache::CachedElement<complex<double> >(const string& integrator, Type type, const Elem* elem, DesignElement::Type dir, PtrCoefFct shadow_coef)
{
  FormData* data = GetFormData(integrator, type, dir, shadow_coef);
  assert(data != NULL);
  assert(data->type == type);
  assert(data->isComplex);

  if(regular_)
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


#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
template bool LocalElementCache::CachedOrgElement<double>(Matrix<double>& out, BiLinearForm* form, const Elem* elem);
template bool LocalElementCache::CachedOrgElement<Complex>(Matrix<Complex>& out,  BiLinearForm* form, const Elem* elem);
#endif


} // end of namespace

