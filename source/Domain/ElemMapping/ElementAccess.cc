#include "ElementAccess.hh"
#include "PDE/SinglePDE.hh"
#include "Driver/Assemble.hh"
#include "Driver/FormsContexts.hh"
#include "Forms/BiLinForms/BDBInt.hh"

using namespace CoupledField;


ElementAccess::ElementAccess(SinglePDE* pde, const std::string& bilinform, RegionIdType reg)
{
  this->pde = pde;
  assert(pde != NULL);
  this->bdb = dynamic_cast<BaseBDBInt*>(pde->GetAssemble()->GetBiLinForm(bilinform, reg, pde)->GetIntegrator());
  assert(bdb != NULL);
}


ElementAccess::ElementAccess(BiLinFormContext* context)
{
  assert(context != NULL);
  this->pde = context->GetFirstPde();
  assert(pde != NULL);
  this->bdb = dynamic_cast<BaseBDBInt*>(context->GetIntegrator());
  assert(bdb != NULL);
}

ElemShapeMap* ElementAccess::SetElem(const Elem* elem)
{
  if(elem != NULL)
  {
    assert(elem->type != Elem::ET_UNDEF);

    // esm = domain->GetGrid()->GetElemShapeMap(elem);
    // esm.reset(elem->ptrShapeMap);  // ???
    // esm = elem->ptrShapeMap;
    esm = (elem)->GetElemShapeMap(domain->GetGrid());

    if(bdb != NULL)
    {
      curr_base_fe = bdb->GetFeSpace1()->GetFe(elem->elemNum);

      // TODO could be cached
      bdb->GetFeSpace1()->GetIntegration(curr_base_fe, elem->regionId, method, order);
      bdb->GetIntScheme()->GetIntPoints(Elem::GetShapeType(elem->type), method, order, intPoints, weights );

      // not every time time needed, but for convenience
      bdb->GetFeSpace1()->GetElemEqns(elem_eqn_nr, elem);

      // copy for 0-based indices
      elem_eqn_idx.Resize(elem_eqn_nr.GetSize());
      for(unsigned int i = 0; i < elem_eqn_nr.GetSize(); i++)
        elem_eqn_idx[i] = elem_eqn_nr[i] -1;
    }
  }
  else
  {
    curr_base_fe = NULL;
    // a little bit slow but might avoid bugs which are difficult to find
    elem_eqn_idx.Init(-1);
    elem_eqn_nr.Init(-1);
  }

  // reset
  curr_ip = NULL;
  curr_elem = elem;

  return curr_elem != NULL ? esm.get() : NULL;
}

LocPointMapped& ElementAccess::SetIP(unsigned int ip)
{
  assert(curr_elem != NULL);

  // Calculate for each integration point the LocPointMapped
  lpm.Set(intPoints[ip], esm, weights[ip]);

  curr_ip = &intPoints[ip];

  return lpm;
}


std::string ElementAccess::ToString(int level) const
{
  std::stringstream ss;

  switch(level)
  {
  case 0:
  case 1:
    ss << "EA: pde=" << pde->GetName() << " bdb=" << bdb->GetName();
    break;

  case 2:
    ss << "EA: e=" << (curr_elem != NULL? curr_elem->elemNum : -1) << " << method=" << method << " order=" << order.ToString() << " ips=" << intPoints.ToString()
       << " eqn_idx=" << elem_eqn_idx.ToString();
    break;

  case 3:
    ss <<  "EA e=" << (curr_elem != NULL? curr_elem->elemNum : -1) << " ip=";
    if(curr_ip != NULL)
      ss << curr_ip->number << "/" << curr_ip->coord.ToString();
    else
      ss << "NULL";
    ss << " jacDet=" << lpm.jacDet;
    break;

  default:
    assert(false);
  }

  return ss.str();
}



