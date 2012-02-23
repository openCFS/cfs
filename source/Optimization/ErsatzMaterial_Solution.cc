#include <map>

#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"
#include "Domain/domain.hh"
#include "Domain/grid.hh"
#include "Driver/assemble.hh"
#include "OLAS/algsys/basesystem.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/ErsatzMaterial.hh"
#include "Optimization/Excitation.hh"

using namespace std;

namespace CoupledField {

DECLARE_LOG(emsol)
DEFINE_LOG(emsol, "ersatzMaterialSolutions")

ErsatzMaterial::Solutions::Solutions()
{
  this->em_ = NULL;
  forward_data_ = NULL;
  isForward = false;
}

ErsatzMaterial::Solutions::~Solutions()
{
  StdVector<Function*> funcs = GetFunctions();
  for(unsigned int fi = 0; fi < funcs.GetSize(); fi++)
  {
    StdVector<StdVector<Unit*> >& list = data_[funcs[fi]];
    for(unsigned int deriv = 0; deriv < list.GetSize(); deriv++){
      for(unsigned int ts = 0; ts < list[deriv].GetSize(); ++ts)
        delete list[deriv][ts];
    }
  }
}

void ErsatzMaterial::Solutions::Init(ErsatzMaterial* em)
{
  this->em_ = em;
}

void ErsatzMaterial::Solutions::Init(Function* f)
{
  assert(em_ != NULL);

  StdVector<StdVector<Unit*> >& list = data_[f];
  
  unsigned int nts = domain->GetDriver()->GetNumSteps();

  // we have to delete the old data before overwriting with new stuff!
  for(unsigned int deriv = 0; deriv < list.GetSize(); deriv++){
    for(unsigned int ts = 0; ts < list[deriv].GetSize(); ++ts)
      delete list[deriv][ts];
    list[deriv].Resize(0); // this is not costly, as reallocation is only done if needed
  }
  
  unsigned int nderiv = SECOND_DERIV;
  if(f != NULL || nts == 0){ // in case of non-transient or adjoint, do not allocate derivative space
    nderiv = NO_DERIVTYPE;
  }
  list.Resize(nderiv+1);
  for(unsigned int deriv = NO_DERIVTYPE; deriv <= nderiv; deriv++){
    list[deriv].Resize(nts);
    for(unsigned int ts = 0; ts < nts; ++ts)
      list[deriv][ts] = new Unit(em_);
  }
}


ErsatzMaterial::Solutions::Unit::Unit(ErsatzMaterial* em)
{
  // we have to delete the old data before overwriting with new stuff!
  for(unsigned int i = 0, s = data.GetSize(); i < s; ++i)
      delete data[i];
  
  data.Clear();
  data.Resize(em->me->excitations.GetSize());
  
  for(unsigned int i = 0, s = data.GetSize(); i < s; ++i)
    data[i] = new Solution(em);
}

ErsatzMaterial::Solutions::Unit::~Unit()
{
  for(unsigned int i = 0, s = data.GetSize(); i < s; ++i)
  {
    delete data[i];
    data[i] = NULL;
  }
}

ErsatzMaterial::Solution* ErsatzMaterial::Solutions::Get(Excitation& excitation, Function* f, unsigned int timestep, const DERIVType derivative)
{
  return Get(excitation.index, f, timestep, derivative);
}

ErsatzMaterial::Solution* ErsatzMaterial::Solutions::Get(int excitation_index, Function* f, unsigned int timestep, const DERIVType derivative)
{
  if(isForward){ // if this is true, f is ignored and forward_data_ is used to avoid one map access
    if(forward_data_ == NULL){
      if(data_.find(NULL) == data_.end()){
        Init((Function*)NULL);
      }
      forward_data_ = &data_[NULL];
    }
    return((*forward_data_)[derivative][timestep]->data[excitation_index]);
  }else{
    // do we have to init first?
    if(data_.find(f) == data_.end())
      Init(f);
    assert(data_.find(f) != data_.end());
    return data_[f][derivative][timestep]->data[excitation_index];
  }
}

StdVector<Function*> ErsatzMaterial::Solutions::GetFunctions() const
{
  StdVector<Function*> result;

  for(map<Function*, StdVector<StdVector<Unit*> > >::const_iterator it = data_.begin(); it != data_.end(); ++it)
  {
    // LOG_DBG2(emsol) << "GetFunctions(): f=" << (it->first == NULL ? "NULL" : it->first->ToString());
    if(it->first != NULL)
      result.Push_back(it->first);
  }

  return result;
}

ErsatzMaterial::Solution::Solution(ErsatzMaterial* em)
{
  this->em_ = em;
  this->raw = NULL;
  this->rhs = NULL;
  this->select = NULL;
}

ErsatzMaterial::Solution::~Solution()
{
  delete raw;
  delete rhs;
  delete select;

  std::map<Application, StdVector<SingleVector* > >::iterator elem_iter;
  for(elem_iter = elem.begin(); elem_iter != elem.end(); elem_iter++)
  {
    StdVector<SingleVector* >& data = elem_iter->second;
    for(unsigned int i = 0; i < data.GetSize(); i++)
    {
      delete data[i];
    }
  }
}

SingleVector* ErsatzMaterial::Solution::Read(StorageType st, StdPDE* pde, Application app, bool save_sol, DERIVType derivative)
{
  if (em_->harmonic)
    return Read<std::complex<double> > (st, pde, app, save_sol, derivative);
  else
    return Read<double> (st, pde, app, save_sol, derivative);
}

/** Writes the solution (raw vector) back to the pde */
void ErsatzMaterial::Solution::Write(StdPDE* pde)
{
  if (em_->harmonic)
    Write<std::complex<double> >(pde);
  else
    Write<double>(pde);
}

template <class T>
void ErsatzMaterial::Solution::Write(StdPDE* pde)
{
  T* ptr = NULL;
  Vector<T>* v = static_cast<Vector<T>*>(raw);
  ptr = v->GetPointer();
  assert(ptr != NULL);
  assert(raw->GetSize() != 0);
  pde->SaveSolution(ptr, raw->GetSize());
}


void ErsatzMaterial::Solution::Write(StdPDE* pde, Solutions& sol, Function* f, int time_step, StdVector<Excitation>& excitations)
{
  if(BasePDE::IsComplex(domain->GetDriver()->GetAnalysisType()))
    Write<complex<double> >(pde, sol, f, time_step, excitations);
  else
    Write<double>(pde, sol, f, time_step, excitations);
}

template <class T>
void ErsatzMaterial::Solution::Write(StdPDE* pde, Solutions& sol, Function* f, int time_step, StdVector<Excitation>& excitations)
{
  if(excitations.GetSize() == 1)
  {
    sol.Get(0, f)->Write(pde);
  }
  else
  {
    assert(f == NULL || sol.data_.find(f) != sol.data_.end()); // if f != NULL it has to be in the map
    Solutions::Unit* unit = sol.data_[f][NO_DERIVTYPE][time_step];
    assert(unit->data.GetSize() == excitations.GetSize());

    Vector<T> sum(unit->data[0]->raw->GetSize());
    sum.Init();

    for(unsigned int ex = 0; ex < excitations.GetSize(); ex ++)
    {
      Solution* s =unit->data[ex];
      sum.Add((T) excitations[ex].normalized_weight, *(s->raw));
    }
    Write(pde, &sum);
  }
}


void ErsatzMaterial::Solution::Write(StdPDE* pde, SingleVector* vec)
{
  // we are static
  bool harmonic = BasePDE::IsComplex(domain->GetDriver()->GetAnalysisType());
  LOG_DBG2(emsol) << "S:W pde=" << pde->GetName() << " h=" << harmonic << " vec=" << vec->ToString();
  if(harmonic)
  {
    Vector<complex<double> >& sol = dynamic_cast<Vector<complex<double> >& >(*vec);
    pde->SaveSolution(sol.GetPointer(), sol.GetSize());
  }
  else
  {
    Vector<double>& sol = dynamic_cast<Vector<double>& >(*vec);
    pde->SaveSolution(sol.GetPointer(), sol.GetSize());
  }
}

template <class T>
SingleVector* ErsatzMaterial::Solution::GetVector(StorageType st, bool create)
{
  switch(st)
  {
  case RAW_VECTOR:
    assert(raw != NULL || create);
    if(raw == NULL && create) // this will crash for !create in release
      raw = new Vector<T>();
    return raw;

  case RHS_VECTOR:
    assert(rhs != NULL || create);
    if(rhs == NULL && create)
      rhs = new Vector<T>();
    return rhs;

  case SEL_VECTOR:
    assert(select != NULL || create);
    if(select == NULL && create)
      select = new Vector<T>();
    return select;

  default:
    assert(false);
  }

  EXCEPTION("false");
}

SingleVector* ErsatzMaterial::Solution::GetVector(StorageType st)
{
  // as create is false, it makes no difference of we use the double or complex variant
  return GetVector<double>(st, false);
}

Vector<double>& ErsatzMaterial::Solution::GetRealVector(StorageType st)
{
  // we know what we want - hence we can create the result on the fly
  return dynamic_cast<Vector<double>& >(*GetVector<double>(st, true));
}

Vector<complex<double> >& ErsatzMaterial::Solution::GetComplexVector(StorageType st)
{
  return dynamic_cast<Vector<complex<double> >& >(*GetVector<complex<double> >(st, true));
}

template <class T>
SingleVector* ErsatzMaterial::Solution::Read(StorageType st, StdPDE* pde, Application app, bool save_sol, DERIVType derivative)
{
  SolutionType solt;
  switch(app)
  {
  case NO_APP: // up to now
    solt = dynamic_cast<SinglePDE*>(pde)->GetNativeSolutionType();
    break;
  case MECH:
    solt = MECH_DISPLACEMENT;
    break;
  case ELEC:
    solt = ELEC_POTENTIAL;
    break;
  case HEAT:
    solt = HEAT_TEMPERATURE;
    break;
  case ACOUSTIC:
    solt = ACOU_POTENTIAL;
    break;
  case LAPLACE:
    app = MECH;
    solt = MECH_DISPLACEMENT;
    break;
  default:
    EXCEPTION("Solution type not implemented");
  }

  shared_ptr<ResultInfo> resinfo = pde->GetResultInfo(solt);
  Grid* grid = domain->GetGrid();

  switch(st)
    {
    case GRIDELEM_VECTORS:
    {
      StdVector<SingleVector*>& elem_vec = gridelem[app];
      int n = grid->GetNumElems();
      if(elem_vec.GetSize() == 0){
        elem_vec.Resize(n);
        for(int ve = 0; ve < n; ve++){
          elem_vec[ve] = new Vector<T>;
        }
      }
      ElemList elemList(grid);
      for(int e = 0; e < n; e++){
        elemList.SetElement(grid->GetElem(e+1)); // GetElem is 1-based
        const EntityIterator& it = elemList.GetIterator();
        pde->GetSolVecOfElement((Vector<T>&) *elem_vec[e], it, resinfo);
      }
      return NULL;
    }
    case ELEMENT_VECTORS:
    {
      if(save_sol)
      {
        // store the solution in the PDE! This is necessary to read the element vector
        // in the adjoint case!
        Vector<T> tmpSol;
        em_->assemble_->GetAlgSys()->GetSolutionVal(tmpSol);
        pde->SaveSolution(tmpSol.GetPointer(), tmpSol.GetSize() );
      }

      // we save the element vectors in elem_vec. Might be empty the first call
      StdVector<SingleVector*>& elem_vec = elem[app];

      int n = em_->design->GetNumberOfElements(); // the standard design elements
      int pn = em_->design->CalcPseudoDesignElements(); // optional (multiple) pseudo design elements

      // check for first call
      if(elem_vec.GetSize() == 0)
      {
        elem_vec.Resize(n + pn);
        for(int ve = 0; ve < (n + pn); ve++)
          elem_vec[ve] = new Vector<T>;
      }

      // create an element list to gain the iterator in the loop
      ElemList elemList(grid);

      // store the results of the standard design elements in our own structure
      for(int e = 0; e < n; e++)
      {
        DesignElement* de = &em_->design->data[e];

        elemList.SetElement(de->elem);
        const EntityIterator& it = elemList.GetIterator();

        assert((int) de->GetElementSolutionIndex() == e);
        
        pde->GetAnyDerivSolVecOfElement((Vector<T>&) *elem_vec[e], it, resinfo, derivative);
      }
      // the pseudo design if we have some
      for(unsigned int r = 0; r < em_->design->GetPseudoDesignRegions().GetSize(); r++)
      {
        StdVector<DesignElement>& data = em_->design->GetPseudoDesignRegions()[r];

        for(unsigned int e = 0; e < data.GetSize(); e++)
        {
          DesignElement* de = &data[e];

          elemList.SetElement(de->elem);
          const EntityIterator& it = elemList.GetIterator();

          assert(de->GetElementSolutionIndex() >= em_->design->GetNumberOfElements());
          pde->GetAnyDerivSolVecOfElement((Vector<T>&) *elem_vec[de->GetElementSolutionIndex()], it, resinfo, derivative);
        }
      }

      return NULL;
    }
    case SEL_VECTOR: // interpreted as RHS_VECTOR
    case RAW_VECTOR:
    case RHS_VECTOR:
    {
      // It is best to get the RHS from the algebraic system (OLAS), the PDE
      // might not check about a further adjont calculation

      // get access to solution
      BaseSystem* bs = em_->assemble_->GetAlgSys();
      Vector<T>** tmp; // tmp is ** Vec, *tmp is *Vec = raw/rhs/select, **tmp is Vec
      switch(st){
      case RAW_VECTOR:
        tmp = (Vector<T>**)&raw; break;
      case RHS_VECTOR:
        tmp = (Vector<T>**)&rhs; break;
      case SEL_VECTOR:
        tmp = (Vector<T>**)&select; break;
      default:
        assert(false);
      }
      if(*tmp == NULL){
        *tmp = new Vector<T>(0);
      }
      if(st == RAW_VECTOR) {
        if(derivative == NO_DERIVTYPE){
          bs->GetSolutionVal(**tmp);
        }else{
          **tmp = pde->getTimeStepping()->GetDeriveMap()[derivative]; // assigning, data is copied
        }
      } else {
        bs->GetRHSVal(**tmp); // RHS_VECTOR and SEL_VECTOR, they do not have derivatives, so this is ignored
      }
      return *tmp;
    }
  }
  assert(false);
  return(NULL);
}

}
