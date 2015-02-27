#include <map>

#include "Optimization/StateSolution.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/Logging/log.hpp"
#include "Domain/Domain.hh"
#include "Domain/Mesh/Grid.hh"
#include "Driver/Assemble.hh"
#include "OLAS/algsys/AlgebraicSys.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/ErsatzMaterial.hh"
#include "Optimization/Excitation.hh"

using namespace std;

namespace CoupledField {

DECLARE_LOG(statesol)
DEFINE_LOG(statesol, "stateSolution")

StateSolutions::StateSolutions()
{
  this->em_ = NULL;
  forward_data_ = NULL;
  isForward = false;
}

StateSolutions::~StateSolutions()
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

void StateSolutions::Init(ErsatzMaterial* em)
{
  this->em_ = em;
}

void StateSolutions::Init(Function* f)
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


StateSolutions::Unit::Unit(ErsatzMaterial* em)
{
  // we have to delete the old data before overwriting with new stuff!
  for(unsigned int i = 0, s = data.GetSize(); i < s; ++i)
      delete data[i];
  
  data.Clear();
  data.Resize(em->me->excitations.GetSize());
  
  for(unsigned int i = 0, s = data.GetSize(); i < s; ++i)
    data[i] = new StateSolution(em);
}

StateSolutions::Unit::~Unit()
{
  for(unsigned int i = 0, s = data.GetSize(); i < s; ++i)
  {
    delete data[i];
    data[i] = NULL;
  }
}

StateSolution* StateSolutions::Get(Excitation& excitation, Function* f, unsigned int timestep, const DERIVType derivative)
{
  StateSolution* sol = Get(excitation.index, f, timestep, derivative);
  LOG_DBG2(statesol) << "S:G: e=" << excitation.index << " f=" << (f != NULL ? f->type.ToString(f->GetType()) : "NULL") << " ts=" << timestep << " d=" << derivative
                  << " -> rhs=" << sol->ToString();
  return sol;
}

StateSolution* StateSolutions::Get(int excitation_index, Function* f, unsigned int timestep, const DERIVType derivative)
{
  LOG_DBG2(statesol) << "S:G: ei=" << excitation_index << " f=" << (f != NULL ? f->type.ToString(f->GetType()) : "NULL") << " ts=" << timestep << " d=" << derivative << " iF=" << isForward;

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

StdVector<Function*> StateSolutions::GetFunctions() const
{
  StdVector<Function*> result;

  for(map<Function*, StdVector<StdVector<Unit*> > >::const_iterator it = data_.begin(); it != data_.end(); ++it)
  {
    // LOG_DBG2(statesol) << "GetFunctions(): f=" << (it->first == NULL ? "NULL" : it->first->ToString());
    if(it->first != NULL)
      result.Push_back(it->first);
  }

  return result;
}

StateSolution::StateSolution(ErsatzMaterial* em)
{
  this->em_ = em;
  this->raw = NULL;
  this->rhs = NULL;
  this->select = NULL;
}

StateSolution::~StateSolution()
{
  delete raw;
  delete rhs;
  delete select;

  std::map<Optimization::Application, StdVector<SingleVector* > >::iterator elem_iter;
  for(elem_iter = elem.begin(); elem_iter != elem.end(); elem_iter++)
  {
    StdVector<SingleVector* >& data = elem_iter->second;
    for(unsigned int i = 0; i < data.GetSize(); i++)
    {
      delete data[i];
    }
  }
}

SolutionType StateSolution::GetSolutionType(SinglePDE* pde, Optimization::Application app)
{
  switch(app)
  {
  case Optimization::NO_APP: // up to now
    assert(pde != NULL);
    return pde->GetNativeSolutionType();
  case Optimization::MECH:
    return MECH_DISPLACEMENT;
  case Optimization::ELEC:
    return ELEC_POTENTIAL;
  case Optimization::HEAT:
    return HEAT_TEMPERATURE;
  case Optimization::ACOUSTIC:
    return ACOU_POTENTIAL;
  case Optimization::LAPLACE:
    assert(false);
    break;
    // app = MECH;
    //solt = MECH_DISPLACEMENT;
  case Optimization::LBM:
    return LBM_PROBABILITY_DISTRIBUTION;
  default:
    EXCEPTION("Solution type not implemented");
  }
  return NO_SOLUTION_TYPE;
}




std::string StateSolution::ToString()
{
  std::stringstream ss;
  ss <<  "raw=" << (raw != NULL ) << " rhs=" << (rhs != NULL) << " select=" << (select != NULL);
  return ss.str();
}

SingleVector* StateSolution::Read(StorageType st, SinglePDE* pde, Optimization::Application app, bool save_sol, DERIVType derivative)
{
  if (em_->complex_)
    return Read<std::complex<double> > (st, pde, app, save_sol, derivative);
  else
    return Read<double> (st, pde, app, save_sol, derivative);
}

/** Writes the solution (raw vector) back to the pde */
void StateSolution::Write(SinglePDE* pde)
{
  if (em_->complex_)
    Write<std::complex<double> >(pde);
  else
    Write<double>(pde);
}

template <class T>
void StateSolution::Write(SinglePDE* pde)
{
  // TODO make robust for LBM
  if(raw != NULL)
  {
    assert(raw->GetSize() != 0);

    SolutionType solt = GetSolutionType(pde);
    shared_ptr<BaseFeFunction> fe = pde->GetFeFunction(solt);
    *(fe->GetSingleVector()) = *raw;
  }
  else
    LOG_DBG2(statesol) << "S:W raw not written as it was not set";
}


void StateSolution::Write(SinglePDE* pde, StateSolutions& sol, Function* f, int time_step, StdVector<Excitation>& excitations)
{
  if(domain->GetDriver()->IsComplex())
    Write<complex<double> >(pde, sol, f, time_step, excitations);
  else
    Write<double>(pde, sol, f, time_step, excitations);
}

template <class T>
void StateSolution::Write(SinglePDE* pde, StateSolutions& sol, Function* f, int time_step, StdVector<Excitation>& excitations)
{
  if(excitations.GetSize() == 1)
  {
    sol.Get(0, f)->Write(pde);
  }
  else
  {
    assert(f == NULL || sol.data_.find(f) != sol.data_.end()); // if f != NULL it has to be in the map
    StateSolutions::Unit* unit = sol.data_[f][NO_DERIVTYPE][time_step];
    assert(unit->data.GetSize() == excitations.GetSize());

    Vector<T> sum(unit->data[0]->raw->GetSize());
    sum.Init();

    for(unsigned int ex = 0; ex < excitations.GetSize(); ex ++)
    {
      StateSolution* s =unit->data[ex];
      sum.Add((T) excitations[ex].normalized_weight, *(s->raw));
    }
    Write(pde, &sum);
  }
}


void StateSolution::Write(SinglePDE* pde, SingleVector* vec)
{
  LOG_DBG2(statesol) << "S:W pde=" << pde->GetName() << " vec=" << vec->ToString(1);

  // get the coefficients from the fefunction
  SingleVector* sys = pde->GetFeFunction(pde->GetNativeSolutionType())->GetSingleVector();
  assert(sys != NULL && sys->GetSize() == vec->GetSize());

  // let the assignment operator do the job
  *sys = *vec;

  LOG_DBG3(statesol) << "S:W sys -> " << sys->ToString(0);
}

template <class T>
SingleVector* StateSolution::GetVector(StorageType st, bool create)
{
  LOG_DBG2(statesol) << "S:GV st=" << st << ToString() << " create=" << create;
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

SingleVector* StateSolution::GetVector(StorageType st)
{
  // as create is false, it makes no difference of we use the double or complex variant
  return GetVector<double>(st, false);
}

Vector<double>& StateSolution::GetRealVector(StorageType st)
{
  // we know what we want - hence we can create the result on the fly
  return dynamic_cast<Vector<double>& >(*GetVector<double>(st, true));
}

Vector<complex<double> >& StateSolution::GetComplexVector(StorageType st)
{
  return dynamic_cast<Vector<complex<double> >& >(*GetVector<complex<double> >(st, true));
}

template <class T>
SingleVector* StateSolution::Read(StorageType st, SinglePDE* pde, Optimization::Application app, bool save_sol, DERIVType derivative)
{
  assert(derivative == NO_DERIVTYPE); // would change solt!

  SolutionType solt = GetSolutionType(pde, app);

  if(app == Optimization::LAPLACE)
  {
    assert(false); // FIXME
    app = Optimization::MECH;
    solt = MECH_DISPLACEMENT;
  }

  shared_ptr<BaseFeFunction> fe = pde->GetFeFunction(solt);

  switch(st)
  {
    case GRIDELEM_VECTORS:
    {
      Grid* grid = domain->GetGrid();

      StdVector<SingleVector*>& elem_vec = gridelem[app];
      int n = grid->GetNumElems();
      if(elem_vec.GetSize() == 0)
      {
        elem_vec.Resize(n);
        for(int ve = 0; ve < n; ve++)
          elem_vec[ve] = new Vector<T>;
      }

      ElemList elemList(grid);
      for(int e = 0; e < n; e++)
      {
        elemList.SetElement(grid->GetElem(e+1)); // GetElem is 1-based
        assert(false);
        // FIXME const EntityIterator& it = elemList.GetIterator();
        // FIXME pde->GetAnyDerivSolVecOfElement((Vector<T>&) *elem_vec[e], it, resinfo, derivative);
      }
      return NULL;
    }
    case ELEMENT_VECTORS:
    {
      if(save_sol)
      {
        // we need to copy the solution from the algebraic system to the feFunction
        LOG_DBG3(statesol) << "S:R: fe sol=" << fe->GetSingleVector()->ToString();
        Vector<T> tmpSol;
        fe->GetSystem()->GetSolutionVal(tmpSol, fe->GetFctId(), true); // set idbc
        LOG_DBG3(statesol) << "S:R: sys sol=" << tmpSol.ToString();
        dynamic_cast<Vector<T>& >(*(fe->GetSingleVector())) = tmpSol;
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

      // store the results of the standard design elements in our own structure
      for(int e = 0; e < n; e++)
      {
        DesignElement* de = &em_->design->data[e];

        fe->GetEntitySolution((Vector<T>&) *elem_vec[e], de->elem);
      }

      // the pseudo design if we have some
      for(unsigned int r = 0; r < em_->design->GetPseudoDesignRegions().GetSize(); r++)
      {
        StdVector<DesignElement>& data = em_->design->GetPseudoDesignRegions()[r];

        for(unsigned int e = 0; e < data.GetSize(); e++)
        {
          DesignElement* de = &data[e];
          assert(de->GetElementSolutionIndex() >= em_->design->GetNumberOfElements());
          fe->GetEntitySolution(*elem_vec[de->GetElementSolutionIndex()], de->elem);
          // FIXME pde->GetAnyDerivSolVecOfElement((Vector<T>&) *elem_vec[de->GetElementSolutionIndex()], it, resinfo, derivative);
        }
      }

      return NULL;
    }
    case SEL_VECTOR: // interpreted as RHS_VECTOR
    case RAW_VECTOR:
    case RHS_VECTOR:
    {
      SingleVector* vec = GetVector<T>(st, true); // create when needed
      if(st == RAW_VECTOR)
      {
        // we need to copy the solution from the algebraic system to the feFunction
        // LOG_DBG3(statesol) << "S:R: fe sol=" << fe->GetSingleVector()->ToString(); // data will be outdated
        fe->GetSystem()->GetSolutionVal(*vec, fe->GetFctId(), true); // set idbc
        assert(derivative == NO_DERIVTYPE);
        // if not, the above might work ?!
        // FIXME **tmp = pde->getTimeStepping()->GetDeriveMap()[derivative]; // assigning, data is copied
      }
      else
        fe->GetSystem()->GetRHSVal(*vec, fe->GetFctId());

      LOG_DBG3(statesol) << "S:R: st=" << st << " vec=" << vec->ToString();

      return vec;
    }
  }
  assert(false);
  return(NULL);
}

}
