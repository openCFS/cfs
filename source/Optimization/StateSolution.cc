#include <map>

#include "Optimization/StateSolution.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Domain/Domain.hh"
#include "Domain/Mesh/Grid.hh"
#include "Driver/Assemble.hh"
#include "OLAS/algsys/AlgebraicSys.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/ErsatzMaterial.hh"
#include "Optimization/Excitation.hh"

using namespace std;

using namespace CoupledField;

DEFINE_LOG(statesol, "stateSolution")

StateContainer::StateContainer()
{
}

StdVector<double> StateContainer::CollectEigenfrequencies(Excitation& ex)
{
  assert(Optimization::manager.context[ex.sequence-1].IsEigenvalue()); // 1-based

  const StdVector<StateSolution*> states = Search(&ex, ex.sequence);

  StdVector<double> efs;
  efs.Reserve(states.GetSize()); // shall be the final size ?!

  for(unsigned int i = 0; i < states.GetSize(); i++)
    if(states[i]->eigenfreq != -1.0) // for degenerate structures we might have negative eigenvalues!!
      efs.Push_back(states[i]->eigenfreq);
    else
      throw Exception("degenerate eigenfrequency for ex=" + ex.GetFullLabel()); // why shall this happen??

  assert(efs.GetSize() == states.GetSize()); // what shall be the scenario?! If one exists just delete the assert

  // the current context is not necessarily the eigenvalue context and we typically don't have constraints for each state
  assert(Optimization::manager.context[ex.sequence-1].num_eigenmodes >= efs.GetSize());

  return efs;
}

Matrix<double> StateContainer::CollectBlochEigenfrequencies(Context* ctxt)
{
  assert(ctxt->DoBloch());
  assert(ctxt->num_eigenmodes >= 1 && ctxt->num_bloch_wave_vectors >= 1);

  Matrix<double> mat(ctxt->num_eigenmodes, ctxt->num_bloch_wave_vectors);
  mat.Init();

  for(unsigned int e = 0; e < ctxt->excitations.GetSize(); e++)
  {
    const Excitation* ex = ctxt->excitations[e];
    const StdVector<StateSolution*> states = Search(ex, ex->sequence);
    for(unsigned int i = 0; i < states.GetSize(); i++)
    {
      assert(states[i]->eigenfreq != -1); // why is this handled in CollectEigenfrequencies()?
      mat[i][e] = states[i]->eigenfreq;
    }
  }
  return mat;
}

StateSolution* StateContainer::Get(const Excitation* ex, const Function* f, int timestep_mode_local, TimeDeriv derivative)
{
  assert(ex != NULL);
  return Get(*ex, f, timestep_mode_local, derivative);
}

StateSolution* StateContainer::Get(const Excitation& ex, const Function* f, int timestep_mode_local, TimeDeriv derivative)
{
  LOG_DBG2(statesol) << "SC:G: ex=" << ex.index << " f=" << (f != NULL ? f->type.ToString(f->GetType()) : "NULL")
                     << " ts=" << timestep_mode_local << " d=" << derivative << " size=" << data_.GetSize();

  assert(timestep_mode_local == -1 || (f != NULL && f->IsLocal()) || Optimization::manager.any().eigenvalue); // add transient

  if(data_.GetCapacity() == 0)
  {
    ErsatzMaterial* em = dynamic_cast<ErsatzMaterial*>(domain->GetOptimization());
    assert(em != NULL);

    unsigned int size = 2; // more does not harm
    for(auto& cond : em->constraints.all)
      if(cond->IsAdjointBased())
        size += cond->IsLocal() ? cond->GetLocal()->virtual_elem_map.GetSize() : 1;

    size *= em->GetMultipleExcitation()->excitations.GetSize();

    for(unsigned int i = 0; i < em->manager.context.GetSize(); i++)
      size *= std::max(em->manager.context[i].num_eigenmodes, (unsigned int) 1);

    data_.Reserve(size);

    LOG_DBG(statesol) << "SC:G reserved " << size;
  }

  StdVector<StateSolution*> found = Search(&ex, ex.sequence, f, timestep_mode_local, derivative);
  if(found.GetSize() > 1)
    EXCEPTION("request non-unique state from StateContainer: ex=" << ex.index << " func=" << (f == NULL ? "NULL" : f->ToString()) << " timestep_mode_local=" << timestep_mode_local << " deriv=" << derivative);

  if(found.GetSize() == 1)
    return found[0];

  // we need to add one
  if(data_.GetSize() >= data_.GetCapacity()) // We may not re-map the data as we might have returned pointers to the content before
    EXCEPTION("StateContainer capacity " << data_.GetCapacity() << " turned out to be too small");

  data_.Resize(data_.GetSize() + 1);
  StateSolution& state = data_.Last();

  state.excitation = &ex;
  state.func = f;
  state.timestep_mode_local = timestep_mode_local;
  state.derivative = derivative;

  LOG_DBG2(statesol) << "SC:G: add state ex=" << ex.index << " f=" << (f == NULL ? "NULL" : f->ToString()) << " tsm=" << timestep_mode_local << " der=" << derivative << " -> " << data_.GetSize();

  return &state;
}


StdVector<StateSolution*> StateContainer::Search(const Excitation* ex, int seq, const Function* f, int timestep_mode_local, const TimeDeriv derivative)
{
  StdVector<StateSolution*> found;

  assert(seq == -1 || (seq >= 1 && seq <= (int) Optimization::manager.context.GetSize())); // 1-based!
  assert(ex == NULL || (seq == -1 || ex->sequence == seq)); // require consistency
  assert(!(ex == NULL && seq == -1)); // one needs to be given

  LOG_DBG2(statesol) << "SC:S search for ex=" << (ex == NULL ? -1 : ex->index) << " seq=" << seq << " f=" << (f == NULL ? "NULL" : f->ToString())
                     << " ts_m=" << timestep_mode_local << " der=" << derivative << " size=" << data_.GetSize();

  for(unsigned int i = 0, n = data_.GetSize(); i < n; i++)
  {
    StateSolution& s = data_[i];
    bool ex_ok = ex == NULL || s.excitation == ex; // if no excitation is given we don't exclude
    bool seq_ok = seq == -1 ? true : s.excitation->sequence == seq;
    bool tsm_ok = timestep_mode_local == -1 ? true : s.timestep_mode_local == timestep_mode_local;

    if(ex_ok && seq_ok && s.func == f && tsm_ok && s.derivative == derivative)
      found.Push_back(&s);

    LOG_DBG2(statesol) << "SC:S cand i=" << i << " ex=" << (s.excitation == NULL ? -1 : s.excitation->index) << " seq="
                       << (s.excitation == NULL ? -1 : s.excitation->sequence) << " f=" << (s.func == NULL ? "NULL" : s.func->ToString())
                       << " ts_m=" << s.timestep_mode_local << " der=" << s.derivative << " set=" << s.ContainsState() << " -> " << found.GetSize();
  }

  return found;
}

void StateContainer::WriteAverage(SinglePDE* pde, int sequence, const Function* f)
{
  assert(sequence >= 1 && sequence <= (int) Optimization::manager.context.GetSize()); // 1-based

  if(Optimization::manager.context[sequence-1].IsComplex())
    WriteAverage<complex<double> >(pde, sequence, f);
  else
    WriteAverage<double>(pde, sequence, f);
}



template <class T>
void StateContainer::WriteAverage(SinglePDE* pde, int sequence, const Function* f)
{
  StdVector<StateSolution*> states = Search(NULL, sequence, f);

  if(states.GetSize() == 0)
    return; // shall be only LBM?

  if(states.GetSize() == 1)
    states[0]->Write(pde);
  else
  {
    Vector<T>* raw = dynamic_cast<Vector<T>* >(states[0]->GetVector(StateSolution::RAW_VECTOR));

    Vector<T> sum(raw->GetSize());

    for(unsigned int i = 0; i < states.GetSize(); i++)
    {
      raw = dynamic_cast<Vector<T>* >(states[i]->GetVector(StateSolution::RAW_VECTOR));
      assert(raw->GetSize() == sum.GetSize());
      assert(states[i]->excitation != NULL);
      sum.Add((T) states[i]->excitation->normalized_weight, *raw);
    }
    StateSolution::Write(pde, &sum);
  }
}


StdVector<const Function*> StateContainer::GetFunctions() const
{
  StdVector<const Function*> result;

  for(unsigned int i = 0; i < data_.GetSize(); i++)
    if(data_[i].func != NULL && !result.Contains(data_[i].func ))
      result.Push_back(data_[i].func);

  return result;
}

void StateContainer::Dump()
{
  for(auto& state : data_)
    std::cout << state.ToString() << std::endl;
}



StateSolution::StateSolution()
{
  this->raw = NULL;
  this->rhs = NULL;
  this->select = NULL;
  this->eigenfreq = -1.0;
  this->set_ = false;

  // by this we identify the solution uniquely
  this->func = NULL;
  this->derivative = NO_DERIVTYPE;
  this->timestep_mode_local = -1; // only set different for eigenmodes or transient
  this->excitation = NULL;
}

StateSolution::~StateSolution()
{
  delete raw;
  delete rhs;
  delete select;

  std::map<App::Type, StdVector<SingleVector* > >::iterator elem_iter;
  for(elem_iter = elem.begin(); elem_iter != elem.end(); elem_iter++)
  {
    StdVector<SingleVector* >& data = elem_iter->second;
    for(unsigned int i = 0; i < data.GetSize(); i++)
    {
      delete data[i];
    }
  }
}

SolutionType StateSolution::GetSolutionType(SinglePDE* pde, App::Type app)
{
  // what is the physical solution type we need to request from cfs to store as state solution?
  switch(app)
  {
  case App::NO_APP: // up to now
    assert(pde != NULL);
    return pde->GetNativeSolutionType();
  case App::MECH:
  case App::BUCKLING:
    return MECH_DISPLACEMENT;
  case App::ELEC:
    return ELEC_POTENTIAL;
  case App::HEAT:
    return HEAT_TEMPERATURE;
  case App::MAG:
	return MAG_POTENTIAL;
  case App::ACOUSTIC:
    return ACOU_PRESSURE;
  case App::LAPLACE:
    assert(false);
    break;
  case App::LBM:
    return LBM_PROBABILITY_DISTRIBUTION;
  default:
    EXCEPTION("Solution type not implemented");
  }
  return NO_SOLUTION_TYPE;
}




std::string StateSolution::ToString()
{
  std::stringstream ss;
  ss << " excite=" << (excitation != NULL ? excitation->GetFullLabel() : "NULL");
  ss << " func=" << (func != NULL ? func->ToString() : "NULL");
  ss << " raw=" << (raw != NULL ? (int) raw->GetSize() : -1 ) << " " << (raw != NULL ? BaseMatrix::entryType.ToString(raw->GetEntryType()) : "");
  ss << " rhs=" << (rhs != NULL ? (int) rhs->GetSize() : -1 ) << " " << (rhs != NULL ? BaseMatrix::entryType.ToString(rhs->GetEntryType()) : "");
  ss << " select=" << (select != NULL ? (int) select->GetSize() : -1) << " " << (select != NULL ? BaseMatrix::entryType.ToString(select->GetEntryType()) : "");
  return ss.str();
}

SingleVector* StateSolution::Read(StorageType st, SinglePDE* pde, App::Type app, bool save_sol, TimeDeriv derivative)
{
  bool buckling_adjoint = Optimization::context->DoBuckling() && app == App::MECH && pde->GetAnalysisType() == BasePDE::STATIC;

  if (Optimization::context->IsComplex() && !buckling_adjoint)
    return Read<std::complex<double> > (st, pde, app, save_sol, derivative);
  else
    return Read<double>(st, pde, app, save_sol, derivative);
}

template <class T>
SingleVector* StateSolution::Read(StorageType st, SinglePDE* pde, App::Type app, bool save_sol, TimeDeriv derivative)
{
  assert(derivative == NO_DERIVTYPE); // would change solt!

  // if called via StateContainer::Get() many parameters are set - this is not the case when called via EigenvalueState::Init()

  // these properties were set in StateContainer::Get() on the creation of the object
  assert(this->derivative == derivative);
  assert(this->excitation == NULL || this->excitation == Optimization::context->GetExcitation());

  this->excitation = Optimization::context->GetExcitation(); // for EigenvalueState::Init()
  set_ = true;  // might be already set in practice

  SolutionType solt = GetSolutionType(pde, app);

  LOG_DBG2(statesol) << "SS:R st=" << st << " org='" << ToString() << "'";

  if(app == App::LAPLACE)
  {
    assert(false); // FIXME
    app = App::MECH;
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
      DesignSpace* design = domain->GetOptimization()->GetDesign();

      if(save_sol)
      {
        // we need to copy the solution from the algebraic system to the feFunction
        LOG_DBG3(statesol) << "SS:R1: fe sol=" << fe->GetSingleVector()->ToString();
        Vector<T> tmpSol;
        if(app == App::MAG)
          fe->GetSystem()->GetSolutionVal(tmpSol, fe->GetFctId(), false);
        else
          fe->GetSystem()->GetSolutionVal(tmpSol, fe->GetFctId(), true); // set idbc

        LOG_DBG3(statesol) << "SS:R2: sys sol=" << tmpSol.ToString();
        dynamic_cast<Vector<T>& >(*(fe->GetSingleVector())) = tmpSol;
        LOG_DBG3(statesol) << "SS:R3: fe sol-nachher=" << fe->GetSingleVector()->ToString();
      }

      // we save the element vectors in elem_vec. Might be empty the first call
      StdVector<SingleVector*>& elem_vec = elem[app];

      int n = design->GetNumberOfElements(); // the standard design elements
      int pn = design->CalcPseudoDesignElements(); // optional (multiple) pseudo design elements

      // check for first call
      if(elem_vec.GetSize() == 0)
      {
        elem_vec.Resize(n + pn); // save resize as the size was checked for zero
        for(int ve = 0; ve < (n + pn); ve++)
          elem_vec[ve] = new Vector<T>; // FIXME: where do we delete
      }

      // store the results of the standard design elements in our own structure
      for(int e = 0; e < n; e++)
      {
        DesignElement* de = &design->data[e];

        fe->GetEntitySolution((Vector<T>&) *elem_vec[e], de->elem);
      }

      // the pseudo design if we have some
      for(unsigned int r = 0; r < design->GetPseudoDesignRegions().GetSize(); r++)
      {
        StdVector<DesignElement>& data = design->GetPseudoDesignRegions()[r];

        for(unsigned int e = 0; e < data.GetSize(); e++)
        {
          DesignElement* de = &data[e];
          assert(de->GetElementSolutionIndex() >= design->GetNumberOfElements());
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
        if(app == App::MAG)
          //fe->GetSystem()->GetSolutionVal(*vec, fe->GetFctId(), false); // set idbc
          *vec = dynamic_cast<Vector<T>& >(*(fe->GetSingleVector()));
        else
          fe->GetSystem()->GetSolutionVal(*vec, fe->GetFctId(), true); // set idbc

        assert(derivative == NO_DERIVTYPE);
        // if not, the above might work ?!
        // FIXME **tmp = pde->getTimeStepping()->GetDeriveMap()[derivative]; // assigning, data is copied
      }
      else
        fe->GetSystem()->GetRHSVal(*vec, fe->GetFctId());


      LOG_DBG3(statesol) << "SS:R: st=" << st << " vec=" << vec->ToString() << " size= " << vec->GetSize();

      return vec;
    }
  }
  assert(false);
  return(NULL);
}


/** Writes the solution (raw vector) back to the pde */
void StateSolution::Write(SinglePDE* pde)
{
  if (Optimization::context->IsComplex())
    Write<std::complex<double> >(pde);
  else
    Write<double>(pde);
}

template <class T>
void StateSolution::Write(SinglePDE* pde)
{
  assert(set_);

  LOG_DBG(statesol) << "SS:W pde=" << pde->GetName() << " ss=" << ToString();

  // TODO make robust for App::LBM
  if(raw != NULL)
  {
    assert(raw->GetSize() != 0);

    SolutionType solt = GetSolutionType(pde);
    shared_ptr<BaseFeFunction> fe = pde->GetFeFunction(solt);
    LOG_DBG3(statesol) << "SS:W raw = " << raw->ToString();
    LOG_DBG3(statesol) << "SS:W fe =  " << fe->GetSingleVector()->ToString();
    *(fe->GetSingleVector()) = *raw; // out of two pointers we make references and then use the assignment operator
    LOG_DBG3(statesol) << "SS:W fe =  " << fe->GetSingleVector()->ToString();
  }
  else
    LOG_DBG2(statesol) << "SS:W raw not written as it was not set";
}


void StateSolution::Write(SinglePDE* pde, SingleVector* vec)
{
  LOG_DBG2(statesol) << "SS:W pde=" << pde->GetName() << " vec=" << vec->ToString();

  // get the coefficients from the fefunction
  SingleVector* sys = pde->GetFeFunction(pde->GetNativeSolutionType())->GetSingleVector();
  assert(sys != NULL && sys->GetSize() == vec->GetSize());

  // let the assignment operator do the job
  *sys = *vec;

  LOG_DBG3(statesol) << "SS:W sys -> " << sys->ToString();
}




template <class T>
SingleVector* StateSolution::GetVector(StorageType st, bool create)
{
  LOG_DBG2(statesol) << "SS:GV st=" << st << " org='" << ToString() << "' create=" << create;
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
