#include "Optimization/Excitation.hh"
#include "Optimization/Optimization.hh"
#include "Driver/assemble.hh"
#include "Driver/basedriver.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"
#include "PDE/StdPDE.hh"
#include "PDE/SinglePDE.hh"
#include "PDE/mechPDE.hh"
#include "PDE/elecPDE.hh" // for polarization matrix, see class TopGrad
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/Logging/cfslog.hh"

using namespace CoupledField;

Enum<MultipleExcitation::Type>  MultipleExcitation::type;

MultipleExcitation::MultipleExcitation(bool multiple, PtrParamNode pn)
{
  // set defaults
  this->stride = 1;
  this->damping = 1.0;
  this->max_gain = 1e4;
  this->multiple_excitation_ = multiple;
  this->type_ = NO_TYPE; // to be eventually overwritten soon

  // if disabled, we don't read anything
  if(!multiple || pn == NULL) return;

  this->type_ = type.Parse(pn->Get("type")->As<std::string>());

  // adjust defaults
  if(pn->Has("adjustWeights")) {
    this->stride   = pn->Get("adjustWeights")->Get("stride")->As<Integer>();
    this->max_gain = pn->Get("adjustWeights")->Get("max_gain")->As<Double>();
    this->damping  = pn->Get("adjustWeights")->Get("damping")->As<Double>();
  }
}

void MultipleExcitation::ToInfo(PtrParamNode in) const
{
  if(!multiple_excitation_) return;

  if(type_ == META_OBJECTIVE)  {
    PtrParamNode in_ = in->Get("metaObjective");
    in_->Get("type")->SetValue("best_value");
    in_->Get("damping")->SetValue(damping);
    in_->Get("stride")->SetValue(stride);
    in_->Get("max_gain")->SetValue(max_gain);
  } else
    in->Get("type")->SetValue(type.ToString(type_));
}


Excitation::Excitation()
{
  this->index = -1; // must be updated
  this->frequency = -1.0;
  this->f_link = NULL;
  this->weight = 1.0;
  this->normalized_weight = 1.0;
  linForms = new StdVector<LinearFormContext*>();
  AddLinFormsFromAssemble();
}

void Excitation::AddLinFormsFromAssemble(){
  // all already set linear Forms
  StdVector<LinearFormContext*>& assLinForms = domain->GetBasePDE()->getPDE_assemble()->GetLinForms();
  for(unsigned int i = 0; i < assLinForms.GetSize(); ++i){
    linForms->Push_back(assLinForms[i]);
  }
}

void Excitation::Apply()
{
  domain->GetOptimization()->applied_excitation = this;
  Assemble* a = domain->GetBasePDE()->getPDE_assemble();
  a->SetLinForms(linForms); // this works always, if not used just a copy of assembles linearForms
  if(! loads.IsEmpty()){
    a->SetLoads(loads);
  }
  // a frequency cannot really be applied but has to be used as parameter
  // in the driver call
}

double Excitation::GetOmega()
{
  if(frequency < 0.0)
    EXCEPTION("No frequency given");
  return 2 * M_PI * frequency;
}

double Excitation::GetFactor(Function* cost)
{
  double factor = 1.0; // default

  if(cost->FactorOmegaOmega())
  {
    if(frequency < 0.0) EXCEPTION("No frequency given");
    factor *= 4 * M_PI * M_PI * frequency * frequency;
  }

  return factor;
}

double Excitation::GetWeightedFactor(Function* cost)
{
  return normalized_weight * GetFactor(cost);
}

void Excitation::ReadTrackings(PtrParamNode ts){
  StdPDE* mech = domain->GetStdPDE("mechanic");
  trackings.Clear();
  ParamNodeList tracking_list = ts->GetChildren();
  for(unsigned int i = 0; i < tracking_list.GetSize(); i++){
    std::string name, dof, value;
    dof = "";
    tracking_list[i]->GetValue("name", name);
    tracking_list[i]->GetValue("dof", dof, ParamNode::PASS);
    tracking_list[i]->GetValue("value", value);
    shared_ptr<LoadBc> actLoad( new LoadBc );
    shared_ptr<EntityList> actList = domain->GetGrid()->GetEntityList( EntityList::NODE_LIST, name, EntityList::NAMED_NODES );
    actLoad->entities = actList;
    actLoad->eqnMap = mech->GetEqnMap();
    if ( dof == "" ) {
      actLoad->dof = 1;
    } else {
      actLoad->dof = mech->GetResultInfo(MECH_DISPLACEMENT)->GetDofIndex(dof);
    }
    actLoad->value = value;
    trackings.Push_back(actLoad);
  }
}

void Excitation::ReadLoads(PtrParamNode ls)
{
  // loads
  loads.Clear();
  MechPDE* mech = (MechPDE*)domain->GetSinglePDE("mechanic");
  mech->ReadLoads(ls->GetList("load"), loads);

  // pressures
  StdVector<shared_ptr<EntityList> > pressSurf;
  StdVector<std::string> pressVals;
  StdVector<std::string> pressPhase;
  mech->ReadPressureLoadsFromXML(ls, pressSurf, pressVals, pressPhase);
  mech->DefinePressureIntegrators(pressSurf, pressVals, pressPhase, linForms);

  // regionLoads
  std::map<RegionIdType, SinglePDE::RegionLoad> regionLoads;
  mech->ReadRegionLoadsFromXML(ls, regionLoads);
  mech->DefineRegionLoadIntegrators(regionLoads, linForms);
}

void Excitation::ReadTestStrain(const Vector<double>& vec)
{
  this->test_strain = vec;

  loads.Clear();
  MechPDE* mech = dynamic_cast<MechPDE*>(domain->GetSinglePDE("mechanic"));
  mech->DefineTestStrainIntegrators(vec, linForms);
}

void Excitation::SetPolarizationMatrixRHS(const Vector<double>& mechp,
    const Vector<double>& elecp, const int num)
{
  loads.Clear();
  MechPDE* mech = dynamic_cast<MechPDE*>(domain->GetSinglePDE("mechanic"));
  mech->DefinePolarizationMatrixIntegrators(mechp, linForms, num);

  ElecPDE* elec = dynamic_cast<ElecPDE*>(domain->GetSinglePDE("electrostatic"));
  elec->DefinePolarizationMatrixIntegrators(elecp, linForms, num);
}


