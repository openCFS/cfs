
#include <assert.h>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <list>
#include <map>
#include <sstream>

#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/XmlReader.hh"
#include "Domain/Domain.hh"
#include "Domain/ElemMapping/Elem.hh"
#include "Domain/Mesh/Grid.hh"
#include "Driver/BucklingDriver.hh"
#include "Driver/EigenFrequencyDriver.hh"
#include "Driver/FormsContexts.hh"
#include "General/Environment.hh"
#include "General/Exception.hh"
#include "MatVec/Matrix.hh"
#include "Optimization/Condition.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Design/ShapeDesign.hh"
#include "Optimization/ErsatzMaterial.hh"
#include "Optimization/Excitation.hh"
#include "Optimization/Optimization.hh"
#include "Utils/ToolsFull.hh"
#include <boost/lexical_cast.hpp>

#if defined(WIN32) && defined(__INTEL_COMPILER)
  #include <ciso646>
#endif

using std::string;
using std::pair;
using std::get;

using namespace CoupledField;
namespace CoupledField {class DesignStructure;}

DEFINE_LOG(conditions, "conditions")

// instantiation of the static elements
Enum<Condition::Bound> Condition::bound;
double Condition::SLACK_VALUE = -45217861;
double Condition::ALPHA_VALUE = -45217858;
double Condition::ALPHA_MINUS_SLACK_VALUE = -45217860;
double Condition::ALPHA_PLUS_SLACK_VALUE = -45217859;

Condition::Condition(PtrParamNode pn) : Function(pn)
{
  volume_fraction = 0.0;
  blown_up_ = false;
  index_ = -1; // to be set by ConditionContainer::Read()
  virtual_base_index_ = -1;

  observation_ = pn->Get("mode")->As<string>() == "observation";

  // the bound value is mandatory when we have a constraint
  if(!observation_ && !pn->Has("bound") && type_ != ISOTROPY && type_ != ISO_ORTHOTROPY && type_ != ORTHOTROPY)
    throw Exception("bound type for constraint '" + type.ToString(type_) + "' mandatory");
  bound_ = pn->Has("bound") ? bound.Parse(pn->Get("bound")->As<string>()) : EQUAL;
  // the bound value is called value in the problem file!
  // there must not  be a value when a homogenization tensor is given
  this->boundValue_ = -1.0;

  // special handling of scaling
  objective_scaling_ = pn->Get("scaling")->As<string>() == "objective";
  manual_scaling_value = objective_scaling_ ? -1.0 : pn->Get("scaling")->As<double>();


  delta_logging_ignored_ = false;
  delta_logging = pn->Get("log_delta")->As<bool>();

  if(pn->Get("log_delta")->As<bool>() && (!pn->Has("value") && !pn->Has("tensor") && !pn->Has("isotropic")))
    delta_logging_ignored_ = true;

  if(pn->Has("coord"))
    ReadCoord(pn);

  penalty = pn->Get("penalty")->As<double>();

  // validated in StressConstraint::GetApplications()
  stressType_ = stressType.Parse(pn->Get("stress")->As<string>());

  // set number of displacement constraints realized by multiple output constraints
  if (pn->Has("output") && pn->Get("output")->Has("displacement") && pn->Get("output")->Get("displacement")->Has("multiple_nodes"))
    output_multiple_nodes = pn->Get("output")->Get("displacement")->Get("multiple_nodes")->As<double>();
  else
    output_multiple_nodes = 0;

  bloch_extremal_ = false; // set in the proper case

  if(pn->Has("value"))
  {
    string v = pn->Get("value")->As<string>();
    if(v == "slack")
      this->boundValue_ = SLACK_VALUE;
    else if (v == "alpha")
      this->boundValue_ = ALPHA_VALUE;
    else if (v == "alpha+slack")
      this->boundValue_ = ALPHA_PLUS_SLACK_VALUE;
    else if (v == "alpha-slack")
      this->boundValue_ = ALPHA_MINUS_SLACK_VALUE;
    else
    {
      // interpret the value as expression to allow "1/nx". Does not evaluate each function evaluation
      // for this the handle needs to be stored in the function and care must be taken for optimizer interface and
      // local function performance
      this->boundValue_ = pn->Get("value")->MathParse<double>();
      LOG_DBG(conditions) << "C: " << type.ToString(type_) << " p=" << penalty << " os=" << objective_scaling_ << " msv=" << manual_scaling_value
                          << " bv=" << boundValue_ << " -> " << (boundValue_ * manual_scaling_value);
      if(!objective_scaling_)
        this->boundValue_ *= manual_scaling_value;
    }

    if((boundValue_ == ALPHA_PLUS_SLACK_VALUE && bound_ == UPPER_BOUND) || (boundValue_ == ALPHA_MINUS_SLACK_VALUE && bound_ == LOWER_BOUND)) {
      std::string msg =  "are you sure about value '" + v + "' and bound '" + bound.ToString(bound_) + "' in constraint '" + ToString() + "'?";
      domain->GetInfoRoot()->Get("optimization")->Get(ParamNode::HEADER)->Get("constraints")->SetWarning(msg, true); // domain->GetOptimization() does not work yet!
    }
  }


  // value is not mandatory for almost all constraints. Check for homogenization later
  if(!observation_)
  {
    switch(type_)
    {
    case HOM_TENSOR:
    case HOM_TRACKING:
      if(!pn->Has("tensor") && !pn->Has("value") && !pn->Has("isotropic"))
        throw Exception("Neither value nor tensor is given for constraint '" + type.ToString(type_) + "'");
      break;
    case HOM_FROBENIUS_PRODUCT:
      if(!pn->Has("tensor"))
        throw Exception("Tensor is mandatory for '" + type.ToString(type_) + "'");
      break;
    // case FMO_POS_DEF:
    case POS_DEF_DET_MINOR_1:
    case POS_DEF_DET_MINOR_2:
    case POS_DEF_DET_MINOR_3:
    case BENSON_VANDERBEI_1:
    case BENSON_VANDERBEI_2:
    case BENSON_VANDERBEI_3:
    if(!pn->Has("parameter"))
        throw Exception("'parameter' (very small value) mandatory for '" + type.ToString(type_) + "'");
      break;
    case ISOTROPY:
    case ISO_ORTHOTROPY:
    case ORTHOTROPY:
      if(pn->Has("value"))
        throw Exception("No value allowed for constraint '" + type.ToString(type_) + "'");
      break; // ok without value
    case EIGENFREQUENCY:
      if(Optimization::context->DoBloch()) {
        if(!pn->Has("bloch"))
           throw Exception("For Bloch optimization constraints '" + type.ToString(type_) + "' require the 'bloch' attribute to be set");
        bloch_extremal_ = pn->Get("bloch")->As<string>() == "extremal";
      }
      break;
    case EXPRESSION:
      if(!pn->Has("parameter"))
        throw Exception("'parameter' mandatory for '" + type.ToString(type_) + "' to formulate e.g. 'parameter' larger alpha-slack");
      // warn about boundValue != alpha +/- slack in ToInfo
     break;

    default:
      if(!pn->Has("value"))
        throw Exception("No value given for constraint '" + type.ToString(type_) + "'");
      break;
    }
  }
}

void Condition::PostProc(DesignSpace* space, DesignStructure* structure, ErsatzMaterial* em)
{
  SetElements(space, region); // before Function::PostProc() because of virtual_elem_map

  if(type_ == DESIGN_TRACKING)
    ReadDesignTrackingPattern(space, structure);

  if(boundValue_ == ALPHA_VALUE|| boundValue_ == ALPHA_MINUS_SLACK_VALUE || boundValue_ == ALPHA_PLUS_SLACK_VALUE) {
    if(!space->HasAlphaVariable())
      throw Exception("design variable 'alpha' is missing.");
    if(!space->HasSlackVariable())
      throw Exception("design variable 'alpha' requires also design variable 'slack'.");
  }


  // shall not be necessary when we register all pdes!
  //if((type_ == STRESS) && stressType_ != App::MECH)
  // {
    // it might be that we do piezo stresses on a pure elastic optimization problem.
    // Then register the App::ELEC PDE such that it is stored for the stress calculation by StressConstraint()
    // if we do PiezoSIMP this is simply redundant
    // Optimization::context->pdes[App::ELEC] = domain->GetSinglePDE("electrostatic");
  // }

  // note, meanwhile we have info_ set! but not yet in the constructor
  Function::PostProc(space, structure, em);

}

bool Condition::ReadCoord(PtrParamNode pn)
{
  string val = pn->Get("coord")->As<string>();
  if(val == "all") return false;

  assert(val.size() == 2);
  coords.Resize(1);

  ParseCoord(pn, coords[0]);

  return true;
}



void Condition::AddCondition(PtrParamNode pn, StdVector<Condition*>& list, int i, std::string entName)
{
  Type t = type.Parse(pn->Get("type")->As<string>());
  list.Push_back(IsLocal(t) ? new LocalCondition(pn) : new Condition(pn));

  // note that the pointer becomes invalid by AddSubCondition()
  Condition* g = list.Last();
  g->index_ = -1; // to be defined by the virtual "all" list in the ConditionContainer::Read() method

  // homogenization are special constraints which constraints which "blow up" to more constraints
  if(g->type_ == HOM_TENSOR)
    AddHomogenizationTensorConstraints(pn, list, g);

  // isotropy is a special constraint which blows up special tensor entry constraints
  if(g->type_ == ISOTROPY || g->type_ == ISO_ORTHOTROPY || g->type_ == ORTHOTROPY)
    AddXtropyConstraints(pn, list, g);

  // if OUTPUT is defined and the multiple_node option is turned on, multiple constraints are added to represent displacement constraints
  if(g->type_ == OUTPUT && i > 0)
    AddOutputConstraints(pn,list,g,i,entName);



  //if(g->type_ == FMO_POS_DEF_MINOR_1 || FMO_POS_DEF_MINOR_2 || POS_DEF_DET_MINOR_3)
  //  AddFMOPosDefConstraints(pn, list, g);
}

// modify ParamNode pn of constraint, add number i to node name of output constraint. Necessary for automatic numbering of displacement constraints with multiple_node option
void Condition::AddOutputConstraints(PtrParamNode pn, StdVector<Condition*>& list, Condition* g,int i,std::string entName) {
  assert(g->GetType() == OUTPUT && i > 0);

  PtrParamNode output;
  ParamNodeList elems;
  if (pn->Has("output"))
    output = pn->Get("output");
  if (output->Has("displacement"))
    elems = output->GetList ("displacement");

  // add number i to node name of output constraint
  assert(elems.GetSize() == 1);
  PtrParamNode xml = elems[0];
  //std::string entName = xml->Get("name")->As<std::string>();
  entName.assign(entName + boost::lexical_cast<std::string>(i));
  xml->Get("name")->SetValue(entName);
}


void Condition::AddXtropyConstraints(PtrParamNode pn, StdVector<Condition*>& list, Condition* g)
{
  // isotropy is a special constraint which blows up special tensor entry constraints
  assert(g->GetType() == ISOTROPY || g->GetType() == ISO_ORTHOTROPY || g->type_ == ORTHOTROPY);

  // we reset the type, therefore keep it
  Type org = g->type_;

  if(pn->Has("coord"))
    throw Exception("don't use attribute 'coord' for constraint 'isotropy'/'iso-orthotropy'/'orthotropy'");

  if(g->bound_ != EQUAL)
    throw Exception("the 'isotropy'/'iso-orthotropy'/'orthotropy' constraint requires equality constraint type");

  // become an HOM_TENSOR constraint!
  g->type_ = HOM_TENSOR;

  // isotropic specific is a condition relating E11, shear and off-diagonal. See Lame-notation in Richter
  // all constraints in isotropic are equal zero.
  //
  // iso-orthotropic is orthotropic with all E, G and all nu equal. E and nu are related but not G
  // notation:
  // E = E(1,1) (= E(2,2) = E(3,3))
  // nu = E(1,2) (= E(1,3) = E(2,3))
  // G = E(4,4) (= E(5,5) = E(6,6))
  // this conditions are lower bound conditions!

  if(domain->GetDim() == 2)
  {
    assert(g->coords.GetSize() == 0);
    g->boundValue_ = 0;

    if(org == ISOTROPY || org == ISO_ORTHOTROPY)
    {
      // E11 - E22 = 0
      g->coords.Push_back(std::make_tuple(1,1,1.0)); // no AppendSubCondtion()/Clear() in first case
      g->coords.Push_back(std::make_tuple(2,2,-1.0));

      g = g->AppendSubCondition(list);
      g->coords.Clear(); // do here, see below
    }

    // E13 = 0
    // AppendSubCondtion()/Clear() not in ISO_ORTHOTROPY case and done otherwise above.
    // this is very ugly, but switching the order makes problems for SCPIP to solve trivial max E11 in the testsuite!!
    g->coords.Push_back(std::make_tuple(1,3,1.0));

    // E23 = 0
    g = g->AppendSubCondition(list);
    g->coords.Clear();
    g->coords.Push_back(std::make_tuple(2,3,1.0));


    if(org == ISOTROPY)
    {
      // E11 - E12 - 2E33 = E11 - E12 - E33 - E33 = 0
      g = g->AppendSubCondition(list);
      g->coords.Clear(); // the copy constructor above copies old stuff
      g->coords.Push_back(std::make_tuple(1,1,1.0));
      g->coords.Push_back(std::make_tuple(1,2,-1.0));
      g->coords.Push_back(std::make_tuple(3,3,-2.0));
    } // else case is common for 2D and 3D
  }
  else
  {
    assert(g->coords.GetSize() == 0);

    g->boundValue_ = 0;

    if(org == ISOTROPY || org == ISO_ORTHOTROPY)
    {
      // non-shear diagonal is constant
      // E11 = E22 = E33 -> E11 - E22 = 0, E22 - E33 = 0
      g->coords.Push_back(std::make_tuple(1,1,1.0));
      g->coords.Push_back(std::make_tuple(2,2,-1.0));

      g = g->AppendSubCondition(list);
      g->coords.Clear();
      g->coords.Push_back(std::make_tuple(2,2,1.0));
      g->coords.Push_back(std::make_tuple(3,3,-1.0));

      // upper non-shear triangle is constant
      // E12 = E13 = E23 -> E12 - E13 = 0, E13 - E23 = 0
      g = g->AppendSubCondition(list);
      g->coords.Clear();
      g->coords.Push_back(std::make_tuple(1,2,1.0));
      g->coords.Push_back(std::make_tuple(1,3,-1.0));

      g = g->AppendSubCondition(list);
      g->coords.Clear();
      g->coords.Push_back(std::make_tuple(1,3,1.0));
      g->coords.Push_back(std::make_tuple(2,3,-1.0));

      // shear diagonal is constant
      // E44 = E55 = E66 -> E44 - E55 = 0, E55 - E66 = 0
      g = g->AppendSubCondition(list);
      g->coords.Clear();
      g->coords.Push_back(std::make_tuple(4,4,1.0));
      g->coords.Push_back(std::make_tuple(5,5,-1.0));

      g = g->AppendSubCondition(list);
      g->coords.Clear();
      g->coords.Push_back(std::make_tuple(5,5,1.0));
      g->coords.Push_back(std::make_tuple(6,6,-1.0));

      g = g->AppendSubCondition(list);
      g->coords.Clear(); // see 2D why we have to do this nonsense
    }


    // the zero entries
    // E14 = E15 = E16 = E24 = E25 = E26 = E34 = E35 = E36 = E45 = E46 = E56 = 0
    g->coords.Push_back(std::make_tuple(1,4,1.0));

    g = g->AppendSubCondition(list);
    g->coords.Clear();
    g->coords.Push_back(std::make_tuple(1,5,1.0));

    g = g->AppendSubCondition(list);
    g->coords.Clear();
    g->coords.Push_back(std::make_tuple(1,6,1.0));

    g = g->AppendSubCondition(list);
    g->coords.Clear();
    g->coords.Push_back(std::make_tuple(2,4,1.0));

    g = g->AppendSubCondition(list);
    g->coords.Clear();
    g->coords.Push_back(std::make_tuple(2,5,1.0));

    g = g->AppendSubCondition(list);
    g->coords.Clear();
    g->coords.Push_back(std::make_tuple(2,6,1.0));

    g = g->AppendSubCondition(list);
    g->coords.Clear();
    g->coords.Push_back(std::make_tuple(3,4,1.0));

    g = g->AppendSubCondition(list);
    g->coords.Clear();
    g->coords.Push_back(std::make_tuple(3,5,1.0));

    g = g->AppendSubCondition(list);
    g->coords.Clear();
    g->coords.Push_back(std::make_tuple(3,6,1.0));

    g = g->AppendSubCondition(list);
    g->coords.Clear();
    g->coords.Push_back(std::make_tuple(4,5,1.0));

    g = g->AppendSubCondition(list);
    g->coords.Clear();
    g->coords.Push_back(std::make_tuple(4,6,1.0));

    g = g->AppendSubCondition(list);
    g->coords.Clear();
    g->coords.Push_back(std::make_tuple(5,6,1.0));

    if(org == ISOTROPY)
    {
      // relationship of the three unique values
      // E11 - E12 - 2E66 = E11 - E12 - E66 - E66 = 0
      g = g->AppendSubCondition(list);
      g->coords.Clear(); // the copy constructor above copies old stuff
      g->coords.Push_back(std::make_tuple(1,1,1.0));
      g->coords.Push_back(std::make_tuple(1,2,-1.0));
      g->coords.Push_back(std::make_tuple(6,6,-2.0));
    }

  }
}

void Condition::AddHomogenizationTensorConstraints(PtrParamNode pn, StdVector<Condition*>& list, Condition* g)
{
  // homogenization are special constraints which constraints which "blow up" to more constraints
  assert(g->GetType() == HOM_TENSOR);
  // there has been a extensive test in the constructor
  // do we need to blow-up?
  if(!g->ReadCoord(pn)) // this is done if coord="all" is given in xml!
  {
    // the first entry is this constraint
    // for the conversion of the indices see Thesis from Ole Sigmund, p 30 and book of Manfred
    // p 42 (first edition)
    g->coords.Resize(1);
    std::tuple<int, int, double>& entry = g->coords[0];
    get<0>(entry) = 1; // 1. ijkl = 1111
    get<1>(entry) = 1;
    get<2>(entry) = 1.0;
    g->boundValue_ = g->tensor_[1-1][1-1]; // one-based!
    g = g->AppendSubCondition(list, 2,2); // 2. 2222
    if(domain->GetDim() == 2)
    {
      g = g->AppendSubCondition(list, 1,2); // 3. 1122

      //g = g->AppendSubCondition(list, 3,1); // no covered by Sigmund
      //g = g->AppendSubCondition(list, 3,2); // no covered by Sigmund

      g = g->AppendSubCondition(list, 3,3); // 4. 1212 // would be 6,6 according to book

      g = g->AppendSubCondition(list, 1,3); // 4. 1122
      g = g->AppendSubCondition(list, 2,3); // 5. 2233
    }
    else
    {
      g = g->AppendSubCondition(list, 3,3); // 3. 3333
      g = g->AppendSubCondition(list, 1,2); // 4. 1122
      g = g->AppendSubCondition(list, 1,3); // 5. 1133
      g = g->AppendSubCondition(list, 2,3); // 6. 2233
      g = g->AppendSubCondition(list, 6,6); // 7. 1212
      g = g->AppendSubCondition(list, 5,5); // 8. 1313
      g = g->AppendSubCondition(list, 4,4); // 9. 2323

      g = g->AppendSubCondition(list, 1,4);
      g = g->AppendSubCondition(list, 1,5);
      g = g->AppendSubCondition(list, 1,6);
      g = g->AppendSubCondition(list, 2,4);
      g = g->AppendSubCondition(list, 2,5);
      g = g->AppendSubCondition(list, 2,6);
      g = g->AppendSubCondition(list, 3,4);
      g = g->AppendSubCondition(list, 3,5);
      g = g->AppendSubCondition(list, 3,6);

      g = g->AppendSubCondition(list, 4,5);
      g = g->AppendSubCondition(list, 4,6);
      g = g->AppendSubCondition(list, 5,6);
    }
  }
}

void Condition::AddExcitationStressConstraints(StdVector<Condition*>& list, MultipleExcitation* me)
{
  // no multiple excitations, no additional stress constraints
  if(!me->IsEnabled())
    return;

  // do we need to blow up the list? This is the case when there is stress constraint with
  // the default excitation attribute all which is Function::excite_ = -1.
  // Otherwise we rely on the opt_unique_constraint xsd:unique constraint
  int blow_up = -1;
  for(unsigned int i = 0; i < list.GetSize(); i++)
  {
    switch(list[i]->GetType())
    {
    case GLOBAL_STRESS:
    case LOCAL_STRESS:
      assert(!Optimization::context->DoMultiSequence());
      if(list[i]->DoEvaluateAlways(1)) // sequence 1
        blow_up = i;
      else
        if(blow_up != -1)
          throw Exception("You cannot mix stress constraints with excitation and with default 'all' excitation");
      break;
    default:
      break;
    }
  }

  // are there stress constraints to blow up?
  if(blow_up == -1) return;

  Condition& g = *(list[blow_up]);
  g.SetExcitation(me, me->excitations[0].index);
  assert(!Optimization::context->DoMultiSequence());
  for(unsigned int e = 1; e < me->excitations.GetSize(); e++)
  {
    Condition* tmp = new Condition(g);
    tmp->SetExcitation(me, me->excitations[e].index);

    list.Insert(blow_up + e, tmp);
  }
}


void Condition::AddBlochEigenConstraints(StdVector<Condition*>& all_cond, MultipleExcitation* me)
{
  // this is a static function
  for(unsigned int c = 0; me->IsEnabled() && c < Optimization::manager.context.GetSize(); c++)
  {
    Context& ctxt = Optimization::manager.context[c];

    if(ctxt.DoBloch())
    {
      // we need to find all eigenvalue constraints. Then extend the ones by excitation which are bloch=full

      // extract all eigenfrequency constraints to become full to full_ev
      StdVector<Condition> full_ev; // instances as list will be enlarged which involves copying
      for(unsigned int i = 0; i < all_cond.GetSize(); i++)
      {
        if(all_cond[i]->GetType() == EIGENFREQUENCY)
        {
          Condition* g = all_cond[i];
          assert(g->ctxt->sequence == ctxt.sequence);
          assert(ctxt.excitations[0]->index >= 0);

          // expand the the constraints to all wave vectors?
          if(g->DoFullBloch())
          {
            // reset excitation to the first wave vector
            g->SetExcitation(me, ctxt.excitations[0]->index);
            full_ev.Push_back(*g);
            LOG_DBG(conditions) << "ABEC: seq=" << ctxt.sequence << " i=" << i << " ev=" << full_ev.GetSize() << " ex=" << ctxt.excitations[0]->index << " -> " <<  full_ev.Last().ToString();
          }
          else
          {
            // we evaluate the function at the very last wave vector as we have to search for the extremals
            g->SetExcitation(me, ctxt.excitations.Last()->index);
            LOG_DBG(conditions) << "ABEC: seq=" << ctxt.sequence << " i=" << i << " g=" << g->ToString();
          }

        }
      }

      // expand only for bloch=full
      if(!full_ev.IsEmpty())
      {

        assert(ctxt.num_bloch_wave_vectors * me->GetNumberRobust(&ctxt, true) == ctxt.excitations.GetSize());
        for(unsigned int e = 1; e < ctxt.excitations.GetSize(); e++) // start from 1!
        {
          LOG_DBG2(conditions) << "ABEC: e=" << e << " -> " << ctxt.excitations[e]->index;
          // note that we traverse ev and not list again!
          for(unsigned int g = 0; g < full_ev.GetSize(); g++)
          {
            assert(full_ev[g].IsExcitationSensitive());
            assert(full_ev[g].GetExcitation()->index >= 0);

            Condition* tmp = new Condition(full_ev[g]);
            tmp->SetExcitation(me, ctxt.excitations[e]->index);
            LOG_DBG2(conditions) << "ABEC: e=" << e << " g=" << g << " -> " << tmp->ToString();
            assert(ctxt.excitations[e]->index >= 0);

            all_cond.Push_back(tmp);
          }
        }
      }
    }
  }
}



Condition* Condition::AppendSubCondition(StdVector<Condition*>& list, bool biisotropy, bool imag)
{
  if(this->IsLocalCondition())
    list.Push_back(new LocalCondition(*dynamic_cast<LocalCondition*>(this)));
  else
    list.Push_back(new Condition(*this)); // make a copy of this element by the (default) copy constructor

  Condition* sub = list.Last(); // copy this entry as reference
  sub->index_ = list.GetSize() - 1;
  sub->blown_up_ = true;
  imag_ = imag;
  biisotropy_ = biisotropy;
  return sub;
}


Condition* Condition::AppendSubCondition(StdVector<Condition*>& list, int pos_x, int pos_y)
{
  Condition* sub = AppendSubCondition(list);
  sub->coords.Resize(1);
  std::tuple<int, int, double>& entry = sub->coords[0];
  get<0>(entry) = pos_x;
  get<1>(entry) = pos_y;
  get<2>(entry) = 1.0; // default
  sub->boundValue_ = sub->tensor_[pos_x-1][pos_y-1];
  sub->blown_up_ = true;
  return sub;
}

void Condition::ReadDesignTrackingPattern(DesignSpace* space, DesignStructure* structure)
{
  assert(type_ == DESIGN_TRACKING);
  assert(elements.GetSize() > 0); // SetElements() needs to be called prior this one

  // if elements is not set by a region we overwrite elements and search for the elements with periodic
  // boundary conditions (the outer frame)
  if(region == ALL_REGIONS)
  {
    // ensure it is initialized
    VicinityElement::Init(space, structure);

    assert(elements[0]->vicinity != NULL); // it shall not be a ghost element
    elements.Resize(0); // capacity is still there so we can push back
    for(unsigned int i = 0; i < space->data.GetSize(); i++)
    {
      DesignElement& de = space->data[i];
      if(de.vicinity->periodic)
        elements.Push_back(&de);
    }

    if(elements.GetSize() == 0)
      throw Exception("Constraint 'designTracking' requires attribute 'region' when there are no periodic boundary conditions");
  }

  // read the pattern file
  if(!pn->Has("designTarget"))
    throw Exception("Attribute 'designTarget' holding a density file name is mandatory of 'designTracking'");
  string file = pn->Get("designTarget")->As<string>();
  PtrParamNode xml = XmlReader::ParseFile(file);

  // check this file
  if (xml->Count("set") == 0)
    throw Exception("There are no design sets in the pattern file " + file);

  // read the target in a huge temporary list such that it is cheap to compare against the design elements
  unsigned int grid_size = domain->GetGrid()->GetNumElems();
  StdVector<double> tmp;
  tmp.Resize(grid_size + 1, 0.0);

  ParamNodeList elems = xml->GetList("set").Last()->GetList("element");
  if(elems.GetSize() > grid_size)
    EXCEPTION("The 'designTarget' file '" << file << "' has " << elems.GetSize() << " elements and the mesh only " << grid_size)
  if(!elems[0]->Has("physical"))
    throw Exception("'designTracking' requires the attribute 'physical' in the 'designTarget' file " + file);

  for(unsigned int i = 0; i < elems.GetSize(); i++)
    tmp[elems[i]->Get("nr")->As<int>()] = elems[i]->Get("physical")->As<double>();

  // copy from pattern what we actually need
  pattern.Resize(elements.GetSize());
  for(unsigned int i = 0, n = elements.GetSize(); i < n; i++)
    pattern[i] = tmp[elements[i]->elem->elemNum];
}


double Condition::CalcFeasibility() const
{
  double diff = GetValue() - boundValue_; // handles also local constraints!
  return diff;
}

bool Condition::IsFeasible() const
{
  double diff = CalcFeasibility();

  switch(bound_)
  {
    case EQUAL:
      return std::abs(diff) < 1e-4;

    case LOWER_BOUND:
      return diff + 1e-4 > 0;

    case UPPER_BOUND:
      return diff - 1e-4 < 0;
  }

  assert(false);
  return false;
}

bool Condition::IsFeasibilityConstraint() const
{
  switch(type_)
  {
  case POS_DEF_DET_MINOR_1:
  case POS_DEF_DET_MINOR_2:
  case POS_DEF_DET_MINOR_3:
  case BENSON_VANDERBEI_1:
  case BENSON_VANDERBEI_2:
  case BENSON_VANDERBEI_3:
  case DESIGN:
    return true;
  default:
    return false;
  }
}


string Condition::ToString() const
{
  std::ostringstream os;
  
  if(delta_logging) os << "delta_";

  os << Function::ToString(); // includes physical

  if(region != ALL_REGIONS)
    os << "_" << domain->GetGrid()->GetRegion().ToString(region);

  if(design_ != DesignElement::DEFAULT)
    os << "_(" << DesignElement::type.ToString(design_) << ")";

  if(type_ == HOM_TENSOR)
    os << "_" << ToString(coords);

  // with multiple output constraints we need to identify
  if((type_ == OUTPUT || type_ == SQUARED_OUTPUT) && !output_forms.IsEmpty())
    os << "_" << output_forms[0]->GetEntities()->GetName();

  // e.g. stresses are extended for every excitation
  if(GetExcitation() != NULL && domain->GetOptimization()->GetMultipleExcitation()->IsEnabled())
  {
    if(type_ == GLOBAL_STRESS || type_ == LOCAL_STRESS)
      os << "_" << GetExcitation()->GetFullLabel(); // change to excite label
    else if(domain->GetOptimization()->GetMultipleExcitation()->DoMetaExcitation(GetExcitation()->sequence))
      os << "_" << GetExcitation()->GetMetaLabel();
  }

  if(type_ == EIGENFREQUENCY)
    os << "_" << eigenvalue_id_;

  if(type_ == EIGENFREQUENCY && GetExcitation() != NULL && GetExcitation()->DoBloch()) // might not be set meantime - e.g. due to early logging
  {
    if(DoFullBloch())
      os << "_wv_" << GetExcitation()->GetWaveNumber();
    else
      os << "_" << (bound_ == Condition::LOWER_BOUND ? "min" : "max");
  }

  // add bound type if multiple unique conditions exist
  if(domain->GetOptimization()->constraints.RequiresBoundForUniqueness(this))
    os << "_" << bound.ToString(bound_);

  return os.str();  
}

string Condition::ToString(const StdVector<std::tuple<int, int, double> >& coords)
{
  assert(coords.GetSize() > 0);
  assert(get<2>(coords[0]) == 1.0); // so we don't have to start with a minus

  std::ostringstream os;

  // 11_p12_2m33
  for(unsigned int i = 0; i < coords.GetSize(); i++)
  {
    const std::tuple<int, int, double>& entry = coords[i];
    assert(std::floor(get<2>(entry)) - get<2>(entry) == 0.0);
    int factor = static_cast<int>(get<2>(entry));
    assert((i == 0 && factor == 1) || i != 0);
    if(std::abs(factor) != 1.0)
      os << factor;
    if(factor < 0)
      os << "m";
    else
      os << "p";

    os << get<0>(entry) << get<1>(entry);

    if(i < coords.GetSize() - 1)
      os << "_";
  }

  return os.str();
}


void Condition::ToInfo(PtrParamNode in)
{
  Function::ToInfo(in);

  if(IsActive())
  {
    if(type_ != HOM_TRACKING)
    {
      if(boundValue_ == SLACK_VALUE)
        in->Get("bound_value")->SetValue("slack");
      else if (boundValue_ == ALPHA_VALUE)
        in->Get("bound_value")->SetValue("alpha");
      else if(boundValue_ == ALPHA_MINUS_SLACK_VALUE)
        in->Get("bound_value")->SetValue("alpha-slack");
      else if(boundValue_ == ALPHA_PLUS_SLACK_VALUE)
        in->Get("bound_value")->SetValue("alpha+slack");
      else
        // FIXME: does not handle objective_scaling. Also the scaling shall not be encoded in the bound value :(
        in->Get("bound_value")->SetValue(boundValue_ / manual_scaling_value);
    }
    in->Get("bound")->SetValue(bound.ToString(bound_));

    if(objective_scaling_)
      in->Get("scaling")->SetValue("objective");
    else if(manual_scaling_value != 1.0)
      in->Get("scaling")->SetValue(manual_scaling_value);
  }
  if(type_ == HOM_TENSOR)
    in->Get("tensor_entry")->SetValue(ToString(coords));

  if(observation_)
      in->Get("mode")->SetValue("observation");

  in->Get("design")->SetValue(DesignElement::type.ToString(design_));

  if(type_ == DESIGN_TRACKING)
    in->Get("elements")->SetValue(elements.GetSize());

  if(type_ == GLOBAL_STRESS || type_ == LOCAL_STRESS)
    in->Get("stress")->SetValue(stressType.ToString(stressType_));

  if(type_ == EIGENFREQUENCY && GetExcitation()->DoBloch())
    in->Get("bloch")->SetValue(bloch_extremal_ ? "extremal" : "full");

  if(type_ == EXPRESSION && (boundValue_ != ALPHA_MINUS_SLACK_VALUE && boundValue_ != ALPHA_PLUS_SLACK_VALUE))
   info_->SetWarning("be sure to know what condition 'expression' with alpha+/-slack bound means");

  if(domain->GetOptimization()->GetMultipleExcitation()->IsEnabled())
  {
    if(DoEvaluateAlways(ctxt->sequence))
      in->Get("excitation")->SetValue(Optimization::context->DoMultiSequence() ? "always within sequence" : "always");
    else
      in->Get("excitation")->SetValue(GetExcitation()->GetFullLabel());
  }

  // TODO somehow scaling does not work ??
  // if(IsHomogenization() && !objective_scaling_ && !blown_up_) // warn only the first time!
  //  in->SetWarning("Doing homogenization without 'objective' scaling constraint '" + type.ToString(type_) + "'");


  if(type_ == VOLUME && IsPhysical() && !observation_)
    info_->SetWarning("a physical volume constraint should make no sense");

  if((type_ == VOLUME || type_ == TENSOR_TRACE) && design_ == DesignElement::MECH_TRACE)
    info_->Get("notation")->SetValue(tensorNotation.ToString(notation_));

  // the bounds are essential as we have to flip sign!
  if(type_ == OVERHANG_HOR && bound_ != LOWER_BOUND)
    throw Exception("overhang constraints for horizontal structures restrict the lower boundary only and this boundary shall be steep enough -> 'lower_bound'");

  if(type_ == OVERHANG_VERT && bound_ != UPPER_BOUND)
    throw Exception("overhang constraints for vertical structures restrict the left boundary for left overhangs and vice versa. -> 'upper_bound'");
}

void Condition::DescribeProperties(StdVector<std::pair<string,string> >& map) const
{
  Function::DescribeProperties(map);
  map.Push_back(std::make_pair("constraint", observation_ ? "observation" : "applied"));
  map.Push_back(std::make_pair("bound", bound.ToString(bound_)));
  map.Push_back(std::make_pair("bound_value", std::to_string(boundValue_)));
}


bool Condition::IsForRegion(RegionIdType regionId)
{
  return(region == ALL_REGIONS || region == regionId);
}


LocalCondition::LocalCondition(PtrParamNode pn) : Condition(pn)
{
}


Function::Local::Identifier& LocalCondition::GetCurrentVirtualContext()
{
  assert(IsLocal());

  unsigned int idx = current_view_index_ - virtual_base_index_;
  return local->virtual_elem_map[idx];
}

const Function::Local::Identifier& LocalCondition::GetCurrentVirtualContext() const
{
  assert(IsLocal());

  unsigned int idx = current_view_index_ - virtual_base_index_;
  return local->virtual_elem_map[idx];
}


/** fixed (or mapped) feature variables are not part of the optimization space and drop out of
 * sparsity patterns: -1, otherwise the opt index (e.g. a 'distance' node on a partially fixed pill) */
static int OptIndexIfVariable(const BaseDesignElement* bde)
{
  const FeatureVariable* fv = dynamic_cast<const FeatureVariable*>(bde);
  return fv != nullptr ? fv->EffectiveOptIndex() : (int) bde->GetOptIndex();
}

unsigned int LocalCondition::GetSparsityPatternSize() const
{
  if(IsAdjointBased())
  {
    assert(type_ == LOCAL_STRESS || type_ == LOCAL_BUCKLING_LOAD_FACTOR); // the only known cases up to now
    assert(!jac_sparsity_.empty());
    return jac_sparsity_.GetSize();
  }

  if(this->isFiltered())
    return ((Condition*) this)->GetSparsityPattern().GetSize();
  else {
    // some local constraints have non-uniform neighbor size, e.g. curvature in shape mapping.
    // the element itself is not a neighbor of itself, therefore i starts at -1
    const Function::Local::Identifier& id = GetCurrentVirtualContext();
    unsigned int nnz = 0;
    for(int i = -1, nn = id.neighbor.GetSize(); i < nn; i++)
      if(OptIndexIfVariable(id.GetElement(i)) >= 0)
        nnz++;
    return nnz;
  }
}

StdVector<unsigned int>& LocalCondition::GetSparsityPattern()
{
  assert(IsLocal());

  // up to now only LOCAL_STRESS and LOCAL_BUCKLING_LOAD_FACTOR have a full gradient
  assert(!IsStateDependent() || type_ == LOCAL_STRESS || type_ == LOCAL_BUCKLING_LOAD_FACTOR);
  if(IsStateDependent())
    return Function::GetSparsityPattern();


  Function::Local::Identifier& id = GetCurrentVirtualContext();

    // we shall sort the indices
  std::list<unsigned int> indices;
  for(int i = -1 ; i < (int) id.neighbor.GetSize(); i++)
  {
    BaseDesignElement* bde = id.GetElement(i);
    assert(bde != NULL);
    // int other_idx = local->space->Find(de); // needs to be fast!
    int other_idx = OptIndexIfVariable(bde);
    if(other_idx < 0)
      continue; // fixed feature variables are constants
    indices.push_back(other_idx);

    if(this->isFiltered())
    {
      DesignElement* de = dynamic_cast<DesignElement*>(bde);
      if(de != NULL && !de->simp->filter.IsEmpty())
      {
        const StdVector<Filter::NeighbourElement> neighborhood = de->simp->filter[de->simp->DetermineFilterIndexNonInlined()].neighborhood;
        for(unsigned int j = 0, n = neighborhood.GetSize(); j < n; j++)
          indices.push_back(neighborhood[j].neighbour->GetOptIndex());
      }
    }
  }

  // sort and copy
  indices.sort();
  jac_sparsity_.Resize(0); // keeps capacity, hence Push_back is cheap
  for(std::list<unsigned int>::const_iterator it = indices.begin(); it != indices.end(); ++it)
    jac_sparsity_.Push_back(*it);

  LOG_DBG2(conditions) << "LC:GSP: " << ToString() << " a=" << access.ToString(this->GetAccess()) << " f=" << isFiltered()
                       << " nnz=" << jac_sparsity_.GetSize() << " -> " << jac_sparsity_.ToString();
  return jac_sparsity_;
}

Matrix<unsigned int>& LocalCondition::GetHessianSparsityPattern()
{
  assert(IsLocal());

  bool elec = GetDesignType() == DesignElement::DIELEC_ALL;

  Function::Local::Identifier& id = GetCurrentVirtualContext();

  switch(type_)
  {
  case POS_DEF_DET_MINOR_2:
  {
    DesignElement::Type t11 = elec ? DesignElement::DIELEC_11 : DesignElement::MECH_11;
    DesignElement::Type t12 = elec ? DesignElement::DIELEC_12 : DesignElement::MECH_12;
    DesignElement::Type t22 = elec ? DesignElement::DIELEC_22 : DesignElement::MECH_22;

    hess_sparsity_.Resize(2, 2);

    hess_sparsity_(0, 0) = id.GetElementByType(t11)->GetOptIndex();
    hess_sparsity_(0, 1) = id.GetElementByType(t22)->GetOptIndex();

    hess_sparsity_(1, 0) = id.GetElementByType(t12)->GetOptIndex();
    hess_sparsity_(1, 1) = id.GetElementByType(t12)->GetOptIndex();

    break;
  }
  case POS_DEF_DET_MINOR_3:
    assert(!elec);
    hess_sparsity_.Resize(12, 2);

    hess_sparsity_(0, 0) = id.GetElementByType(DesignElement::MECH_11)->GetOptIndex();
    hess_sparsity_(0, 1) = id.GetElementByType(DesignElement::MECH_22)->GetOptIndex();

    hess_sparsity_(1, 0) = id.GetElementByType(DesignElement::MECH_11)->GetOptIndex();
    hess_sparsity_(1, 1) = id.GetElementByType(DesignElement::MECH_23)->GetOptIndex();

    hess_sparsity_(2, 0) = id.GetElementByType(DesignElement::MECH_11)->GetOptIndex();
    hess_sparsity_(2, 1) = id.GetElementByType(DesignElement::MECH_33)->GetOptIndex();

    hess_sparsity_(3, 0) = id.GetElementByType(DesignElement::MECH_12)->GetOptIndex();
    hess_sparsity_(3, 1) = id.GetElementByType(DesignElement::MECH_12)->GetOptIndex();

    hess_sparsity_(4, 0) = id.GetElementByType(DesignElement::MECH_12)->GetOptIndex();
    hess_sparsity_(4, 1) = id.GetElementByType(DesignElement::MECH_13)->GetOptIndex();

    hess_sparsity_(5, 0) = id.GetElementByType(DesignElement::MECH_12)->GetOptIndex();
    hess_sparsity_(5, 1) = id.GetElementByType(DesignElement::MECH_23)->GetOptIndex();

    hess_sparsity_(6, 0) = id.GetElementByType(DesignElement::MECH_12)->GetOptIndex();
    hess_sparsity_(6, 1) = id.GetElementByType(DesignElement::MECH_33)->GetOptIndex();

    hess_sparsity_(7, 0) = id.GetElementByType(DesignElement::MECH_22)->GetOptIndex();
    hess_sparsity_(7, 1) = id.GetElementByType(DesignElement::MECH_13)->GetOptIndex();

    hess_sparsity_(8, 0) = id.GetElementByType(DesignElement::MECH_22)->GetOptIndex();
    hess_sparsity_(8, 1) = id.GetElementByType(DesignElement::MECH_33)->GetOptIndex();

    hess_sparsity_(9, 0) = id.GetElementByType(DesignElement::MECH_13)->GetOptIndex();
    hess_sparsity_(9, 1) = id.GetElementByType(DesignElement::MECH_13)->GetOptIndex();

    hess_sparsity_(10, 0) = id.GetElementByType(DesignElement::MECH_13)->GetOptIndex();
    hess_sparsity_(10, 1) = id.GetElementByType(DesignElement::MECH_23)->GetOptIndex();

    hess_sparsity_(11, 0) = id.GetElementByType(DesignElement::MECH_23)->GetOptIndex();
    hess_sparsity_(11, 1) = id.GetElementByType(DesignElement::MECH_23)->GetOptIndex();

    break;
  case DISTANCE:
  {
    // c = ||Q-P|| over the node coordinates [px py (pz) qx qy (qz)] (neighbor order -1,0,..., see
    // CalcDistance); fixed nodes are constants and drop out, the (row, col) opt-index pairs of the
    // free nodes remain, profile p is absent. CalcHessian fills the matching values.
    int nc = 2 * domain->GetDim(); // number of node coordinates
    int o[6];
    int nf = 0;
    for(int i = -1; i <= nc - 2; i++)
    {
      int oi = OptIndexIfVariable(id.GetElement(i));
      if(oi >= 0)
        o[nf++] = oi;
    }
    hess_sparsity_.Resize(nf * nf, 2);
    for(int i = 0, k = 0; i < nf; i++)
      for(int j = 0; j < nf; j++, k++)
      {
        hess_sparsity_(k, 0) = o[i];
        hess_sparsity_(k, 1) = o[j];
      }
    break;
  }
  default:
    hess_sparsity_.Resize(0, 0);
    break;
  }

  LOG_DBG3(conditions) << "LC:GHSP: g=" << ToString() << " -> " << hess_sparsity_.ToString() << " n=" << BaseDesignElement::ToString(id.neighbor, true);
  return hess_sparsity_;
}

void LocalCondition::CalcHessian(StdVector<double>& out, double factor)
{
  assert(IsLocal());


  assert(out.GetSize() == GetHessianSparsityPattern().GetNumRows());

  switch(type_)
  {
  case POS_DEF_DET_MINOR_2:
    // (6.69) in the diss of Sonja Lehmann
    out[0] = factor * 1;
    out[1] = factor * -2;
    break;
  case POS_DEF_DET_MINOR_3:
  {
    Matrix<double> E;
    // (E - v*I) >= gamma
    double v = GetParameter();
    double eps = 1.0 * GetBoundValue();

    //Function::Local::Identifier& id = GetCurrentVirtualContext();
    assert(false);
    // local->space->GetErsatzMaterialTensor(E, PLANE_STRAIN, dynamic_cast<DesignElement*>(id.element)->elem, DesignElement::NO_DERIVATIVE, HILL_MANDEL); // the sub-tensor-type does'nt matter
    double e11 = E[0][0]; // 1
    double e12 = E[0][1]; // 2
    double e22 = E[1][1]; // 3
    double e13 = E[0][2]; // 4
    double e23 = E[1][2]; // 5
    double e33 = E[2][2]; // 6

    // (6.72) in the diss of Sonja Lehmann. For the eps see Function::CalcPosDefDeterminant()!
    out[0]  = factor * (e33 - v);      // 1
    out[1]  = factor * (-2.0 * e23);   // 2
    out[2]  = factor * (e22 - v);      // 3
    out[3]  = factor * (-2.0*(e33-v)); // 4
    out[4]  = factor * (2.0*e23);      // 5
    out[5]  = factor * (2.0*e13);      // 6
    out[6]  = factor * (-2.0*e12);     // 7
    out[7]  = factor * (-2.0*e13);     // 8
    out[8]  = factor * (e11-v-eps);    // 9
    //out[8]  = factor * (e11-v);    // 9
    out[9]  = factor * (-2.0*(e22-v)); // 10
    out[10] = factor * (2.0*e12);      // 11
    out[11] = factor * (-2.0*(e11-v-eps)); // 12
    //out[11] = factor * (-2.0*(e11-v)); // 12
    break;
  }
  case DISTANCE:
  {
    // exact Hessian of c = ||Q-P|| w.r.t. [px py (pz) qx qy (qz)]: the [[Hv,-Hv],[-Hv,Hv]] block with
    // Hv = (I - dd^T/r^2)/r, d = Q-P, r = |d| (same as FeatureMappingDesign::Pill::HessLength). The
    // values match the (row,col) order of GetHessianSparsityPattern().
    Function::Local::Identifier& id = GetCurrentVirtualContext();
    int dim = domain->GetDim();
    int nc = 2 * dim; // number of node coordinates
    double H[6][6];
    if(dim == 2)
    {
      double px = id.GetElement(-1)->GetPlainDesignValue();
      double py = id.GetElement(0)->GetPlainDesignValue();
      double qx = id.GetElement(1)->GetPlainDesignValue();
      double qy = id.GetElement(2)->GetPlainDesignValue();
      double dx = qx - px, dy = qy - py;
      double r2 = dx*dx + dy*dy, r = std::sqrt(r2);
      assert(r > 0.0);
      double a = (1.0 - dx*dx/r2) / r;
      double b = -dx*dy / (r2*r);
      double c = (1.0 - dy*dy/r2) / r;
      const double H2[4][4] = {{ a,  b, -a, -b},
                               { b,  c, -b, -c},
                               {-a, -b,  a,  b},
                               {-b, -c,  b,  c}};
      for(int i = 0; i < 4; i++)
        std::copy(H2[i], H2[i]+4, H[i]);
    }
    else
    {
      double C[6]; // P and Q as in CalcDistance
      for(int i = 0; i < 6; i++)
        C[i] = id.GetElement(i-1)->GetPlainDesignValue();
      double d[3] = {C[3]-C[0], C[4]-C[1], C[5]-C[2]};
      double r2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
      double r = std::sqrt(r2);
      assert(r > 0.0);
      for(int i = 0; i < 6; i++)
        for(int j = 0; j < 6; j++)
        {
          int di = i < 3 ? i : i - 3;
          int dj = j < 3 ? j : j - 3;
          H[i][j] = (i/3 == j/3 ? 1.0 : -1.0) * ((di == dj ? 1.0 : 0.0) - d[di]*d[dj]/r2)/r;
        }
    }
    // fixed nodes drop out, matching GetHessianSparsityPattern()
    int pos[6];
    int nf = 0;
    for(int i = -1; i <= nc - 2; i++)
      if(OptIndexIfVariable(id.GetElement(i)) >= 0)
        pos[nf++] = i + 1;
    for(int i = 0, k = 0; i < nf; i++)
      for(int j = 0; j < nf; j++, k++)
        out[k] = factor * H[pos[i]][pos[j]];
    break;
  }
  default:
    assert(out.GetSize() == 0);
    break;
  }
}

double LocalCondition::CalcMeanAbsValue() const
{
  assert(local->local_values.GetSize() > 0);

  double sum = 0.0;

  // we cannot use GetValue() as we have no current_local_index
  for(double val : local->local_values)
    sum += std::abs(val);

  return sum / local->local_values.GetSize();
}

double LocalCondition::CalcMinMaxAbsValue() const
{
  assert(local->virtual_elem_map.GetSize() == local->local_values.GetSize());
  double minmax;
  if(bound_ == LOWER_BOUND)
    minmax = std::numeric_limits<double>::infinity();
  else
    minmax = 0.0;
  for(double val : local->local_values)
    if(bound_ == LOWER_BOUND)
      minmax = std::min(minmax, std::abs(val));
    else
      minmax = std::max(minmax, std::abs(val));
  return minmax;
}


int LocalCondition::CountInfeasibles() const
{
  assert(local->virtual_elem_map.GetSize() == local->local_values.GetSize());

  int cnt = 0;
  for(unsigned int i = 0; i < local->local_values.GetSize(); i++)
  {
    double v = local->local_values[i];
    double d = v - GetBoundValue();
    LOG_DBG2(conditions) << "LC:CI check f=" << ToString() << " i=" << i << " b=" << Condition::bound.ToString(bound_) << " v=" << v << " bv=" << GetBoundValue() << " d=" << d << " cnt=" << cnt;

    // upper_bound: 0 < 3 -> -3 < 0 (ok)   4 < 3 -> 1 < 0 (false)
    // lower_bound: 3 > 2 ->  1 > 0 (ok)   3 > 4 -> -1 > 0 (false)

    if((bound_ == Condition::EQUAL && std::abs(d) > 1e-5) ||
       (bound_ == Condition::LOWER_BOUND && d <= 1e-6) ||
       (bound_ == Condition::UPPER_BOUND && d >= -1e-6))
    {
      cnt++;
      LOG_DBG(conditions) << "LC:CI -> count f=" << ToString() << " i=" << i << " v=" << v << " bv=" << GetBoundValue() << " d=" << d << " cnt=" << cnt;
    }
  }

  return cnt;
}

double LocalCondition::GetValue() const
{
  if(IsLocal())
    return local->local_values[current_view_index_ - virtual_base_index_];
  else
    return value_;
}

void LocalCondition::SetValue(double val)
{
  if(IsLocal())
  {
    local->local_values[current_view_index_ - virtual_base_index_] = val;
    value_ = -1.0; // invalidated
  }
  else
    value_ = val; // only for Done() allowed!!
}


string LocalCondition::ToString(MultipleExcitation* me) const
{
  std::stringstream ss;

  ss << Condition::ToString();

   if(IsLocal())
      ss << " cvi=" << current_view_index_;

  return ss.str();
}


ConditionContainer::ConditionContainer()
{
  view = NULL;
  space_ = NULL;
  // set in Read()
}

ConditionContainer::~ConditionContainer()
{
  delete(view);

  for(unsigned int i = 0; i < active.GetSize(); i++) {
    delete active[i];
    active[i] = NULL;
  }

  for(unsigned int i = 0; i < observe.GetSize(); i++) {
    delete observe[i];
    observe[i] = NULL;
  }
}

void ConditionContainer::Read(ParamNodeList pn_list)
{
  // call only once
  assert(all.IsEmpty());

  // slope constraints need to be post processed in ErsatzMaterial
  bool displacement_constr;
  for(unsigned int i = 0; i < pn_list.GetSize(); i++)
  {
    PtrParamNode pn = pn_list[i];
    bool act = pn->Get("mode")->As<string>() == "constraint";

    // Add multiple displacement constraints using output displacement constraint on multiple numbered nodes
    PtrParamNode output;
    ParamNodeList elems;
    displacement_constr = false;
    if (pn->Has("output")) {
      output = pn->Get("output");
      if (output->Has("displacement")) {
        elems = output->GetList ("displacement");
        assert(elems.GetSize() == 1);
        PtrParamNode xml = elems[0];
       if (xml->Has("multiple_nodes")) {
          UInt end = xml->Get("multiple_nodes")->As<UInt>();
          displacement_constr = true;
          std::string entName = xml->Get("name")->As<std::string>();

          for (UInt j = 0; j < end; j++) {
            Condition::AddCondition(pn, act ? active : observe,j+1,entName);
          }
        }
      }
    }
    // General constraint (Non displacement constraint case)
    if (!displacement_constr)
      Condition::AddCondition(pn, act ? active : observe);
  }

  // process the virtual containers
  Refresh(); // Save for first time call
}

void ConditionContainer::Refresh()
{
  // process the virtual containers
  all.Resize(0); // when we do it again.
  all.Reserve(active.GetSize() + observe.GetSize());

  for(unsigned int i = 0; i < active.GetSize(); i++)
    all.Push_back(active[i]);

  for(unsigned int i = 0; i < observe.GetSize(); i++)
    all.Push_back(observe[i]);

  // set index
  for(unsigned int i = 0; i < all.GetSize(); i++)
    all[i]->index_ = i;

  if(view == NULL)
    view = new VirtualView(this);

  view->Refresh();
}



void ConditionContainer::PostProc(DesignSpace* space, DesignStructure* structure, MultipleExcitation* me, ErsatzMaterial* em)
{
  this->space_ = space;

  // the conditions have no space
  for(unsigned int i = 0; i < all.GetSize(); i++)
    if(all[i]->HasDenseJacobian())
      all[i]->SetDenseSparsityPattern(space);

  for(unsigned int i = 0; i < all.GetSize(); i++)
  {
    all[i]->PostProc(space, structure, em);
    all[i]->SetExcitation(me);
  }

  // check for uniqueness of the eigenvalue id
  if(Optimization::context->IsEigenvalue())
  {
    unsigned int max = Optimization::context->GetEigenFrequencyDriver()->GetNumSteps();

    StdVector<unsigned int> ids;

    for(unsigned int i = 0; i < all.GetSize(); i++)
    {
      if(all[i]->GetType() == Function::EIGENFREQUENCY)
      {
        unsigned int id = all[i]->GetEigenValueID();
        if(all[i]->ctxt->DoBuckling())
        {
          StdVector<unsigned int> order = all[i]->ctxt->GetBucklingDriver()->GetModeOrder();
          id = order[id];
        }

        assert(id > 0); // ensured by xml schema
        if(id > max)
          EXCEPTION("eigenvalue id 'ev'" << id << " larger than the " << max << " calculated eigenfrequencies");

        if(ids.Contains(id))
          EXCEPTION("eigenvalue id 'ev'" << id << " is not unique");

        ids.Push_back(id);
      }
    }
  }

  // for stress constraints with multiple excitation we insert additional stress constraints for
  // the specific excitations. This cannot be done in all.
  Condition::AddExcitationStressConstraints(active, me);
  Condition::AddExcitationStressConstraints(observe, me);


  // in the bloch mode optimization case we either have a constraint for every wave vector or search for the extremals
  Condition::AddBlochEigenConstraints(active, me);
  Condition::AddBlochEigenConstraints(observe, me);

  Refresh(); // inform about the news if the slopes created a lot of virtual objectives!
}

void ConditionContainer::ToInfo(PtrParamNode in)
{
  for(unsigned int i = 0; i < all.GetSize(); i++)
    all[i]->ToInfo(in->Get("constraint", ParamNode::APPEND));
}


Condition* ConditionContainer::Get(Condition::Type type, DesignElement::Type design, Condition::Bound bound, bool throw_exception)
{
  assert(design != DesignElement::NO_TYPE);

  for(unsigned int i = 0; i < active.GetSize(); i++)
    if(active[i]->GetType() == type && active[i]->design_ == design && active[i]->bound_ == bound)
      return active[i]; // shall be unique

  if(throw_exception)
    throw Exception("have no active constraint " + Condition::type.ToString(type) + " with design " + DesignElement::type.ToString(design) + " and bound " + Condition::bound.ToString(bound));

  return NULL;
}

Condition* ConditionContainer::Get(const std::string name, bool throw_exception)
{
  for(Condition* g : all)
    if(g->ToString() == name)
      return g;

  if(throw_exception)
    throw Exception("condition '" + name + "' not known within our " + std::to_string(all.GetSize()) + " constraints/observes");
  else
    return nullptr;
}

StdVector<Condition*> ConditionContainer::GetList(Condition::Type type, DesignElement::Type design, bool only_active, Function::Access access)
{
  StdVector<Condition*> result;

  for(unsigned int i = 0, n = active.GetSize() + (only_active ? 0 : observe.GetSize()); i < n; i++)
  {
    assert(!(only_active && i >= active.GetSize()));
    Condition* g = i < active.GetSize() ? active[i] : observe[i-active.GetSize()];

    if(g->GetType() != type)
      continue;

    if(design != DesignElement::NO_TYPE && g->design_ != design)
      continue;

    if(access != Function::NO_ACCESS && g->GetAccess() != access)
      continue;

    result.Push_back(g);
  }
  return result;
}


bool ConditionContainer::HasUniqueBounds(const StdVector<Condition*>& list)
{
  bool lower = false;
  bool upper = false;
  bool equal = false;

  for(unsigned int i = 0;i < list.GetSize();i++) {
    switch (list[i]->GetBound()) {
     case Condition::LOWER_BOUND:
       if (lower) {
         return false;
       } else {
         lower = true;
       }
       break;
     case Condition::UPPER_BOUND:
       if (upper) {
         return false;
       } else {
         upper = true;
       }
       break;
     case Condition::EQUAL:
       if (equal) {
         return false;
       } else {
         equal = true;
       }
       break;
    }
  }
  return true;
}

bool ConditionContainer::RequiresBoundForUniqueness(const Condition* g) {
  const StdVector<Condition*> list = GetList(g->GetType(),g->GetDesignType(),false, g->GetAccess());
  return list.GetSize() > 1 && HasUniqueBounds(list);
}


bool ConditionContainer::IsAllStateDependent() const
{
  for(auto c : active)
    if(c->IsStateDependent())
      return true;

  return false;
}


bool ConditionContainer::Has(Condition::Type type, DesignElement::Type design, bool only_active)
{
  // be save and check for uniqueness!
  StdVector<Condition*> list = GetList(type, design, only_active);

  return !list.IsEmpty();
}


Condition* ConditionContainer::Get(Condition::Type type, DesignElement::Type design, bool only_active, bool throw_exception)
{
  // be save and check for uniqueness!
  StdVector<Condition*> list = GetList(type, design, only_active);

  if(list.GetSize() == 0)
  {
    if(throw_exception)
      throw Exception("no active constraint '" + Condition::type.ToString(type) + "' found");
    else
      return NULL;
  }

  if(list.GetSize() > 1 && !HasUniqueBounds(list) && throw_exception)
    throw Exception("constraint " + Condition::type.ToString(type) + " is not unique");

  return list[0];
}

void ConditionContainer::PushBackHistory()
{
  for(Condition* g : all)
    g->history.Push_back(g->IsLocal() ? 0.0 : g->value_);
}


ConditionContainer::VirtualView::VirtualView(ConditionContainer* constraints)
{
  container_ = constraints;
  Refresh();
}

void ConditionContainer::VirtualView::Refresh()
{
  // find the indices of LocalConditions, they need to be sorted.
  std::list<unsigned int> tmp;

  // search also for observe conditions!
  for(int i = 0; i < 100; i++) // assume to get all valid function types
  {
    Function::Type cand = (Function::Type) i;
    if(Function::IsLocal(cand))
    {
      const StdVector<Condition*>& c = container_->GetList(cand, DesignElement::NO_TYPE, false);
      for(UInt i=0; i<c.GetSize(); ++i) {
        tmp.push_back(c[i]->GetIndex());
        LOG_DBG2(conditions) << "CC:VV:R: add " << Function::type.ToString(cand) << " c=" << c.GetSize() << " idx=" << c[i]->GetIndex();
      }
    }
  }
  tmp.sort();

  // copy sorted
  local_cond_index_.Resize(0);
  for(std::list<unsigned int>::iterator it = tmp.begin(); it != tmp.end(); ++it)
    local_cond_index_.Push_back(*it);

  // determine the virtual sizes for the container
  virtual_active_size_ = container_->active.GetSize();
  virtual_total_size_ = container_->all.GetSize();

  for(unsigned int i = 0; i < local_cond_index_.GetSize(); i++)
  {
    LocalCondition* lc = dynamic_cast<LocalCondition*>(container_->all[local_cond_index_[i]]);

    // refresh is called multiple times, maybe we are too early
    if(lc == NULL || lc->GetLocal() == NULL) continue;

    // replace the global slope by many local slopes -> if it is initialized!
    if(lc->IsActive() && lc->GetConstraintSize() > 0)
      virtual_active_size_ += lc->GetConstraintSize() -1; // replace means remove ourselves

    if(lc->GetConstraintSize() > 0)
      virtual_total_size_ += lc->GetConstraintSize() -1;

    LOG_DBG2(conditions) << "CC:VV:R i=" << i << " lci=" << local_cond_index_[i] << " lc=" << lc->ToString() << " vas=" << virtual_active_size_ << " vts=" << virtual_total_size_;
  }

  // set the virtual base indices
  int curr = 0;
  for(unsigned int i = 0; i < container_->all.GetSize(); i++)
  {
    Condition* g = container_->all[i];
    g->virtual_base_index_ = curr;
    if(g->IsLocalCondition() && g->GetLocal() != NULL) // does not need to be initialized yet
      curr += std::max((int) dynamic_cast<LocalCondition*>(g)->GetConstraintSize(), 1);
    else
      curr++;
    LOG_DBG2(conditions) << "CC:VV:R g=" << Condition::type.ToString(g->GetType()) << " vbi=" << g->virtual_base_index_ << " new curr=" << curr;
  }
  assert(curr == virtual_total_size_);
}

Condition* ConditionContainer::VirtualView::Get(int view_index)
{
  StdVector<Condition*>& all = container_->all;;
  assert(all.GetSize() > 0);

  // traverse the conditions, if we are above virtual_base_index we gone one too far
  Condition* g = NULL;
  for(unsigned int i = 0; g == NULL && i < all.GetSize()-1; i++)
    if(all[i+1]->virtual_base_index_ > view_index) // we are right if the next is too far
      g = all[i];
  if(g == NULL) g = all.Last();
  assert(g->virtual_base_index_ <= view_index);

  if(g->IsLocalCondition())
  {
    dynamic_cast<LocalCondition*>(g)->SetCurrentViewIndex(view_index);
  }

  return g;
}

void ConditionContainer::VirtualView::Done()
{
  for(unsigned int g = 0; g < local_cond_index_.GetSize(); g++) // no local conditions, nothing to do
  {
    LocalCondition* lc = dynamic_cast<LocalCondition*>(container_->all[local_cond_index_[g]]);
    lc->SetCurrentViewIndex(-1); // reset to global mode

    // shall we give the values as special result?
    DesignElement::ValueSpecifier vs = DesignElement::MAX_OSCILLATION; // overwrite if necessary
    if(lc->GetType() == Function::SLOPE)        vs = DesignElement::MAX_SLOPE;
    if(lc->GetType() == Function::MOLE)         vs = DesignElement::MAX_MOLE;
    if(lc->GetType() == Function::JUMP)         vs = DesignElement::MAX_JUMP;

    int idx = container_->space_->GetSpecialResultIndex(DesignElement::DEFAULT, vs);
    if(idx >= 0)
    {
      // we add up the max value and not elements have a slope constraint, therefore reset
      for(unsigned int e = 0; e < lc->elements.GetSize(); e++)
        lc->elements[e]->specialResult[idx] = 0.0; // initialize

      StdVector<Function::Local::Identifier>& vem = lc->GetLocal()->virtual_elem_map;

      for(unsigned int i = 0; i < vem.GetSize(); i++)
      {
        Function::Local::Identifier& id = vem[i];
        assert(lc->GetType() !=  Function::SHAPE_INF);
        DesignElement* de =  dynamic_cast<DesignElement*>(id.element);
        double sv = id.EvalFunction(lc->local);

        // in checkerboard we must not use abs
        double corr = lc->GetType() == Function::OSCILLATION ? sv : std::abs(sv);


        de->specialResult[idx] = std::max(de->specialResult[idx], corr);
      }
    }
  }
}

int ConditionContainer::VirtualView::CalcNumberOfJacobianNonZeros()
{
  int nnz = 0;
  for(int c = 0; c < GetNumberOfActiveConstraints(); c++)
    nnz += Get(c)->GetSparsityPatternSize();
  Done(); // reset a potential slope constraint back to global mode, like EvalGradConstraints()
  return nnz;
}



