#include <assert.h>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>

#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/Logging/log.hpp"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Domain/BCs.hh"
#include "Domain/Domain.hh"
#include "Domain/ElemMapping/Elem.hh"
#include "Domain/ElemMapping/EntityLists.hh"
#include "Domain/Mesh/Grid.hh"
#include "Domain/ElemMapping/SurfElem.hh"
#include "Driver/Assemble.hh"
#include "Driver/BaseDriver.hh"
#include "Driver/FormsContexts.hh"
#include "Driver/HarmonicDriver.hh"
#include "Driver/SolveSteps/BaseSolveStep.hh"
#include "Driver/TransientDriver.hh"
#include "Driver/EigenFrequencyDriver.hh"
#include "FeBasis/BaseFE.hh"
#include "FeBasis/FeSpace.hh"
#include "Forms/LinForms/LinearForm.hh"
#include "Forms/LinForms/SingleEntryInt.hh"
#include "Domain/CoefFunction/CoefXpr.hh"
#include "General/Enum.hh"
#include "General/defs.hh"
#include "General/Environment.hh"
#include "General/Exception.hh"
#include "MatVec/SingleVector.hh"
#include "MatVec/exprt/xpr1.hh"
#include "MatVec/exprt/xpr2.hh"
#include "MatVec/Matrix.hh"
#include "MatVec/Vector.hh"
#include "MatVec/SBM_Matrix.hh"
#include "Materials/MechanicMaterial.hh"
#include "Optimization/Condition.hh"
#include "Optimization/Design/DensityFile.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Design/AuxDesign.hh"
#include "Optimization/Design/DesignStructure.hh"
#include "Optimization/ErsatzMaterial.hh"
#include "Optimization/Excitation.hh"
#include "Optimization/Function.hh"
#include "Optimization/Objective.hh"
#include "Optimization/Optimization.hh"
#include "Optimization/Optimizer/BaseOptimizer.hh"
#include "Optimization/OptimizationMaterial.hh"
#include "Optimization/SIMP.hh"
#include "Optimization/StressConstraint.hh"
#include "Optimization/TransferFunction.hh"
#include "Optimization/Context.hh"
#include "PDE/SinglePDE.hh"
#include "PDE/StdPDE.hh"
#include "PDE/BasePDE.hh"
#include "PDE/MechPDE.hh"
#include "PDE/HeatPDE.hh"
#include "PDE/LatticeBoltzmannPDE.hh"
#include "Utils/Point.hh"
#include "Utils/StdVector.hh"
#include "Utils/mathParser/mathParser.hh"
#include "Utils/tools.hh"
#include "Utils/Timer.hh"

namespace CoupledField {
class BaseMaterial;
class CoordSystem;
class DenseMatrix;
struct ResultInfo;
}  // namespace CoupledField

using namespace std;

DECLARE_LOG(conditions)
DEFINE_LOG(conditions, "conditions")

DECLARE_LOG(em)
DEFINE_LOG(em, "ersatzMaterial")


namespace CoupledField {

Enum<ErsatzMaterial::Method> ErsatzMaterial::method;

ErsatzMaterial::ErsatzMaterial() :
  Optimization(),
  dim(grid->GetDim())
{
  /** We store here the solution */
  volume_fraction_ = 0.0;
  structure_ = NULL;
  densityFile = NULL;
  bitensor_ = false;
  trackingFunc_ = NULL;

  interfaceDrivenGradCalc_ = false;

  pn = domain->GetParamRoot()->Get("optimization/ersatzMaterial");

  method_ = method.Parse(pn->Get("method")->As<std::string>());

  // we set the calc_u1ku2_timer_ only for non-regular meshes but then as sub-timer
  calc_u1ku2_timer_ = grid->IsGridRegular() ? boost::shared_ptr<Timer>() : boost::shared_ptr<Timer>(new Timer("calc_U1KU2", true));

  // region stuff - we might have the attribute region or a list in region but not both or none
  if(!pn->Has("region") && !pn->Has("regions") && (method_ != SHAPE_OPT && method_ != SHAPE_PARAM_MAT))
    throw Exception("give a region as 'ersatzMaterial' attribute or as 'regions' list");

  if(pn->Has("region") && pn->Has("regions") && (method_ != SHAPE_OPT && method_ != SHAPE_PARAM_MAT))
    throw Exception("you may not give regions via 'ersatzMaterial' attribute or as 'regions' list concurrently");

  ParamNodeList region_list;
  if(pn->Has("regions"))
    region_list = pn->Get("regions")->GetList("region");
  else
    region_list = pn->GetList("region");

  if(region_list.IsEmpty() && (method_ != SHAPE_OPT && method_ != SHAPE_PARAM_MAT))
    EXCEPTION("no region given!");

  StdVector<RegionIdType> regions;

  if(method_ != SHAPE_OPT){ // we want no regions here, if only shape is optimized, even if there would be any in xml
    for(unsigned int i=0; i < region_list.GetSize(); i++){
      // we are compatible with the region attribute and unbounded region elements
      std::string reg = region_list[i]->Has("name") ? region_list[i]->Get("name")->As<std::string>() : region_list[i]->As<std::string>();
      std::string bimat = region_list[i]->Has("bimaterial") ? region_list[i]->Get("bimaterial")->As<std::string>() : "";
      if(!grid->GetRegion().IsValid(reg))
        throw Exception("region given in ersatzMaterial is invalid");
      if (std::count (regions.Begin(), regions.End(), grid->GetRegion().Parse(reg)) > 0)
        throw Exception("region "+ reg + " is given multiple times in ErsatzMaterial");
      regions.Push_back(grid->GetRegion().Parse(reg));
    }
  }

  design = DesignSpace::CreateInstance(regions, pn, method_);

  // the L-mesh of the stress constraint benchmark is meshed by gid with different positions of
  // element nodes, such that one cannot use the same element matrix, even if the grid is regular
  // therefore the attribute enforce_unstructured

  // optionally write the densities to a density.xml file
  if(pn->Has("export"))
  {
    ParamNodeList design_list = pn->GetList("design");
    ParamNodeList transfer_list = pn->GetList("transferFunction");
    densityFile = new DensityFile(design, pn->Get("export"), design_list, transfer_list, pn->Get("filters", ParamNode::PASS));
  }

  // check our constraints, the shall have only valid designs
  for(unsigned int i = 0; i < constraints.all.GetSize();  i++)
  {
    Condition* g = constraints.all[i];
    DesignElement::Type dt = g->GetDesignType();
    if(dt == DesignElement::UNITY && (method_ == SHAPE_OPT || method_ == SHAPE_PARAM_MAT))
      continue;

    if((   dt == DesignElement::MECH_TRACE || dt == DesignElement::MECH_ALL
        || dt == DesignElement::DIELEC_TRACE || dt == DesignElement::DIELEC_ALL
        || dt == DesignElement::PIEZO_ALL || dt == DesignElement::ALL_DESIGNS)
       && (method_ == PARAM_MAT || method_ == SHAPE_PARAM_MAT))
      continue;

    if(dt != DesignElement::DEFAULT && design->FindDesign(g->GetDesignType(), false) == -1)
      throw Exception("constraint '" + Function::type.ToString(g->GetType()) + "' operates on invalid design variable"); // ToString() may trigger too much within constructor

  }

  // give the domain this data, s.th. the ersatz material approach is applied
  domain->SetDesign(design);
}

ErsatzMaterial::~ErsatzMaterial()
{
  // if write to file close the xml envelope and the file
  if(densityFile != NULL) { delete densityFile; densityFile = NULL; }

  delete structure_;

  // "remove" the ersatzMaterial (=data) from the domain
  domain->SetDesign(NULL);
}

void ErsatzMaterial::PostInit()
{
  // updates context which we need for the filters (pde)
  Optimization::PostInit();

  // from the filters we detect robustness which we need for multiple excitations
  if(pn->Has("filters"))
  {
    ParamNodeList list = pn->Get("filters")->GetList("filter");
    // this is save for design=polarization
    for(unsigned int i = 0; i < list.GetSize(); i++)
    {
      if(structure_ == NULL)
        structure_ = new DesignStructure(this);
      structure_->SetFilter(list[i], this->optInfoNode);
    }
  }

  // check for multiple load cases (might be frequencies)
  me->InitializeMultipleExcitations(this, &manager);
  for(unsigned int i = 0; i < manager.context.GetSize(); i++)
    me->PrepareMultipleExcitations(this, &(manager.context[i]));
  me->FinalizeMultipleExcitations(this, &manager, optimizer_ == EVALUATE_INITIAL_DESIGN);
  me->excitations.First().Apply(false); // this sets the first excitation but does not switch context. This is done below

  // add optimization results to the pde
  for(unsigned int c = 0; c < manager.context.GetSize(); c++) {
    assert(manager.context[c].pdes.size() == 1); // extend!
    design->AppendOptimizationResults(manager.context[c].pde, !context->DoMultiSequence()); // don't warn in multi-sequence case
  }

  // might be constructed in SIMP::PostInit() or ParamMat::PostInit()
  if(structure_ == NULL)
    structure_ = new DesignStructure(this);

  // post init slope constraints when the design is there
  constraints.PostProc(design, structure_, GetMultipleExcitation(), this);
  // same for the objectives
  objectives.PostProc(design, structure_, GetMultipleExcitation());

  // the constraints size is only now known and the shapeDesign constructor is finished -> PostInit design
  design->PostInit(objectives.data.GetSize(), constraints.all.GetSize());

  unsigned int total = objectives.data.GetSize() + constraints.active.GetSize() + constraints.observe.GetSize();

  for(unsigned int i = 0; i < total; i++)
  {
    Function* f = i < objectives.data.GetSize() ? dynamic_cast<Function*>(objectives.data[i]) : dynamic_cast<Function*>(constraints.all[i - objectives.data.GetSize()]);

    std::string func = "'" + f->type.ToString(f->GetType()) + "'";

    // checks
    switch(f->GetType())
    {
    case Function::COMPLIANCE:
    case Function::OUTPUT: // it would work but is saver not to allow
    case Function::SQUARED_OUTPUT:
    case Function::TRACKING:
    case Function::HOM_TENSOR:
    case Function::HOM_TRACKING:
    case Function::HOM_FROBENIUS_PRODUCT:
    case Function::POISSONS_RATIO:
    case Function::YOUNGS_MODULUS:
    case Function::YOUNGS_MODULUS_E1:
    case Function::YOUNGS_MODULUS_E2:
    case Function::ELEC_ENERGY: // it simply does not work yet in the harmonics
      if(f->ctxt->IsComplex())
        throw Exception(func + " is only for static state problems");
      break;

    case Function::DYNAMIC_OUTPUT:
    case Function::GLOBAL_DYNAMIC_COMPLIANCE:
      if(!f->ctxt->IsComplex())
        throw Exception(func + " is only for harmonic state problems");
      break;

    case Function::EIGENFREQUENCY:
      if(!f->ctxt->IsEigenvalue())
        throw Exception(func + " is only for eigenvalue state problems");
      break;

    default:
      break;
    }

    // actions
    switch(f->GetType())
    {
    case Function::OUTPUT:
    case Function::SQUARED_OUTPUT:
    case Function::DYNAMIC_OUTPUT:
    case Function::ABS_OUTPUT:
    {
      if(!f->pn->Has("output"))
        throw Exception("element 'output' is mandatory for function 'output'");
      PtrParamNode output = f->pn->Get("output");

      StdVector<shared_ptr<EntityList> > ent;
      StdVector<PtrCoefFct > coef;
      bool geo = false;
      assert(!context->DoMultiSequence()); // the pdes are not know yet!
      SinglePDE* pde = context->pde;

      if(output->Has("displacement"))
        pde->ReadRhsExcitation("displacement", pde->GetDofNames(MECH_DISPLACEMENT), ResultInfo::VECTOR, context->IsComplex(), ent, coef, geo, output);

      if(output->Has("elecPotential"))
        pde->ReadRhsExcitation("elecPotential", pde->GetDofNames(ELEC_POTENTIAL), ResultInfo::SCALAR, context->IsComplex(), ent, coef, geo, output);

      if(output->Has("acoustic"))
        assert(false);

        //domain->GetSinglePDE("acoustic")->ReadLoads(output->GetList("acoustic"), f->output_nodes);

      // we store the loads in forms of linear forms
      for(unsigned int i = 0; i < ent.GetSize(); ++i )
      {
        assert(ent[i]->GetType() == EntityList::NODE_LIST);

        if(ent[i]->GetSize() > 1 ) { // MechPDE.cc -> "force"
          Global::ComplexPart part = context->IsComplex() ? Global::COMPLEX : Global::REAL;
          coef[i] = CoefFunction::Generate(domain->GetMathParser(), part, CoefXprVecScalOp(domain->GetMathParser(), coef[i], boost::lexical_cast<std::string>(ent[i]->GetSize()), CoefXpr::OP_DIV));
        }

        LinearForm* lin = new SingleEntryInt(coef[i]);
        lin->SetName("NodalForceInt");
        LinearFormContext *ctx = new LinearFormContext( lin );
        ctx->SetEntities( ent[i] );
        ctx->SetFeFunction(pde->GetFeFunction(pde->GetNativeSolutionType()));

        LOG_DBG2(em) << "PI: o:" << f->IsObjective() << " add output form " << ent[i]->GetName() << " size=" << ent[i]->GetSize();

        f->output_forms.Push_back(ctx);
      }
      LOG_DBG2(em) << "PI: size of output_forms: " <<f->output_forms.GetSize() ;
      if(f->output_forms.GetSize() == 0)
        throw Exception("no output optimization targets given");
      break;
    }
    // we do energy flux in realtime
    default:
      break;
    }
  }

  // if loadErsatzMaterial is used with optimization specifying a starting point,
  // we have to load it here, before scaling is done.
  if(DensityFile::NeedLoadErsatzMaterial())
    DensityFile::ReadErsatzMaterial(design);

  // plausibility check for homogenization
  assert(manager.context.GetSize() >= 1);
  assert(manager.IsInitialized());

  for(unsigned int i = 0; i < objectives.data.GetSize(); i++)
    if(objectives.data[i]->IsHomogenization())
      objectives.data[i]->ctxt->homogenization = true;

  for(unsigned int i = 0; i < constraints.all.GetSize(); i++)
    if(constraints.all[i]->IsHomogenization())
      constraints.all[i]->ctxt->homogenization = true;

  for(unsigned int i = 0; i < manager.context.GetSize(); i++)
    manager.context[i].infoNode->Get("homogenization")->SetValue(manager.context[i].homogenization);

  for(unsigned int i = 0; i < manager.context.GetSize(); i++)
  {
    const Context& c = manager.context[i];

    if(c.homogenization && (!me->IsEnabled(c.sequence) || !(me->DoHomogenization())))
      throw Exception("A homogenization objective/constraint is set but no homogenization test strain excitation");
    if(me->IsEnabled(c.sequence) && me->GetSequence() == c.sequence && me->DoHomogenization() && !c.homogenization)
      throw Exception("No homogenization objective/constraint for homogenization test strain excitation");
  }

  // for transformations we might have more than only one tensor
  Context* ctxt = manager.GetHomogenization();
  if(ctxt) {
    homogenizedTensor.Resize(me->GetNumberMeta(ctxt, true));
    for(unsigned int i = 0; i < homogenizedTensor.GetSize(); i++)
      homogenizedTensor[i].Resize(dim == 2 ? 3 : 6, dim == 2 ? 3 : 6);
  }

  if(manager.any().eigenvalue && optParamNode->Has("eigenvalue/sort"))
  {
    for(unsigned int i = 0; i < constraints.all.GetSize(); i++)
    {
      int idx = constraints.all[i]->GetEigenValueID();
      if(idx > 0)
        log.AddToHeader("mode_" + boost::lexical_cast<std::string>(idx));
    }
  }

  // make basic logging
  design->ToInfo(this);

  if(calc_u1ku2_timer_)
    optInfoNode->Get(ParamNode::SUMMARY)->Get("calcUKU/timer")->SetValue(calc_u1ku2_timer_);
}


void ErsatzMaterial::StoreResults(double step_val)
{
  CommitMode cm = commitMode_;

  double real_step = step_val != -1 ? step_val : currentIteration;

  if(cm == EACH_FORWARD)
  {
    for(unsigned int e = 0; e < me->excitations.GetSize(); e++)
    {
      Excitation& ex = me->excitations[e];
      Context& ctxt = manager.GetContext(&ex);
      // in case of eigenvalues we have for each excitation many modes but write only the first one here
      // in that case we need mode number 0 instead of the default -1
      int mode = ctxt.IsEigenvalue() ? 0 : -1;
      forward.Get(ex, NULL, mode)->Write(ctxt.pde); // forward is function NULL
      // call real implementation in Optimization. sum up in excitation fractions up to smaller 0.5
      Optimization::StoreResults(real_step + (0.5 / me->excitations.GetSize()) * e);
    }
  }

  if(cm == FORWARD || cm == BOTH) {
//    assert(!context->DoMultiSequence()); // extend!
    forward.WriteAverage(context->pde, context->sequence); // func = NULL
    Optimization::StoreResults(step_val);
  }

  // over all functions, mostly only one or none
  StdVector<const Function*> funcs = adjoint.GetFunctions();
  for(unsigned int fi = 0; fi < funcs.GetSize(); fi++)
  {
    LOG_DBG(em) << "StoreResults(" << step_val << ") rs=" << real_step << " adjoint_function=" << funcs[fi]->ToString();
    if(cm == EACH_ADJOINT)
    {
      for(unsigned int e = 0; me->excitations.GetSize(); e++)
      {
        assert(me->excitations[e].sequence == context->sequence); // will fail for multiple sequence!
        adjoint.Get(me->excitations[e], funcs[fi])->Write(context->pde);
        // call real implementation in Optimization. sum up in excitation fractions up to smaller 0.5
        double index = (me->excitations.GetSize() * funcs.GetSize()) * (fi * funcs.GetSize()) * e;
        Optimization::StoreResults(real_step + 0.5 + (0.5 / index));
      }
    }

    if(cm == ADJOINT || cm == BOTH) {
      // sum up if there are more excitations
      assert(!context->DoMultiSequence()); // extend!)
      adjoint.WriteAverage(context->pde, context->sequence, funcs[fi]); // time step??
      Optimization::StoreResults(real_step + 0.5 + (0.5 / funcs.GetSize()) * fi);
    }
  }
}


void ErsatzMaterial::LogFileLine(std::ofstream* out, PtrParamNode iteration)
{
  // first the basic stuff to have iteration and such stuff first
  Optimization::LogFileLine(out, iteration);

  // now the stuff prepared in ErsatzMaterial::CommitIteration()
  // in case of bloch_info
  for(unsigned int i = 0; i < log.bloch_info.GetSize(); i++)
    iteration->Get(boost::get<0>(log.bloch_info[i]))->SetValue(boost::get<1>(log.bloch_info[i]));

  // add our multiple excitation stuff here (only in info.xml, this would be to complex for dat
  if(me->IsEnabled())
  {
    for(unsigned int e = 0; e < me->excitations.GetSize(); e++)
    {
      Excitation& excite = me->excitations[e];
      PtrParamNode info = iteration->Get("excitation", ParamNode::APPEND);
      info->Get("index")->SetValue(excite.index);
      info->Get("objective")->SetValue(excite.cost);
      info->Get("objective_weight")->SetValue(excite.normalized_weight);
      for(unsigned int c = 0; c < constraints.all.GetSize(); c++)
      {
        Condition* g = constraints.all[c];
        if(g->IsExcitationSensitive() && g->DoEvaluate(&excite))
        {
          info->Get(g->ToString())->SetValue(g->GetValue());
          if(g->GetType() == Function::EIGENFREQUENCY && g->GetExcitation()->DoBloch() && !g->DoFullBloch()) {
            string label = "ef_" + lexical_cast<string>(g->GetEigenValueID()) + "_wv";
            info->Get(label)->SetValue(g->bloch.col);
          }
        }
      }
    }
  }

  for(unsigned int ci = 0; ci < manager.context.GetSize(); ci++)
  {
    Context* ctxt = &(manager.context[ci]);

    if(ctxt->homogenization)
    {
      for(unsigned int t = 0; t < homogenizedTensor.GetSize(); t++)
      {
        PtrParamNode in = iteration->Get("homogenizedTensor", ParamNode::APPEND);

        // assert(!(context->DoMultiSequence() && me->DoMetaExcitation(ctxt))); // check the base_index below!
        if(me->DoMetaExcitation(ctxt))
          in->Get("case")->SetValue(ctxt->GetExcitation(0, t)->GetMetaLabel());

        Matrix<double>& ht = homogenizedTensor[t];

        in->Get("norm_L2")->SetValue(ht.NormL2());
        in->Get("trace")->SetValue(ht.Trace());

        PtrParamNode iso = in->Get("isotropy");
        StdVector<std::pair<string, double> > isop = MechanicMaterial::CalcIsotropicProperties(ht, ctxt->stt);
        for(unsigned int p = 0; p < isop.GetSize(); p++)
          iso->Get(isop[p].first)->SetValue(isop[p].second);

        PtrParamNode orth = in->Get("orthotropy");
        // for the orthotropic case we need the design. This might be excitation dependent on the robust case
        assert(me->DoMetaExcitation(ctxt) || (ctxt->excitations.GetSize() == 3 || ctxt->excitations.GetSize() == 6)); // no robust!
        Excitation* ex = ctxt->GetExcitation(0, t);
        LOG_DBG2(em) << "CI hom t=" << t << " ex=" << ex->GetFullLabel() << " ht=" << ht.ToString();
        StdVector<std::pair<string, double> > ortho = GetOrthotropeProperties(ht, ex);
        for(unsigned int p = 0; p < ortho.GetSize(); p++)
          orth->Get(ortho[p].first)->SetValue(ortho[p].second);

        LOG_DBG(em) << "CI t=" << t << " ortho:" << ortho[0].first << "=" << ortho[0].second << " ht=" << ht.ToString();

        in->Get("tensor")->SetValue(ht);
      }
    }
  }
}

PtrParamNode ErsatzMaterial::CommitIteration()
{

  // to this before CommitIteration such that we can store the info in log.bloch_info. Later we use it to add it to the ParamNode
  // in case we do bloch and have eigenvalue with bloch=full (alpha+/-slack formulation) we additionally print here min max frequencies
  StdVector<Condition*> ev = constraints.GetList(Condition::EIGENFREQUENCY);
  for(int c = 0; c < (int) ev.GetSize()-1; c++) // allow always a next
  {
    Condition* g = ev[c];
    Condition* n = ev[c+1]; // see loop!
    assert(g->GetType() == Condition::EIGENFREQUENCY);

    if(g->GetExcitation()->DoBloch() && g->DoFullBloch() && g->GetBound() == Condition::UPPER_BOUND && n->GetBound() == Condition::LOWER_BOUND)
    {
      assert(n->GetExcitation()->DoBloch() && n->DoFullBloch());
      //assert(n->GetEigenValueID() == g->GetEigenValueID() + 1); // who knows what happens else?!
      const Matrix<double> mat = forward.CollectBlochEigenfrequencies(&(manager.GetContext(g->GetExcitation())));
      double lower, upper; // see ErsatzMaterial::CalcEigenFrequency()
      SearchMinMax(mat, (unsigned int) g->GetEigenValueID()-1, false, &lower);
      SearchMinMax(mat, (unsigned int) n->GetEigenValueID()-1, true, &upper);

      // replace the key, we have only "bandgap" in log
      //iter->Get("bandgap_" + lexical_cast<string>(g->GetEigenValueID()) + "_" + lexical_cast<string>(n->GetEigenValueID()))->SetValue(upper - lower);
      boost::get<0>(log.bloch_info.First()) ="bandgap_" + lexical_cast<string>(g->GetEigenValueID()) + "_" + lexical_cast<string>(n->GetEigenValueID());
      boost::get<1>(log.bloch_info.First()) = upper - lower;

      LOG_DBG(em) << "CI g=" << g->ToString() << "/" << g->GetEigenValueID() << " n=" << n->ToString() << "/" << n->GetEigenValueID();
      break; // assume only one lower/upper constraint gap
    }
  }

  // now add the min/max for the ev such we can analyse possible nonsmoothneses. However the gap is the information to search manually
  std::map<unsigned int, bool> ev_done; // what we did
  for(unsigned int c = 0; c < ev.GetSize(); c++)
  {
    Condition* g = ev[c];

    if(g->GetExcitation()->DoBloch() && g->DoFullBloch())
    {
      // we have the ev-constraint for every wave vector as we are full
      if(!ev_done[g->GetEigenValueID()]) // bool is false by default
      {
         const Matrix<double> mat = forward.CollectBlochEigenfrequencies(&(manager.GetContext(g->GetExcitation())));
         double freq; // see ErsatzMaterial::CalcEigenFrequency()
         SearchMinMax(mat, (unsigned int) g->GetEigenValueID()-1, g->GetBound() == Condition::LOWER_BOUND, &freq);

         //iter->Get(Condition::type.ToString(g->GetType()) + "_" + lexical_cast<string>(g->GetEigenValueID())+ (g->GetBound() == Condition::LOWER_BOUND ? "_min" : "_max"))->SetValue(freq);
         boost::get<1>(log.bloch_info[c+1]) = freq; // first entry is gap
         ev_done[g->GetEigenValueID()] = true;
         LOG_DBG(em) << "CI g=" << g->ToString() << " evid=" << g->GetEigenValueID() << " b=" << g->GetBound() << " mm=" << freq;
      }
    }
  }

  // will write the cfs results and the log file using possibly set log.bloch_info
  // by calling virtual LogFileLine()
  PtrParamNode iter = Optimization::CommitIteration();



  // write the current info file, if the writing frequency is not too high.
  domain->GetInfoRoot()->ToFile();

  if(densityFile != NULL)
    densityFile->SetAndWriteCurrent(currentIteration - 1); // already written in DesignSpace::ReadDesignFromExtern()

  return iter;
}

  StdVector<std::pair<string,double> > ErsatzMaterial::GetOrthotropeProperties(const Matrix<double>& tensor, Excitation* ex)
  {
    if(design->regions.GetSize() != 1)
    {
      StdVector<std::pair<string, double> > result;
      return result; // empty result
    }
    else
    {
      LOG_DBG2(em) << "GOP tensor=" << tensor.ToString();
      assert(ex != NULL);
      Context& ctxt = manager.GetContext(ex);
      assert(ex->sequence == ctxt.sequence);
      ex->Apply(false); // we read the design. When we do robust, this must match the filter associated to the tensor

      BaseMaterial* bm = NULL;
      // this happens when doing shape optimization with homTracking!
      // we then have no design region and need to skip GetForm
      if(design->GetRegionId() != -1)
        bm = ctxt.pde->GetMaterialData()[design->GetRegionId()];
      Objective vf(Function::VOLUME, 0.0, Function::PHYSICAL); // physical!
      assert(design->GetRegionIds().GetSize() ==1);
      vf.SetElements(design, design->GetRegionId());
      double vol = CalcVolume(&vf, NULL, false, true);
      StdVector<std::pair<string, double> > ortho = MechanicMaterial::CalcOrthotropeProperties(tensor, bm, ctxt.stt, vol);
      return ortho;
    }
  }

  string ErsatzMaterial::GetIterationFrequency()
  {
    if (!manager.any().harmonic)
      return "";

    // make clear, when doing *real* multiple excitations, that this is no single
    // frequency result
    if (me->IsEnabled() && me->excitations.GetSize() > 1)
      return "(mult)";

    // search frequency, assume only one context has a frequency
    for(unsigned int i = 0; i < manager.context.GetSize(); i++)
    {
      Context& ctxt = manager.context[i];
      if(ctxt.IsHarmonic())
      {
        double frequency = ctxt.GetHarmonicDriver()->GetActFreq();
        // as we control the fractional digits, we do not use lexical_cast<string>
        stringstream ss;
        ss << fixed << std::setprecision(1) << frequency;
        return ss.str();
      }
    }
    assert(false); // there shall be a harmonic driver!
    return "no freq found";
  }

  int ErsatzMaterial::GetSpecialResultIndex(App::Type app1, App::Type app2, CalcMode calcMode, Condition* constraint)
  {
    stringstream label;
    label << application.ToString(app1) << "_" << application.ToString(app2);
    if (calcMode == CONJ_QUAD)
      label << "_quad";

    DesignElement& de = design->data[0];
    int index = -1;
    DesignElement::Detail detail = de.NONE;
    if(de.detail.IsValid(label.str()))
    {
      detail = de.detail.Parse(label.str());
      DesignElement::ValueSpecifier vs = constraint != NULL ? de.CONSTRAINT_GRADIENT : de.COST_GRADIENT;
      index = design->GetSpecialResultIndex(de.DEFAULT, vs, detail, de.PLAIN);
    }
    LOG_DBG2(em) << "GetSpecialResultIndex: label='" << label.str() << "' detail=" << de.detail.ToString(detail) << " index=" << index;
    return index;
  }

  void ErsatzMaterial::CalcNewmarkDerivative(Excitation& excite, StateContainer& forward, StateContainer& adjoint, double factor, Objective* c, Condition* g)
  {
    assert(!context->DoMultiSequence());
    Assemble* assemble = context->pde->GetAssemble();

    // this calculates p^T (dF - dA) u
    // where p is solution of adjoint, dF is derivative of newmark update, dA is derivative of system matrix, u is solution of forward problem
    Function* f = Function::Cast(c, g);
    UInt timesteps = context->GetDriver()->GetNumSteps();
    MathParser* parser = domain->GetMathParser();
    unsigned int mathParserHandle = parser->GetNewHandle();
    // FIXME assert(domain->HasErsatzMaterialTensor());
    Matrix<double> dK(1, 1), dM(1, 1);
    Vector<double> dKp(0), dMp(0), dDp(0);
    // this is only caching, access using the double maps for every get slows down this procedure by about 50% for 200 timesteps
    StdVector<StdVector<SingleVector*> *> forwards;
    StdVector<StdVector<SingleVector*> *> forwarddt;
    StdVector<StdVector<SingleVector*> *> forwarddtt;
    StdVector<StdVector<SingleVector*> *> adjoints;
    forwards.Resize(timesteps);
    forwarddt.Resize(timesteps);
    forwarddtt.Resize(timesteps);
    adjoints.Resize(timesteps);
    for(unsigned int t = 0; t < timesteps; ++t)
    {
      forwards[t] = (&forward.Get(excite, NULL, t)->elem[App::MECH]);
      forwarddt[t] = (&forward.Get(excite, NULL, t, FIRST_DERIV)->elem[App::MECH]);
      forwarddtt[t] = (&forward.Get(excite, NULL, t, SECOND_DERIV)->elem[App::MECH]);
      adjoints[t] = (&adjoint.Get(excite, f, t)->elem[App::MECH]);
    }
    const TransferFunction* ktf = design->GetTransferFunction(DesignElement::DENSITY, App::MECH);
    const TransferFunction* mtf = design->GetTransferFunction(DesignElement::DENSITY, App::MASS);
    // the outer most loop is over all elements, so element matrices can be reused as much as possible
    int upper = design->data.GetSize();
    int elements = design->GetNumberOfElements();
    for (int base = 0;base < upper;base += elements)
    {
      // loop over all designs
      for (int e = 0;e < elements;++e)
      {
        // loop over all elements
        DesignElement* de = &design->data[base + e];
        bool notDampingElement = de->GetType() != DesignElement::DAMPINGALPHA && de->GetType() != DesignElement::DAMPINGBETA;
        SetElementK(f->ctxt, de, ktf, App::MECH, dynamic_cast<DenseMatrix*>(&dK), notDampingElement);
        SetElementK(f->ctxt, de, mtf, App::MASS, dynamic_cast<DenseMatrix*>(&dM), notDampingElement);
        // The damping matrix is alpha * Mass + beta * Stiffness, so it's derivative is also alpha * dMass + beta * dStiffness
        // We need to get alpha and beta, from the integrators
        // if we get Damping Information from the DesignSpace, we use that, else we use the "traditional" one
        double dampingAlpha = 0.0;
        double dampingBeta = 0.0;
        //if(!design->designMaterial->GetMaterialDamping(dampingAlpha, dampingBeta, de->elem, notDampingElement ? DesignElement::NO_DERIVATIVE : de->GetType()))
        assert(false); // FIXME
        if(false) // FIMXE
        {
          RegionIdType regionId = de->elem->regionId;
          SinglePDE* pde = context->pde;
          BiLinFormContext* linElastIntCtxt = assemble->GetBiLinForm("LinElastInt", regionId, pde, pde, false);
          BiLinFormContext* linMassIntCtxt = assemble->GetBiLinForm("MassInt", regionId, pde, pde, false);
          if (linElastIntCtxt->GetSecDestMat() != NOTYPE)
          {
            parser->SetExpr(mathParserHandle, linElastIntCtxt->GetSecMatFac());
            dampingBeta = parser->Eval(mathParserHandle);
          }
          if (linMassIntCtxt->GetSecDestMat() != NOTYPE)
          {
            parser->SetExpr(mathParserHandle, linMassIntCtxt->GetSecMatFac());
            dampingAlpha = parser->Eval(mathParserHandle);
          }
        }

        const bool damping = dampingAlpha > 0.0 || dampingBeta > 0.0;
        double v = 0.0;
        for (unsigned int t = 0;t < timesteps;++t)
        {
          // loop over all time steps
          Vector<double>& p_vec = dynamic_cast<Vector<double>&>(*(*adjoints[t])[e]);
          dKp = dK * p_vec;
          if (notDampingElement)
          {
            // K/dDamp = M/dDamp = 0
            Vector<double>& u_vec = dynamic_cast<Vector<double>&>(*(*forwards[t])[e]);
            v -= u_vec * dKp;
          }
          if (t > 0 || !IsFirstTransientStepStatic()) // all transiently calculated timesteps
          {
            dMp = dM * p_vec;
            if (notDampingElement)
            {
              Vector<double>& utt_vec = dynamic_cast<Vector<double>&>(*(*forwarddtt[t])[e]);
              v -= utt_vec * dMp;
            }
            if (damping)
            {
              // just performance; if alpha = beta = 0, we can omit one vector product
              Vector<double>& ut_vec = dynamic_cast<Vector<double>&>(*(*forwarddt[t])[e]);
              assert(false);
              // FIXME dDp = dampingAlpha * dMp + dampingBeta * dKp;
              v -= ut_vec * dDp;
            }
          }

        } // loop over timesteps

        de->AddGradient(c, g, factor * v);
      }

    }

    parser->ReleaseHandle(mathParserHandle);
  }

  double ErsatzMaterial::CalcU1KU2(TransferFunction* tf, StdVector<SingleVector*>& u1, App::Type k, StdVector<SingleVector*>& u2, DesignDependentRHS* rhs, double factor, CalcMode calcMode, Function* f, int res_idx, double ev)
  {
    //Special case when doing mapping optimization
    assert(f != NULL);
    if(f->ctxt->IsComplex())
      return CalcU1KU2<std::complex<double> >(tf, u1, k, u2, rhs, factor, calcMode, f, res_idx, ev);
    else
      return CalcU1KU2<double>(tf, u1, k, u2, rhs, factor, calcMode, f, res_idx, ev);
  }

  template<class T>
  double ErsatzMaterial::CalcU1KU2(TransferFunction* tf, StdVector<SingleVector*>& u1, App::Type app, StdVector<SingleVector*>& u2, DesignDependentRHS* rhs, double factor, CalcMode calcMode, Function* f, int res_idx, double ev)
  {
    if(calc_u1ku2_timer_)
      calc_u1ku2_timer_->Start();
    // LOG_DBG2(em) << "CalcU1KU2: tf=" << (tf ? tf->ToString() : "NULL") << " app=" << application.ToString(app) << "(" << app << ")"
    //              << " #u1=" << u1.GetSize() << " #u2=" << u2.GetSize() << " calcMode=" << calcMode << " factor=" << factor << " rhs=" << (rhs == NULL && rhs->vec == NULL ? "NULL" : rhs->ToString(1)) << " ev=" << ev;
    // This solves <l,K'*u-f'> or <u1, K' * u2 - f'> for all elements and adds it up to the element gradients
    // Note to perform "<f',u>" from <f',u> + <l,K'*u-f'> manually
    assert(u1.GetSize() != 0);
    assert(u1.GetSize() == u2.GetSize());
    assert(f != NULL); // for context or relax

    double sum = 0.0;
    // mat will be filled by SetElementK where also the derivative form most cases is built in
    // the dimensions of our matrix is determined by u1_vec and u2_vec.
    Matrix<T> mat(u1[0]->GetSize(), u2[0]->GetSize());//NOTE: SetElementK (In PiezoSimp) relies on the matrix already having the right size!!!
    Vector<T> mat_vec(u1[0]->GetSize());
    TransferFunction* rtf = rhs != NULL && rhs->valid ? design->GetTransferFunction(tf->GetDesign(), rhs->app) : NULL;

    // the context->GetExcitation() is now the last one as we solve and store all excitations first before calculating the gradients
    Transform* trans = f != NULL && f->GetExcitation() != NULL ? f->GetExcitation()->transform : NULL; // even ->transform might be NULL

    // traverse over our elements
    // in ErsatzMaterialTensor case we loop over all elements, else only over the elements belonging to this design
    int elements = design->GetNumberOfElements();
    int base_lower = 0;
    int base_upper = design->data.GetSize(); // ErsatzMatzerialTensor and MultiMaterial
    if(design->designMaterial == NULL && !design->HasMultiMaterial())
    {
      base_lower = design->FindDesign(tf->GetDesign()) * elements;
      base_upper = base_lower + elements;
    }
    LOG_DBG2(em) << "elements=" << elements << " base=" << base_lower << " base_upper=" << base_upper;
    // create an element list to gain the iterator in the loop
    ElemList elemList(grid);

    // for ParamMat we need the derivative w.r.t. every designvariable, else the base loop is only run once
    for(int base = base_lower; base < base_upper; base += elements)
    {
      for(int e = 0; e < elements; e++)
      {
        // if we do transformation, the physical u is calculated based on transformed elements
        Vector<T>& u1_vec = dynamic_cast<Vector<T>& >(*u1[e]);
        Vector<T>& u2_vec = dynamic_cast<Vector<T>& >(*u2[e]);

        DesignElement* org = &design->data[e + base];

        // de is the potentially transformed stuff. Note, that we also store the stuff for the transformed element!
        // the general idea about gradients of transformation is the following
        // - in the forward problem the state is calculated for the transformed (rotated) element
        // - for the gradient the state is already transformed, we no need ONLY to transform the element index for
        //   - dK/drho
        //   - storing the gradient
        //   - do NOT use the element state for the transformed element -> this has already been done for the forward problem!
        DesignElement* de = design->ApplyTransformations(org, org, trans); // fallback to design if there is no transformation

        LOG_DBG3(em) << "nodes:" << e << ": " << de->elem->connect.ToString() << " dt=" << de->type.ToString(de->GetType()) << " e=" << de->elem->elemNum;
        LOG_DBG3(em) << "u1:" << e << ": " << u1_vec.ToString();
        LOG_DBG3(em) << "u2:" << e << ": " << u2_vec.ToString();

        // u1^T (K' u2 - f') -> find "K'"
        SetElementK(f->ctxt, de, tf, app, dynamic_cast<DenseMatrix*>(&mat), true, calcMode, ev); // derivative = true
        LOG_DBG3(em) << "mat: " << mat.ToString();

        // We generally solve u1^T (K' u2 - f')
        // u1^T (K' u2 - f') -> calc "K' u2"
        mat_vec = mat * u2_vec;
        LOG_DBG3(em) << "mat * u2: " << mat_vec.ToString();

        // u1^T (K' u2 - f') -> calc "- f'"
        assert(!(calcMode == CONJ_QUAD && rtf != NULL));// no sensitive rhs here!
        assert(!(rtf != NULL && f->ctxt->IsStrainExcitedSystem()));

        if(rtf != NULL) {
          if (rhs->isInterfaceDriven_)
            SubstractInterfaceDrivenGradRHS(f, tf, de, mat_vec);
          else
            SubtractGradSurfaceRHS(de, rtf, rhs, mat_vec);
        }

        if(f->ctxt->IsStrainExcitedSystem())
          SubtractGradStrainRHS(de, tf, rhs, mat_vec);

        LOG_DBG3(em) << "-f': " << mat_vec.ToString();

        // u1^T(K' u2 - f') -> calc "u1^T *" or <u1, *>
        // the difference is the conjugate complex in the harmonic inner product case!
        T sp;
        if(calcMode == CONJ_QUAD || calcMode == EIGENFREQ)
          mat_vec.Inner(u1_vec, sp);// u1 = u2 = u!
        else
          sp = mat_vec * u1_vec;

        // when doing complex Jensen 22.07.07 shows that we always have 2 * Re(lamda * grad S * u)
        // the factor gives the negative sign
        // in real case it is simple value = factor * sp.
        // factor shall be +/- 1!
        double this_value = factor;
        if(f->ctxt->IsHarmonic() && calcMode == STANDARD)
          this_value *= 2 * ((complex<double>) sp).real();// 2 * Re{...}
        else
          this_value *= ((complex<double>) sp).real();// CONJ_QUAD, EIGENFREQ or real STANDARD

        de->AddGradient(f, this_value);

        LOG_DBG3(em) << "CalcU1KU2: e=" << e << "->" << de->GetIndex() << " de=" << de->ToString() << " <l,K'*u-f'>  = " << sp << " -> " << this_value << " sum = " << de->GetPlainGradient(f);

        sum += this_value;

        if(res_idx != -1) de->specialResult[res_idx] = this_value;
      }
    }
    if(calc_u1ku2_timer_)
      calc_u1ku2_timer_->Stop();
    return sum;
  }


  void ErsatzMaterial::AddMassToStiffness(Context* ctxt, const TransferFunction* mtf, DesignElement* de, Matrix<complex<double> >& K_in_S_out, bool derivative, bool bimaterial, CalcMode mode, double ev)
  {
    // The result matrix is
    // S = K + i*omega*C - omega^2*M
    // with purely imaginary C = alpha_k*K+alpha_m*M
    // S = K + i*omega*alpha_k*K + i*omega*alpha_m*M - omega^2*M
    // with m_factor which is the transfer function. The K transfer function is already in S
    // S = K + i*omega*alpha_k*K + i*omega*alpha_m*m_factor*M - omega^2*m_factor*M
    // With S = k_factor K in the beginning this is: (k_factor is transfer function or its derivative)
    // S += i*alpha_k*S + (i*alpha_m*m_factor-omega^2*m_factor)*M
    //
    // in case we have pamping (e.g. Sigmund; Morhology; 2007) there is to add
    // j*omega*pamping*rho*(1-rho)*M and for the derivative case
    // j*omega*pamping*rho'*M - j*2*omega*pamping*rho*rho'*M = j*omega*pamping*rho'(1-2*rho)
    //
    // the eigenvalue derivative is u^T (K' - ev M') u
    Assemble* assemble = ctxt->pde->GetAssemble(); // shall work even if current context != ctxt

    double mtv(0.0), mdv(0.0), m_factor(1.0);

    if(this->method_ != ErsatzMaterial::PARAM_MAT) // density is treated in Mass(...) function in case of ParamMat
    {
      mtv =  mtf->Transform(de, DesignElement::SMART, bimaterial);
      mdv =  mtf->Derivative(de, DesignElement::SMART, bimaterial);
      m_factor = derivative ? mdv : mtv;
    }
    assert(mode != EIGENFREQ || (derivative == true && ev > 0)); // EIGENVALUE only for derivative
    if(mode == EIGENFREQ)
      m_factor *= ev;
    // change name only
    Matrix<complex<double> >& S = K_in_S_out;
    LOG_DBG3(em) << "AMTS: 1. e=" << de->elem->elemNum << " ev=" << ev << " m_factor=" << m_factor << " K_in_S_out=" << S.ToString();
    // find alpha, beta and omega. We have no omega for the eigenvalue case and 1.0 eliminates it
    double omega = mode != EIGENFREQ ? 2.0 * M_PI * ctxt->pde->GetSolveStep()->GetActFreq() : 1.0 ;  // todo: check with multiple excitation frequencies!
    double alpha_k = 0.0;
    double alpha_m  = 0.0;
    double pamping_m = 0.0; // add on without omega
    // do we have damping (C = alpha*M+beta*K) -> this is pure imaginary!
    RegionIdType regionId = de->elem->regionId;
    if(ctxt->pde->GetDamping(regionId) == RAYLEIGH)
    {
      assert(mode != EIGENFREQ);
      SinglePDE* pde = ctxt->pde;
      assert(pde != NULL);
      // the alpha and beta might be calculated and adjusted, get them
      // from the integrators in the form as they are used for the state problem!
      alpha_k = assemble->GetBiLinForm("LinElastInt", regionId, pde, pde)->EvalSecMatFac();
      // now alpha_m
      alpha_m = assemble->GetBiLinForm("MassInt", regionId, pde, pde)->EvalSecMatFac();
      assert(alpha_k > 0 && alpha_m > 0&& omega > 0);
      // pamping stuff without omega
      double pamping = design->GetPampingValue(); // 0 if not applicable
      if(!derivative)
        pamping_m = pamping * mtv * (1.0 - mtv);
      else // pamping*rho'(1-2*rho)
        pamping_m = pamping * mdv * (1.0 - 2.0 * mtv);
      assert(this->method_ != ErsatzMaterial::PARAM_MAT || pamping == 0.0);
    }
    assert(mode != EIGENFREQ || (omega == 1.0 && m_factor != 0 && alpha_m == 0.0 && pamping_m == 0.0)); // note that we might have very_small negative eigenvalues!
    const unsigned int srows = S.GetNumRows();
    const unsigned int scols = S.GetNumCols();
    // we first add the K part of C (= pure imaginary). E.G. in the bloch case S=K might already have an imaginary part
    for(unsigned int r = 0; r < srows; r++)
      for(unsigned int c = 0; c < scols; c++)
        S[r][c] = complex<double>(S[r][c].real(), S[r][c].imag() + omega * alpha_k * S[r][c].real());
    LOG_DBG3(em) << "AMTS: 2. e=" << de->elem->elemNum << " add K o=" << omega << " a_K=" << alpha_k << " S=" << S.ToString();
    // we the add the M part of C and the real mass part
    complex<double> damp_mass = complex<double>(-1.0 *omega*omega*m_factor, omega*(alpha_m*m_factor  + pamping_m));
    // multimaterial stuff
    int index = de->multimaterial != NULL ? de->multimaterial->index : -1;
    LOG_DBG3(em) << "AMTS: e=" << de->elem->elemNum << " S=" << S.ToString();
    if(ctxt->mat->ComplexElementMatrix(de->elem->regionId))
    {
      // only accessed as derivative in ParamMat case
      assert(this->method_ != ErsatzMaterial::PARAM_MAT);
      const Matrix<Complex>& M = dynamic_cast<const Matrix<Complex>&>(ctxt->mat->Mass(de->elem, bimaterial, index, (this->method_ == ErsatzMaterial::PARAM_MAT) ? de->GetType() : DesignElement::NO_DERIVATIVE));
      assert(S.GetNumRows() == M.GetNumRows() && S.GetNumCols() == M.GetNumCols());
      Add<Complex, Complex>(S, damp_mass, M);
      LOG_DBG3(em) << "AMTS: 3. complex e=" << de->elem->elemNum << " damp_mass=" << damp_mass << " S=" << S.ToString();
    }
    else
    {
      // only accessed as derivative in ParamMat case
      const Matrix<double>& M = dynamic_cast<const Matrix<double>&>(ctxt->mat->Mass(de->elem, bimaterial, index, (this->method_ == ErsatzMaterial::PARAM_MAT) ? de->GetType() : DesignElement::NO_DERIVATIVE));
      assert(S.GetNumRows() == M.GetNumRows() && S.GetNumCols() == M.GetNumCols());
      Add<Complex, double>(S, damp_mass, M);
      LOG_DBG3(em) << "AMTS: 3. real e=" << de->elem->elemNum << " damp_mass=" << damp_mass << " S=" << S.ToString();
    }
    LOG_DBG2(em) << "AddMassToStiffness: d=" << de->elem->elemNum << " der=" << derivative << " bm=" << bimaterial << " mode="  << mode << " ev=" << ev
                   << " m_factor:" << m_factor << " alpha_k: " << alpha_k << " alpha_m: " << alpha_m << " pamping_m:" << pamping_m
                   << " omega: " << omega << " K_img: " << (omega * alpha_k) << " damp_mass: " << damp_mass;
    LOG_DBG3(em) << "AMTS: 4. e=" << de->elem->elemNum << " S -> " << S.ToString();
  }


  template<class T>
  void ErsatzMaterial::SubtractGradStrainRHS(DesignElement* de, TransferFunction* tf, DesignDependentRHS* rhs, Vector<T>& in_out)
  {
    assert(!context->DoMultiSequence());
    assert(rhs == NULL || rhs->app == App::STRESS);
    MechPDE::TestStrain ts = rhs != NULL ? rhs->test_strain : MechPDE::NOT_SET;
    // OptMechMat is base for any further child!
    const Vector<double>& vec = dynamic_cast<MechMat*>(context->mat)->MechStrainRHS(de->elem, ts);
    double factor = tf->Derivative(de, DesignElement::SMART);
    // LOG_DBG3(em) << "SGSR: de=" << de->elem->elemNum << " in_out=" << in_out.ToString();
    // LOG_DBG3(em) << "SGSR: de=" << de->elem->elemNum << "    vec=" << vec.ToString();
    // LOG_DBG3(em) << "SGSR: de=" << de->elem->elemNum << "    val=" << de->GetDesign(DesignElement::PLAIN) << " drho=" << factor;
    in_out.Add(-1.0 * factor, vec);// -1.0 as we want to subtract!
    // LOG_DBG3(em) << "SGSR: de=" << de->elem->elemNum << "     ->=" << in_out.ToString();
  }

  template<class T>
  void ErsatzMaterial::CalcAndStoreInterfaceDrivenGrad(Function* f, TransferFunction* tf)
  {
    if (interfaceDrivenGradCalc_)
      return;

    // get nodes where homogeneous Dirichlet BC is enforced
    shared_ptr<BaseFeFunction> fe = f->ctxt->pde->GetFeFunction(f->ctxt->pde->GetNativeSolutionType());
    StdVector<unsigned int> idBcNodes;
    IdBcList& idBcs = fe->GetInHomDirichletBCs();

    // find indices of nodes with a hom Dirichlet bc
    for (unsigned int i = 0; i < idBcs.GetSize(); i++) {
      EntityIterator entIt = idBcs[i]->entities->GetIterator();
      for ( ; !entIt.IsEnd(); entIt++)
        idBcNodes.Push_back(entIt.GetNode());
    }
    assert(f != NULL);
    //FIXME Assume design elements are all of the same type and application is HEAT

    for (unsigned int id = 0; id < design->data.GetSize(); id++) {
      DesignElement* de = &design->data[id];
      // differentiation factor when using filter
      // for each node we have f' =  4* ds_i /drho_i * (1-2*s_i)
      // except for bc nodes, there f' is 0
      StdVector<unsigned int>& nodes = de->elem->connect;
      for(unsigned int n = 0; n < nodes.GetSize(); n++) {
        unsigned int node = nodes[n];
        if(!idBcNodes.Contains(node)) // gradient is 0 at bc nodes
        {
          double factor = 0.0;
          if (f->ctxt->pde->GetParamNode()->Has("bcsAndLoads/designDependentHeatSource/value"))
            f->ctxt->pde->GetParamNode()->GetValue("bcsAndLoads/designDependentHeatSource/value",factor);
          else
            f->ctxt->pde->GetParamNode()->GetValue("bcsAndLoads/heatSource",factor);

          de->interfaceDrivenLoadGrad_[n] = design->EvalInterfaceFunction(node, true) / design->data.GetSize() * factor * tf->Derivative(de, DesignElement::SMART,false);

        } //if
      } // node
    } // elem

  } // function

  template<class T>
  void ErsatzMaterial::SubstractInterfaceDrivenGradRHS(Function* f, TransferFunction* tf, const DesignElement* de, Vector<T>& in_out)
  {
      if (!interfaceDrivenGradCalc_) {
        CalcAndStoreInterfaceDrivenGrad<double>(f,tf);
        interfaceDrivenGradCalc_ = true;
      }

      assert(in_out.GetSize() > 0);
      assert(!f->ctxt->IsComplex());
      assert(in_out.GetEntryType() == BaseMatrix::DOUBLE);

      Vector<double>& vec = dynamic_cast<Vector<double>& >(in_out);
      assert(vec.GetSize() == de->interfaceDrivenLoadGrad_.GetSize());
      vec -= de->interfaceDrivenLoadGrad_;
  }

  template<class T>
  void ErsatzMaterial::SubtractGradSurfaceRHS(DesignElement* de, TransferFunction* tf, DesignDependentRHS* ref, Vector<T>& in_out)
  {
    // we have to find the nodes which are common between de->elem
    // and the surface element which is one dimension smaller
    // not all elements do necessary lay on a surface and then not all nodes
    assert(ref != NULL && ref->valid);
    // nodes (numbers) of our design element
    StdVector<unsigned int>& de_nodes = de->elem->connect;
    Vector<T>& rhs = dynamic_cast<Vector<T>&>(*(ref->vec));
    // in_out is scalar (potential) or vectorial (x,y,(z))
    assert(in_out.GetSize() >= de_nodes.GetSize());
    int dof = in_out.GetSize() / de_nodes.GetSize();
    assert(dof == 1 || dof == 2 || dof == 3);
    assert(dof == (int) rhs.GetSize());
    // all node numbers of the surface are in a set
    std::set<unsigned int>::iterator it;
    // compare with the node numbers of our design element
    for(unsigned int n = 0; n < de_nodes.GetSize(); n++)
    {
      it = ref->nodes.find(de_nodes[n]);
      if(it != ref->nodes.end())
      {
        LOG_DBG3(em) << "SubtractGradSurfaceRHS : node " << n << " is common with elem "
        << *it << " in surface: K'u = " << in_out.ToString();

        // find the the sensitivity of the rhs w.r.t the design volume element!
        double factor = tf->Derivative(de, DesignElement::SMART);

        // we do not really construct a rhs vector (with some/many zeros) but substract
        // for all design nodes common with the surface directly the entries
        for(int d = 0; d < dof; d++)
        in_out[n*dof + d] -= factor * rhs[d];
        LOG_DBG3(em) << "... -" << factor << "*" << rhs.ToString() << " -> " << in_out.ToString();
      }
    }
  }

  double ErsatzMaterial::CalcHomTensor(Objective* c, Condition* g, bool derivative)
  {
    Function* f = Function::GetFunction(c, g);
    if(c != NULL && derivative && c->HasHomogenizationEntry())
    {
      // if there s no "coord" set it is only meant for evaluate for forward homogenization
      StdVector<double> tmp;
      CalcHomogenizedTensorEntry(f->ctxt, c->coord, true, tmp, f->GetExcitation()->meta_index);
      for(unsigned int e = 0, ne = design->GetNumberOfElements(); e < ne; e++)
      design->data[e].AddGradient(c, NULL, tmp[e]);
      return 0.0;
    }
    if(c != NULL && !derivative)
    {
      Matrix<double> hom_tensor = CalcHomogenizedTensor(f);
      if(c->HasHomogenizationEntry())
      {
        return hom_tensor[boost::get<0>(c->coord)-1][boost::get<1>(c->coord)-1];
      }
      else
      {
        std::cout << "Homogenized Tensor: " << std::endl << hom_tensor.ToString(0, true);

        std::cout << "Isotrope properties: ";
        SubTensorType stt = f->ctxt->stt;
        std::cout << " E=" << MechanicMaterial::CalcIsotropicYoungsModulus(hom_tensor, stt);
        std::cout << " v=" << MechanicMaterial::CalcIsotropicPoissonsRatio(hom_tensor, stt);
        std::cout << " err=" << MechanicMaterial::CalcIsotropyError(hom_tensor, stt) << "\n";

        StdVector<std::pair<string, double> > ortho = GetOrthotropeProperties(hom_tensor, f->GetExcitation());
        std::cout << "Orthotrope properties: ";
        if(ortho.GetSize() == 0)
        std::cout << " in 2D only for single region\n";

        for(unsigned int i = 0; i < ortho.GetSize(); i++)
        std::cout << " " << ortho[i].first << "=" << ortho[i].second;
        std::cout << "\n";

        return hom_tensor.NormL2();
      }
    }
    if(g != NULL)
    {
      return CalcHomogenizedTensorConstraint(g, derivative);
    }
    return 0.0;
  }

  double ErsatzMaterial::CalcFunction(Excitation& excite, Function* f, bool derivative)
  {
    interfaceDrivenGradCalc_ = false;

    assert(f != NULL);
    // for legacy reasions there is also the difference between Objective and Condition, to be replaced once
    Objective* c = f->IsObjective() ? dynamic_cast<Objective*>(f) : NULL;
    Condition* g = f->IsObjective() ? NULL : dynamic_cast<Condition*>(f);
    double result = 0.0;// stays  for the derivative
    switch(f->GetType())
    {
      case Function::VOLUME:
      case Function::PENALIZED_VOLUME:
      case Function::GAP:
      result = CalcVolume(c, g, derivative, true);
      break;

      case Function::REALVOLUME:
      result = CalcVolume(c, g, derivative, false);
      break;

      case Objective::TYCHONOFF:
      result = IntegrateDesignVariable(c, g, derivative, DesignElement::NO_TYPE, true, true, 2.0);
      break;

      case Function::COMPLIANCE:
      result = CalcCompliance(excite, c, g, derivative);
      break;

      case Objective::TRACKING:
      result = CalcTracking(excite, c, g, derivative);
      break;

      case Objective::TEMP_TRACKING_AT_INTERFACE:
      {
        Vector<double> res;
        result = CalcStateTrackingAtInterface(excite, f, derivative, f->GetParameter());
        break;
      }

      case Function::GREYNESS:
      assert(c == NULL);
      result = CalcGreyness(g, derivative);
      break;

      case Function::FILTERING_GAP:
      result = CalcFilteringGap(g,derivative);
      break;

      case Objective::STRESS:
      case Objective::STRESS_DENSITY:
      {
        // copy data for element von Mises stress
        Vector<double> data;
        if(f->ctxt->IsComplex())
        {
          StressConstraint<complex<double> > sc(&excite, f, this, &forward);
          sc.CalcStresses(data);
        }
        else
        {
          StressConstraint<double> sc(&excite, f, this, &forward);
          sc.CalcStresses(data);
        }
        result = CalcGlobalFunction(f, false, &data);
        break;
      }

      case Function::GLOBAL_SLOPE:
      case Function::GLOBAL_MOLE:
      case Function::GLOBAL_OSCILLATION:
      case Function::GLOBAL_JUMP:
      case Function::GLOBAL_SUM_MODULI:
      case Function::GLOBAL_TWO_SCALE_VOL:
      case Function::GLOBAL_TENSOR_TRACE:
      case Function::GLOBAL_ORTHOTROPIC_TENSOR_TRACE:
      case Function::GLOBAL_CURVATURE:
      case Function::GLOBAL_DESIGN:
      case Function::PERIMETER:
        result = CalcGlobalFunction(f, derivative);
      break;

      case Function::SLOPE:
      case Function::MOLE:
      case Function::OSCILLATION:
      case Function::JUMP:
      case Function::BUMP:
      case Function::CURVATURE:
      case Function::OVERHANG_VERT:
      case Function::OVERHANG_HOR:
      case Function::PERIODIC:
      case Function::SUM_MODULI:
      case Function::TWO_SCALE_VOL:
      case Function::ORTHOTROPIC_TENSOR_TRACE:
      case Function::TENSOR_TRACE:
      case Function::TENSOR_NORM:
      case Function::PARAM_PS_POS_DEF:
      case Function::POS_DEF_DET_MINOR_1:
      case Function::POS_DEF_DET_MINOR_2:
      case Function::POS_DEF_DET_MINOR_3:
      case Function::BENSON_VANDERBEI_1:
      case Function::BENSON_VANDERBEI_2:
      case Function::BENSON_VANDERBEI_3:
      case Function::DESIGN:
      case Function::MULTIMATERIAL_SUM:
      case Function::SHAPE_INF:
      assert(c == NULL);
      result = CalcLocalConstraint(g, derivative);
      break;

      case Function::DESIGN_TRACKING:
      assert(c == NULL);
      result = CalcDesignTracking(g, derivative);
      break;

      case Function::HOM_TENSOR:
      result = CalcHomTensor(c, g, derivative);
      break;

      case Function::HOM_TRACKING:
      if(derivative)
      {
        CalcHomogenizedTrackingGradient(f->GetTensor(), CalcHomogenizedTensor(f), f);
      }
      else
      {
        double diff = f->GetTensor().DiffNormL2(CalcHomogenizedTensor(f));
        result = 0.5 * diff * diff;
      }
      break;

      case Function::HOM_FROBENIUS_PRODUCT:
        if(derivative)
          CalcHomFrobeniusProductGradient(f->GetTensor(), CalcHomogenizedTensor(f), f);
        else
          return f->GetTensor().FrobeniusProduct(CalcHomogenizedTensor(f));
        break;

      case Function::POISSONS_RATIO:
      case Function::YOUNGS_MODULUS:
      case Function::YOUNGS_MODULUS_E1:
      case Function::YOUNGS_MODULUS_E2:
      result = CalcPoissonsRatioAndYoungsModulus(f, derivative);
      break;

      case Function::GLOBAL_DYNAMIC_COMPLIANCE:
      assert(!derivative); // SIMP!
      result = CalcGlobalDynamicCompliance(excite, f);
      break;

      case Function::OUTPUT:
      case Function::SQUARED_OUTPUT:
      case Function::DYNAMIC_OUTPUT:
      case Function::CONJUGATE_COMPLIANCE:
      case Function::ABS_OUTPUT:
      if(f->ctxt->IsComplex())
      {
        if (derivative)
        {
          App::Type app = Context::ToApp(f->ctxt->pde);
          // synthesis of compliant mechanism: As our adjoint PDE
          // c' = l K' u
          TransferFunction* tf = design->GetTransferFunction(DesignElement::Default(f->ctxt), TransferFunction::Default(f->ctxt), true, true); // excpetion and use_single
          double weight = excite.GetWeightedFactor(f);
//          LOG_DBG(simp) << "CalcFunction(idx=" << excite.index << ") norm_weight= " <<  excite.normalized_weight  << " factor=" << excite.GetFactor(f) << " weight=" << weight;
          CalcU1KU2(tf, adjoint.Get(excite, f)->elem[app], app, forward.Get(excite)->elem[app], NULL, weight, STANDARD, f);
          return 0.0;
        }
        else
          result = CalcOutput<complex<double> >(excite, f);
      }
      else
        result = CalcOutput<double>(excite, f);
      break;

      case Function::ENERGY_FLUX:
      assert(!derivative);
      result = CalcEnergyFlux(excite, c);
      break;

      case Function::TEMPERATURE:
      break;// FIXMEHEAT

      case Function::HEAT_ENEGRY:
      result = CalcHeatEnergy(excite, c, g, derivative);
      break;

      case Function::ELEC_ENERGY:
        assert(false);// shall be handled before
        break;

      case Function::PRESSURE_DROP:
        if (!derivative)
          result = f->ctxt->GetLatticeBoltzmannPDE()->CalcPressureDrop();
        else
          f->ctxt->GetLatticeBoltzmannPDE()->SensitivityAnalysis(design->GetTransferFunction(f->elements[0]), f, design);
        break;

      case Function::SLACK:
        if(!derivative)
          result = design->GetSlackVariable();
        else
          dynamic_cast<AuxDesign*>(design)->GetSlackDesign()->AddGradient(f, 1.0);
        break;
      case Function::EIGENFREQUENCY:
        result = CalcEigenfrequency(excite, f, derivative);
        break;

      case Function::BANDGAP:
        result = CalcBandGap(excite, f, derivative);
        break;

      case Function::SLACK_FNCT:
        result = CalcSlackFunction(f, derivative);
        break;

      case Function::EXPRESSION:
        result = CalcExpression(g, derivative);
        break;

      case Function::ISOTROPY:
      case Function::ISO_ORTHOTROPY:
      case Function::ORTHOTROPY:
      case Function::MULTI_OBJECTIVE:
        assert(false);// no valid function
        break;
      // no default, gcc warns
    }
    LOG_DBG2(em) << "CalcFunction " << f->ToString() << " cost=" << f->IsObjective() << " -> " << (derivative ? "derivative" : lexical_cast<std::string>(result));
    return result;
  }

  double ErsatzMaterial::IntegrateDesignVariable(Objective* c, Condition* g, bool derivative, DesignElement::Type dtype, bool normalized, bool scale, double exponent)
  {
    // this replaces and enhances calculation of volume, it is used by regularization
    // when not assuming a regular grid, computation of Volume is not as simple
    // we also consider working only on a given region, when used as constraint
    // use dtype == NO_TYPE to iterate over all designs, but do not calculate tensor trace even if available
    // do we want the physical value? Don't make GTF() fault tolerant as we assume the physical value!
    Grid* grid = domain->GetGrid();
    Function* f = Function::GetFunction(c, g);
    SubTensorType stt = f->ctxt->stt;
    TransferFunction* tf = Function::GetFunction(c, g)->IsPhysical() ? design->GetTransferFunction(dtype, App::MECH) : NULL;

    double fraction = c != NULL ? volume_fraction_ : g->volume_fraction;
    bool allDesignsRelevant = dtype == DesignElement::MECH_TRACE  || dtype == DesignElement::DIELEC_TRACE || dtype == DesignElement::DEFAULT || dtype == DesignElement::NO_TYPE;
    // tensor trace is calculated if dtype == DEFAULT or TENSOR_TRACE and a tensor available
    bool calculateTensorTrace = design->designMaterial != NULL && (dtype == DesignElement::MECH_TRACE || dtype == DesignElement::DIELEC_TRACE || dtype == DesignElement::DEFAULT);
    if (calculateTensorTrace && scale)
    {
      throw Exception("Cannot calculate Tensor Trace Volume on scaled design variables!");
    }
    // check whether we have to calculate the full volume
    if(normalized)
    {
      if(g != NULL && g->GetDesignType() == DesignElement::UNITY)
      { // this will always return 1, does only make sense for unnormalized (real) volume
        // the gradient is not set, it is really 0 on every element but should not be used
        //TODO: Do We need a warning here?
        return(derivative ? 0.0 : 1.0);
      }
    }
    else // no normalization needed, we set factor
    {
      if(design->IsRegular())
      { // we use 1/(the volume of the first element) as fraction
        fraction = 1.0 / grid->GetElemShapeMap(design->data[0].elem, false)->CalcVolume();
      }
      else
      {
        fraction = 1.0; // if no normalization needed, fraction is simply 1.0
      }
    }
    const unsigned int nd = design->regions.GetSize();
    if( fraction == 0.0 )
    {
      for(unsigned int d = 0; d < nd; d++)
      { // tensortrace breaks out after the first design
        const unsigned int nr = design->regions[d].GetSize();
        if(allDesignsRelevant || design->design[d].design == dtype)
        {
          for(unsigned int r = 0; r < nr; r++)
          {
            DesignSpace::DesignRegion& cur_reg = design->regions[d][r];
            if(c != NULL || g->IsForRegion(cur_reg.regionId))
            {
              if(design->IsRegular())
              {
                fraction += cur_reg.elements;
              }
              else
              {
                const unsigned int u = cur_reg.base + cur_reg.elements;
                for(unsigned int i = cur_reg.base; i < u; i++)
                {
                  DesignElement* de = &design->data[i];
                  fraction += grid->GetElemShapeMap(de->elem, false)->CalcVolume();
                }
              }
            }
          }
        }
        if(calculateTensorTrace)
        {
          break;
        }
      }
      fraction = 1.0 / fraction;
      if(g == NULL)
      {
        volume_fraction_ = fraction;
      }
      else
      {
        g->volume_fraction = fraction;
      }
    }
    double sum = 0.0;
    for(unsigned int d = 0; d < nd; d++)
    {
      if(allDesignsRelevant || design->design[d].design == dtype)
      { // tensortrace breaks out after the first design
        const unsigned int nr = design->regions[d].GetSize();
        for(unsigned int r = 0; r < nr; r++)
        {
          DesignSpace::DesignRegion& cur_reg = design->regions[d][r];
          if(c != NULL || g->IsForRegion(cur_reg.regionId))
          {
            const double scaling = scale ? cur_reg.scale_design : 1.0;
            const double rscaling = 1.0 / scaling;
            const double translation = scale ? cur_reg.translate_design : 0.0;
            const unsigned int u = cur_reg.base + cur_reg.elements;
            for(unsigned int i = cur_reg.base; i < u; i++)
            {
              DesignElement* de = &design->data[i];
              // standard or derivative case?
              if(derivative)
              {
                double val = 0.0;
                if(calculateTensorTrace)
                {
                  Matrix<double> material;
                  design->designMaterial->GetTensor(material, dtype, stt, de->elem, de->GetType(), f->GetNotation());
                  val = material.Trace();
                  if(exponent != 1.0)
                  {
                    // chain rule, original, non derived tensor
                    design->designMaterial->GetTensor(material, dtype, stt, de->elem, DesignElement::NO_DERIVATIVE, f->GetNotation());
                    double des = material.Trace();
                    val *= exponent * std::pow(des, exponent - 1.0);
                    LOG_DBG(em) << " material = "<< material.ToString();
                    LOG_DBG(em) << " trace = "<< des;
                  }
                }
                else
                {
                  if(exponent != 1.0)
                  {
                    // by tf we mark if we want the physical stuff.
                    double des_val = tf != NULL ? tf->Derivative(de, DesignElement::SMART) : de->GetDesign(DesignElement::SMART);
                    val = exponent * std::pow(((des_val - translation) * rscaling), exponent - 1.0);
                  }
                  else
                  val = 1.0;
                }
                val *= rscaling; // the gradient will be multiplied by scaling later again, but if it is calculated on the scaled designs here (as in regularization) this should not happen
                val *= fraction;
                if(!design->IsRegular())
                {
                  const double vol = grid->GetElemShapeMap(de->elem, false)->CalcVolume();
                  val *= vol;
                }
                de->AddGradient(c, g, val);
              }
              else
              {
                // no derivative
                double des;
                if(calculateTensorTrace)
                { // use the trace of the stiffness Tensor as "volume"
                  Matrix<double> material;
                  design->designMaterial->GetTensor(material, dtype, stt, de->elem, DesignElement::NO_DERIVATIVE, f->GetNotation());
                  des = material.Trace();
                }
                else
                {
                  // tf marks the physical function
                  double des_val = tf != NULL ? tf->Transform(de, DesignElement::SMART) : de->GetDesign(DesignElement::SMART);
                  des = (des_val - translation) * rscaling;
                }
                des = std::pow(des, exponent);
                if(design->IsRegular())
                {
                  sum += des;
                }
                else
                {
                  const double vol = grid->GetElemShapeMap(de->elem, false)->CalcVolume();
                  sum += des * vol;
                }
              } // if derivative
            } // for element (i)
          } // right region
        } // for region
        if(calculateTensorTrace && !derivative)
        { // the derivative has to be calculated for all designs but the tensor itself only once
          break;
        }
      } // if relevant
    }
    return derivative ? -1.0 : sum * fraction;
  }

  double ErsatzMaterial::CalcVolume(Objective* f, Condition* g, bool derivative, bool normalized)
  {
    Function* func = Function::GetFunction(f, g);
    DesignElement::Type des = func->GetDesignType();
    switch(func->GetType())
    {
      case Function::VOLUME:
      // Bastian's stuff needs IntegrateDesignVariable(). The SIMP stuff is better with CalcTrivialVolume() as
      // we cannot assume SIMP_TYPE transfer functions here!
      if(method_ == SIMP_METHOD || GetDesign()->FindDesign(DesignElement::DENSITY, false) > -1)
        return CalcTrivialVolume(func, derivative, normalized);
      else// FIXME check if it is ok not to give an exponent in the physical case!
        return IntegrateDesignVariable(f, g, derivative, des, normalized, false, 1.0);// no scaling, exponent=1

      case Function::PENALIZED_VOLUME:
      case Function::REALVOLUME:
      return IntegrateDesignVariable(f, g, derivative, des, normalized, false, func->GetParameter());

      case Function::GAP:
      if(!derivative)
      {
        double vol = IntegrateDesignVariable(f, g, false, des, normalized, false, 1.0);
        double pen = IntegrateDesignVariable(f, g, false, des, normalized, false, func->GetParameter());
        return vol - pen;
      }
      else
      {
        if(!design->IsRegular())
        throw Exception(Function::type.ToString(func->GetType()) + " only implemented for regular grids");

        CalcRegularGapConstraint(func, des);
        return -1.0;
      }

      default:
      assert(false);
      break;
    }
    return -1.0; // cannot happen due to assert
  }

  double ErsatzMaterial::CalcTrivialVolume(Function* f, bool derivative, bool normalized)
  {
    assert(!f->elements.IsEmpty());
    // In CalcOrthotropeMaterialProperties() we construct a dummy function, this has Function::elements not set :(
    // only for physical
    // TODO: assumes a single transfer function for all regions!
    // TODO: App::MECH ist stupid when we do App::LBM
    TransferFunction* tf = f->IsPhysical() ? design->GetTransferFunction(f->GetDesignType(), App::MECH) : NULL;
    bool regular = design->IsRegular();
    unsigned int numEls = f->elements.GetSize();
    unsigned int base;
    if(design->HasMultiMaterial())
      base = 0;
    else
      base = design->FindDesign(DesignElement::DENSITY)*numEls;
    double sum = 0.0;
    // we need the total volume in the non-regular case
    double total_vol = 0.0;
    if(!normalized)
    {
      total_vol = 1.0;
    }
    else
    {
      if (!regular)
      {
        for (unsigned int i = base; i < base+numEls; i++)
          total_vol += f->elements[i]->CalcVolume();
      }
      else
      total_vol = numEls;
    }
    // in the multimaterial case we have to consider, the multiple design element case
    if(design->HasMultiMaterial())
      total_vol /= design->GetMultiMaterials().GetSize();

    LOG_DBG(em) << "CTV: d=" << derivative << " p=" << f->IsPhysical() << " n=" << normalized << " tv=" << total_vol << " ex=" << context->GetExcitation()->GetFullLabel();
    for (unsigned int i = base; i < base+numEls; i++)
    {
      DesignElement* de = f->elements[i];

      double val = f->IsPhysical() ? tf->Transform(de, DesignElement::SMART) : de->GetDesign(DesignElement::SMART);
      double vol = (regular ? 1.0 : de->CalcVolume())/total_vol;
      sum += vol * val;
      if(derivative)
        de->AddGradient(f, vol);

      LOG_DBG2(em) << "CTV de=" << de->elem->elemNum << " val=" << val << " vol=" << vol << " -> " << vol*val;
    }
    return derivative ? -1.0 : sum;
  }



  double ErsatzMaterial::CalcDesignTracking(Condition* g, bool derivative)
  {
    assert(g->elements.GetSize() > 0 && g->elements.GetSize() == g->pattern.GetSize());
    // g = 1/N * sum (rho - rho^*)^2 where rho is the physical rho
    // g' = 2/N * (rho - rho^*) * rho'  and the derivative of the filter if any
    int res_idx = design->GetSpecialResultIndex(DesignElement::DEFAULT, DesignElement::DESIGN_TRACKING);
    double result = 0.0;
    TransferFunction* tf = design->GetTransferFunction(ToDesign(g->ctxt->pde), g->ctxt->ToApp());
    for (unsigned int i = 0, n = g->elements.GetSize();i < n;i++)
    {
      DesignElement* de = g->elements[i];
      double rho = tf->Transform(de, DesignElement::SMART); // physical design
      double target = g->pattern[i];
      if (derivative)
      {
        // for the non-element designs nothing is added and it stays 0
        de->AddGradient(g, (2.0 / n) * (rho - target) * tf->Derivative(de, DesignElement::SMART));
      }
      else
      {
        double val = (rho - target) * (rho - target);
        result += val;
        if (res_idx > 0)
        de->specialResult[res_idx] = val;
      }
    }

    return result / g->elements.GetSize();
  }

  void ErsatzMaterial::CalcRegularGapConstraint(Function* f, DesignElement::Type dt)
  {
    assert(design->IsRegular());
    unsigned int des = design->FindDesign(dt);
    unsigned int ele = design->GetNumberOfElements();
    // exponent for penalized volume
    const double exp = f->GetParameter();
    for (unsigned int i = des * ele;i < (des + 1) * ele;i++)
    {
      DesignElement& de = design->data[i];
      // derivative of vol is 1
      // derivative of pen (x^p) is p*x^(p-1)
      double pen_grad = exp * std::pow(de.GetDesign(DesignElement::SMART), exp - 1.0);
      // normalize
      double grad = (1.0 - pen_grad) / (double)(ele);
      de.AddGradient(f, grad);
    }
  }

  double ErsatzMaterial::CalcGlobalDynamicCompliance(Excitation& excite, Function* f)
  {
    //GLOBAL_DYNAMIC_COMPLIANCE
    // c = u^T conj(u) -> "A note on sensitivity analysis of linear dynamic systems with
    //                          harmonic excitation"; Jakob S. Jensen; June 22, 2007
    Vector<complex<double> >& u = forward.Get(excite)->GetComplexVector(StateSolution::RAW_VECTOR);
    assert(u.GetSize() != 0);
    complex<double> csp;
    u.Inner(u, csp);// the inner product is sum over u_i * conj(u_i);
    assert(csp.imag() == 0);
    return csp.real() * excite.GetFactor(f);
  }

  double ErsatzMaterial::CalcCompliance(Excitation& excite, Objective* f, Condition* g, bool derivative)
  {
    // note for the derivative case gradients are summed up (with weights)
    assert(f != NULL || g != NULL);
    assert(f == NULL || g == NULL);
    Function* func = Function::GetFunction(f, g);
    double result = 0.0;
    if(derivative)
    {
      // calculate the compliance which is according to
      // "A 99 line topology optimization code written in Matlab"; O.Sigmund, 2001
      // -> dc/dx_e = -p * x_e ^(p-1) u_e^T k_0 u_e
      TransferFunction* tf = design->GetTransferFunction(func->GetDesignType() , App::MECH, true);
      double factor = excite.GetWeightedFactor(func);

      if(IsTransient())
      {
        // this computes the complete derivative of the Newmark scheme, up to now, all objectives/constraints can be handled like that
        // as the derivative of all objectives/constraint is calculated p^T (dF - dA) u
        // where p is solution of adjoint, dF is derivative of newmark update, dA is derivative of system matrix, u is solution of forward problem
        CalcNewmarkDerivative(excite, forward, adjoint, factor, f, g);
      }
      else
      {
        CalcU1KU2(tf, forward.Get(excite)->elem[App::MECH], App::MECH, forward.Get(excite)->elem[App::MECH], NULL, -factor, STANDARD, func);
      }
    }
    else // now comes not derivative
    {
      UInt timesteps = func->ctxt->GetDriver()->GetNumSteps();
      result = 0.0;
      for(unsigned int ts = 0; ts < timesteps; ++ts)
      {
        // this formulation works for transient as well as static cases, integral over time
        // compliance is easier computed using f^T u on nodes with force
        // to avoid any work for assembling force again, we simply calculate solution times rhs from the system

        // when not transient or eigenvalue we dont't store timestep_mode as the default 0 is valid there as timestep-nr/mode-nr.
        int corr_ts = timesteps == 1 ? -1 : (int) ts;
        LOG_DBG2(em) << "CC: ts=" << ts << " corr_ts=" << corr_ts;
        Vector<double>& u = forward.Get(excite, NULL, corr_ts)->GetRealVector(StateSolution::RAW_VECTOR);
        Vector<double>& rhs = forward.Get(excite, NULL, corr_ts)->GetRealVector(StateSolution::RHS_VECTOR);
        double sp = 0.0;
        u.Inner(rhs, sp);
        result += sp * excite.GetFactor(func) * GetStepWeight(ts);
        LOG_DBG(em) << "CC: result=" << result << " sp=" << sp << " u=" << u.ToString() << " func=" << func->ToString();
      }
    }
    return result;
  }

  double ErsatzMaterial::CalcHeatEnergy(Excitation& excite, Objective* c, Condition* g, bool derivative)
  {
    assert(c != NULL || g != NULL);
    assert(c == NULL || g == NULL);
    Function* f = Function::GetFunction(c, g);

    double res = 0.0;

    if(derivative)
    {
      TransferFunction* tf = design->GetTransferFunction(f->GetDesignType() , App::HEAT, true);

      if (!interfaceDrivenGradCalc_) {
        CalcAndStoreInterfaceDrivenGrad<double>(f,tf);
        interfaceDrivenGradCalc_ = true;
      }

      double factor = excite.GetWeightedFactor(f);
      HeatPDE* heat = dynamic_cast<HeatPDE*>(f->ctxt->pde);
      assert(heat != NULL);
      DesignDependentRHS* rhs = NULL;
      if (heat->HasInterfaceDrivenRHS())
      {
        rhs = new DesignDependentRHS();
        rhs->Init<double>(design,App::HEAT);
        // f'^Tu de->AddGradient(f, this_value);
        StdVector<SingleVector*>& stateSol = forward.Get(excite)->elem[App::HEAT];
        for (unsigned int id = 0; id < design->data.GetSize(); id++) {
          DesignElement* de = &design->data[id];
          Vector<double> gradRHS = de->interfaceDrivenLoadGrad_; // f'
          double val = gradRHS.Inner(*stateSol[id]);
          de->AddGradient(f,val);
        }
      }
      CalcU1KU2(tf, forward.Get(excite)->elem[App::HEAT], App::HEAT, forward.Get(excite)->elem[App::HEAT], rhs, -factor, STANDARD, f);
    }
    else {
      Vector<double>& u = forward.Get(excite, NULL)->GetRealVector(StateSolution::RAW_VECTOR);
      Vector<double>& rhs = forward.Get(excite, NULL)->GetRealVector(StateSolution::RHS_VECTOR);
      u.Inner(rhs,res);
      res *= excite.GetFactor(f);
    }
    return res;
  }

  template<class T>
  double ErsatzMaterial::CalcOutput(Excitation& excite, Function* f)
  {
    // The output is <l,u> where l is -1* rhs of the adjoint pde
    // Here the rhs of the adoint pde is not -1 -> hence we do -1*<l,u>
    //
    // We always take l from rhs and don't store it explicitly.
    // Note that one has to use the algsys RHS! The PDE RHS is still from the
    // forward simulation!
    Vector<T>& u = dynamic_cast<Vector<T> & >(*(forward.Get(excite, NULL)->GetVector(StateSolution::RAW_VECTOR)));
    Vector<T>& l = dynamic_cast<Vector<T>&>(*(adjoint.Get(excite, f)->GetVector(StateSolution::SEL_VECTOR)));

    // temporary vector for displacement used for squared_output
    Vector<T> u_square(u.GetSize());

    assert(u.GetSize() == l.GetSize());
    LOG_DBG2(em) << "CO: f=o: " << f->IsObjective() << " adjoint sel (l): " << l.ToString(1);
    LOG_DBG2(em) << "CO: forward sol (u): " << u.ToString(0);
    double result = 0.0;
    switch(f->GetType())
    {
      case Function::OUTPUT:
      case Function::SQUARED_OUTPUT:
      {
        T inner;
        if (f->GetType() == Objective::SQUARED_OUTPUT) {
          // this is <l, u o u>
          for(unsigned int i = 0;i<u.GetSize();i++) {
            u_square[i] = u[i] * u[i];
            //l[i] *= l[i];
          }
          inner = u_square.Inner(l);
        } else {
          // this is <l, u> which is for complex not really defined as it might be non-real
          inner = u.Inner(l);
        }
        result = ((complex<double>) inner).real();
        //if (f->GetType() == Objective::SQUARED_OUTPUT)
        //  result *= result;
        result *= excite.GetFactor(f);
        LOG_DBG2(em) << "CO: <l,u>: " << inner << " * " << excite.GetFactor(f) << " -> " << result;
        break;
      }

      case Function::ABS_OUTPUT:
      {
        // |<u,l>|
        T ul = u.Inner(l);
        result = std::abs(ul);
        LOG_DBG2(em) << "CO: |<u,l>| = |" << ul << "| -> " << result;
        break;
      }

      case Objective::DYNAMIC_OUTPUT:
      case Objective::CONJUGATE_COMPLIANCE:
      {
        // this is <u,L conj(u)> and only defined for the harmonic case!
        if(!f->ctxt->IsComplex()) throw Exception("'" + f->type.ToString(f->GetType()) + "' is only defined for harmonic!");

        // we loop over the vectors and do the scalar product by hand as we have
        // no diagonal matrix version of l
        for(unsigned int i = 0; i < u.GetSize(); i++ )
        {
          if(l[i] == 0.0) continue; // we skip this so we can make output for the real cases

          complex<double> u_val = (complex<double>) u[i];
          // make sure we have no penalization stuff!
          assert(std::abs(u_val) < 1e15);

          double sp = std::real(std::conj(u_val) * l[i] * u_val);
          LOG_DBG2(em) << "CO: CalcObjective: " << std::conj(u_val) << " * " << l[i] << " * " << u_val << " -> " << sp;
          result += sp;
        }
        LOG_DBG2(em) << "CO: <u,L u*>: " << result << " * " << excite.GetFactor(f) << " -> " << result * excite.GetFactor(f);
        result *= excite.GetFactor(f);
        break;
      }

      default: EXCEPTION("Not handled");
      break;
    }
    return result;
  }


  /* tracking and transient
  void ErsatzMaterial::SetAdjointRhs(AdjointParameters* adjointParams)
  {
    Assemble* assemble = context->pde->GetAssemble();
    int ts = context->GetDriver()->GetActStep("mech") - 1; // drivers count timesteps starting with 1
    Excitation& excite = *(adjointParams->GetExcitation());
    switch(adjointParams->GetFunction()->GetType())
    {
      case Objective::TRACKING:
      SetTrackingAdjointRhs(excite, ts);
      break;
      case Objective::COMPLIANCE: // adjoint has the original rhs (scaled with weight per timestep), but time walks backwards
      {
        assemble->AssembleLinRHS();
        Vector<Double> rhs;
        assert(false);
        // FIXME assemble->GetAlgSys()->GetRHSVal(rhs);
        double w = GetStepWeight(ts);
        for(unsigned int i = 0; i < rhs.GetSize(); ++i)
        {
          rhs[i] *= w;
        }
        assert(false);
        // FIXME assemble->GetAlgSys()->InitRHS(rhs);
      }
      break;

      default:
      EXCEPTION("adjoint RHS not known for Objective " + adjointParams->GetFunction()->ToString());
    }
    // if the first step is static, we have to readjust the right hand side for ts == 0. Note that the transientDriver/solveStep does no rhs update any more
    if (IsFirstTransientStepStatic() && ts == 0)
    {
      assert(false);
       FIXME FE-Space
      double dt = dynamic_cast<TransientDriver*>(domain->GetDriver())->GetDeltaT();
      double gamma =pde->getTimeStepping()->GetNewmarkGamma();
      double beta = pde->getTimeStepping()->GetNewmarkBeta();

      Vector<Double>& pp = pde->getTimeStepping()->GetDeriveMap()[FIRST_DERIV];
      Vector<Double>& ppp = pde->getTimeStepping()->GetDeriveMap()[SECOND_DERIV];

      // note, that these point eveluations are divided by dt, as our integral is missing multiplication by dt, it is normed by the number of timesteps using GetStepWeight
      Vector<Double> coeffMass = pp;
      coeffMass.ScalarMult(1.0 / dt);
      coeffMass.Add(1.0 - gamma, ppp);
      assemble->GetAlgSys()->UpdateRHS(CoupledField::App::MASS, coeffMass);

      // look up, whether the damping matrix exists
      std::set<FEMatrixType> matTypes;
      assemble->GetAlgSys()->GetFEMatrixTypes(matTypes);
      if(matTypes.find(CoupledField::DAMPING) != matTypes.end()){
        Vector<Double> coeffDamping(0);
        assemble->GetAlgSys()->GetSolutionVal(coeffDamping);
        coeffDamping.ScalarMult(1.0 / dt);
        coeffDamping.Add(0.5, pp);
        coeffDamping.Add(0.5 * (gamma - 2*beta) * dt, ppp);
        assemble->GetAlgSys()->UpdateRHS(CoupledField::DAMPING, coeffDamping);
      }

    }

    // in case of contact, we have to inform the solver, that an adjoint system is solved
    assert(false);
    // FIXME assemble->GetAlgSys()->PrepareForAdjoint(forward.Get(excite, NULL, ts)->GetRealVector(Solution::RAW_VECTOR));
  }

  */

  double ErsatzMaterial::CalcEnergyFlux(Excitation& excite, Objective* f)
  {
    // calc 1/2 Re{j omega psi^R Q psi^*} where Q is the grad_n matrix B applied at some points L
    // this is the global element solution, indexed by equation number
    Vector<complex<double> > u_glob = forward.Get(excite)->GetComplexVector(StateSolution::RAW_VECTOR);
    // here we store Q*u^* as we have to determine the nodal entries by count
    Vector<complex<double> > q_u_glob(u_glob.GetSize());
    SetEnergyFluxVector(f, u_glob, false, q_u_glob);
    // calculate the product Re(j*omega*u*(Q*u^*)) which is -Im(u*(Q*u^*))
    // the (Q*u^*) term is normalized
    double result = 0;
    for (unsigned int i = 0, in = q_u_glob.GetSize();i < in;i++)
      result += std::real(excite.GetOmega() * complex<double>(0, 0.5) * u_glob[i] * q_u_glob[i]);
    return result;
  }

  void ErsatzMaterial::SetEnergyFluxVector(Function* f, const Vector<complex<double> >& u_glob, bool adjoint, Vector<complex<double> >& q_u_glob)
  {
    // determine the global vector Q*u^* or (Q - Q^T)^T*u^* in the adjoint case
    // Q is the grad_z vector (selected B matrix) applied to selection points L defined by
    // a surface region in the "energyFlux" element.
    // always determined again, no caching.
    // can be easily extended to other pdes!
    SinglePDE* mypde = f->ctxt->ToPDE(App::ACOUSTIC);
    shared_ptr<ResultInfo> res_info = mypde->GetResultInfo(ACOU_POTENTIAL);
    assert(false);
    shared_ptr<EqnMap> eqnMap; // FIMXE = mypde->GetEqnMap();
    // get the surface region description
    PtrParamNode srpn = f->pn->Get("output/energyFlux/surfRegion");
    // surface element list
    assert(false);
    shared_ptr<EntityList> sel; // FIXME = grid->GetEntityList(EntityList::SURF_ELEM_LIST, srpn->Get("name")->As<string>(), EntityList::REGION);
    // neighbor region is mandatory to define the volume elements
    assert(false);
    // FIXME RegionIdType vol_neigh = grid->GetRegion().Parse(srpn->Get("neighborRegion")->As<string>());
    // reset output
    q_u_glob.Init(complex<double>(0.0, 0.0));
    // here we count the entries to q_u_glob to normalize in the end
    StdVector<int> count(u_glob.GetSize());
    count.Init(0);
    // an element solution vector -> we need a 1 dof solution up to now!
    Vector<complex<double> > elem_sol;
    // node numbers common on a surface element and the matching volume element
    StdVector<unsigned int> common_nodes;
    // this contains an element B (grad_n) matrix to be applied with the solution.
    Matrix<complex<double> > q_mat;
    // surface normal
    Vector<double> normal;
    // traverse our surface elements
    EntityIterator it = sel->GetIterator();
    for(it.Begin(); !it.IsEnd(); it++)
    {
      assert(false);
      /* FIXME
      const SurfElem* se = it.GetSurfElem();
      const Elem* vol = se->ptVolElem1 != NULL && se->ptVolElem1->regionId == vol_neigh ? se->ptVolElem1 : se->ptVolElem2;
      assert(vol->regionId == vol_neigh);

      // determine grad_n u^* on the se nodes coinciding with the vol nodes
      FindCommonNodes(se, vol, common_nodes);

      // get the element solution
      ElemList se_it(grid);//single element iterator
      se_it.SetElement(vol);
      mypde->GetSolVecOfElement(elem_sol, se_it.GetIterator(), res_info);

      // determine selected element grad_n matrix (includes selection defined by surface elements)
      // the matrix is squared but for non-common nodes entries it is zero
      grid->CalcSurfNormal(normal, *se);
      LOG_DBG2(em) << "SEFV: se=" << se->elemNum << " normal=" << normal.ToString() << " relevant_direction=" << Point::GetCartesianOrientation(&normal) + 1;
      mypde->CalcElemGradMatrix(vol, common_nodes, Point::GetCartesianOrientation(&normal) + 1, vol->ptElem->GetAnsatzFct(), q_mat);// x, y, z-direction -> 1-based!! :(

      for(unsigned int n = 0; n < vol->connect.GetSize(); n++)
      {
        unsigned int node = vol->connect[n];

        for(unsigned int t = 0; t < vol->connect.GetSize(); t++)
        assert(common_nodes.Contains(node) || abs(q_mat[n][t]) < 1e-16);

        // sum up grad_N * u^*
        complex<double> sum = 0;
        for(unsigned idx = 0; idx < vol->connect.GetSize(); idx++)
        {
          if(adjoint)
          {
            //(Q -Q^T)^T = (Q^T - Q)
            sum += (q_mat[idx][n] - q_mat[n][idx]) * std::conj(elem_sol[idx]);
          }
          else
          {
            // simple Q u^*
            sum += q_mat[n][idx] * std::conj(elem_sol[idx]);
          }
        }

        int eqn_nr = eqnMap->GetNodeEqn(node,1); // ACOUSTIC!
        // we shall not be in the boundary conditions!
        assert(eqn_nr >= 0);
        int eqn_idx = eqn_nr -1;// FUCK!!!
        if(!close(sum, complex<double>(0.0, 0.0)))
        {
          q_u_glob[eqn_idx] += sum;
          count[eqn_idx]++;
        }
        LOG_DBG3(em) << "SEFV: vol=" << vol->elemNum << " node=" << node << " eqn_nr=" << eqn_nr << " count=" << count[eqn_idx] << " sum=" << sum;
      }
      */
    }
    LOG_DBG2(em) << "SEFV: q_u_glob=" << q_u_glob.ToString(1);
    LOG_DBG2(em) << "SEFV: count=" << count.ToString(1);
// normalize Q*u^*
    for (unsigned int i = 0, in = q_u_glob.GetSize();i < in;i++)
    if (count[i] != 0)
    q_u_glob[i] /= (double)(count[i]);
  }

  void ErsatzMaterial::FindCommonNodes(const SurfElem* se, const Elem* vol, StdVector<unsigned int>& common_nodes) const
  {
    const StdVector<unsigned int>& se_nodes = se->connect;
    const StdVector<unsigned int>& vol_nodes = vol->connect;
// for higher order elements the center node is not on the vol, for lin common_nodes = se nodes
    common_nodes.Reserve(se_nodes.GetSize());
    common_nodes.Clear();
    for (unsigned int i = 0;i < vol_nodes.GetSize();i++)
    {
      unsigned int n = vol_nodes[i];
      if (se_nodes.Contains(n))
        common_nodes.Push_back(n);
    }
    LOG_DBG3(em) << "FCN se=" << se->elemNum << " vol=" << vol->elemNum << " common=" << common_nodes.ToString();
  }

  /*
  void ErsatzMaterial::SetTrackingAdjointRhs(Excitation& excite, int ts)
  {
    assert(false);

    // this is for the static and for the transient case.
    Vector<double>& u = forward.Get(excite, NULL, ts)->GetRealVector(Solution::RAW_VECTOR);
    LOG_DBG3(em) << "SolveTrackingProblem: displacement vector: (" << u.GetSize() << ") " << u.ToString();
    assert(false);
    shared_ptr<EqnMap> eqnMap; // FIXME = pde->GetEqnMap();
    MathParser* parser = domain->GetMathParser();
    unsigned int mathParserHandle = parser->GetNewHandle();
    CoordSystem* coosy = domain->GetCoordSystem();
    Vector<Double> rhs;
    assert(false)
    // FIXME assemble->GetAlgSys()->GetRHSVal(rhs);
// set rhs to 0
    rhs.Init();
    double w = GetStepWeight(ts);
// set rhs for tracking nodes
    Vector<double> GlobCoord;
    for(unsigned int l=0; l < excite.trackings.GetSize(); l++)
    {
      TrackingBc const & actTrack = *(excite.trackings[l]);
      EntityIterator it = actTrack.entities->GetIterator();
      const int dof = actTrack.dof;
      for(it.Begin(); !it.IsEnd(); it++)
      {
        grid->GetNodeCoordinate(GlobCoord, it.GetNode());
        parser->SetCoordinates(mathParserHandle, *coosy, GlobCoord);
        parser->SetExpr(mathParserHandle, actTrack.value);
        const double uref = parser->Eval(mathParserHandle);
        const int eqnr = eqnMap->GetNodeEqn(it.GetNode(), dof) - 1; // equation numbers are 1 based but vector u is 0 based
        if(eqnr>=0)
        {
          double v = u[eqnr] - uref;
          rhs[eqnr] = v * w;
          LOG_DBG2(em) << "SolveTrackingProblem: tracking setting RHS equation " << eqnr << " (Node: " << it.GetNode() << ", dof: " << dof << ") to " << rhs[eqnr];
        }
      }
    }
    assemble->GetAlgSys()->InitRHS(rhs);
    parser->ReleaseHandle(mathParserHandle);
  }
  */


  double ErsatzMaterial::CalcExpression(Condition* g, bool derivative)
  {
    assert(g->GetBound() == Condition::ALPHA_MINUS_SLACK_VALUE || g->GetBound() == Condition::ALPHA_PLUS_SLACK_VALUE);

    // the comparison with alpha and slack  is done in BaseOptimizer(). We give back simply param

    // in the derivative case the slack derivatives are handled in AuxDesign().
    // The derivatives w.r.t standard design is zero, we simply don't add to the design gradient.

    return g->GetParameter(); // no harm to return on derivative
   }


  double ErsatzMaterial::CalcEigenfrequency(Excitation& org_excite, Function* f, bool derivative)
  {
    // for the bloch mode case this might be complex!
    // the eigenvalues lambda = (2*pi*ef)^2 !!
    Condition* g = dynamic_cast<Condition*>(f); // NULL when f is objective

    assert(f->GetEigenValueID() >= 1);
    unsigned int mode_idx = f->GetEigenValueID() - 1; // 0-based

    // each mode is encoded in forward as timestep_mode and in the bloch mode case the excitation idx is the wave number

    // we have the bloch=full case, then org_excite is exatly the wave number constraint
    // for bloch=extremal we need to search for the wave number (=excitation) whic is minimal/maximal

    double freq=-1.0;
    Excitation* ex = NULL;
    if(f->IsObjective() || g->DoFullBloch())
    {
      StdVector<double> efs = forward.CollectEigenfrequencies(org_excite);
      freq = efs[mode_idx];
      ex = &org_excite;
    }
    else
    {
      Context& ctxt = manager.GetContext(&org_excite);

      Matrix<double> mat = forward.CollectBlochEigenfrequencies(&ctxt);
      assert(g->GetBound() == Condition::LOWER_BOUND || g->GetBound() == Condition::UPPER_BOUND);
      // when we are lower bounded we search for the minimum. Also set freq
      SearchMinMax(mat, mode_idx, g->GetBound() == Condition::LOWER_BOUND, &freq, &(g->bloch));
      ex = ctxt.excitations[g->bloch.col]; // freq set above
      LOG_DBG2(em) << "CE: f=" << f->ToString() << " mat=" << mat.ToString(2);
      LOG_DBG(em) << "CE: mode_idx=" << mode_idx << " col_idx=" << g->bloch.col << " min=" << (g->GetBound() == Condition::LOWER_BOUND) << " f=" << freq;
    }
    LOG_DBG(em) << "CE: mode_idx=" << mode_idx << " f=" << freq;

    if(derivative)
    {
      // the sensitivity for a single modal eigenvalue is (Bendsoe & Sigmund (2.2) page 73)
      // d_ev = u^T (d_K - ev d_M) u with u the mode for the i-th eigenvalue and ev the i-th eigenvalue
      // hence this is the standard CalcU1KU2() scheme
      //
      // with ef = sqrt(ev)/(2*pi)
      // d_ef = 1/(4*pi)*ev^-0.5 * d_ev = 1/(8*pi^2*ef) * d_ev
      //
      // the modes are not stores via StorePDESolution() but held in the ArpackSolver
      TransferFunction* tf = design->GetTransferFunction(f->GetDesignType() , App::MECH, true);

      double factor = 1.0 / ( 8.0 * M_PI * M_PI * freq);
      // our eigenvalue
      double ev = std::pow(2.0 * M_PI * freq, 2);

      StateSolution* sol = forward.Get(*ex, NULL, mode_idx); // never give a function for forward!
      assert(sol->ContainsState());

      LOG_DBG2(em) << "CE mode_idx=" << mode_idx << " f=" << freq << " sol=" << sol->GetVector(StateSolution::RAW_VECTOR)->ToString();

      // we need to set the current wave_vector such that SetElementK determines the right stiffness matrices!
      if(f->ctxt->DoBloch())
        f->ctxt->GetEigenFrequencyDriver()->SetCurrentWaveVector(ex->GetWaveNumber());

      CalcU1KU2(tf, sol->elem[App::MECH], App::MECH, sol->elem[App::MECH], NULL, factor, EIGENFREQ, f, -1, ev);
    }
    return freq;
  }


  /** is a lot of copy and paste from CalcEigenfrquency :( */
  double ErsatzMaterial::CalcBandGap(Excitation& gap_excite, Function* f, bool derivative)
  {
    Context& ctxt = manager.GetContext(&gap_excite);
    Matrix<double> mat = forward.CollectBlochEigenfrequencies(&ctxt);

    double lower_freq = 0.0;
    double upper_freq = 0.0;
    f->bandgap.lower.col = SearchMinMax(mat, (unsigned int) f->bandgap.lower_ev - 1, false, &lower_freq); // maximum
    f->bandgap.upper.col = SearchMinMax(mat, (unsigned int) f->bandgap.upper_ev - 1, true, &upper_freq);  // minimum

    if(derivative)
    {
      // see CalcEigenfrequency()
      TransferFunction* tf = design->GetTransferFunction(f->GetDesignType() , App::MECH, true);

      double factor = 1.0 / ( 8.0 * M_PI * M_PI * upper_freq);
      StateSolution* sol_upper = forward.Get(*ctxt.excitations[f->bandgap.upper.col], NULL, f->bandgap.upper_ev - 1);
      f->ctxt->GetEigenFrequencyDriver()->SetCurrentWaveVector(ctxt.excitations[f->bandgap.upper.col]->GetWaveNumber());
      CalcU1KU2(tf, sol_upper->elem[App::MECH], App::MECH, sol_upper->elem[App::MECH], NULL, factor, EIGENFREQ, f, -1, std::pow(2.0 * M_PI * upper_freq, 2));

      factor = -1.0 / ( 8.0 * M_PI * M_PI * lower_freq); // we substract!!
      StateSolution* sol_lower = forward.Get(*ctxt.excitations[f->bandgap.lower.col], NULL, f->bandgap.lower_ev - 1);
      f->ctxt->GetEigenFrequencyDriver()->SetCurrentWaveVector(ctxt.excitations[f->bandgap.lower.col]->GetWaveNumber());
      CalcU1KU2(tf, sol_lower->elem[App::MECH], App::MECH, sol_lower->elem[App::MECH], NULL, factor, EIGENFREQ, f, -1, std::pow(2.0 * M_PI * lower_freq, 2));
    }
    return upper_freq - lower_freq;
  }


  double ErsatzMaterial::CalcStateTrackingAtInterface(Excitation& excite, Function* f, bool derivative, double trackVal)
  {
    assert(Context::ToApp(f->ctxt->pde) == App::HEAT);
    assert(f->GetType() == Condition::TEMP_TRACKING_AT_INTERFACE || f->GetType() == Objective::TEMP_TRACKING_AT_INTERFACE);
    trackingFunc_ = f;
    double res = 0.0;

    double sourceVal = 0.0;
    f->ctxt->pde->GetParamNode()->GetValue("bcsAndLoads/designDependentHeatSource/value",sourceVal);

    if (derivative)
    { // (u - u_)^T * F'(u - u_), where u_ is tracked temperature and F diag(f)
      TransferFunction* tf = design->GetTransferFunction(f->GetDesignType() , App::HEAT, true);
      double factor = excite.GetWeightedFactor(f);
      HeatPDE* heat = dynamic_cast<HeatPDE*>(f->ctxt->pde);
      assert(heat != NULL);
      DesignDependentRHS* rhs = NULL;
      if (heat->HasInterfaceDrivenRHS())
      {
        rhs = new DesignDependentRHS();
        rhs->Init<double>(design,App::HEAT);
        StdVector<SingleVector* >& all_u_elem = forward.Get(excite)->elem[App::HEAT];

        if (!interfaceDrivenGradCalc_) {
          CalcAndStoreInterfaceDrivenGrad<double>(f,tf);
          interfaceDrivenGradCalc_ = true;
        }

        for (unsigned int e = 0; e < design->data.GetSize(); e++) { // (u_i - u_ref)^2, element based
          Vector<double>& u_elem = dynamic_cast<Vector<double>& >(*(all_u_elem[e]));
          DesignElement* de = &design->data[e];
          Vector<double> gradLoad = de->interfaceDrivenLoadGrad_; // f'
          double val = 0.0;
          // f'_i * (u_i - u_track)^2
          for (unsigned int n = 0; n < gradLoad.GetSize(); n++)
            val += gradLoad[n] * (u_elem[n] - trackVal) * (u_elem[n] - trackVal);

//          de->AddGradient(f,val*design->data.GetSize() / sourceVal / sqrt(domain->GetGrid()->GetNumElems()));
          de->AddGradient(f,val*design->data.GetSize() / sourceVal);
        }
      }
      CalcU1KU2(tf, adjoint.Get(excite,f)->elem[App::HEAT], App::HEAT, forward.Get(excite)->elem[App::HEAT], rhs, factor, STANDARD, f);
    }
    else
    {
      StdVector<unsigned int> nodeList;
      domain->GetGrid()->GetNodesByRegion(nodeList,design->GetRegionId());

      for (unsigned int i = 0; i < nodeList.GetSize(); i++)
//        res += CalcStateTrackingAtNode(nodeList[i]) / sqrt(domain->GetGrid()->GetNumElems());
        res += CalcStateTrackingAtNode(nodeList[i]);

    } // if-else

    return res;
  }

  double ErsatzMaterial::CalcSlackFunction(Function* f, bool derivative)
  {
    assert(f->GetType() == Function::SLACK_FNCT);
    assert(f->GetSlackFnct() != Function::NO_FUNCTION);
    assert(design->HasAlphaVariable() && design->HasSlackVariable()); // shall be checked already

    double a = design->GetAlphaVariable();
    double s = design->GetSlackVariable();

    double result = 0.0;

    if(!derivative)
    {
      switch(f->GetSlackFnct())
      {
      case Function::ALPHA_SLACK_QUOTIENT:
        if(IsNoise(a) && IsNoise(s))
          optInfoNode->SetWarning("'alpha' and 'slack' are both close to zero. Best adjust your bounds to avoid division by zero.");
        result = a / s;
        break;

      case Function::REL_BANDGAP:
        if(IsNoise(a-s))
          optInfoNode->SetWarning("denominator of '" + Function::slackFnct.ToString(Function::REL_BANDGAP) + "' is close to zero. Adjust bounds or use " + Function::slackFnct.ToString(Function::NORM_BANDGAP));
        assert(std::abs(a-s) > 1e-8);
        result = (2*s)/(a-s);
        break;

      case Function::NORM_BANDGAP:
        if(IsNoise(a))
          optInfoNode->SetWarning("'alpha' is close to zero for function " + Function::slackFnct.ToString(Function::NORM_BANDGAP));
        assert(std::abs(a) > 1e-8);
        result = (2*s)/a;
        break;

      case Function::ALPHA_MINUS_SLACK:
        result = a-s;
        break;

      case Function::NO_FUNCTION:
        assert(false);
        break;
      }
    }
    else // derivative case
    {
      AuxDesign* ad = dynamic_cast<AuxDesign*>(design);
      assert(ad != NULL);

      double da = 0;
      double ds = 0;

      switch(f->GetSlackFnct())
      {
      case Function::ALPHA_SLACK_QUOTIENT:
        da = 1/s ;
        ds = -a/(s*s);
        break;

      case Function::REL_BANDGAP:
        da = -2*s / ((a-s)*(a-s));
        ds = 2*s / ((a-s)*(a-s)) + 2/(a-s);
        break;

      case Function::NORM_BANDGAP:
        da = -2*s / (a*a);
        ds = 2/a;
        break;

      case Function::ALPHA_MINUS_SLACK:
        da = 1;
        ds = -1;
        break;

      case Function::NO_FUNCTION:
        assert(false);
        break;
      }

      ad->GetAlphaDesign()->AddGradient(f, da);
      ad->GetSlackDesign()->AddGradient(f, ds);
    }

    return result;
  }


  double ErsatzMaterial::CalcStateTrackingAtNode(int node)
  {
    assert(node > 0);
    assert(trackingFunc_ != NULL);

    if (trackingFunc_ == NULL) // in case tracking result should be written out, but tracking is actually not an active function
      return 0.0;

    NodeList nodeList(domain->GetGrid());
    StdVector<UInt> nodeId(1);
    nodeId[0] = node;
    nodeList.SetNodes(nodeId);

    shared_ptr<BaseFeFunction> fe = trackingFunc_->ctxt->pde->GetFeFunction(trackingFunc_->ctxt->pde->GetNativeSolutionType());

    Vector<double> stateSol(1); // we get one scalar
    fe->GetEntitySolution(stateSol,nodeList.GetIterator()); // state solution at node 'node'

    shared_ptr<BaseFeFunction> rhsFe = trackingFunc_->ctxt->pde->GetRhsFeFunctions()[HEAT_TEMPERATURE];
    Vector<double> load(1); // scalar
    rhsFe->GetEntitySolution(load,nodeList.GetIterator()); // load at node 'node'

    double trackVal = trackingFunc_->GetParameter();
    double factor = 0.0;
    trackingFunc_->ctxt->pde->GetParamNode()->GetValue("bcsAndLoads/designDependentHeatSource/value",factor);

    LOG_DBG3(em) << "CSTAN node=" << node << " u=" << stateSol[0];
    assert(stateSol[0] > -1e15);
    // assert((stateSol[0] - trackVal) * (stateSol[0] - trackVal));
    assert(factor > 0);

    return load[0] * (stateSol[0] - trackVal) * (stateSol[0] - trackVal) * design->data.GetSize() / factor;
  }

  double ErsatzMaterial::CalcTempAtInterface(int node)
  {
    if (trackingFunc_ == NULL)
      return 0.0;

    // tempAtInterface(node) = load(node) * temp(node), here load is normed to 1
    NodeList nodeList(domain->GetGrid());
    StdVector<UInt> nodeId(1);
    nodeId[0] = node;
    nodeList.SetNodes(nodeId);

    shared_ptr<BaseFeFunction> fe = trackingFunc_->ctxt->pde->GetFeFunction(trackingFunc_->ctxt->pde->GetNativeSolutionType());

    Vector<double> stateSol(1); // we get one scalar
    fe->GetEntitySolution(stateSol,nodeList.GetIterator()); // state solution at node 'node'

    shared_ptr<BaseFeFunction> rhsFe = trackingFunc_->ctxt->pde->GetRhsFeFunctions()[HEAT_TEMPERATURE];
    Vector<double> load(1);
    rhsFe->GetEntitySolution(load,nodeList.GetIterator()); // load at node 'node'

    double factor = 0.0;
    trackingFunc_->ctxt->pde->GetParamNode()->GetValue("bcsAndLoads/designDependentHeatSource/value",factor);

    return load[0] * stateSol[0] * design->data.GetSize() / factor;
  }

  void ErsatzMaterial::CalcAdjointRHSStateTracking(Excitation& excite, Function* f, double trackVal, Vector<double>& out)
  {
    Vector<double> stateSol = forward.Get(excite,NULL)->GetRealVector(StateSolution::RAW_VECTOR);
    out.Resize(stateSol.GetSize(),0.0);

    Vector<double> loads = forward.Get(excite, NULL)->GetRealVector(StateSolution::RHS_VECTOR);
    double factor = 0.0;
    f->ctxt->pde->GetParamNode()->GetValue("bcsAndLoads/designDependentHeatSource/value",factor);

    for (unsigned int i = 0; i < stateSol.GetSize(); i++)
      out[i] = - 2.0 * loads[i] * (stateSol[i] - trackVal) * design->data.GetSize() / factor;
  }


  double ErsatzMaterial::CalcTracking(Excitation& excite, Objective* c, Condition* g, bool derivative)
  {
    /* FIXME
    Function* f = Function::Cast(c, g);
    UInt timesteps = context->GetDriver()->GetNumSteps();
    if(derivative)
    {
      // calculate the tracking functional gradient, which is z^T k_i u,
      // where Kz = ut
      // where ut = M^T M (u-u0) = I_\Gamma (u-u0)

      // calculate gradient z^T k_i u
      // note that in multiple excitations case (this is always now), we do sum this up
      TransferFunction* tf = design->GetTransferFunction(DesignElement::DENSITY, App::MECH, false);
      double factor = excite.GetWeightedFactor(c);

      if(IsTransient())
      {
        // this computes the complete derivative of the Newmark scheme, up to now, all objectives/constraints can be handled like that
        // as the derivative of all objectives/constraint is calculated p^T (dF - dA) u
        // where p is solution of adjoint, dF is derivative of newmark update, dA is derivative of system matrix, u is solution of forward problem
        CalcNewmarkDerivative(excite, forward, adjoint, factor, c, g);
      }
      else
      {
        CalcU1KU2(tf, adjoint.Get(excite, f)->elem[App::MECH], App::MECH, forward.Get(excite)->elem[App::MECH], NULL, -factor, STANDARD, f);
      }
      return 0.0;
    }
    else
    {
      // prepare for calculation
      shared_ptr<EqnMap> eqnMap = pde->GetEqnMap();
      MathParser * parser = domain->GetMathParser();
      unsigned int mathParserHandle = parser->GetNewHandle();
      CoordSystem * coosy = domain->GetCoordSystem();

      // the tracking functional is ut^T ut (ut as above)
      double val = 0.0;
      Vector<double> GlobCoord;

      double dt = 0.0;
      if(IsTransient())
      {
        dt = dynamic_cast<TransientDriver*>(domain->GetDriver())->GetDeltaT();
      }
      parser->SetValue(MathParser::GLOB_HANDLER, "dt", dt);

      for(unsigned int ts = 0; ts < timesteps; ++ts)
      { // this formulation works for transient as well as static cases, integral over time
        Vector<double>& u = forward.Get(excite, NULL, ts)->GetRealVector(Solution::RAW_VECTOR);// Tracking is only implemented for non-harmonic
        LOG_DBG3(em) << "CalcTracking: displacement vector: (" << u.GetSize() << ") " << u.ToString();
        parser->SetValue(MathParser::GLOB_HANDLER, "t", dt*(ts+1));
        parser->SetValue(MathParser::GLOB_HANDLER, "step", ts+1);

        double v = 0.0;
        for(unsigned int l=0; l < excite.trackings.GetSize(); l++)
        {
          TrackingBc const & actTrack = *(excite.trackings[l]);
          EntityIterator it = actTrack.entities->GetIterator();
          const int dof = actTrack.dof;
          for(it.Begin(); !it.IsEnd(); it++)
          {
            grid->GetNodeCoordinate(GlobCoord, it.GetNode(), true); // get Updated Coordinate
            parser->SetCoordinates(mathParserHandle, *coosy, GlobCoord);
            parser->SetExpr(mathParserHandle, actTrack.value);
            const double uref = parser->Eval(mathParserHandle);
            const int eqnr = eqnMap->GetNodeEqn(it.GetNode(), dof) - 1;// equation numbers are 1 based, the vector u is 0 based
            const double uact = eqnr>=0 ? u[eqnr] : 0.0;
            v += (uact - uref) * (uact - uref);
          }
        }
        val += v * GetStepWeight(ts);
      }
      parser->ReleaseHandle(mathParserHandle);
      return 0.5 * val;
    }
    */
    return -1;
  }

  double ErsatzMaterial::CalcPoissonsRatioAndYoungsModulus(Function* f, bool derivative)
  {
    Function::Type ft = f->GetType();
    assert(ft == f->POISSONS_RATIO || ft == f->YOUNGS_MODULUS || ft == f->YOUNGS_MODULUS_E1 || ft == f->YOUNGS_MODULUS_E2);
    SubTensorType stt = f->ctxt->stt;
    assert(stt == PLANE_STRAIN || stt == PLANE_STRESS || stt == FULL);
    Matrix<double> hom_tensor = CalcHomogenizedTensor(f);
    LOG_DBG(em) << "CPRAYM der=" << derivative << " ht=" << hom_tensor.ToString();
    const double E11 = hom_tensor[1 - 1][1 - 1];
    const double E12 = hom_tensor[1 - 1][2 - 1];
    const double E22 = hom_tensor[2 - 1][2 - 1];
    double result = 0.0;
    if(derivative)
    {
      // for iso-orthotropy we use the formulas for isotropy
      // Poisson's Ratio:
      // see MechanicMaterial::CalcIsotropicPoissonsRatio()
      // FULL + PLANE_STRAIN: v = E12 / (E11 + E12)
      // PLANE_STRESS: v = E12 / E11
      //
      // Young's Modulus
      // see MechanicMaterial::CalcIsotropicYoungsModulus()
      // FULL + PLANE_STRAIN: E = E11 * (1+v) * (1-2v)/(1-v)
      // PLANE_STRESS: E = E11 * (1-v^2)

      StdVector<double> dE11;
      CalcHomogenizedTensorEntry(f->ctxt, boost::make_tuple(1,1,1.0), true, dE11, f->GetExcitation()->meta_index);
      StdVector<double> dE12;
      CalcHomogenizedTensorEntry(f->ctxt, boost::make_tuple(1,2,1.0), true, dE12, f->GetExcitation()->meta_index);
      StdVector<double> dE22;
      CalcHomogenizedTensorEntry(f->ctxt, boost::make_tuple(2,2,1.0), true, dE22, f->GetExcitation()->meta_index);

      double grad(0.0);

      Transform* trans = f->GetExcitation()->transform; // transform might be NULL

      for(unsigned int o = 0, ne = design->GetNumberOfElements(); o < ne; o++)
      {
        // in case of transformation, the state is already transformed for the forward simulation
        // we need the transformation for the design and to store the result transformed
        
        DesignElement* de = design->ApplyTransformations(&design->data[o], &design->data[o], trans);
        unsigned int e = de->GetIndex();

        if(ft== f->POISSONS_RATIO && (stt == FULL || stt == PLANE_STRAIN))
        {
          grad = (dE12[e] * E11 - E12 * dE11[e]) / ((E11 + E12) * (E11 + E12));
        }
        if(ft == f->POISSONS_RATIO && stt == PLANE_STRESS)
        {
          grad = (dE12[e] * E11 - E12 * dE11[e]) / (E11 * E11);
        }
        if(ft == f->YOUNGS_MODULUS && (stt == FULL || stt == PLANE_STRAIN))
        {
          grad = (E11 * E11 + 2.0 * E11 * E12 + 3.0 * E12 * E12) * dE11[e];
          grad -= (4.0 * E11 * E12 + 2.0 * E12 * E12) * dE12[e];
          grad /= (E11 + E12) * (E11 + E12);
        }
        if(ft == f->YOUNGS_MODULUS && stt == PLANE_STRESS)
        {
          grad = (E11 * E11 + E12 * E12) * dE11[e];
          grad -= 2.0 * E11 * E12 * dE12[e];
          grad /= E11 * E11;
        }
        if((ft == f->YOUNGS_MODULUS_E1 || ft == f->YOUNGS_MODULUS_E2) && stt == PLANE_STRESS)
        {
          double t1 = dE11[e] * E22 + E11 * dE22[e] - 2.0 * E12 * dE12[e];
          double t2 = E11 * E22 - E12 * E12;

          if(ft == f->YOUNGS_MODULUS_E1)
          grad = (t1 * E22 - t2 * dE22[e]) / (E22 * E22);
          else
          grad = (t1 * E11 - t2 * dE11[e]) / (E11 * E11);
        }
        if((ft == f->YOUNGS_MODULUS_E1 || ft == f->YOUNGS_MODULUS_E2) && (stt == FULL || stt == PLANE_STRAIN))
        {
          throw Exception("youngsModulusE1/2 only implemented for plane stress");
        }
        LOG_DBG2(em) << "CPRAYM f=" << f->ToString() << " deriv o=" << o << " e=" << e << " de=" << de->ToString() << " gr=" << grad;
        //design->data[o].AddGradient(f, grad);
        de->AddGradient(f, grad);
      }
    }
    else
    {
      switch(ft)
      {
        case Function::POISSONS_RATIO:
        result = MechanicMaterial::CalcIsotropicPoissonsRatio(hom_tensor, stt);
        break;
        case Function::YOUNGS_MODULUS:
        result = MechanicMaterial::CalcIsotropicYoungsModulus(hom_tensor, stt);
        break;
        case Function::YOUNGS_MODULUS_E1:
        result = (E11 * E22 - E12 * E12) / E22;
        break;
        case Function::YOUNGS_MODULUS_E2:
        result = (E11 * E22 - E12 * E12) / E11;
        break;
        default:
        assert(false);
      }
      LOG_DBG(em) << "CPRAYM f=" << f->ToString() << " r=" << result << " ht=" << hom_tensor.ToString();
    }
    return result;
  }

  void ErsatzMaterial::SetTestStrainMatrix(Matrix<double>& matrix, const Vector<double>& vec)
  {
    assert(matrix.GetNumCols() == dim);
    assert(matrix.GetNumRows() == dim);
    matrix[0][0] = vec[0];
    matrix[1][1] = vec[1];
    matrix[0][1] = vec[5]; // Voigt notation!
    matrix[1][0] = 0.0;// because of symmetry we need factor 0.5 : FIXME why not vec[5]!!
    if (dim == 3)
    {
      matrix[2][2] = vec[2];
      matrix[1][2] = vec[3];
      matrix[2][1] = 0.0; // symmetry again FIXME why not vec[3]!!
      matrix[0][2] = vec[4];
      matrix[2][0] = 0.0;// symmetry again FIXME why not vec[4]!!
    }
  }

  Matrix<double> ErsatzMaterial::CalcHomogenizedTensor(Function* f)
  {
    const double cube_vol = grid->CalcHullVolume();
    unsigned int ex_size = me->GetNumberHomogenization(); // also ok when we do transform or robust

    assert((dim == 2 && ex_size == 3) || (dim == 3 && ex_size == 6));

    LOG_DBG(em) << "CHT f=" << f->ToString() << " ctxt=" << f->ctxt->GetExcitation()->robust_filter_idx << " f=" << f->GetExcitation()->robust_filter_idx;
    assert(f->ctxt->GetExcitation()->robust_filter_idx == f->GetExcitation()->robust_filter_idx);

    Matrix<double> test_strain_matrix_ij(dim, dim);
    Matrix<double> test_strain_matrix_kl(dim, dim);
    Matrix<double> result(ex_size, ex_size);
    result.Init();

    // we might have transformations
    Transform* trans = f->GetExcitation() != NULL ? f->GetExcitation()->transform : NULL;
    unsigned int meta = f->GetExcitation()->meta_index; // for transformations but also for robust
    for (unsigned int ij = 0;ij < ex_size;++ij)
    {
      // we need the transformation here to have the proper forward solution when we have multiple meta excitations
      // -> more than one rotation or robust
      SetTestStrainMatrix(test_strain_matrix_ij, f->ctxt->GetExcitation(ij, meta)->test_strain);
      StdVector<SingleVector*>& u1 = forward.Get(f->ctxt->GetExcitation(ij, f))->elem[App::MECH]; // equal to \chi^{ij}
      for (unsigned int kl = 0;kl < ex_size;++kl)
      {
        if (ij > kl) // already computed this entry!
        {
          // by construction, the resulting matrix must be symmetric
          // so we can skip the calculation and do a simple assignment instead
          result[ij][kl] = result[kl][ij];
          continue;
        }
        SetTestStrainMatrix(test_strain_matrix_kl, f->ctxt->GetExcitation(kl, meta)->test_strain);
        StdVector<SingleVector*>& u2 = forward.Get(f->ctxt->GetExcitation(kl, f))->elem[App::MECH]; // equal to \chi^{kl}
        // loop over elements. In the gradient case not summed up

        for (int e = 0, ne = design->GetNumberOfElements(); e < ne; ++e)
        {
          // When we rotate, the state u is based on a transformation of e, hence we need here to transform the element
          // BUT do NOT transform the state (they match already).
          DesignElement* de = design->ApplyTransformations(&design->data[e], &design->data[e], trans);
          LOG_DBG2(em) << "CHT: trans e=" << e << " -> " << de->GetIndex();

          Vector<double>& u1_vec = dynamic_cast<Vector<double>&>(*u1[e]);
          Vector<double>& u2_vec = dynamic_cast<Vector<double>&>(*u2[e]);
          // prepare for calculation
          LOG_DBG3(em) << "CHT f=" << f->ToString() << " ij=" << ij << " kl=" << kl << " e=" << e << " u1=" << u1_vec.ToString();
          LOG_DBG3(em) << "CHT f=" << f->ToString() << " ij=" << ij << " kl=" << kl << " e=" << e << " u2=" << u2_vec.ToString();

          // transformed de
          double p = CalcHomogenizedElementProduct(this, f->ctxt, de, false, u1_vec, u2_vec, test_strain_matrix_ij, test_strain_matrix_kl);

          assert(p < 1e100);

          LOG_DBG3(em) << "CHT f=" << f->ToString() << " ij=" << ij << " kl=" << kl << " e=" << e << " p=" << p;


          result[ij][kl] += p / cube_vol;// normalize for volume
        }
      } // end of kl loop
    } // end of ij loop

    LOG_DBG(em) << "CHT f=" << f->ToString() << " ex=" << f->GetExcitation()->GetFullLabel() << " mi=" << f->GetExcitation()->meta_index << " -> " << result.ToString();
    // save e.g. for CommitIteration()
    homogenizedTensor[f->GetExcitation()->meta_index].Assign(result, 1.0);
    return result;
  }

  void ErsatzMaterial::CalcHomogenizedTrackingGradient(const Matrix<double>& target, const Matrix<double>& hom, Function* f)
  {
    const double cube_vol = grid->CalcHullVolume();
    Context* ctxt = f->ctxt;

    Matrix<double> diff_tensor;
    diff_tensor = target - hom;
    const unsigned int ex_size(me->excitations.GetSize());
    assert((dim == 2 && ex_size == 3) || (dim == 3 && ex_size == 6));
    Matrix<double> test_strain_matrix_ij(dim, dim);
    Matrix<double> test_strain_matrix_kl(dim, dim);

    // our derivative tensor
    Matrix<double> hom_tensor_deriv(ex_size, ex_size);
    hom_tensor_deriv.Init();// we set and do not add - hence one init is enough

    // we might have transformations
    assert(f->GetExcitation() != NULL);

    unsigned int meta = f->GetExcitation()->meta_index;
    // loop over elements.
    for (int e = 0, ne = design->GetNumberOfElements();e < ne;++e)
    {
      DesignElement* de = &design->data[e];
      for (unsigned int ij = 0;ij < ex_size;++ij)
      {
        SetTestStrainMatrix(test_strain_matrix_ij, ctxt->GetExcitation(ij, meta)->test_strain);
        StdVector<SingleVector*>& u1 = forward.Get(ctxt->GetExcitation(ij, f))->elem[App::MECH]; // equal to \chi^{ij}
        Vector<double>& u1_vec = dynamic_cast<Vector<double>&>(*u1[e]);
        for (unsigned int kl = 0;kl < ex_size;++kl)
        {
          if (ij > kl) // already computed this entry!
          {
            // by construction, the resulting matrix must be symmetric
            // so we can skip the calculation and do a simple assignment instead
            hom_tensor_deriv[ij][kl] = hom_tensor_deriv[kl][ij];
            continue;
          }
          SetTestStrainMatrix(test_strain_matrix_kl, ctxt->GetExcitation(kl, meta)->test_strain);
          StdVector<SingleVector*>& u2 = forward.Get(ctxt->GetExcitation(kl, f))->elem[App::MECH]; // equal to \chi^{kl}
          Vector<double>& u2_vec = dynamic_cast<Vector<double>&>(*u2[e]);
          // prepare for calculation
          double p = CalcHomogenizedElementProduct(this, f->ctxt, de, true, u1_vec, u2_vec, test_strain_matrix_ij, test_strain_matrix_kl);
          hom_tensor_deriv[ij][kl] = p / cube_vol;// normalize for volume
        } // end of kl loop

      } // end of ij loop

      // hom_tensor_deriv is completely set.
      // (E^* - E^H) * - d(E^H)/d(rho_e) -> therefore the minus !
      double grad = -1.0 * diff_tensor.FrobeniusProduct(hom_tensor_deriv);
      de->AddGradient(f, grad);
    } // element loop

  }

  void ErsatzMaterial::CalcHomFrobeniusProductGradient(const Matrix<double>& par, const Matrix<double>& hom, Function* f)
  {
    // J  = sum_ij E_ij*D_ij
    // dJ = sum_ij dE_ij*D_ij
    // CalcHomogenizedTensorEntry((i, j), derivative = true, tmp_grad_out) sets the dE_ij in tmp_grad_out
    StdVector<double> tmp_grad_out;
    unsigned int meta = f->GetExcitation()->meta_index;
    for (unsigned int y = 0;y < par.GetNumRows();y++)
    {
      for (unsigned int x = 0;x < par.GetNumCols();x++)
      {
        boost::tuple<int,int,double> entry = boost::make_tuple(x + 1, y + 1, 0.0);
        tmp_grad_out.Init(0.0);
        CalcHomogenizedTensorEntry(f->ctxt, entry, true, tmp_grad_out, meta);
        double d_ij = par[y][x];
        for (int e = 0, ne = design->GetNumberOfElements();e < ne;++e)
        {
          DesignElement* de = &design->data[e];
          de->AddGradient(f, tmp_grad_out[e] * d_ij);
        }
      }

    }

  }

  double ErsatzMaterial::CalcHomogenizedTensorConstraint(Condition* g, bool derivative)
  {
    // me make use of the multi-purpose CalcHomogenizedTensorEntry()
    StdVector<double> grad;
    double result = 0.0;

    Transform* trans = g->GetExcitation()->transform;
    unsigned int meta = g->GetExcitation()->meta_index;

    // we have a list of int,int,double tuples which are added with the double factor.
    // E11 = <0,0,x>
    for(unsigned int i = 0; i < g->coords.GetSize(); i++)
    {
      boost::tuple<int, int, double>& entry = g->coords[i];
      double t = CalcHomogenizedTensorEntry(g->ctxt, entry, derivative, grad, meta);
      double factor = boost::get<2>(entry);

      if(derivative)
      {
        for(int o = 0, ne = design->GetNumberOfElements(); o < ne; ++o)
        {
          DesignElement* de = design->ApplyTransformations(&design->data[o], &design->data[o], trans);
          de->AddGradient(NULL, g, factor * grad[de->GetIndex()]);
          LOG_DBG2(em) << "CHTC: g=" << g->ToString() << " deriv " << i << " t=" << t << " fac=" << factor << " o=" << o << " de=" << de->ToString()
                       << " dg=" << grad[de->GetIndex()];
        }
      }
      else
      {
        result += factor * t;

        Matrix<double>& ht = homogenizedTensor[g->GetExcitation()->meta_index];
        LOG_DBG(em) << "CHTC: g=" << g->ToString() <<  " ht_idx=" << g->GetExcitation()->meta_index;

        ht[boost::get<0>(entry)-1][boost::get<1>(entry)-1] = t;
        // all tensors are symmetric. Makes reading easier!
        ht[boost::get<1>(entry)-1][boost::get<0>(entry)-1] = t;

        LOG_DBG(em) << "CHTC: g=" << g->ToString() << " coord=" << i << " [" << boost::get<0>(entry)-1
                     << "][" << boost::get<1>(entry)-1 << "] = " << t;
      }
    }
    return result;
  }

  double ErsatzMaterial::CalcHomogenizedTensorEntry(Context* ctxt, const boost::tuple<int,int,double> entry, bool derivative, StdVector<double>& grad_out, unsigned int meta)
  {
    const double cube_vol = grid->CalcHullVolume();

    assert((dim == 2 && ctxt->excitations.GetSize() >= 3) || (dim == 3 && ctxt->excitations.GetSize() >= 6)); // for meta exctiations it is more
    Matrix<double> test_strain_matrix_ij(dim, dim);
    Matrix<double> test_strain_matrix_kl(dim, dim);
    const unsigned int ij = boost::get<0>(entry) - 1;
    const unsigned int kl = boost::get<1>(entry) - 1;

    // the test strain itself shall be independent of the meta excitation
    assert(ctxt->excitations[ij]->test_strain == ctxt->GetExcitation(ij, meta)->test_strain);
    SetTestStrainMatrix(test_strain_matrix_ij, ctxt->excitations[ij]->test_strain);

    // for multiple meta excitations (rotations, robustness) take the corresponding state
    StdVector<SingleVector*>& u1 = forward.Get(ctxt->GetExcitation(ij, meta))->elem[App::MECH]; // equal to \chi^{ij}

    SetTestStrainMatrix(test_strain_matrix_kl, ctxt->excitations[kl]->test_strain);
    StdVector<SingleVector*>& u2 = forward.Get(ctxt->GetExcitation(kl, meta))->elem[App::MECH];// equal to \chi^{kl}

    double result = 0.0;

    if (derivative)
      grad_out.Resize(design->GetNumberOfElements());

    Transform* trans = ctxt->GetExcitation(0, meta)->transform; // the base 0 is absolutely ok as this is the fast excitation index with meta the slow index

    // loop over elements. In the gradient case not summed up
    for (int e = 0, ne = design->GetNumberOfElements();e < ne;++e)
    {
      Vector<double>& u1_vec = dynamic_cast<Vector<double>&>(*u1[e]);
      Vector<double>& u2_vec = dynamic_cast<Vector<double>&>(*u2[e]);

      // for transformation we transform the element design for CHEP and the index for storing the gradient
      // but we do NOT transform the state as this has already been done for the state problem
      DesignElement* de = design->ApplyTransformations(&design->data[e], &design->data[e], trans);

      // prepare for calculation
      double p = CalcHomogenizedElementProduct(this, ctxt, de, derivative, u1_vec, u2_vec, test_strain_matrix_ij, test_strain_matrix_kl);
      result += p / cube_vol;// normalize for volume

      LOG_DBG2(em) << "CHTE ij=" << ij << " kl=" << kl << " der=" << derivative << " meta=" << meta << " e=" << e << "de=" << de->ToString() << " p=" << p << " re=" << result;
      if (derivative)
      {
        grad_out[de->GetIndex()] = result;
        result = 0.0; // reset such that we do not sum up for the next case!
      }
    }

    LOG_DBG(em) << "CHTE ij=" << ij << " kl=" << kl << " der=" << derivative << " meta=" << meta << " re=" << result;
    return result; // in the non-derivative case this is the sum.
  }


  double ErsatzMaterial::CalcHomogenizedElementProduct(ErsatzMaterial* obj, Context* ctxt, DesignElement* de, bool derivative, Vector<double>& u1_vec, Vector<double>& u2_vec, Matrix<double>& test_strain_matrix_ij, Matrix<double>& test_strain_matrix_kl)
  {
    assert(u1_vec.NormL2() > 0);
    assert(u2_vec.NormL2() > 0);
    assert(test_strain_matrix_ij.NormL2() > 0);
    assert(test_strain_matrix_kl.NormL2() > 0);
    assert(u1_vec.GetSize() == u2_vec.GetSize());

    const unsigned int dim_ = obj->dim;

    LOG_DBG3(em) << "CHEP: de=" << de->ToString() << " der=" << derivative; // << " u1_vec=" << u1_vec.ToString() << " u2_vec=" << u2_vec.ToString();
    // TODO too much temporary matrices and vectors!
    // from the coordinates of this element we build a "test displacement" vector
    // u1(,2)_0. it contains linear strains which are assembled in the following lines
    // these strains are not unique! an arbitrary constant can be added without change
    // coordinates of "this" element
    Matrix<double> tmp_mat;
    // coordinates of current element
    domain->GetGrid()->GetElemNodesCoord(tmp_mat, de->elem->connect, true);
    Matrix<double> u1_tmp;
    u1_tmp = test_strain_matrix_ij * tmp_mat;
    Matrix<double> u2_tmp;
    u2_tmp = test_strain_matrix_kl * tmp_mat;
    assert(u1_tmp.GetNumCols() == u2_tmp.GetNumCols());
    assert(u1_tmp.GetNumRows() == u2_tmp.GetNumRows());
    assert(u1_tmp.GetNumRows() == dim_);
    assert(u1_tmp.GetNumCols() * dim_ == u2_vec.GetSize());
    assert(u1_tmp.GetNumRows() * u1_tmp.GetNumCols() == u1_vec.GetSize());
    // u1_tmp, u2_tmp have to be transformed into vectors
    // 2D: from 2x4 to 8x1 on quad elems, 2x3 to 6x1 on triangles
    Vector<double> u1_0(u1_vec.GetSize());
    Vector<double> u2_0(u2_vec.GetSize());// u1 and u2 have the same size

    for (unsigned int out = 0, cols = u1_tmp.GetNumCols();out < cols;++out)
    {
      for (unsigned int in = 0;in < dim_;++in)
      {
        u1_0[out * dim_ + in] = u1_tmp[in][out]; // equal to \chi^{0(ij)}
        u2_0[out * dim_ + in] = u2_tmp[in][out]; // equal to \chi^{0(kl)}
      }
    }
    u1_0 -= u1_vec;
    u2_0 -= u2_vec;

    // any transformation is done outside by getting the tranformed e. The state is transformed from the state problem

    // reuse tmp_mat as elementK-Matrix
    // Matrix<double> k_mat;
    TransferFunction* tf = obj->design->GetTransferFunction(DesignElement::DENSITY, App::MECH);
    obj->SetElementK(ctxt, de, tf, App::MECH, &tmp_mat, derivative);

    assert(tmp_mat.GetNumRows() == tmp_mat.GetNumCols() && tmp_mat.GetNumCols() == u1_0.GetSize());

    // Vector<double> mat_vec = tmp_mat * u1_0;
    Vector<double> mat_vec(u1_0.GetSize());
    tmp_mat.Mult(u1_0, mat_vec);

    // LOG_DBG3(em) << "CHEP de=" << de->ToString() << " tmp_mat=" << tmp_mat.ToString();
    // LOG_DBG3(em) << "CHEP de=" << de->ToString() << " mat_vec=" << mat_vec.ToString();
    // LOG_DBG3(em) << "CHEP de=" << de->ToString() << " u2_0=" << u2_0.ToString();

    assert(mat_vec.GetSize() == u2_0.GetSize());
    // double result = mat_vec * u2_0;

    double result = mat_vec.Inner(u2_0);
    LOG_DBG3(em) << "CHEP de=" << de->ToString() << " result=" << result;
    assert(result < 1e100);
    return result;
  }

  double ErsatzMaterial::CalcFilteringGap(Condition* g, bool derivative) {
    /* Calculates squared difference between filtered and non-filtered tensor E*/
    //TODO: asserts
    //assert(g->GetDesignType() != )

    double result = 0, grad = 0;
    unsigned int n_elem = design->GetNumberOfElements();
    unsigned int dtype = design->FindDesign(g->GetDesignType());
    for(unsigned int i = 0; i < n_elem; i++)
    {
      DesignElement* de = dynamic_cast<DesignElement*>(design->GetDesignElement(dtype*n_elem+i));
      if (!derivative) {
        // (E_(ij) - filtered(E_ij))^2
        double error = de->GetDesign(DesignElement::PLAIN)- de->GetDesign(DesignElement::SMART) ;
        result += error * error;
      } else {
        // calculate derivative
        // We filter over this element and the neighbors.
        assert(de->simp != NULL);
        unsigned int fix = de->simp->DetermineFilterIndexNonInlined();
        const Filter& f = de->simp->filter[fix];

        assert(f.GetType() == Filter::DENSITY);
        //assert(de == DesignElement::COST_GRADIENT || de == DesignElement::CONSTRAINT_GRADIENT);
        //assert((g == NULL || (g->IsObjective() && de == DesignElement::COST_GRADIENT)) || (g == NULL || (!g->IsObjective() && de == DesignElement::CONSTRAINT_GRADIENT)) || (g == NULL || (g->IsObjective())));
        // projection has density filtering only in the fake filter problem but not in the original problem (which should not be density filtered anyway)
        //assert(g == NULL || g->ForDensityFiltering());

        // Density filtering for gradient is (Sigmund; Morphology-based black and white filters for topology optimization; 2007; eqn (35). (36)
        // p is rho and P is rho filtered! d f/d p_e = sum_i(in N_e) d f/d P_i * d P_i/d p_e with d P_i/d p_e = w(x_e)/ sum_j(in N_i) w(x_j)
        // note, that the stored value is already v = d f/d P_i
        grad = 0.0;
        if(f.density_ == Filter::STANDARD)
        {
          for(int j = -1, nj = (int) f.neighborhood.GetSize(); j < nj; j++)
          {
            const Filter::NeighbourElement* ne = j == -1 ? NULL : &f.neighborhood[j];
            const DesignElement* de_iter = j == -1 ? de : ne->neighbour;
            double v = 2 * (de_iter->GetDesign(DesignElement::PLAIN)- de_iter->GetDesign(DesignElement::SMART));
            double w = j == -1 ? f.weight : ne->weight;
            double var = j == -1 ? 1.:0.;

            if (de_iter->simp->filter[fix].weight_sum < 0.0)
                de_iter->simp->filter[fix].weight_sum = de_iter->simp->filter[fix].CalcWeightSum(true);

            double summand = v * (var -(w / de_iter->simp->filter[fix].weight_sum));
            grad += summand;

            // LOG_DBG3(desel) << "GDFG: el=" << de_->elem->elemNum << ": curr=" << de->elem->elemNum
            //                << " v= " << v  << " h=" << h << " w=" << w << " x_n=" << x_n << " w_sum=" << w_sum
            //                << " summand=" << summand << " sum=" << sum;
          }
        }
        design->data[dtype*n_elem+i].AddGradient(NULL, g, grad);
      }
      LOG_DBG2(em) << "GDFG: el=" << de->elem->elemNum << " de = "<< de->ToString() << " filtering_gap = "<< result <<" derivative ="<<derivative << " grad = "<< grad;
    }
    return result;
  }


  double ErsatzMaterial::CalcGreyness(Condition* g, bool derivative)
  {
    double greyness = 0.0; // element greyness
    int counter = 0;// to make it sure for different design variables!
    double lb, ub, span, org_value, value, grad, eval;
    lb = ub = value = grad = eval = span = 0.0;
    // we have to divide the gradients by their relative volume = fraction
    double fraction = g->GetDesignType() == DesignElement::DEFAULT ? design->data.GetSize() : design->GetNumberOfElements();
    // do we want the physical value?
    TransferFunction* tf = g->IsPhysical() ? design->GetTransferFunction(g->GetDesignType(), App::MECH) : NULL;
    // go over the complete design space to set gradients of other types to 0
    for(unsigned int i = 0; i < g->elements.GetSize(); i++)
    {
      DesignElement* de = g->elements[i];
      bool relevant = g->GetDesignType() == DesignElement::DEFAULT || g->GetDesignType() == de->GetType();
      if(relevant)
      {
        if(g->IsPhysical())
        {
          lb = tf->Transform(de->GetLowerBound());
          ub = tf->Transform(de->GetUpperBound());
          org_value = derivative ? tf->Derivative(de, DesignElement::SMART) : tf->Transform(de, DesignElement::SMART);
        }
        else
        {
          lb = de->GetLowerBound();
          ub = de->GetUpperBound();
          org_value = de->GetDesign(DesignElement::PLAIN);
        }

        span = ub-lb;

        if(span < 0.01) EXCEPTION("cannot calculate grayness with almost equal design bounds");

        // We normalize for a design variable from [0;1]
        // this has minor effect on density [0.001;1] but is important
        // for polarization[-1;1]
        value = (org_value - lb) / span;
      }

      if(derivative)
      {
        if(relevant)
        {
          // standard greyness without parameters! times 4 so we have 1 for maximum greyness
          // Note that we transformed to [0,1]
          // f(x)=4*(1-x)x * 4 = 4*(-x^2+x+1)
          // f'(x)= 4*(-2x+1)
          grad = (-8.0 * value + 4.0)/span;
        }
        // divide by fraction
        // the gradient of non-relevant is 0
        grad = relevant ? grad / fraction : 0.0;

        // set 0.0 if not relevant
        design->data[i].AddGradient(NULL, g, grad);
      }
      else // not derivative but function evaluation
      {
        // we normalized the greyness to value from [0;1]
        if(relevant)
        {
          eval = 4.0* (1-value) * (value);
          greyness += eval;
          counter++;
        }
      }
      LOG_DBG3(conditions) << Condition::type.ToString(g->GetType())
      << " derive=" << derivative << " relevant=" << relevant
      << " elem " << de->elem->elemNum << " des_value: " << de->GetDesign(DesignElement::PLAIN)
      << " value = " << org_value
      << " -> " << value << " grad=" << grad << " eval=" << eval
      << " fraction=" << fraction << " counter=" << counter;

    }
    return greyness / (double)(counter);
  }

  double ErsatzMaterial::CalcLocalConstraint(Condition* g, bool derivative)
  {
    LocalCondition* loc_cond = dynamic_cast<LocalCondition*>(g);
// take care, similar logic in SlopeCondition::GetSparsityPattern() !
    assert(loc_cond->IsLocal());
// The neighborhood is already determined
    Function::Local* local = g->GetLocal();
    assert(local != NULL);
    Function::Local::Identifier& id = loc_cond->GetCurrentVirtualContext();
    double res = -1.0;
    if (derivative)
    {
      id.EvalGradient(local);
    }
    else
    {
      res = id.EvalFunction(local);
    }
    return res;
  }

  double ErsatzMaterial::CalcGlobalFunction(Function* f, bool derivative, const Vector<double>* von_mises_stress)
  {
    LOG_DBG(em) << "CGF c=" << f->type.ToString(f->GetType()) << " derivative=" << derivative;
    Function::Local* local = f->GetLocal();
    assert(local != NULL);
    // the neighborhoods are already defined by Local!
    StdVector<Function::Local::Identifier>& vem = local->virtual_elem_map;
    if(!derivative)
    {
      // evaluate the function values, which is
      // max(0, x_i - x_i+1 - c) and max(0,x_i+1 - x_i - c)
      double res = 0.0;

      assert(von_mises_stress == NULL || (von_mises_stress->GetSize() == vem.GetSize()));
      for(unsigned int i = 0; i < vem.GetSize(); i++)
      {
        Function::Local::Identifier& id = vem[i];
        double fv = id.EvalFunction(local, false, von_mises_stress != NULL ? (*von_mises_stress)[i] : -1.0);
        res += fv;
        LOG_DBG2(em) << "CGF: !d c=" << f->type.ToString(f->GetType()) << " i=" << i << " de="
                     << ( typeid(id.element) == typeid(DesignElement*) ? (int)dynamic_cast<DesignElement*>(id.element)->elem->elemNum : -1 ) << " sign=" << id.sign
                     << " fv=" << fv << " -> " << res;
      }

      return res;
    }
    else
    {
      // the gradient g/x_i = 0 or 2 * (x_i+1 - x_i - c) * -1 or 2 * (x_i - x_i+1 - c) * 1
      // in the non-periodic case is the number of functions per design variable not constant,
      // e.g. the most upper right design has no slope constraint.
      for(int j = 0; j < (int) vem.GetSize(); j++)
      {
        Function::Local::Identifier& id = vem[j];
        id.EvalGradient(local);
      }

      return 0.0; // gradient case has no information
    }
  }

  void ErsatzMaterial::SolveStateProblem(Excitation* ev_only_excite)
  {
    // if ev_only_excite is set we use the given excitation
    // -> it shall not coincide
    assert(!(ev_only_excite != NULL && me->IsEnabled()));
    // shall we normalize afterwards?
    bool normalize = false;

    // we have to check objectives and active (non local) constraints

    // We traverse all excitations and conditionally perform a context switch. Because of the context switch
    // we need to solve the adjoints within the same context

    StdVector<Function*> funcs = GetFunctions(false);
    // when ev_only_excite is given we evaluate for one excitation only
    for(unsigned int e = 0; e < (ev_only_excite != NULL ? 1 : me->excitations.GetSize()); e++)
    {
      Excitation& excite = ev_only_excite != NULL ? *ev_only_excite : me->excitations[e];
      // ! sets the context!!
      bool switched = excite.Apply(true); // make the multiple sequence switch if necessary
      assert(excite.sequence == context->sequence);

      if(context->DoBloch() && (e == 0 || switched)) // handle no multiple sequence case and multiple sequence case with bloch not first
        context->GetEigenFrequencyDriver()->SetupBlochPlot(); // the plot is written for each iteration and contains all modes for all wave numbers

      if(context->DoLBM()) {
        // in autoscale case we are still in the BaseOptimizer constructor
        boost::shared_ptr<Timer> eval_timer = baseOptimizer_ != NULL ? baseOptimizer_->GetRunnungEvalTimer() : boost::shared_ptr<Timer>();
        if(eval_timer)
          eval_timer->Stop();

        LatticeBoltzmannPDE* lbmPde = context->GetLatticeBoltzmannPDE();
        assert(lbmPde != NULL);
        lbmPde->Solve();

        if(eval_timer)
          eval_timer->Start();
      }
      else
        Optimization::SolveStateProblem(&excite); // this is true for all problem types

      if(!context->DoLBM() && !IsTransient()) // transient solutions are read per timestep
      {
        // in the eigenvalue case we store the modes separately, similar to timesteps
        if(!context->IsEigenvalue())
          StorePDESolution(forward, excite, NULL, -1, true, true, false, NO_DERIVTYPE, "forward"); // no solution and mode is -1 as it is not set
        else
          for(int m = 0; m < (int) context->GetEigenFrequencyDriver()->eigenFreqs->GetSize() ; m++)
            StorePDESolution(forward, excite, NULL, m, true, true, true, NO_DERIVTYPE, "forward"); // only in the ev case we need to save the solution
      }

      // check for each excitation all functions if we shall solve the adjoint - take care about the context!
      for(unsigned fi = 0; fi < funcs.GetSize(); fi++)
      {
        Function* f = funcs[fi];
        assert(f != NULL);
        // some functions need the selection vector for function evaluation, e.g. output
        if(f->NeedsSelectionVector())
          ConstructSelection(excite, f, false);// don't change the rhs of the system but restore

        // in the harmonic case the system matrix depends on the frequency. Hence we have to
        // use the current assembly and factorization to solve the adjoint problem.
        if(f->IsAdjointBased() && DoSolveAdjointWithState() && f->ctxt == context) // the context is set properly
          SolveAdjointProblem(&excite, f); // not called in a standard case

        // when we do multiple excitations with adjusted weights we calculate the objective here
        // to find the best weights. CalcObjective is so cheap, it is done later again by
        // the optimizer but then the weights are set
        if(me->DoAdjustWeights())
        {
          // in the first iteration we adjust the weights
          // stride = 0 is only the first time
          if((me->stride < 1 && GetCurrentIteration() == 0) || ((GetCurrentIteration() % me->stride) == 0))
          {
            excite.cost = CalcFunction(excite, f, false); // to be normalized
            normalize = true;
          }
        }
      }
    }
    // we need only to normalize when the properties have changed. This also reflects the stride
    if(normalize)
      me->NormalizeMultipleExcitations(&objectives);
  }
  void ErsatzMaterial::SolveAdjointProblems(Excitation* ev_only_excite)
  {
    // solve all adjoints needed for gradient calculation
    assert(!(ev_only_excite != NULL && me->IsEnabled()));

    // check ErsatzMaterial::SolveStateProblem() and DoSolveAdjointWithState().
    if(!DoSolveAdjointWithState()) // was it already computed in ErsatzMaterial::SolveStateProblem() ?
    {
      for(unsigned int e = 0; e < me->excitations.GetSize(); ++e)
      {
        Excitation* excite = ev_only_excite != NULL ? ev_only_excite : &me->excitations[e];

        // it seems an excite.Apply() is missing for robustness ?!

        // will check for all functions and call SolveAdjointProblem() if the function requires it
        Optimization::SolveAdjointProblems(excite);
      }
    }
  }

  void ErsatzMaterial::StorePDESolution(StateContainer& solutions, Excitation& excite, Function* f, int timestep_mode, bool read_sol, bool read_rhs, bool save_sol, TimeDeriv derivative, const std::string& comment)
  {
    assert(context->pde != NULL);
    assert(context->sequence == excite.sequence);
    assert(f == NULL || f != NULL); // f is NULL for forward problem

    Assemble* assemble = context->pde->GetAssemble(); // context is valid as it was switched by Excitation::Apply(true)
    StateSolution& sol = *(solutions.Get(excite, f, timestep_mode, derivative));
    SingleVector* raw = NULL;

    // in the eigenvalue case we have not only one solution vector but one for each mode.
    // Stores the mode in the the solution part of algsys
    if(context->IsEigenvalue())
      assemble->GetAlgSys()->GetEigenMode(timestep_mode);

    // store solution element wise for gradient and raw vector for objective.
    // This is redundant as currently the solution is the global one!
    for(map<App::Type, SinglePDE*>::iterator it = context->pdes.begin(); it != context->pdes.end(); ++it)
    {
      stringstream ss;
      ss << "SPDES: prob=" << comment << " excite=" << excite.index << " pde: " << it->first << " timestep_mode=" << timestep_mode;

      if(read_sol)
      {
        sol.Read(StateSolution::ELEMENT_VECTORS, it->second, it->first, save_sol, derivative);
        raw = sol.Read(StateSolution::RAW_VECTOR, it->second, it->first, save_sol, derivative);

        LOG_DBG2(em) << ss.str() << " sol: " << raw->ToString();
      }

      if(read_rhs)
      {
        sol.Read(StateSolution::RHS_VECTOR, it->second, it->first, save_sol, derivative);
        LOG_DBG2(em) << ss.str() << " rhs: " << sol.GetVector(StateSolution::RHS_VECTOR)->ToString();
      }

      if(context->IsEigenvalue())
      {
        assert(timestep_mode >= 0);
        SingleVector* ef = context->GetEigenFrequencyDriver()->eigenFreqs;
        sol.eigenfreq = ef->GetEntryType() == BaseMatrix::DOUBLE ? ef->GetDoubleEntry(timestep_mode) : ef->GetComplexEntry(timestep_mode).real();
        LOG_DBG(em) << ss.str() << " store freq " << f;
      }
    }
  }
/* tracking and transient
  void ErsatzMaterial::TimeStepCalculated(UInt timeStep, AdjointParameters* adjParams)
  {
    // is only called in transient case
    assert(IsTransient());

    // drivers start counting steps with 1
    if(adjParams == NULL)
    {
      StorePDESolution(forward, *context->GetExcitation(), NULL, timeStep-1, true, false, false, NO_DERIVTYPE, "forward");
      StorePDESolution(forward, *context->GetExcitation(), NULL, timeStep-1, true, false, false, FIRST_DERIV, "forward-derivative");
      StorePDESolution(forward, *context->GetExcitation(), NULL, timeStep-1, true, false, false, SECOND_DERIV, "forward-second-derivative");
    }
    else
    {
      StorePDESolution(adjoint, *context->GetExcitation(), adjParams->GetFunction(), timeStep-1, true, false, false, NO_DERIVTYPE, "adjoint");
    }
  }
*/
  /* tracking and transient
  void ErsatzMaterial::RhsCalculated(AdjointParameters* adjParams)
  {
    if(IsTransient())
    { // only in transient case this is needed
      if(adjParams == NULL)
      {
        StorePDESolution(forward, *context->GetExcitation(), NULL, context->GetDriver()->GetActStep("mech")-1, false, true, false, NO_DERIVTYPE, "forward");
      }
      else
      {
        StorePDESolution(adjoint, *context->GetExcitation(), adjParams->GetFunction(), context->GetDriver()->GetActStep("mech")-1, false, true, false, NO_DERIVTYPE, "adjoint");
      }
    }
  }
  */
/*
  StdVector<pair<SinglePDE*,IdBcList> > ErsatzMaterial::SetHDBC()
  {
    // IDBC values in the forward problem are homogeneous in the adjoint PDE. Consider elec and mech
    // store org value as idbc_list per pde
    StdVector<pair<SinglePDE*,IdBcList> > org_idbc;
    assert(false);


    org_idbc.Reserve(pdes.size());
    for(map<App::Type, SinglePDE*>::iterator it = pdes.begin(); it != pdes.end(); ++it)
    {
      // we need a deep copy. Create new list for the PDE
      org_idbc.Push_back(make_pair(it->second, IdBcList()));

      // get the idbc list of the current pde
      shared_ptr<BaseFeFunction> fe = it->second->GetFeFunction(it->second->GetNativeSolutionType());

      IdBcList idbc_list = fe->GetInHomDirichletBCs();
      for(unsigned int bc = 0; bc < idbc_list.GetSize(); bc++)
      {
        // store original value -> we have to do a deep copy!
        org_idbc.Last().second.Push_back(shared_ptr<InhomDirichletBc>(new InhomDirichletBc(*(idbc_list[bc]))));
        // make homogeneous values! -> this is only stored in the PDE, we have to call SetBCs() later!
        assert(false);
        (*idbc_list[bc]).value->Get = "0.0";
        (*idbc_list[bc]).phase = "0.0";
        LOG_DBG(em) << "Set IDBC to HDBC (" << it->second->GetName() << ") -> value "
                     << (*(org_idbc.Last().second[bc])).value << " -> " << (*(it->second->GetIDBCList()[bc])).value;
      }

      // apply the new boundary condition
      it->second->SetBCs();
    }
    return org_idbc;
  }

  void ErsatzMaterial::ResetHDBC(StdVector<pair<SinglePDE*,IdBcList> >& org_idbc)
  {
    // reset the original IDBC which were HDBC for the adjoint PDE
    assert(false);
    for(unsigned int p = 0; p < org_idbc.GetSize(); p++)
    {
      // get the idbc list of the current pde
      IdBcList idbc_list = org_idbc[p].first->GetIDBCList();
      IdBcList org_list = org_idbc[p].second;
      assert(idbc_list.GetSize() == org_list.GetSize());
      for(unsigned int bc = 0; bc < idbc_list.GetSize(); bc++)
      {
        (*idbc_list[bc]).value = (*org_list[bc]).value;
        (*idbc_list[bc]).phase = (*org_list[bc]).phase;
        LOG_DBG(em) << "Reset HDBC (" << org_idbc[p].first->GetName() << ") -> " << (*(org_idbc[p].first->GetIDBCList()[bc])).value;
      }
      org_idbc[p].first->SetBCs();
    }

  }
*/

  void ErsatzMaterial::SolveAdjointProblem(Excitation* excite, Function* f)
  {
    if(context->IsComplex())
      SolveAdjointProblem<std::complex<double> >(excite, f);
    else
      SolveAdjointProblem<double>(excite, f);
  }

  template<class T>
  void ErsatzMaterial::SolveAdjointProblem(Excitation* excite, Function* f)
  {
    assert(baseOptimizer_ != NULL || !baseOptimizer_->GetOptimierTimer()->IsRunning()); // https://cfs.mdmt.tuwien.ac.at/trac/ticket/263#ticket
    boost::shared_ptr<Timer> eval_timer = baseOptimizer_ != NULL ? baseOptimizer_->GetRunnungEvalTimer() : boost::shared_ptr<Timer>();
    if(eval_timer)
      eval_timer->Stop();

    assert(excite != NULL);
    excite->Apply(false); // the context shall be already switched
    assert(excite->sequence == context->sequence);
    assert(context->pde != NULL);
    assert(f->ctxt == context);

    switch(f->GetType())
    {
      // these objectives need their adjoint problems only for gradient calculation
      case Function::COMPLIANCE:
      if(IsTransient())
      { // in transient case, everything has an adjoint
        // Optimization::SolveAdjointProblem(excite, f);
      }
      break;
      case Function::TRACKING:
      // these objectives need their adjoint problems only for gradient calculation
      // Optimization::SolveAdjointProblem(excite, f);

      if(!IsTransient())
      { // transient solutions are read every timestep
        StorePDESolution(adjoint, *excite, f, 0, true, false, false, NO_DERIVTYPE, "adjoint");
      }

      // write back the solution s.th. CommitIteration() makes StoreResults() properly.
      forward.Get(excite)->Write(context->pde);
      break;

      case Function::OUTPUT:
      case Function::SQUARED_OUTPUT:
      case Function::CONJUGATE_COMPLIANCE:
      case Function::ABS_OUTPUT:
      case Function::GLOBAL_DYNAMIC_COMPLIANCE:
      case Function::DYNAMIC_OUTPUT:
      case Function::ELEC_ENERGY:
      case Function::ENERGY_FLUX:
      case Function::STRESS:
      case Function::TEMP_TRACKING_AT_INTERFACE: // track boundary driven load
      case Function::STRESS_DENSITY:
      {
        // these objectives need their adjoint problems for the calculation of the objective value
        // they are directly solved after the StateProblem
        // they are no more solved for gradient calculation (this only works if the optimizer always evaluates the function before the gradient)
        // when doing linesearch, this slows down optimization if solution is only needed for gradient

        // the forward problem was already solved and stored !!

        // Set the rhs and solve for it
        SetAndSolveAdjointRHS<T>(*excite, f);

        // store the stuff -> no rhs but special handling of element results
        StorePDESolution(adjoint, *excite, f, -1, true, false, true, NO_DERIVTYPE, "adjoint");

        // write back the solution s.th. CommitIteraion() makes StoreResults() properly.
        forward.Get(excite)->Write(context->pde);
        break;
      }

      default:
      assert(false);
    }
    if(eval_timer)
      eval_timer->Start();
  }

  template<class T>
  void ErsatzMaterial::SetAndSolveAdjointRHS(Excitation& excite, Function* f)
  {
    assert(context->sequence == excite.sequence);
    assert(f->ctxt == context);
    Assemble* assemble = context->pde->GetAssemble();

    // the adjoint RHS might be an output stuff, then the loads are changed.
    // save and restore them in any case.
    StdVector<LinearFormContext*> org_forms = assemble->GetLinForms();
    // set pseudo loads (if there are output nodes)
    if (f->NeedsSelectionVector()) // TODO: rhs? no, since selection vector is assembled automatically
      ConstructSelection(excite, f, true);// is actually already set for the forward calculation - who cares?

    // any adjoint PDE has HDBC instead of IDBC. We Store the IDBC, add the BC as HDBC, solve, reset the IDBC and remove the additional HDBC
    shared_ptr<BaseFeFunction> fe = context->pde->GetFeFunction(context->pde->GetNativeSolutionType()); // no reference but copy constructor
    IdBcList  org_idbc = fe->GetInHomDirichletBCs();
    fe->GetHomDirichletBCs().Reserve(org_idbc.GetSize()); // what will be added temporarily
    for(unsigned int i = 0; i < org_idbc.GetSize(); i++)
    {
      shared_ptr<HomDirichletBc> hdbc(new HomDirichletBc);
      hdbc->dofs = org_idbc[i]->dofs;
      hdbc->entities = org_idbc[i]->entities;
      hdbc->result = org_idbc[i]->result;
      fe->GetHomDirichletBCs().Push_back(hdbc);
    }
    fe->GetInHomDirichletBCs().Resize(0);
    fe->ApplyBC();

    // set the adjoint rhs
    ConstructAdjointRHS(excite, f);

    assert(context->GetDriver()->GetAnalysisId().adjoint == false);
    context->GetDriver()->GetAnalysisId().adjoint = true;

    // calculate adjoint problem
    assemble->GetAlgSys()->Solve();

    context->GetDriver()->GetAnalysisId().adjoint = false;

    // reset the boundary conditions
    fe->GetInHomDirichletBCs() = org_idbc; // I love copy constructors
    fe->GetHomDirichletBCs().Resize(fe->GetHomDirichletBCs().GetSize() - org_idbc.GetSize()); // remove "artificial" hdbc
    fe->ApplyBC();

    // reset the original loads, they have been changed in the output case
    assemble->GetLinForms() = org_forms;
  }

  void ErsatzMaterial::ConstructSelection(Excitation& excite, Function* f, bool alter_rhs)
  {
    assert(context->sequence == excite.sequence);
    assert(f->ctxt == context);
    Assemble* assemble = context->pde->GetAssemble();
    // in SolveStateProblem() the clean variant for the objective
    StdVector<LinearFormContext*> org_forms;
    if (!alter_rhs)
      org_forms = assemble->GetLinForms();

    if (f->GetType() != Objective::CONJUGATE_COMPLIANCE){
      // overwrite the assemble loads with "pseudo loads"s loads
      assemble->GetLinForms() = f->output_forms;
    }
    // set our own RHS but delete first as Assemble adds
    assemble->GetAlgSys()->InitRHS();
    // assemble the output nodes
    assemble->AssembleLinRHS();

    // save the "pseudo loading" which is the selection as the rhs for the adjoint
    // This is exactly what has been constructed. Not that for an adjoint RHS it needs
    // post processing.
    adjoint.Get(excite, f)->Read(StateSolution::SEL_VECTOR, context->pde);
    if (!alter_rhs)
      assemble->GetLinForms() = org_forms;

    LOG_DBG2(em) << "ConstructSelection: excite=" << excite.index << " f=" << f->ToString() << " alter=" << alter_rhs
        << " sel=" << adjoint.Get(excite, f)->GetVector(StateSolution::SEL_VECTOR)->ToString(1);

  }
  void ErsatzMaterial::ConstructRealAdjointRHS(Excitation& excite, Function* f)
  {
    assert(context->sequence == excite.sequence);
    assert(f->ctxt == context);
    Vector<double> rhs; // own OLAS vector
    switch(f->GetType())
    {
      case Function::OUTPUT:
      case Function::SQUARED_OUTPUT:
      {
        Vector<double>& l = adjoint.Get(excite, f)->GetRealVector(StateSolution::SEL_VECTOR);
        rhs.Resize(l.GetSize());
        rhs = l * -1.0;
        break;
      }
      case Function::STRESS:
      case Function::STRESS_DENSITY:
      {
        StressConstraint<double> sc(&excite, f, this, &forward);
        sc.CalcAdjointRHS(rhs);
        break;
      }
      case Function::TEMP_TRACKING_AT_INTERFACE:
      {
        CalcAdjointRHSStateTracking(excite, f, f->GetParameter(), rhs);
        break;
      }
      default:
      assert(false);
      break;
    }

    shared_ptr<BaseFeFunction> fe = context->pde->GetFeFunction(context->pde->GetNativeSolutionType());
    // we cannot easily set the rhs. Therefore we set it to 0 and add our own rhs
    fe->GetSystem()->InitRHS(fe->GetFctId());
    fe->GetSystem()->SetFncRHS(rhs, fe->GetFctId());

    LOG_DBG2(em) << "CARHS<double>: f=" << f->ToString() << " obj=" << f->IsObjective() << " rhs before solving: " << rhs.ToString(1) << " max_norm=" << rhs.NormMax();
    assert(rhs.NormMax() != 0.0);
  }

  void ErsatzMaterial::ConstructComplexAdjointRHS(Excitation& excite, Function* f)
  {
    // we handle only complex cases
    assert(context->sequence == excite.sequence);
    assert(f->ctxt == context);

    Vector<complex<double> >& u = forward.Get(excite)->GetComplexVector(StateSolution::RAW_VECTOR);
    Vector<complex<double> >& l = adjoint.Get(excite, f)->GetComplexVector(StateSolution::SEL_VECTOR);
    LOG_DBG2(em) << "AdjustComplexAdjointRHS: u = " << u.ToString();
    LOG_DBG2(em) << "AdjustComplexAdjointRHS: l = " << l.ToString();

    // create a OLAS vector
    Vector<complex<double> >& rhs = adjoint.Get(excite, f)->GetComplexVector(StateSolution::RHS_VECTOR);
    rhs.Resize(u.GetSize());
    switch(f->GetType())
    {
      case Function::DYNAMIC_OUTPUT: // rhs is from "output loads" and set in adjoint...rhs
      // the correct conjugate_output case is -L * u*, always complex!
      for(unsigned int i = 0; i < rhs.GetSize(); i++)
      rhs[i] = -1.0 * l[i] * std::conj(u[i]);
      break;

      case Function::ABS_OUTPUT:
      {
        // J = |u^T l| = sqrt( <u_R,l>^2 + <u_I,l>^2)
        // S lambda = - <u^*, l>/ 2*J * l = alpha * l

        complex<double> lu = l.Inner(u);// (u^*)^T l = <l, u>
        complex<double> ul = u.Inner(l);// J = |u^T l|
        complex<double> alpha = -0.5 * lu / std::abs(ul);

        LOG_DBG2(em) << "ACARHS: <u,l>=" << ul << " <l,u>=" << lu << " alpha=" << alpha;

        for(unsigned int i = 0; i < rhs.GetSize(); i++)
        rhs[i] = alpha * l[i];
        break;
      }

      case Function::CONJUGATE_COMPLIANCE: // rhs is from original excitation, we stored it in forward...rhs
      {
        forward.Get(excite)->Read(StateSolution::RHS_VECTOR, context->pde); // set
        Vector<complex<double> >& org_rhs = forward.Get(excite)->GetComplexVector(StateSolution::RHS_VECTOR);// read
        assert(org_rhs.GetSize() == u.GetSize());

        // the actual rhs for the adjoint pde is org_rhs * conj(u) -> this is only stored in OLAS!!!
        for(unsigned int i = 0, n = rhs.GetSize(); i < n; i++)
        rhs[i] = -1.0 * org_rhs[i] * std::conj(u[i]);
        break;
      }

      case Function::GLOBAL_DYNAMIC_COMPLIANCE:
      // S lambda = -conj(u);
      for(int i = 0, n = u.GetSize(); i < n; ++i)
      rhs[i] = -1.0 * std::conj(u[i]);
      break;

      case Function::ENERGY_FLUX:
      // reuse the RHS vector space
      SetEnergyFluxVector(f, u, true, rhs);
      // 0.25 * j* omega * (Q-Q^T)^T * u^* where 0.25 is 0.5 from standard adjoint formulation and 0.5 from objective
      for(unsigned int i = 0, n = rhs.GetSize(); i < n; i++)
      rhs[i] = complex<double>(0, 0.25) * excite.GetOmega() * rhs[i];
      break;

      case Function::STRESS:
      case Function::STRESS_DENSITY:
      {
        StressConstraint<complex<double> > sc(&excite, f, this, &forward);
        sc.CalcAdjointRHS(rhs);
        break;
      }

      default:
      assert(true); // e.g. for ELEC_ENERGY the rhs is set in PiezoSIMP::ConstructAdjointRHS()
    }
    shared_ptr<BaseFeFunction> fe = context->pde->GetFeFunction(context->pde->GetNativeSolutionType());
    // we cannot easily set the rhs. Therefore we set it to 0 and add our own rhs
    fe->GetSystem()->InitRHS(fe->GetFctId());
    fe->GetSystem()->SetFncRHS(rhs, fe->GetFctId());
    assert(!(rhs.NormMax() == 0.0 && f->GetLocal() != NULL)); // globalized stuff might have zero adjoint!
    LOG_DBG2(em) << "CARHS<complex>: f=" << context->GetDriver()->GetActStep(context->pde->GetName()) << " rhs before solving: " << rhs.ToString(1);
  }

  void ErsatzMaterial::ConstructAdjointRHS(Excitation& excite, Function* f)
  {
    assert(context->sequence == excite.sequence);
    assert(f->ctxt == context);

    // cannot be inlined due to linker problems
    if(context->IsComplex())
      ConstructComplexAdjointRHS(excite, f);
    else
      ConstructRealAdjointRHS(excite, f);
  }

// template instantiation stuff
template double ErsatzMaterial::CalcU1KU2<double>(TransferFunction* tf, StdVector<SingleVector*>& u1,  App::Type app, StdVector<SingleVector*>& u2,
    DesignDependentRHS* rhs, double factor, CalcMode calcMode, Function* f, int res_idx, double ev);

template double ErsatzMaterial::CalcU1KU2<complex<double> >(TransferFunction* tf, StdVector<SingleVector*>& u1,
    App::Type app, StdVector<SingleVector*>& u2,
    DesignDependentRHS* rhs, double factor, CalcMode calcMode, Function* f,  int res_idxm, double ev);

} // end of namespace
