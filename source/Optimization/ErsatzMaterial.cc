#include <string>
#include <cmath>
#include <limits>
#include <iomanip>
#include <fstream>

#include "Optimization/SIMP.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/DesignStructure.hh"
#include "Optimization/Design/DensityFile.hh"
#include "Optimization/OptimizationMaterial.hh"
#include "Domain/domain.hh"
#include "Domain/surfElem.hh"
#include "General/exception.hh"
#include "Utils/basenodestoresol.hh"
#include "PDE/SinglePDE.hh"
#include "PDE/mechPDE.hh"
#include "PDE/heatCondPDE.hh"
#include "PDE/acousticPDE.hh"
#include "Forms/baseForm.hh"
#include "Forms/linSurfForm.hh"
#include "Elements/basefe.hh"
#include "CoupledPDE/DirectCoupledPDE.hh"
#include "Utils/StdVector.hh"
#include "MatVec/vector.hh"
#include "Utils/result.hh"
#include "Utils/mathParser/mathParser.hh"
#include "Driver/assemble.hh"
#include "Driver/baseSolveStep.hh"
#include "Driver/harmonicDriver.hh"
#include "Driver/transientdriver.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/programOptions.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "OLAS/algsys/basesystem.hh"
#include "OLAS/algsys/baseentrymanipulator.hh"
#include "MatVec/basematrix.hh"
#include "MatVec/stdmatrix.hh"
#include "Materials/mechanicMaterial.hh"

using namespace CoupledField;
using namespace std;

DECLARE_LOG(conditions)
DEFINE_LOG(conditions, "conditions")

DECLARE_LOG(em)
DEFINE_LOG(em, "ersatzMaterial")

Enum<ErsatzMaterial::Method> ErsatzMaterial::method;

ErsatzMaterial::ErsatzMaterial() :
  Optimization(),
  material(NULL), // to be set in PostInit()
  dim(grid->GetDim())
{
  /** We store here the solution */
  volume_fraction_ = 0.0;
  structure_ = NULL;
  pde = NULL;
  densityFile = NULL;

  pn = param->Get("optimization")->Get("ersatzMaterial");

  method_ = method.Parse(pn->Get("method")->As<std::string>());

  homogenization_ = false;
  for(unsigned int i = 0; i < objectives.data.GetSize(); i++)
    if(objectives.data[i]->IsHomogenization()) homogenization_ = true;

  for(unsigned int i = 0; i < constraints.all.GetSize(); i++)
    if(constraints.all[i]->IsHomogenization()) homogenization_ = true;

  optInfoNode->Get(ParamNode::HEADER)->Get("homogenization")->SetValue(homogenization_);

  homogenizedTensor.Resize(dim == 2 ? 3 : 6, dim == 2 ? 3 : 6);
  homogenizedTensor.Init();

  // region stuff
  ParamNodeList region_list = pn->GetList("region");
  if(region_list.IsEmpty() && (method_ != SHAPE_OPT && method_ != SHAPE_PARAM_MAT))
    EXCEPTION("no region given!");

  for(unsigned int i=0; i < region_list.GetSize(); i++){
    std::string reg = region_list[i]->Has("name") ? region_list[i]->Get("name")->As<std::string>() : region_list[i]->As<std::string>();
    if(!grid->GetRegion().IsValid(reg))
      throw Exception("region given in ersatzMaterial is invalid");
    regionIds.Push_back(grid->GetRegion().Parse(reg));
  }


  // set up the design space elements, note PiezoSIMP have only POLARIZATION
  // this includes the transfer functions!
  ParamNodeList design_list = pn->GetList("design");
  ParamNodeList transfer_list = pn->GetList("transferFunction");
  ParamNodeList result = pn->GetList("result");
  design = DesignSpace::CreateInstance(regionIds, design_list, transfer_list, result, method_);
  // make basic loggings
  design->ToInfo(optInfoNode->Get(ParamNode::HEADER)->Get("designSpace"));

  // the L-mesh of the stress constraint benchmark is meshed by gid with different positions of
  // element nodes, such that one cannot use the same element materix, even if the grid is regular
  // therefore the attribute enforce_unstructured
  assume_constant_element_matrices_ = design->IsRegular() && method_ != ErsatzMaterial::PARAM_MAT
                                      && !pn->Get("enforce_unstructured")->As<bool>();
  LOG_TRACE2(em) << "EM:EM: const_mat=" << assume_constant_element_matrices_ << " reg=" << design->IsRegular()
                 << " PARAM_MAT=" << (method_ == ErsatzMaterial::PARAM_MAT) << " enforce_unstr=" << pn->Get("enforce_unstructured")->As<bool>();

  // optionally write the densities to an xml file
  if(pn->Has("export"))
    densityFile = new DensityFile(design, pn->Get("export"),
                                  design_list,
                                  transfer_list,
                                  pn->Get("SIMP/regularization", ParamNode::PASS));

  // we assume always to have either a mechanic or heatcond pde.
  SetPDEs();

  harmonic = BasePDE::IsComplex(pde->GetAnalysisType());

  if(homogenization_ && !pde->HasPeriodicBC())
    EXCEPTION("homogenization requires periodic boundary conditions");

  // Get the assemble class
  assemble_ = pde->getPDE_assemble();

  // check our constraints, the shall have only valid designs
  for(unsigned int i = 0; i < constraints.all.GetSize();  i++)
  {
    Condition* g = constraints.all[i];
    if(g->design == DesignElement::UNITY && (method_ == SHAPE_OPT || method_ == SHAPE_PARAM_MAT))
      continue;

    if(g->design == DesignElement::TENSOR_TRACE && (method_ == PARAM_MAT || method_ == SHAPE_PARAM_MAT))
      continue;

    if(design->FindDesign(g->design, false) == -1)
      throw Exception("constraint " + g->ToString() + " operates on invalid design variable");

  }

  // give the domain this data, s.th. the ersatz material approach is applied
  domain->SetErsatzMaterial(design);

  // add optimization results to the pde
  design->AppendOptimizationResults(pde);

  // forward and adjoint are initialized in PostInit()
  forward.Init(this);
  forward.SetIsForward(true);
  adjoint.Init(this);

  // check for multiple loadcases (might be frequencies)
  me->PrepareMultipleExcitations(pde, optInfoNode, harmonic, optimizer_ == EVALUATE_INITIAL_DESIGN);
}

ErsatzMaterial::~ErsatzMaterial()
{
  // if write to file close the xml envelope and the file
  if(densityFile != NULL) { delete densityFile; densityFile = NULL; }

  delete material;

  delete structure_;

  // "remove" the ersatzMaterial (=data) from the domain
  domain->SetErsatzMaterial(NULL);
}

void ErsatzMaterial::PostInit()
{
  // might be constructed in SIMP::PostInit()
  if(structure_ == NULL)
    structure_ = new DesignStructure(this);

  // post init slope constraints when the design is there
  constraints.PostProc(design, structure_, GetMultipleExcitation());
  // same for the objectives
  objectives.PostProc(design, structure_, GetMultipleExcitation());

  // the constraints size is only now known and the shapeDesign constructor is finished -> PostInit design
  design->PostInit(objectives.data.GetSize(), constraints.all.GetSize());


  // number of design variables in shape optimization is only known here
  // make a copy of the old iteration to calculate the move
  last_iteration.Resize(design->GetNumberOfVariables());

  // note the difference between function evaluations (line search) and iterations!
  last_evaluation.Resize(design->GetNumberOfVariables());

  for(unsigned int i = 0; i < objectives.data.GetSize(); i++)
  {
    Objective* cost = objectives.data[i];

    std::string objective = "'" + cost->type.ToString(cost->GetType()) + "'";

    // checks
    switch(cost->GetType())
    {
      case Objective::COMPLIANCE:
      case Objective::OUTPUT: // it would work but is saver not to allow
      case Objective::TRACKING:
      case Objective::HOMOGENIZATION_TENSOR:
      case Objective::HOMOGENIZATION_TRACKING:
      case Objective::POISSONS_RATIO:
      case Objective::YOUNGS_MODULUS:
      case Objective::ELEC_ENERGY: // it simply does not work yet in the harmonics
        if(harmonic) throw Exception(objective + " is only for static state problems");
        break;

      case Objective::DYNAMIC_OUTPUT:
      case Objective::GLOBAL_DYNAMIC_COMPLIANCE:
      case Objective::ABS_DYN_OUTPUT_SQUARED:
      case Objective::ENERGY_FLUX:
        if(!harmonic) throw Exception(objective + " is only for harmonic state problems");
        break;

      default:
        break;
    }

    // actions
    switch(cost->GetType())
    {
      case Objective::OUTPUT:
      case Objective::DYNAMIC_OUTPUT:
      case Objective::ABS_DYN_OUTPUT_SQUARED:
      {
        PtrParamNode output = cost->pn->Get("output");

        if(output->Has("displacement"))
          pde->ReadLoads(output->GetList("displacement"), output_nodes_);
        if(output->Has("elecPotential"))
          domain->GetSinglePDE("electrostatic")->ReadLoads(output->GetList("elecPotential"), output_nodes_);
        if(output->Has("acoustic"))
          domain->GetSinglePDE("acoustic")->ReadLoads(output->GetList("acoustic"), output_nodes_);

        if(output_nodes_.GetSize() == 0)
          throw Exception("no output optimization targets given");
        break;
      }
      // we do energy flux in realtime
      default:
        break;
    }
  }

  // create Material class
  if(pde != NULL)
  {
    OptimizationMaterial::System system = 
        OptimizationMaterial::system.Parse(pn->Get("material")->As<std::string>());
    switch(system)
    {
    case OptimizationMaterial::PIEZO:
      material = new OptPiezoMat(this);
      break;

    case OptimizationMaterial::MECHANIC:
      material = new OptMechMat(this);
      break;

    case OptimizationMaterial::HEAT:
      material = new HeatMat(this);
      break;

    default: assert(false);
    }
    optInfoNode->Get(ParamNode::HEADER)->Get("material")->SetValue(OptimizationMaterial::system.ToString(material->GetSystem()));
  }

  // if loadErsatzMaterial is used with optimization specifying a starting point,
  // we have to load it here, before scaling is done.
  if(DensityFile::NeedLoadErsatzMaterial())
    DensityFile::ReadErsatzMaterial(design);

  // plausibility check for homogenization
  if(homogenization_ && (!me->IsEnabled() || !me->DoHomogenization()))
    throw Exception("A homogenization objective/constraint is set but no homogenization test strain excitation");
  if(me->IsEnabled() && me->DoHomogenization() && !homogenization_)
    throw Exception("No homogenization objective/constraint for homogenization test strain excitation");

}




void ErsatzMaterial::StoreResults(double step_val)
{
  CommitMode cm = commitMode_;

  double real_step = step_val != -1 ? step_val : currentIteration;

  if(cm == EACH_FORWARD)
  {
    for(unsigned int e = 0; e < me->excitations.GetSize(); e++) {
      forward.Get(e)->Write(pde);
      // call real implementation in Optimization. sum up in excitation fractions up to smaller 0.5
      Optimization::StoreResults(real_step + (0.5 / me->excitations.GetSize()) * e);
    }
  }

  if(cm == FORWARD || cm == BOTH) {
    Solution::Write(pde, forward, NULL, 0, me->excitations); // no function for forward, no time step
    Optimization::StoreResults(step_val);
  }

  // over all functions, mostly only one or none
  StdVector<Function*> funcs = adjoint.GetFunctions();
  for(unsigned int fi = 0; fi < funcs.GetSize(); fi++)
  {
    LOG_DBG(em) << "StoreResults(" << step_val << ") rs=" << real_step << " adjoint_function=" << funcs[fi]->ToString();
    if(cm == EACH_ADJOINT)
    {
      for(unsigned int e = 0; me->excitations.GetSize(); e++) {
        adjoint.Get(e, funcs[fi])->Write(pde);
        // call real implementation in Optimization. sum up in excitation fractions up to smaller 0.5
        double index = (me->excitations.GetSize() * funcs.GetSize()) * (fi * funcs.GetSize()) * e;
        Optimization::StoreResults(real_step + 0.5 + (0.5 / index));
      }
    }

    if(cm == ADJOINT || cm == BOTH) {
      // sum up if there are more excitations
      Solution::Write(pde, adjoint, funcs[fi], 0, me->excitations); // TODO no time step set!
      Optimization::StoreResults(real_step + 0.5 + (0.5 / funcs.GetSize()) * fi);
    }
  }
}


PtrParamNode ErsatzMaterial::CommitIteration(bool keep_iteration_number)
{
  // will write the cfs results and the log file
  PtrParamNode iter = Optimization::CommitIteration(keep_iteration_number);
  // add our multiple excitation stuff here (only in info.xml, this would be to complex for dat
  if(me->IsEnabled())
  {
    for(unsigned int e = 0; e < me->excitations.GetSize(); e++)
    {
      Excitation& excite = me->excitations[e];
      PtrParamNode info = iter->Get("excitation", ParamNode::APPEND);
      info->Get("index")->SetValue(excite.index);
      info->Get("objective")->SetValue(excite.cost);
      info->Get("normalized_weight")->SetValue(excite.normalized_weight);
      for(unsigned int c = 0; c < constraints.all.GetSize(); c++)
      {
        Condition* g = constraints.all[c];
        if(g->IsExcitationSensitive() && g->DoEvaluate(&excite))
          info->Get(g->ToString(me))->SetValue(g->GetValue());
      }
    }
  }

  if(homogenization_)
  {
    PtrParamNode in = iter->Get("homogenizedTensor");
    in->Get("norm_L2")->SetValue(homogenizedTensor.NormL2());

    PtrParamNode iso = in->Get("isotropy");
    iso->Get("err")->SetValue(MechanicMaterial::CalcIsotropyError(homogenizedTensor, false));
    iso->Get("err_normed")->SetValue(MechanicMaterial::CalcIsotropyError(homogenizedTensor, true));
    iso->Get("poissons_ratio")->SetValue(MechanicMaterial::CalcIsotropicPoissonsRatio(homogenizedTensor));
    iso->Get("E")->SetValue(MechanicMaterial::CalcIsotropicYoungsModulus(homogenizedTensor));

    PtrParamNode orth = in->Get("orthotropy");
    StdVector<std::pair<string, double> > ortho = MechanicMaterial::CalcOrthotropeProperties(homogenizedTensor);
    for(unsigned int i = 0; i < ortho.GetSize(); i++)
      orth->Get(ortho[i].first)->SetValue(ortho[i].second);

    in->Get("tensor")->SetValue(new Matrix<double>(homogenizedTensor));
  }

  if(densityFile != NULL)
    densityFile->SetCurrent(currentIteration - 1); // the counter was already incremented!

  return iter;
}



string ErsatzMaterial::GetIterationFrequency()
{
  if(!harmonic) return "";

  // make clear, when doing *real* multiple excitations, that this is no single
  // frequency result
  if(me->IsEnabled() && me->excitations.GetSize() > 1) return "(mult)";

  double frequency = dynamic_cast<HarmonicDriver*>(domain->GetDriver())->GetActFreq();
  // as we control the fractional digits, we do not use lexical_cast<string>
  stringstream ss;
  ss << fixed << std::setprecision(1) << frequency;
  return ss.str();
}


BiLinFormContext* ErsatzMaterial::GetFormContext(const RegionIdType regionId, StdPDE* pde1, StdPDE* pde2, const std::string& integrator)
{
  return(assemble_->GetBiLinForm(regionId, pde1, pde2, integrator));
}

BaseForm* ErsatzMaterial::GetForm(const RegionIdType regionId, StdPDE* pde1, StdPDE* pde2, const std::string& integrator)
{
  return(GetFormContext(regionId, pde1, pde2, integrator)->GetIntegrator());
}

int ErsatzMaterial::GetSpecialResultIndex(Application app1, Application app2, CalcMode calcMode, Condition* constraint)
{
  stringstream label;
  label << application.ToString(app1) << "_" << application.ToString(app2);
  if(calcMode == CONJ_QUAD) label << "_quad";

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

void ErsatzMaterial::CalcNewmarkDerivative(Excitation& excite, Solutions& forward, Solutions& adjoint, double factor, Objective* c, Condition* g){
  // this calculates p^T (dF - dA) u
  // where p is solution of adjoint, dF is derivative of newmark update, dA is derivative of system matrix, u is solution of forward problem
  
  Function* f = Function::Cast(c, g);

  UInt timesteps = domain->GetDriver()->GetNumSteps();
  double dt = dynamic_cast<TransientDriver*>(domain->GetDriver())->GetDeltaT();
  double gamma = pde->getTimeStepping()->GetNewmarkGamma();
  double beta = pde->getTimeStepping()->GetNewmarkBeta();
  
  MathParser * parser = domain->GetMathParser();
  unsigned int mathParserHandle = parser->GetNewHandle();

  assert(domain->HasErsatzMaterialTensor()); // this is not implemented for SIMP
  
  Matrix<double> dK(1, 1), dM(1, 1);
  Vector<double> dKu(0), dMu(0);

  // the outer most loop is over all elements, so element matrices can be reused as much as possible
  int upper = design->data.GetSize();
  int elements = design->GetNumberOfElements();
  for(int base = 0; base < upper; base += elements) { // loop over all designs
    for(int e = 0 ; e < elements; ++e) { // loop over all elements
      DesignElement* de = &design->data[base + e];
      SetElementK(de, MECH, dynamic_cast<DenseMatrix*>(&dK), STANDARD, true);
      SetElementK(de, MASS, dynamic_cast<DenseMatrix*>(&dM), STANDARD, true);
      
      // The damping matrix is alpha * Mass + beta * Stiffness, so it's derivative is also alpha * dMass + beta * dStiffness
      // We need to get alpha and beta, from the integrators
      RegionIdType regionId = de->elem->regionId;
      BiLinFormContext* linElastIntCtxt = GetFormContext(regionId, pde, pde, "linElastInt");
      BiLinFormContext* linMassIntCtxt = GetFormContext(regionId, pde, pde, "MassInt");
      double dampingAlpha = 0.0; double dampingBeta = 0.0;
      if(linElastIntCtxt->GetSecDestMat() != NOTYPE){
        parser->SetExpr(mathParserHandle, linElastIntCtxt->GetSecMatFac());
        dampingBeta = parser->Eval(mathParserHandle);
      }
      if(linMassIntCtxt->GetSecDestMat() != NOTYPE){
        parser->SetExpr(mathParserHandle, linMassIntCtxt->GetSecMatFac());
        dampingAlpha = parser->Eval(mathParserHandle);
      }
      
      double vK = 0.0;
      double vM = 0.0;
      double vC = 0.0;
      for(unsigned int t = 0; t < timesteps; ++t){ // loop over all time steps in u
        Vector<double>& u_vec = dynamic_cast<Vector<double>& >(*forward.Get(excite, NULL, t)->elem[MECH][e]);
        dKu = dK * u_vec;
        dMu = dM * u_vec;
        // the dA part
        Vector<double>& p_vecd = dynamic_cast<Vector<double>& >(*adjoint.Get(excite, f, t)->elem[MECH][e]);
        double dvK = p_vecd * dKu;
        vK -= dvK;
        double dvM = p_vecd * dMu;
        vM -= dvM;
        double dvC = (gamma / (beta * dt) ) * ( dampingAlpha * dvM + dampingBeta * dvK);
        vC -= dvC;
        double u = 1.0; double upp = 1.0 / (beta*dt*dt); double up = upp * gamma * dt;
        if(t == 0 && IsFirstTransientStepStatic()){ // reset up and upp as in simulation (StdSolveStep)
          upp = 0.0;
          up = 0.0;
          vM = 0.0; // reset the mass part for the first step
          vC = 0.0;
        }
        for(unsigned int tp = t+1; tp < timesteps; ++tp){ // loop over all time steps in p
          Vector<double>& p_vec = dynamic_cast<Vector<double>& >(*adjoint.Get(excite, f, tp)->elem[MECH][e]);
          // these are the scalar factors for the parts of the derivative of F
          double ut = u +  (up + 0.5*upp*(1.0-2.0*beta)*dt ) *dt;
          double upt = up + (1.0 - gamma) * dt * upp;
          double pdMu = p_vec * dMu;
          double tvM = ut * pdMu;
          vM += tvM;
          double pdKu = p_vec * dKu;
          double tvC = ( gamma * ut / (beta * dt) - upt ) * (dampingAlpha * pdMu + dampingBeta * pdKu);
          vC += tvC;
          u = 0.0;
          upp = (u - ut) / (beta * dt * dt);
          up = (upt + upp * gamma * dt);
        }
      }
      de->AddGradient(c, g, factor * (vK + vM / (beta * dt * dt) + vC) );
    }
  }
}

double ErsatzMaterial::CalcU1KU2(TransferFunction* tf, StdVector<SingleVector*>& u1,
    Application k, StdVector<SingleVector*>& u2, SurfaceRef* rhs,
    double factor, CalcMode calcMode, Function* f, int res_idx)
{
  if (harmonic)
    return CalcU1KU2<std::complex<double> > (tf, u1, k, u2, rhs, factor, calcMode, f, res_idx);
  else
    return CalcU1KU2<double> (tf, u1, k, u2, rhs, factor, calcMode, f, res_idx);
}

template <class T>
double ErsatzMaterial::CalcU1KU2(TransferFunction* tf, StdVector<SingleVector*>& u1,
                       Application app, StdVector<SingleVector*>& u2,
                       SurfaceRef* rhs, double factor, CalcMode calcMode, Function* f, int res_idx)
{
  LOG_DBG2(em) << "CalcU1KU2(): tf=" << (tf ? tf->ToString() : "NULL") << " #u1=" << u1.GetSize()
                 << " app=" << application.ToString(app) << " #u2=" << u2.GetSize() << " calcMode="
                 << calcMode << " factor=" << factor << " rhs=" << (rhs == NULL ? "NULL" : rhs->ToString(1));

  // This solves <l,K'*u-f'> or <u1, K' * u2 - f'> for all elements and adds it up to the element gradients

  assert(u1.GetSize() != 0);
  assert(u1.GetSize() == u2.GetSize());

  double sum = 0.0;
  // mat will be filled by SetElementK where also the derivative form most cases is built in
  // the dimensions of our matrix is determined by u1_vec and u2_vec.
  Matrix<T> mat(u1[0]->GetSize(), u2[0]->GetSize()); //NOTE: SetElementK (In PiezoSimp) relies on the matrix already having the right size!!!
  Vector<T> mat_vec(u1[0]->GetSize());

  TransferFunction* rtf = rhs != NULL && rhs->valid ? design->GetTransferFunction(tf->GetDesign(), rhs->app) : NULL;

  // traverse over our elements
  // in ErsatzMaterialTensor case we loop over all elements, else only over the elements belonging to this design
  int elements = design->GetNumberOfElements();
  int base_lower = 0;
  int base_upper = 0;
  if(domain->HasErsatzMaterialTensor()){
    base_upper = design->data.GetSize();
  }else{
    base_lower = design->FindDesign(tf->GetDesign()) * elements;
    base_upper = base_lower + elements;
  }
  LOG_DBG2(em) << "elements=" << elements << " base=" << base_lower << " base_upper=" << base_upper;

  // create an element list to gain the iterator in the loop
  ElemList elemList(grid);

  // for ParamMat we need the derivative w.r.t. every designvariable, else the base loop is only run once
  for(int base = base_lower; base < base_upper; base += elements)
  {
    for(int e = 0 ; e < elements; e++)
    {
      Vector<T>& u1_vec = dynamic_cast<Vector<T>& >(*u1[e]);
      Vector<T>& u2_vec = dynamic_cast<Vector<T>& >(*u2[e]);

      DesignElement* de = &design->data[e + base];

      LOG_DBG3(em) << "nodes:" << e << ": " << de->elem->connect.ToString();
      LOG_DBG3(em) << "u1:" << e << ": " << u1_vec.ToString();
      LOG_DBG3(em) << "u2:" << e << ": " << u2_vec.ToString();

      // u1^T (K' u2 - f') -> find "K'"
      SetElementK(de, app, dynamic_cast<DenseMatrix*>(&mat), calcMode);
      LOG_DBG3(em) << "mat: " << mat.ToString();

      // We generally solve u1^T (K' u2 - f') as u1^T (K' u2 - f')
      // u1^T (K' u2 - f') -> calc K' u2"
      mat_vec = mat * u2_vec;
      LOG_DBG3(em) << "mat * u2: " << mat_vec.ToString();

      // u1^T (K' u2 - f') -> calc "- f'"
      assert(!(calcMode == CONJ_QUAD && rtf != NULL)); // no sensitive rhs here!
      if(rtf != NULL) SubstractGradSurfaceRHS(de, rtf, rhs, mat_vec);
      LOG_DBG3(em) << "-f': " << mat_vec.ToString();

      // u1^T(K' u2 - f') -> calc "u1^T *" or <u1, *> 
      // the difference is the conjugate complex in the harmonic inner product case!
      T sp;
      if(calcMode == CONJ_QUAD) mat_vec.Inner(u1_vec, sp); // u1 = u2 = u! 
                           else sp = mat_vec * u1_vec ;

      // when doing complex Jensen 22.07.07 shows that we always have 2 * Re(lamda * grad S * u)
      // the factor gives the negative sign
      // in real case it is simple value = factor * sp.
      // factor shall be +/- 1!
      double this_value = factor;
      if(harmonic && calcMode != CONJ_QUAD) this_value *= 2 * ((complex<double>) sp).real(); // 2 * Re{...}
                                      else this_value *= ((complex<double>) sp).real(); // CONJ_QUAD or real STANDARD

      de->AddGradient(f, this_value);

      LOG_DBG3(em) << "CU1Ku2:" << de->elem->elemNum << " <l,K'*u-f'>  = "
                   << sp << " -> " << this_value << " sum = " << de->GetPlainGradient(f);

      sum += this_value;

      if(res_idx != -1) de->specialResult[res_idx] = this_value;
    }
  }

  return sum;
}


template <class T>
void ErsatzMaterial::SubstractGradSurfaceRHS(DesignElement* de, TransferFunction* tf, SurfaceRef* ref, Vector<T>& in_out)
{
  // we have to find the nodes which are common between de->elem
  // and the surface element which is one dimension smaller

  // not all elements do necessary lay on a surface and then not all nodes
  assert(ref != NULL && ref->valid);


  // nodes (numbers) of our design element
  StdVector<unsigned int>& de_nodes = de->elem->connect;
  Vector<T>& rhs = dynamic_cast<Vector<T>& >(*(ref->vec));

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
      LOG_DBG3(em) << "SubstractGradSurfaceRHS : node " << n << " is common with elem "
                     << *it << " in surface: K'u = " << in_out.ToString();

      // find the the sensitivity of the rhs w.r.t the design volume element!
      double factor = tf->Derivative(de);

      // we do not really construct a rhs vector (with some/many zeros) but substract
      // for all design nodes common with the surface directly the entries
      for(int d = 0; d < dof; d++)
        in_out[n*dof + d] -= factor * rhs[d];
      LOG_DBG3(em) << "... -" << factor << "*" << rhs.ToString()
                     << " -> " << in_out.ToString();
    }
  }
}

double ErsatzMaterial::CalcObjective()
{
  design->WriteDesignToExtern(last_evaluation.GetPointer());

  // in objective.value_ we store the sum over all excitations w/o penalty but with normalization
  // in excitation.cost we store the sum over all objectives with penalty but w/o normalization

  // reset the objective values such that we can sum up normalized but unpenalized values
  for(unsigned int o = 0; o < objectives.data.GetSize(); o++)
    objectives.data[o]->ResetValue();

  double result = 0.0;

  // the multiple excitation case is a special case - for all other cases this is executed once
  for(unsigned int e = 0; e < me->excitations.GetSize(); e++)
  {
    Excitation& excite = me->excitations[e];
    excite.cost = 0.0;

    for(unsigned int o = 0; o < objectives.data.GetSize(); o++)
    {
      Objective* f = objectives.data[o];

      // some objectives are only to be evaluated for the last excitation
      if(!f->DoEvaluate(&excite)) continue;

      //double ov = CalcObjective(excite, f); // this is virtual!
      double ov = CalcFunction(excite, f, false); // this is virtual!
      excite.cost += ov * f->GetPenalty();

      // we ignore the weight if the evaluation happens only once! TODO why not omega*omega? - Fabian
      double weight = !f->DoEvaluateAlways() ? 1.0 : excite.normalized_weight;

      f->AddValue(ov * weight);

      result += ov * f->GetPenalty() * weight;
      LOG_DBG(em) << "CalcObjective: ex=" << e << " obj=" << f->type.ToString(f->GetType()) << " ov=" << ov
                  << " penalty" << f->GetPenalty() << " ex.cost=" << excite.cost << " nw=" << excite.normalized_weight
                  << " wei=" << weight << " f->val=" << f->GetValue() << " result=" << result;
    }
  }

  return result;
}

void ErsatzMaterial::CalcObjectiveGradient(StdVector<double>* grad_out)
{
  // reset the cost gradients in the design elements and sum them up in a weighted way
  // to perform multiple loads
  design->Reset(DesignElement::COST_GRADIENT);

  for(unsigned int obj = 0; obj < objectives.data.GetSize(); obj++)
  {
    Objective* cost = objectives.data[obj];
    // the multiple excitation case is a special case - for all other cases this is executed once
    for(unsigned int idx = 0; idx < me->excitations.GetSize(); idx++)
    {
      Excitation& excite = me->excitations[idx];
      // some objectives are only to be evaluated for the last excitation
      if(!cost->DoEvaluate(&excite)) continue;

      CalcFunction(excite, cost, true);
    }
  }

  if(grad_out != NULL)
    design->WriteGradientToExtern(*grad_out, DesignElement::COST_GRADIENT, DesignElement::SMART);
}



double ErsatzMaterial::CalcConstraint(Condition* g)
{
  // assume when we have only one constraint which is not explicitly given, this is not the stress constraint!
  assert((g == NULL && constraints.active.GetSize() == 1 && constraints.active[0]->DoEvaluateAlways()) || g != NULL);

  if(g == NULL)
    g = constraints.active[0];

  double result = 0.0;

  for(unsigned int e = 0; e < me->excitations.GetSize(); e++)
  {
    Excitation& excite = me->excitations[e];
    // in the evaluate once case only the last excitation
    double v = g->DoEvaluate(&excite) ? CalcFunction(excite, g, false) : 0.0;
    double w = g->DoEvaluateAlways() ? excite.GetFactor(g) : 1.0;
    result += v * w;
    LOG_DBG(em) << "CC ex=" << e << " eval=" << g->DoEvaluate(&excite) << " v=" << v << " w=" << w << " -> " << result;
  }

  g->SetValue(result);
  return result;
}

void ErsatzMaterial::CalcConstraintGradient(Condition* g, StdVector<double>* grad_out)
{
  // assume when we have only one constraint which is not explicitly given, this is not the stress constraint!
  assert((g == NULL && constraints.active.GetSize() == 1 && !constraints.active[0]->DoEvaluateAlways()) || g != NULL);

  if(g == NULL)
    g = constraints.active[0];

  for(unsigned int i = 0; i < me->excitations.GetSize(); i++)
  {
    if(g->DoEvaluate(&me->excitations[i]))
      CalcFunction(me->excitations[i], g, true);
  }

  // copies from the design element gradient data to a memory array for external optimizers
  if(grad_out != NULL)
    design->WriteGradientToExtern(*grad_out, DesignElement::CONSTRAINT_GRADIENT, DesignElement::SMART, g);

  // if there is a <result ... value="constraintGradient" detail="penalizedVolume/*"
  if(g->special_result_idx != -1)
  {
    int base = design->FindDesign(g->design);
    int n    = design->GetNumberOfElements();
    for(int i = n * base; i < n * (base + 1); i++) // TODO add access!
      design->data[i].specialResult[g->special_result_idx] = design->data[i].GetPlainGradient(NULL, g);
  }
}



double ErsatzMaterial::CalcFunction(Excitation& excite, Function* f, bool derivative)
{
  assert(f != NULL);

  // for legacy reasions there is also the difference between Objective and Condition, to be replaced once
  Objective* c = f->IsObjective() ? dynamic_cast<Objective*>(f) : NULL;
  Condition* g = f->IsObjective() ? NULL : dynamic_cast<Condition*>(f);

  double result = 0.0; // stays  for the derivative

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

    case Function::GREYNESS:
         assert(c == NULL);
         result = CalcGreyness(g, derivative);
         break;

    case Objective::STRESS: {
           const Vector<double> data = CalcVonMisesStressVector(excite, f, false, false); // copy data for element von Mises stress
           result = CalcGlobalFunction(f, false, &data);
           break;
         }

    case Function::GLOBAL_SLOPE:
    case Function::GLOBAL_MOLE:
    case Function::GLOBAL_OSCILLATION:
    case Function::GLOBAL_JUMP:
         result = CalcGlobalFunction(f, derivative);
         break;

    case Function::SLOPE:
    case Function::MOLE:
    case Function::OSCILLATION:
    case Function::JUMP:
         assert(c == NULL);
         result = CalcLocalConstraint(g, derivative);
         break;

    case Function::HOMOGENIZATION_TENSOR:
          if(c != NULL && derivative && c->HasHomogenizationEntry())
          {
            // if there s no "coord" set it is only meant for evaluateInitialDesign for forward homogenization
            StdVector<double> tmp;
            CalcHomogenizedTensorEntry(c->coord, true, tmp);
            for(unsigned int e = 0, ne = design->GetNumberOfElements(); e < ne; e++)
              design->data[e].AddGradient(c, NULL, tmp[e]);
          }
          if(c != NULL && !derivative)
          {
            Matrix<double> hom_tensor = CalcHomogenizedTensor();
            if(c->HasHomogenizationEntry())
            {
              result = hom_tensor[c->coord.first -1][c->coord.second - 1];
            }
            else
            {
              std::cout << "Homogenized Tensor: " << std::endl << hom_tensor.ToString(0, true);

              std::cout << "Orthotrope properties: ";
              StdVector<std::pair<string, double> > ortho = MechanicMaterial::CalcOrthotropeProperties(hom_tensor);
              for(unsigned int i = 0; i < ortho.GetSize(); i++)
                std::cout << " " << ortho[i].first << "=" << ortho[i].second;
              std::cout << "\n";

              result = hom_tensor.NormL2();
            }
          }
          if(g != NULL)
          {
            result = CalcHomogenizedTensorConstraint(g, derivative);
          }
          break;

    case Function::HOMOGENIZATION_TRACKING:
         if(derivative)
         {
           CalcHomogenizedTrackingGradient(f->GetTensor(), CalcHomogenizedTensor(), c, g);
         }
         else
         {
           double diff = f->GetTensor().DiffNormL2(CalcHomogenizedTensor());
           result = 0.5 * diff * diff;
         }
         break;

    case Function::POISSONS_RATIO:
    case Function::YOUNGS_MODULUS:
          result = CalcPoissonsRatioAndYoungsModulus(c, g, derivative);
          break;

    case Function::GLOBAL_DYNAMIC_COMPLIANCE:
          assert(!derivative); // SIMP!
          result = CalcGlobalDynamicCompliance(excite, c);
          break;

     case Function::OUTPUT:
     case Function::DYNAMIC_OUTPUT:
     case Function::CONJUGATE_COMPLIANCE:
     case Function::ABS_DYN_OUTPUT_SQUARED:
          assert(!derivative); // SIMP!
          if(harmonic)
            result = CalcOutputObjective<complex<double> >(excite, c);
          else
            result = CalcOutputObjective<double>(excite, c);
          break;

     case Function::ENERGY_FLUX:
           assert(!derivative);
           result = CalcEnergyFlux(excite, c);
           break;

     case Function::TEMPERATURE:
           break; // FIXMEHEAT

     case Function::ELEC_ENERGY:
       assert(false); // shall be handled before

     case Function::ISOTROPY:
     case Function::MULTI_OBJECTIVE:
       assert(false); // no valid function
    // no default, gcc warns
  }

  LOG_DBG2(em) << "CalcFunction " << f->ToString() << " cost=" << f->IsObjective() << " -> " << (derivative ? "derivative" : lexical_cast<std::string>(result));
  return result;
}


double ErsatzMaterial::IntegrateDesignVariable(Objective* f, Condition* g, bool derivative,
    DesignElement::Type dtype, bool normalized, bool scale, double exponent)
{
  // this replaces and enhances calculation of volume, it is used by regularization
  // when not assuming a regular grid, computation of Volume is not as simple
  // we also consider working only on a given region, when used as constrain
  // use dtype == NO_TYPE to iterate over all designs, but do not calculate tensor trace even if available

  // do we want the physical value? Don't make GTF() fault tolerant as we assume the physcial value!
  TransferFunction* tf = Function::GetFunction(f, g)->IsPhysical() ? design->GetTransferFunction(dtype, MECH) : NULL;

  Matrix<Double> cornerCoords;
  double fraction = f != NULL ? volume_fraction_ : g->volume_fraction;
  bool allDesignsRelevant = dtype == DesignElement::TENSOR_TRACE || dtype == DesignElement::DEFAULT || dtype == DesignElement::NO_TYPE;
  bool calculateTensorTrace = domain->HasErsatzMaterialTensor() && (dtype == DesignElement::TENSOR_TRACE || dtype == DesignElement::DEFAULT); // tensor trace is calculated if dtype == DEFAULT or TENSOR_TRACE and a tensor available
  if(calculateTensorTrace && scale){
    throw Exception("Cannot calculate Tensor Trace Volume on scaled design variables!");
  }
  // check whether we have to calculate the full volume
  if(normalized){
    if(g != NULL && g->design == DesignElement::UNITY){ // this will always return 1, does only make sense for unnormalized (real) volume
      // the gradient is not set, it is really 0 on every element but should not be used
      //TODO: Do We need a warning here? 
      return(derivative ? 0.0 : 1.0);
    }
  }else{
    if(design->IsRegular()){ // we use 1/(the volume of the first element) as fraction
      grid->GetElemNodesCoord(cornerCoords, design->data[0].elem->connect, true);
      fraction = 1.0 / design->data[0].elem->ptElem->CalcVolume(cornerCoords, false);
    }else{
      fraction = 1.0; // if no normalization needed, fraction is simply 1.0
    }
  }
  if( fraction == 0.0 ){
    for(unsigned int d = 0; d < design->design.GetSize(); d++){
      if(allDesignsRelevant || design->design[d] == dtype){
        const unsigned int base = d * design->elements;
        for(unsigned int r = 0; r < design->regions.GetSize(); r++){
          if(f != NULL || g->IsForRegion(design->regions[r].regionId)){
            if(design->IsRegular()){
              fraction += design->regions[r].elements;
            }else{
              const unsigned int u = base + design->regions[r].base + design->regions[r].elements;
              for(unsigned int i = base + design->regions[r].base; i < u; i++){
                DesignElement* de = &design->data[i];
                grid->GetElemNodesCoord(cornerCoords, de->elem->connect, true);
                fraction += de->elem->ptElem->CalcVolume(cornerCoords, false);
              }
            }
          }
        }
      }
    }
    fraction = 1.0 / fraction;
    if(g == NULL){
      volume_fraction_ = fraction;
    }else{
      g->volume_fraction = fraction;
    }
  }

  double sum = 0.0;

  for(unsigned int d = 0; d < design->design.GetSize(); d++){
    if(allDesignsRelevant || design->design[d] == dtype){
      const unsigned int base = d * design->elements;
      for(unsigned int r = 0; r < design->regions.GetSize(); r++){
        if(f != NULL || g->IsForRegion(design->regions[r].regionId)){
          const double scaling = scale ? design->scale_design[d][r] : 1.0;
          const double rscaling = 1.0 / scaling;
          const double translation = scale ? design->translate_design[d][r] : 0.0;
          const unsigned int u = base + design->regions[r].base + design->regions[r].elements;
          for(unsigned int i = base + design->regions[r].base; i < u; i++)
          {
            DesignElement* de = &design->data[i];
            // standard or derivative case?
            if(derivative)
            {
              double val = 0.0;
              if(calculateTensorTrace){
                Matrix<double> material;
                GetErsatzMaterialTensor(material, de->elem, de->GetType());
                material.Trace(val);
                if(exponent != 1.0){
                  double des;
                  GetErsatzMaterialTensor(material, de->elem);
                  material.Trace(des);
                  val *= exponent * std::pow(des, exponent - 1.0);
                }
              }else{
                if(exponent != 1.0){
                  // by tf we mark if we want the physical stuff.
                  double des_val = tf != NULL ? tf->Derivative(de, DesignElement::SMART) : de->GetDesign(DesignElement::SMART);
                  val = exponent * std::pow(((des_val - translation) * rscaling), exponent - 1.0);
                }
                else
                  val = 1.0;
              }
              val *= rscaling;  // the gradient will be multiplied by scaling later again, but if it is calculated on the scaled designs here (as in regularization) this should not happen
              val *= fraction;
              if(!design->IsRegular()){
                grid->GetElemNodesCoord(cornerCoords, de->elem->connect, true);
                const double vol = de->elem->ptElem->CalcVolume(cornerCoords, false);
                val *= vol;
              }
              de->AddGradient(f, g, val);
            }
            else
            {
              // no derivative
              double des;
              if(calculateTensorTrace){ // use the trace of the stiffness Tensor as "volume"
                Matrix<double> material;
                GetErsatzMaterialTensor(material, de->elem);
                material.Trace(des);
              }
              else
              {
                // tf marks the physical function
                double des_val = tf != NULL ? tf->Transform(de, DesignElement::SMART) : de->GetDesign(DesignElement::SMART);
                des = (des_val - translation) * rscaling;
              }
              des = std::pow(des, exponent);
              if(design->IsRegular()){
                sum += des;
              }else{
                grid->GetElemNodesCoord(cornerCoords, de->elem->connect, true);
                const double vol = de->elem->ptElem->CalcVolume(cornerCoords, false);
                sum += des * vol;
              }
            } // if derivative
          } // for element (i)
        } // right region
      } // for region
      if(calculateTensorTrace && !derivative){  // the derivative has to be calculated for all designs but the tensor itself only once
        break;
      }
    } // if relevant
  } // for design
  return derivative ? -1.0 : sum * fraction;
}

double ErsatzMaterial::CalcVolume(Objective* f, Condition* g, bool derivative, bool normalized)
{
  // obly a constraint has a design type set (if it has)
  DesignElement::Type des  = g == NULL ? DesignElement::DEFAULT : g->design;

  // parameter is form the base function
  Function::Type func = f != NULL ? f->GetType() : g->GetType();
  double        exp  = f != NULL ? f->GetParameter() : g->GetParameter();

  switch(func)
  {
  case Function::VOLUME:
    return IntegrateDesignVariable(f, g, derivative, des, normalized, false, 1.0); // no scaling, exponent=1

  case Function::PENALIZED_VOLUME:
  case Function::REALVOLUME:
    return IntegrateDesignVariable(f, g, derivative, des, normalized, false, exp);

  case Function::GAP:
  {
    if(!derivative)
    {
      double vol = IntegrateDesignVariable(f, g, false, des, normalized, false, 1.0);
      double pen = IntegrateDesignVariable(f, g, false, des, normalized, false, exp);
      return vol - pen;
    }
    else
    {
      if(!design->IsRegular())
        throw Exception(Function::type.ToString(func) + " only implemented for regular grids");

      CalcRegularGapConstraint(f, g, des);
      return -1.0;
    }
  }

  default: assert(false);
  }
  return -1.0; // cannot happen due to assert
}

void ErsatzMaterial::CalcRegularGapConstraint(Objective* f, Condition* g, DesignElement::Type dt)
{
  assert(design->IsRegular());

  unsigned int des = design->FindDesign(dt);
  unsigned int ele = design->GetNumberOfElements();

  // exponent for penalized volume
  const double exp = f != NULL ? g->GetParameter() : g->GetParameter();

  for(unsigned int i = des * ele; i < (des + 1) * ele; i++)
  {
    DesignElement& de = design->data[i];

    // derivative of vol is 1
    // derivative of pen (x^p) is p*x^(p-1)
    double pen_grad = exp * std::pow(de.GetDesign(DesignElement::SMART), exp - 1.0);
    // normalize
    double grad = (1.0 - pen_grad) / (double) ele;

    de.AddGradient(f, g, grad);
  }
}


double ErsatzMaterial::CalcGlobalDynamicCompliance(Excitation& excite, Objective* f)
{
  //GLOBAL_DYNAMIC_COMPLIANCE
  // c = u^T conj(u) -> "A note on sensitivity analysis of linear dynamic systems with
  //                          harmonic excitation"; Jakob S. Jensen; June 22, 2007
  Vector<complex<double> >& u = forward.Get(excite)->GetComplexVector(Solution::RAW_VECTOR);
  assert(u.GetSize() != 0);

  complex<double> csp;
  u.Inner(u, csp);  // the inner product is sum over u_i * conj(u_i);
  assert(csp.imag() == 0);
  return csp.real() * excite.GetFactor(f);
}

double ErsatzMaterial::CalcCompliance(Excitation& excite, Objective* f, Condition* g, bool derivative)
{
  // note for the derivative case gradients are summed up (with weights)
  assert(f != NULL || g != NULL);
  assert(f == NULL || g == NULL);
  Function* func = f != NULL ? dynamic_cast<Function*>(f) : dynamic_cast<Function*>(g);
  double result = 0.0;
  if(derivative)
  {
    // calculate the compliance which is according to
    // "A 99 line topology optimization code written in Matlab"; O.Sigmund, 2001
    // -> dc/dx_e = -p * x_e ^(p-1) u_e^T k_0 u_e
    TransferFunction* tf = design->GetTransferFunction(DesignElement::DENSITY, MECH, true);
    double factor = excite.GetWeightedFactor(func);
    
    if(IsTransient()){
      // this computes the complete derivative of the Newmark scheme, up to now, all objectives/constraints can be handled like that
      // as the derivative of all objectives/constraint is calculated p^T (dF - dA) u 
      // where p is solution of adjoint, dF is derivative of newmark update, dA is derivative of system matrix, u is solution of forward problem
      CalcNewmarkDerivative(excite, forward, adjoint, factor, f, g);
    }else{
      CalcU1KU2(tf, forward.Get(excite)->elem[MECH], MECH, forward.Get(excite)->elem[MECH], NULL, -factor, STANDARD, func);
    }
  }
  else
  {
    UInt timesteps = domain->GetDriver()->GetNumSteps();
    result = 0.0;
    for(unsigned int ts = 0; ts < timesteps; ++ts){ // this formulation works for transient as well as static cases, integral over time
      // compliance is easier computed using f^T u on nodes with force
      // to avoid any work for assembling force again, we simply calculate solution times rhs from the system
      Vector<double>& u = forward.Get(excite, NULL, ts)->GetRealVector(Solution::RAW_VECTOR);
      Vector<double>& rhs = forward.Get(excite, NULL, ts)->GetRealVector(Solution::RHS_VECTOR);
      double sp = 0.0;
      u.Inner(rhs, sp);
      result += sp * excite.GetFactor(func) * GetStepWeight(ts);
      LOG_DBG(em) << "CalcCompliance(): result=" << result << " sp=" << sp << " u=" << u.ToString() << " func=" << func->ToString();
    }
  }
  return result;
}

template <class T>
double ErsatzMaterial::CalcOutputObjective(Excitation& excite, Objective* cost)
{
  // The output is <l,u> where l is -1* rhs of the adjoint pde
  // Here the rhs of the adoint pde is not -1 -> hence we do -1*<l,u>
  //
  // We always take l from rhs and don't store it explicitly.
  // Note that one has to use the algsys RHS! The PDE RHS is still from the
  // forward simulation!
  Vector<T>& u = dynamic_cast<Vector<T> & >(*(forward.Get(excite, NULL)->GetVector(Solution::RAW_VECTOR)));
  Vector<T>& l = dynamic_cast<Vector<T> & >(*(adjoint.Get(excite, cost)->GetVector(Solution::SEL_VECTOR)));
  assert(u.GetSize() == l.GetSize());

  LOG_DBG2(em) << "OUPTUT: adjoint sel (l): " << l.ToString(1);
  LOG_DBG2(em) << "OUPTUT: forward sol (u): " << u.ToString(0);

  double result = 0.0;

  switch(cost->GetType())
  {
    case Objective::OUTPUT:
    {
      // this is <l, u> which is for complex not really defined as it might be non-real
      T inner = u.Inner(l);
      result = ((complex<double>) inner).real();
      result *= excite.GetFactor(cost);
      LOG_DBG2(em) << "output <l,u>: " << inner << " * " << excite.GetFactor(cost) << " -> " << result;
      break;
    }

    case Objective::ABS_DYN_OUTPUT_SQUARED:
    {
      // |<u,l>|^2 = Re{<u,l>}^2 + Im{<u,l>}^2 
      T ul = u.Inner(l);
      double re_ul = ((complex<double>) ul).real();
      double im_ul = ((complex<double>) ul).imag();
      result = re_ul * re_ul + im_ul * im_ul;
      result *= excite.GetFactor(cost);
      LOG_DBG2(em) << "output |<u,l>|^2 = Re{<u,l>}^2 + Im{<u,l>}^2: " << re_ul << "^2  + " << im_ul << "^2 -> " << result;
      break;
    }


    case Objective::DYNAMIC_OUTPUT:
    case Objective::CONJUGATE_COMPLIANCE:
    {
      // this is <u,L conj(u)> and only defined for the harmonic case!
      if(!harmonic) throw Exception("'" + cost->type.ToString(cost->GetType()) + "' is only defined for harmonic!");

      // we loop over the vectors and do the scalar product by hand as we have
      // no diagonal matrix version of l
      for(unsigned int i = 0; i < u.GetSize(); i++ )
      {
        if(l[i] == 0.0) continue; // we skip this so we can make output for the real cases

        complex<double> u_val = (complex<double>) u[i];
        // make sure we have no penalization stuff!
        assert(std::abs(u_val) < 1e15);

        double sp = std::real(std::conj(u_val) * l[i] * u_val);
        LOG_DBG2(em) << "CalcObjective: " << std::conj(u_val) << " * " << l[i] << " * " << u_val << " -> " << sp;
        result += sp;
      }
      LOG_DBG2(em) << "output <u,L u*>: " << result << " * " << excite.GetFactor(cost) << " -> " << result * excite.GetFactor(cost);
      result *= excite.GetFactor(cost);
      break;
    }

    default: EXCEPTION("Not handled");
  }

  return result;
}

void ErsatzMaterial::SetAdjointRhs(AdjointParameters* adjointParams){
  int ts = domain->GetDriver()->GetActStep("mech") - 1; // drivers count timesteps starting with 1
  unsigned int nts = domain->GetDriver()->GetNumSteps();

  Excitation& excite = *(adjointParams->GetExcitation());
  
  switch(adjointParams->GetFunction()->GetType()){
  case Objective::TRACKING:
    SetTrackingAdjointRhs(excite, ts);
    break;
  case Objective::COMPLIANCE: // adjoint has the original rhs (scaled with weight per timestep), but time walks backwards
  {
    assemble_->AssembleLinRHS(NULL);
    Vector<Double> rhs;
    assemble_->GetAlgSys()->GetRHSVal(rhs);
    double w = GetStepWeight(ts);
    for(unsigned int i = 0; i < rhs.GetSize(); ++i){
      rhs[i] *= w;
    }
    assemble_->GetAlgSys()->InitRHS(rhs);
  }
    break;
  
  default:
    EXCEPTION("adjoint RHS not known for Objective " + adjointParams->GetFunction()->ToString());
  }
  
  // if the first step is static, we have to readjust the right hand side for ts == 0. Note that the transientDriver/solveStep does no rhs update any more
  if(IsFirstTransientStepStatic() && ts == 0){
    double dt = dynamic_cast<TransientDriver*>(domain->GetDriver())->GetDeltaT();
    double gamma = pde->getTimeStepping()->GetNewmarkGamma();
    double beta = pde->getTimeStepping()->GetNewmarkBeta();
    Vector<Double> coeffMass;
    Vector<Double> coeffDamping;
    // look up, whether the damping matrix exists
    std::set<FEMatrixType> matTypes;    
    assemble_->GetAlgSys()->GetFEMatrixTypes(matTypes);
    bool damping = ( matTypes.find(CoupledField::DAMPING) != matTypes.end() );
    double u = 1.0; double upp = 1.0 / (beta*dt*dt); double up = upp * gamma * dt;
    upp = 0.0;
    up = 0.0;
    for(unsigned int t = 1; t < nts; ++t){
      Vector<Double>& p_vec = adjoint.Get(excite, adjointParams->GetFunction(), t)->GetRealVector(Solution::RAW_VECTOR);
      double ut = u +  (up + 0.5*upp*(1.0-2.0*beta)*dt ) *dt;
      double upt = up + (1.0 - gamma) * dt * upp;
      // now we have the factor ut / (beta * dt * dt) for the mass matrix
      coeffMass = p_vec * (ut / (beta * dt * dt));
      assemble_->GetAlgSys()->UpdateRHS(CoupledField::MASS, coeffMass);
      if(damping){
        coeffDamping = p_vec * (ut * gamma / (beta * dt) - upt);
        assemble_->GetAlgSys()->UpdateRHS(CoupledField::DAMPING, coeffDamping);
      }
      u = 0.0;
      upp = (u - ut) / (beta * dt * dt);
      up = (upt + upp * gamma * dt);      
    }
  }
  
}

Vector<double> ErsatzMaterial::CalcVonMisesStressVector(Excitation& excite, Function* f, bool adjoint_rhs, bool grad_contrib)
{
  // Kocvara and Stingl; 2007
  // stress = rho^p E_0 B(ip) * u
  // von Mises stress per element: sum_ip stress^T * M * stress
  // rhs w.r.p to one element: 2* sum_ip strees^T * M * rho^p E_0 B(ip)
  //            if i=j:      + 2 * stress^T * M * (rho^p)' E_0 B(ip) u

  // we follow the qp-relaxation: Bruggi 2008, actually going back to Duysinx and Bendsoe 1998
  // The idea of the qp-relaxion is, that the stress is calculated based on the standard simp penalization (p).
  // To have vanishing constraints for vanishing stress we divide by rho^q (stress).
  // In this implementation we make a new penalization and apply it to the stress, this is easier as we
  // need no chain rule for the derivative. Clearly works only if the mech and stress transfer functions are simp.

  LOG_DBG2(em) << "CVMSV: excite=" << excite.index << " adjont=" << adjoint_rhs << " grad_contrib=" << grad_contrib;
  assert((adjoint_rhs && !grad_contrib) || !adjoint_rhs);
  // adjoint_rhs:   result is raw vector size with summed up 2 * stress^T * M * (rho^p * E_0 * B)
  // !grad_contrib: result is design size with stress^T * M stress (element von Mises stress)
  // grad_contrib:  result is design size with 2 * stress^T * M * ( (rho^p)' * E_0 * B * u)

  // see comment in header for locality of stress!

  assert(design->design.GetSize() == 1);

  StdVector<SingleVector*>& all_u_elem = forward.Get(excite)->elem[MECH];

  // the result vector is either per design element or the rhs for the adjoint problem
  Vector<double> result(adjoint_rhs ? forward.Get(excite)->GetVector(Solution::RAW_VECTOR)->GetSize() : design->data.GetSize());
  result.Init();

  Matrix<double> M = dynamic_cast<MechPDE*>(ToPDE(MECH))->GetVonMisesMatrix(dim);
  Matrix<double> B;
  Matrix<double> E; // the material matrix with applied pseudo density
  Matrix<double> coords;

  assert(design->design.GetSize() == 1); // easy to extend
  assert(design->regions.GetSize() == 1); // easy to extend

  // the are special transfer functions for the stress (q) as relaxation and the fem simulation (p). See above!
  TransferFunction* tf_p = design->GetTransferFunction(DesignElement::DENSITY, MECH);
  TransferFunction* tf_q = design->GetTransferFunction(DesignElement::DENSITY, STRESS);

  if(tf_q->GetType() != TransferFunction::SIMP_TYPE || tf_p->GetType() != TransferFunction::SIMP_TYPE)
    throw Exception("Stress constraints are currently only implemented for simp type transfer functions in 'mech' and 'stress'.");
  TransferFunction qp(NO_APP, TransferFunction::SIMP_TYPE, tf_p->GetParam() - tf_q->GetParam(), tf_p->GetDesign());
  LOG_DBG2(em) << "CVMSV: p=" << tf_p->GetParam() << " q=" << tf_q->GetParam() << " qp=" << qp.GetParam();

  linElastInt* form = dynamic_cast<linElastInt*>(GetForm(design->GetRegionId(), pde, pde, "linElastInt"));

  Vector<double> stress(dim == 2 ? 3 : 6); // elem stress
  Vector<double> strain;
  Vector<double> M_stress;
  Matrix<double> E_B;
  Matrix<double> M_E_B;
  Matrix<double> stress_transp(1, stress.GetSize());
  Matrix<double> rhs_transp;
  Vector<double> intPoint;
  shared_ptr<EqnMap> eqnMap = ToPDE(MECH)->GetEqnMap();

  ElemList elemList(grid);

  // in the adjoint case we need alpha which is the gradient of the globalization function
  Vector<double> alpha;
  if(adjoint_rhs) alpha = CalcVonMisesStressGlobalizationFactor(excite, f);

  // output of penalized von mises stresses? Only used in !adjoint_rhs and !grad_contrib
  int res_idx = design->GetSpecialResultIndex(DesignElement::DEFAULT, DesignElement::PENALIZED_STRESS, DesignElement::NONE, DesignElement::PLAIN, excite.label);

  for(unsigned int e = 0, en = design->data.GetSize(); e < en; e++)
  {
    DesignElement* de = &design->data[e];

    grid->GetElemNodesCoord(coords, de->elem->connect, false); // geometric non-lin

    Vector<double>& u_elem = dynamic_cast<Vector<double>& >(*(all_u_elem[e]));

    // apply our own physical densities!
    form->calcDMat(E, NULL, DesignElement::DENSITY, qp.Transform(de, DesignElement::SMART));

    elemList.SetElement(de->elem);
    EntityIterator it = elemList.GetIterator();
    form->ExtractElemInfo(it);
    de->elem->ptElem->GetCoordMidPoint(intPoint);
    form->SetIntPoint(intPoint);

    // same numerical results as with  for(int ip = 1, nip = de->elem->ptElem->GetNumIntPoints(); ip <= nip; ip++)
    form->calcBMatOnly(B, 1, de->elem->ptElem, coords);

    // calc stress
    strain = B * u_elem; // strain
    stress = E * strain; // pure stress

    if(adjoint_rhs)
    {
      // this is actually 2* von Mises stress with one missing u_elem: 2 * stress^T * M * rho^p * E_0 * B, note the rho^p is already in E
      E_B = E * B;
      M_E_B = M * E_B;

      // we have to transpose the stress (3*1 or 6*1) manually :(
      for(unsigned int si = 0; si < stress.GetSize(); si++)
        stress_transp[0][si] = stress[si];

      rhs_transp = stress_transp * M_E_B;
      rhs_transp *= -2.0;

      // there is a factor from the globalization function, which is the gradient of the glob function(func_val)
      rhs_transp *= alpha[e];


      assert(rhs_transp.GetNumCols() == u_elem.GetSize());
      assert(rhs_transp.GetNumRows() == 1);

      LOG_DBG2(em) << "CVMSV de=" << de->ToString() << " ar=" << adjoint_rhs  << " alpha=" << alpha[e]
                   << " rhs_transp=" << rhs_transp.ToString();

      // sum it up to the global rhs vector
      assert(de->elem->connect.GetSize() * dim == rhs_transp.GetNumCols());
      for(unsigned int n = 0; n < de->elem->connect.GetSize(); n++)
      {
        unsigned int node =  de->elem->connect[n];

        for(unsigned int dof = 1; dof <= dim; dof++)
        {
          // fuck 1-based in GetNodeEqn() !!
          int eqn_nr = std::abs(eqnMap->GetNodeEqn(node,dof)); // map periodic bc to the master equation
          assert(eqn_nr >= 0);
          int eqn_idx = eqn_nr -1; // fuck 1-based!!
          // don't set the homogeneous dirichlet boundary conditions :)
          if(eqn_idx >= 0)
          {
            result[eqn_idx] += rhs_transp[0][dim * n + (dof-1)];

            LOG_DBG3(em) << "CVMSV de=" << de->ToString() << " ar=" << adjoint_rhs  << " n=" << n << " node=" << node
                << " dof=" << dof << " eqn_idx=" << eqn_idx << " idx=" << (dim * n + (dof-1)) << " val="
                << rhs_transp[0][dim * n + (dof-1)] << " -> " << result[eqn_idx];
          }
        }
      }

    }
    else
    {
      // normal and "corrected" von Mises stress element results
      double correct = 1.0;
      M_stress = M * stress;
      result[e] = stress.Inner(M_stress);

      if(!grad_contrib)
      {
        // output von mises stress? Note, that this is excitation specific!
        if(res_idx != -1)
        {
          de->specialResult[res_idx] = result[e];
          LOG_DBG3(em) << "CVMSV:sr de=" << de->ToString() << " si=" << res_idx << " v=" << result[e];
        }
      }
      else
      {
        // additional factor for grad_contrib. We replace one rho^p by (rho^p)'
        correct = 2.0 * qp.Derivative(de, DesignElement::SMART) / qp.Transform(de, DesignElement::SMART);
      }

      LOG_DBG2(em) << "CVMSV de=" << de->ToString() << " gv= " << grad_contrib << " rho=" << de->GetDesign(DesignElement::SMART) << " sMs=" << result[e] << " deriv="
                   << qp.Derivative(de, DesignElement::SMART) << " trans=" <<  qp.Transform(de, DesignElement::SMART) << " correct=" << correct << " -> " << (correct * result[e]);

      result[e] *= correct;
    }
  }

  return result;
}


Vector<double> ErsatzMaterial::CalcVonMisesStressGlobalizationFactor(Excitation& excite, Function* f)
{
  assert(f->GetType() == Function::STRESS);

  const Function::Local* local = f->GetLocal();
  assert(local != NULL);

  const Vector<double> stress = CalcVonMisesStressVector(excite, f, false, false);
  Vector<double> result(stress.GetSize());

  // this is a complicated way to calc element wise power*max(0,stress)^(power-1) but that way
  // we are open for further globalization functions and it is not that expensive

  const StdVector<Function::Local::Identifier>& vem = local->virtual_elem_map;

  for(unsigned int i = 0; i < vem.GetSize(); i++)
  {
    const Function::Local::Identifier& id = vem[i];
    assert(id.neighbor.GetSize() == 0);
    double gfv = id.EvalFunction(local, true, &stress);
    int idx = design->Find(id.element);
    result[idx] = gfv;
  }

  LOG_DBG2(em) << "CVMSGF ex=" << excite.index << " stress=" << stress.ToString();
  LOG_DBG2(em) << "CVMSGF ex=" << excite.index << " alpha=" << result.ToString();

  return result;
}


double ErsatzMaterial::CalcEnergyFlux(Excitation& excite, Objective* f)
{
  // calc 1/2 Re{j omega psi^R Q psi^*} where Q is the grad_n matrix B applied at some points L

    // this is the global element solution, indexed by equation number
  Vector<complex<double> > u_glob = forward.Get(excite)->GetComplexVector(Solution::RAW_VECTOR);
  // here we store Q*u^* as we have to determine the nodal entries by count
  Vector<complex<double> > q_u_glob(u_glob.GetSize());

  SetEnergyFluxVector(f, u_glob, false, q_u_glob);

  // calculate the product Re(j*omega*u*(Q*u^*)) which is -Im(u*(Q*u^*))
  // the (Q*u^*) term is normalized
  double result = 0;
  for(unsigned int i = 0, in = q_u_glob.GetSize(); i < in; i++)
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
  SinglePDE* mypde = ToPDE(ACOUSTIC);
  shared_ptr<ResultInfo> res_info = mypde->GetResultInfo(ACOU_POTENTIAL);
  shared_ptr<EqnMap> eqnMap = mypde->GetEqnMap();

  // get the surface region description
  PtrParamNode srpn = f->pn->Get("output/energyFlux/surfRegion");
  // surface element list
  shared_ptr<EntityList> sel = grid->GetEntityList(EntityList::SURF_ELEM_LIST, srpn->Get("name")->As<string>(), EntityList::REGION);
  // neighbor region is mandatory to define the volume elements
  RegionIdType vol_neigh = grid->GetRegion().Parse(srpn->Get("neighborRegion")->As<string>());

  // reset output
  q_u_glob.Init(complex<double>(0.0, 0.0));
  // here we count the entries to q_u_glob to normalize in the end
  Vector<int> count(u_glob.GetSize());
  count.Init(0);

  // an element solution vector -> we need a 1 dof solution up to now!
  Vector<complex<double> > elem_sol;

  // node numbers common on a surface element and the matching volume element
  StdVector<unsigned int> common_nodes;

  // this contains an element B (grad_n) matrix to be applied with the solution.
  Matrix<complex<double> > q_mat;

  // traverse our surface elements
  EntityIterator it = sel->GetIterator();
  for(it.Begin(); !it.IsEnd(); it++)
  {
    const SurfElem* se = it.GetSurfElem();
    const Elem*     vol = se->ptVolElem1 != NULL && se->ptVolElem1->regionId == vol_neigh ? se->ptVolElem1 : se->ptVolElem2;
    assert(vol->regionId == vol_neigh);

    // determine grad_n u^* on the se nodes coinciding with the vol nodes
    FindCommonNodes(se, vol, common_nodes);

    // get the element solution
    ElemList se_it(grid); //single element iterator
    se_it.SetElement(vol);
    mypde->GetSolVecOfElement(elem_sol, se_it.GetIterator(), res_info);

    // determine selected element grad_n matrix (includes selection defined by surface elements)
    // the matrix is squared but for non-common nodes entries it is zero
    mypde->CalcElemGradMatrix(vol, common_nodes, 3, vol->ptElem->GetAnsatzFct(), q_mat); // z-direction!

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
      int eqn_idx = eqn_nr -1; // FUCK!!!
      if(!close(sum, complex<double>(0.0, 0.0)))
      {
        q_u_glob[eqn_idx] += sum;
        count[eqn_idx]++;
      }
      LOG_DBG3(em) << "SEFV: vol=" << vol->elemNum << " node=" << node << " eqn_nr=" << eqn_nr << " count=" << count[eqn_idx] << " sum=" << sum;
    }
  }

  LOG_DBG2(em) << "SEFV: q_u_glob=" << q_u_glob.ToString(1);
  LOG_DBG2(em) << "SEFV: count="    << count.ToString(1);

  // normalize Q*u^*
  for(unsigned int i = 0, in = q_u_glob.GetSize(); i < in; i++)
    if(count[i] != 0)
      q_u_glob[i] /= (double) count[i];
}


void ErsatzMaterial::FindCommonNodes(const SurfElem* se, const Elem* vol, StdVector<unsigned int>& common_nodes) const
{
  const StdVector<unsigned int>& se_nodes = se->connect;
  const StdVector<unsigned int>& vol_nodes = vol->connect;

  // for higher order elements the center node is not on the vol, for lin common_nodes = se nodes
  common_nodes.Reserve(se_nodes.GetSize());
  common_nodes.Clear();

  for(unsigned int i = 0; i < vol_nodes.GetSize(); i++)
  {
    unsigned int n = vol_nodes[i];
    if(se_nodes.Contains(n))
      common_nodes.Push_back(n);
  }

  LOG_DBG3(em) << "FCN se=" << se->elemNum << " vol=" << vol->elemNum << " common=" << common_nodes.ToString();
}

void ErsatzMaterial::SetTrackingAdjointRhs(Excitation& excite, int ts){
  // this is for the static and for the transient case.
  Vector<double>& u = forward.Get(excite, NULL, ts)->GetRealVector(Solution::RAW_VECTOR); // Tracking is only implemented for non-harmonic
  LOG_DBG3(em) << "SolveTrackingProblem: displacement vector: (" << u.GetSize() << ") " << u.ToString();

  shared_ptr<EqnMap> eqnMap = pde->GetEqnMap();
  MathParser * parser = domain->GetMathParser();
  unsigned int mathParserHandle = parser->GetNewHandle();
  CoordSystem * coosy = domain->GetCoordSystem();

  Vector<Double> rhs;
  assemble_->GetAlgSys()->GetRHSVal(rhs);

  // set rhs to 0
  rhs.Init();
  
  double w = GetStepWeight(ts);

  // set rhs for tracking nodes
  Vector<double> GlobCoord;
  for(unsigned int l=0; l < excite.trackings.GetSize(); l++){
    TrackingBc const & actTrack = *(excite.trackings[l]);
    EntityIterator it = actTrack.entities->GetIterator();
    const int dof = actTrack.dof;
    for(it.Begin(); !it.IsEnd(); it++){
      grid->GetNodeCoordinate(GlobCoord, it.GetNode());
      parser->SetCoordinates(mathParserHandle, *coosy, GlobCoord);
      parser->SetExpr(mathParserHandle, actTrack.value);
      const double uref = parser->Eval(mathParserHandle);
      const int eqnr = eqnMap->GetNodeEqn(it.GetNode(), dof) - 1; // equation numbers are 1 based but vector u is 0 based
      if(eqnr>=0){
        double v = u[eqnr] - uref;
        rhs[eqnr] = v * w;
        LOG_DBG2(em) << "SolveTrackingProblem: tracking setting RHS equation " << eqnr << " (Node: " << it.GetNode() << ", dof: " << dof << ") to " << rhs[eqnr];
      }
    }
  }
  assemble_->GetAlgSys()->InitRHS(rhs);
}

double ErsatzMaterial::CalcTracking(Excitation& excite, Objective* c, Condition* g, bool derivative)
{
  Function* f = Function::Cast(c, g);

  UInt timesteps = domain->GetDriver()->GetNumSteps();
  if(derivative){
    // calculate the tracking functional gradient, which is z^T k_i u,
    // where Kz = ut
    // where ut = M^T M (u-u0) = I_\Gamma (u-u0)

    // calculate gradient z^T k_i u
    // note that in multiple excitations case (this is always now), we do sum this up
    TransferFunction* tf = design->GetTransferFunction(DesignElement::DENSITY, MECH, false);
    double factor = excite.GetWeightedFactor(c);
    
    if(IsTransient()){
      // this computes the complete derivative of the Newmark scheme, up to now, all objectives/constraints can be handled like that
      // as the derivative of all objectives/constraint is calculated p^T (dF - dA) u 
      // where p is solution of adjoint, dF is derivative of newmark update, dA is derivative of system matrix, u is solution of forward problem
      CalcNewmarkDerivative(excite, forward, adjoint, factor, c, g);
    }else{
      CalcU1KU2(tf, adjoint.Get(excite, f)->elem[MECH], MECH, forward.Get(excite)->elem[MECH], NULL, -factor, STANDARD, f);
    }
    return 0.0;
  }else{
    // prepare for calculation
    shared_ptr<EqnMap> eqnMap = pde->GetEqnMap();
    MathParser * parser = domain->GetMathParser();
    unsigned int mathParserHandle = parser->GetNewHandle();
    CoordSystem * coosy = domain->GetCoordSystem();

    // the tracking functional is ut^T ut (ut as above)
    double val = 0.0;
    Vector<double> GlobCoord;
    
    double dt = 0.0;
    if(IsTransient()){
      dt = dynamic_cast<TransientDriver*>(domain->GetDriver())->GetDeltaT();
    }
    parser->SetValue(MathParser::GLOB_HANDLER, "dt", dt);
    
    for(unsigned int ts = 0; ts < timesteps; ++ts){ // this formulation works for transient as well as static cases, integral over time
      Vector<double>& u = forward.Get(excite, NULL, ts)->GetRealVector(Solution::RAW_VECTOR); // Tracking is only implemented for non-harmonic
      LOG_DBG3(em) << "CalcTracking: displacement vector: (" << u.GetSize() << ") " << u.ToString();
      parser->SetValue(MathParser::GLOB_HANDLER, "t", dt*(ts+1));
      parser->SetValue(MathParser::GLOB_HANDLER, "step", ts+1);

      double v = 0.0;
      for(unsigned int l=0; l < excite.trackings.GetSize(); l++){
        TrackingBc const & actTrack = *(excite.trackings[l]);
        EntityIterator it = actTrack.entities->GetIterator();
        const int dof = actTrack.dof;
        for(it.Begin(); !it.IsEnd(); it++){
          grid->GetNodeCoordinate(GlobCoord, it.GetNode(), true); // get Updated Coordinate
          parser->SetCoordinates(mathParserHandle, *coosy, GlobCoord);
          parser->SetExpr(mathParserHandle, actTrack.value);
          const double uref = parser->Eval(mathParserHandle);
          const int eqnr = eqnMap->GetNodeEqn(it.GetNode(), dof) - 1; // equation numbers are 1 based, the vector u is 0 based
          const double uact = eqnr>=0 ? u[eqnr] : 0.0;
          v += (uact - uref) * (uact - uref);
        }
      }
      val += v * GetStepWeight(ts);
    }
    return 0.5 * val;
  }
}

double ErsatzMaterial::CalcPoissonsRatioAndYoungsModulus(Objective* cost, Condition* g, bool derivative)
{
  assert(cost != NULL || g != NULL); // at least one of them must be given!
  
  bool poisson(false);
  if(cost != NULL)
  {
    assert(cost->GetType() == Objective::POISSONS_RATIO || cost->GetType() == Objective::YOUNGS_MODULUS);
    poisson = cost->GetType() == Objective::POISSONS_RATIO;
  }
  else
  {
    assert(g->GetType() == Condition::POISSONS_RATIO || g->GetType() == Condition::YOUNGS_MODULUS);
    poisson = g->GetType() == Condition::POISSONS_RATIO;
  }

  Matrix<double> hom_tensor = CalcHomogenizedTensor();

  double result = 0.0;

  if(derivative)
  {
    StdVector<double> dE11;
    CalcHomogenizedTensorEntry(make_pair(1,1), true, dE11);
    StdVector<double> dE12;
    CalcHomogenizedTensorEntry(make_pair(1,2), true, dE12);

    const double E11 = hom_tensor[0][0];
    const double E12 = hom_tensor[0][1];
    double grad(0.0);
    const double denom = (E11 + E12) * (E11 + E12);
    
    for(unsigned int e = 0, ne = design->GetNumberOfElements(); e < ne; e++)
    {      
      if(poisson) 
        grad = (dE12[e] * E11 - E12 * dE11[e]);
      else
      {
        grad  = (E11 * E11 + 2.0 * E11 * E12 + 3.0 * E12 * E12) * dE11[e];
        grad -= (4.0 * E11 * E12 + 2.0 * E12 * E12) * dE12[e];
      }
      
      grad /= denom;

      design->data[e].AddGradient(cost, g, grad);
    }
  }
  else
  {
    if(poisson)
      result = MechanicMaterial::CalcIsotropicPoissonsRatio(hom_tensor);
    else
      result = MechanicMaterial::CalcIsotropicYoungsModulus(hom_tensor);
  }

  return result;
}


void ErsatzMaterial::SetTestStrainMatrix(Matrix<double> &matrix, const Vector<double> &vec)
{
  assert(matrix.GetNumCols() == dim);
  assert(matrix.GetNumRows() == dim);
  matrix[0][0] = vec[0];
  matrix[1][1] = vec[1];
  matrix[0][1] = vec[5]; // Voigt notation!
  matrix[1][0] = 0.0; // because of symmetry we need factor 0.5

  if(dim == 3)
  {
    matrix[2][2] = vec[2];
    matrix[1][2] = vec[3];
    matrix[2][1] = 0.0; // symmetry again

    matrix[0][2] = vec[4];
    matrix[2][0] = 0.0; // symmetry again
  }
}

Matrix<double> ErsatzMaterial::CalcHomogenizedTensor()
{
  const double cube_vol(grid->CalcVolumeSpannedByNamedNodes());

  unsigned int ex_size = me->excitations.GetSize();
  assert((dim == 2 && ex_size == 3) || (dim == 3 && ex_size == 6));

  Matrix<double> test_strain_matrix_ij(dim, dim);
  Matrix<double> test_strain_matrix_kl(dim, dim);

  Matrix<double> result(ex_size, ex_size);
  result.Init();

  for(unsigned int ij = 0; ij < ex_size; ++ij)
  {
    SetTestStrainMatrix(test_strain_matrix_ij, me->excitations[ij].test_strain);
    StdVector<SingleVector*> &u1 = forward.Get(ij)->elem[MECH]; // equal to \chi^{ij}

    for(unsigned int kl = 0; kl < ex_size; ++kl)
    {
      if(ij > kl) // already computed this entry!
      {
        // by construction, the resulting matrix must be symmetric
        // so we can skip the calculation and do a simple assignment instead
        result[ij][kl] = result[kl][ij];
        continue;
      }

      SetTestStrainMatrix(test_strain_matrix_kl, me->excitations[kl].test_strain);
      StdVector<SingleVector*> &u2 = forward.Get(kl)->elem[MECH]; // equal to \chi^{kl}

      // loop over elements. In the gradient case not summed up
      for(int e = 0, ne = design->GetNumberOfElements(); e < ne; ++e)
      {
        DesignElement* de = &design->data[e];
        Vector<double>& u1_vec = dynamic_cast<Vector<double>& >(*u1[e]);
        Vector<double>& u2_vec = dynamic_cast<Vector<double>& >(*u2[e]);

        // prepare for calculation
        double p = CalcHomogenizedElementProduct(this, de, false, u1_vec, u2_vec, test_strain_matrix_ij, test_strain_matrix_kl);
        result[ij][kl] += p/cube_vol; // normalize for volume
      }
    } // end of kl loop
  } // end of ij loop

  // save e.g. for CommitIteration()
  homogenizedTensor.Assign(result, 1.0);

  return result;
}

void ErsatzMaterial::CalcHomogenizedTrackingGradient(const Matrix<double>& target, const Matrix<double>& hom, Objective* f, Condition* g)
{
  const double cube_vol(grid->CalcVolumeSpannedByNamedNodes());

  // TODO Would be E^* - E^H if expression templates would work
  Matrix<double> diff_tensor;
  diff_tensor = target - hom;

  const unsigned int ex_size(me->excitations.GetSize());
  assert((dim == 2 && ex_size == 3) || (dim == 3 && ex_size == 6));

  Matrix<double> test_strain_matrix_ij(dim, dim);
  Matrix<double> test_strain_matrix_kl(dim, dim);

  // our derivative tensor
  Matrix<double> hom_tensor_deriv(ex_size, ex_size);
  hom_tensor_deriv.Init(); // we set and do not add - hence one init is enough

  // loop over elements.
  for(int e = 0, ne = design->GetNumberOfElements(); e < ne; ++e)
  {
    DesignElement* de = &design->data[e];

    for(unsigned int ij = 0; ij < ex_size; ++ij)
    {
      SetTestStrainMatrix(test_strain_matrix_ij, me->excitations[ij].test_strain);
      StdVector<SingleVector*> &u1 = forward.Get(ij)->elem[MECH]; // equal to \chi^{ij}
      Vector<double>& u1_vec = dynamic_cast<Vector<double>& >(*u1[e]);

      for(unsigned int kl = 0; kl < ex_size; ++kl)
      {
        if(ij > kl) // already computed this entry!
        {
          // by construction, the resulting matrix must be symmetric
          // so we can skip the calculation and do a simple assignment instead
          hom_tensor_deriv[ij][kl] = hom_tensor_deriv[kl][ij];
          continue;
        }

        SetTestStrainMatrix(test_strain_matrix_kl, me->excitations[kl].test_strain);
        StdVector<SingleVector*> &u2 = forward.Get(kl)->elem[MECH]; // equal to \chi^{kl}
        Vector<double>& u2_vec = dynamic_cast<Vector<double>& >(*u2[e]);

        // prepare for calculation
        double p = CalcHomogenizedElementProduct(this, de, true, u1_vec, u2_vec, test_strain_matrix_ij, test_strain_matrix_kl);
        hom_tensor_deriv[ij][kl] = p/cube_vol; // normalize for volume
      } // end of kl loop
    } // end of ij loop

    // hom_tensor_deriv is completely set.
    // (E^* - E^H) * - d(E^H)/d(rho_e) -> therefore the minus !
    double grad = -1.0 * diff_tensor.ScalarProduct(hom_tensor_deriv);

    de->AddGradient(f, g, grad);
  } // element loop
}


double ErsatzMaterial::CalcHomogenizedTensorConstraint(Condition* g, bool derivative)
{
  // me make use of the multi-purpose CalcHomogenizedTensorEntry()
  StdVector<double> grad;

  double result = 0.0;

  // we have the form Exx or Exx-Eyy or Exx-Eyy-Ezz-Ezz for isotropy constraints

  for(unsigned int i = 0; i < g->coords.GetSize(); i++)
  {
    double t = CalcHomogenizedTensorEntry(g->coords[i], derivative, grad);

    if(derivative)
    {
      for(int e = 0, ne = design->GetNumberOfElements(); e < ne; ++e)
        design->data[e].AddGradient(NULL, g, (i == 0 ? 1.0 : -1.0) * grad[e]);
    }
    else
    {
      result = (i == 0 ? t : result-t);

      homogenizedTensor[g->coords[i].first - 1][g->coords[i].second - 1] = t;
      // all tensors are symmetric. Makes reading easier!
      homogenizedTensor[g->coords[i].second - 1][g->coords[i].first - 1] = t;

      //LOG_DBG(em) << "CHTC: g=" << g->ToString() << " coord=" << i << " ["
      //    << g->coords[i].first << "-1][" << g->coords[i].second << "-1] = " << t;
    }
  }

  return result;
}

double ErsatzMaterial::CalcHomogenizedTensorEntry(const std::pair<int, int> entry, bool derivative, StdVector<double>& grad_out)
{
  const double cube_vol(grid->CalcVolumeSpannedByNamedNodes());

  assert((dim == 2 && me->excitations.GetSize() == 3) || (dim == 3 && me->excitations.GetSize() == 6));

  Matrix<double> test_strain_matrix_ij(dim, dim);
  Matrix<double> test_strain_matrix_kl(dim, dim);

  const unsigned int ij = entry.first - 1;
  const unsigned int kl = entry.second - 1;

  SetTestStrainMatrix(test_strain_matrix_ij, me->excitations[ij].test_strain);
  StdVector<SingleVector*> &u1 = forward.Get(ij)->elem[MECH]; // equal to \chi^{ij}

  SetTestStrainMatrix(test_strain_matrix_kl, me->excitations[kl].test_strain);
  StdVector<SingleVector*> &u2 = forward.Get(kl)->elem[MECH]; // equal to \chi^{kl}

  double result = 0.0;

  if(derivative) grad_out.Resize(design->GetNumberOfElements());

  // loop over elements. In the gradient case not summed up
  for(int e = 0, ne = design->GetNumberOfElements(); e < ne; ++e)
  {
    DesignElement* de = &design->data[e];
    Vector<double>& u1_vec = dynamic_cast<Vector<double>& >(*u1[e]);
    Vector<double>& u2_vec = dynamic_cast<Vector<double>& >(*u2[e]);

    // prepare for calculation
    double p = CalcHomogenizedElementProduct(this, de, derivative, u1_vec, u2_vec, test_strain_matrix_ij, test_strain_matrix_kl);
    result += p/cube_vol; // normalize for volume

    if(derivative)
    {
      grad_out[e] = result;
      result = 0.0; // reset such that we do not sum up for the next case!
    }
  }

  return result; // in the non-derivative case this is the sum.
}


double ErsatzMaterial::CalcHomogenizedElementProduct(ErsatzMaterial* obj, DesignElement* de, bool derivative,
    Vector<double>& u1_vec, Vector<double>& u2_vec,
    Matrix<double>& test_strain_matrix_ij, Matrix<double>& test_strain_matrix_kl)
{
  assert(u1_vec.GetSize() == u2_vec.GetSize());

  const unsigned int dim_ = obj->dim;

  LOG_DBG3(em) << "CHEP: de=" << de->ToString() << " u1_vec=" << u1_vec.ToString() << " u2_vec=" << u2_vec.ToString();

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
  assert(u1_tmp.GetNumRows() * u1_tmp.GetNumCols() == u1_vec.GetSize());

  // u1_tmp, u2_tmp have to be transformed into vectors
  // 2D: from 2x4 to 8x1 on quad elems, 2x3 to 6x1 on triangles
  Vector<double> u1_0(u1_vec.GetSize());
  Vector<double> u2_0(u2_vec.GetSize()); // u1 and u2 have the same size
  for(unsigned int out = 0, cols = u1_tmp.GetNumCols(); out < cols; ++out)
  {
    for(unsigned int in = 0; in < dim_; ++in)
    {
      u1_0[out*dim_ + in] = u1_tmp[in][out]; // equal to \chi^{0(ij)}
      u2_0[out*dim_ + in] = u2_tmp[in][out]; // equal to \chi^{0(kl)}
    }
  }

  u1_0 -= u1_vec;
  u2_0 -= u2_vec;

  // reuse tmp_mat as elementK-Matrix
  // Matrix<double> k_mat;
  obj->SetElementK(de, MECH, &tmp_mat, STANDARD, derivative);

  Vector<double> mat_vec;
  mat_vec = tmp_mat * u1_0;
  double result = mat_vec * u2_0;

  LOG_DBG3(em) << "CHEP de=" << de->ToString() << " result=" << result;

  return result;
}

double ErsatzMaterial::CalcGreyness(Condition* g, bool derivative)
{
  double greyness = 0.0; // element greyness
  int counter = 0; // to make it sure for different design variables!

  double lb, ub, span, org_value, value, grad, eval;
  lb = ub = value = grad = eval = span = 0.0;

  // we have to divide the gradients by their relative volume = fraction
  double fraction = g->design == DesignElement::DEFAULT ?
                    design->data.GetSize() : design->GetNumberOfElements();

  // do we want the physical value?
  TransferFunction* tf = g->IsPhysical() ? design->GetTransferFunction(g->design, MECH) : NULL;

  // go over the complete design space to set gradients of other types to 0
  for(unsigned int i = 0; i < design->data.GetSize(); i++)
  {
    DesignElement* de = &design->data[i];
    bool relevant = g->design == DesignElement::DEFAULT || g->design == de->GetType();
    if(relevant)
    {
      if(g->IsPhysical())
      {
        lb = tf->Transform(de, DesignElement::PLAIN, de->GetLowerBound());
        ub = tf->Transform(de, DesignElement::PLAIN, de->GetUpperBound());
        org_value = derivative ? tf->Derivative(de, DesignElement::SMART) : tf->Transform(de, DesignElement::SMART);
      }
      else
      {
        lb = de->GetLowerBound();
        ub = de->GetUpperBound();
        org_value = de->GetDesign(DesignElement::PLAIN);
      }

      span = ub-lb;


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

  return greyness / (double) counter;
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
  
  if(derivative)
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
    local->infeasible = 0;

    for(unsigned int i = 0; i < vem.GetSize(); i++)
    {
      Function::Local::Identifier& id = vem[i];
      double fv = id.EvalFunction(local, false, von_mises_stress);
      res += fv;
      if(fv > 0) local->infeasible++;
      LOG_DBG2(em) << "CGF: !d c=" << f->type.ToString(f->GetType()) << " i=" << i << " de="
                   << id.element->elem->elemNum << " sign=" << id.sign << " fv=" << fv  << " infeasible=" << local->infeasible << " -> " << res;
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



void ErsatzMaterial::SolveStateProblem(Excitation* ev_only_exite)
{
  // if ev_only_exite is set we use the given excitation
  // -> it shall not coincide
  assert(!(ev_only_exite != NULL && me->IsEnabled()));

  // shall we normalize afterwards?
  bool normalize = false;

  // we have to check objectives and active (non local) constraints
  StdVector<Function*> funcs = GetActiveFunctions();

  for(unsigned int e = 0; e < me->excitations.GetSize(); e++)
  {
    Excitation& excite = ev_only_exite != NULL ? *ev_only_exite : me->excitations[e];
    excite.Apply();

    // this is true for all problem types
    Optimization::SolveStateProblem(&excite);
    if(!IsTransient()){ // transient solutions are read per timestep
      StorePDESolution(forward, excite, NULL, 0, true, true, false, "forward");
    }

    for(unsigned fi = 0; fi < funcs.GetSize(); fi++)
    {
      // some functions need the selection vector for function evaluation, e.g. output
      if(funcs[fi]->NeedsSelectionVector())
        ConstructSelection(excite, funcs[fi], false); // don't change the rhs of the system but restore

      // in the harmonic case the system matrix depends on the frequency. Hence we have to
      // use the current assembly and factorization to solve the adjoint problem.
      // Note, that SolveAdointProblem*s*() must not be called, it would overwrite the adjoints with wrong results
      if(harmonic && me->excitations.GetSize() > 1)
        SolveAdjointProblem(&excite, funcs[fi]);

      // when we do multiple excitations with adjusted weights we calculate the objective here
      // to find the best weights. CalcObjective is so cheap, it is done later again by
      // the optimizer but then the weights are set
      if(me->DoAdjustWeights())
      {
        // in the first iteration we adjust the weights
        // stride = 0 is only the first time
        if((me->stride < 1 && GetCurrentIteration() == 0) ||
            ((GetCurrentIteration() % me->stride) == 0))
        {
          excite.cost = CalcFunction(excite, funcs[fi], false); // to be normalized
          normalize = true;
        }
      }
    }
  }

  // we need only to normalize when the properties have changed. This also reflects the stride
  if(normalize)
    me->NormalizeMultipleExcitations(&objectives);
}

void ErsatzMaterial::SolveAdjointProblems(Excitation* ev_only_exite)
{
  // solve all adjoints needed for gradient calculation
  assert(!(ev_only_exite != NULL && me->IsEnabled()));

  // in the harmonic multiple frequency case me must not solve for the adjoints, as
  // only in the forward problems the matrix is reassembled and the system factorized
  if(!harmonic || me->excitations.GetSize() == 1)
  {
    for(unsigned int e = 0; e < me->excitations.GetSize(); ++e)
    {
      Excitation* excite = ev_only_exite != NULL ? ev_only_exite : &me->excitations[e];

      // processs all functions and call SolveAdjointProblem() if the function requires it
      Optimization::SolveAdjointProblems(excite);
    }
  }
}

void ErsatzMaterial::StorePDESolution(Solutions& solutions, Excitation &excite, Function* f, unsigned int timestep,
                                      bool read_sol, bool read_rhs, bool save_sol, const std::string& comment)
{
  Solution& sol = *(solutions.Get(excite, f, timestep));

  SingleVector* raw = NULL; // last result = every result for multiple solutions

  // store solution element wise for gradient and raw vector for objective.
  // This is redundant as currently the solution is the global one!
  for(map<Application, SinglePDE*>::iterator it = pdes.begin(); it != pdes.end(); ++it)
  {
    stringstream ss;
    ss << "StorePDESolution: prob=" << comment << " excite=" << excite.index << " pde: " << it->first
       << " timestep=" << timestep;
    
    if(read_sol){
      sol.Read(Solution::ELEMENT_VECTORS, it->second, it->first, save_sol);
      raw = sol.Read(Solution::RAW_VECTOR, it->second, it->first, save_sol);

      LOG_DBG2(em) << ss.str() << " sol: " << raw->ToString();
    }


    if(read_rhs)
    {
      sol.Read(Solution::RHS_VECTOR, it->second, it->first, save_sol);
      LOG_DBG2(em) << ss.str() << " rhs: " << sol.GetVector(Solution::RHS_VECTOR)->ToString();
    }
  }
}

void ErsatzMaterial::TimeStepCalculated(UInt timeStep, AdjointParameters* adjParams){
  // is only called in transient case
  assert(IsTransient());
  
  // drivers start counting steps with 1
  if(adjParams == NULL){
    StorePDESolution(forward, *applied_excitation, NULL, timeStep-1, true, false, false, "forward");
  }else{
    StorePDESolution(adjoint, *applied_excitation, adjParams->GetFunction(), timeStep-1, true, false, false, "adjoint");
  }
  
}

void ErsatzMaterial::RhsCalculated(AdjointParameters* adjParams){
  if(IsTransient()){ // only in transient case this is needed
    if(adjParams == NULL){
      StorePDESolution(forward, *applied_excitation, NULL, domain->GetDriver()->GetActStep("mech")-1, false, true, false, "forward");
    }else{
      StorePDESolution(adjoint, *applied_excitation, adjParams->GetFunction(), domain->GetDriver()->GetActStep("mech")-1, false, true, false, "adjoint");
    }
  }
}


StdVector<pair<SinglePDE*, IdBcList> > ErsatzMaterial::SetHDBC()
{
  // IDBC values in the forward problem are homogeneous in the adjoint PDE. Consider elec and mech
  // store org value as idbc_list per pde
  StdVector<pair<SinglePDE*, IdBcList> > org_idbc;
  org_idbc.Reserve(pdes.size());

  for(map<Application, SinglePDE*>::iterator it = pdes.begin(); it != pdes.end(); ++it)
  {
    // we need a deep copy. Create new list for the PDE
    org_idbc.Push_back(make_pair(it->second, IdBcList()));

    // get the idbc list of the current pde
    IdBcList  idbc_list = it->second->GetIDBCList();
    for(unsigned int bc = 0; bc < idbc_list.GetSize(); bc++)
    {
      // store original value -> we have to do a deep copy!
      org_idbc.Last().second.Push_back(shared_ptr<InhomDirichletBc>(new InhomDirichletBc(*(idbc_list[bc]))));
      // make homogeneous values! -> this is only stored in the PDE, we have to call SetBCs() later!
      (*idbc_list[bc]).value = "0.0";
      (*idbc_list[bc]).phase = "0.0";
      LOG_DBG(em) << "Set IDBC to HDBC (" << it->second->GetName() << ") -> value "
                    << (*(org_idbc.Last().second[bc])).value << " -> " << (*(it->second->GetIDBCList()[bc])).value;
    }

    // apply the new boundary condition
    it->second->SetBCs();
  }

  return org_idbc;
}

void ErsatzMaterial::ResetHDBC(StdVector<pair<SinglePDE*, IdBcList> >& org_idbc)
{
  // reset the original IDBC which were HDBC for the adjoint PDE
  for(unsigned int p = 0; p < org_idbc.GetSize(); p++)
  {
    // get the idbc list of the current pde
    IdBcList idbc_list = org_idbc[p].first->GetIDBCList();
    IdBcList org_list  = org_idbc[p].second;
    assert(idbc_list.GetSize() == org_list.GetSize());
    for(unsigned int bc = 0; bc < idbc_list.GetSize(); bc++) {
      (*idbc_list[bc]).value = (*org_list[bc]).value;
      (*idbc_list[bc]).phase = (*org_list[bc]).phase;
      LOG_DBG(em) << "Reset HDBC (" << org_idbc[p].first->GetName() << ") -> " << (*(org_idbc[p].first->GetIDBCList()[bc])).value;
    }
    org_idbc[p].first->SetBCs();
  }
}

void ErsatzMaterial::SolveAdjointProblem(Excitation* excite, Function* f)
{
  if (harmonic)
    SolveAdjointProblem<std::complex<double> > (excite, f);
  else
    SolveAdjointProblem<double> (excite, f);
}


template <class T>
void ErsatzMaterial::SolveAdjointProblem(Excitation* excite, Function* f)
{
  excite->Apply();
  switch(f->GetType())
  {
    // these objectives need their adjoint problems only for gradient calculation
    case Function::COMPLIANCE:
      if(IsTransient()){ // in transient case, everything has an adjoint
        Optimization::SolveAdjointProblem(excite, f);
      }
      break;
    case Function::TRACKING:
      // these objectives need their adjoint problems only for gradient calculation
      Optimization::SolveAdjointProblem(excite, f);

      if(!IsTransient()){ // transient solutions are read every timestep
        StorePDESolution(adjoint, *excite, f, 0, true, false, false, "adjoint");
      }

      // write back the solution s.th. CommitIteraion() makes StoreResults() properly.
      forward.Get(*excite)->Write(pde);
      break;

    case Function::OUTPUT:
    case Function::CONJUGATE_COMPLIANCE:
    case Function::ABS_DYN_OUTPUT_SQUARED:
    case Function::GLOBAL_DYNAMIC_COMPLIANCE:
    case Function::DYNAMIC_OUTPUT:
    case Function::ELEC_ENERGY:
    case Function::ENERGY_FLUX:
    case Function::STRESS:
    {
      // these objectives need their adjoint problems for the calculation of the objective value
      // they are directly solved after the StateProblem
      // they are no more solved for gradient calculation (this only works if the optimizer always evaluates the function before the gradient)
      // when doing linesearch, this slows down optimization if solution is only needed for gradient

      // the forward problem was already solved and stored !!

      // Set the rhs and solve for it
      SetAndSolveAdjointRHS<T>(*excite, f);

      // store the stuff -> no rhs but special handling of element results
      StorePDESolution(adjoint, *excite, f, 0, true, false, true, "adjoint");

      // write back the solution s.th. CommitIteraion() makes StoreResults() properly.
      forward.Get(*excite)->Write(pde);
      break;
    }

    default:
      assert(false);
  }
}



template <class T>
void ErsatzMaterial::SetAndSolveAdjointRHS(Excitation& excite, Function* f)
{
  // the adjoint RHS might be an output stuff, then the loads are changed.
  // save and restore them in any case.
  LoadList org_loads = assemble_->GetLoads();

  // set pseudo loads (if there are output nodes)
  if(f->NeedsSelectionVector())
    ConstructSelection(excite, f, true); // is actually already set for the forward calculation - who cares?

  // any adjoint PDE has HDBC instead of IDBC -> will be undone later ResetHDBC()
  StdVector<pair<SinglePDE*, IdBcList> > org_idbc = SetHDBC();

  // set the adjoint rhs
  ConstructAdjointRHS(excite, f);

  // calculate adjoint problem
  assemble_->GetAlgSys()->Solve(CreateAdjointAnalysisIdNode());

  ResetHDBC(org_idbc);

  // reset the original loads, they have been changed in the output case
  assemble_->SetLoads(org_loads);
}

void ErsatzMaterial::ConstructSelection(Excitation& excite, Function* f, bool alter_rhs)
{
  // in SolveStateProblem() the clean variant for the objective
  LoadList org_loads;
  if(!alter_rhs)
    org_loads = assemble_->GetLoads();

  // overwrite the assemble loads with "pseudo loads"s loads
  assemble_->SetLoads(output_nodes_);

  // set our own RHS but delete first as Assemble adds
  assemble_->GetAlgSys()->InitRHS();

  // assemble the output nodes
  assemble_->AssembleRHSLoads();

  // save the "pseudo loading" which is the selection as the rhs for the adjoint
  // This is exactly what has been constructed. Not that for an adjoint RHS it needs
  // post processing.
  adjoint.Get(excite, f)->Read(Solution::SEL_VECTOR, pde);

  if(!alter_rhs)
    assemble_->SetLoads(org_loads);

  LOG_DBG2(em) << "ConstructSelection: excite=" << excite.index << " f=" << f->ToString() << " alter=" << alter_rhs
               << " sel=" << adjoint.Get(excite, f)->GetVector(Solution::SEL_VECTOR)->ToString(1);
}


void ErsatzMaterial::ConstructRealAdjointRHS(Excitation& excite, Function* f)
{
  Vector<double> rhs; // own OLAS vector

  switch(f->GetType())
  {
  case Function::OUTPUT:
  {
    Vector<double>& l = adjoint.Get(excite, f)->GetRealVector(Solution::SEL_VECTOR);
    rhs.Resize(l.GetSize());
    rhs = l * -1.0;
    break;
  }
  case Function::STRESS:
  {
    rhs = CalcVonMisesStressVector(excite, f, true, false); // gives us exactly what we want
    break;
  }
  default:
    assert(false);
    break;
  }

  assemble_->GetAlgSys()->InitRHS(rhs);
  // assert(rhs.NormMax() != 0.0);
  LOG_DBG2(em) << "CARHS<double>: f=" << f->ToString() << " rhs before solving: " << rhs.ToString(1) << " max_norm=" << rhs.NormMax();
}

void ErsatzMaterial::ConstructComplexAdjointRHS(Excitation& excite, Function* f)
{
  // we handle only complex cases
  Vector<complex<double> >& u = forward.Get(excite)->GetComplexVector(Solution::RAW_VECTOR);
  Vector<complex<double> >& l = adjoint.Get(excite, f)->GetComplexVector(Solution::SEL_VECTOR);

  // create a OLAS vector
  Vector<complex<double> >& rhs = adjoint.Get(excite, f)->GetComplexVector(Solution::RHS_VECTOR);
  rhs.Resize(u.GetSize());

  switch(f->GetType())
  {
  case Function::DYNAMIC_OUTPUT:       // rhs is from "output loads" and set in adjoint...rhs
    // the correct conjugate_output case is -L * u*, always complex!
    for(unsigned int i = 0; i < rhs.GetSize(); i++)
      rhs[i] = -1.0 * l[i] * std::conj(u[i]);
    break;

  case Function::ABS_DYN_OUTPUT_SQUARED: // J = Re{<u,l>}^2 + Im{<u,l>}^2 .
  {
    // RHS = -0.5 * (2*<u_R,l> l - j 2*<u_I,l> l)
    //     = - <u_R,l> l + j <u_I,l> l
    //     = - Re{<u,l>} l + j Im{<u,l>} l
    // find the complex scalar product ul <u,l>
    complex<double> ul = u.Inner(l);
    LOG_DBG2(em) << "AdjustComplexAdjointRHS: <u,l> = " << ul;

    for(unsigned int i = 0; i < rhs.GetSize(); i++)
      rhs[i] = complex<double>(-1.0 * ul.real() * l[i].real(), ul.imag() * l[i].real());

  }
  case Function::CONJUGATE_COMPLIANCE: // rhs is from original excitation, we stored it in forward...rhs
  {
    forward.Get(excite)->Read(Solution::RHS_VECTOR, pde); // set
    Vector<complex<double> >& org_rhs = forward.Get(excite)->GetComplexVector(Solution::RHS_VECTOR); // read
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
    // 0.25 * j* omega * (Q-Q^T)^T * u^* where 0.25 is 0.5 from standard ajoint formulation and 0.5 from objective
    for(unsigned int i = 0, n = rhs.GetSize(); i < n; i++)
      rhs[i] = complex<double>(0, 0.25) * excite.GetOmega() * rhs[i];
    break;

  default:
    assert(true); // e.g. for ELEC_ENERGY the rhs is set in PiezoSIMP::ConstructAdjointRHS()
  }

  assemble_->GetAlgSys()->InitRHS(rhs);
  assert(rhs.NormMax() != 0.0);
  LOG_DBG2(em) << "CARHS<complex>: f=" << domain->GetDriver()->GetActStep(pde->GetName())
               << " rhs before solving: " << rhs.ToString(1);
}


void ErsatzMaterial::ConstructAdjointRHS(Excitation& excite, Function* f)
{
  // cannot be inlined due to linker problems
  if(harmonic) ConstructComplexAdjointRHS(excite, f);
          else ConstructRealAdjointRHS(excite, f);
}

void ErsatzMaterial::SetPDEs()
{
  // we either set MECH (almost always) or HEAT -> set the common pde attribute!
  pde = domain->GetSinglePDE("mechanic", false);
  if(pde != NULL)
    pdes[MECH] = pde;
  else
  {
    pde = domain->GetSinglePDE("heatConduction", true);
    pdes[HEAT] = pde;
  }

  // make it more smart when using energy flux for other pdes
  if(objectives.Has(Function::ENERGY_FLUX))
    pdes[ACOUSTIC] = domain->GetSinglePDE("acoustic", true);


  // ELEC is set in PiezoSIMP()
}



SinglePDE* ErsatzMaterial::ToPDE(Application app, bool throw_exception) const
{
  map<Application, SinglePDE*>::const_iterator it = pdes.find(app);
  if(it != pdes.end())
    return it->second;

  // nothing found
  if(throw_exception)
      EXCEPTION("No PDE '" << app << "' stored");

  return NULL;
}


Optimization::Application ErsatzMaterial::ToApp(DesignElement::Type dt)
{
  switch(dt)
  {
  case DesignElement::DENSITY:
    return MECH;
  case DesignElement::POLARIZATION:
    return ELEC;
  default:
    EXCEPTION("DesignType " << DesignElement::type.ToString(dt) << " doesn't map to Application");
  }
}


Optimization::Application ErsatzMaterial::ToApp(SinglePDE* pde) const
{
  if(pde->GetName() == "electrostatic") return ELEC;
  if(pde->GetName() == "mechanic") return MECH;
  if(pde->GetName() == "heatConduction") return HEAT;
  if(pde->GetName() == "acoustic") return ACOUSTIC;

  throw Exception("invalid");
}


void ErsatzMaterial::GetErsatzMaterialTensor(Matrix<double>& mat, Elem* elem, DesignElement::Type direction){
  linElastInt* form = (linElastInt*)GetForm(elem->regionId, pde, pde, "linElastInt");
  form->calcDMat(mat, elem, direction);
}

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
    StdVector<Unit*>& list = data_[funcs[fi]];
    for(unsigned int ts = 0; ts < list.GetSize(); ++ts)
      delete list[ts];
  }
}

void ErsatzMaterial::Solutions::Init(ErsatzMaterial* em)
{
  this->em_ = em;
}

void ErsatzMaterial::Solutions::Init(Function* f)
{
  assert(em_ != NULL);

  StdVector<Unit*>& list = data_[f];
  list.Resize(domain->GetDriver()->GetNumSteps());

  for(unsigned int ts = 0; ts < domain->GetDriver()->GetNumSteps(); ++ts)
    list[ts] = new Unit(em_);
}


ErsatzMaterial::Solutions::Unit::Unit(ErsatzMaterial* em)
{
  data.Resize(em->me->excitations.GetSize());
  for(unsigned int i=0; i < data.GetSize(); i++)
    data[i] = new Solution(em);
}

ErsatzMaterial::Solutions::Unit::~Unit()
{
  for(unsigned int i = 0; i < data.GetSize(); i++)
  {
    delete data[i];
    data[i] = NULL;
  }
}

ErsatzMaterial::Solution* ErsatzMaterial::Solutions::Get(Excitation& excitation, Function* f, unsigned int timestep)
{
  return Get(excitation.index, f, timestep);
}

ErsatzMaterial::Solution* ErsatzMaterial::Solutions::Get(int excitation_index, Function* f, unsigned int timestep)
{
  if(isForward){ // if this is true, f is ignored and forward_data_ is used to avoid one map access
    if(forward_data_ == NULL){
      if(data_.find(NULL) == data_.end()){
        Init((Function*)NULL);
      }
      forward_data_ = &data_[NULL];
    }
    return((*forward_data_)[timestep]->data[excitation_index]);
  }else{
    // do we have to init first?
    if(data_.find(f) == data_.end())
      Init(f);
    assert(data_.find(f) != data_.end());
    return data_[f][timestep]->data[excitation_index];
  }
}

StdVector<Function*> ErsatzMaterial::Solutions::GetFunctions() const
{
  StdVector<Function*> result;

  for(map<Function*, StdVector<Unit*> >::const_iterator it = data_.begin(); it != data_.end(); ++it)
  {
    LOG_DBG2(em) << "GetFunctions(): f=" << (it->first == NULL ? "NULL" : it->first->ToString());
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

SingleVector* ErsatzMaterial::Solution::Read(StorageType st, StdPDE* pde, Application app, bool save_sol)
{
  if (em_->harmonic)
    return Read<std::complex<double> > (st, pde, app, save_sol);
  else
    return Read<double> (st, pde, app, save_sol);
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
    Solutions::Unit* unit = sol.data_[f][time_step];
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
  LOG_DBG2(em) << "S:W pde=" << pde->GetName() << " h=" << harmonic << " vec=" << vec->ToString();
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
SingleVector* ErsatzMaterial::Solution::Read(StorageType st, StdPDE* pde, Application app, bool save_sol)
{
  SolutionType solt;
  switch(app)
  {
  case NO_APP: // up to now
    solt = MECH_DISPLACEMENT;
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
      int n = em_->design->GetNumberOfElements();

      // check for first call
      if(elem_vec.GetSize() == 0)
      {
        elem_vec.Resize(n);
        for(int ve = 0; ve < n; ve++)
          elem_vec[ve] = new Vector<T>;
      }

      // create an element list to gain the iterator in the loop
      ElemList elemList(grid);

      // store the results in our own structure
      for(int e = 0; e < n; e++)
      {
        DesignElement* de = &em_->design->data[e];

        elemList.SetElement(de->elem);
        const EntityIterator& it = elemList.GetIterator();

        pde->GetSolVecOfElement((Vector<T>&) *elem_vec[e], it, resinfo);
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
      Vector<T> tmp;

      if(st == RAW_VECTOR)
        bs->GetSolutionVal(tmp);
      else
        bs->GetRHSVal(tmp); // RHS_VECTOR and SEL_VECTOR

      SingleVector* ptr = NULL;
      switch(st)
      {
      case RAW_VECTOR:
        if(raw == NULL) raw = new Vector<T>(tmp.GetSize());
        ptr = raw;
        break;
      case RHS_VECTOR:
        if(rhs == NULL) rhs = new Vector<T>(tmp.GetSize());
        ptr = rhs;
        break;
      case SEL_VECTOR:
        if(select == NULL) select = new Vector<T>(tmp.GetSize());
        ptr = select;
        break;
      default:
        assert(false);
      }

      *ptr = tmp; // copy constructor

      return ptr;
    }
  }

  throw Exception("false");
}

// template instantiation stuff
template double ErsatzMaterial::CalcU1KU2<double>(TransferFunction* tf, StdVector<SingleVector*>& u1,
    Application app, StdVector<SingleVector*>& u2,
    SurfaceRef* rhs, double factor, CalcMode calcMode, Function* f, int res_idx);

template double ErsatzMaterial::CalcU1KU2<std::complex<double> >(TransferFunction* tf, StdVector<SingleVector*>& u1,
    Application app, StdVector<SingleVector*>& u2,
    SurfaceRef* rhs, double factor, CalcMode calcMode, Function* f,  int res_idx);
