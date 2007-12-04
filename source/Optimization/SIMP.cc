#include "Optimization/SIMP.hh"
#include "Optimization/DesignSpace.hh"
#include "Optimization/DesignElement.hh"
#include "Domain/domain.hh"
#include "Domain/surfElem.hh"
#include "General/exception.hh"
#include "Utils/basenodestoresol.hh"
#include "PDE/SinglePDE.hh"
#include "PDE/mechPDE.hh"
#include "Forms/baseForm.hh"
#include "Forms/SurfaceNormalInt.hh"
#include "CoupledPDE/DirectCoupledPDE.hh"
#include "Utils/StdVector.hh"
#include "Utils/vector.hh"
#include "Utils/result.hh"
#include "Utils/mathParser/mathParser.hh"
#include "Driver/assemble.hh"
#include "Driver/baseSolveStep.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "OLAS/olas.hh"
#include "OLAS/algsys/baseentrymanipulator.hh"
#include "OLAS/matvec/basematrix.hh"
#include "OLAS/matvec/stdmatrix.hh"

#include <string>

using namespace CoupledField;

using OLAS::StdMatrix;
using std::complex;

DECLARE_LOG(conditions)
DEFINE_LOG(conditions, "conditions")

DECLARE_LOG(simp)
DEFINE_LOG(simp, "simp")

Enum<SIMP::System> SIMP::system;

SIMP::SIMP() : Optimization()
{
  /** We store here the solution */
  forward_ = NULL;
  adjoint_ = NULL;
  output_vector_ = NULL;
  mathParserHandle_ = 0;

  pn = param->Get("optimization")->Get("SIMP");  
  system_ = system.Parse(pn->Get("system"));

  // region stuff
  regionId = domain->GetGrid()->RegionNameToId(pn->Get("region")->AsString());

  // find the right PDE which is the mechanic from the piezo
  mech = dynamic_cast<MechPDE*>(domain->GetSinglePDE("mechanic"));

  harmonic = BasePDE::IsComplex(mech->GetAnalysisType());

  // Get the assemble class
  assemble_ = mech->getPDE_assemble();

  // manipulate the PDE (the Integrator of its bilinear form) such that 
  // it is set up with the new densitiy values in the optimization loops
  GetForm(mech, mech, "linElastInt")->SetSolDependent(true);
  if(harmonic)
    GetForm(mech, mech, "MassInt")->SetSolDependent(true);

  // set up the design space elements, note PiezoSIMP have only POLARIZATION
  // this includes the transfer functions!
  StdVector<ParamNode*> design_list = pn->GetList("design");
  StdVector<ParamNode*> transfer_list = pn->GetList("transferFunction");
  StdVector<ParamNode*> result = pn->GetList("result");
  design = new DesignSpace((int)regionId, design_list, transfer_list, result);

  // There might be a filter regularization based on the design element.
  if(pn->Has("regularization", "type", "filter"))
  {
    StdVector<ParamNode*> list = pn->Get("regularization")->GetList("filter");
    // this is save for design=polarization
    for(unsigned int i = 0; i < list.GetSize(); i++)
      DesignElement::InitFilter(design->data, list[i]);
  }
  else 
  {
    if(pn->Has("regularization"))
      throw Exception("regularization not implemented"); 
  }

  // check our constraints, the shall have only valid designs
  for(unsigned int i = 0; i < constraints.GetSize();  i++)
    if(design->FindDesign(constraints[i].design, false) == -1)
      throw Exception("constraint " + constraints[i].ToString() + " operates on invalid design variable");
  for(unsigned int i = 0; i < outputs.GetSize();  i++)
    if(design->FindDesign(outputs[i].design, false) == -1)
      throw Exception("output constraint " + outputs[i].ToString() + " operates on invalid design variable");

  // now we know the constraints size -> PostInit design
  design->PostInit(constraints.GetSize());

  // give the domain this data, s.th. the ersatz material approch is appield
  domain->SetErsatzMaterial(design);      

  // add optimization results to the pde
  design->AppendOptimizationResults(mech); 

  // optionally write the densities to an xml file
  if(pn->Has("ersatzMaterial"))
    CreateErsatzMaterialFile(pn->Get("ersatzMaterial")->Get("file")->AsString(), design_list, transfer_list);
  else ersatzMaterialFile = NULL;  

  // has to be done after we have mech and we have design
  GetElementMatrix(GetForm(mech, mech, "linElastInt"), mechStiffness);
  if(harmonic)
    GetElementMatrix(GetForm(mech, mech, "MassInt"), mechMass);
  // handle surfaceNormal in PostInit()
  
  // we store the forward solution always 
  forward_ = new Solution(this);
  // it's cheap - even if not used
  adjoint_ = new Solution(this);

  // make a copy of the old iteration to calculate the move
  last_iteration.Resize(design->data.GetSize());

  // note the difference between function evaluations (line search) and iterations!   
  last_evaluation.Resize(design->data.GetSize());

  // this is needed for move calculation!
  design_tmp_.Resize(design->data.GetSize());
}

SIMP::~SIMP()
{
  if(forward_ != NULL)       { delete forward_; forward_ = NULL; }
  if(adjoint_ != NULL)       { delete adjoint_; adjoint_ = NULL; }
  if(output_vector_ != NULL) { delete output_vector_; output_vector_ = NULL; }
  
  // if write to file close the xml envelope and the file
  if(ersatzMaterialFile != NULL)
  {
    (*ersatzMaterialFile) << "</cfsErsatzMaterial>" << std::endl;
    ersatzMaterialFile->close();
    delete ersatzMaterialFile;
    ersatzMaterialFile = NULL;
  }

  // we do not necessary realy use one
  if(mathParserHandle_ != 0) domain->GetMathParser()->ReleaseHandle(mathParserHandle_);
  
  
  // "remove" the ersatzMaterial (=data) from the domain
  domain->SetErsatzMaterial(NULL);
}

void SIMP::PostInit()
{
  std::string objective = "'" + objectiveType.ToString(cost->type) + "'"; 

  // checks
  switch(cost->type)
  {
  case COMPLIANCE:
    if(harmonic) throw Exception(objective + " is only for static state problems");
    break;
    
  case GLOBAL_DYNAMIC_COMPLIANCE:
  case RADIATION:  
    if(!harmonic) throw Exception(objective + " is only for harmonic state problems");
    break;
    
  default:
    break;
  }
  
  // actions
  // note, that Radiation has its bilinear form already set in MechPDE::DefineIntegrators()
  switch(cost->type)
  {
  case OUTPUT:
  case CONJUGATE_OUTPUT:
  {
    // optimizing for nodes is mechanic, optimizing for loads is electrostatic, even if we SIMP mechanic!
    ParamNode* output = pn->Get("output");
    if(output->Has("nodes") && output->Has("loads")) 
      throw Exception("current implementation cannot optimize nodes and loads concurrently");

    if(output->Has("nodes"))
      mech->ReadLoads(output->GetList("nodes"), output_nodes_);
    if(output->Has("loads"))
      domain->GetSinglePDE("electrostatic")->ReadLoads(output->GetList("loads"), output_nodes_);

    if(output_nodes_.GetSize() == 0)
      throw Exception("no output optimization nodes/loads given");
    break;
  }
  case RADIATION:
  {
    // note, that the surface is not the optimization region with is all material
    std::string region = pn->Get("radiation")->Get("surfRegion")->Get("name")->AsString();
    RegionIdType sreg = domain->GetGrid()->RegionNameToId(region);
    BaseForm* form = assemble_->GetBiLinForm(sreg, mech, mech, "SurfaceNormalInt")->GetIntegrator();
    // the surface elements for the matrix are not design elements! 
    // Therefore provide an explicit element
    StdVector<Elem*> elems;
    domain->GetGrid()->GetElems(elems, sreg);
    assert(elems.GetSize() != 0);
    GetElementMatrix(form, surfaceNormal, elems[0]);
    break;
  }
  default:
    break;
  }
  
  Optimization::PostInit();  
}


void SIMP::CreateErsatzMaterialFile(const std::string& filename, StdVector<ParamNode*>& des, StdVector<ParamNode*>& tfs)
{
   ersatzMaterialFile = new std::ofstream(filename.c_str());
   if(ersatzMaterialFile == NULL) 
     throw Exception("cannot open file " + filename + " for writing");
     
   (*ersatzMaterialFile) << "<cfsErsatzMaterial>" << std::endl;
   (*ersatzMaterialFile) << "  <header>" << std::endl;
   for(unsigned int i = 0; i < des.GetSize(); i++)
   {
     (*ersatzMaterialFile) << "    ";
     des[i]->ToXML(*ersatzMaterialFile);
   } 
   for(unsigned int i = 0; i < tfs.GetSize(); i++)
   {
     (*ersatzMaterialFile) << "    ";
     tfs[i]->ToXML(*ersatzMaterialFile);
   } 
   (*ersatzMaterialFile) << "  </header>" << std::endl;
}

void SIMP::CommitIteration()
{
  // will write the cfs results and the log file
  Optimization::CommitIteration();

  if(ersatzMaterialFile != NULL)
  { 
    // add the entry, not tha the itteratation counter was incremented in base implementation
    *ersatzMaterialFile << "  <set id=\"" << (currentIteration-1) << "\">" << std::endl;

    for(unsigned int i = 0; i < design->data.GetSize(); i++)    
    {
      DesignElement* de = &design->data[i];
      *ersatzMaterialFile << "    <element nr=\"" << de->elem->elemNum << "\""
                          << " type=\"" << DesignElement::type.ToString(de->GetType()) << "\"" 
                          << " design=\"" << de->GetDesign(DesignElement::PLAIN) << "\""
                          << " gradient=\"" << de->GetObjectiveGradient(DesignElement::PLAIN) << "\""
                          << " filt_grad=\"" << de->GetObjectiveGradient(DesignElement::SMART) << "\""                                
                          << "/>" << std::endl; 
    }

    *ersatzMaterialFile << "  </set>" << std::endl;
  }
}

BaseForm* SIMP::GetForm(StdPDE* pde1, StdPDE* pde2, const std::string& integrator)
{
  return assemble_->GetBiLinForm(regionId, pde1, pde2, integrator)->GetIntegrator();
}




BiLinFormContext* SIMP::CreateSurfaceNormalMatrix(SinglePDE* pde, BaseMaterial* baseMat, shared_ptr<ResultInfo> result)
{
  // killme! this is called in mechPDE, it would be cooler if one could really add this
  // bilinear form in the optimization PostInit() ;(
  
  // check if we have this setting, otherwise return NULL
  if(!param->Has("optimization")) return NULL;
  // using the Enum stuff would be cooler but yet there is not prior SetEnums() calling
  if(param->Get("optimization")->Get("costFunction")->Get("type")->AsString() != "radiation") return NULL;

  // read what to optimize
  Grid* grid = domain->GetGrid();

  // killme works only for one region
  std::string region = param->Get("optimization")->Get("SIMP")->Get("radiation")->Get("surfRegion")->Get("name")->AsString();
  shared_ptr<EntityList> list = 
    grid->GetEntityList(EntityList::SURF_ELEM_LIST, region, EntityList::REGION);
  
  // will be automatically deleted in the BiLinFormContext destructor by assemble
  SurfaceNormalInt* sni = new SurfaceNormalInt(baseMat); 

  // this bilinear form is to be assembled in an extra matrix. it is
  // only used to calculate the rhs for radiation optimization
  BiLinFormContext * bilifc = new BiLinFormContext(sni, AUXILIARY);
  bilifc->SetPtPdes(pde, pde); // actually we don't need the pde
  
  bilifc->SetResults(result, result, list, list);

  return bilifc;
}


void SIMP::GetElementMatrix(BaseForm* form, Matrix<double>& out, Elem* elem)
{
  // create an element list to gain the iterator in the loop
  ElemList elemList(domain->GetGrid());
  
  // when there is no element given, we use the first from our design element.
  // For this process we disable all transfer functions. - they shoul be enabled!
  // This works only for isotropic material.

  // temporarily disables the transfer functions
  design->DisableTransferFunctions();
  
  // We have to do this because CalcElementMatrix asks the density via domain and
  // multiplies it with it's matrix!
  
  if(elem == NULL) elemList.SetElement(design->data[0].elem);
              else elemList.SetElement(elem);
  
  const EntityIterator& it = elemList.GetIterator();

  // get our element stiffness matrix -> it could be saved over the iteration loops!  
  form->CalcElementMatrix(out,const_cast<EntityIterator&>(it), const_cast<EntityIterator&>(it));
  LOG_DBG3(simp) << "CalcElemMatrix for " << form->GetName() << " -> " << out.ToString();

  // enable again our tranfer functions
  design->EnableTransferFunctions();
}

template <class T>
double SIMP::CalcObjective()
{
  design->WriteDesignToExtern(last_evaluation.GetPointer());
  
  switch(cost->type)
  {
     case COMPLIANCE:
     {
       // The compliance is sum over elements \rho^p u K u 
       // where u is the displacement and K the local stiffness matrix=elasticity tensor

       // there is a feature in MechPDE (mech->CalcEnergy<double>(regionId)) which is slower
       // but should be used when using non-uniform cells
       TransferFunction* tf = design->GetTransferFunction(DesignElement::DENSITY, MECH, true);
       // no derivative, set it directly, (+) u K u
       double val = CalcU1KU2(tf, forward_->elem[MECH], MECH, forward_->elem[MECH], false, false, 1.0);
       //double val = CalcU1KU2(tf, MECH, forward_, mechStiffness, MECH, forward_, false, false, 1.0);
       cost->SetValue(val);
       break;
     }
     case GLOBAL_DYNAMIC_COMPLIANCE:
     {
       // c = u^T transpose(u) -> "A note on sensitivity analysis of linear dynamic systems with
       //                          harmonic excitation"; Jakob S. Jensen; June 22, 2007
       CFSVector* cfsvec = forward_->raw[MECH];
       assert(cfsvec != NULL);
       assert(cfsvec->GetSize() != 0);
       Vector<complex<double> >& u = dynamic_cast<Vector<complex<double> >& >(*cfsvec);
       assert(u.GetSize() != 0);
       
       complex<double> csp = u * u; // the inner product is sum over u_i * conj(u_i);
       cost->SetValue(csp.real()); // handles also the case of a wrong inner product!
       break;
     }
     
     case OUTPUT:
     case CONJUGATE_OUTPUT:
     {
       // The output is <l,u> where l is -1* rhs of the adjoint pde
       // Here the rhs of the adoint pde is not -1 -> hence we do -1*<l,u> 
       //
       // We always take l from rhs and don't store it explicitly.
       // Note that one has to use the algsys RHS! The PDE RHS is still from the
       // forward simulation!
       assert(output_vector_ != NULL);
       Vector<T>& u = dynamic_cast<Vector<T> & >(*(forward_->raw[MECH]));
       Vector<T>& l = dynamic_cast<Vector<T> & >(*output_vector_);
       assert(u.GetSize() == l.GetSize());
       
       LOG_DBG3(simp) << "OUPTUT: adjoint rhs (l): " << l.ToString();
       LOG_DBG3(simp) << "OUPTUT: forward sol (u): " << u.ToString();            

       if(cost->type == OUTPUT)
       {
         // this is <l, u> which is for complex not really defined!
         T inner;
         u.Inner(l, inner);
         cost->SetValue(((complex<double>) inner).real());
         LOG_DBG2(simp) << "output <l,u>: " << inner << " -> " << cost->GetValue();
       }
       else
       {
         // this is <u,L conj(u)> and only defined for the harmonic case!
         if(!harmonic) throw Exception("'conjugateOutput' is only defined for harmonic!");
         
         double result = 0.0;
         // we loop over the vectors and do the scalar product by hand as we have
         // no diagonal matrix version of l
         for(unsigned int i = 0; i < u.GetSize(); i++ )
         {
           if(l[i] == 0.0) continue; // we skip this so we can make output for the real cases 

           complex<double> u_val = (complex<double>) u[i];
           // make sure we have no penalization stuff!
           assert(std::abs(u_val) < 1e15);
               
           double sp = std::real(std::conj(u_val) * l[i] * u_val);
           LOG_DBG2(simp) << "CalcObjective: " << std::conj(u_val) << " * " << l[i] << " * " << u_val << " -> " << sp;  
           result += sp;
         }
         LOG_DBG2(simp) << "output <u,L u*>: " << result;
       }
       break;            
     }
     case RADIATION:
     {
       // get surface normal matrix
       StdMatrix* snm = assemble_->GetAlgSys()->GetSysMat(AUXILIARY);
       assert(snm != NULL);
       // assert(out.G)assert(snm->GetMaxDiag() != 0.0);

       // Make an OLAS Vector and feed it with the solution vector
       Vector<T>& sol   = dynamic_cast<Vector<T>& >(*forward_->raw[MECH]);
       OLAS::Vector<T> olas_sol;
       // OLAS is one based!!
       olas_sol.Replace(sol.GetSize(), sol.GetPointer()-1, false);

       // The result is also an olas vector
       OLAS::Vector<T> olas_prod;
       olas_prod.Resize(olas_sol.GetSize());
       
       // multiply (S_n U)
       snm->Mult(olas_sol, olas_prod);
       LOG_DBG(simp) << "radiation objective: ||S_n*U||=" << olas_prod.NormEuclid();

       // scalar product U (S_n U)
       T tsp; 
       olas_sol.Inner(olas_prod, tsp);
       LOG_DBG(simp) << "radiation objective: U*S_n*U=" << tsp;
       // check that the scalar product is (numerical) real
       assert(((complex<double>) tsp).imag() < 1e-10 * ((complex<double>) tsp).real());
       double sp = ((complex<double>) tsp).real();
        
       // the objective (Du, Olhoff 2007) is pi=0.5 * gamma * c * omega^2* U*S_n*U
       double omega = mech->GetSolveStep()->GetActFreq();
       double c     = pn->Get("radiation")->Get("soundSpeed")->AsDouble();
       double gamma = pn->Get("radiation")->Get("specificMass")->AsDouble();
       
       double result = 0.5 * gamma * c * omega * omega * sp;
       LOG_DBG(simp) << " radiation objective: 0.5 * " << gamma << " * " << c  
                    << " *  " << omega << "^2 * " << sp << " -> " << result;
       cost->SetValue(result);
       break;
     }
          
     default: throw Exception("objective no handled");
  }         
       
  return cost->GetValue();
}

template <class T>
void SIMP::CalcSurfaceNormalTimesSolution(OLAS::Vector<T>& olas_prod)
{
  // get surface normal matrix
  StdMatrix* snm = assemble_->GetAlgSys()->GetSysMat(AUXILIARY);
  assert(snm != NULL);
  assert(snm->GetMaxDiag() != 0.0);

  // Make an OLAS Vector and feed it with the solution vector
  Vector<T>& sol   = dynamic_cast<Vector<T>& >(*forward_->raw[MECH]);
  OLAS::Vector<T> olas_sol;
  // OLAS is one based!!
  olas_sol.Replace(sol.GetSize(), sol.GetPointer()-1, false);
  LOG_DBG3(simp) << "Solution 'backcasted' to OLAS::Vector -> " << olas_sol.ToString(); 

  // The result is also an olas vector
  olas_prod.Resize(olas_sol.GetSize());
  
  // multiply (S_n U)
  snm->Mult(olas_sol, olas_prod);
  LOG_DBG(simp) << "||S_n*U||=" << olas_prod.NormEuclid();
  LOG_DBG3(simp) << "S_n * U -> " << olas_prod.ToString();
}

void SIMP::CalcObjectiveGradient(double* grad_out)
{
  TransferFunction* tf = design->GetTransferFunction(DesignElement::DENSITY, MECH, true);
  
  switch(cost->type)
  {
  case COMPLIANCE:
    // calculate the compliance which is according to 
    // "A 99 line topology optimization code written in Matlab"; O.Sigmund, 2001
    // -> dc/dx_e = -p * x_e ^(p-1) u_e^T k_0 u_e
    // Here we use a constant element stiffness matrix
    CalcU1KU2(tf, forward_->elem[MECH], MECH, forward_->elem[MECH], true);
    // CalcU1KU2(tf, MECH, forward_, mechStiffness, MECH, forward_, true);
    break;

  case GLOBAL_DYNAMIC_COMPLIANCE:
    CalcU1KU2(tf, forward_->elem[MECH], MECH, forward_->elem[MECH], true);
    // CalcU1KU2(tf, MECH, forward_, mechStiffness, MECH, forward_, true);
    break;

  case OUTPUT:
    // synthesis of compliant mechanism: As our adjoint PDE
    // c' = l K' u
    CalcU1KU2(tf, adjoint_->elem[MECH], MECH, forward_->elem[MECH], true);
    // CalcU1KU2(tf, MECH, adjoint_, mechStiffness, MECH, forward_, true);
    break;
    
  case RADIATION:
  {
    // copy & paste from CalcObjective
    double omega = mech->GetSolveStep()->GetActFreq();
    double c     = pn->Get("radiation")->Get("soundSpeed")->AsDouble();
    double gamma = pn->Get("radiation")->Get("specificMass")->AsDouble();
    
    // for constant rhs: grad = gamma * c * omega^2 * -1 * U_s^T * S' * U
    // U_s is adjoint solution, S' is derivative of S = stiffness, masss (, damping)
    double factor = -1.0 * gamma * c * omega * omega; 
    
    CalcU1KU2(tf, adjoint_->elem[MECH], MECH, forward_->elem[MECH], true, false, factor);
    break;
  }
  default: throw Exception("objective gradient no handled");
  }

  if(grad_out != NULL) design->WriteGradientToExtern(grad_out, DesignElement::DENSITY, 
                               DesignElement::COST_GRADIENT, DesignElement::SMART);
}


template <class T>
double SIMP::CalcU1KU2(TransferFunction* tf, StdVector<CFSVector*>& u1, 
                       Application k, StdVector<CFSVector*>& u2, 
                       bool derivative, bool add, double factor)
{
  LOG_DBG2(simp) << "CalcU1KU2(): tf=" << tf->ToString() << " #u1=" << u1.GetSize()
                << " #u2=" << u2.GetSize() << " derivative=" << derivative
                << " add=" << add << " factor=" << factor;
  
  double sum = 0.0;
  int elements = design->GetNumberOfElements();
  int base = design->FindDesign(tf->GetDesign()) * elements; 
  Vector<T> tmp;
  assert(u1.GetSize() != 0);
  assert(u1.GetSize() == u2.GetSize());

  LOG_DBG3(simp) << "elements=" << elements << " base=" << base; 
 
  
  // create an element list to gain the iterator in the loop
  ElemList elemList(domain->GetGrid());
  
  for(int i = base; i < base+elements; i++)    
  {
    Vector<T>& u1_vec = dynamic_cast<Vector<T>& >(*u1[i]);
    Vector<T>& u2_vec = dynamic_cast<Vector<T>& >(*u2[i]);
    
    DesignElement* de = &design->data[i];

    LOG_DBG3(simp) << "u1:" << i << ": " << u1_vec.ToString();
    LOG_DBG3(simp) << "u2:" << i << ": " << u2_vec.ToString();

    // compliance: -u_e^T K_0 u_e
    // mechanism:(-)l_e^T K_0 u_e
    
    // this handles the (derivative) of the transfer function 
    CalcElementKU2(&u2_vec, k, de, derivative, &tmp);
    T sp = tmp * u1_vec; // here without transfer function
    
    // when doing complex Jensen 22.07.07 shows that we always have 2 * Re(lamda * grad S * u)
    // the factor gives the negative sign
    // in real case it is simple value = factor * sp. 
    // factor shall be +/- 1!
    double value = (harmonic ? 2 : 1) * ((complex<double>) sp).real() * factor;

    if(derivative) 
    {
      // do we have to sum up (-> piezo) or simply set?
      if(add) value += de->GetObjectiveGradient(DesignElement::PLAIN);
      de->SetObjectiveGradient(value);
    } 

    LOG_DBG3(simp) << "u1*K0*u2 = " << sp << " -> " << value;

    sum += value;
  }
  
  return sum;
}

template <class T>
void SIMP::CalcElementKU2(const CFSVector* cfs_in, Application app, DesignElement* de, bool derivative, CFSVector* cfs_out)
{
  assert(cfs_in->GetSize() != 0);  
  
  TransferFunction* tf = NULL;
  double k_factor = 0.0; // comes from the transfer function, might be derivative 
  double m_factor = 0.0; 
  
  const Vector<T>& in  = dynamic_cast<const Vector<T>& >(*cfs_in);
  Vector<T>& out = dynamic_cast<Vector<T>& >(*cfs_out);
  Vector<T> result;
  
  switch(app)
  {
  case SURFACE_NORMAL:
    // this does not depend on a transfer fuction!
    assert(surfaceNormal.GetSizeCol() == in.GetSize());
    // only with the operators we can multiply a Matrix<double> with Vector<complex>
    out = surfaceNormal * in;
    break;
    
  case MECH:
    // check the stiffness matrix is in proper size
    assert(mechStiffness.GetSizeCol() == in.GetSize());
    // get the transfer function for our design element (e.g. DENSITY, MECH)
    tf = design->GetTransferFunction(de->GetType(), app);
    // with simp this is e.g. de->design^3 or 3*de->design^2 as factor for K
    k_factor = derivative ? tf->Derivative(de) : tf->Transform(de);
    result = mechStiffness * in;
    result *= k_factor;
    out = result;
    LOG_DBG3(simp) << "CalcElementKU2 MECH: " << k_factor << "*K*u -> " << out.ToString();

    // In the harmonic case we have (K + i*omega*C - omega^2*M)*u
    if(harmonic)
    {
      // We do it component wise and already have K*u, now do -omega^2*M*u
      assert(mechMass.GetSizeCol() == in.GetSize());
      // the mass is mostly w/o penalization
      tf = design->GetTransferFunction(de->GetType(), app);
      m_factor = derivative ? tf->Derivative(de) : tf->Transform(de);
      // another factor comes from the omage 
      double omega = mech->GetSolveStep()->GetActFreq();
      // alltogether, all real
      result = mechMass * in;
      result *= (m_factor * -1.0 * omega * omega);
      LOG_DBG3(simp) << m_factor << "*-omega^2*M*u -> " << result.ToString();
      out += result;
      LOG_DBG3(simp) << "K*u-x*M*u -> " << out.ToString();
      
      // do we have damping (C = alpha*M+beta*K) -> this is pure imaginary!
      if(mech->GetDampingList().find(regionId) != mech->GetDampingList().end()
          && mech->GetDampingList().find(regionId)->second == RAYLEIGH)
      {
        // we need a math parser (is an unsigned int)
        if(mathParserHandle_ == 0) mathParserHandle_ = domain->GetMathParser()->GetNewHandle();
        
        // the alpha and beta might be calculated and adjusted, get them
        // from the integrators in the form as they are used for the state problem!
        std::string txt = assemble_->GetBiLinForm(regionId, mech, mech, "linElastInt")->GetSecMatFac();
        domain->GetMathParser()->SetExpr(mathParserHandle_, txt);
        double alpha_k = domain->GetMathParser()->Eval(mathParserHandle_);
        LOG_DBG2(simp) << "alpha_k: '" << txt << "' -> " << alpha_k;
        assert(alpha_k > 0.0);
        
        // actually, all element matrices are <double>, but in the harmonic case the
        // result and and in are Vector<complex> (therefore we need the overloaded operators)
        result = mechStiffness * in;
        result *= (alpha_k * omega * k_factor);
        LOG_DBG3(simp) << "alpha*omega*K*u -> " << result.ToString();
        // make it really complex
        complex<double> unit(0,1);
        for(unsigned int i = 0; i < result.GetSize(); i++)
          ((complex<double>) result[i]) *= unit;
        out += result;
        LOG_DBG3(simp) << "(K-omega^2*M+i*omega*alpha_k*K)*u -> " << out.ToString();
        
        // now alpha_m
        txt = assemble_->GetBiLinForm(regionId, mech, mech, "MassInt")->GetSecMatFac();
        domain->GetMathParser()->SetExpr(mathParserHandle_, txt);
        double alpha_m = domain->GetMathParser()->Eval(mathParserHandle_);
        LOG_DBG2(simp) << "alpha_m: '" << txt << "' -> " << alpha_m;
        assert(alpha_m > 0.0);
        
        result = mechMass * in;
        result *= (alpha_m * omega * m_factor);        

        LOG_DBG3(simp) << "alpha*omega*M*u -> " << result.ToString();
        // make it really complex
        for(unsigned int i = 0; i < result.GetSize(); i++)
          ((complex<double>) result[i]) *= unit;
        out += result;        
        LOG_DBG3(simp) << "(K-omega^2*M+i*omega*(alpha_k*K+alpha_m*M))*u -> " << out.ToString();
      }
    }
    break;
  default:
    assert(false); // other cases should be handled in PiezoSIMP
  } // end switch 
}


double SIMP::CalcConstraint(Condition* constraint, bool derivative, double* grad_out)
{
  // handle default parameter
  if(constraint == NULL)
  {
    if(constraints.GetSize() != 1)
      throw Exception("default constraint only valid with exactly one constraint");
    constraint = &constraints[0];  
  }
  double result = 0.0;
  
  switch(constraint->GetName())
  {
    case Condition::VOLUME:
         result = CalcVolume(derivative, constraint, grad_out);
         break; 
         
    case Condition::GREYNESS:
         result = CalcGreyness(derivative, constraint, grad_out);
         break;     

    case Condition::GAUSS_GREYNESS:
         result = CalcGreyness(derivative, constraint, grad_out);
         break;     
      
    default: throw Exception("Constraint notimplemented", __FILE__, __LINE__);  
  }
  
  LOG_TRACE2(simp) << "CalcConstraint " << constraint->ToString() 
                  << " derivative=" << derivative << " -> " << result;
  return result;                 
}

double SIMP::CalcVolume(bool derivative, Condition* constraint, double* grad_out)
{
  double sum = 0.0;
  int counter = 0;

  double fraction = 1.0 / (double) design->GetNumberOfElements();

  // loop over all elements to set 0.0 for not relevant gradients  
  for(unsigned int i = 0; i < design->data.GetSize(); i++)    
  {
    DesignElement* de = &design->data[i];
    bool relevant = constraint->design == DesignElement::DEFAULT || constraint->design == de->GetType();    
    
    if(derivative) 
    {
      double val = relevant ? fraction : 0.0;
      de->SetConstraintGradient(constraint, val);
      if(grad_out != NULL) grad_out[i] = val;
    }
    else
    {
      if(relevant)
      {
        sum += de->GetDesign(DesignElement::PLAIN);
        counter++; 
      }
    }
  }
 
  return !derivative ? sum / (double) counter : -1.0;
}



double SIMP::CalcGreyness(bool derivative, Condition* constraint, double* grad_out)
{
  double greyness = 0.0; // element greyness
  int counter = 0; // to make it sure for different design variables!
  
  double lb, ub, span, org_value, value, grad, eval;
  lb = ub = value = grad = eval = span = 0.0;
  bool gauss = constraint->GetName() == Condition::GAUSS_GREYNESS;
  // only used in gauss case
  double h = constraint->parameter;
  
  // we have to divide the gradients by their relalive volume = fraction
  double fraction = constraint->design == DesignElement::DEFAULT ?
                    design->data.GetSize() : design->GetNumberOfElements(); 
  
  // go over the complemte design space to set gradients of other types to 0
  for(unsigned int i = 0; i < design->data.GetSize(); i++)    
  {
    DesignElement* de = &design->data[i];
    bool relevant = constraint->design == DesignElement::DEFAULT || constraint->design == de->GetType();
    if(relevant)
    {
      lb = de->GetLowerBound();
      ub = de->GetUpperBound();
      span = ub-lb;
      org_value = design->data[i].GetDesign(DesignElement::PLAIN);

      if(gauss)
      {
        // We normalize for a design variable to [-1;1] ! 
        // nothing to do for polarization [-1:1] but an issue for density [0.001:1]
        value = (org_value - lb) * (2.0/span) - 1.0;  
      } 
      else
      {
        // We normalize for a design variable from [0;1]
        // this has minor effect on density [0.001;1] but is important
        // for polarization[-1;1]
        value = (org_value - lb) / span;
      }
    }
    
    if(derivative) 
    {
      if(relevant && gauss)
      {
        // a gauss type greyness (own idea, but most probably not unique :))
        // y = is a normalized version of x for the range [-1:1] (see above!)
        // f(y)=(e^(-1 * (y^2)/h) - e^(-1*1/h))/(1-e^(-1*1/h))
        // f'(y)=-2*y/h * e^(-1 * (y^2)/h)/(1-e^(-1*1/h))
        grad = (2.0/span) * -2.0 * (value/h) * exp(-1.0 * value * value / h) / (1-exp(-1.0/h));
      }
      if(relevant && !gauss)
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
      design->data[i].SetConstraintGradient(constraint, grad);
      if(grad_out != NULL) grad_out[i] = grad;
    }
    else // not derivative but function evaluation
    {
       // we normalized the greyness to value from [0;1]
       if(relevant)
       {
         if(gauss) eval = (exp(-1.0 * value * value / h) - exp(-1.0/h))
                         / (1-exp(-1.0/h));
              else eval = 4.0* (1-value) * (value);       
         greyness += eval;
         counter++;
       }
    }
    LOG_DBG3(conditions) << Condition::name.ToString(constraint->GetName()) 
                         << " derive=" << derivative << " relevant=" << relevant
                         << " elem " << de->elem->elemNum << " value: " << org_value
                         << " -> " << value << " grad=" << grad << " eval=" << eval
                         << " fraction=" << fraction << " counter=" << counter 
                         << " h=" << h;
 
  }
 
  return greyness / (double) counter;
}


void SIMP::SolveStateProblem()
{
  
  switch(cost->type)
  {
    // In the compliance case we store the solution as forward solution for CalcU1KU2
    case COMPLIANCE:
      Optimization::SolveStateProblem();

      // store solution element wise for gradient and vector for objective
      forward_->ReadSolution(Solution::ELEMENT_VECTORS, mech, MECH);
      break;       
  
    // when our objective is output we have to solve also the adjoint problem
    case GLOBAL_DYNAMIC_COMPLIANCE:
    case RADIATION:
    case OUTPUT:
      SolveAdjointProblem(mech, NULL);
      break;

    default:
      Optimization::SolveStateProblem();
      break;
  }
}


template <class T>
void SIMP::SetAndSolveAdjointRHS()
{
  double time = mech->GetSolveStep()->GetActTime();
  
  switch(cost->type)
  {
  case OUTPUT:
  case CONJUGATE_OUTPUT:
  {
    // overwrite the assemble loads with our adjoint rhs loads
    // Note, that in case of inhomogeneous dirichlet values there might be the
    // penalty term on the rhs!
    LoadList org_loads = assemble_->GetLoads();
    assemble_->SetLoads(output_nodes_);

    // set our own RHS but delete first as Assemble adds
    assemble_->GetAlgSys()->InitRHS();

    assemble_->AssembleRHSLoads(time);

    // save this rhs, it has no idbc and penalization set yet
    T* rhs_ptr;
    int rhs_size = assemble_->GetAlgSys()->GetRHSVal(rhs_ptr);
    if(output_vector_ == NULL) output_vector_ = new Vector<T>(rhs_size);
    dynamic_cast<Vector<T> *>(output_vector_)->Replace(rhs_size, rhs_ptr, false);
    
    if(cost->type == CONJUGATE_OUTPUT)
    {
      // the correct conjugate_output case is -L * u* -- always complex!
      Vector<complex<double> >& rhs = dynamic_cast<Vector<complex<double> >& >(*output_vector_);
      Vector<complex<double> >& u   = dynamic_cast<Vector<complex<double> >& >(*forward_->raw[MECH]);

      LOG_DBG2(simp) << "SetAndSolveAdjointRHS: pure output vector w/o idbc: " << rhs.ToString(); 
      
      for(int i = 0; i < rhs_size; i++)
      {
        // shall be pure real
        assert( ((complex<double>) rhs[i]).imag() == 0);
        // should not have penalization!
        assert( ((complex<double>) rhs[i]).real() < 1e15);
        rhs[i] *= std::conj(u[i]);
      }
      
      LOG_DBG2(simp) << "SetAndSolveAdjointRHS: conjugate output adjoint rhs before solving: " << rhs.ToString();
      assemble_->GetAlgSys()->InitRHS((double*) rhs.GetPointer());
    }
    
    // calculate adjoint problem
    assemble_->GetAlgSys()->Solve(GetSolveComment() + "_adjoint");

    // reset the original loads
    assemble_->SetLoads(org_loads);
    
    break;
  }
  
  case GLOBAL_DYNAMIC_COMPLIANCE:
  {
    // S lambda = -conj(u);

    // GLOBAL_DYNAMIC_COMPLIANCE can only be complex!
    complex<double>* rhs = NULL;
    int size = assemble_->GetAlgSys()->GetRHSVal(rhs);
    assert(size != 0);
    
    Vector<complex<double> >& u = dynamic_cast<Vector<complex<double> >& >(*forward_->raw[MECH]);

    for(int i = 0; i < size; i++)
      rhs[i] = -1.0 * std::conj(u[i]);
    
    assemble_->GetAlgSys()->Solve(GetSolveComment() + "_adjoint");
    break;
  }
  
  case RADIATION:
  {
    // the rhs for radiation is S_n * U
    
    // We have to use an OLAS Vector here
    OLAS::Vector<std::complex<double> > rhs;
    CalcSurfaceNormalTimesSolution(rhs); // does S_n * U as we also need it in the objective
    
    // rhs is OLAS 1-based but handles this internally
    assemble_->GetAlgSys()->InitRHS(const_cast<const OLAS::Vector<std::complex<double> >* >(&rhs));

    // calculate adjoint problem
    assemble_->GetAlgSys()->Solve(GetSolveComment() + "_adjoint");
    break;
  }
  default:
    assert(false);
  }
}

template <class T>
void SIMP::SolveAdjointProblem(StdPDE* mech, StdPDE* elec)
{
  // mech is always present only elec is an option (PiezoSIMP=
  assert(mech != NULL);
  
  // calculate forward problem.
  Optimization::SolveStateProblem();
  
  // store solution element wise for gradient and vector for objective
  forward_->ReadSolution(Solution::ELEMENT_VECTORS, mech, MECH);
  forward_->ReadSolution(Solution::RAW_VECTOR, mech, MECH);

  LOG_DBG3(simp) << "forward mech solution: " << forward_->raw[MECH]->ToString();
  
  if(elec != NULL)  {
    forward_->ReadSolution(Solution::ELEMENT_VECTORS, elec, ELEC);
    forward_->ReadSolution(Solution::RAW_VECTOR, elec, ELEC);
    LOG_DBG2(simp) << "forward elec solution: " << forward_->raw[ELEC]->ToString();
  }
    
  // Set the rhs!
  SetAndSolveAdjointRHS<T>();

  // store the solution in the PDE! This is necessary to read the element vector
  T* ptr;
  int length = assemble_->GetAlgSys()->GetSolutionVal(ptr);
  mech->SaveSolution(ptr, length);
  // now store own version
  adjoint_->ReadSolution(Solution::ELEMENT_VECTORS, mech, MECH);
  if(elec != NULL)  {
    elec->SaveSolution(ptr, length);
    adjoint_->ReadSolution(Solution::ELEMENT_VECTORS, elec, ELEC);
  }
  
  // write back the solution s.th. CommitIteraion() makes StoreResults() properly.
  assert(length == (int) forward_->raw[MECH]->GetSize());
  forward_->WriteSolution(mech, MECH);
  if(elec != NULL) {
    assert(length == (int) forward_->raw[ELEC]->GetSize());
    forward_->WriteSolution(elec, ELEC);
  }
}


SIMP::Solution::Solution(SIMP* simp)
{
  this->simp = simp;
}

SIMP::Solution::~Solution()
{
  std::map<Application, CFSVector* >::iterator raw_iter;
  for(raw_iter = raw.begin(); raw_iter != raw.end(); raw_iter++)
    delete raw_iter->second;

  std::map<Application, StdVector<CFSVector* > >::iterator elem_iter;
  for(elem_iter = elem.begin(); elem_iter != elem.end(); elem_iter++)
  {
    StdVector<CFSVector* >& data = elem_iter->second;
    for(unsigned int i = 0; i < data.GetSize(); i++)
    {
      delete data[i];
    }
  }

  /** This is the algsys solution vector */
  std::map<Application, CFSVector* > raw;
  
}

template <class T>
void SIMP::Solution::WriteGenericSolution(StdPDE* pde, Application app)
{
  T* ptr = NULL; 
  raw[app]->GetPointer(ptr);
  assert(ptr != NULL);
  assert(raw[app]->GetSize() != 0);
  pde->SaveSolution(ptr, raw[app]->GetSize());
}


template <class T>
void SIMP::Solution::ReadGenericSolution(StorageType st, StdPDE* pde, Application app)
{
  switch(st)
  {
    case ELEMENT_VECTORS:
    {
      // we save the element vectors in elem_vec. Might be empty the first call
      StdVector<CFSVector*>& elem_vec = elem[app];
      int n = simp->design->GetNumberOfElements();
      
      // check for first call
      if(elem_vec.GetSize() == 0)
      { 
        elem_vec.Resize(n);
        for(int ve = 0; ve < n; ve++)
          elem_vec[ve] = new Vector<T>;
      }

      // create an element list to gain the iterator in the loop
      ElemList elemList(domain->GetGrid());

      // store the results in our own structure
      for(int e = 0; e < n; e++)
      {
        DesignElement* de = &simp->design->data[e];

        elemList.SetElement(de->elem);
        const EntityIterator& it = elemList.GetIterator();

        SolutionType solt = app == MECH ? MECH_DISPLACEMENT : ELEC_POTENTIAL;
        pde->GetSolVecOfElement((Vector<T>&) *elem_vec[e], it, pde->GetResultInfo(solt));
      }
    }
    case RAW_VECTOR:
    {
      // It is best to get the RHS from the algebraic system (OLAS), the PDE
      // might not check about a further adjont calculation

      // get access to solution
      T* ptr;
      int size = simp->assemble_->GetAlgSys()->GetSolutionVal(ptr);
      
      // check for first call
      if(raw.find(app) == raw.end())
        raw[app] = new Vector<T>(size);
      
      CFSVector* raw_vec = raw[app];

      for(int i = 0; i < size; i++)
        raw_vec->SetEntry(i, ptr[i]);
    }
  }
}
