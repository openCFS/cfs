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
#include "Forms/linSurfForm.hh"
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

DECLARE_LOG(em)
DEFINE_LOG(em, "ersatzMaterial")

Enum<ErsatzMaterial::Method> ErsatzMaterial::method;

ErsatzMaterial::ErsatzMaterial() : Optimization()
{
  /** We store here the solution */
  forward_ = NULL;
  adjoint_ = NULL;
  output_vector_ = NULL;

  pn = param->Get("optimization")->Get("ersatzMaterial");
  
  method_ = method.Parse(pn->Get("method"));

  // region stuff -> to be extended by Bastian :)
  regionId = domain->GetGrid()->RegionNameToId(pn->Get("region")->AsString());

  // we assume always to have a mechanic pde. in pde[] we might have more pdes.
  // historically "mech" is a pointer to pde[0]
  // when there is no more mech, then "mech" might be replaced by pde[0]
  mech = dynamic_cast<MechPDE*>(domain->GetSinglePDE("mechanic"));
  assert(pde.GetSize() == 0);
  pde.Push_back(mech);

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
  design = new DesignSpace((int)regionId, design_list, transfer_list, result, method_);

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
  if(pn->Has("export"))
    CreateErsatzMaterialFile(pn->Get("export")->Get("file")->AsString(), design_list, transfer_list);
  else ersatzMaterialFile = NULL;  
  
  // we store the forward solution always 
  forward_ = new Solution(this);
  // it's cheap - even if not used
  adjoint_ = new Solution(this);

  // make a copy of the old iteration to calculate the move
  last_iteration.Resize(design->data.GetSize());

  // note the difference between function evaluations (line search) and iterations!   
  last_evaluation.Resize(design->data.GetSize());
}

ErsatzMaterial::~ErsatzMaterial()
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
  
  // "remove" the ersatzMaterial (=data) from the domain
  domain->SetErsatzMaterial(NULL);
}

void ErsatzMaterial::PostInit()
{
  std::string objective = "'" + objectiveType.ToString(cost->type) + "'"; 

  // checks
  switch(cost->type)
  {
  case OUTPUT:
  case COMPLIANCE:
    if(harmonic) throw Exception(objective + " is only for static state problems");
    break;
    
  case CONJUGATE_OUTPUT:
  case GLOBAL_DYNAMIC_COMPLIANCE:
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
    ParamNode* output = param->Get("optimization")->Get("costFunction")->Get("output");
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
  // radiation is done in SIMP!
  default:
    break;
  }
  
  Optimization::PostInit();  
}


void ErsatzMaterial::CreateErsatzMaterialFile(const std::string& filename, StdVector<ParamNode*>& des, StdVector<ParamNode*>& tfs)
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

void ErsatzMaterial::StoreResults(double step_val)
{
  // check if we have to play with the results, forward is default mode!
  if(pn->Has("commit") && pn->Get("commit")->Get("mode")->AsString() != "forward")
  {
    if(pn->Get("commit")->Get("mode")->AsString() == "both_cases")
    {
      for(unsigned int i = 0; i < pde.GetSize(); i++)
        forward_->WriteSolution(pde[i], ToApp(pde[i]));
      
      Optimization::StoreResults(step_val);
    }
      
    // "forward" is not possible -> hence it is both_cases or adjoint
    if(adjoint_->raw.empty()) EXCEPTION("Ajdoint solution cannot be written on 'commit' because there is none");

    for(unsigned int i = 0; i < pde.GetSize(); i++)
      adjoint_->WriteSolution(pde[i], ToApp(pde[i]));
    
    step_val = step_val == -1 ? currentIteration + 0.5 : step_val + 0.5;
  }
  
  // call real implementation in Optimization
  Optimization::StoreResults(step_val);  
}


void ErsatzMaterial::CommitIteration()
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

BaseForm* ErsatzMaterial::GetForm(StdPDE* pde1, StdPDE* pde2, const std::string& integrator)
{
  return assemble_->GetBiLinForm(regionId, pde1, pde2, integrator)->GetIntegrator();
}


void ErsatzMaterial::GetElementMatrix(BaseForm* form, Matrix<double>& out, Elem* elem)
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
  LOG_DBG3(em) << "CalcElemMatrix for " << form->GetName() << " -> " << out.ToString();

  // enable again our tranfer functions
  design->EnableTransferFunctions();
}




template <class T>
double ErsatzMaterial::CalcObjective()
{
  design->WriteDesignToExtern(last_evaluation.GetPointer());
  
  switch(cost->type)
  {
     case COMPLIANCE:
     {
       // compliance is easier computed using f^T u on nodes with force 
       // to avoid any work for assembling force again, we simply calculate solution times rhs from the system
       double val = 0.0;
       Vector<double>& u = dynamic_cast< Vector<double>& > (*forward_->raw[MECH]); 
       double* f; 
       const int size = assemble_->GetAlgSys()->GetRHSVal(f);
       
       assert(size == (int) u.GetSize());

       for(int i= 0; i < size; i++) 
         val += u[i] * f[i]; 
       
       
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
       
       LOG_DBG3(em) << "OUPTUT: adjoint rhs (l): " << l.ToString();
       LOG_DBG3(em) << "OUPTUT: forward sol (u): " << u.ToString();            

       if(cost->type == OUTPUT)
       {
         // this is <l, u> which is for complex not really defined!
         T inner;
         u.Inner(l, inner);
         cost->SetValue(((complex<double>) inner).real());
         LOG_DBG2(em) << "output <l,u>: " << inner << " -> " << cost->GetValue();
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
           LOG_DBG2(em) << "CalcObjective: " << std::conj(u_val) << " * " << l[i] << " * " << u_val << " -> " << sp;  
           result += sp;
         }
         LOG_DBG2(em) << "output <u,L u*>: " << result;
       }
       break;            
     }
          
     default: throw Exception("objective no handled");
  }         
       
  return cost->GetValue();
}


template <class T>
double ErsatzMaterial::CalcU1KU2(TransferFunction* tf, StdVector<CFSVector*>& u1, 
                       Application app, StdVector<CFSVector*>& u2, 
                       SurfaceRef* rhs, bool add, double factor)
{
  LOG_DBG2(em) << "CalcU1KU2(): tf=" << tf->ToString() << " #u1=" << u1.GetSize()
                 << " app=" << application.ToString(app) << " #u2=" << u2.GetSize() << " add=" 
                 << add << " factor=" << factor << " rhs=" << (rhs == NULL ? "NULL" : rhs->ToString(1));

  // This solves <l,K'*u-f'> or <u1, K' * u2 - f'> for all elements
  
  assert(u1.GetSize() != 0);
  assert(u1.GetSize() == u2.GetSize());
  
  double sum = 0.0;
  // the dimenions of our matrix is determinded by u1_vec and u2_vec.
  // mat will be filled by SetElementK where also the derivative form most cases is built in
  Matrix<T> mat(u1[0]->GetSize(), u2[0]->GetSize());
  Vector<T> k_u2(u1[0]->GetSize());
  
  TransferFunction* rtf = rhs != NULL && rhs->valid ? design->GetTransferFunction(tf->GetDesign(), rhs->app) : NULL;

  // traverse over our lements 
  int elements = design->GetNumberOfElements();
  int base     = design->FindDesign(tf->GetDesign()) * elements;
  LOG_DBG3(em) << "elements=" << elements << " base=" << base;
  
  // create an element list to gain the iterator in the loop
  ElemList elemList(domain->GetGrid());
              
  for(int e = 0 ; e < elements; e++)    
  {
    Vector<T>& u1_vec = dynamic_cast<Vector<T>& >(*u1[e]);
    Vector<T>& u2_vec = dynamic_cast<Vector<T>& >(*u2[e]);
        
    DesignElement* de = &design->data[e + base];

    LOG_DBG3(em) << "u1:" << e << ": " << u1_vec.ToString();
    LOG_DBG3(em) << "u2:" << e << ": " << u2_vec.ToString();

    // <u1, K' * u2 - f'> -> find "K'"
    SetElementK(de, app, dynamic_cast<CFSMatrix*>(&mat));
    LOG_DBG3(em) << "mat: " << mat.ToString();
    
    // <u1, K' * u2 - f'> -> calc "K' * u2"
    k_u2 = mat * u2_vec;
    LOG_DBG3(em) << "mat * u2: " << k_u2.ToString();
    
    // <u1, K' * u2 - f'> -> calc "- f'"
    if(rtf != NULL) SubstractGradSurfaceRHS(de, rtf, rhs, k_u2);    
    LOG_DBG3(em) << "-f': " << k_u2.ToString();
   
    // <u1, K' * u2 - f'> -> calc "<u1, *>"
    T sp = u1_vec * k_u2; 

    // when doing complex Jensen 22.07.07 shows that we always have 2 * Re(lamda * grad S * u)
    // the factor gives the negative sign
    // in real case it is simple value = factor * sp. 
    // factor shall be +/- 1!
    double value = (harmonic ? 2 : 1) * ((complex<double>) sp).real() * factor;

    // do we have to sum up (-> piezo) or simply set?
    if(add)
    {
      LOG_DBG3(em) << "CalcU1KU2():de=" << de->elem->elemNum << " (" << de->type.ToString(de->GetType()) << ") "
                   << " app=" << tf->GetApplication()
                   << " old=" << de->GetObjectiveGradient(DesignElement::PLAIN)
                   << " new=" << de->GetObjectiveGradient(DesignElement::PLAIN) + value;

      value += de->GetObjectiveGradient(DesignElement::PLAIN);
    }
    de->SetObjectiveGradient(value);

    LOG_DBG3(em) << "<l,K'*u-f'> = " << sp << " -> " << value;

    sum += value;
  }
  
  return sum;
}

template double ErsatzMaterial::CalcU1KU2<double>(TransferFunction* tf, StdVector<CFSVector*>& u1, 
                       Application app, StdVector<CFSVector*>& u2, 
                       SurfaceRef* rhs, bool add, double factor);

template double ErsatzMaterial::CalcU1KU2<std::complex<double> >(TransferFunction* tf, StdVector<CFSVector*>& u1, 
                       Application app, StdVector<CFSVector*>& u2, 
                       SurfaceRef* rhs, bool add, double factor);

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


double ErsatzMaterial::CalcConstraint(Condition* constraint, bool derivative, double* grad_out)
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
  
  LOG_TRACE2(em) << "CalcConstraint " << constraint->ToString() 
                  << " derivative=" << derivative << " -> " << result;
  return result;                 
}

double ErsatzMaterial::CalcVolume(bool derivative, Condition* constraint, double* grad_out)
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



double ErsatzMaterial::CalcGreyness(bool derivative, Condition* constraint, double* grad_out)
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


void ErsatzMaterial::SolveStateProblem()
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
      SolveAdjointProblem();
      break;

    default:
      Optimization::SolveStateProblem();
      break;
  }
}

StdVector<IdBcList> ErsatzMaterial::IDBC_2_HDBC()
{
  // IDBC values are homogeneous in the adjoint PDE. Consider elec and mech 
  // store org value as idbc_list per pde
  StdVector<IdBcList> org_idbc(pde.GetSize());
  for(unsigned int p = 0; p < pde.GetSize(); p++) {
    // get the idbc list of the current pde
    IdBcList idbc_list = pde[p]->GetIDBCList();
    for(unsigned int bc = 0; bc < idbc_list.GetSize(); bc++) {
      // store original value -> we have to do a deep copy!
      org_idbc[p].Push_back(shared_ptr<InhomDirichletBc>(new InhomDirichletBc(*(idbc_list[bc]))));
      // make homogeneous values! -> this is only stored in the PDE, we have to call SetBCs() later!
      (*idbc_list[bc]).value = "0.0";
      (*idbc_list[bc]).phase = "0.0";
      LOG_DBG(em) << "Set IDBC to HDBC (" << pde[p]->GetName() << ") -> value " 
                    << (*(org_idbc[p][bc])).value << " -> " << (*(pde[p]->GetIDBCList()[bc])).value;
    }
  }
  
  return org_idbc;
}

void ErsatzMaterial::ResetHDBC(StdVector<IdBcList> org_idbc)
{
  // reset the original IDBC which were IDBC for the adjoint PDE
  for(unsigned int p = 0; p < pde.GetSize(); p++) {
    // get the idbc list of the current pde
    IdBcList idbc_list = pde[p]->GetIDBCList();
    IdBcList org_list  = org_idbc[p];
    assert(idbc_list.GetSize() == org_list.GetSize());
    for(unsigned int bc = 0; bc < idbc_list.GetSize(); bc++) {
      (*idbc_list[bc]).value = (*org_list[bc]).value;
      (*idbc_list[bc]).phase = (*org_list[bc]).phase;
      LOG_DBG(em) << "Reset HDBC (" << pde[p]->GetName() << ") -> " << (*(pde[p]->GetIDBCList()[bc])).value;
    }
    pde[p]->SetBCs();
  }
}  
  
template <class T>
LoadList ErsatzMaterial::SetOutputRHS()
{
  // overwrite the assemble loads with our adjoint rhs loads
  // Note, that in case of inhomogeneous dirichlet values there might be the
  // penalty term on the rhs!
  LoadList org_loads = assemble_->GetLoads();

  assemble_->SetLoads(output_nodes_);

  // assemble the output nodes
  assemble_->AssembleRHSLoads(mech->GetSolveStep()->GetActTime());

  // save this rhs, it has no idbc and penalization set yet
  T* rhs_ptr;
  int rhs_size = assemble_->GetAlgSys()->GetRHSVal(rhs_ptr);
  if(output_vector_ == NULL) output_vector_ = new Vector<T>(rhs_size);
  dynamic_cast<Vector<T> *>(output_vector_)->Replace(rhs_size, rhs_ptr, false);

  LOG_DBG2(em) << "SetAndSolveAdjointRHS: pure output vector w/o idbc: " << output_vector_->ToString(); 
 
  return org_loads;
}


template <class T>
void ErsatzMaterial::CalcAdjointRHS()
{
  // calculate the RHS
  switch(cost->type)
  {
  case OUTPUT:
       // real output set with SetOutputRHS()
       break;
  
  case CONJUGATE_OUTPUT:
  {
    // the correct conjugate_output case is -L * u* -- always complex!
    Vector<complex<double> >& rhs = dynamic_cast<Vector<complex<double> >& >(*output_vector_);
    Vector<complex<double> >& u   = dynamic_cast<Vector<complex<double> >& >(*forward_->raw[MECH]);

    for(unsigned int i = 0; i < rhs.GetSize(); i++)
    {
      // shall be pure real
      assert( ((complex<double>) rhs[i]).imag() == 0);
      // should not have penalization!
      assert( ((complex<double>) rhs[i]).real() < 1e15);
      rhs[i] *= std::conj(u[i]);
    }

    LOG_DBG2(em) << "SetAndSolveAdjointRHS: conjugate output adjoint rhs before solving: " << rhs.ToString();
    assemble_->GetAlgSys()->InitRHS((double*) rhs.GetPointer());
  
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

    break;
  }
  
  default:
    assert(false);
  }
}

template <class T>
void ErsatzMaterial::SetAndSolveAdjointRHS()
{
  // set our own RHS but delete first as Assemble adds
  assemble_->GetAlgSys()->InitRHS();

  // check if we have to set an output RHS
  LoadList org_loads;
  if(cost->type == OUTPUT || cost->type == CONJUGATE_OUTPUT) org_loads = SetOutputRHS<T>();

  // any adjoint PDE has HDBC instead of IDBC
  StdVector<IdBcList> org_idbc = IDBC_2_HDBC();
  
  // calculate the RHS
  CalcAdjointRHS<T>();

  // make IDBC to HDBC -> moves the stuff to OLAS
  for(unsigned int p = 0; p < pde.GetSize(); p++) 
    pde[p]->SetBCs();
  
  // calculate adjoint problem
  assemble_->GetAlgSys()->Solve(GetSolveComment() + "_adjoint");

  ResetHDBC(org_idbc);
  
  // reset the original loads
  if(cost->type == OUTPUT || cost->type == CONJUGATE_OUTPUT) assemble_->SetLoads(org_loads);
}

template <class T>
void ErsatzMaterial::SolveAdjointProblem()
{
  // mech is always present only elec is an option (PiezoSIMP=
  assert(mech != NULL);
  
  // calculate forward problem.
  Optimization::SolveStateProblem();
  
  // store solution element wise for gradient and vector for objective
  // aktually this is redundant as currently the solution is the global one!
  for(unsigned int i=0; i < pde.GetSize(); i++) {
    forward_->ReadSolution(Solution::ELEMENT_VECTORS, pde[i], ToApp(pde[i]));
    forward_->ReadSolution(Solution::RAW_VECTOR, pde[i], ToApp(pde[i]));
  }
  LOG_DBG3(em) << "forward solution: " << forward_->raw[MECH]->ToString(); 
    
  // Set the rhs!
  SetAndSolveAdjointRHS<T>();

  // store the solution in the PDE! This is necessary to read the element vector
  T* ptr;
  int length = assemble_->GetAlgSys()->GetSolutionVal(ptr);
  for(unsigned int i=0; i < pde.GetSize(); i++) {
    pde[i]->SaveSolution(ptr, length);
    // now store own version
    adjoint_->ReadSolution(Solution::ELEMENT_VECTORS, pde[i], ToApp(pde[i]));
  }
  
  // write back the solution s.th. CommitIteraion() makes StoreResults() properly.
  for(unsigned int i=0; i < pde.GetSize(); i++) {
    assert(length == (int) forward_->raw[MECH]->GetSize());
    forward_->WriteSolution( pde[i], ToApp(pde[i]));
  }
}

SinglePDE* ErsatzMaterial::ToPDE(Optimization::Application app, bool throw_exception) const
{
  switch(app)
  {
  case MECH: return pde[0];
  case ELEC: return pde[1];
  default: if(throw_exception) throw Exception("cannot map to pde");
                          else return NULL;
  }
}


Optimization::Application ErsatzMaterial::ToApp(SinglePDE* pde) const
{
  if(pde->GetName() == "electrostatic") return ELEC;
  if(pde->GetName() == "mechanic") return MECH;
  EXCEPTION("pdetype unknown;")
}



ErsatzMaterial::Solution::Solution(ErsatzMaterial* em)
{
  this->em = em;
}

ErsatzMaterial::Solution::~Solution()
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
void ErsatzMaterial::Solution::WriteSolution(StdPDE* pde, Application app)
{
  T* ptr = NULL; 
  raw[app]->GetPointer(ptr);
  assert(ptr != NULL);
  assert(raw[app]->GetSize() != 0);
  pde->SaveSolution(ptr, raw[app]->GetSize());
}


template <class T>
void ErsatzMaterial::Solution::ReadSolution(StorageType st, StdPDE* pde, Application app)
{
  switch(st)
  {
    case ELEMENT_VECTORS:
    {
      // we save the element vectors in elem_vec. Might be empty the first call
      StdVector<CFSVector*>& elem_vec = elem[app];
      int n = em->design->GetNumberOfElements();
      
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
        DesignElement* de = &em->design->data[e];

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
      int size = em->assemble_->GetAlgSys()->GetSolutionVal(ptr);
      
      // check for first call
      if(raw.find(app) == raw.end())
        raw[app] = new Vector<T>(size);
      
      CFSVector* raw_vec = raw[app];

      for(int i = 0; i < size; i++)
        raw_vec->SetEntry(i, ptr[i]);
    }
  }
}

