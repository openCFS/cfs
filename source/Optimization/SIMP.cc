#include "Optimization/SIMP.hh"
#include "Optimization/DesignSpace.hh"
#include "Optimization/DesignElement.hh"
#include "Domain/domain.hh"
#include "General/exception.hh"
#include "Utils/basenodestoresol.hh"
#include "PDE/SinglePDE.hh"
#include "PDE/mechPDE.hh"
#include "Forms/baseForm.hh"
#include "CoupledPDE/DirectCoupledPDE.hh"
#include "Driver/stdSolveStep.hh"
#include "Utils/StdVector.hh"
#include "Utils/vector.hh"
#include "Utils/result.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/Logging/cfslog.hh"

#include <string>

using namespace CoupledField;

DECLARE_LOG(conditions)
DEFINE_LOG(conditions, "conditions")

DECLARE_LOG(opt)
DEFINE_LOG(opt, "opt")


Enum SIMP::system;

SIMP::SIMP() : Optimization()
{
   // first some plausibility about optimality condition
   if(optimizer == OPTIMALITY_CONDITION)
   {
     if(constraints.GetSize() > 1) 
       throw Exception("optimality condition is only possible with volume constraint");
   }

   // base for bisection search
   lambda_ = 1.0e3;

   ParamNode* pn = param->Get("optimization")->Get("SIMP");

   system_ = (System) system.Parse(pn->Get("system")->AsString());
   
   // region stuff
   regionId = domain->GetGrid()->RegionNameToId(pn->Get("region")->AsString());
   
   // find the right PDE which is the mechanic from the piezo
   mech = dynamic_cast<MechPDE*>(domain->GetSinglePDE("mechanic"));
   
   // manipulate the PDE (the Integrator of its bilinear form) such that 
   // it is set up with the new densitiy values in the optimization loops
   mech->getPDE_assemble()->GetBiLinForm(regionId, mech, mech)->GetIntegrator()->SetSolDependent(true);
   
  
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
   SetElementStiffness(mech, mech, mechStiffness);

   // in case of optimality condition we log the number of lambda iterations
   if(optimizer == OPTIMALITY_CONDITION) logFileHeader += "\tlambda\tlambda iterations";
   
   // check the average and standard deviation
   for(unsigned int i = 0; i < design->design.GetSize(); i++)
   {
     std::string str = DesignElement::type.ToString(design->design[i]);
     logFileHeader += "\tx_avg_" + str + "\tx_std_dev_" + str;
   }
   
   // make a copy of the old iteration to calculate the move
   last_iteration.Resize(design->data.GetSize());
   
   // note the difference between function evaluations (line search) and iterations!   
   last_evaluation.Resize(design->data.GetSize());
   
   // this is needed for optimality condition AND move calculation!
   design_tmp_.Resize(design->data.GetSize());
}

SIMP::~SIMP()
{
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
    
double SIMP::CalcObjective()
{
  design->WriteDesignToExtern(last_evaluation.GetPointer());
  
  switch(cost->type)
  {
     case COMPLIANCE:
          cost->SetValue(CalcCompliance());
          break;

     default: throw Exception("objective no handled");
  }         
       
  return cost->GetValue();
}

void SIMP::SetElementStiffness(StdPDE* pde1, StdPDE* pde2, Matrix<double>& out)
{
    // we manipulate the first element, set the density to 1 and calc the element stiffness
    // matrix. This works only for isotropic material.

    // We have to do this because CalcElementMatrix asks the density via domain and
    // multiplies it with it's matrix!
    
    DesignElement* element = &design->data[0]; 
    double org_density = element->GetDesign(DesignElement::PLAIN);
    element->SetDesign(1.0);
    
    // temporarily disables the transfer functions -> set to IDENTITY
    design->DisableTransferFunctions();
    
    // get the bilinear form from the assemble class, there is the integrator.
    BaseForm* bf = mech->getPDE_assemble()->GetBiLinForm(regionId, pde1, pde2)->GetIntegrator();

    // create an element list to gain the iterator in the loop
    ElemList elemList(domain->GetGrid());

    // we need an iterator to calc the element matrix
    elemList.SetElement(element->elem);
    const EntityIterator& it = elemList.GetIterator();

    // get our element stiffness matrix -> it could be saved over the iteration loops!  
    bf->CalcElementMatrix(out,const_cast<EntityIterator&>(it), const_cast<EntityIterator&>(it));

    element->SetDesign(org_density);
    
    // enable again our tranfer functions
    design->EnableTransferFunctions();
}


double SIMP::CalcCompliance()
{
  // The compliance is sum over elements \rho^p u K u 
  // where u is the displacement and K the local stiffness matrix=elasticity tensor
  
  // there is a feature in MechPDE (mech->CalcEnergy<double>(regionId)) which is slower
  // but should be used when using non-uniform cells
  return CalcUKU(false); 
}


void SIMP::CalcObjectiveGradient(double* grad_out)
{
  if(cost->type != COMPLIANCE) throw Exception("cost function not implemented", __FILE__, __LINE__);
  CalcComplianceGradient(grad_out);
}

void SIMP::CalcComplianceGradient(double* grad_out)
{
  // calculate the compliance which is according to 
  // "A 99 line topology optimization code written in Matlab"; O.Sigmund, 2001
  // -> dc/dx_e = -p * x_e ^(p-1) u_e^T k_0 u_e
  // Here we use a constant element stiffness matrix
  CalcUKU(true);

  // if we output for an externel optimizer, we use the filtered element if filtering is
  // enabled 
  for(unsigned int i = 0; grad_out != NULL && i < design->data.GetSize(); i++)
  {
    DesignElement* de = &design->data[i]; 

    // the gradient is zero for not DENSITY
    bool relevant = de->GetType() == DesignElement::DENSITY;
    double grad = relevant ? de->GetObjectiveGradient(DesignElement::SMART) : 0.0;
    grad_out[i] = grad;
  }    
}

double SIMP::CalcUKU(bool derivative)
{
  Vector<double> vec;
  Vector<double> tmp;
  double sum = 0.0;
 
  // create an element list to gain the iterator in the loop
  ElemList elemList(domain->GetGrid());
  
  // Get the transfer function
  TransferFunction* tf = design->GetTransferFunction(DesignElement::DENSITY, MECH, true);
 
  int elements = design->GetNumberOfElements();
  int base = design->FindDesign(DesignElement::DENSITY) * elements; 
 
  for(int i = base; i < base+elements; i++)    
  {
     DesignElement* de = &design->data[i];

     elemList.SetElement(de->elem);
     const EntityIterator& it = elemList.GetIterator(); 
 
     // read the displacements of the nodes of our element
     mech->GetSolVecOfElement(vec, it, mech->GetResultInfo(MECH_DISPLACEMENT)); 

     //  u_e^T k_0 u_e 
     tmp = mechStiffness * vec;  
     double value = tmp * vec; // here without transfer function

     if(derivative) 
     {
       value *= -1.0 * tf->Derivative(de);
       de->SetObjectiveGradient(value);
     } 
     else
     {
       // We don't need GetErsatzMaterialFactor() as we do only with MECH 
       value *= tf->Transform(de);
     }   
      
     sum += value;
  }
  
  return sum;
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
  
  LOG_TRACE2(opt) << "CalcConstraint " << constraint->ToString() 
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
 
  return sum / (double) counter;
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



void SIMP::LogFileLine(std::ofstream* out)
{
  Optimization::LogFileLine(out);
  
  if(optimizer == OPTIMALITY_CONDITION) 
    *out << "\t" << lambda_ << "\t" << lambda_iterations_;

  for(unsigned int i = 0; i < design->design.GetSize(); i++)
  {
    int n = design->GetNumberOfElements();
    double* ptr = last_iteration.GetPointer() + i * n;
    *out << "\t" << Average(ptr, n) << "\t" << StandardDeviation(ptr, n);
  }
}


void SIMP::CalcNextDensity()
{
  // we store the current densities in the temp variable. Otherwise we cannot
  // find the proper lambda
  design->WriteDesignToExtern(design_tmp_.GetPointer());
  
  // do a standard bisection - we restart from the last try. This 
  // does not work if inital lambda_ is too low or if a real next
  // lamba_ is twice or more than the old -> see stuff in code  
  double lower = 0;
  double upper = lambda_ * 2;

  // we count to handle the problem of above
  lambda_iterations_ = 0;
  int count = 0;
  double err = 25; // initial value so whe run the loop at least once
   
  while(std::abs(err) > 1e-3)
  {
    // for first loop lambda has the same value as on function entry
    lambda_iterations_++; 
    lambda_ = 0.5 * (upper + lower);

    // restore original density from temp so we always start the calculation 
    // on the same base but with different lambda
    design->ReadDesignFromExtern(design_tmp_.GetPointer());

    OptimalityConditionStep(lambda_);
    double vol = CalcConstraint();
    err = GetConstraint(Condition::VOLUME).value - vol;
    
    // std::cout << "lower = " << lower << " upper = " << upper << " lambda = " << lambda_ << " constraint = " << vol << " err = " << err << std::endl;     

    if(err > 0) upper = lambda_; // = center
           else lower = lambda_;
     
    // check if we have a problem, that our inital upper is too large
    if(++count > 50)
    {
      std::cout << "Problems in SIMP::CalcNextDensityenlarge with bisection -> enlarge upper limit" << std::endl;
      upper *= 1e3;
      count = 0;
    }
   }
}


void SIMP::OptimalityConditionStep(double lambda)
{
   // we assume DensityElement.objective_gradient to be set
   // we assume DensityElement.constraint_gradient to be set
  
   //double next_sum = 0.0;
   //double next_real = 0.0;

   Condition& condition = GetConstraint(Condition::VOLUME);

   for(unsigned int i = 0; i < design->data.GetSize(); i++)    
   {
     DesignElement* de = &design->data[i];
     // rho_e is the old rho
     double rho_e = de->GetDesign(DesignElement::PLAIN);   
    
     // if filter is enabled we use the filtered value otherwise the plain one
     double b_e = -1.0 * de->GetObjectiveGradient(DesignElement::SMART);

     if(isnan(b_e))
      throw Exception("b_e is nan");

     b_e /= (lambda * de->GetConstraintGradient(&condition)); 
     
     // maximization is tricky
     // if(cost->task == MAXIMIZE) b_e = 2-b_e;
     
     // next is density times b_e which is compared with box constraints and move limit
     double next = rho_e * std::pow(b_e, oc_damping);        
                     
     double lower = std::max(de->GetLowerBound(), rho_e - move_limit);
     double upper = std::min(de->GetUpperBound(), rho_e + move_limit);            

     //next_real += next;

     de->SetDesign(next);  
     if(next <= lower) de->SetDesign(lower);
     if(upper <= next) de->SetDesign(upper);
     
     //next_sum += data[i].GetDesign(DesignElement::PLAIN);
   }
   
   // std::cout << " real next would be " << (next_real / (double) data.GetSize()) << " limited next = " << (next_sum / (double) data.GetSize()) << std::endl; 
}

void SIMP::GetElasticityTensor(Matrix<double>& matrix_out)
{
    // the elasticity tenser of the mech PDE
    std::map<RegionIdType, BaseMaterial*> materials = mech->getPDEMaterialData();
    
    BaseMaterial* bm = materials[regionId];
    if(bm == NULL) 
       throw Exception("did not find material for optimized region", __FILE__, __LINE__); 
    bm->GetTensor(matrix_out,MECH_STIFFNESS_TENSOR,REAL,FULL);
}


void SIMP::ExecuteOptimizationStep()
{
    // calc gradients to store the results in data[element]...
    CalcObjectiveGradient(NULL);
    CalcConstraintGradient(NULL, NULL);
    // do a SIMP Optimality Condition step
    CalcNextDensity();
    // calc the objective for the logging in CommitIteration(),
    // for the optimality condition it is not required.
    CalcObjective();
}


