#include "Optimization/PiezoSIMP.hh"
#include "Optimization/DesignSpace.hh"
#include "Optimization/DesignElement.hh"
#include "Driver/basedriver.hh"
#include "Driver/singleDriver.hh"
#include "Forms/baseForm.hh"
#include "Forms/linearForm.hh"
#include "Forms/linNeumannInt.hh"
#include "Forms/linPressureInt.hh"
#include "Domain/domain.hh"
#include "Domain/bcs.hh"
#include "PDE/mechPDE.hh"
#include "PDE/elecPDE.hh"
#include "Driver/stdSolveStep.hh"
#include "General/exception.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/Logging/cfslog.hh"

#include "boost/lexical_cast.hpp"

using namespace CoupledField; 

DECLARE_LOG(simp)


PiezoSIMP::Transduction::Transduction()
{ 
  isLoad = true;
  value = 0.0;
  found = false;
}

// <mech type="load" value="-1e0"/>
void PiezoSIMP::Transduction::Set(ParamNode* pn)
{
  isLoad = pn->Get("type")->AsString() == "load";
  value  = pn->Get("value")->AsDouble();
}


PiezoSIMP::PiezoSIMP()
{
   elec = dynamic_cast<ElecPDE*>(domain->GetSinglePDE("electrostatic"));

   case_elec_ = cost->type == TRANSDUCTION ? new Solution(this) : NULL;
   case_mech_ = cost->type == TRANSDUCTION ? new Solution(this) : NULL;
   
   // set the stiffness matrices. mechStiffnes is set in SIMP.
   GetElementMatrix(GetForm(elec, elec, "linElecInt"), elecStiffness);
   GetElementMatrix(GetForm(mech, elec, "linPiezoCoupling"), coupledStiffness);    
   coupledStiffness.Transpose(coupledStiffnessTransposed);

   // set the corresponding systems to nonlinear -> mech, mech was done in SIMP
   GetForm(mech, elec, "linPiezoCoupling")->SetSolDependent(true);      
   GetForm(elec, elec, "linElecInt")->SetSolDependent(true);
   
   // validate the transfer functions
   if(design->design.Find(DesignElement::DENSITY) >= 0)
   {
      if(design->GetTransferFunction(DesignElement::DENSITY, MECH, false) == NULL)
        throw Exception("miss transfer function for densitiy and mechanic");
        
      if(design->GetTransferFunction(DesignElement::DENSITY, ELEC, false) == NULL &&
         design->GetTransferFunction(DesignElement::POLARIZATION, ELEC, false) == NULL)
        throw Exception("miss transfer function for densitiy/polarization and electrostatic");
        
      if(design->GetTransferFunction(DesignElement::DENSITY, PIEZO_COUPLING, false) == NULL)
        throw Exception("miss transfer function for densitiy and coupling");
   }
   if(design->design.Find(DesignElement::POLARIZATION) >= 0)
   {
      if(design->GetTransferFunction(DesignElement::POLARIZATION, PIEZO_COUPLING, false) == NULL)
        throw Exception("miss transfer function for polarization and coupling");
   }
   
   // read additional data - transduction in piezo is mandatory if we have this objective
   ParamNode* pn = param->Get("optimization")->Get("SIMP");
   if(cost->type == TRANSDUCTION)
   {
     if(!pn->Has("piezo")) 
       throw Exception("Doing simp with trandsduction requires 'piezo' element in SIMP");
     
     if(!pn->Get("piezo")->Has("transduction"))
       throw Exception("'transduction' element in 'SIMP->piezo' is missing");   

     // set the transductions via default copy constructor 
     mech_.Set(pn->Get("piezo")->Get("transduction")->Get("mech"));
     elec_.Set(pn->Get("piezo")->Get("transduction")->Get("elec"));
   } 
   
   // the storage is optional
   if(pn->Get("piezo")->Has("storage"))
     storage_ = storage.Parse(pn->Get("piezo")->Get("storage")->Get("save"));
   else
     storage_ = COMMIT;  
   
   if(storage_ != COMMIT && cost->type != TRANSDUCTION)
     throw Exception("Please choose storage to 'COMMIT' when not optimizing transduction.");
   
   // we determine the transduction in more detail
   if(cost->type == TRANSDUCTION)
   {
     logFileHeader += "\ttransd_coupl\ttc_simp\ttransd_elec\tte_simp";
   }
} 

PiezoSIMP::~PiezoSIMP()
{
  if(case_elec_ != NULL) { delete case_elec_; case_elec_ = NULL; }
  if(case_mech_ != NULL) { delete case_mech_; case_mech_ = NULL; }
}

double PiezoSIMP::CalcObjective()
{
   design->WriteDesignToExtern(last_evaluation.GetPointer());
   
   switch(cost->type)
   {
      case TRANSDUCTION:        
           cost->SetValue(CalcTransduction());
           break;
          
      default: 
           SIMP::CalcObjective();
   }         
   return cost->GetValue(); 
}


void PiezoSIMP::CalcObjectiveGradient(double* grad_out)
{
  switch(cost->type)
  {
  case OUTPUT:
    // we calculate lambda^T K sol where 
    // sol = displacement and potential and K = K_uu, K_up, K_up^T, K_pp
    // we calcualate the 4 individual vec Mat vec via CalcU1KU2() and sum it up
    for(unsigned int i = 0; i < design->design.GetSize(); i++)
    {
      DesignElement::Type dt = design->design[i];
      TransferFunction*   tf = NULL;

      // reset as we sum up with CalcU1KU2() 
      design->Reset(dt, DesignElement::COST_GRADIENT);

      // we allow NULL for the transfer functions, 
      // then the gradient is 0 what is done via Reset! 
      // lambda_u * K_uu' * u
      tf = design->GetTransferFunction(dt, MECH, false); // we allow NULL
      if(tf != NULL)
        CalcU1KU2(tf, adjoint_->elem[MECH], MECH, forward_->elem[MECH], true, true);
        // CalcU1KU2(tf, MECH, adjoint_, mechStiffness, MECH, forward_, true, true); 
     
      // lambda_u * K_up' * p
      tf = design->GetTransferFunction(dt, PIEZO_COUPLING, false); 
      if(tf != NULL)
        CalcU1KU2(tf, adjoint_->elem[MECH], PIEZO_COUPLING, forward_->elem[ELEC], true, true);
        //CalcU1KU2(tf, MECH, adjoint_, coupledStiffness, ELEC, forward_, true, true); 
      
      // lambda_p * (K_up^T)' * u
      tf = design->GetTransferFunction(dt, PIEZO_COUPLING, false); 
      if(tf != NULL)
        CalcU1KU2(tf, adjoint_->elem[ELEC], PIEZO_COUPLING, forward_->elem[MECH], true, true);
        //CalcU1KU2(tf, ELEC, adjoint_, coupledStiffnessTransposed, MECH, forward_, true, true); 
      
      // lambda_p * K_pp' * p
      tf = design->GetTransferFunction(dt, ELEC, false); 
      if(tf != NULL)
        CalcU1KU2(tf, adjoint_->elem[ELEC], ELEC, forward_->elem[ELEC], true, true);
        //CalcU1KU2(tf, ELEC, adjoint_, elecStiffness, ELEC, forward_, true, true);
    }
    break;
    
  case TRANSDUCTION:        
    CalcTransductionGradient(grad_out); 
    break;

  default:
    SIMP::CalcObjectiveGradient(grad_out);
  }  
  
  if(grad_out != NULL) design->WriteGradientToExtern(grad_out, DesignElement::NO_TYPE, 
                               DesignElement::COST_GRADIENT, DesignElement::SMART);
  
}


template <class T>
void PiezoSIMP::CalcElementKU2(const CFSVector* cfs_in, Application app, DesignElement* de, bool derivative, CFSVector* cfs_out)
{
  const Vector<T>& in  = dynamic_cast<const Vector<T>& >(*cfs_in);
  Vector<T>& out = dynamic_cast<Vector<T>& >(*cfs_out);
  
  // all cases a are a multiplication with an Matrix<double>
  Matrix<double>* stiff = NULL;
  switch(app)
  {
  case ELEC:
    stiff = &elecStiffness;
    break;
    
  case PIEZO_COUPLING:
    // we see on the size of in if we cave to be transposed!
    if(in.GetSize() == coupledStiffness.GetSizeCol())
      stiff = &coupledStiffness;
    else
      stiff = &coupledStiffnessTransposed;
    break;
    
  default: 
    // mech and surface normal matrix are handled in SIMP
    SIMP::CalcElementKU2(cfs_in, app, de, derivative, cfs_out);
    return; // all calculation done there (or assert!)
  }
  
  assert(in.GetSize() == stiff->GetSizeCol());  
  
  TransferFunction* tf = design->GetTransferFunction(de->GetType(), app);
  double factor = derivative ? tf->Derivative(de) : tf->Transform(de);

  out = *stiff * in;
  out *= factor;
  
  LOG_DBG3(simp) << "PiezoSIMP::CalcElementKU2 MECH: " << factor << "*K*u -> " << out.ToString();
}

double PiezoSIMP::CalcTransduction()
{
  double result = 0.0;

  log_coupled_ = log_coupled_simp_ = log_elec_ = log_elec_simp_ = 0.0;
  
  // two temporary vectors as they have different size
  CFSVector* tmp1 = harmonic ? (CFSVector*) new Vector<std::complex<double> > : (CFSVector*) new Vector<double>;
  CFSVector* tmp2 = harmonic ? (CFSVector*) new Vector<std::complex<double> > : (CFSVector*) new Vector<double>;

  StdVector<CFSVector*>& mech_sol_1 = case_elec_->elem[MECH];           
  StdVector<CFSVector*>& elec_sol_1 = case_elec_->elem[ELEC]; 
  StdVector<CFSVector*>& elec_sol_2 = case_mech_->elem[ELEC];
  
  // for the case of special results, here are the indices. -1 is not registered
  int objective_idx = design->GetSpecialResultIndex(DesignElement::DENSITY, DesignElement::OBJECTIVE, DesignElement::NONE);
  int pKu_idx = design->GetSpecialResultIndex(DesignElement::DENSITY, DesignElement::OBJECTIVE, DesignElement::PKU);
  int pKp_idx = design->GetSpecialResultIndex(DesignElement::DENSITY, DesignElement::OBJECTIVE, DesignElement::PKP);

  for(unsigned int i = 0; i < design->GetNumberOfElements(); i++)
  {
    // The mean transduction is given in [1] as 
    // $L_{21} = \phi_2^T K_{u \phi}^T u_1 + \phi_2^T K_{\phi \phi} \phi_1$

    // note, that GetErsatzMaterialFactor() is able to multiply the design and polarization 
    // transfer functions!
    double factor = design->GetErsatzMaterialFactor(i, PIEZO_COUPLING);
    // tmp1 = coupledStiffnessTransposed * mech_sol_1[i];
    coupledStiffnessTransposed.Mult(*mech_sol_1[i], *tmp1);
    double pKu; // = elec_sol_2[i] * tmp1;
    elec_sol_2[i]->Inner(*tmp1, pKu);
    
    log_coupled_ += pKu;
    log_coupled_simp_ += factor * pKu;

    LOG_DBG3(simp) << "calcTransduction elem " << design->data[i].elem->elemNum << " coupling: factors " << pKu << " with " << factor;
    pKu *= factor;

    factor = design->GetErsatzMaterialFactor(i, ELEC);
    // tmp2 = elecStiffness * elec_sol_1[i];
    elecStiffness.Mult(*elec_sol_1[i], *tmp2);
    double pKp; // = elec_sol_2[i] * tmp2;
    elec_sol_2[i]->Inner(*tmp2, pKp);
    
    log_elec_ += pKp;
    log_elec_simp_ += factor * pKp;

    LOG_DBG3(simp) << "calcTransduction elem " << design->data[i].elem->elemNum << " elec: factors " << pKp << " with " << factor;
    pKp *= factor;
    
    // see if we have stuff special registered -> is normally not the case and can hence be slow
    if(pKu_idx != -1) design->data[i].specialResult[pKu_idx] = pKu;
    if(pKp_idx != -1) design->data[i].specialResult[pKp_idx] = pKp;
    if(objective_idx != -1) design->data[i].specialResult[objective_idx] = pKu + pKp;
    
    // sum up the total result for the objective 
    result += pKu + pKp;
  }
  
  delete tmp1;
  delete tmp2;
  return result;
}

void PiezoSIMP::LogFileLine(std::ofstream* out)
{
  Optimization::LogFileLine(out);
  
  if(cost->type == TRANSDUCTION)
  {
    *out << "\t" << log_coupled_ << "\t" << log_coupled_simp_  
         << "\t" << log_elec_ << "\t" << log_elec_simp_;
  }
}


 
void PiezoSIMP::CalcTransductionGradient(double* grad_out)
{
  // The gradient ist given in [1] as 
  // $\frac{\partial L_{21}}{\partial \rho} = - U_1^T K_\rho U_2$
  // where $U_1 = (u_1 \phi_1)^T$ and K is the set of 
  // $K_{uu}, K_{u \phi}, K_{u \phi}^T, K_{\phi \phi}$
  // Note, the simp penaltisation is not given in this notation!!
  //
  // here we split this in four terms of the form
  // $-p_c \rho^{p-1} u_1^T K_{uu} u_2, u_1^T K_{uu}$, here as
  // variabe uKu. p is here used as $\phi$   
  CFSVector* tmp_uKu = harmonic ? (CFSVector*) new Vector<std::complex<double> > : (CFSVector*) new Vector<double>;
  CFSVector* tmp_uKp = harmonic ? (CFSVector*) new Vector<std::complex<double> > : (CFSVector*) new Vector<double>;
  CFSVector* tmp_pKu = harmonic ? (CFSVector*) new Vector<std::complex<double> > : (CFSVector*) new Vector<double>;
  CFSVector* tmp_pKp = harmonic ? (CFSVector*) new Vector<std::complex<double> > : (CFSVector*) new Vector<double>;
    
  StdVector<CFSVector*>& mech_sol_1 = case_elec_->elem[MECH];           
  StdVector<CFSVector*>& mech_sol_2 = case_mech_->elem[MECH];           
  StdVector<CFSVector*>& elec_sol_1 = case_elec_->elem[ELEC]; 
  StdVector<CFSVector*>& elec_sol_2 = case_mech_->elem[ELEC];
  
  
  // the element contributions.    
  double uKu, uKp, pKu, pKp;

  // for the case of special results where we store element contributions, here are the indices. -1 is not registered
  int uKu_idx = design->GetSpecialResultIndex(DesignElement::DENSITY, DesignElement::COST_GRADIENT, DesignElement::UKU);
  int pKu_idx = design->GetSpecialResultIndex(DesignElement::DENSITY, DesignElement::COST_GRADIENT, DesignElement::PKU);
  int uKp_idx = design->GetSpecialResultIndex(DesignElement::DENSITY, DesignElement::COST_GRADIENT, DesignElement::UKP);
  int pKp_idx = design->GetSpecialResultIndex(DesignElement::DENSITY, DesignElement::COST_GRADIENT, DesignElement::PKP);  
    
  // run over all elements and set the values to zero where the design type does not match
  // Do this via runnign over the all design elements in block of elements.
  // Note, that we store all elements for density and then all elements for polarization
  int elements = design->GetNumberOfElements();
  for(unsigned int d = 0; d < design->design.GetSize(); d++)
  {
    DesignElement* de = &design->data[d * elements];
    // The transfer functions. Might be NULL fpr polarization!
    TransferFunction* tf_u  = design->GetTransferFunction(de->GetType(), MECH, false); // false -> might be NULL! 
    TransferFunction* tf_up = design->GetTransferFunction(de->GetType(), PIEZO_COUPLING, false);
    TransferFunction* tf_p  = design->GetTransferFunction(de->GetType(), ELEC, false);       
      
    for(int e = 0; e < elements; e++)
    {
      int idx = d * elements + e;
      de = &design->data[idx];
      // initialization is important for polarization or other dt case
      uKu = uKp = pKu = pKp = 0.0;
      
      #ifdef CHECK_INDEX
        if(de->GetType() != tf_up->GetDesign()) throw Exception("index/type mixed up");
      #endif
      
      if(tf_u != NULL)
      {
        // tmp_uKu = mechStiffness * mech_sol_2[e];
        mechStiffness.Mult(*mech_sol_2[e], *tmp_uKu);
        // uKu = mech_sol_1[e] * tmp_uKu;
        mech_sol_1[e]->Inner(*tmp_uKu, uKu);
        LOG_DBG3(simp) << "calcTransductionGradient elem [" << idx << "]=" << de->elem->elemNum << " uKu = " 
                            << uKu << " * (-1) " << tf_u->Derivative(de) << " = " << uKu * -1.0 * tf_u->Derivative(de);   
        uKu *= -1.0 * tf_u->Derivative(de);
      }  
  
      if(tf_up != NULL)
      {
        // tmp_uKp = coupledStiffness * elec_sol_2[e];
        coupledStiffness.Mult(*elec_sol_2[e], *tmp_uKp);
        // uKp = mech_sol_1[e] * tmp_uKp;
        mech_sol_1[e]->Inner(*tmp_uKp, uKp);
        LOG_DBG3(simp) << "calcTransductionGradient elem [" << idx << "]=" << de->elem->elemNum << " uKp = " 
                            << uKp << " * (-1)" << tf_up->Derivative(de) << " = " << uKp * -1.0 * tf_up->Derivative(de);   
        uKp *= -1.0 * tf_up->Derivative(de);
          
        // tmp_pKu = coupledStiffnessTransposed * mech_sol_2[e];
        coupledStiffnessTransposed.Mult(*mech_sol_2[e], *tmp_pKu);
        // pKu = elec_sol_1[e] * tmp_pKu;
        elec_sol_1[e]->Inner(*tmp_pKu, pKu);
        LOG_DBG3(simp) << "calcTransductionGradient elem [" << idx << "]=" << de->elem->elemNum << " pKu = "
                            << pKu << " * (-1)" << tf_up->Derivative(de) << " = " << pKu * -1.0 * tf_up->Derivative(de);
        pKu *= -1.0 * tf_up->Derivative(de);
      }      
      
      if(tf_p != NULL)
      {
        // tmp_pKp = elecStiffness * elec_sol_2[e];
        elecStiffness.Mult(*elec_sol_2[e], *tmp_pKp);
        // pKp = elec_sol_1[e] * tmp_pKp;
        elec_sol_1[e]->Inner(*tmp_pKp, pKp);
        LOG_DBG3(simp) << "calcTransductionGradient elem [" << idx << "]=" << de->elem->elemNum << " pKp = " 
                            << pKp << " * (-1)" << tf_p->Derivative(de) << " = " << pKp * -1.0 * tf_p->Derivative(de);
        pKp *= -1.0 * tf_p->Derivative(de);
      }
  
      // the gradient for polarization is 0 + uKp + pKu + 0  
      double gradient = uKu + uKp + pKu + pKp;
      
      design->data[idx].SetObjectiveGradient(gradient);
      
      // eventually store the special results
      if(uKu_idx != -1) design->data[idx].specialResult[uKu_idx] = uKu;
      if(pKu_idx != -1) design->data[idx].specialResult[pKu_idx] = pKu;
      if(uKp_idx != -1) design->data[idx].specialResult[uKp_idx] = uKp;
      if(pKp_idx != -1) design->data[idx].specialResult[pKp_idx] = pKp;                  
    }
  }
 
  delete tmp_uKu;
  delete tmp_uKp;
  delete tmp_pKu;
  delete tmp_pKp;  
  
  design->WriteGradientToExtern(grad_out, DesignElement::NO_TYPE, DesignElement::COST_GRADIENT, DesignElement::SMART);
} 

void PiezoSIMP::SolveTransductionSubProblem(Application subProblem)
{
  StdSolveStep* sss = dynamic_cast<StdSolveStep*>(mech->GetSolveStep());
  Assemble* assemble = sss->GetAssemble();
  
  // first check for loads
  LoadList& loads = assemble->GetLoads();
  for(unsigned int i=0; i < loads.GetSize(); i++ ) 
  {
     LoadBc& bc = *loads[i];
            
     // dof=1 is charge, this is zero in the mech case
     if(bc.dof == 1)
     {
       if(!elec_.isLoad) throw Exception("you define a charge load in the PDE but have an integrator in 'transduction'.");
       bc.value = boost::lexical_cast<std::string>(subProblem == ELEC ? elec_.value : 0.0);
       elec_.found = true;
     }

     if(bc.dof == 3)
     {
       if(!mech_.isLoad) throw Exception("you define a force load in the PDE but have an integrator in 'transduction'.");      
       bc.value = boost::lexical_cast<std::string>(subProblem == MECH ? mech_.value : 0.0);
       mech_.found = true;
     }
  }

  std::set<LinearFormContext*>* linForms = assemble->GetLinForms();
  std::set<LinearFormContext*>::iterator it;
  for(it = linForms->begin(); it != linForms->end(); it++)
  { 
      // get integrator
      LinearFormContext& actContext = **it;

      // Check, if lin/non-lin type of Context matches parameter nonLin
      if(actContext.IsNonLin()) continue;

      LinearForm* form = actContext.GetIntegrator();
      
      if(form->GetName() == "LinNeumannInt")
      {
        std::string value = boost::lexical_cast<std::string>(subProblem == ELEC ? elec_.value : 0.0);
        if(elec_.isLoad) throw Exception("you define a surface density in the PDE but have a load in 'transduction'.");
        (dynamic_cast<LinNeumannInt*>(form))->SetAmplitude(value);
        elec_.found = true;
      }
      if(form->GetName() == "PressureLinForm")
      {
        std::string value = boost::lexical_cast<std::string>(subProblem == MECH ? mech_.value : 0.0);
        if(mech_.isLoad) throw Exception("you define a pressure in the PDE but have a load in 'transduction'.");
        (dynamic_cast<PressureLinForm*>(form))->SetValue(value);
        mech_.found = true;
      }      
  }  
  
  
  if(!mech_.found || !elec_.found)
    throw Exception("Define load or pressure/charge density for the mechanic and electrostatic PDEs."); 

  // solve this transduction problem
  Optimization::SolveStateProblem();
  
  // create an element list to gain the iterator in the loop
  ElemList elemList(domain->GetGrid());

  // store the results in our own structure
  if(subProblem == ELEC)
  {
    // ELEC is by definition the first case
    case_elec_->ReadSolution(Solution::ELEMENT_VECTORS, mech, MECH);
    case_elec_->ReadSolution(Solution::ELEMENT_VECTORS, elec, ELEC);
  }
  else
  {
    case_mech_->ReadSolution(Solution::ELEMENT_VECTORS, mech, MECH);
    case_mech_->ReadSolution(Solution::ELEMENT_VECTORS, elec, ELEC);
  }
}

void PiezoSIMP::StoreResults(double step_val)
{
  // called by Optimization::CommitIteration(). We store the last solved problem
  // (which is the mech case) of the current iteration. In SCPIP one iteration
  // is often one problem, in IPOPT there might be much more line searches
  if(storage_ == COMMIT) Optimization::StoreResults(step_val);
  // no COMMIT -> solve explicitly in SolveStateProblem().    
}

void PiezoSIMP::SolveStateProblem()
{
  // when our objective is tranasduction we have to calculate 2 cases
  switch(cost->type)
  {
    // when our objective is output we have to solve also the adjoint problem
    case OUTPUT:
      SolveAdjointProblem(mech, elec);
      break;

    case TRANSDUCTION:
    {
      // the resuls are stored in <mech/elec>_sol_<1/2>
      SolveTransductionSubProblem(ELEC); // by definition first
      
      // format: <iteration>.wwcc where ww is promen within iteration and cc is the case 01 or 02
      // we divide problemsWithinIteration as we have two cases. 
      if(problemWithinIteration > 200) 
        throw Exception("> 100 problems within iteration, set storagage save to commit");
      double step_val = (double) currentIteration + 0.005 * (double) (problemWithinIteration+1) + 0.001;
      
      if(storage_ == ELEC_CASE || storage_ == BOTH_CASES)
         Optimization::StoreResults(step_val); 
         // PiezoSIMP::StoreResults() does nothing when not COMMIT
      
      SolveTransductionSubProblem(MECH);
      if(storage_ == MECH_CASE || storage_ == BOTH_CASES)
         Optimization::StoreResults(step_val + 0.001);
      break;   
    }
    default: // compliance
      SIMP::SolveStateProblem();
      break;
  }
}

  
  
