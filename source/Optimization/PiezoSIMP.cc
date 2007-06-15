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

// declare class specific logging stream
DECLARE_LOG(piezoSimp)
DEFINE_LOG(piezoSimp, "piezoSimp")


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

   // here we store the results of the two mean transduction cases      
   mech_sol_1.Resize(design->GetNumberOfElements());
   mech_sol_2.Resize(design->GetNumberOfElements());
   elec_sol_1.Resize(design->GetNumberOfElements());
   elec_sol_2.Resize(design->GetNumberOfElements());
   
   // set the stiffness matrices. mechStiffnes is set in SIMP                    
   SetElementStiffness(elec, elec, elecStiffness);
   SetElementStiffness(mech, elec, coupledStiffness);    
   coupledStiffness.Transpose(coupledStiffnessTransposed);

   // set the corresponding systems to nonlinear -> mech, mech was done in SIMP
   mech->getPDE_assemble()->GetBiLinForm(regionId, mech, elec)->GetIntegrator()->SetSolDependent(true);      
   elec->getPDE_assemble()->GetBiLinForm(regionId, elec, elec)->GetIntegrator()->SetSolDependent(true);
   
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
     storage_ = (Storage) storage.Parse(pn->Get("piezo")->Get("storage")->Get("save")->AsString());
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

double PiezoSIMP::CalcObjective()
{
   design->WriteDesignToExtern(last_evaluation.GetPointer());
   
   switch(cost->type)
   {
      case COMPLIANCE:
           cost->SetValue(CalcCompliance());
           break;
           
      case TRANSDUCTION:        
           cost->SetValue(CalcTransduction());
           break;
          
      default: throw Exception("objective no handled");
   }         
   return cost->GetValue(); 
}


void PiezoSIMP::CalcObjectiveGradient(double* grad_out)
{
   switch(cost->type)
   {
      case COMPLIANCE:
           CalcComplianceGradient(grad_out);
           break;
           
      case TRANSDUCTION:        
           CalcTransductionGradient(grad_out); 
           break;
          
      default: throw Exception("objective no handled");
   }         
}

double PiezoSIMP::CalcTransduction()
{
  double result = 0.0;
  
  log_coupled_ = log_coupled_simp_ = log_elec_ = log_elec_simp_ = 0.0;
  
  // two temporary vectors as they have different size
  Vector<double> tmp1;
  Vector<double> tmp2;

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
    tmp1 = coupledStiffnessTransposed * mech_sol_1[i];
    double pKu = elec_sol_2[i] * tmp1;
    
    log_coupled_ += pKu;
    log_coupled_simp_ += factor * pKu;

    LOG_DBG3(piezoSimp) << "calcTransduction elem " << design->data[i].elem->elemNum << " coupling: factors " << pKu << " with " << factor;
    pKu *= factor;

    factor = design->GetErsatzMaterialFactor(i, ELEC);
    tmp2 = elecStiffness * elec_sol_1[i];
    double pKp = elec_sol_2[i] * tmp2;
    
    log_elec_ += pKp;
    log_elec_simp_ += factor * pKp;

    LOG_DBG3(piezoSimp) << "calcTransduction elem " << design->data[i].elem->elemNum << " elec: factors " << pKp << " with " << factor;
    pKp *= factor;
    
    // see if we have stuff special registered -> is normally not the case and can hence be slow
    if(pKu_idx != -1) design->data[i].specialResult[pKu_idx] = pKu;
    if(pKp_idx != -1) design->data[i].specialResult[pKp_idx] = pKp;
    if(objective_idx != -1) design->data[i].specialResult[objective_idx] = pKu + pKp;
    
    // sum up the total result for the objective 
    result += pKu + pKp;
  }
  
  return result;
}

void PiezoSIMP::LogFileLine(std::ofstream* out)
{
  SIMP::LogFileLine(out);
  
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
  Vector<double> tmp_uKu;
  Vector<double> tmp_uKp;
  Vector<double> tmp_pKu;
  Vector<double> tmp_pKp;
    
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
        tmp_uKu = mechStiffness * mech_sol_2[e];
        uKu = mech_sol_1[e] * tmp_uKu;
        LOG_DBG3(piezoSimp) << "calcTransductionGradient elem [" << idx << "]=" << de->elem->elemNum << " uKu = " 
                            << uKu << " * (-1) " << tf_u->Derivative(de) << " = " << uKu * -1.0 * tf_u->Derivative(de);   
        uKu *= -1.0 * tf_u->Derivative(de);
      }  
  
      if(tf_up != NULL)
      {
        tmp_uKp = coupledStiffness * elec_sol_2[e];
        uKp = mech_sol_1[e] * tmp_uKp;
        LOG_DBG3(piezoSimp) << "calcTransductionGradient elem [" << idx << "]=" << de->elem->elemNum << " uKp = " 
                            << uKp << " * (-1)" << tf_up->Derivative(de) << " = " << uKp * -1.0 * tf_up->Derivative(de);   
        uKp *= -1.0 * tf_up->Derivative(de);
          
        tmp_pKu = coupledStiffnessTransposed * mech_sol_2[e];
        pKu = elec_sol_1[e] * tmp_pKu;
        LOG_DBG3(piezoSimp) << "calcTransductionGradient elem [" << idx << "]=" << de->elem->elemNum << " pKu = "
                            << pKu << " * (-1)" << tf_up->Derivative(de) << " = " << pKu * -1.0 * tf_up->Derivative(de);
        pKu *= -1.0 * tf_up->Derivative(de);
      }      
      
      if(tf_p != NULL)
      {
        tmp_pKp = elecStiffness * elec_sol_2[e];
        pKp = elec_sol_1[e] * tmp_pKp;
        LOG_DBG3(piezoSimp) << "calcTransductionGradient elem [" << idx << "]=" << de->elem->elemNum << " pKp = " 
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
 
    // give the output for the external optimizer most probably filtered   
  for(unsigned int i = 0; grad_out != NULL && i < design->data.GetSize(); i++)
  {
    // see the documentation of GetObjectiveGradient() to understand design times gradient!
    grad_out[i] = design->data[i].GetObjectiveGradient(DesignElement::SMART);
  }    
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
  for(unsigned int i = 0; i < design->GetNumberOfElements(); i++)
  {
    DesignElement* de = &design->data[i];
 
    elemList.SetElement(de->elem);
    const EntityIterator& it = elemList.GetIterator();

    if(subProblem == ELEC)
    {
      // ELEC is by definition the first case 
      mech->GetSolVecOfElement(mech_sol_1[i], it, mech->GetResultInfo(MECH_DISPLACEMENT));
      elec->GetSolVecOfElement(elec_sol_1[i], it, mech->GetResultInfo(ELEC_POTENTIAL));
    }
    else
    {
      mech->GetSolVecOfElement(mech_sol_2[i], it, mech->GetResultInfo(MECH_DISPLACEMENT));
      elec->GetSolVecOfElement(elec_sol_2[i], it, mech->GetResultInfo(ELEC_POTENTIAL));
    }
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
      Optimization::SolveStateProblem();
      break;
  }
}

  
  
