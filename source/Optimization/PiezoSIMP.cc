#include "Optimization/PiezoSIMP.hh"
#include "Optimization/SIMP.hh"
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

using std::complex;

DECLARE_LOG(simp)


PiezoSIMP::PiezoSIMP()
{
   elec = dynamic_cast<ElecPDE*>(domain->GetSinglePDE("electrostatic"));
   assert(pde.GetSize() == 1);
   pde.Push_back(elec);
   
   // set the stiffness matrices. mechStiffnes is set in SIMP.
   GetElementMatrix(GetForm(elec, elec, "linElecInt"), elecStiffness);
   GetElementMatrix(GetForm(mech, elec, "linPiezoCoupling"), coupledStiffness);    
   coupledStiffness.Transpose(coupledStiffnessTransposed);

   // set the corresponding systems to nonlinear -> mech, mech was done in SIMP
   GetForm(mech, elec, "linPiezoCoupling")->SetSolDependent(true);      
   GetForm(elec, elec, "linElecInt")->SetSolDependent(true);
   // The linear forms (pressure, charge density) are set in SoluctionRef::Init()
   
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
} 

PiezoSIMP::~PiezoSIMP()
{
  if(elecRHS.vec != NULL) { delete elecRHS.vec;   elecRHS.vec = NULL; }
}

void PiezoSIMP::PostInit()
{
  if(harmonic) elecRHS.Init<complex<double> >(design, CHARGE_DENSITY); // mechRHS in SIMP!
          else elecRHS.Init<double>(design, CHARGE_DENSITY);

  SIMP::PostInit();
}


void PiezoSIMP::CalcObjectiveGradient(double* grad_out)
{
  switch(cost->type)
  {
  case OUTPUT:
  {
    // we see at mechRHS and elecRHS if we have such an excitation
    SurfaceRef* mech_rhs = mechRHS.valid ? &mechRHS : NULL;
    SurfaceRef* elec_rhs = elecRHS.valid ? &elecRHS : NULL;
    
    LOG_TRACE2(simp) << "CalcObjectiveGradient: output -> sensitive rhs's: "
                     << mechRHS.ToString(1) << ", " << elecRHS.ToString(1);

    
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
        CalcU1KU2(tf, adjoint_->elem[MECH], MECH, forward_->elem[MECH], mech_rhs, true);
     
      // lambda_u * K_up' * p
      tf = design->GetTransferFunction(dt, PIEZO_COUPLING, false); 
      if(tf != NULL)
        CalcU1KU2(tf, adjoint_->elem[MECH], PIEZO_COUPLING, forward_->elem[ELEC], mech_rhs, true);
      
      // lambda_p * (K_up^T)' * u
      tf = design->GetTransferFunction(dt, PIEZO_COUPLING, false); 
      if(tf != NULL)
        CalcU1KU2(tf, adjoint_->elem[ELEC], PIEZO_COUPLING, forward_->elem[MECH], elec_rhs, true);
      
      // lambda_p * K_pp' * p
      tf = design->GetTransferFunction(dt, ELEC, false); 
      if(tf != NULL)
        CalcU1KU2(tf, adjoint_->elem[ELEC], ELEC, forward_->elem[ELEC], elec_rhs, true);
    }
    break;
  }

  default:
    SIMP::CalcObjectiveGradient(grad_out);
  }  
  
  if(grad_out != NULL) design->WriteGradientToExtern(grad_out, DesignElement::NO_TYPE, 
                               DesignElement::COST_GRADIENT, DesignElement::SMART);
  
}


template <class T>
void PiezoSIMP::SetElementK(DesignElement* de, Application app, CFSMatrix* mat_out)
{
  TransferFunction* tf = design->GetTransferFunction(de->GetType(), app);
  double factor = tf->Derivative(de);
  
  Matrix<T>& out = dynamic_cast<Matrix<T>& >(*mat_out);
  
  switch(app)
  {
  case ELEC:
    Assign(out, elecStiffness, factor);
    break;
    
  case PIEZO_COUPLING:
    // we see on the size of in if we cave to be transposed!
    if(out.GetSizeCol() == coupledStiffness.GetSizeCol())
      Assign(out, coupledStiffness, factor);
    else
      Assign(out, coupledStiffnessTransposed, factor);      
    break;
    
  default: 
    // mech and surface normal matrix are handled in SIMP
    SIMP::SetElementK(de, app, mat_out);
    return; // all calculation done there (or assert!)
  }

  LOG_DBG2(simp) << "PiezoSIMP::SetElementK elem: " << de->elem->elemNum << " app: " << application.ToString(app) << " factor: " << factor;
}



  
  
