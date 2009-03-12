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

   for(unsigned int r = 0; r < regionIds.GetSize(); r++){
     GetElementMatrix(GetForm(regionIds[r], elec, elec, "linElecInt"), elecStiffness_map[regionIds[r]], NULL, DesignElement::NO_DERIVATIVE, -1.0);
     GetElementMatrix(GetForm(regionIds[r], mech, elec, "linPiezoCoupling"), coupledStiffness_map[regionIds[r]]);
     coupledStiffness_map[regionIds[r]].Transpose(coupledStiffnessTransposed_map[regionIds[r]]);
     GetForm(regionIds[r], mech, elec, "linPiezoCoupling")->SetSolDependent(true);
     GetForm(regionIds[r], elec, elec, "linElecInt")->SetSolDependent(true);
   }
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


void PiezoSIMP::CalcObjectiveGradient(Excitation& excite)
{
  switch(cost->type)
  {
  case OUTPUT:
  case DYNAMIC_OUTPUT:
  case GLOBAL_DYNAMIC_COMPLIANCE:
  {
    unsigned int idx = excite.index;

    double factor = excite.normalized_weight;
    if(cost->FactorOmegaOmega()) factor *= excite.GetOmegaOmega();
    LOG_DBG2(simp) << "PS::CalcObjectiveGradient(idx=" << excite.index << ") factor = "
                   << excite.normalized_weight << " * "
                   << (cost->FactorOmegaOmega() ? excite.GetOmegaOmega() : 1.0) << " -> " << factor;

    // we see at mechRHS and elecRHS if we have such an excitation
    SurfaceRef* mech_rhs = mechRHS.valid ? &mechRHS : NULL;
    SurfaceRef* elec_rhs = elecRHS.valid ? &elecRHS : NULL;

    LOG_TRACE2(simp) << "PS:CalcOutputGradient: output -> sensitive rhs's: "
                     << mechRHS.ToString(1) << ", " << elecRHS.ToString(1);

    // we calculate lambda^T K sol where
    // sol = displacement and potential and K = K_uu, K_up, K_up^T, K_pp
    // we calcualate the 4 individual vec Mat vec via CalcU1KU2() and sum it up
    for(unsigned int i = 0; i < design->design.GetSize(); i++)
    {
      DesignElement::Type dt = design->design[i];
      TransferFunction*   tf = NULL;

      // we allow NULL for the transfer functions,
      // then the gradient is 0 what is done via Reset!
      // lambda_u * K_uu' * u
      tf = design->GetTransferFunction(dt, MECH, false); // we allow NULL
      if(tf != NULL)
        CalcU1KU2(tf, adjoint.data[idx]->elem[MECH], MECH, forward.data[idx]->elem[MECH], mech_rhs, true, factor);

      // lambda_u * K_up' * p
      tf = design->GetTransferFunction(dt, PIEZO_COUPLING, false);
      if(tf != NULL)
        CalcU1KU2(tf, adjoint.data[idx]->elem[MECH], PIEZO_COUPLING, forward.data[idx]->elem[ELEC], mech_rhs, true, factor);

      // lambda_p * (K_up^T)' * u
      tf = design->GetTransferFunction(dt, PIEZO_COUPLING, false);
      if(tf != NULL)
        CalcU1KU2(tf, adjoint.data[idx]->elem[ELEC], PIEZO_COUPLING, forward.data[idx]->elem[MECH], elec_rhs, true, factor);

      // lambda_p * K_pp' * p
      tf = design->GetTransferFunction(dt, ELEC, false);
      if(tf != NULL)
        CalcU1KU2(tf, adjoint.data[idx]->elem[ELEC], ELEC, forward.data[idx]->elem[ELEC], elec_rhs, true, factor);
    }
  }
  break;

  default:
    SIMP::CalcObjectiveGradient(excite);
  }
}

template <class T>
void PiezoSIMP::SetElementK(DesignElement* de, Application app, DenseMatrix* mat_out)
{
  TransferFunction* tf = design->GetTransferFunction(de->GetType(), app);
  double factor = tf->Derivative(de);

  Matrix<T>& out = dynamic_cast<Matrix<T>& >(*mat_out);

  switch(app)
  {
  case ELEC:
    Assign(out, ElecStiffness(de->elem), factor);
    break;

  case PIEZO_COUPLING:
  {
    const Matrix<double>& coupledStiffness = CoupledStiffness(de->elem);
    // we see on the size of in if we cave to be transposed!
    if(out.GetNumCols() == coupledStiffness.GetNumCols())
      Assign(out, coupledStiffness, factor);
    else
      Assign(out, CoupledStiffnessTransposed(de->elem), factor);
    break;
  }

  default:
    // mech and surface normal matrix are handled in SIMP
    SIMP::SetElementK(de, app, mat_out);
    return; // all calculation done there (or assert!)
  }

  LOG_DBG2(simp) << "PiezoSIMP::SetElementK elem: " << de->elem->elemNum << " app: " << application.ToString(app) << " factor: " << factor;
}


const Matrix<double>& PiezoSIMP::ElecStiffness(Elem* elem){
  return elecStiffness_map[elem->regionId];
}

const Matrix<double>& PiezoSIMP::CoupledStiffness(Elem* elem){
  return coupledStiffness_map[elem->regionId];
}

const Matrix<double>& PiezoSIMP::CoupledStiffnessTransposed(Elem* elem){
  return coupledStiffnessTransposed_map[elem->regionId];
}



