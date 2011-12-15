// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <assert.h>
#include <stddef.h>
#include <string>

#include "CoupledPDE/BasePairCoupling.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "Domain/domain.hh"
#include "Domain/elem.hh"
#include "Domain/entityList.hh"
#include "Domain/resultInfo.hh"
#include "Driver/assemble.hh"
#include "Driver/formsContext.hh"
#include "Forms/PiezoStressStrain.hh"
#include "Forms/baseForm.hh"
#include "Forms/gradfieldop.hh"
#include "Forms/linPiezoCoupling.hh"
#include "Forms/nLinPiezoCoupling.hh"
#include "PDE/SinglePDE.hh"
#include "PDE/StdPDE.hh"
#include "Utils/nodestoresol.hh" // IWYU pragma: keep

namespace CoupledField {
class BaseNodeStoreSol;
class Grid;
}  // namespace CoupledField

DECLARE_LOG(ss)

namespace CoupledField
{

  // =============================================================================
  // base class for mechanical stress computation
  // =============================================================================

  // base class for calculation of mechanical stresses
  template <class TYPE>
  PiezoStressStrain<TYPE>::PiezoStressStrain(BaseMaterial* mechMat, SubTensorType type, BasePairCoupling* bpc)
    : MechStressStrain<TYPE>(mechMat, type)
  {
    this->bpc_ = bpc;
    this->name_ = "PiezoStressStrain";

    StdPDE* mech = bpc->GetPde1();
    StdPDE* elec = bpc->GetPde2();

    assert(mech != NULL && mech->GetName() == "mechanic");
    assert(elec != NULL && elec->GetName() == "electrostatic");

    BaseNodeStoreSol * elec_bnss = elec->getPDESolution();
    NodeStoreSol<TYPE> * elec_nss = dynamic_cast<NodeStoreSol<TYPE>*>(elec_bnss);

    Grid* grid = domain->GetGrid();

    // Determines gradient of electric potential, i.e. E=\grad \phi
    gradOp = new GradientFieldOp<TYPE>(grid, elec, bpc->GetEqnMap2(), *elec_nss, bpc_->GetResultInfos2()[0]->fctType, mech->GetIsaxi());


  }

  template <class TYPE>
  PiezoStressStrain<TYPE>::~PiezoStressStrain()
  {
    delete gradOp; gradOp = NULL;
  }



  template <class TYPE>
  void PiezoStressStrain<TYPE>::CalcStressVec(Vector<TYPE>& stressVec, UInt ip, EntityIterator& ent)
  {
    MechStressStrain<TYPE>::CalcStressVec(stressVec, ip, ent); // base class call

    // this is clearly slow, but fast enough up to the FE-Space changing everything!
    StdPDE* mech = bpc_->GetPde1();
    StdPDE* elec = bpc_->GetPde2();

    RegionIdType reg = ent.GetElem()->regionId;

    // TODO check for imaginary coupling as in harmAxiQC
    BiLinFormContext* bc = mech->getPDE_assemble()->GetBiLinForm(reg, mech, elec, "linPiezoCoupling", true, Global::REAL);
    if(bc != NULL)
    {
      linPiezoCoupling* couplingForm = dynamic_cast<linPiezoCoupling*>(bc->GetIntegrator());

      // proper transposed coupling matrix with applied pseudo density
      couplingForm->calcDMat(piezoCouplingMatT, ent.GetElem());

    }
    else
    {
      bc = mech->getPDE_assemble()->GetBiLinForm(reg, mech, elec, "nLinPiezoCoupling", false); // now needs to be nonlin
      nLinPiezoCoupling* couplingForm = dynamic_cast<nLinPiezoCoupling*>(bc->GetIntegrator());

      couplingForm->calcDMat(piezoCouplingMatT, ip, this->ptCoord_);
    }
    // determine the electric field at integration point
    gradOp->CalcElemGradField(grad, ent, this->intPoint_, 1.0);
    tmpCouplStress = piezoCouplingMatT * grad;

    stressVec -= tmpCouplStress;
  }


#ifdef __GNUC__
  // Explicite template instantiation
  template class PiezoStressStrain<double>;
  template class PiezoStressStrain<Complex>;
#endif

} // end of namespace
