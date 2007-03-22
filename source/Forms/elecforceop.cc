// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "elecforceop.hh"

#include <string>

#include "Domain/elem.hh"
#include "Domain/grid.hh"
#include "Utils/vector.hh"
#include "Matrix/matrix.hh"
#include "PDE/StdPDE.hh"
#include "Utils/elemstoresol.hh"
#include "Materials/baseMaterial.hh"

namespace CoupledField
{
 
  ElecForceOp::ElecForceOp(Grid * ptGrid,
                           StdPDE * ptPDE,
                           shared_ptr<EqnMap> eqnMap,
                           NodeStoreSol<Double> & sol,
                           UInt dim,
                           std::map<RegionIdType,BaseMaterial*>& matData,
                           shared_ptr<ResultInfo> result,
                           bool isaxi, bool coordUpdate ) 
    : BaseForceOp(ptGrid, ptPDE, eqnMap, sol, dim, 
                  matData, isaxi, coordUpdate)
  {
    ENTER_FCN( "ElecForceOp::ElecForceOp" );

    gradFieldOp_ = new GradientFieldOp<Double>(ptGrid, ptPDE, eqnMap, 
                                               sol, ELEC_POTENTIAL, result,isaxi,
                                               coordUpdate);
    solType_ = ELEC_FORCE_VWP;
    sign_    = 1.0;
  }

  ElecForceOp::~ElecForceOp()
  {
    ENTER_FCN( "ElecForceOp::~ElecForceOp" );

    if (gradFieldOp_) 
      delete gradFieldOp_;
  }


  void ElecForceOp::ComputeField(Vector<Double> & Field, 
                                 const EntityIterator& ent,
                                 const Vector<Double> & lCoord)
  {
    ENTER_FCN( "ElecForceOp::ComputeField" );

    gradFieldOp_->CalcElemGradField(Field, ent, lCoord, 1);

  } 

  Double ElecForceOp::GetMatVal(RegionIdType actRegion)
  {
    ENTER_FCN( "ElecForceOp::GetMatVal" );

    Double epsilon;
    materials_[actRegion]->GetScalar(epsilon,ELEC_PERMITTIVITY,REAL);
    return epsilon;
  } 

}// end of namespace CoupledField
