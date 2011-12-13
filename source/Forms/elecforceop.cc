// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "Domain/resultInfo.hh"
#include "Forms/gradfieldop.hh"
#include "MatVec/vector.hh"
#include "Materials/baseMaterial.hh"
#include "elecforceop.hh"

namespace CoupledField {
class EntityIterator;
class EqnMap;
class Grid;
class StdPDE;
template <typename T> class NodeStoreSol;
}  // namespace CoupledField

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

    gradFieldOp_ = new GradientFieldOp<Double>(ptGrid, ptPDE, eqnMap, sol,
                                               result->fctType, isaxi, coordUpdate);
    solType_ = ELEC_FORCE_VWP;
    sign_    = 1.0;
  }

  ElecForceOp::~ElecForceOp()
  {

    if (gradFieldOp_) 
      delete gradFieldOp_;
  }


  void ElecForceOp::ComputeField(Vector<Double> & Field, 
                                 const EntityIterator& ent,
                                 const Vector<Double> & lCoord)
  {

    gradFieldOp_->CalcElemGradField(Field, ent, lCoord, 1);

  } 

  Double ElecForceOp::GetMatVal(RegionIdType actRegion)
  {

    Double epsilon;
    materials_[actRegion]->GetScalar(epsilon,ELEC_PERMITTIVITY,Global::REAL);
    return epsilon;
  } 

}// end of namespace CoupledField
