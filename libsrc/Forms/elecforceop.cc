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
                           NodeEQN * ptEQN,
                           NodeStoreSol<Double> & sol,
                           UInt dim,
                           StdVector<BaseMaterial*> &matData,
                           StdVector<RegionIdType> & allSubdoms,
                           Boolean isaxi) 
    : BaseForceOp(ptGrid, ptPDE, ptEQN, sol, dim, matData, allSubdoms, isaxi)
  {
    ENTER_FCN( "ElecForceOp::ElecForceOp" );

    gradFieldOp_ = new GradientFieldOp<Double>(ptGrid, ptPDE, ptEQN, 
                                               sol, ELEC_POTENTIAL, isaxi);
    solType_ = ELEC_FORCE_VWP;
    sign_    = 1.0;
  }

  ElecForceOp::~ElecForceOp()
  {
    ENTER_FCN( "ElecForceOp::~ElecForceOp" );

    if (gradFieldOp_) 
      delete gradFieldOp_;
  }


  void ElecForceOp::ComputeField(Vector<Double> & Field, const Elem * ptElement,
                                 const Vector<Double> & lCoord)
  {
    ENTER_FCN( "ElecForceOp::ComputeField" );

    gradFieldOp_->CalcElemGradField(Field, ptElement, lCoord, 1);

  } 

  Double ElecForceOp::GetMatVal(UInt actSD)
  {
    ENTER_FCN( "ElecForceOp::GetMatVal" );

    Double epsilon;
    materialData_[actSD]->GetScalar(epsilon,ELEC_PERMITTIVITY,REAL);
    return epsilon;
  } 

}// end of namespace CoupledField
