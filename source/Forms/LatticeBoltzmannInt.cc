#include <assert.h>
#include <string>

#include "Forms/LatticeBoltzmannInt.hh"
#include "MatVec/matrix.hh"
#include "Materials/baseMaterial.hh"
#include "Domain/elem.hh"
#include "Domain/entityList.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"

DECLARE_LOG(lbm_int)
DEFINE_LOG(lbm_int, "lbm_int")

namespace CoupledField
{

LatticeBoltzmannInt::LatticeBoltzmannInt(BaseMaterial* matData, SubTensorType type,  bool coordUpdate)
  : BaseForm(matData, type, coordUpdate)
{
  name_ = "LatticeBoltzmannInt";
  baseType_ = STIFFNESS;
}

void LatticeBoltzmannInt::CalcElementMatrix(Matrix<Double>& elemMat, EntityIterator& ent1, EntityIterator& ent2)
{
   elemMat.Resize(getNrDofs(), getNrDofs());
   elemMat.InitValue(0.0);
   // this is dummy data, it will be set in LatticeBoltzmannPDE as the LBM problem is solved afterwards
   for(unsigned int i = 0; i < getNrDofs(); i++)
     elemMat[i][i] = 1.0;
}

} // end of namespace

