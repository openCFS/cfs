#ifndef PIEZOSTRESSSTRAIN_HH_
#define PIEZOSTRESSSTRAIN_HH_

#include "Forms/mechStressStrain.hh"
#include "Forms/gradfieldop.hh"

namespace CoupledField
{
class BasePairCoupling;
class BaseMaterial;

//! class for calculation of mechanical stresses and strains
template <class TYPE>
class PiezoStressStrain : public MechStressStrain<TYPE>
{
public:
  /** @param coupleMat the PiecoCoupling material, not the mech material */
  PiezoStressStrain(BaseMaterial* coupleMat, SubTensorType type, BasePairCoupling* bpc);

  virtual ~PiezoStressStrain();

  /** overwrites MechStressStrain::CalcStressVec() and calculates the coupled stress  */
  void CalcStressVec(Vector<TYPE>& stressVec, UInt ip, EntityIterator& ent);

private:

  BasePairCoupling* bpc_;
  BaseMaterial*     coupleMat_;
  GradientFieldOp<TYPE>* gradOp;

  Vector<TYPE> grad;
  Vector<TYPE> tmpCouplStress;
  Vector<double> pont_location_; // the elec grad field operator wants the location
  Matrix<double> piezoCouplingMatT; // we consider only real valued material!
};

} // end namespace

#endif /* PIEZOSTRESSSTRAIN_HH_ */
