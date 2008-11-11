#ifndef FILE_FLUIDMECHOPTMASSINT
#define FILE_FLUIDMECHOPTMASSINT

#include "fluidMechInt.hh"

namespace CoupledField
{

//**************************************************************************************************************************
//***************** Stabilized FEM ****************************************************************************************
//**************************************************************************************************************************

  /// Class for calculation the contribution of the convective term in non conservative form to the element stiffness matrix
  /// in 2D Plane
  class FluidMechPlaneMassInt_UV : public FluidMechInt {
public:
  FluidMechPlaneMassInt_UV(Double density, Double kinematicViscosity,Matrix<Double> & stabilParams);
  virtual ~FluidMechPlaneMassInt_UV();
  void CalcElementMatrix( Matrix<Double>& elemMat, EntityIterator& ent1, EntityIterator& ent2 );
protected:    
  Matrix<Double> & stabilParams_;
};
//***************************************************
class FluidMechPlaneMassInt_UQ : public FluidMechInt {
public:
  FluidMechPlaneMassInt_UQ(Double density,Double kinematicViscosity,Matrix<Double> & stabilParams);
  virtual ~FluidMechPlaneMassInt_UQ();
  void CalcElementMatrix( Matrix<Double>& elemMat, EntityIterator& ent1, EntityIterator& ent2 );
protected:    
  Matrix<Double> & stabilParams_;
};


//**************************************************************************************************************************
//***************** mixed FEM **********************************************************************************************
//**************************************************************************************************************************
class FluidMechPlaneMixedMassInt_UV : public FluidMechInt {
public:
FluidMechPlaneMixedMassInt_UV(Double density, Double kinematicViscosity,Matrix<Double> & stabilParams);
virtual ~FluidMechPlaneMixedMassInt_UV();
void CalcElementMatrix( Matrix<Double>& elemMat, EntityIterator& ent1, EntityIterator& ent2 );
protected:    
Matrix<Double> & stabilParams_;
};


}

#endif // FILE_FLUIDMECHOPTMASSINT
