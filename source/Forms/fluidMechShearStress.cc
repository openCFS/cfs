#include <iostream>
#include <fstream>

#include "fluidMechShearStress.hh"

namespace CoupledField
{
  
  // =============================================================================
  // base class for fluid dynamical shear stress computation
  // =============================================================================

  // base class for calculation of shear stresses
  template <class TYPE>
  FluidMechShearStress<TYPE>::FluidMechShearStress(BaseMaterial* matData, 
                                                   SubTensorType type) 
    : linElastInt(matData, type)
  {
    name_="FluidMechShearStress";

    coordUpdate_ = true;
    
    matData->GetScalar(density_,DENSITY,REAL);
    matData->GetScalar(dynamicViscosity_,DYNAMIC_VISCOSITY,REAL);
    matData->GetScalar(kinematicViscosity_,KINEMATIC_VISCOSITY,REAL);

  }
 
  template <class TYPE>
  FluidMechShearStress<TYPE>::~FluidMechShearStress()
  {
  }


  /// calculates shear stresses T in vector notation
  // T = mu * d v_i / d x_j with mu the dynamic viscosity

  template <class TYPE>
  void FluidMechShearStress<TYPE>::
  CalcShearStressVec(Vector<TYPE>& stressVec, UInt ip, EntityIterator& ent)
  {

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );

    // convert velocities of all elem nodes into one vector: 
    // (vNode1X, vNode1Y, vNode2X, vNode2Y, ...)
    Vector<TYPE> veloVec;
    elemVelo_.ConvertToVec_AppendCols(veloVec);

    // linear differential operator B_lin
    Matrix<Double> linBMat;    
    calcBMat( linBMat, ip, ptCoord_);

    Vector<TYPE> linVelocityStrain(linBMat.GetSizeRow());
    linVelocityStrain.Init();
    stressVec.Init();
    linVelocityStrain = linBMat * veloVec;
    stressVec = linVelocityStrain*dynamicViscosity_;
  }

  // returns linear B - matrix

  template <class TYPE>
  void FluidMechShearStress<TYPE>::
  calcBMat(Matrix<Double> & bMat, UInt ip, Matrix<Double> & ptCoord)
  {

    // linear differential operator B_lin
    linElastInt::calcBMat(bMat, ip, ptCoord);
  }




#ifdef __GNUC__
  // Explicite template instantiation
  template class FluidMechShearStress<Double>;
  template class FluidMechShearStress<Complex>;
#endif

} // end namespace CoupledField
