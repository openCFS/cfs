#include "Forms/linElastInt.hh"
#include "General/environment.hh"
#include "MatVec/vector.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"
#include "Elements/basefe.hh"
#include "Materials/baseMaterial.hh" // IWYU pragma: keep
#include "fluidMechShearStress.hh"

namespace CoupledField {
class EntityIterator;
}  // namespace CoupledField

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
    
    matData->GetScalar(density_,DENSITY,Global::REAL);
    matData->GetScalar(dynamicViscosity_,DYNAMIC_VISCOSITY,Global::REAL);
    matData->GetScalar(kinematicViscosity_,KINEMATIC_VISCOSITY,Global::REAL);

  }
 
  template <class TYPE>
  FluidMechShearStress<TYPE>::~FluidMechShearStress()
  {
  }


  template <class TYPE>
  void FluidMechShearStress<TYPE>::CalcShearStressVec( \
                                   Matrix<TYPE>& viscStress , \
                                   EntityIterator& volEnt, \
                                   const EntityIterator& surfEnt, \
                                   const Vector<Double> outNormal)
  {
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( volEnt );
    const Elem* volElem = volEnt.GetElem();
    const Elem* surfElem = surfEnt.GetElem();

    Grid* ptGrid = domain->GetGrid();
    const UInt spaceDim = ptGrid->GetDim();
    if (spaceDim != 2)
    {
      EXCEPTION("Surface viscous stress only implemented for 2D");
    }
    Double jacDet, auxFactor;

    Vector<Double> *intPoints;
    Vector<Double> intWeights;
    Vector<Double> point;
    Matrix<Double> derivShpFncVolume;
    Vector<Double> DxXi, DyXi, xi;
    Vector<Double> DxXi_vel, DyXi_vel;
    Matrix<Double> globalIp;
    Matrix<Double> localVolIp;
    Matrix<Double> viscStressTensor;
    Vector<Double> viscStressVector;
    Matrix<Double> viscStressTmp;

    const UInt numIntPoints = surfElem->ptElem->GetNumIntPoints();

    const UInt nrVolFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const UInt nrSurfFncs = surfElem->ptElem->GetNumFncs( ansatzFct1_ );
    viscStress.Resize(nrSurfFncs, spaceDim);
    viscStress.Init();
    viscStressTmp.Resize(nrSurfFncs, spaceDim);
    viscStressTensor.Resize(spaceDim, spaceDim);
    viscStressVector.Resize(spaceDim);
    DxXi.Resize(nrVolFncs);
    DyXi.Resize(nrVolFncs);
    DxXi_vel.Resize(spaceDim);
    DyXi_vel.Resize(spaceDim);
    xi.Resize(nrSurfFncs);

    intPoints = surfElem->ptElem->GetIntPoints();
    intWeights = surfElem->ptElem->GetIntWeights();

    globalIp.Resize(spaceDim, numIntPoints);
    localVolIp.Resize(spaceDim, numIntPoints);
    point.Resize(spaceDim);

    Matrix<Double> surfElemCoords;
    ptGrid->GetElemNodesCoord( surfElemCoords, \
                               surfElem->connect, \
                               true );

    for (UInt actIntPt = 0; actIntPt < numIntPoints; ++actIntPt)
    {
      surfElem->ptElem->Local2GlobalCoord(point, intPoints[actIntPt], \
                                          surfElemCoords, NULL);

      for (UInt j = 0; j < spaceDim; ++j)
      {
        globalIp[j][actIntPt] = point[j];
      }
    }
    try {
      // transform global to local coordinates of volume element
      volElem->ptElem->Global2LocalCoords(localVolIp, globalIp, \
                                            ptCoord_);
    } catch (Exception& ex) {
      EXCEPTION("global to local threw an error in Elem: "
          << volElem->elemNum);
    }

    jacDet = surfElem->ptElem->CalcJacobianDet(intPoints[0], \
                                     surfElemCoords, surfElem);
    for (UInt actIntPt = 0; actIntPt < numIntPoints; ++actIntPt)
    {
      for (UInt j = 0; j < spaceDim; ++j)
      {
        point[j] = localVolIp[j][actIntPt];
      }
      try {
        volElem->ptElem->GetGlobDerivShFnc(derivShpFncVolume, \
                                           point, ptCoord_, volElem);
      } catch (Exception& ex) {
        WARN(ex.GetMsg() << "\n volume element: " << volElem->elemNum
             << "\tsurface element: " << surfElem->elemNum);
      }
      surfElem->ptElem->GetShFncAtIp(xi, actIntPt +1, surfElem);

      for (UInt i = 0; i < nrVolFncs; ++i)
      {
        DxXi[i] = derivShpFncVolume[i][0];
        DyXi[i] = derivShpFncVolume[i][1];
      }

      DxXi_vel = elemVelo_ * DxXi;
      DyXi_vel = elemVelo_ * DyXi;

      // only 2D !
      viscStressTensor[0][0] = 2 * DxXi_vel[0];
      viscStressTensor[0][1] = DxXi_vel[1] + DyXi_vel[0];
      viscStressTensor[1][0] = DyXi_vel[0] + DxXi_vel[1];
      viscStressTensor[1][1] = 2 * DyXi_vel[1];
      viscStressVector = viscStressTensor * outNormal;

      auxFactor = intWeights[actIntPt] * jacDet * dynamicViscosity_;
      viscStressVector *= auxFactor;

      for (UInt i = 0; i < nrSurfFncs; ++i)
      {
        for (UInt j = 0; j < spaceDim; ++j)
        {
          viscStressTmp[i][j]  = viscStressVector[j];
          viscStressTmp[i][j] *= xi[i];
        }
      }
      viscStress += viscStressTmp;
    }
  }

  template <class TYPE>
  void FluidMechShearStress<TYPE>::CalcPressureForce( \
                                   Matrix<TYPE>& pressureForce, \
                                   EntityIterator& volEnt, \
                                   const EntityIterator& surfEnt, \
                                   const Vector<Double> outNormal)
  {
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( volEnt );
    const Elem* volElem = volEnt.GetElem();
    const Elem* surfElem = surfEnt.GetElem();

    Grid* ptGrid = domain->GetGrid();
    const UInt spaceDim = ptGrid->GetDim();
    if (spaceDim != 2)
    {
      EXCEPTION("Surface viscous stress only implemented for 2D");
    }

    Matrix<Double> pressureForceTmp;
    Vector<Double> *intPoints;
    Vector<Double> intWeights;
    Vector<Double> point;
    Vector<Double> shpFncVolume;
    Vector<Double> xi;
    Matrix<Double> globalIp;
    Matrix<Double> localVolIp;
    Vector<Double> pressVector;
    Double pressScalar;
    Double jacDet;

    const UInt numIntPoints = surfElem->ptElem->GetNumIntPoints();

    //const UInt nrVolFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const UInt nrSurfFncs = surfElem->ptElem->GetNumFncs( ansatzFct1_ );
    pressureForce.Resize(nrSurfFncs, spaceDim);
    pressureForce.Init();
    pressureForceTmp.Resize(nrSurfFncs, spaceDim);
    pressVector.Resize(spaceDim);
    xi.Resize(nrSurfFncs);

    intPoints = surfElem->ptElem->GetIntPoints();
    intWeights = surfElem->ptElem->GetIntWeights();

    globalIp.Resize(spaceDim, numIntPoints);
    localVolIp.Resize(spaceDim, numIntPoints);
    point.Resize(spaceDim);

    Matrix<Double> surfElemCoords;
    ptGrid->GetElemNodesCoord( surfElemCoords, \
                               surfElem->connect, \
                               true );

    jacDet = surfElem->ptElem->CalcJacobianDetAtIp(1, \
                                      surfElemCoords, surfElem );
    for (UInt actIntPt = 0; actIntPt < numIntPoints; ++actIntPt)
    {
      surfElem->ptElem->Local2GlobalCoord(point, intPoints[actIntPt], \
                                          surfElemCoords, NULL);

      for (UInt j = 0; j < spaceDim; ++j)
      {
        globalIp[j][actIntPt] = point[j];
      }
    }
    try {
      // transform global to local coordinates of volume element
      volElem->ptElem->Global2LocalCoords(localVolIp, globalIp, \
                                            ptCoord_);
    } catch (Exception& ex) {
      EXCEPTION("global to local threw an error in Elem: "
          << volElem->elemNum);
    }

    for (UInt actIntPt = 0; actIntPt < numIntPoints; ++actIntPt)
    {
      for (UInt j = 0; j < spaceDim; ++j)
      {
        point[j] = localVolIp[j][actIntPt];
      }
      try {
        volElem->ptElem->GetShFnc(shpFncVolume, point, volElem);
      } catch (Exception& ex) {
        WARN(ex.GetMsg() << "\n volume element: " << volElem->elemNum
             << "\tsurface element: " << surfElem->elemNum);
      }
      surfElem->ptElem->GetShFncAtIp(xi, actIntPt +1, surfElem);

      pressScalar  = elemPres_ * shpFncVolume;
      pressScalar *= -1.0 * density_;
      pressScalar *= intWeights[actIntPt];
      pressScalar *= jacDet;

      // only 2D !
      pressVector[0] = pressScalar * outNormal[0];
      pressVector[1] = pressScalar * outNormal[1];

      for (UInt i = 0; i < nrSurfFncs; ++i)
      {
        for (UInt j = 0; j < spaceDim; ++j)
        {
          pressureForceTmp[i][j]  = pressVector[j];
          pressureForceTmp[i][j] *= xi[i];
        }
      }
      pressureForce += pressureForceTmp;
    }
  }

#ifdef __GNUC__
  // Explicite template instantiation
  template class FluidMechShearStress<Double>;
  //template class FluidMechShearStress<Complex>;
#endif

} // end namespace CoupledField
