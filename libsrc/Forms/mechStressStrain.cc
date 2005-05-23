#include <iostream>
#include <fstream>

#include "mechStressStrain.hh"

namespace CoupledField
{

  // =============================================================================
  // base class for mechanical stress computation
  // =============================================================================

  // base class for calculation of mechanical stresses
  MechStressStrain::MechStressStrain(BaseFE * aptelem, MaterialData & matData) 
    : linElastInt(aptelem, matData)

  {
    ENTER_FCN( "MechStressStrain::MechStressStrain" );
  }



  MechStressStrain::MechStressStrain(MaterialData & matData) 
    : linElastInt(matData)
  {
    ENTER_FCN( "MechStressStrain::MechStressStrain" );
  }
 

  MechStressStrain::~MechStressStrain()
  {
    ENTER_FCN( "MechStressStrain::~MechStressStrain" );
  }


  /// calculates of stresses T (vector notation)
  // T = c . S with c the material tensor snd S the linear strains
  // S = B_lin * u 
  // see Habil. M. Kaltenbacher 
  void MechStressStrain::
  CalcStressVec(Vector<Double>& stressVec, Integer ip, Matrix<Double> & ptCoord)
  {
    ENTER_FCN( "MechStressStrain::calcPiolaStressTensor" );

    Matrix<Double> dMat;
    calcDMat(dMat);
  
    // convert displacement of all elem nodes into one vector: 
    // (uNode1X, uNode1Y, uNode2X, uNode2Y, ...)
    Vector<Double> displVec;
    elemDisp_.ConvertToVec_AppendCols(displVec);

    // linear differential operator B_lin
    Matrix<Double> linBMat;    
    calcBMat( linBMat, ip, ptCoord);

    Vector<Double> linStrain(linBMat * displVec );
  
    stressVec = dMat * linStrain;
  }


  // returns linear B - matrix
  void MechStressStrain::
  calcBMat(Matrix<Double> & bMat, Integer ip, Matrix<Double> & ptCoord)
  {
    ENTER_FCN( "MechStressStrain::calcLinBMat" );

    // linear differential operator B_lin
    linElastInt::calcBMat(bMat, ip, ptCoord);
  }



  // =============================================================================
  // class for 3d stresses
  // =============================================================================
  MechStressStrain3D::MechStressStrain3D(BaseFE * aptelem, MaterialData & matData) 
    :MechStressStrain(aptelem, matData)
  {
    ENTER_FCN( "mechStressStrain3D::mechStressStrain3D" );
  }

  MechStressStrain3D::MechStressStrain3D(MaterialData & matData) 
    :MechStressStrain(matData)
  {
    ENTER_FCN( "mechStressStrain3D::mechStressStrain3D");
  }
 

  MechStressStrain3D::~MechStressStrain3D()
  {
    ENTER_FCN( "mechStressStrain3D::~mechStressStrain3D" );

  }

  void MechStressStrain3D::calcDMat(Matrix<Double> & dMat)
  {
    ENTER_FCN( "MechStressStrain3D::calcMaterialDMat" );

    linElastInt::Calc3DMaterialMat(dMat);
  }


  // =============================================================================
  // 2D axi case
  // =============================================================================

  MechStressStrainAxi::MechStressStrainAxi(BaseFE * aptelem, MaterialData & matData) 
    : MechStressStrain(aptelem, matData)
  {
    ENTER_FCN( "MechStressStrainAxi::MechStressStrainAxi" );

    isaxi_ = TRUE;
  }

  MechStressStrainAxi::MechStressStrainAxi(MaterialData & matData) 
    : MechStressStrain(matData)
  {
    ENTER_FCN( "MechStressStrainAxi::MechStressStrainAxi" );
    isaxi_ = TRUE;
  }

  MechStressStrainAxi::~MechStressStrainAxi()
  {
    ENTER_FCN( "MechStressStrainAxi::~MechStressStrainAxi" );
  }

  void MechStressStrainAxi::calcDMat(Matrix<Double> & dMat)
  {
    ENTER_FCN( "MechStressStrainAxi::calcMaterialDMat" );

    linElastInt::CalcAxiMaterialMat(dMat,actOrientation);
  }


  // ===================================================================================
  // plane strain case
  // ===================================================================================

  MechStressStrainPlaneStrain::MechStressStrainPlaneStrain(BaseFE * aptelem, MaterialData & matData) 
    : MechStressStrain(aptelem, matData)
  {
    ENTER_FCN( "MechStressStrainPlaneStrain::MechStressStrainPlaneStrain" );
  }


  MechStressStrainPlaneStrain::MechStressStrainPlaneStrain(MaterialData & matData) 
    : MechStressStrain(matData)
  {
    ENTER_FCN( "MechStressStrainPlaneStrain::MechStressStrainPlaneStrain" );
  }
 
  MechStressStrainPlaneStrain::~MechStressStrainPlaneStrain()
  {
    ENTER_FCN( "MechStressStrainPlaneStrain::~MechStressStrainPlaneStrain" );
  }

  void MechStressStrainPlaneStrain::calcDMat(Matrix<Double> & dMat)
  {
    ENTER_FCN( "MechStressStrainPlaneStrain::calcMaterialDMat" );

    linElastInt::CalcPlaneStrainMaterialMat(dMat,actOrientation);
  }


} // end namespace CoupledField
