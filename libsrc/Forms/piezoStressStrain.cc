#include <iostream>
#include <fstream>

#include "piezoStressStrain.hh"

namespace CoupledField
{

  // =============================================================================
  // base class for mechanical stress computation
  // =============================================================================

  // base class for calculation of mechanical stresses
  PiezoStressStrain::PiezoStressStrain(BaseFE * aptelem, MaterialData & matData) 
    : linPiezoInt(aptelem, matData)

  {
    ENTER_FCN( "PiezoStressStrain::PiezoStressStrain" );
  }



   PiezoStressStrain::PiezoStressStrain(MaterialData & matData) 
    : linPiezoInt(matData)
  {
    ENTER_FCN( "PiezoStressStrain::PiezoStressStrain" );
  }
 

   PiezoStressStrain::~PiezoStressStrain()
  {
    ENTER_FCN( "PiezoStressStrain::~PiezoStressStrain" );
  }


  /// calculates of stresses T (vector notation)
  // T = c . S - e^T E with c tensor of mechanical moduli and e the piezoelectric tensor
  // S = Bmech * u 
  // E = Belec V
  // see Habil. M. Kaltenbacher 
  void PiezoStressStrain::
  CalcStressVec(Vector<Double>& stressElecVec, Integer ip, Matrix<Double> & ptCoord)
  {
    ENTER_FCN( "PiezoStressStrain::CalcStressVec" );

    Matrix<Double> dMat;
    calcDMat(dMat);
 
    // convert displacement of all elem nodes into one vector: 
    // (uNode1X, uNode1Y, VNode1, uNode2X, uNode2Y,VNode2,  ...)
    Vector<Double> solVec;
    elemDisp_.ConvertToVec_AppendCols(solVec);

    // linear differential operator B_lin
    Matrix<Double> linBMat;    
    calcBMat( linBMat, ip, ptCoord);

    Vector<Double> linStrainElec(linBMat * solVec );

    // | c Bmech u - e^T Belec V |
    // | e Bmech u + eps Belec V |
    //    Vector<Double> stressElecVec = dMat * linStrainElec;
    stressElecVec = dMat * linStrainElec;

    //get out the stresses!
    //    for (Integer comp=0; comp < getStressComp(); comp++)
    //      stressVec[comp] = stressElecVec[comp];
  }


  // returns linear B - matrix
  void PiezoStressStrain::
  calcBMat(Matrix<Double> & bMat, Integer ip, Matrix<Double> & ptCoord)
  {
    ENTER_FCN( "PiezoStressStrain::calcLinBMat" );

    // linear differential operator B_lin
    linPiezoInt::calcBMat(bMat, ip, ptCoord);
  }



  // =============================================================================
  // class for 3d stresses
  // =============================================================================
  PiezoStressStrain3D::PiezoStressStrain3D(BaseFE * aptelem, MaterialData & matData) 
    :PiezoStressStrain(aptelem, matData)
  {
    ENTER_FCN( "mechStressStrain3D::mechStressStrain3D" );
  }

  PiezoStressStrain3D::PiezoStressStrain3D(MaterialData & matData) 
    :PiezoStressStrain(matData)
  {
    ENTER_FCN( "mechStressStrain3D::mechStressStrain3D");
  }
 

  PiezoStressStrain3D::~PiezoStressStrain3D()
  {
    ENTER_FCN( "mechStressStrain3D::~mechStressStrain3D" );

  }

  void PiezoStressStrain3D::calcDMat(Matrix<Double> & dMat)
  {
    ENTER_FCN( "PiezoStressStrain3D::calcDMat" );

    linPiezoInt::Calc3DMaterialMat(dMat);
  }


  // =============================================================================
  // 2D axi case
  // =============================================================================

  PiezoStressStrainAxi::PiezoStressStrainAxi(BaseFE * aptelem, MaterialData & matData) 
    : PiezoStressStrain(aptelem, matData)
  {
    ENTER_FCN( "PiezoStressStrainAxi::PiezoStressStrainAxi" );

    isaxi_ = TRUE;
  }

  PiezoStressStrainAxi::PiezoStressStrainAxi(MaterialData & matData) 
    : PiezoStressStrain(matData)
  {
    ENTER_FCN( "PiezoStressStrainAxi::PiezoStressStrainAxi" );
    isaxi_ = TRUE;
  }

  PiezoStressStrainAxi::~PiezoStressStrainAxi()
  {
    ENTER_FCN( "PiezoStressStrainAxi::~PiezoStressStrainAxi" );
  }

  void PiezoStressStrainAxi::calcDMat(Matrix<Double> & dMat)
  {
    ENTER_FCN( "PiezoStressStrainAxi::calcMaterialDMat" );

    linPiezoInt::CalcAxiMaterialMat(dMat);
  }


  // ===================================================================================
  // plane strain case
  // ===================================================================================

  PiezoStressStrainPlaneStrain::PiezoStressStrainPlaneStrain(BaseFE * aptelem, MaterialData & matData) 
    : PiezoStressStrain(aptelem, matData)
  {
    ENTER_FCN( "PiezoStressStrainPlaneStrain::PiezoStressStrainPlaneStrain" );
  }


  PiezoStressStrainPlaneStrain::PiezoStressStrainPlaneStrain(MaterialData & matData) 
    : PiezoStressStrain(matData)
  {
    ENTER_FCN( "PiezoStressStrainPlaneStrain::PiezoStressStrainPlaneStrain" );
  }
 
  PiezoStressStrainPlaneStrain::~PiezoStressStrainPlaneStrain()
  {
    ENTER_FCN( "PiezoStressStrainPlaneStrain::~PiezoStressStrainPlaneStrain" );
  }

  void PiezoStressStrainPlaneStrain::calcDMat(Matrix<Double> & dMat)
  {
    ENTER_FCN( "PiezoStressStrainPlaneStrain::calcMaterialDMat" );

    linPiezoInt::CalcPlaneStrainMaterialMat(dMat);
  }


} // end namespace CoupledField
