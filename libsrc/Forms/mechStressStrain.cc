#include <iostream>
#include <fstream>

#include "mechStressStrain.hh"

namespace CoupledField
{

  // =============================================================================
  // base class for mechanical stress computation
  // =============================================================================

  // base class for calculation of mechanical stresses
  template <class TYPE>
  MechStressStrain<TYPE>::MechStressStrain(BaseFE * aptelem, MaterialData & matData) 
    : linElastInt(aptelem, matData)

  {
    ENTER_FCN( "MechStressStrain::MechStressStrain" );
  }


  template <class TYPE>
  MechStressStrain<TYPE>::MechStressStrain(MaterialData & matData) 
    : linElastInt(matData)
  {
    ENTER_FCN( "MechStressStrain::MechStressStrain" );
  }
 
  template <class TYPE>
  MechStressStrain<TYPE>::~MechStressStrain()
  {
    ENTER_FCN( "MechStressStrain::~MechStressStrain" );
  }


  /// calculates of stresses T (vector notation)
  // T = c . S with c the material tensor snd S the linear strains
  // S = B_lin * u 
  // see Habil. M. Kaltenbacher 

  template <class TYPE>
  void MechStressStrain<TYPE>::
  CalcStressVec(Vector<TYPE>& stressVec, UInt ip, Matrix<Double> & ptCoord)
  {
    ENTER_FCN( "MechStressStrain::calcPiolaStressTensor" );

    Matrix<Double> dMat;
    calcDMat(dMat);
  
    // convert displacement of all elem nodes into one vector: 
    // (uNode1X, uNode1Y, uNode2X, uNode2Y, ...)
    Vector<TYPE> displVec;
    elemDisp_.ConvertToVec_AppendCols(displVec);

    // linear differential operator B_lin
    Matrix<Double> linBMat;    
    calcBMat( linBMat, ip, ptCoord);

    Vector<TYPE> linStrain(linBMat.GetSizeRow());
    //Matrix<TYPE>(linBMat).Mult(displVec,linStrain);
    //Matrix<TYPE>(dMat).Mult(linStrain,stressVec);
    linStrain = linBMat * displVec;
    stressVec = dMat * linStrain;

  }

  /// calculates green-lagrangian strains (linear part, vector notation)
  // see Habil. M. Kaltenbacher 

  template <class TYPE>
  void MechStressStrain<TYPE>::
  CalcStrainVec(Vector<TYPE>& strainVec, UInt ip, Matrix<Double> & ptCoord)
  {
    ENTER_FCN( "MechStressStrain::CalcStrainVec" );

    Matrix<Double> dMat;
    calcDMat(dMat);
  
    // convert displacement of all elem nodes into one vector: 
    // (uNode1X, uNode1Y, uNode2X, uNode2Y, ...)
    Vector<TYPE> displVec;
    elemDisp_.ConvertToVec_AppendCols(displVec);

    // linear differential operator B_lin
    Matrix<Double> linBMat;    
    calcBMat( linBMat, ip, ptCoord);

    //strainVec = linBMat * displVec;
    strainVec.Resize(linBMat.GetSizeRow());
    Matrix<TYPE>(linBMat).Mult(displVec,strainVec);
  
  }


  // returns linear B - matrix

  template <class TYPE>
  void MechStressStrain<TYPE>::
  calcBMat(Matrix<Double> & bMat, UInt ip, Matrix<Double> & ptCoord)
  {
    ENTER_FCN( "MechStressStrain::calcLinBMat" );

    // linear differential operator B_lin
    linElastInt::calcBMat(bMat, ip, ptCoord);
  }



  // =============================================================================
  // class for 3d stresses
  // =============================================================================

  template <class TYPE> 
  MechStressStrain3D<TYPE>::MechStressStrain3D(BaseFE * aptelem, MaterialData & matData) 
    : MechStressStrain<TYPE>(aptelem, matData)
  {
    ENTER_FCN( "mechStressStrain3D::mechStressStrain3D" );
  }

  template <class TYPE>
  MechStressStrain3D<TYPE>::MechStressStrain3D(MaterialData & matData) 
    :MechStressStrain<TYPE>(matData)
  {
    ENTER_FCN( "mechStressStrain3D::mechStressStrain3D");
  }
 
  template <class TYPE>
  MechStressStrain3D<TYPE>::~MechStressStrain3D()
  {
    ENTER_FCN( "mechStressStrain3D::~mechStressStrain3D" );

  }

  template <class TYPE>
  void MechStressStrain3D<TYPE>::calcDMat(Matrix<Double> & dMat)
  {
    ENTER_FCN( "MechStressStrain3D::calcMaterialDMat" );

    linElastInt::Calc3DMaterialMat(dMat);
  }


  // =============================================================================
  // 2D axi case
  // =============================================================================

  template <class TYPE>
  MechStressStrainAxi<TYPE>::MechStressStrainAxi(BaseFE * aptelem, MaterialData & matData) 
    : MechStressStrain<TYPE>(aptelem, matData)
  {
    ENTER_FCN( "MechStressStrainAxi::MechStressStrainAxi" );

    isaxi_ = TRUE;
  }

  template <class TYPE>
  MechStressStrainAxi<TYPE>::MechStressStrainAxi(MaterialData & matData) 
    : MechStressStrain<TYPE>(matData)
  {
    ENTER_FCN( "MechStressStrainAxi::MechStressStrainAxi" );
    isaxi_ = TRUE;
  }

  template <class TYPE>
  MechStressStrainAxi<TYPE>::~MechStressStrainAxi()
  {
    ENTER_FCN( "MechStressStrainAxi::~MechStressStrainAxi" );
  }

  template <class TYPE>
  void MechStressStrainAxi<TYPE>::calcDMat(Matrix<Double> & dMat)
  {
    ENTER_FCN( "MechStressStrainAxi::calcMaterialDMat" );

    linElastInt::CalcAxiMaterialMat(dMat,actOrientation);
  }


  // ===================================================================================
  // plane strain case
  // ===================================================================================

  template <class TYPE> MechStressStrainPlaneStrain<TYPE>::
  MechStressStrainPlaneStrain(BaseFE * aptelem, MaterialData & matData) 
    : MechStressStrain<TYPE>(aptelem, matData)
  {
    ENTER_FCN( "MechStressStrainPlaneStrain::MechStressStrainPlaneStrain" );
  }

  template <class TYPE>  MechStressStrainPlaneStrain<TYPE>::
  MechStressStrainPlaneStrain(MaterialData & matData) 
    : MechStressStrain<TYPE>(matData)
  {
    ENTER_FCN( "MechStressStrainPlaneStrain::MechStressStrainPlaneStrain" );
  }

  template <class TYPE>
  MechStressStrainPlaneStrain<TYPE>::~MechStressStrainPlaneStrain()
  {
    ENTER_FCN( "MechStressStrainPlaneStrain::~MechStressStrainPlaneStrain" );
  }

  template <class TYPE>
  void MechStressStrainPlaneStrain<TYPE>::calcDMat(Matrix<Double> & dMat)
  {
    ENTER_FCN( "MechStressStrainPlaneStrain::calcMaterialDMat" );

    linElastInt::CalcPlaneStrainMaterialMat(dMat,actOrientation);
  }


} // end namespace CoupledField
